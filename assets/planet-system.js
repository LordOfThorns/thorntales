/* planet-system.js — "Planetary system" block of the map generator (terrain.html).
   PROTOTYPE, in development.

   The reader picks the star (and an optional close companion or a distant one), the
   distance of the planet and the shape of its orbit, its axial tilt, how fast it turns,
   its atmosphere and whether it is tidally locked. Here that becomes what terrain-map.js
   needs (PlanetSystem.planet()):
     flux   = year-mean sunlight on the planet, Earth = 1
              (luminosity / distance^2, more on a stretched orbit: / sqrt(1 - e^2),
              plus the light of a distant companion)
     tilt   = axial tilt in degrees           cells = wind cells per hemisphere (1, 3, 5)
     atmo   = thin | earth | thick | venus    locked = one side always faces the star
     ecc    = eccentricity of the orbit (for the seasons)
   Things that need no new map go through PlanetSystem.onLive(cb): the moment of the year
   (PlanetSystem.season(): solar longitude in radians, or null for the whole year) and the
   colour of the starlight (PlanetSystem.tint()).

   The orrery shows the star(s), the habitable zone, the orbit (an ellipse with the star in
   its focus) and the planet with its axis (it keeps its direction in space all year: that
   is what makes seasons) and its air. Dragging on the orrery moves the planet closer to the
   star or further away. Distances are drawn on a log scale, otherwise a red dwarf's planet
   and a blue star's planet would not both fit.
*/
(function () {
  'use strict';

  const $ = (id) => document.getElementById(id);
  const box = $('planetBox');
  const orrery = $('orrery');
  if (!box || !orrery) return;

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);

  // Main sequence stars, rough: luminosity and mass in Suns; tint = colour of the light
  // on the map (multiplied over it), null for a Sun-like star
  const STARS = {
    red:    { L: 0.04, M: 0.3,  color: '#ff6a48', r: 5,    tint: [255, 176, 140], k: 0.32 },
    orange: { L: 0.35, M: 0.75, color: '#ffae5e', r: 6.5,  tint: [255, 214, 170], k: 0.22 },
    yellow: { L: 1,    M: 1,    color: '#ffe58a', r: 7.5,  tint: null, k: 0 },
    white:  { L: 4,    M: 1.5,  color: '#f1f4ff', r: 8.5,  tint: [228, 236, 255], k: 0.18 },
    blue:   { L: 50,   M: 6,    color: '#9cc0ff', r: 10.5, tint: [196, 214, 255], k: 0.3 },
  };
  const D_MIN = 0.02, D_MAX = 20;     // AU, the distance slider (log scale)
  const C_MIN = 5, C_MAX = 300;       // AU, the distant companion
  const SLIDER = 1000;
  const PAIR_SEP = 0.05;              // AU between the stars of a close pair
  const PERIHELION = 283;             // deg of the year where the planet is closest (as Earth)

  const els = {
    star: $('plStar'), star2: $('plStar2'), dist: $('plDist'), distOut: $('plDistOut'),
    ecc: $('plEcc'), eccOut: $('plEccOut'),
    tilt: $('plTilt'), tiltOut: $('plTiltOut'), day: $('plDay'), hz: $('plToHz'),
    atmo: $('plAtmo'), locked: $('plLocked'),
    comp: $('plComp'), compDist: $('plCompDist'), compOut: $('plCompOut'), compBox: $('plCompBox'),
    seasonOn: $('plSeasonOn'), season: $('plSeason'), seasonOut: $('plSeasonOut'),
    play: $('plPlay'), seasonBox: $('plSeasonBox'),
    summary: $('plSummary'),
  };

  const logFrom = (v, lo, hi) => lo * Math.pow(hi / lo, v / SLIDER);
  const logTo = (d, lo, hi) => Math.round(SLIDER * Math.log(d / lo) / Math.log(hi / lo));
  const clampN = (x, lo, hi) => Math.max(lo, Math.min(hi, x));

  function state() {
    const s1 = STARS[els.star.value] || STARS.yellow;
    const s2 = STARS[els.star2.value] || null;
    const comp = STARS[els.comp.value] || null;
    const L = s1.L + (s2 ? s2.L : 0);
    const M = s1.M + (s2 ? s2.M : 0);
    const d = logFrom(Number(els.dist.value), D_MIN, D_MAX);
    const e = clampN(Number(els.ecc.value) / 100, 0, 0.6);
    const cd = logFrom(Number(els.compDist.value), C_MIN, C_MAX);
    const compFlux = comp ? comp.L / (cd * cd) : 0;
    return {
      s1, s2, comp, L, M, d, e, cd, compFlux,
      flux: L / (d * d) / Math.sqrt(1 - e * e) + compFlux,
      tilt: Number(els.tilt.value),
      cells: Number(els.day.value),
      atmo: els.atmo.value || 'earth',
      locked: els.locked.checked,
      year: Math.sqrt(d * d * d / M),            // Kepler, in Earth years
      hzIn: Math.sqrt(L / 1.1), hzOut: Math.sqrt(L / 0.53), // habitable zone, AU
    };
  }

  // What render() gets
  function planet() {
    const s = state();
    // the sliders cannot hit Earth exactly: near enough counts as Earth (the old maps)
    const flux = Math.abs(s.flux - 1) < 0.02 ? 1 : s.flux;
    const tilt = s.tilt === 23 ? 23.44 : s.tilt;
    return { flux, tilt, cells: s.cells, atmo: s.atmo, locked: s.locked, ecc: s.e };
  }

  // The moment of the year shown on the map (radians), or null for the year as a whole
  function season() {
    if (!els.seasonOn.checked || els.locked.checked) return null;
    return Number(els.season.value) * Math.PI / 180;
  }

  // Colour of the starlight over the map: { color: [r, g, b], k } or null
  function tint() {
    const s = state();
    const star = s.s2 && s.s2.L > s.s1.L ? s.s2 : s.s1; // the brighter one colours the light
    return star.tint ? { color: star.tint, k: star.k } : null;
  }

  const live = [];
  let liveQueued = false;
  function fireLive() {
    if (liveQueued) return;
    liveQueued = true;
    requestAnimationFrame(() => { liveQueued = false; for (const cb of live) cb(); });
  }

  const fmt = (x, digits = 2) => {
    const lang = (window.I18n && I18n.lang) || 'ru';
    return x.toLocaleString(lang, { maximumSignificantDigits: digits });
  };

  function seasonName(deg) {
    const a = ((deg % 360) + 360) % 360;
    const key = a < 45 || a >= 315 ? 'spring' : a < 135 ? 'summer' : a < 225 ? 'autumn' : 'winter';
    let text = t('planet.season.' + key);
    const e = Number(els.ecc.value) / 100;
    if (e >= 0.05) {
      const off = Math.abs(((a - PERIHELION + 540) % 360) - 180); // 0 at the nearest point
      if (off < 40) text += t('planet.season.near');
      else if (off > 140) text += t('planet.season.far');
    }
    return text;
  }

  function summary() {
    const s = state();
    els.distOut.textContent = t('planet.au', { d: fmt(s.d) });
    els.tiltOut.textContent = `${s.tilt}°`;
    els.eccOut.textContent = fmt(s.e, 2);
    els.compOut.textContent = t('planet.au', { d: fmt(s.cd) });
    els.compBox.hidden = !s.comp;
    els.seasonBox.classList.toggle('is-off', s.locked);
    els.season.disabled = els.seasonOn.disabled = els.play.disabled = s.locked;
    els.seasonOut.textContent = s.locked ? t('planet.seasonLocked') : seasonName(Number(els.season.value));

    const f = s.flux;
    const warmth = f < 0.55 ? 'frozen' : f < 0.85 ? 'cold' : f < 1.2 ? 'earth' : f < 1.7 ? 'warm' : 'hot';
    const tilt = s.tilt < 10 ? 'none' : s.tilt < 35 ? 'earth' : s.tilt < 54 ? 'strong' : 'inverted';
    const year = s.year < 0.1
      ? t('planet.yearDays', { n: fmt(s.year * 365.25) })
      : t('planet.yearYears', { n: fmt(s.year) });
    const lines = [
      t('planet.light', { x: fmt(f) }) + ' ' + t('planet.warmth.' + warmth),
      t('planet.year', { year }),
    ];
    if (s.locked) lines.push(t('planet.lockedSum'));
    else lines.push(t('planet.tilt.' + tilt));
    lines.push(t('planet.cells', { n: s.cells }));
    if (s.atmo !== 'earth') lines.push(t('planet.atmo.' + s.atmo + 'Sum'));
    if (s.e >= 0.05) lines.push(t('planet.eccSum', { x: fmt(Math.pow((1 + s.e) / (1 - s.e), 2)) }));
    if (s.comp && s.compFlux >= 0.005) lines.push(t('planet.compSum', { x: fmt(s.compFlux) }));
    if (s.s2 && s.d < PAIR_SEP * 4) lines.push(t('planet.warnPair'));
    if (els.star.value === 'red' && s.d < 0.4 && !s.locked) lines.push(t('planet.warnLocked'));
    els.summary.innerHTML = lines.map((l) => `<span>${escapeHtml(l)}</span>`).join('');
  }

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }

  // ---------- Orrery ----------
  const ctx = orrery.getContext('2d');
  let size = 0, dpr = 1;
  function resize() {
    dpr = window.devicePixelRatio || 1;
    size = orrery.clientWidth || 260;
    orrery.width = Math.round(size * dpr);
    orrery.height = Math.round(size * dpr);
  }
  const R0 = 18; // px from the centre where the distance scale starts
  const rOf = (d) => R0 + (size / 2 - 16 - R0) * Math.log(d / D_MIN) / Math.log(D_MAX / D_MIN);
  const dOf = (r) => D_MIN * Math.pow(D_MAX / D_MIN, (r - R0) / (size / 2 - 16 - R0));
  // distance from the star on the orbit at true anomaly nu (periapsis to the right)
  const orbitR = (a, e, nu) => a * (1 - e * e) / (1 + e * Math.cos(nu));

  function glow(x, y, r, color) {
    const g = ctx.createRadialGradient(x, y, 0, x, y, r * 3.2);
    g.addColorStop(0, color);
    g.addColorStop(0.3, color + 'aa');
    g.addColorStop(1, color + '00');
    ctx.fillStyle = g;
    ctx.beginPath(); ctx.arc(x, y, r * 3.2, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = color;
    ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.fill();
  }

  const AIR = { thin: null, earth: [1, 0.3, '150, 200, 255'], thick: [2.4, 0.45, '170, 210, 255'], venus: [4, 0.6, '240, 220, 160'] };

  let nu = 0, spin = 0, pairAngle = 0, last = 0;
  function frame(now) {
    const dt = last ? Math.min(0.05, (now - last) / 1000) : 0;
    last = now;
    if (!size) resize();
    const s = state();
    const lam = season();
    if (lam != null) {
      nu = lam - PERIHELION * Math.PI / 180; // the planet stands where the season says
    } else {
      // one Earth year takes ~10 s on screen; faster near the star on a stretched orbit
      const w = (2 * Math.PI / 10) / clampN(s.year, 0.15, 6);
      nu += dt * w * Math.pow(1 + s.e * Math.cos(nu), 2) / Math.pow(1 - s.e * s.e, 1.5);
    }
    spin += dt * (s.cells === 1 ? 0.6 : s.cells === 5 ? 9 : 3.5);
    pairAngle += dt * 2.2;

    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, size, size);
    const c = size / 2;

    // habitable zone
    const a = rOf(Math.max(D_MIN, s.hzIn)), b = rOf(Math.min(D_MAX, s.hzOut));
    if (b > R0) {
      ctx.beginPath();
      ctx.arc(c, c, Math.max(R0, b), 0, Math.PI * 2);
      ctx.arc(c, c, Math.max(R0, a), 0, Math.PI * 2, true);
      ctx.fillStyle = 'rgba(96, 176, 118, .16)';
      ctx.fill();
      ctx.strokeStyle = 'rgba(120, 200, 140, .35)';
      ctx.setLineDash([3, 4]);
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(c, c, Math.max(R0, a), 0, Math.PI * 2); ctx.stroke();
      ctx.beginPath(); ctx.arc(c, c, Math.max(R0, b), 0, Math.PI * 2); ctx.stroke();
      ctx.setLineDash([]);
    }

    // orbit: an ellipse with the star in its focus (radii on the log scale)
    ctx.strokeStyle = 'rgba(214, 178, 92, .4)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let k = 0; k <= 96; k++) {
      const q = (k / 96) * Math.PI * 2;
      const rr = rOf(clampN(orbitR(s.d, s.e, q), D_MIN, D_MAX * 2));
      const x = c + Math.cos(q) * rr, y = c + Math.sin(q) * rr;
      if (k) ctx.lineTo(x, y); else ctx.moveTo(x, y);
    }
    ctx.stroke();

    // a distant companion: at the edge, up and to the right, with its distance
    if (s.comp) {
      const cx = c + (size / 2 - 14) * Math.cos(-0.8), cy = c + (size / 2 - 14) * Math.sin(-0.8);
      glow(cx, cy, s.comp.r * 0.45, s.comp.color);
    }

    // star(s)
    if (s.s2) {
      const k = s.s2.M / (s.s1.M + s.s2.M); // the lighter star swings wider
      const sep = 9;
      const dx = Math.cos(pairAngle) * sep, dy = Math.sin(pairAngle) * sep;
      glow(c - dx * k, c - dy * k, s.s1.r * 0.8, s.s1.color);
      glow(c + dx * (1 - k), c + dy * (1 - k), s.s2.r * 0.8, s.s2.color);
    } else {
      glow(c, c, s.s1.r, s.s1.color);
    }

    // planet: lit on the star side; the axis keeps its direction all year
    const r = rOf(clampN(orbitR(s.d, s.e, nu), D_MIN, D_MAX * 2));
    const px = c + Math.cos(nu) * r, py = c + Math.sin(nu) * r;
    const pr = 7;
    const air = AIR[s.atmo];
    if (air) {
      ctx.strokeStyle = `rgba(${air[2]}, ${air[1]})`;
      ctx.lineWidth = air[0];
      ctx.beginPath(); ctx.arc(px, py, pr + air[0] / 2 + 0.5, 0, Math.PI * 2); ctx.stroke();
    }
    const sunDir = Math.atan2(c - py, c - px);
    const g = ctx.createLinearGradient(
      px + Math.cos(sunDir) * pr, py + Math.sin(sunDir) * pr,
      px - Math.cos(sunDir) * pr, py - Math.sin(sunDir) * pr);
    g.addColorStop(0, '#7fb2d8');
    g.addColorStop(0.55, '#3b6b8c');
    g.addColorStop(1, '#0d1a26');
    ctx.fillStyle = g;
    ctx.beginPath(); ctx.arc(px, py, pr, 0, Math.PI * 2); ctx.fill();
    if (s.locked) {
      // the same face always to the star: a mark that never leaves the day side
      ctx.fillStyle = 'rgba(214, 178, 92, .95)';
      ctx.beginPath();
      ctx.arc(px + Math.cos(sunDir) * pr * 0.55, py + Math.sin(sunDir) * pr * 0.55, 1.8, 0, Math.PI * 2);
      ctx.fill();
    } else {
      const tilt = s.tilt * Math.PI / 180;
      const ax = Math.sin(tilt), ay = -Math.cos(tilt);
      ctx.strokeStyle = 'rgba(230, 220, 190, .85)';
      ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(px - ax * 12, py - ay * 12); ctx.lineTo(px + ax * 12, py + ay * 12); ctx.stroke();
      // a mark on the equator turns with the day
      const ex = Math.cos(spin) * pr * 0.85;
      ctx.fillStyle = 'rgba(214, 178, 92, .95)';
      ctx.beginPath();
      ctx.arc(px + ex * Math.cos(tilt), py + ex * Math.sin(tilt), 1.6, 0, Math.PI * 2);
      ctx.fill();
    }

    if (visible) requestAnimationFrame(frame); else last = 0;
  }

  // animate only while the orrery is on screen
  let visible = true;
  if (window.IntersectionObserver) {
    new IntersectionObserver((entries) => {
      const was = visible;
      visible = entries[0].isIntersecting;
      if (visible && !was) requestAnimationFrame(frame);
    }).observe(orrery);
  }

  // drag on the orrery: the planet goes where the pointer is (its distance from the star)
  let dragging = false;
  function setFromPointer(e) {
    const rect = orrery.getBoundingClientRect();
    const x = e.clientX - rect.left - rect.width / 2, y = e.clientY - rect.top - rect.height / 2;
    const d = clampN(dOf(Math.hypot(x, y) * size / rect.width), D_MIN, D_MAX);
    els.dist.value = logTo(d, D_MIN, D_MAX);
    if (season() == null) nu = Math.atan2(y, x);
    summary();
  }
  orrery.addEventListener('pointerdown', (e) => {
    if (e.button !== 0) return;
    dragging = true;
    try { orrery.setPointerCapture(e.pointerId); } catch (err) { /* not a real pointer */ }
    setFromPointer(e);
  });
  orrery.addEventListener('pointermove', (e) => { if (dragging) setFromPointer(e); });
  const stop = () => { dragging = false; };
  orrery.addEventListener('pointerup', stop);
  orrery.addEventListener('pointercancel', stop);

  // to the middle of the habitable zone
  els.hz.addEventListener('click', () => {
    const s = state();
    els.dist.value = logTo(clampN(Math.sqrt(s.L), D_MIN, D_MAX), D_MIN, D_MAX);
    summary();
  });

  // the year in motion: the season slider walks round, the map follows
  let playing = 0;
  function setPlaying(on) {
    if (on && !playing) {
      els.seasonOn.checked = true;
      playing = setInterval(() => {
        els.season.value = (Number(els.season.value) + 6) % 360;
        summary();
        fireLive();
      }, 140);
    } else if (!on && playing) {
      clearInterval(playing);
      playing = 0;
    }
    els.play.textContent = t(playing ? 'planet.pause' : 'planet.play');
    els.play.dataset.i18n = playing ? 'planet.pause' : 'planet.play';
  }
  els.play.addEventListener('click', () => { setPlaying(!playing); fireLive(); });

  box.addEventListener('input', (e) => {
    summary();
    // the season and the colour of the light need no new map
    if (e.target === els.season || e.target === els.seasonOn) fireLive();
  });
  box.addEventListener('change', (e) => {
    summary();
    if (e.target === els.locked && els.locked.checked) setPlaying(false);
    if ([els.season, els.seasonOn, els.star, els.star2, els.locked].includes(e.target)) fireLive();
  });
  addEventListener('resize', () => { resize(); });

  for (const el of [els.dist, els.compDist]) el.max = SLIDER;
  els.dist.value = logTo(1, D_MIN, D_MAX);
  els.compDist.value = logTo(50, C_MIN, C_MAX);
  resize();
  summary();
  requestAnimationFrame(frame);

  window.PlanetSystem = {
    planet, season, tint,
    onLive: (cb) => live.push(cb),
    stopYear: () => setPlaying(false),
  };
})();
