/* planet-system.js — "Planetary system" block of the map generator (terrain.html).
   PROTOTYPE, in development.

   The reader picks the star (and an optional close companion), the distance of the planet,
   its axial tilt and how fast it turns. Here that becomes what terrain-map.js needs:
     flux  = sunlight on the planet, Earth = 1  (luminosity of the stars / distance^2)
     tilt  = axial tilt in degrees              (where the sun shines most over a year)
     cells = wind cells per hemisphere          (slow 1, Earth 3, fast 5)
   PlanetSystem.planet() returns { flux, tilt, cells }; terrain-page.js passes it to render.

   The orrery shows the star(s), the habitable zone, the orbit and the planet with its axis
   (it keeps its direction in space all year: that is what makes seasons). Dragging on the
   orrery moves the planet closer to the star or further away. Distances are drawn on a
   log scale, otherwise a red dwarf's planet and a blue star's planet would not both fit.
*/
(function () {
  'use strict';

  const $ = (id) => document.getElementById(id);
  const box = $('planetBox');
  const orrery = $('orrery');
  if (!box || !orrery) return;

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);

  // Main sequence stars, rough: luminosity and mass in Suns
  const STARS = {
    red:    { L: 0.04, M: 0.3,  color: '#ff6a48', r: 5 },
    orange: { L: 0.35, M: 0.75, color: '#ffae5e', r: 6.5 },
    yellow: { L: 1,    M: 1,    color: '#ffe58a', r: 7.5 },
    white:  { L: 4,    M: 1.5,  color: '#f1f4ff', r: 8.5 },
    blue:   { L: 50,   M: 6,    color: '#9cc0ff', r: 10.5 },
  };
  const D_MIN = 0.02, D_MAX = 20;     // AU, the distance slider (log scale)
  const SLIDER = 1000;
  const PAIR_SEP = 0.05;              // AU between the stars of a close pair

  const els = {
    star: $('plStar'), star2: $('plStar2'), dist: $('plDist'), distOut: $('plDistOut'),
    tilt: $('plTilt'), tiltOut: $('plTiltOut'), day: $('plDay'), hz: $('plToHz'),
    summary: $('plSummary'),
  };

  const distFromSlider = (v) => D_MIN * Math.pow(D_MAX / D_MIN, v / SLIDER);
  const sliderFromDist = (d) => Math.round(SLIDER * Math.log(d / D_MIN) / Math.log(D_MAX / D_MIN));

  function state() {
    const s1 = STARS[els.star.value] || STARS.yellow;
    const s2 = STARS[els.star2.value] || null;
    const L = s1.L + (s2 ? s2.L : 0);
    const M = s1.M + (s2 ? s2.M : 0);
    const d = distFromSlider(Number(els.dist.value));
    return {
      s1, s2, L, M, d,
      flux: L / (d * d),
      tilt: Number(els.tilt.value),
      cells: Number(els.day.value),
      year: Math.sqrt(d * d * d / M),            // Kepler, in Earth years
      hzIn: Math.sqrt(L / 1.1), hzOut: Math.sqrt(L / 0.53), // habitable zone, AU
    };
  }

  // What render() gets
  function planet() {
    const s = state();
    // the slider cannot hit Earth exactly: near enough counts as Earth (the old maps)
    const flux = Math.abs(s.flux - 1) < 0.02 ? 1 : s.flux;
    const tilt = s.tilt === 23 ? 23.44 : s.tilt;
    return { flux, tilt, cells: s.cells };
  }

  const fmt = (x, digits = 2) => {
    const lang = (window.I18n && I18n.lang) || 'ru';
    return x.toLocaleString(lang, { maximumSignificantDigits: digits });
  };

  function summary() {
    const s = state();
    els.distOut.textContent = t('planet.au', { d: fmt(s.d) });
    els.tiltOut.textContent = `${s.tilt}°`;

    const f = s.flux;
    const warmth = f < 0.55 ? 'frozen' : f < 0.85 ? 'cold' : f < 1.2 ? 'earth' : f < 1.7 ? 'warm' : 'hot';
    const tilt = s.tilt < 10 ? 'none' : s.tilt < 35 ? 'earth' : s.tilt < 54 ? 'strong' : 'inverted';
    const year = s.year < 0.1
      ? t('planet.yearDays', { n: fmt(s.year * 365.25) })
      : t('planet.yearYears', { n: fmt(s.year) });
    const lines = [
      t('planet.light', { x: fmt(f) }) + ' ' + t('planet.warmth.' + warmth),
      t('planet.year', { year }),
      t('planet.tilt.' + tilt),
      t('planet.cells', { n: s.cells }),
    ];
    if (s.s2 && s.d < PAIR_SEP * 4) lines.push(t('planet.warnPair'));
    if (els.star.value === 'red' && s.d < 0.4) lines.push(t('planet.warnLocked'));
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

  let angle = 0, spin = 0, pairAngle = 0, last = 0;
  function frame(now) {
    const dt = last ? Math.min(0.05, (now - last) / 1000) : 0;
    last = now;
    if (!size) resize();
    const s = state();
    // one Earth year takes ~10 s on screen; very short or long years are slowed / sped up
    angle += dt * (2 * Math.PI / 10) / clampN(s.year, 0.15, 6);
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

    // orbit
    const r = rOf(s.d);
    ctx.strokeStyle = 'rgba(214, 178, 92, .4)';
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(c, c, r, 0, Math.PI * 2); ctx.stroke();

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
    const px = c + Math.cos(angle) * r, py = c + Math.sin(angle) * r;
    const pr = 7;
    const sunDir = Math.atan2(c - py, c - px);
    const g = ctx.createLinearGradient(
      px + Math.cos(sunDir) * pr, py + Math.sin(sunDir) * pr,
      px - Math.cos(sunDir) * pr, py - Math.sin(sunDir) * pr);
    g.addColorStop(0, '#7fb2d8');
    g.addColorStop(0.55, '#3b6b8c');
    g.addColorStop(1, '#0d1a26');
    ctx.fillStyle = g;
    ctx.beginPath(); ctx.arc(px, py, pr, 0, Math.PI * 2); ctx.fill();
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

    if (visible) requestAnimationFrame(frame); else last = 0;
  }
  const clampN = (x, lo, hi) => Math.max(lo, Math.min(hi, x));

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
    els.dist.value = sliderFromDist(d);
    angle = Math.atan2(y, x);
    summary();
  }
  orrery.addEventListener('pointerdown', (e) => {
    if (e.button !== 0) return;
    dragging = true;
    orrery.setPointerCapture?.(e.pointerId);
    setFromPointer(e);
  });
  orrery.addEventListener('pointermove', (e) => { if (dragging) setFromPointer(e); });
  const stop = () => { dragging = false; };
  orrery.addEventListener('pointerup', stop);
  orrery.addEventListener('pointercancel', stop);

  // to the middle of the habitable zone
  els.hz.addEventListener('click', () => {
    const s = state();
    els.dist.value = sliderFromDist(clampN(Math.sqrt(s.L), D_MIN, D_MAX));
    summary();
  });

  box.addEventListener('input', summary);
  box.addEventListener('change', summary);
  addEventListener('resize', () => { resize(); });

  els.dist.max = SLIDER;
  els.dist.value = sliderFromDist(1);
  resize();
  summary();
  requestAnimationFrame(frame);

  window.PlanetSystem = { planet };
})();
