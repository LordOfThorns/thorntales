/* terrain-page.js — the settings panel of the map generator (terrain.html).

   The panel speaks the reader's language (i18n.js) and shows the settings as plain choices;
   here they become the numbers TerrainMap.render() expects:
     sea       5..95 %             -> seaLevel 0.05..0.95
     rivers    0 (few)..100 (many) -> riverThreshold 2000..200 (a lower threshold, more rivers)
     mountains 0..100 %            -> mountainsStrength 0..0.6
     deserts   0..100 %            -> desertStrength 0..1

   Drawing a map takes one to four seconds and blocks the page, so a spinner is shown
   first (two frames are let through so it is painted; it turns with a CSS transform, which
   keeps going while the script works), the button is off until the map is ready, and clicks
   made meanwhile are dropped.
*/
(function () {
  'use strict';

  const $ = (id) => document.getElementById(id);
  const form = $('mapForm');
  const canvas = $('map');
  const button = $('gen');
  const loading = $('mapLoading');
  if (!form || !canvas || !window.TerrainMap) return;

  const int = (el, lo, hi, def) => {
    const v = Math.round(Number(el.value));
    return Number.isFinite(v) ? Math.max(lo, Math.min(hi, v)) : def;
  };

  // values shown next to the sliders
  const outputs = [
    ['sea', 'seaOut'], ['rivers', 'riversOut'], ['mount', 'mountOut'], ['desert', 'desertOut'],
  ];
  function showValues() {
    for (const [inp, out] of outputs) $(out).textContent = `${$(inp).value}%`;
  }

  function settings() {
    let contMin = int($('contMin'), 1, 10, 2);
    let contMax = int($('contMax'), 1, 10, 6);
    if (contMax < contMin) [contMin, contMax] = [contMax, contMin];
    $('contMin').value = contMin;
    $('contMax').value = contMax;
    const seed = int($('seed'), -2147483648, 2147483647, 6000);
    $('seed').value = seed;
    return {
      canvas: work, // drawn off screen, then shown through the window (show())
      seed,
      size: int($('size'), 512, 2048, 1280), // width; the height is half
      mode: $('mode').value || 'biomes',
      continentsMin: contMin,
      continentsMax: contMax,
      seaLevel: int($('sea'), 5, 95, 30) / 100,
      riverThreshold: Math.round(2000 - int($('rivers'), 0, 100, 78) * 18),
      mountainsStrength: (int($('mount'), 0, 100, 30) / 100) * 0.6,
      desertStrength: int($('desert'), 0, 100, 35) / 100,
      plates: 18,
      windDir: 'W',
      planet: window.PlanetSystem ? PlanetSystem.planet() : null, // prototype (planet-system.js)
    };
  }

  let busy = false;
  function setBusy(on) {
    busy = on;
    button.disabled = on;
    $('seedRandom').disabled = on;
    loading.hidden = !on;
    form.setAttribute('aria-busy', String(on));
    canvas.classList.toggle('is-busy', on);
  }

  function generate() {
    if (busy) return; // a map is being drawn: impatient clicks change nothing
    const opts = settings();
    setBusy(true);
    // let the browser paint the spinner before the heavy work blocks the page
    requestAnimationFrame(() => requestAnimationFrame(() => setTimeout(() => {
      try {
        TerrainMap.render(opts);
        if (work.width !== shownN) scale = 0; // another size: choose the scale again
        shownN = work.width;
        compose();
        draw();
        if (window.MapLegend) MapLegend.render(opts.mode);
      } catch (e) {
        console.error('map failed', e);
      }
      // clicks made while the page was busy arrive right after: they meet a disabled button
      setTimeout(() => setBusy(false), 150);
    }, 0)));
  }

  form.addEventListener('submit', (e) => {
    e.preventDefault();
    generate();
  });
  form.addEventListener('input', (e) => {
    if (e.target.classList.contains('map-slider')) showValues();
  });
  $('seedRandom').addEventListener('click', () => {
    if (busy) return;
    $('seed').value = 1 + Math.floor(Math.random() * 999999);
    generate();
  });

  // ---------- The window onto the world ----------
  // The world is twice as wide as high and wraps around (terrain-map.js). It is drawn off
  // screen (work) and shown through a window (.map-stage): as wide as the page and about as
  // tall as the old square map. The reader can drag the map inside it (east-west it never
  // ends, north-south it stops at the poles), resize the window by its bottom-right corner
  // (CSS resize) and switch to the zoomed-out view, where the whole world fits the old width.
  // TerrainMap.view = { ox, oy, scale } (map pixel at the top-left corner, css px per map
  // pixel) is read by map-hover.js.
  const work = document.createElement('canvas');
  const stage = canvas.parentElement;
  const zoomBtn = $('zoomToggle');
  const COLUMN = 948;   // width of the page column, the old width of the map
  let base = null;      // what is shown: the drawn map (work) with the time of year and the starlight
  const shown = document.createElement('canvas');
  let ox = 0, oy = 0;   // map pixel shown at the top-left corner of the window
  let scale = 0;        // css px per map pixel; 0 = choose for the current layout
  let zoomedOut = false;
  let shownN = 0;       // width of the map now shown

  // Window size for the current mode (forgets a size the reader set by hand)
  function layout() {
    const full = stage.parentElement.clientWidth;
    if (zoomedOut) {
      const w = Math.min(full, COLUMN);
      stage.style.width = `${w}px`;
      stage.style.height = `${Math.round(w / 2)}px`;
    } else {
      // tall enough to look at, within the screen, no taller than wide (phones get a square)
      // and never so low that the world would be seen twice
      const h = Math.max(Math.min(800, full, Math.max(260, innerHeight * 0.75)), full / 2);
      stage.style.width = '';
      stage.style.height = `${Math.round(h)}px`;
    }
    scale = 0;
    if (zoomBtn) {
      zoomBtn.dataset.i18n = zoomedOut ? 'terrain.zoomIn' : 'terrain.zoomOut';
      if (window.I18n) zoomBtn.textContent = I18n.t(zoomBtn.dataset.i18n);
    }
    draw();
  }

  // The drawn map plus what needs no new map (planet-system.js): snow, ice and dry seasons
  // for the chosen moment of the year (or ice that stays all year), and the colour of
  // the starlight. Not in the wind view.
  function compose() {
    if (!work.width) return;
    shown.width = work.width;
    shown.height = work.height;
    const g = shown.getContext('2d');
    g.drawImage(work, 0, 0);
    const PS = window.PlanetSystem;
    if (TerrainMap.lastMode !== 'wind' && TerrainMap.seasonLayer) {
      const layer = TerrainMap.seasonLayer(PS ? PS.season() : null);
      if (layer) g.drawImage(layer, 0, 0);
      const tint = PS && PS.tint();
      if (tint) {
        g.globalCompositeOperation = 'multiply';
        g.globalAlpha = tint.k;
        g.fillStyle = `rgb(${tint.color.join(',')})`;
        g.fillRect(0, 0, shown.width, shown.height);
        g.globalCompositeOperation = 'source-over';
        g.globalAlpha = 1;
      }
    } else {
      TerrainMap.lastIce = null;
    }
    base = shown;
    texDirty = true;
    draw();
  }
  if (window.PlanetSystem) PlanetSystem.onLive(() => { if (!busy && base) compose(); });

  function draw() {
    const w = stage.clientWidth, h = stage.clientHeight;
    const dpr = window.devicePixelRatio || 1;
    const cw = Math.max(1, Math.round(w * dpr)), ch = Math.max(1, Math.round(h * dpr));
    if (canvas.width !== cw || canvas.height !== ch) { canvas.width = cw; canvas.height = ch; }
    if (!base || !w || !h) return;
    if (globe) { drawGlobe(); return; }
    TerrainMap.pick = null;
    const n = base.width, nh = base.height;
    if (!scale) scale = zoomedOut ? w / n : h / nh;
    scale = Math.max(scale, w / n); // never show the world twice side by side
    // a window taller than the map shrinks to it
    if (h > nh * scale + 1) { stage.style.height = `${Math.floor(nh * scale)}px`; return; }
    const viewW = w / scale, viewH = h / scale;
    ox = ((ox % n) + n) % n;
    oy = Math.max(0, Math.min(nh - viewH, oy));
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr * scale, 0, 0, dpr * scale, 0, 0);
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(base, -ox, -oy);
    if (n - ox < viewW) ctx.drawImage(base, n - ox, -oy); // the part past the seam
    TerrainMap.view = { ox, oy, scale };
  }

  // ---------- The globe ----------
  // "Show as a globe": the same map (with the time of year and the starlight) wrapped on a
  // sphere, orthographic view. Every pixel of the disc is turned back into latitude and
  // longitude and takes its colour from the map; a soft light from the upper left and a
  // blue haze at the rim (the air). Drag to turn it; left alone it slowly turns by itself.
  // The disc is drawn at most GLOBE_RES px across and scaled up, to stay fast.
  // TerrainMap.pick(x, y) (css px in the canvas -> map pixel) is read by map-hover.js.
  const globeBtn = $('globeToggle');
  const GLOBE_RES = 560;
  let globe = false;
  let gLon = 0, gTilt = 0.35;   // longitude at the centre of the disc, tilt towards the viewer
  let lastTouch = 0;            // the globe turns by itself again 2.5 s after the reader lets go
  let tex = null, texDirty = true;
  let geo = null;               // per-pixel point on the sphere, light and rim for size geo.S
  const gCanvas = document.createElement('canvas');
  let gImg = null;

  function makeGeo(S) {
    const N = S * S;
    const x = new Float32Array(N), y = new Float32Array(N), z = new Float32Array(N);
    const light = new Float32Array(N), rim = new Float32Array(N), alpha = new Float32Array(N);
    const L = [-0.48, 0.42, 0.77]; // light from the upper left, a bit in front
    for (let j = 0; j < S; j++) {
      for (let i = 0; i < S; i++) {
        const p = j * S + i;
        const nx = ((i + 0.5) / S) * 2 - 1, ny = ((j + 0.5) / S) * 2 - 1;
        const r2 = nx * nx + ny * ny;
        if (r2 > 1) continue;
        const nz = Math.sqrt(1 - r2);
        x[p] = nx; y[p] = -ny; z[p] = nz;
        const d = Math.max(0, nx * L[0] - ny * L[1] + nz * L[2]);
        light[p] = Math.min(1.08, 0.38 + 0.74 * d);
        rim[p] = Math.pow(1 - nz, 2.2) * 0.6;
        alpha[p] = Math.min(1, (1 - Math.sqrt(r2)) * S / 2); // soft edge, one pixel wide
      }
    }
    gCanvas.width = gCanvas.height = S;
    gImg = gCanvas.getContext('2d').createImageData(S, S);
    return { S, x, y, z, light, rim, alpha };
  }

  // For the current tilt: map row start (row * width) and longitude 0..1 of each pixel
  function tiltTables(th) {
    const { S, x, y, z, alpha } = geo;
    const N = S * S;
    geo.rows = geo.rows || new Int32Array(N);
    geo.lon = geo.lon || new Float32Array(N);
    const ct = Math.cos(gTilt), st = Math.sin(gTilt);
    const tw = tex.width;
    for (let p = 0; p < N; p++) {
      if (alpha[p] <= 0) continue;
      // from the view back to the globe: undo the tilt
      const yw = y[p] * ct + z[p] * st;
      const zw = -y[p] * st + z[p] * ct;
      const lat = Math.asin(yw > 1 ? 1 : yw < -1 ? -1 : yw);
      const ty = Math.min(th - 1, Math.max(0, ((0.5 - lat / Math.PI) * th) | 0));
      geo.rows[p] = ty * tw;
      geo.lon[p] = Math.atan2(x[p], zw) / (Math.PI * 2) + 0.5; // 0..1
    }
    geo.tilt = gTilt;
    geo.th = th;
  }

  // disc size and place in the canvas (device px)
  function globeBox() {
    const W = canvas.width, H = canvas.height;
    const D = Math.min(W, H) * 0.9;
    return { W, H, D, left: (W - D) / 2, top: (H - D) / 2 };
  }

  function drawGlobe() {
    const ctx = canvas.getContext('2d');
    const { W, H, D, left, top } = globeBox();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, W, H);
    if (texDirty || !tex) {
      tex = base.getContext('2d').getImageData(0, 0, base.width, base.height);
      texDirty = false;
    }
    const S = Math.max(64, Math.min(GLOBE_RES, Math.round(D)));
    if (!geo || geo.S !== S) geo = makeGeo(S);
    const { light, rim, alpha } = geo;
    const tw = tex.width, th = tex.height, td = tex.data, od = gImg.data;
    const TWO_PI = Math.PI * 2;
    // latitude (map row) and longitude of every pixel depend only on the tilt: worked out
    // once per tilt; turning the globe then only shifts the longitude
    if (!geo.rows || geo.tilt !== gTilt || geo.th !== th) tiltTables(th);
    const { rows, lon } = geo;
    let shift = gLon / TWO_PI; // lon[] already holds the +0.5
    shift -= Math.floor(shift);
    for (let p = 0, N = S * S; p < N; p++) {
      const o = p * 4;
      const a = alpha[p];
      if (a <= 0) { od[o + 3] = 0; continue; }
      let u = lon[p] + shift;
      if (u >= 1) u -= 1;
      const tx = (u * tw) | 0;
      const ti = (rows[p] + (tx < tw ? tx : tw - 1)) * 4;
      const l = light[p], r = rim[p];
      od[o] = td[ti] * l * (1 - r) + 150 * r;
      od[o + 1] = td[ti + 1] * l * (1 - r) + 192 * r;
      od[o + 2] = td[ti + 2] * l * (1 - r) + 255 * r;
      od[o + 3] = a * 255;
    }
    gCanvas.getContext('2d').putImageData(gImg, 0, 0);
    // the air: a faint glow around the disc
    const cx = W / 2, cy = H / 2, R = D / 2;
    const halo = ctx.createRadialGradient(cx, cy, R * 0.96, cx, cy, R * 1.09);
    halo.addColorStop(0, 'rgba(130, 180, 255, 0.35)');
    halo.addColorStop(1, 'rgba(130, 180, 255, 0)');
    ctx.fillStyle = halo;
    ctx.beginPath(); ctx.arc(cx, cy, R * 1.09, 0, TWO_PI); ctx.fill();
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(gCanvas, left, top, D, D);
    TerrainMap.view = null;
    TerrainMap.pick = pickGlobe;
  }

  // css px in the canvas -> map pixel { x, y } on the globe, or null off the disc
  function pickGlobe(cx, cy) {
    if (!globe || !base) return null;
    const dpr = canvas.width / (canvas.getBoundingClientRect().width || 1);
    const { W, H, D } = globeBox();
    const nx = (cx * dpr - W / 2) / (D / 2), ny = (cy * dpr - H / 2) / (D / 2);
    const r2 = nx * nx + ny * ny;
    if (r2 > 1) return null;
    const vy = -ny, vz = Math.sqrt(1 - r2);
    const yw = vy * Math.cos(gTilt) + vz * Math.sin(gTilt);
    const zw = -vy * Math.sin(gTilt) + vz * Math.cos(gTilt);
    let u = (Math.atan2(nx, zw) + gLon) / (Math.PI * 2) + 0.5;
    u -= Math.floor(u);
    const n = base.width, nh = base.height;
    return {
      x: Math.min(n - 1, Math.floor(u * n)),
      y: Math.min(nh - 1, Math.max(0, Math.floor((0.5 - Math.asin(clampN(yw, -1, 1)) / Math.PI) * nh))),
    };
  }
  const clampN = (v, lo, hi) => Math.max(lo, Math.min(hi, v));

  // slow turning while nobody touches it, only while it is on screen
  let onScreen = true, spinRaf = 0, spinLast = 0;
  function spin(now) {
    spinRaf = 0;
    if (!globe) return;
    const dt = spinLast ? Math.min(0.1, (now - spinLast) / 1000) : 0;
    spinLast = now;
    if (!drag && !busy && now - lastTouch > 2500 && onScreen && base) {
      gLon += dt * 0.12;
      draw();
    }
    spinRaf = requestAnimationFrame(spin);
  }
  if (window.IntersectionObserver) {
    new IntersectionObserver((en) => { onScreen = en[0].isIntersecting; }).observe(stage);
  }

  function setGlobe(on) {
    globe = on;
    stage.classList.toggle('is-globe', on);
    if (globeBtn) {
      globeBtn.dataset.i18n = on ? 'terrain.flat' : 'terrain.globe';
      if (window.I18n) globeBtn.textContent = I18n.t(globeBtn.dataset.i18n);
    }
    const hint = $('mapHint');
    if (hint) {
      hint.dataset.i18n = on ? 'terrain.globeHint' : 'terrain.dragHint';
      if (window.I18n) hint.textContent = I18n.t(hint.dataset.i18n);
    }
    if (on) {
      // start with the part of the map that was in the middle of the window
      if (base && scale) gLon = ((ox + stage.clientWidth / scale / 2) / base.width - 0.5) * Math.PI * 2;
      spinLast = 0;
      if (!spinRaf) spinRaf = requestAnimationFrame(spin);
    } else {
      // back to the flat map, the same longitude in the middle
      if (base && scale) ox = (gLon / (Math.PI * 2) + 0.5) * base.width - stage.clientWidth / scale / 2;
      TerrainMap.pick = null;
    }
    draw();
  }
  if (globeBtn) globeBtn.addEventListener('click', () => setGlobe(!globe));

  // the window follows the page width, the reader's resizing and its own height changes
  if (window.ResizeObserver) new ResizeObserver(() => draw()).observe(stage);
  let resizeTimer = 0;
  addEventListener('resize', () => {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(() => { if (zoomedOut) layout(); else draw(); }, 100);
  });

  if (zoomBtn) zoomBtn.addEventListener('click', () => { zoomedOut = !zoomedOut; layout(); });

  // the corner grip: the reader makes the window smaller or larger (the map keeps its scale)
  const grip = document.createElement('span');
  grip.className = 'map-resize';
  grip.setAttribute('aria-hidden', 'true');
  stage.appendChild(grip);
  let sizing = null;
  grip.addEventListener('pointerdown', (e) => {
    if (e.button !== 0) return;
    e.preventDefault();
    sizing = { x: e.clientX, y: e.clientY, w: stage.offsetWidth, h: stage.offsetHeight };
    grip.setPointerCapture?.(e.pointerId);
  });
  grip.addEventListener('pointermove', (e) => {
    if (!sizing) return;
    const full = stage.parentElement.clientWidth;
    // the window stays centred, so the width changes twice as fast as the pointer moves
    const w = Math.max(200, Math.min(full, sizing.w + (e.clientX - sizing.x) * 2));
    const h = Math.max(140, sizing.h + (e.clientY - sizing.y));
    stage.style.width = w >= full ? '' : `${Math.round(w)}px`;
    stage.style.height = `${Math.round(h)}px`;
    draw();
  });
  const endSizing = () => { sizing = null; };
  grip.addEventListener('pointerup', endSizing);
  grip.addEventListener('pointercancel', endSizing);

  let drag = null;
  canvas.addEventListener('pointerdown', (e) => {
    if (busy || !base || e.button !== 0) return;
    drag = { x: e.clientX, y: e.clientY, ox, oy, lon: gLon, tilt: gTilt };
    try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* not a real pointer */ }
    canvas.classList.add('is-dragging');
  });
  canvas.addEventListener('pointermove', (e) => {
    if (!drag) return;
    if (globe) {
      // the point under the pointer follows it: one disc radius of travel = one radian
      const R = globeBox().D / 2 / (canvas.width / (canvas.getBoundingClientRect().width || 1));
      gLon = drag.lon - (e.clientX - drag.x) / R;
      gTilt = clampN(drag.tilt + (e.clientY - drag.y) / R, -1.4, 1.4);
      lastTouch = performance.now();
    } else {
      ox = drag.ox - (e.clientX - drag.x) / scale;
      oy = drag.oy - (e.clientY - drag.y) / scale;
    }
    draw();
  });
  const stop = () => { drag = null; canvas.classList.remove('is-dragging'); };
  canvas.addEventListener('pointerup', stop);
  canvas.addEventListener('pointercancel', stop);
  // keyboard: the map can be focused and moved with the arrow keys
  canvas.tabIndex = 0;
  const KEYS = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1] };
  canvas.addEventListener('keydown', (e) => {
    const k = KEYS[e.key];
    if (!base || !k) return;
    e.preventDefault();
    if (globe) {
      gLon += k[0] * Math.PI / 12;
      gTilt = clampN(gTilt - k[1] * 0.15, -1.4, 1.4);
      lastTouch = performance.now();
    } else {
      ox += k[0] * base.width / 24;
      oy += k[1] * base.height / 12;
    }
    draw();
  });

  showValues();
  layout();
  generate();
})();
