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
        base = work;
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
  let base = null;      // the whole map as drawn by TerrainMap.render
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

  function draw() {
    const w = stage.clientWidth, h = stage.clientHeight;
    const dpr = window.devicePixelRatio || 1;
    const cw = Math.max(1, Math.round(w * dpr)), ch = Math.max(1, Math.round(h * dpr));
    if (canvas.width !== cw || canvas.height !== ch) { canvas.width = cw; canvas.height = ch; }
    if (!base || !w || !h) return;
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
    drag = { x: e.clientX, y: e.clientY, ox, oy };
    canvas.setPointerCapture?.(e.pointerId);
    canvas.classList.add('is-dragging');
  });
  canvas.addEventListener('pointermove', (e) => {
    if (!drag) return;
    ox = drag.ox - (e.clientX - drag.x) / scale;
    oy = drag.oy - (e.clientY - drag.y) / scale;
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
    ox += k[0] * base.width / 24;
    oy += k[1] * base.height / 12;
    draw();
  });

  showValues();
  layout();
  generate();
})();
