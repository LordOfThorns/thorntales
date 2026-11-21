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
      canvas,
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
        keepMap();
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

  // ---------- Dragging the map east-west ----------
  // The world wraps around (terrain-map.js), so the map can be turned like a globe: the
  // drawn map is kept in a copy and shown from viewOffset, the part that leaves one edge
  // comes back at the other. map-hover.js reads TerrainMap.viewOffset.
  let base = null;    // the whole map as drawn by TerrainMap.render
  let offset = 0;     // map pixel shown at the left edge of the canvas

  function keepMap() {
    base = base || document.createElement('canvas');
    base.width = canvas.width;
    base.height = canvas.height;
    base.getContext('2d').drawImage(canvas, 0, 0);
    offset %= canvas.width; // the same world turned the same way after a new drawing
    show();
  }

  function show() {
    if (!base) return;
    const n = base.width;
    offset = ((offset % n) + n) % n;
    const ctx = canvas.getContext('2d');
    const o = Math.round(offset);
    ctx.drawImage(base, -o, 0);
    if (o) ctx.drawImage(base, n - o, 0);
    TerrainMap.viewOffset = o;
  }

  let drag = null;
  canvas.addEventListener('pointerdown', (e) => {
    if (busy || !base || e.button !== 0) return;
    drag = { x: e.clientX, offset };
    canvas.setPointerCapture?.(e.pointerId);
    canvas.classList.add('is-dragging');
  });
  canvas.addEventListener('pointermove', (e) => {
    if (!drag) return;
    const scale = canvas.width / canvas.getBoundingClientRect().width; // css px -> map px
    offset = drag.offset - (e.clientX - drag.x) * scale;
    show();
  });
  const stop = () => { drag = null; canvas.classList.remove('is-dragging'); };
  canvas.addEventListener('pointerup', stop);
  canvas.addEventListener('pointercancel', stop);
  // keyboard: the map can be focused and turned with the arrow keys
  canvas.tabIndex = 0;
  canvas.addEventListener('keydown', (e) => {
    if (!base || (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight')) return;
    e.preventDefault();
    offset += (e.key === 'ArrowRight' ? 1 : -1) * canvas.width / 24;
    show();
  });

  showValues();
  generate();
})();
