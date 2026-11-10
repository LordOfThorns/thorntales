/* map-hover.js — "Zone under cursor" on the terrain map (terrain.html).
   When the toggle is on, a small label with the zone's colour and name follows the
   cursor (or the finger on touch screens). Needs terrain-map.js and i18n.js. */
(function () {
  'use strict';

  const STORE_KEY = 'thorntales.mapHover';
  const canvas = document.getElementById('map');
  const tip = document.getElementById('mapTip');
  const toggle = document.getElementById('hoverZones');
  if (!canvas || !tip || !toggle) return;

  const t = (key, fallback) => (window.I18n && I18n.has(key) ? I18n.t(key) : (fallback ?? key));
  const rgb = (c) => `rgb(${c[0] | 0}, ${c[1] | 0}, ${c[2] | 0})`;

  // Remember the choice (storage can be blocked: then it simply starts off)
  try { toggle.checked = localStorage.getItem(STORE_KEY) === '1'; } catch (e) { /* ignore */ }
  toggle.addEventListener('change', () => {
    try { localStorage.setItem(STORE_KEY, toggle.checked ? '1' : '0'); } catch (e) { /* ignore */ }
    if (!toggle.checked) hide();
  });

  // What is at map cell (x, y): { color, label } or null
  function zoneAt(x, y) {
    const TM = window.TerrainMap;
    const ids = TM && TM.lastBiomeIds;
    const n = TM && TM.lastSize;
    if (!ids || !n) return null;
    const i = y * n + x;
    if (TM.lastLake && TM.lastLake[i]) return { color: TM.atlas.lake, label: t('legend.lake') };
    const id = ids[i];
    if (id === 0) return { color: TM.atlas.depthBands[1].c, label: t('legend.sea') };
    const color = TM.biomeColors[id];
    return color ? { color, label: t('biome.' + id) } : null;
  }

  function show(e) {
    const mode = document.getElementById('mode');
    if (!toggle.checked || (mode && mode.value === 'wind')) return hide();

    const rect = canvas.getBoundingClientRect();
    const n = canvas.width;
    // the canvas may be stretched by CSS, so map x and y separately
    const x = Math.floor(((e.clientX - rect.left) / rect.width) * n);
    const y = Math.floor(((e.clientY - rect.top) / rect.height) * n);
    if (x < 0 || y < 0 || x >= n || y >= n) return hide();

    const z = zoneAt(x, y);
    if (!z) return hide();

    tip.innerHTML = `<span class="map-tip__swatch" style="background:${rgb(z.color)}"></span>${escapeHtml(z.label)}`;
    tip.hidden = false;

    // position next to the cursor inside .map-block, flipped near the right edge
    const host = tip.offsetParent || canvas.parentElement;
    const hostRect = host.getBoundingClientRect();
    let left = e.clientX - hostRect.left + 14;
    const top = e.clientY - hostRect.top + 14;
    if (left + tip.offsetWidth > hostRect.width) left = e.clientX - hostRect.left - tip.offsetWidth - 10;
    tip.style.left = `${left}px`;
    tip.style.top = `${top}px`;
  }

  function hide() {
    tip.hidden = true;
  }

  canvas.addEventListener('pointermove', show);
  canvas.addEventListener('pointerdown', show);
  canvas.addEventListener('pointerleave', hide);

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }
})();
