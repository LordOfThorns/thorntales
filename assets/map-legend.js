/* map-legend.js — legend under the terrain map (terrain.html).
   Lists every zone the generator can draw. Needs terrain-map.js and i18n.js.
   The generator calls MapLegend.render(mode) after every render. */
(function () {
  'use strict';

  // Zones in the order of climate belts, so neighbours on the map are neighbours here:
  // high mountains and cold -> temperate -> dry -> tropical
  const ZONE_ORDER = [
    6,  // snow and glaciers
    14, // bare rock
    13, // alpine meadows
    5,  // tundra
    4,  // taiga
    11, // temperate rainforest
    9,  // deciduous forest
    15, // wetlands
    7,  // steppe, grassland
    12, // shrubland
    2,  // cold desert, semi-desert
    1,  // hot desert
    3,  // savanna
    8,  // seasonal tropical forest
    10, // tropical rainforest
  ];

  const t = (key, fallback) => (window.I18n && I18n.has(key) ? I18n.t(key) : (fallback ?? key));
  const rgb = (c) => `rgb(${c[0] | 0}, ${c[1] | 0}, ${c[2] | 0})`;
  const mix = (a, b, k) => [a[0] + (b[0] - a[0]) * k, a[1] + (b[1] - a[1]) * k, a[2] + (b[2] - a[2]) * k];

  function item(swatchStyle, label, extraClass) {
    // Size inline too, so colours show even with a stale cached stylesheet
    const size = extraClass ? 'width:22px;height:4px' : 'width:22px;height:14px';
    return `<li class="legend__item">
      <span class="legend__swatch${extraClass ? ' ' + extraClass : ''}" style="display:inline-block;vertical-align:middle;${size};${swatchStyle}"></span>
      <span class="legend__label">${escapeHtml(label)}</span>
    </li>`;
  }

  // Stepped colour bar with labels at both ends
  function ramp(colors, left, right, hint) {
    const gradient = colors
      .map((c, k) => `${rgb(c)} ${(k / colors.length) * 100}% ${((k + 1) / colors.length) * 100}%`)
      .join(', ');
    return `
      <div class="legend__relief">
        <span class="legend__ramp" style="display:block;height:14px;background:linear-gradient(90deg, ${gradient})"></span>
        <span class="legend__ramp-labels"><span>${escapeHtml(left)}</span><span>${escapeHtml(right)}</span></span>
        ${hint ? `<span class="legend__hint">${escapeHtml(hint)}</span>` : ''}
      </div>`;
  }

  // Lowland colour with n ink layers on top, as the map draws elevation
  function reliefSteps(atlas, base) {
    const steps = [];
    let c = base;
    steps.push(c);
    for (let k = 0; k < atlas.layers.length; k++) {
      c = mix(c, atlas.ink, atlas.inkAlpha);
      steps.push(c);
    }
    return steps;
  }

  function render(mode) {
    const host = document.getElementById('legend');
    const TM = window.TerrainMap;
    if (!host || !TM) return;
    mode = mode || TM.lastMode || 'biomes';

    if (mode === 'wind') {
      host.hidden = true;
      return;
    }
    host.hidden = false;
    const A = TM.atlas;

    const water = [
      item(`background:${rgb(A.lake)}`, t('legend.lake')),
      item(`background:${rgb(A.river)}`, t('legend.river'), 'legend__swatch--line'),
    ].join('');
    const bands = A.depthBands.map((b) => b.c);
    const depthRamp = ramp(bands, t('legend.shelf'), t('legend.ocean'));

    const steps = reliefSteps(A, mode === 'height' ? A.lowland : [196, 200, 170]);
    const relief = ramp(steps, t('legend.low'), t('legend.high'), t('legend.reliefHint'));

    let zones = '';
    if (mode === 'biomes') {
      zones = ZONE_ORDER
        .filter((id) => TM.biomeColors[id])
        .map((id) => item(`background:${rgb(TM.biomeColors[id])}`, t('biome.' + id)))
        .join('');
    }

    host.innerHTML = `
      <h2 class="legend__title">${escapeHtml(t('legend.title'))}</h2>
      ${zones ? `<ul class="legend__list">${zones}</ul>` : ''}
      <h3 class="legend__subtitle">${escapeHtml(t('legend.water'))}</h3>
      ${depthRamp}
      <ul class="legend__list legend__list--after-ramp">${water}</ul>
      <h3 class="legend__subtitle">${escapeHtml(t('legend.relief'))}</h3>
      ${relief}`;
  }

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }

  window.MapLegend = { render };
  // The first map is generated before this file loads
  render();
})();
