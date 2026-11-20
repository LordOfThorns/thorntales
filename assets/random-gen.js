/* random-gen.js
   Custom random generator of everything (random.html).
   Idea from the old Android app "Magic Of Random" (LordOfThorns/chameleonparty):
   a thing has characteristics, each characteristic has values with integer weights
   (weights 3 and 2 = 60/40), a characteristic may give several values ("from N to M").
   New here: links between values of different characteristics drawn on a map.
   A link "if A then B ×f" multiplies the weight of B by f once A is rolled
   (0 = never, 0.5 = half as often, 2 = twice as often); "always" forces B.
   Weights are applied first, factors second (as planned in the old blog post).
   Characteristics are rolled in dependency order; links that close a cycle are ignored.
   Storage: the current generator only (localStorage); keep others as JSON files; share link (#g=...).
*/
(function () {
  'use strict';

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const LANG = (window.I18n && I18n.lang) || 'ru';

  const STORE_KEY = 'thorn-rng-v1';
  const NODE_W = 232;
  const GAP_X = 72;
  const GAP_Y = 36;
  const LEFT = 64; // room on the left for loops of links inside a characteristic
  const MAX_COUNT = 50;
  const MAX_HISTORY = 10;

  const root = document.getElementById('rng');
  if (!root) return;

  // ---------- small helpers ----------

  function h(tag, attrs, ...kids) {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v == null || v === false) continue;
      if (k === 'class') el.className = v;
      else if (k === 'text') el.textContent = v;
      else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
      else if (k in el && k !== 'list') el[k] = v;
      else el.setAttribute(k, v === true ? '' : v);
    }
    for (const kid of kids.flat(Infinity)) {
      if (kid == null || kid === false) continue;
      el.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
    }
    return el;
  }

  const uid = () => Math.random().toString(36).slice(2, 9);
  const clampInt = (v, lo, hi) => Math.max(lo, Math.min(hi, Math.round(Number(v) || 0)));

  // ---------- model ----------

  function newValue(text = '') { return { id: uid(), text, w: 1 }; }
  function newChar(name = '') {
    return { id: uid(), name, x: 0, y: 0, multi: false, min: 1, max: 2, values: [newValue(), newValue()] };
  }
  function newGen(name = '') {
    return { id: uid(), name, count: 1, chars: [newChar()], links: [], needsLayout: true };
  }

  // Make any parsed object a valid generator (files and links come from outside)
  function normalize(raw) {
    if (!raw || typeof raw !== 'object') throw new Error('bad');
    const str = (v, n = 200) => String(v == null ? '' : v).slice(0, n);
    const g = {
      id: uid(),
      name: str(raw.name, 120),
      count: clampInt(raw.count || 1, 1, MAX_COUNT),
      chars: [],
      links: [],
    };
    const valIds = new Map();
    for (const c of Array.isArray(raw.chars) ? raw.chars.slice(0, 60) : []) {
      if (!c || typeof c !== 'object') continue;
      const ch = {
        id: str(c.id, 20) || uid(),
        name: str(c.name, 120),
        x: Math.max(0, Number(c.x) || 0),
        y: Math.max(0, Number(c.y) || 0),
        multi: !!c.multi,
        min: clampInt(c.min || 1, 1, 99),
        max: clampInt(c.max || 2, 1, 99),
        values: [],
      };
      if (g.chars.some((o) => o.id === ch.id)) ch.id = uid();
      for (const v of Array.isArray(c.values) ? c.values.slice(0, 200) : []) {
        if (!v || typeof v !== 'object') continue;
        const val = { id: str(v.id, 20) || uid(), text: str(v.text, 300), w: clampInt(v.w ?? 1, 0, 1e6) };
        if (ch.values.some((o) => o.id === val.id)) val.id = uid();
        ch.values.push(val);
      }
      valIds.set(ch.id, new Set(ch.values.map((v) => v.id)));
      g.chars.push(ch);
    }
    for (const l of Array.isArray(raw.links) ? raw.links.slice(0, 500) : []) {
      if (!l || !l.a || !l.b) continue;
      const ok = (e) => valIds.get(e.c) && valIds.get(e.c).has(e.v);
      if (!ok(l.a) || !ok(l.b) || l.a.v === l.b.v) continue;
      g.links.push({
        id: uid(),
        a: { c: l.a.c, v: l.a.v },
        b: { c: l.b.c, v: l.b.v },
        f: Math.max(0, Math.min(1000, Number(l.f) >= 0 ? Number(l.f) : 2)),
        always: !!l.always,
        ox: Math.max(-5000, Math.min(5000, Math.round(Number(l.ox) || 0))),
        oy: Math.max(-5000, Math.min(5000, Math.round(Number(l.oy) || 0))),
      });
    }
    return g;
  }

  // Only what matters for the generator: used for files and share links
  function exportable(g) {
    return {
      name: g.name,
      count: g.count,
      chars: g.chars.map((c) => ({
        id: c.id, name: c.name, x: Math.round(c.x), y: Math.round(c.y),
        multi: c.multi, min: c.min, max: c.max,
        values: c.values.map((v) => ({ id: v.id, text: v.text, w: v.w })),
      })),
      links: g.links.map((l) => Object.assign({ a: l.a, b: l.b, f: l.f, always: l.always }, l.ox || l.oy ? { ox: l.ox, oy: l.oy } : {})),
    };
  }

  // ---------- examples (ru | en | de) ----------

  // chars: [name, multi ? [min, max] : null, [[text, weight], ...]]
  // links: [charIndex, valueIndex, charIndex, valueIndex, factor or 'always']
  const EXAMPLES = [
    {
      name: 'Кот|Cat|Katze',
      count: 3,
      chars: [
        ['Окрас|Coat|Fell', null, [['чёрный|black|schwarz', 2], ['белый|white|weiß', 1], ['рыжий|ginger|rot', 2], ['серый|grey|grau', 2], ['трёхцветный|calico|dreifarbig', 1]]],
        ['Глаза|Eyes|Augen', null, [['зелёные|green|grün', 2], ['жёлтые|yellow|gelb', 2], ['голубые|blue|blau', 1], ['медные|copper|kupferfarben', 1]]],
        ['Характер|Temper|Wesen', [1, 2], [['ласковый|affectionate|verschmust', 1], ['вредный|spiteful|boshaft', 1], ['ленивый|lazy|faul', 1], ['любопытный|curious|neugierig', 1], ['пугливый|timid|scheu', 1], ['наглый|cheeky|frech', 1], ['смелый|brave|mutig', 1]]],
        ['Где спит|Sleeps|Schläft', null, [['на батарее|on the radiator|auf der Heizung', 2], ['в шкафу|in the wardrobe|im Schrank', 1], ['на клавиатуре|on the keyboard|auf der Tastatur', 2], ['на голове хозяина|on the owner’s head|auf dem Kopf des Menschen', 1]]],
      ],
      links: [[0, 3, 1, 0, 3], [0, 1, 1, 2, 3], [0, 2, 1, 3, 2], [2, 2, 3, 0, 3], [2, 4, 3, 1, 4], [2, 4, 3, 3, 0], [2, 5, 3, 3, 4], [2, 4, 2, 6, 0], [2, 4, 2, 5, 0.3], [2, 2, 2, 3, 0.5], [2, 6, 3, 1, 0.5]],
    },
    {
      name: 'Персонаж|Character|Figur',
      count: 3,
      chars: [
        ['Кто|Who|Wer', null, [['картограф|cartographer|Kartograf', 1], ['смотритель маяка|lighthouse keeper|Leuchtturmwärter', 1], ['переводчик|translator|Übersetzer', 1], ['беглый солдат|deserter|Deserteur', 1], ['аптекарь|apothecary|Apotheker', 1]]],
        ['Возраст|Age|Alter', null, [['юный|young|jung', 1], ['средних лет|middle-aged|mittleren Alters', 2], ['старый|old|alt', 1]]],
        ['Тайна|Secret|Geheimnis', null, [['должен денег|owes money|hat Schulden', 1], ['видит чужие сны|sees other people’s dreams|sieht fremde Träume', 1], ['не помнит своего имени|does not remember his name|erinnert sich nicht an seinen Namen', 1], ['служил не той стороне|served the wrong side|diente der falschen Seite', 1], ['никакой тайны|no secret at all|gar kein Geheimnis', 1]]],
        ['При себе|Carries|Trägt bei sich', [1, 2], [['ржавый ключ|a rusty key|einen rostigen Schlüssel', 1], ['карту без названий|a map without names|eine Karte ohne Namen', 1], ['неотправленное письмо|an unsent letter|einen nie abgeschickten Brief', 1], ['фонарь|a lantern|eine Laterne', 1], ['чужое кольцо|someone else’s ring|einen fremden Ring', 1]]],
      ],
      links: [[0, 0, 3, 1, 4], [0, 1, 3, 3, 4], [0, 3, 2, 3, 3], [1, 2, 2, 2, 2], [1, 0, 2, 3, 0.5]],
    },
    {
      name: 'Сон|Dream|Traum',
      count: 1,
      chars: [
        ['Где|Where|Wo', null, [['вокзал|a railway station|ein Bahnhof', 1], ['старая школа|the old school|die alte Schule', 1], ['дом без дверей|a house without doors|ein Haus ohne Türen', 1], ['берег северного моря|the shore of a northern sea|das Ufer eines nördlichen Meeres', 1], ['башня|a tower|ein Turm', 1]]],
        ['Кто там|Who is there|Wer ist da', null, [['никого|nobody|niemand', 1], ['мама|mother|Mutter', 1], ['незнакомец в пальто|a stranger in a coat|ein Fremder im Mantel', 1], ['давний друг|an old friend|ein alter Freund', 1], ['ты сам, но старше|you, but older|du selbst, aber älter', 1]]],
        ['Что не так|What is wrong|Was stimmt nicht', null, [['часы идут назад|clocks run backwards|die Uhren laufen rückwärts', 1], ['нет теней|there are no shadows|es gibt keine Schatten', 1], ['все говорят на чужом языке|everyone speaks a foreign language|alle sprechen eine fremde Sprache', 1], ['вода стоит в воздухе|water hangs in the air|Wasser steht in der Luft', 1], ['каждая дверь ведёт в ту же комнату|every door leads to the same room|jede Tür führt in dasselbe Zimmer', 1]]],
        ['Чем кончается|How it ends|Wie es endet', null, [['просыпаешься|you wake up|du wachst auf', 2], ['падаешь|you fall|du fällst', 1], ['уходит поезд|the train leaves|der Zug fährt ab', 1], ['тебя зовут, а ты не откликаешься|someone calls you and you do not answer|jemand ruft dich, und du antwortest nicht', 1]]],
      ],
      links: [[0, 0, 3, 2, 4], [0, 2, 2, 4, 0], [1, 0, 2, 2, 0], [0, 4, 3, 1, 3], [1, 4, 3, 3, 3]],
    },
  ];

  const LI = { ru: 0, en: 1, de: 2 }[LANG] ?? 0;
  const pickLang = (s) => { const p = s.split('|'); return p[LI] || p[0]; };

  function exampleGen(ex) {
    const g = { id: uid(), name: pickLang(ex.name), count: ex.count, chars: [], links: [] };
    for (const [name, multi, values] of ex.chars) {
      g.chars.push({
        id: uid(), name: pickLang(name), x: 0, y: 0,
        multi: !!multi, min: multi ? multi[0] : 1, max: multi ? multi[1] : 2,
        values: values.map(([text, w]) => ({ id: uid(), text: pickLang(text), w })),
      });
    }
    for (const [ci, vi, cj, vj, f] of ex.links) {
      g.links.push({
        id: uid(),
        a: { c: g.chars[ci].id, v: g.chars[ci].values[vi].id },
        b: { c: g.chars[cj].id, v: g.chars[cj].values[vj].id },
        f: f === 'always' ? 2 : f,
        always: f === 'always',
      });
    }
    g.needsLayout = true;
    return g;
  }

  // ---------- storage ----------

  // Only the current generator lives in the browser; whatever is worth keeping goes to a file.
  let gen = load();

  function load() {
    try {
      const raw = JSON.parse(localStorage.getItem(STORE_KEY) || 'null');
      const g = raw && (raw.gen || (Array.isArray(raw.gens) && (raw.gens.find((x) => x.id === raw.cur) || raw.gens[0])));
      if (g) return normalize(g);
    } catch (e) { /* fall through: nothing chosen yet */ }
    return null; // the page waits for "Generate" or an example
  }

  let saveTimer = 0;
  function save() {
    if (!gen) return;
    clearTimeout(saveTimer);
    saveTimer = setTimeout(() => {
      try {
        localStorage.setItem(STORE_KEY, JSON.stringify({ gen: exportable(gen) }));
      } catch (e) { /* private mode or full storage: the page still works */ }
    }, 250);
  }


  // ---------- generation ----------

  // Characteristics in an order where link sources come before targets.
  // Kahn's algorithm; what is left in a cycle keeps the order of the list.
  function rollOrder(g) {
    const ids = g.chars.map((c) => c.id);
    const indeg = new Map(ids.map((id) => [id, 0]));
    const edges = new Map(ids.map((id) => [id, new Set()]));
    for (const l of g.links) {
      if (l.a.c === l.b.c) continue; // links inside one characteristic do not order anything
      if (!edges.get(l.a.c).has(l.b.c)) { edges.get(l.a.c).add(l.b.c); indeg.set(l.b.c, indeg.get(l.b.c) + 1); }
    }
    const order = [];
    const done = new Set();
    while (order.length < ids.length) {
      let next = ids.find((id) => !done.has(id) && indeg.get(id) === 0);
      if (!next) next = ids.find((id) => !done.has(id)); // cycle: break it in list order
      done.add(next);
      order.push(next);
      for (const to of edges.get(next)) indeg.set(to, indeg.get(to) - 1);
    }
    return order;
  }

  const isInner = (l) => l.a.c === l.b.c;

  // A link between characteristics works only if its source is rolled before its target;
  // a link inside one characteristic works only if that characteristic gives several values
  function activeLinks(g) {
    const pos = new Map(rollOrder(g).map((id, i) => [id, i]));
    const multi = new Set(g.chars.filter((c) => c.multi).map((c) => c.id));
    return new Set(g.links.filter((l) => (isInner(l) ? multi.has(l.a.c) : pos.get(l.a.c) < pos.get(l.b.c))).map((l) => l.id));
  }

  function offReason(l) { return t(isInner(l) ? 'rng.innerOff' : 'rng.cycle'); }

  function weightedPick(items) {
    const total = items.reduce((s, it) => s + it.w, 0);
    if (total <= 0) return -1;
    let r = Math.random() * total;
    for (let i = 0; i < items.length; i++) {
      r -= items[i].w;
      if (r < 0) return i;
    }
    return items.length - 1;
  }

  function rollOnce(g, active) {
    const byId = new Map(g.chars.map((c) => [c.id, c]));
    const picked = new Map(); // charId -> Set of valueIds
    const out = [];
    for (const cid of rollOrder(g)) {
      const c = byId.get(cid);
      const vals = c.values.filter((v) => v.text.trim());
      if (!vals.length) continue;
      const w = new Map(vals.map((v) => [v.id, v.w]));
      const forced = [];
      for (const l of g.links) {
        if (l.b.c !== cid || isInner(l) || !active.has(l.id) || !w.has(l.b.v)) continue;
        if (!(picked.get(l.a.c) || new Set()).has(l.a.v)) continue;
        if (l.always) { if (!forced.includes(l.b.v)) forced.push(l.b.v); }
        else w.set(l.b.v, w.get(l.b.v) * l.f);
      }
      let n = 1;
      if (c.multi) {
        const lo = Math.min(c.min, c.max);
        const hi = Math.max(c.min, c.max);
        n = lo + Math.floor(Math.random() * (hi - lo + 1));
      }
      n = Math.min(n, vals.length);
      // Links inside the characteristic work both ways (values of one roll have no order):
      // "timid — brave: never" also means "brave — timid: never".
      const inner = g.links.filter((l) => isInner(l) && l.a.c === cid && active.has(l.id));
      const chosen = [];
      const take = (id) => {
        chosen.push(id);
        for (const l of inner) {
          const other = l.a.v === id ? l.b.v : l.b.v === id ? l.a.v : null;
          if (!other || chosen.includes(other) || !w.has(other)) continue;
          if (l.always) { if (!forced.includes(other)) forced.push(other); }
          else w.set(other, w.get(other) * l.f);
        }
      };
      // forced values first; if there are more than fit, pick among them by weight
      while (chosen.length < n) {
        const queue = forced.filter((id) => !chosen.includes(id));
        const pool = queue.length
          ? queue.map((id) => ({ id, w: Math.max(w.get(id), 1e-9) }))
          : vals.filter((v) => !chosen.includes(v.id)).map((v) => ({ id: v.id, w: w.get(v.id) }));
        const i = weightedPick(pool);
        if (i < 0) break; // everything left has weight 0
        take(pool[i].id);
      }
      picked.set(cid, new Set(chosen));
      out.push({
        name: c.name.trim() || t('rng.untitledChar'),
        values: vals.filter((v) => chosen.includes(v.id)).map((v) => v.text.trim()),
        cid,
      });
    }
    // show in the order of the list, not the order of rolling
    const idx = new Map(g.chars.map((c, i) => [c.id, i]));
    return out.sort((a, b) => idx.get(a.cid) - idx.get(b.cid));
  }

  // ---------- page skeleton ----------

  const el = {};
  let past = []; // earlier batches [{ title, rolls }]
  let lastBatch = null;
  let linkFrom = null; // { c, v } while choosing the target of a new link
  let editLink = null; // id of the link whose editor is open on the map

  function build() {
    el.file = h('input', { type: 'file', accept: '.json,application/json', hidden: true, onchange: importFile });
    // Before a generator exists the name field is "what to create" and Enter / the button
    // creates it; afterwards the field renames it and the button starts over.
    el.name = h('input', {
      class: 'input rng-name', type: 'text', maxLength: 120, placeholder: t('rng.namePh'),
      'aria-label': t('rng.name'),
      oninput: () => { el.name.classList.remove('is-missing'); if (gen) { gen.name = el.name.value; save(); } },
      onkeydown: (e) => { if (e.key === 'Enter' && !gen) createGen(); },
    });
    el.toast = h('p', { class: 'rng-toast', role: 'status', 'aria-live': 'polite' });
    el.mainBtn = h('button', { class: 'btn btn--primary rng-main__btn', type: 'button', onclick: () => (gen ? startOver() : createGen()) });

    const files = h('div', { class: 'rng-files' },
      h('button', { class: 'rng-small', type: 'button', onclick: () => el.file.click() }, t('rng.import')),
      el.exportBtn = h('button', { class: 'rng-small', type: 'button', onclick: exportFile }, t('rng.export')),
      el.file,
    );
    const main = h('div', { class: 'rng-main' },
      h('label', { class: 'rng-main__label', text: t('rng.name') }),
      el.name,
      el.mainBtn,
      menu(t('rng.examples'), fillExamples, 'rng-main__ex'),
    );
    main.firstChild.htmlFor = 'rng-name';
    el.name.id = 'rng-name';
    el.start = h('p', { class: 'rng-help rng-start' }, t('rng.start'));
    const toolbar = h('div', { class: 'rng-toolbar' }, files, main, el.start, el.toast);

    // sharing comes after the work is done: under the map of traits
    el.shareRow = h('div', { class: 'rng-share' },
      h('button', { class: 'btn', type: 'button', onclick: shareLink }, t('rng.share')),
    );

    el.canvas = h('div', { class: 'rng-canvas' });
    el.svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    el.svg.setAttribute('class', 'rng-edges');
    el.svg.setAttribute('aria-hidden', 'true');
    el.board = h('div', { class: 'rng-board' }, el.canvas);
    el.banner = h('div', { class: 'rng-banner', hidden: true });

    const mapHead = h('div', { class: 'rng-map-head' },
      h('h2', { class: 'panel__title' }, t('rng.map')),
      h('div', { class: 'rng-map-head__btns' },
        h('button', { class: 'btn', type: 'button', onclick: addChar }, t('rng.addChar')),
        h('button', { class: 'btn', type: 'button', onclick: () => { layout(gen, true); gen.links.forEach((l) => { l.ox = 0; l.oy = 0; }); renderMap(); save(); } }, t('rng.arrange')),
      ),
    );

    el.linksBody = h('div', { class: 'rng-links' });
    el.linksSum = h('summary', {});
    const links = h('details', { class: 'rng-links-wrap', open: !matchMedia('(max-width: 560px)').matches }, el.linksSum, el.linksBody);

    el.count = h('input', {
      class: 'input rng-count', type: 'number', min: 1, max: MAX_COUNT, step: 1,
      onchange: () => { gen.count = clampInt(el.count.value, 1, MAX_COUNT); el.count.value = gen.count; save(); },
    });
    el.results = h('div', { class: 'rng-results', 'aria-live': 'polite' });
    el.history = h('div', { class: 'rng-history' });
    const rollPanel = el.rollPanel = h('section', { class: 'panel rng-roll' },
      h('div', { class: 'rng-roll__row' },
        h('button', { class: 'btn btn--primary rng-roll__btn', type: 'button', onclick: roll }, t('rng.roll')),
        h('label', { class: 'rng-roll__count' }, t('rng.count'), el.count),
        el.copyBtn = h('button', { class: 'btn btn--ghost', type: 'button', hidden: true, onclick: copyResults }, t('rng.copy')), // shown once there is a result
      ),
      el.results,
      el.history,
    );

    root.append(
      h('section', { class: 'panel rng-top' }, toolbar),
      el.mapPanel = h('section', { class: 'panel rng-map' },
        mapHead,
        h('p', { class: 'rng-help' }, t('rng.mapHelp')),
        el.banner,
        el.board,
        links,
        el.shareRow,
      ),
      rollPanel,
    );

    el.board.addEventListener('pointerdown', onBoardPointerDown, true);
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape') { stopLinking(); closeEditor(); closeMenus(); } });
    document.addEventListener('click', (e) => { if (!e.target.closest('.rng-menu__panel')) closeMenus(); });
    window.addEventListener('resize', () => renderEdges());
  }

  // ---------- top menus ----------

  // A button that opens a small list under it; only one list is open at a time
  function menu(label, fill, cls = '', aria) {
    const panel = h('div', { class: 'rng-menu__panel', hidden: true });
    const btn = h('button', {
      class: `btn rng-menu__btn ${cls}`, type: 'button', 'aria-haspopup': 'true', 'aria-expanded': 'false',
      'aria-label': aria,
      onclick: (e) => {
        e.stopPropagation();
        const opening = panel.hidden;
        closeMenus();
        if (!opening) return;
        fill(panel);
        panel.hidden = false;
        btn.setAttribute('aria-expanded', 'true');
      },
    }, label);
    const wrap = h('div', { class: 'rng-menu' }, btn, panel);
    wrap.btn = btn;
    return wrap;
  }

  function closeMenus() {
    root.querySelectorAll('.rng-menu__panel').forEach((p) => { p.hidden = true; });
    root.querySelectorAll('.rng-menu__btn').forEach((b) => b.setAttribute('aria-expanded', 'false'));
  }

  function menuItem(title, desc, onclick, cls = '') {
    return h('button', {
      class: `rng-menu__item ${cls}`, type: 'button',
      onclick: () => { closeMenus(); onclick(); },
    }, h('span', { class: 'rng-menu__title' }, title), desc ? h('span', { class: 'rng-menu__desc' }, desc) : null);
  }

  const traitList = (names) => names.filter(Boolean).join(', ').toLowerCase();

  function fillExamples(panel) {
    panel.replaceChildren(
      ...EXAMPLES.map((ex) => menuItem(
        pickLang(ex.name),
        traitList(ex.chars.map((c) => pickLang(c[0]))),
        () => addExample(ex),
      )),
    );
  }

  function renderAll() {
    // until a generator is chosen, settings would not belong to anything: show only the start
    const none = !gen;
    for (const n of [el.mapPanel, el.rollPanel, el.exportBtn]) n.hidden = none;
    el.start.hidden = !none;
    el.mainBtn.textContent = none ? t('rng.create') : `+ ${t('rng.new')}`;
    el.mainBtn.classList.toggle('btn--primary', none);
    if (none) { el.name.value = ''; return; }
    el.name.value = gen.name;
    el.count.value = gen.count;
    past = []; lastBatch = null; el.copyBtn.hidden = true;
    el.results.replaceChildren();
    el.history.replaceChildren();
    stopLinking();
    editLink = null;
    if (gen.needsLayout) { delete gen.needsLayout; layout(gen, true); }
    renderMap();
  }

  // ---------- the current generator: replace, files, links ----------

  // New, example, file and link all replace the current generator at once (null = back to the start)
  function replaceGen(g, note) {
    gen = g;
    if (g) save();
    else { clearTimeout(saveTimer); try { localStorage.removeItem(STORE_KEY); } catch (e) { /* ignore */ } }
    renderAll();
    if (note) toast(note);
  }

  // "Generate": a new generator named after the field (a name is needed to start)
  function createGen() {
    const name = el.name.value.trim();
    if (!name) { el.name.classList.add('is-missing'); el.name.focus(); return; }
    replaceGen(newGen(name));
    el.canvas.querySelector('.rng-node__name')?.focus({ preventScroll: true });
  }

  // "+ New generator": back to the empty start
  function startOver() {
    replaceGen(null);
    el.name.focus();
  }

  function addExample(ex) { replaceGen(exampleGen(ex)); }

  function fileName() {
    const base = (gen.name.trim() || 'generator').replace(/[\\/:*?"<>|\s]+/g, '-').slice(0, 60);
    return `${base}.json`;
  }

  function exportFile() {
    const blob = new Blob([JSON.stringify(exportable(gen), null, 2) + '\n'], { type: 'application/json' });
    const a = h('a', { href: URL.createObjectURL(blob), download: fileName() });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }

  function importFile() {
    const f = el.file.files[0];
    el.file.value = '';
    if (!f) return;
    f.text().then((txt) => {
      const g = normalize(JSON.parse(txt));
      replaceGen(g, t('rng.imported', { name: g.name.trim() || t('rng.untitled') }));
    }).catch(() => toast(t('rng.importBad')));
  }

  // share link: #g=z<base64url(deflate-raw(json))>, or #g=j<base64url(json)> without CompressionStream
  const b64u = (bytes) => {
    let s = '';
    for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
    return btoa(s).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
  };
  const unb64u = (s) => {
    const bin = atob(s.replace(/-/g, '+').replace(/_/g, '/'));
    return Uint8Array.from(bin, (ch) => ch.charCodeAt(0));
  };
  async function pipe(bytes, stream) {
    const out = new Blob([bytes]).stream().pipeThrough(stream);
    return new Uint8Array(await new Response(out).arrayBuffer());
  }

  async function shareLink() {
    const json = new TextEncoder().encode(JSON.stringify(exportable(gen)));
    let code;
    try { code = 'z' + b64u(await pipe(json, new CompressionStream('deflate-raw'))); }
    catch (e) { code = 'j' + b64u(json); }
    const url = `${location.origin}${location.pathname}#g=${code}`;
    try {
      await navigator.clipboard.writeText(url);
      toast(t('rng.shareDone'));
    } catch (e) {
      toast(t('rng.shareFail'), url);
    }
  }

  async function openSharedLink() {
    const m = /^#g=([zj])([A-Za-z0-9_-]+)$/.exec(location.hash);
    if (!m) return;
    window.history.replaceState(null, '', location.pathname + location.search);
    try {
      let bytes = unb64u(m[2]);
      if (m[1] === 'z') bytes = await pipe(bytes, new DecompressionStream('deflate-raw'));
      const g = normalize(JSON.parse(new TextDecoder().decode(bytes)));
      replaceGen(g, t('rng.imported', { name: g.name.trim() || t('rng.untitled') }));
    } catch (e) {
      toast(t('rng.linkBad'));
    }
  }

  let toastTimer = 0;
  function toast(msg, url) {
    el.toast.replaceChildren(msg, url ? h('input', { class: 'input rng-toast__url', value: url, readOnly: true, onfocus: (e) => e.target.select() }) : '');
    el.toast.classList.add('is-on');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => el.toast.classList.remove('is-on'), url ? 20000 : 3500);
  }

  // ---------- map: layout ----------

  function nodeHeight(c) { return 118 + c.values.length * 42 + (c.multi ? 34 : 0); }

  // Simple column layout: fill as many columns as fit the board (at least 2)
  function layout(g, all) {
    const width = (el.board && el.board.clientWidth) || 900;
    const cols = Math.max(2, Math.floor((width - LEFT) / (NODE_W + GAP_X)));
    const colY = new Array(cols).fill(16);
    g.chars.forEach((c, i) => {
      if (!all && (c.x || c.y)) return;
      const col = all ? i % cols : colY.indexOf(Math.min(...colY));
      c.x = LEFT + col * (NODE_W + GAP_X);
      c.y = colY[col];
      colY[col] += nodeHeight(c) + GAP_Y;
    });
  }

  function placeNew(c) {
    // below the lowest card in the column that is least filled
    const width = el.board.clientWidth || 900;
    const cols = Math.max(2, Math.floor((width - LEFT) / (NODE_W + GAP_X)));
    let best = { col: 0, y: Infinity };
    for (let col = 0; col < cols; col++) {
      const x = LEFT + col * (NODE_W + GAP_X);
      const bottom = gen.chars.filter((o) => o !== c && Math.abs(o.x - x) < NODE_W)
        .reduce((m, o) => Math.max(m, o.y + nodeHeight(o) + GAP_Y), 16);
      if (bottom < best.y) best = { col, y: bottom };
    }
    c.x = LEFT + best.col * (NODE_W + GAP_X);
    c.y = best.y;
  }

  // ---------- map: rendering ----------

  function totalWeight(c) { return c.values.reduce((s, v) => s + (v.text.trim() ? v.w : 0), 0); }
  function pct(c, v) {
    const tw = totalWeight(c);
    if (!v.text.trim() || !tw) return '';
    const p = (v.w / tw) * 100;
    return `${p >= 10 || p === 0 ? Math.round(p) : p.toFixed(1)}%`;
  }

  function renderMap() {
    const keepX = el.board.scrollLeft;
    const keepY = el.board.scrollTop;
    el.canvas.replaceChildren(el.svg);
    for (const c of gen.chars) el.canvas.append(renderNode(c));
    sizeCanvas();
    renderEdges();
    renderLinkList();
    el.board.scrollLeft = keepX;
    el.board.scrollTop = keepY;
  }

  function sizeCanvas() {
    let w = 0;
    let hgt = 0;
    for (const n of el.canvas.querySelectorAll('.rng-node')) {
      w = Math.max(w, n.offsetLeft + n.offsetWidth);
      hgt = Math.max(hgt, n.offsetTop + n.offsetHeight);
    }
    el.canvas.style.width = `${Math.max(w + 140, el.board.clientWidth)}px`;
    el.canvas.style.height = `${Math.max(hgt + 120, 360)}px`;
  }

  function renderNode(c) {
    const node = h('div', { class: 'rng-node', 'data-c': c.id });
    node.style.left = `${c.x}px`;
    node.style.top = `${c.y}px`;

    const nameIn = h('input', {
      class: 'rng-node__name', type: 'text', maxLength: 120, value: c.name, placeholder: t('rng.charPh'),
      'aria-label': t('rng.charPh'),
      oninput: () => { c.name = nameIn.value; save(); },
      onchange: () => renderLinkList(),
    });
    const head = h('div', { class: 'rng-node__head', title: t('rng.drag') },
      h('span', { class: 'rng-node__grip', 'aria-hidden': 'true' }, '⠿'),
      nameIn,
      h('button', {
        class: 'rng-x', type: 'button', title: t('rng.delChar'), 'aria-label': t('rng.delChar'),
        onclick: () => {
          gen.chars = gen.chars.filter((o) => o !== c);
          gen.links = gen.links.filter((l) => l.a.c !== c.id && l.b.c !== c.id);
          save();
          renderMap();
        },
      }, '×'),
    );
    head.addEventListener('pointerdown', (e) => startDrag(e, c, node));

    const list = h('div', { class: 'rng-vals' });
    c.values.forEach((v) => list.append(renderValue(c, v, node)));

    const multiBox = h('input', {
      type: 'checkbox', checked: c.multi,
      onchange: () => { c.multi = multiBox.checked; save(); renderMap(); },
    });
    const minIn = h('input', {
      class: 'rng-num', type: 'number', min: 1, max: 99, value: c.min, 'aria-label': t('rng.from'),
      onchange: () => { c.min = clampInt(minIn.value, 1, 99); minIn.value = c.min; save(); },
    });
    const maxIn = h('input', {
      class: 'rng-num', type: 'number', min: 1, max: 99, value: c.max, 'aria-label': t('rng.to'),
      onchange: () => { c.max = clampInt(maxIn.value, 1, 99); maxIn.value = c.max; save(); },
    });

    const foot = h('div', { class: 'rng-node__foot' },
      h('button', {
        class: 'rng-add', type: 'button',
        onclick: () => {
          const v = newValue();
          c.values.push(v);
          save();
          renderMap();
          el.canvas.querySelector(`.rng-val[data-v="${v.id}"] .rng-val__text`)?.focus();
        },
      }, t('rng.addVal')),
      h('label', { class: 'rng-multi' }, multiBox, t('rng.multi')),
      c.multi ? h('div', { class: 'rng-range' }, t('rng.from'), ' ', minIn, ' ', t('rng.to'), ' ', maxIn) : null,
    );

    node.append(head, h('div', { class: 'rng-vals-head' }, h('span', {}, ''), h('span', {}, t('rng.weight'))), list, foot);
    return node;
  }

  function renderValue(c, v, node) {
    const share = h('span', { class: 'rng-val__pct' }, pct(c, v));
    const refreshPct = () => {
      for (const row of node.querySelectorAll('.rng-val')) {
        const val = c.values.find((x) => x.id === row.dataset.v);
        if (val) row.querySelector('.rng-val__pct').textContent = pct(c, val);
      }
    };
    const textIn = h('input', {
      class: 'rng-val__text', type: 'text', maxLength: 300, value: v.text, placeholder: t('rng.valPh'),
      'aria-label': t('rng.valPh'),
      oninput: () => { v.text = textIn.value; save(); refreshPct(); },
      onchange: () => renderLinkList(),
    });
    const wIn = h('input', {
      class: 'rng-num rng-val__w', type: 'number', min: 0, step: 1, value: v.w, 'aria-label': t('rng.weight'),
      oninput: () => { v.w = clampInt(wIn.value, 0, 1e6); save(); refreshPct(); },
      onchange: () => { wIn.value = v.w; },
    });
    const row = h('div', { class: 'rng-val', 'data-c': c.id, 'data-v': v.id },
      textIn,
      wIn,
      share,
      h('button', {
        class: 'rng-link-btn', type: 'button', title: t('rng.linkBtn'), 'aria-label': t('rng.linkBtn'),
        onclick: () => startLinking(c, v),
      }, '⟶'),
      h('button', {
        class: 'rng-x', type: 'button', title: t('rng.delVal'), 'aria-label': t('rng.delVal'),
        onclick: () => {
          c.values = c.values.filter((o) => o !== v);
          gen.links = gen.links.filter((l) => l.a.v !== v.id && l.b.v !== v.id);
          save();
          renderMap();
        },
      }, '×'),
    );
    return row;
  }

  // ---------- map: links ----------

  function linkKind(l) {
    if (l.always) return 'always';
    if (l.f === 0) return 'never';
    if (l.f < 1) return 'less';
    if (l.f > 1) return 'more';
    return 'same';
  }
  function linkLabel(l) {
    if (l.always) return t('rng.always');
    if (l.f === 0) return t('rng.never');
    return `×${String(+l.f.toFixed(3)).replace('.', LANG === 'en' ? '.' : ',')}`;
  }

  function findVal(e) {
    const c = gen.chars.find((x) => x.id === e.c);
    const v = c && c.values.find((x) => x.id === e.v);
    return { c, v };
  }
  function endText(e) {
    const { c, v } = findVal(e);
    return `${(c && c.name.trim()) || t('rng.untitledChar')}: ${(v && v.text.trim()) || '…'}`;
  }
  // "if Coat: grey → Eyes: green" / "Temper: timid ↔ brave" (inner links work both ways)
  function linkText(l) {
    if (!isInner(l)) return `${t('rng.if')} ${endText(l.a)} → ${endText(l.b)}`;
    const b = findVal(l.b).v;
    return `${endText(l.a)} ↔ ${(b && b.text.trim()) || '…'}`;
  }

  function renderEdges() {
    if (!gen) return;
    const svg = el.svg;
    const W = el.canvas.offsetWidth;
    const H = el.canvas.offsetHeight;
    svg.setAttribute('width', W);
    svg.setAttribute('height', H);
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.innerHTML = '<defs>' + ['more', 'less', 'never', 'always', 'same'].map((k) =>
      `<marker id="rng-arrow-${k}" class="rng-arrow rng-arrow--${k}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 L3 5 Z"/></marker>`).join('') + '</defs>';
    el.canvas.querySelectorAll('.rng-edge-label, .rng-editor').forEach((n) => n.remove());

    const active = activeLinks(gen);
    const box = el.canvas.getBoundingClientRect();
    const rect = (sel) => {
      const n = el.canvas.querySelector(sel);
      if (!n) return null;
      const r = n.getBoundingClientRect();
      return { l: r.left - box.left, r: r.right - box.left, t: r.top - box.top, b: r.bottom - box.top };
    };

    for (const l of gen.links) {
      const na = rect(`.rng-node[data-c="${l.a.c}"]`);
      const nb = rect(`.rng-node[data-c="${l.b.c}"]`);
      const ra = rect(`.rng-val[data-v="${l.a.v}"]`);
      const rb = rect(`.rng-val[data-v="${l.b.v}"]`);
      if (!na || !nb || !ra || !rb) continue;
      const y1 = (ra.t + ra.b) / 2;
      const y2 = (rb.t + rb.b) / 2;
      let x1, x2, c1, c2;
      if (nb.l > na.r) { x1 = na.r; x2 = nb.l; const d = Math.max(40, (x2 - x1) / 2); c1 = x1 + d; c2 = x2 - d; }
      else if (nb.r < na.l) { x1 = na.l; x2 = nb.r; const d = Math.max(40, (x1 - x2) / 2); c1 = x1 - d; c2 = x2 + d; }
      else if (isInner(l)) { x1 = na.l; x2 = na.l; const d = Math.min(56, 22 + Math.abs(y2 - y1) * 0.3); c1 = x1 - d; c2 = c1; } // inner links loop on the left
      else { x1 = na.r; x2 = nb.r; const d = 50 + Math.abs(y2 - y1) * 0.15; c1 = Math.max(x1, x2) + d; c2 = c1; }
      const kind = linkKind(l);
      const on = active.has(l.id);
      const geo = { x1, y1, x2, y2, c1, c2 };
      const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      path.setAttribute('class', `rng-edge rng-edge--${kind}${on ? '' : ' is-off'}`);
      path.setAttribute('marker-end', `url(#rng-arrow-${kind})`);
      if (isInner(l)) path.setAttribute('marker-start', `url(#rng-arrow-${kind})`);
      path.dataset.l = l.id;
      svg.append(path);

      const lab = h('button', {
        class: `rng-edge-label rng-edge-label--${kind}${on ? '' : ' is-off'}${editLink === l.id ? ' is-open' : ''}`,
        type: 'button',
        title: `${linkText(l)}${on ? '' : ` (${offReason(l)})`}\n${t('rng.labelDrag')}`,
        onclick: (e) => e.stopPropagation(), // opening is done in dragLabel (a click that did not move)
        ondblclick: (e) => { e.stopPropagation(); l.ox = 0; l.oy = 0; placeEdge(l, geo, path, lab); save(); },
        onmouseenter: () => markEnds(l, true),
        onmouseleave: () => markEnds(l, false),
      }, linkLabel(l));
      lab.addEventListener('pointerdown', (e) => dragLabel(e, l, geo, path, lab));
      el.canvas.append(lab);
      const { mx, my } = placeEdge(l, geo, path, lab);

      if (editLink === l.id) {
        const ed = h('div', { class: 'rng-editor', onpointerdown: (e) => e.stopPropagation() },
          h('p', { class: 'rng-editor__what' }, linkText(l)),
          linkControls(l),
          on ? null : h('p', { class: 'rng-editor__warn' }, offReason(l)),
          isInner(l) ? h('p', { class: 'rng-editor__help' }, t('rng.innerHelp')) : null,
          h('p', { class: 'rng-editor__help' }, t('rng.factorHelp')),
        );
        ed.style.left = `${Math.min(mx, W - 270)}px`;
        ed.style.top = `${my + 18}px`;
        el.canvas.append(ed);
      }
    }
  }

  // The label sits at the middle of the bezier, moved by the link's own offset (ox, oy).
  // B(0.5) = (P0 + 3·C1 + 3·C2 + P3) / 8, so shifting both control points by δ moves
  // the middle by 0.75·δ: the curve keeps passing through the dragged label.
  function placeEdge(l, g, path, lab) {
    const ox = l.ox || 0;
    const oy = l.oy || 0;
    const dx = ox / 0.75;
    const dy = oy / 0.75;
    path.setAttribute('d', `M${g.x1} ${g.y1} C${g.c1 + dx} ${g.y1 + dy}, ${g.c2 + dx} ${g.y2 + dy}, ${g.x2} ${g.y2}`);
    const mx = (g.x1 + 3 * g.c1 + 3 * g.c2 + g.x2) / 8 + ox;
    const my = (g.y1 + 3 * g.y1 + 3 * g.y2 + g.y2) / 8 + oy;
    lab.style.left = `${mx}px`;
    lab.style.top = `${my}px`;
    return { mx, my };
  }

  // Drag a label to bend its arrow out of the crowd; a click without moving opens the editor
  function dragLabel(e, l, geo, path, lab) {
    if (e.button !== 0 || linkFrom) return;
    e.preventDefault();
    e.stopPropagation();
    const sx = e.clientX;
    const sy = e.clientY;
    const ox = l.ox || 0;
    const oy = l.oy || 0;
    let moved = false;
    lab.setPointerCapture?.(e.pointerId);
    const move = (ev) => {
      if (!moved && Math.hypot(ev.clientX - sx, ev.clientY - sy) < 4) return;
      if (!moved) { moved = true; lab.classList.add('is-dragging'); el.canvas.querySelector('.rng-editor')?.remove(); }
      l.ox = ox + ev.clientX - sx;
      l.oy = oy + ev.clientY - sy;
      placeEdge(l, geo, path, lab);
    };
    const up = () => {
      lab.removeEventListener('pointermove', move);
      lab.removeEventListener('pointerup', up);
      lab.removeEventListener('pointercancel', up);
      lab.classList.remove('is-dragging');
      if (moved) { l.ox = Math.round(l.ox); l.oy = Math.round(l.oy); save(); if (editLink === l.id) renderEdges(); return; }
      editLink = editLink === l.id ? null : l.id;
      renderEdges();
    };
    lab.addEventListener('pointermove', move);
    lab.addEventListener('pointerup', up);
    lab.addEventListener('pointercancel', up);
  }

  function markEnds(l, on) {
    for (const id of [l.a.v, l.b.v]) el.canvas.querySelector(`.rng-val[data-v="${id}"]`)?.classList.toggle('is-lit', on);
  }

  // factor input + presets + delete; used on the map and in the list
  function linkControls(l) {
    const num = h('input', {
      class: 'rng-num rng-factor', type: 'number', min: 0, step: 0.1, value: l.f, disabled: l.always,
      'aria-label': t('rng.factor'),
      onchange: () => { l.f = Math.max(0, Math.min(1000, Number(String(num.value).replace(',', '.')) || 0)); l.always = false; changed(); },
    });
    const preset = (label, f, always) => h('button', {
      class: `rng-chip${(always ? l.always : !l.always && l.f === f) ? ' is-on' : ''}`, type: 'button',
      onclick: () => { l.always = always; if (!always) l.f = f; changed(); },
    }, label);
    return h('div', { class: 'rng-ctl' },
      preset(t('rng.never'), 0, false),
      preset(t('rng.rarer'), 0.5, false),
      preset(t('rng.often'), 2, false),
      preset(t('rng.always'), l.f, true),
      h('label', { class: 'rng-ctl__f' }, '×', num),
      h('button', {
        class: 'rng-x rng-x--link', type: 'button', title: t('rng.delLink'), 'aria-label': t('rng.delLink'),
        onclick: () => { gen.links = gen.links.filter((o) => o !== l); if (editLink === l.id) editLink = null; changed(); },
      }, '✕'),
    );
  }

  function changed() {
    save();
    renderEdges();
    renderLinkList();
  }

  function renderLinkList() {
    const active = activeLinks(gen);
    el.linksSum.textContent = `${t('rng.links')} (${gen.links.length})`;
    if (!gen.links.length) {
      el.linksBody.replaceChildren(h('p', { class: 'rng-help' }, t('rng.noLinks')));
      return;
    }
    el.linksBody.replaceChildren(...gen.links.map((l) => h('div', { class: `rng-link-row${active.has(l.id) ? '' : ' is-off'}` },
      h('p', { class: 'rng-link-row__what' },
        linkText(l),
        active.has(l.id) ? null : h('small', {}, ` (${offReason(l)})`)),
      linkControls(l),
    )));
  }

  function startLinking(c, v) {
    if (linkFrom && linkFrom.v === v.id) { stopLinking(); return; }
    closeEditor();
    linkFrom = { c: c.id, v: v.id };
    el.canvas.classList.add('is-linking');
    el.canvas.querySelectorAll('.rng-val').forEach((row) => {
      row.classList.toggle('is-source', row.dataset.v === v.id);
      row.classList.toggle('is-target', row.dataset.c !== c.id || (c.multi && row.dataset.v !== v.id));
    });
    el.banner.replaceChildren(
      h('span', {}, t('rng.linking', { v: v.text.trim() || '…' })),
      h('button', { class: 'btn btn--ghost', type: 'button', onclick: stopLinking }, t('rng.cancel')),
    );
    el.banner.hidden = false;
  }

  function stopLinking() {
    linkFrom = null;
    if (!el.canvas) return;
    el.canvas.classList.remove('is-linking');
    el.canvas.querySelectorAll('.is-source, .is-target').forEach((n) => n.classList.remove('is-source', 'is-target'));
    el.banner.hidden = true;
  }

  function closeEditor() {
    if (editLink) { editLink = null; renderEdges(); }
  }

  // A link may go to a value of another characteristic, or to another value
  // of the same characteristic if it gives several values
  function canTarget(row) {
    if (row.dataset.c !== linkFrom.c) return true;
    const c = gen.chars.find((x) => x.id === linkFrom.c);
    return !!(c && c.multi && row.dataset.v !== linkFrom.v);
  }

  // While linking, a tap on a value that can be linked finishes the link.
  // Otherwise a tap outside the open editor closes it.
  function onBoardPointerDown(e) {
    if (linkFrom) {
      const row = e.target.closest('.rng-val');
      if (row && canTarget(row)) {
        e.preventDefault();
        e.stopPropagation();
        const b = { c: row.dataset.c, v: row.dataset.v };
        let l = gen.links.find((x) => (x.a.v === linkFrom.v && x.b.v === b.v) || (x.a.c === x.b.c && x.b.v === linkFrom.v && x.a.v === b.v));
        if (!l) {
          l = { id: uid(), a: { ...linkFrom }, b, f: 2, always: false };
          gen.links.push(l);
        }
        stopLinking();
        editLink = l.id;
        changed();
      }
      return;
    }
    if (editLink && !e.target.closest('.rng-editor, .rng-edge-label')) closeEditor();
  }

  // ---------- map: dragging cards ----------

  function startDrag(e, c, node) {
    if (e.button !== 0 || linkFrom || e.target.closest('button')) return;
    // the name field takes most of the header: dragging it moves the card,
    // a plain click (no movement) puts the cursor into it
    const nameIn = e.target.closest('.rng-node__name');
    if (nameIn && document.activeElement === nameIn) return;
    e.preventDefault();
    const sx = e.clientX;
    const sy = e.clientY;
    const ox = c.x;
    const oy = c.y;
    let moved = false;
    if (node.nextElementSibling) el.canvas.append(node); // the card you hold lies on top
    node.setPointerCapture?.(e.pointerId);
    let raf = 0;
    const move = (ev) => {
      if (!moved && Math.hypot(ev.clientX - sx, ev.clientY - sy) < 4) return;
      if (!moved) { moved = true; node.classList.add('is-dragging'); }
      c.x = Math.max(0, ox + ev.clientX - sx);
      c.y = Math.max(0, oy + ev.clientY - sy);
      node.style.left = `${c.x}px`;
      node.style.top = `${c.y}px`;
      if (!raf) raf = requestAnimationFrame(() => { raf = 0; sizeCanvas(); renderEdges(); });
    };
    const up = () => {
      node.classList.remove('is-dragging');
      node.removeEventListener('pointermove', move);
      node.removeEventListener('pointerup', up);
      node.removeEventListener('pointercancel', up);
      if (moved) save();
      else if (nameIn) nameIn.focus();
    };
    node.addEventListener('pointermove', move);
    node.addEventListener('pointerup', up);
    node.addEventListener('pointercancel', up);
  }

  function addChar() {
    const c = newChar();
    gen.chars.push(c);
    placeNew(c);
    save();
    renderMap();
    const n = el.canvas.querySelector(`.rng-node[data-c="${c.id}"]`);
    if (n) {
      el.board.scrollTo({ left: Math.max(0, c.x - 16), top: Math.max(0, c.y - 16), behavior: 'smooth' });
      n.querySelector('.rng-node__name').focus({ preventScroll: true });
    }
  }

  // ---------- rolling and results ----------

  function title() { return gen.name.trim() || t('rng.untitled'); }

  function roll() {
    gen.count = clampInt(el.count.value, 1, MAX_COUNT);
    el.count.value = gen.count;
    save();
    const usable = gen.chars.some((c) => c.values.some((v) => v.text.trim()));
    if (!usable) {
      el.results.replaceChildren(h('p', { class: 'rng-help' }, t('rng.empty')));
      return;
    }
    const active = activeLinks(gen);
    const rolls = Array.from({ length: gen.count }, () => rollOnce(gen, active));
    const prev = el.results.querySelector('.rng-cards');
    if (prev && lastBatch) {
      past.unshift(lastBatch);
      past.length = Math.min(past.length, MAX_HISTORY);
    }
    lastBatch = { title: title(), rolls };
    el.results.replaceChildren(cards(lastBatch, true));
    el.copyBtn.hidden = false;
    renderHistory();
  }

  function cards(batch, fresh) {
    return h('div', { class: `rng-cards${fresh ? ' is-fresh' : ''}` }, batch.rolls.map((r, i) =>
      h('article', { class: 'rng-card' },
        h('h3', { class: 'rng-card__title' }, batch.rolls.length > 1 ? `${batch.title} · ${i + 1}` : batch.title),
        h('dl', { class: 'rng-card__list' }, r.map((x) => [
          h('dt', {}, x.name),
          h('dd', {}, x.values.length ? x.values.join(', ') : t('rng.none')),
        ])),
      )));
  }

  function renderHistory() {
    if (!past.length) { el.history.replaceChildren(); return; }
    el.history.replaceChildren(h('details', { class: 'rng-history__box' },
      h('summary', {}, `${t('rng.history')} (${past.length})`),
      past.map((b) => cards(b, false)),
      h('button', {
        class: 'btn btn--ghost', type: 'button',
        onclick: () => { past.length = 0; renderHistory(); },
      }, t('rng.clear')),
    ));
  }

  function copyResults() {
    const b = lastBatch;
    if (!b) return;
    const txt = b.rolls.map((r, i) => [
      b.rolls.length > 1 ? `${b.title} ${i + 1}` : b.title,
      ...r.map((x) => `${x.name}: ${x.values.join(', ') || '—'}`),
    ].join('\n')).join('\n\n');
    navigator.clipboard?.writeText(txt).then(() => toast(t('rng.copied')), () => toast(t('rng.shareFail'), txt));
  }

  // ---------- start ----------

  build();
  renderAll();
  openSharedLink();
  window.addEventListener('hashchange', openSharedLink);
})();
