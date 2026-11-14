/* texts.js — texts.html: a desk with untidy stacks of books.

   Every section is a book with its name on the spine; nameless books are mixed in
   to make a heap. The heap is laid out anew on every visit (and on resize), a bit
   differently each time. Hover pulls a book out of its stack a little; a click
   pulls it out and opens the section: texts.html?s=<section id>.
   A section page is a stub for now ("this book is still empty").
*/
(function () {
  'use strict';

  const t = (key, fallback) => (window.I18n && I18n.has(key) ? I18n.t(key) : fallback);

  // To add a section, add it here (and its name to i18n.js as texts.<id>)
  const SECTIONS = [
    { id: 'stories', label: t('texts.stories', 'Рассказы') },
    { id: 'tales', label: t('texts.tales', 'Сказки') },
    { id: 'repatriate', label: t('texts.repatriate', 'Записки репатрианта') },
    { id: 'dreams', label: t('texts.dreams', 'Сны') },
    { id: 'chatter', label: t('texts.chatter', 'Болтовня') },
    { id: 'other', label: t('texts.other', 'Другое') },
  ];

  const COLORS = ['#4a1f24', '#1f3a2c', '#1d2742', '#2e2233', '#3d3020', '#16171d', '#40282f', '#24333a', '#3a2a1a', '#2b1c2e'];
  const FALL_MS = 700;        // the book falls off the desk
  const SPILL_AT_MS = 320;    // the book has tipped over the edge: now its pages fly out
  const TAKE_DELAY_MS = 1350; // the section opens while the leaves are still in the air
  const SPRITE_DIR = "/assets/img/texts/"; // <section id>/1.png, 2.png, ... (see README.md there)
  const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');

  const host = document.getElementById('desk');
  if (!host) return;

  const sectionId = new URLSearchParams(location.search).get('s');
  const section = SECTIONS.find((s) => s.id === sectionId);
  if (section) {
    showSection(section);
    return;
  }

  const rand = (a, b) => a + Math.random() * (b - a);
  const pick = (arr) => arr[Math.floor(Math.random() * arr.length)];

  // Text width for a spine title, so a labelled book is long enough for its name
  const measureCtx = document.createElement('canvas').getContext('2d');
  function textWidth(text, px) {
    const family = getComputedStyle(document.body).fontFamily;
    measureCtx.font = `600 ${px}px ${family}`;
    return measureCtx.measureText(text.toUpperCase()).width + text.length * px * 0.14; // + letter-spacing
  }

  let lastWidth = 0;
  let scale = 1;
  // titles are measured with the site font, so wait for it (but not for long)
  Promise.race([document.fonts ? document.fonts.ready : null, new Promise((r) => setTimeout(r, 600))])
    .then(() => layout());
  window.addEventListener('resize', () => {
    if (Math.abs(host.clientWidth - lastWidth) > 40) layout();
  });

  function layout() {
    const W = host.clientWidth;
    lastWidth = W;
    host.innerHTML = '';

    const k = Math.max(0.62, Math.min(1, W / 1000)); // size of everything
    scale = k;
    const titlePx = Math.round(13 * Math.max(0.85, k));
    const stacksN = W < 520 ? 2 : W < 820 ? 3 : 4;
    const deskH = Math.round(46 * k);

    // books: every section once, plus nameless ones for the heap
    const books = SECTIONS.map((s) => ({ section: s }));
    const decor = stacksN * 3 + Math.floor(rand(0, 3));
    for (let i = 0; i < decor; i++) books.push({ section: null });
    books.sort(() => Math.random() - 0.5);

    // one or two nameless ones stand leaning against a stack (behind it, so they
    // never cover a title)
    const leaning = [];
    const leanN = stacksN >= 3 ? 2 : 1;
    for (let i = books.length - 1; i >= 0 && leaning.length < leanN; i--) {
      if (!books[i].section) leaning.push(...books.splice(i, 1));
    }

    // spread the rest over the stacks
    const stacks = Array.from({ length: stacksN }, () => []);
    books.forEach((b, i) => stacks[i % stacksN].push(b));

    const slot = W / stacksN;
    let maxTop = 0;
    const els = [];

    // Height of everything already lying in a stack, every 2px across the desk (the
    // "skyline"). A new book is lowered until it touches it, so a tilted book rests
    // on the highest point under it instead of sinking into its neighbour. Every stack
    // has its own skyline: a book that pokes out over the next stack does not lie on it.
    const STEP = 2;
    const newSky = () => new Float32Array(Math.ceil(W / STEP) + 2).fill(deskH);
    const skyAt = (sky, x) => sky[Math.max(0, Math.min(sky.length - 1, Math.round(x / STEP)))];
    // Corners and edge points of a w x h box turned by rot degrees (CSS: clockwise),
    // relative to its centre, in "up is positive" coordinates.
    const edge = (w, h, rot, top) => {
      const th = -rot * Math.PI / 180, c = Math.cos(th), s = Math.sin(th);
      const y = top ? h / 2 : -h / 2;
      const pts = [];
      for (let i = 0; i <= 24; i++) {
        const x = -w / 2 + (w * i) / 24;
        pts.push([x * c - y * s, x * s + y * c]);
      }
      return pts;
    };
    // Lowest CSS `bottom` at which the book rests on the skyline; then raise the
    // skyline by the book's top edge.
    function rest(sky, left, w, h, rot) {
      const cx = left + w / 2;
      let b = deskH;
      for (const [x, y] of edge(w, h, rot, false)) b = Math.max(b, skyAt(sky, cx + x) - h / 2 - y);
      const top = edge(w, h, rot, true);
      for (let i = 0; i < top.length - 1; i++) {
        const [x0, y0] = top[i], [x1, y1] = top[i + 1];
        for (let x = cx + x0; x <= cx + x1; x += STEP) {
          const yy = b + h / 2 + y0 + ((x - cx - x0) / (x1 - x0 || 1)) * (y1 - y0);
          const j = Math.round(x / STEP);
          if (j >= 0 && j < sky.length) sky[j] = Math.max(sky[j], yy);
        }
        maxTop = Math.max(maxTop, b + h / 2 + y0);
      }
      return b;
    }

    stacks.forEach((stack, si) => {
      const cx = slot * (si + 0.5) + rand(-0.12, 0.12) * slot;
      // the stack keeps to its own lane of the desk: its books never reach into the
      // next stack, so no two stacks seem to pass through each other
      const laneL = slot * si + 5 * k, laneR = slot * (si + 1) - 5 * k;
      let widest = 0;
      const sky = (stack.sky = newSky());
      stack.forEach((b, level) => {
        const minLen = b.section ? textWidth(b.section.label, titlePx) + 64 * k : 0;
        const len = Math.min(laneR - laneL, Math.max(minLen, rand(170, 280) * k));
        const th = Math.round((b.section ? rand(44, 62) : rand(34, 64)) * k); // thick old volumes
        // now and then a book sticks far out of its stack
        const shift = (Math.random() < 0.18 ? rand(30, 48) * (Math.random() < 0.5 ? -1 : 1) : rand(-20, 20)) * k;
        const rot = rand(-3.5, 3.5) * (level === 0 ? 0.3 : 1);
        const left = Math.max(laneL, Math.min(laneR - len, cx - len / 2 + shift));
        const bottom = rest(sky, left, len, th, rot);
        els.push(bookEl(b, {
          left, bottom, w: len, h: th, rot, z: 10 + level,
          color: pick(COLORS), titlePx, stack: si, level, drop: th,
        }));
        widest = Math.max(widest, Math.abs(left - cx), Math.abs(left + len - cx));
      });
      stack.cx = cx;
      stack.half = widest;
    });

    // a nameless book thrown across the tops of two neighbouring stacks: it lies on
    // the highest points of both, tilted towards the lower one
    if (stacksN >= 2) {
      const i = Math.floor(rand(0, stacksN - 1));
      const a = stacks[i], b = stacks[i + 1];
      const x1 = a.cx + rand(-0.15, 0.15) * a.half, x2 = b.cx + rand(-0.15, 0.15) * b.half;
      // it lies on both stacks: their skylines together
      const both = a.sky.map((v, j) => Math.max(v, b.sky[j]));
      const dx = x2 - x1, dy = skyAt(both, x2) - skyAt(both, x1);
      const len = Math.hypot(dx, dy) + rand(40, 90) * k;
      const th = Math.round(rand(36, 50) * k);
      const rot = -Math.atan2(dy, dx) * 180 / Math.PI; // CSS rotation, y up
      const left = Math.max(2, Math.min(W - len - 2, (x1 + x2) / 2 - len / 2));
      const bottom = rest(both, left, len, th, rot);
      els.push(bookEl({ section: null }, {
        left, bottom, w: len, h: th, rot, z: 25, color: pick(COLORS), titlePx,
      }));
    }

    // leaning books: standing on the desk against the side of a stack
    leaning.forEach((b, i) => {
      const stack = stacks[(i * 2 + 1) % stacksN];
      const minLen = b.section ? textWidth(b.section.label, titlePx) + 64 * k : 0;
      const len = Math.max(minLen, rand(200, 250) * k);
      const th = Math.round(rand(44, 58) * k);
      const lean = rand(7, 14);
      // The book is turned up around its bottom-left corner (the pivot): standing,
      // it takes [pivot - th, pivot] on the desk; then it tips towards the stack.
      // It goes on the side of the stack where it fits on the desk.
      const pivotFor = (sd) => (sd > 0 ? stack.cx + stack.half + 6 * k + th : stack.cx - stack.half - 6 * k);
      let side = i % 2 ? -1 : 1; // +1: right of the stack, -1: left of it
      const fits = (sd) => pivotFor(sd) - th >= 0 && pivotFor(sd) + len * Math.sin(lean * Math.PI / 180) <= W;
      if (!fits(side) && fits(-side)) side = -side;
      const pivot = pivotFor(side);
      els.push(bookEl(b, {
        // it stands behind the stack, deeper on the desk: seen from a little above, its foot
        // is higher up and it is in shade; the stack in front covers it
        left: pivot, bottom: deskH + Math.round(rand(8, 14) * k), w: len, h: th, color: pick(COLORS), titlePx, z: 5, back: true,
        upright: true, rot: side > 0 ? -90 - lean : -90 + lean,
      }));
      maxTop = Math.max(maxTop, deskH + len);
    });

    const desk = document.createElement('div');
    desk.className = 'desk__top';
    desk.style.height = `${deskH}px`;
    host.appendChild(desk);
    // the top of the desk seen from a little above: a strip of wood going back, where
    // books further away stand (behind everything else)
    const surface = document.createElement('div');
    surface.className = 'desk__surface';
    surface.style.bottom = `${deskH}px`;
    surface.style.height = `${Math.round(22 * k)}px`;
    host.appendChild(surface);
    els.forEach((el) => host.appendChild(el));
    host.style.height = `${Math.ceil(maxTop + 40 * k)}px`;
    addWeb(W, k);
  }

  function bookEl(b, o) {
    const el = document.createElement(b.section ? 'a' : 'span');
    el.className = 'book' + (b.section ? ' book--section' : '') + (o.upright ? ' book--upright' : '') + (o.back ? ' book--back' : '');
    if (o.stack !== undefined) {
      // place in a stack, so the books above can settle when this one falls out
      el.dataset.stack = o.stack;
      el.dataset.level = o.level;
      el.dataset.drop = o.drop;
    }
    el.style.left = `${o.left.toFixed(1)}px`;
    el.style.bottom = `${o.bottom.toFixed(1)}px`;
    el.style.width = `${o.w.toFixed(1)}px`;
    el.style.height = `${o.h}px`;
    el.style.zIndex = o.z;
    el.style.setProperty('--c', o.color);
    el.style.setProperty('--r', `${o.rot.toFixed(2)}deg`);
    el.style.setProperty('--title', `${o.titlePx}px`);
    // wear: faded colour, its own patch of stains and scuffs, unevenly rubbed corners
    el.style.setProperty('--sat', rand(0.5, 0.88).toFixed(2));
    el.style.setProperty('--br', rand(0.7, 0.96).toFixed(2));
    el.style.setProperty('--bx', `${Math.round(rand(0, 240))}px`);
    el.style.setProperty('--by', `${Math.round(rand(0, 160))}px`);
    // Every book has its own past: most are a bit worn, some almost pristine,
    // a few badly knocked about. That sets the faded ends and the scratches.
    const life = Math.random();
    const kind = life < 0.3 ? 'kept' : life < 0.82 ? 'worn' : 'battered';
    const wear = kind === 'kept' ? rand(0.1, 0.3) : kind === 'worn' ? rand(0.4, 0.7) : rand(0.75, 1);
    el.style.setProperty('--wear', wear.toFixed(2));
    el.style.setProperty('--we', `${(kind === 'battered' ? rand(10, 18) : rand(5, 11)).toFixed(1)}%`);
    el.style.setProperty('--scratch', scratches(o.w, o.h, kind));
    // one rubbed shiny spot
    el.style.setProperty('--wx', `${rand(15, 85).toFixed(0)}%`);
    el.style.setProperty('--wy', `${rand(25, 75).toFixed(0)}%`);
    const c = () => Math.round(rand(1, 4)); // small: a flat spine, not a scroll
    el.style.setProperty('--rad', `${c()}px ${c()}px ${c()}px ${c()}px / ${c()}px ${c()}px ${c()}px ${c()}px`);
    if (o.upright) el.style.transformOrigin = '0 100%';
    if (b.section) {
      el.href = `/texts.html?s=${encodeURIComponent(b.section.id)}`;
      el.dataset.section = b.section.id;
      el.innerHTML = `<span class="book__title">${escapeHtml(b.section.label)}</span>`;
      el.addEventListener('click', (e) => {
        if (e.button !== 0 || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) return;
        if (reducedMotion.matches) return;
        e.preventDefault();
        if (el.classList.contains('is-falling')) return;
        fall(el);
        setTimeout(() => { location.href = el.href; }, TAKE_DELAY_MS);
      });
    } else {
      el.setAttribute('aria-hidden', 'true');
      // a nameless book may still have a worn band of gold in the middle
      if (Math.random() < 0.35) el.classList.add('book--plain-band');
    }
    if (!o.upright) {
      // the cover seen from a little above: a thin strip over the spine, slightly
      // narrower at the back. The book lying on top hides it, as it would.
      const top = document.createElement('span');
      top.className = 'book__top';
      top.style.setProperty('--top', `${Math.round(rand(5, 8) * scale)}px`);
      el.prepend(top);
    }
    return el;
  }

  // Scratches of one book, drawn for its size as a small SVG (a CSS url() value).
  // A scratch is a short straight nick (sometimes broken into two or three), pale with a
  // dark groove beside it, tapered at the ends and kept inside the spine; most run along
  // the spine (books slide against each other), some go any way. Badly worn books also
  // get scuffed spots: bunches of tiny parallel marks.
  function scratches(w, h, kind) {
    const n = kind === 'kept' ? Math.floor(rand(0, 2)) : kind === 'worn' ? Math.round(rand(2, 6)) : Math.round(rand(7, 14));
    const M = 6; // keep every mark this far inside the spine (its corners are rounded)
    const inside = (x, y) => x > M && x < w - M && y > M * 0.6 && y < h - M * 0.6;
    let d = '';
    // One straight nick: a thin sliver, widest in the middle and sharp at both ends,
    // pale on top with a dark groove just below. (Even-width strokes look like hairs.)
    const nick = (x0, y0, x1, y1, hw, a) => {
      const dx = x1 - x0, dy = y1 - y0, len = Math.hypot(dx, dy) || 1;
      const nx = -dy / len * hw, ny = dx / len * hw;
      const mx = (x0 + x1) / 2, my = (y0 + y1) / 2;
      const f = (v) => v.toFixed(1);
      const path = `M${f(x0)} ${f(y0)}Q${f(mx + nx)} ${f(my + ny)} ${f(x1)} ${f(y1)}Q${f(mx - nx)} ${f(my - ny)} ${f(x0)} ${f(y0)}Z`;
      d += `<path d='${path}' fill='rgba(0,0,0,${(a * 1.2).toFixed(2)})' transform='translate(.4 .7)'/>`;
      d += `<path d='${path}' fill='rgba(232,222,200,${a.toFixed(2)})'/>`;
    };
    // A scratch: one to three nicks on one straight line with small gaps between them
    const scratch = (x0, y0, ang, len, hw, a) => {
      const ux = Math.cos(ang), uy = Math.sin(ang);
      if (!inside(x0, y0) || !inside(x0 + ux * len, y0 + uy * len)) return false; // would leave the spine
      const parts = len > 14 && Math.random() < 0.6 ? Math.round(rand(2, 3)) : 1;
      let t = 0;
      for (let p = 0; p < parts; p++) {
        const seg = (len - t) / (parts - p) * rand(0.6, 0.95);
        nick(x0 + ux * t, y0 + uy * t, x0 + ux * (t + seg), y0 + uy * (t + seg), hw * rand(0.7, 1.1), a * rand(0.7, 1));
        t += seg + rand(1.5, 4);
        if (t >= len) break;
      }
      return true;
    };
    for (let i = 0, made = 0; i < n * 4 && made < n; i++) {
      // mostly short; now and then a longer one
      const len = Math.random() < 0.8 ? rand(4, 18) : rand(18, Math.max(20, Math.min(w * 0.3, 44)));
      const along = Math.random() < 0.72; // books slide against each other: most run along the spine
      const ang = (along ? rand(-12, 12) : rand(-70, 70)) * Math.PI / 180;
      if (scratch(rand(M, w - M), rand(M * 0.6, h - M * 0.6), ang, len, rand(0.35, 0.8), rand(0.14, 0.32))) made++;
    }
    if (kind === 'battered') {
      for (let c = 0; c < Math.round(rand(1, 3)); c++) {
        // scuffed spots near the ends, where a book gets knocked: bunches of tiny nicks
        const cx = Math.random() < 0.5 ? rand(M + 6, w * 0.2) : rand(w * 0.8, w - M - 6);
        const cy = rand(M + 2, h - M - 2);
        const ang = rand(-30, 30) * Math.PI / 180;
        for (let j = 0; j < Math.round(rand(6, 12)); j++) {
          scratch(cx + rand(-6, 6), cy + rand(-4, 4), ang + rand(-0.2, 0.2), rand(2, 6), rand(0.25, 0.45), rand(0.1, 0.22));
        }
      }
    }
    if (!d) return 'none';
    const svg = `<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 ${Math.round(w)} ${Math.round(h)}' preserveAspectRatio='none'>${d}</svg>`;
    return `url("data:image/svg+xml,${encodeURIComponent(svg)}")`;
  }

  // ---------- A cobweb under the desk, and its spider ----------
  //
  // The web hangs from the underside of the desk: threads from a hub up to the board,
  // rings between them, a couple of broken threads drooping. The spider hangs on its
  // own thread from the hub and sways slowly. Clicked, it raises a front leg and shakes
  // it at you like a fist.
  function addWeb(W, k) {
    const w = Math.round(150 * k), h = Math.round(86 * k);
    const x = Math.round(W * rand(0.58, 0.82) - w / 2);
    const hub = [w * rand(0.42, 0.58), h * rand(0.5, 0.62)];
    const f = (v) => v.toFixed(1);

    // spokes: to points along the board and a few down to the sides
    const ends = [];
    const nTop = 5;
    for (let i = 0; i < nTop; i++) ends.push([w * (0.04 + (0.92 * i) / (nTop - 1)) + rand(-4, 4), 0]);
    ends.push([0, h * rand(0.15, 0.35)], [w, h * rand(0.15, 0.35)]);
    ends.sort((a, b) => Math.atan2(a[1] - hub[1], a[0] - hub[0]) - Math.atan2(b[1] - hub[1], b[0] - hub[0]));
    let d = ends.map(([ex, ey]) => `M${f(hub[0])} ${f(hub[1])}L${f(ex)} ${f(ey)}`).join('');
    // rings: between neighbouring spokes, sagging a little towards the hub
    for (let r = 0.22; r < 0.95; r += rand(0.14, 0.2)) {
      for (let i = 0; i < ends.length - 1; i++) {
        if (Math.random() < 0.12) continue; // torn here
        const p = (e) => [hub[0] + (e[0] - hub[0]) * r, hub[1] + (e[1] - hub[1]) * r];
        const a = p(ends[i]), b = p(ends[i + 1]);
        const mx = (a[0] + b[0]) / 2 + (hub[0] - (a[0] + b[0]) / 2) * 0.12;
        const my = (a[1] + b[1]) / 2 + (hub[1] - (a[1] + b[1]) / 2) * 0.12 + 1.5;
        d += `M${f(a[0])} ${f(a[1])}Q${f(mx)} ${f(my)} ${f(b[0])} ${f(b[1])}`;
      }
    }
    // broken threads drooping from the web
    for (let i = 0; i < 2; i++) {
      const sx = hub[0] + rand(-0.35, 0.35) * w, sy = hub[1] * rand(0.6, 1);
      d += `M${f(sx)} ${f(sy)}q${f(rand(-6, 6))} ${f(h * 0.25)} ${f(rand(-10, 10))} ${f(h * rand(0.35, 0.6))}`;
    }

    const web = document.createElement('div');
    web.className = 'web';
    web.style.left = `${x}px`;
    web.style.width = `${w}px`;
    web.style.height = `${h}px`;
    web.innerHTML = `<svg viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" aria-hidden="true"><path d="${d}"/></svg>`;
    host.appendChild(web);

    // the spider on its thread, hanging from the hub
    const len = Math.round(rand(60, 90) * k);
    const size = Math.round(46 * Math.max(0.85, k));
    const hang = document.createElement('div');
    hang.className = 'spider-hang';
    hang.style.left = `${Math.round(x + hub[0])}px`;
    hang.style.top = `calc(100% + ${Math.round(hub[1])}px)`;
    hang.style.setProperty('--len', `${len}px`);
    hang.style.setProperty('--size', `${size}px`);
    hang.style.animationDuration = `${rand(3.2, 4.4).toFixed(1)}s`;
    hang.innerHTML = `
      <span class="spider-hang__thread"></span>
      <button class="spider" type="button" aria-label="${escapeHtml(t('texts.spider', 'Паук'))}">
        <svg viewBox="-20 -20 40 40" width="${size}" height="${size}" aria-hidden="true">
          <g class="spider__legs">
            <path d="M-3 -1 Q-10 -9 -15 -4"/><path d="M-3 1 Q-12 -2 -17 4"/>
            <path d="M-3 3 Q-11 6 -15 13"/><path d="M-2 5 Q-7 11 -9 17"/>
            <path d="M3 1 Q12 -2 17 4"/>
            <path d="M3 3 Q11 6 15 13"/><path d="M2 5 Q7 11 9 17"/>
          </g>
          <g class="spider__fist"><path d="M3 -1 Q10 -9 15 -4"/></g>
          <g class="spider__arm">
            <path d="M3 -1 Q8 .5 11.5 -1.5"/>
            <g class="spider__forearm"><path d="M11.5 -1.5 Q12.6 -6 12.2 -11"/><g transform="translate(12.2 -13.6)"><g class="spider__knuckles"><circle r="3.1"/><path d="M-2 -1.1 q1 -1.1 2 0 M-.4 .5 q1 -1.1 2 0"/></g></g></g>
          </g>
          <ellipse class="spider__body" cx="0" cy="6" rx="6.5" ry="7.5"/>
          <circle class="spider__body" cx="0" cy="-3.5" r="4"/>
          <circle class="spider__eye" cx="-1.5" cy="-4.5" r=".9"/><circle class="spider__eye" cx="1.5" cy="-4.5" r=".9"/>
          <g class="spider__brows"><path d="M-3.6 -7.2 L-.7 -5.8"/><path d="M3.6 -7.2 L.7 -5.8"/></g>
        </svg>
      </button>`;
    host.appendChild(hang);
    // room under the desk for the web and the dangling spider
    host.style.marginBottom = `${Math.round(hub[1] + len + size + 24 * k)}px`;

    const spider = hang.querySelector('.spider');
    spider.addEventListener('click', () => {
      spider.classList.remove('is-angry');
      void spider.offsetWidth; // restart the animation on every click
      spider.classList.add('is-angry');
      grumble(spider);
    });
    spider.addEventListener('animationend', (e) => {
      if (e.animationName === 'spider-arm') spider.classList.remove('is-angry');
    });
  }

  // The spider's displeasure, comic-strip style: muttered curses (#@!), a throbbing
  // anger mark and a couple of puffs pop out around it and drift off.
  const ANGER = '<svg viewBox="-10 -10 20 20"><path d="M-2 -8 Q-2 -2 -8 -2 M2 -8 Q2 -2 8 -2 M-2 8 Q-2 2 -8 2 M2 8 Q2 2 8 2"/></svg>';
  function grumble(spider) {
    const k = Math.max(0.8, scale);
    const hostBox = host.getBoundingClientRect();
    const r = spider.getBoundingClientRect();
    const cx = r.left + r.width / 2 - hostBox.left, cy = r.top + r.height / 2 - hostBox.top;
    const words = ['#@!', '%&!', '@#$!', '!!', '*#%', '?!'];
    const marks = [];
    marks.push({ html: ANGER, cls: 'grumble grumble--anger', ang: rand(-150, -110), dist: 26 });
    for (let i = 0; i < 3; i++) {
      marks.push({ text: words[Math.floor(rand(0, words.length))], cls: 'grumble', ang: rand(-170, -10), dist: rand(28, 40) });
    }
    marks.push({ html: '', cls: 'grumble grumble--puff', ang: rand(-60, -20), dist: 22 });
    marks.forEach((m, i) => {
      const el = document.createElement('span');
      el.className = m.cls;
      if (m.html) el.innerHTML = m.html; else el.textContent = m.text;
      el.style.left = `${cx.toFixed(1)}px`;
      el.style.top = `${cy.toFixed(1)}px`;
      el.style.fontSize = `${Math.round(11 * k)}px`;
      host.appendChild(el);
      const a = (m.ang * Math.PI) / 180, d = m.dist * k;
      const x = Math.cos(a) * d, y = Math.sin(a) * d;
      const tilt = rand(-18, 18);
      el.animate([
        { transform: `translate(-50%,-50%) translate(${(x * 0.4).toFixed(1)}px, ${(y * 0.4).toFixed(1)}px) scale(.2) rotate(${tilt}deg)`, opacity: 0 },
        { offset: 0.18, transform: `translate(-50%,-50%) translate(${x.toFixed(1)}px, ${y.toFixed(1)}px) scale(1.15) rotate(${tilt}deg)`, opacity: 1 },
        { offset: 0.3, transform: `translate(-50%,-50%) translate(${x.toFixed(1)}px, ${y.toFixed(1)}px) scale(.95) rotate(${tilt}deg)`, opacity: 1 },
        { offset: 0.42, transform: `translate(-50%,-50%) translate(${x.toFixed(1)}px, ${y.toFixed(1)}px) scale(1.08) rotate(${tilt}deg)`, opacity: 1 },
        { transform: `translate(-50%,-50%) translate(${(x * 1.3).toFixed(1)}px, ${(y * 1.3 - 14 * k).toFixed(1)}px) scale(1) rotate(${tilt}deg)`, opacity: 0 },
      ], { duration: 1300, delay: i * 110, easing: 'ease-out', fill: 'both' }).finished.then(() => el.remove(), () => {});
    });
  }

  // ---------- Dust ----------
  //
  // Old books are dusty: when one moves, a soft cloud puffs up around it and specks
  // of dust drift slowly down, spreading, and fade. Specks are small radial gradients
  // (no blur filter, so they stay cheap); only transform and opacity are animated.
  function dust(from, count, strength) {
    const k = scale;
    const hostBox = host.getBoundingClientRect();
    const r = from.getBoundingClientRect();
    const frag = document.createDocumentFragment();
    const jobs = [];

    // a cloud or two: grows and thins out
    for (let i = 0; i < 2; i++) {
      const c = document.createElement('span');
      c.className = 'dust dust--cloud';
      const size = (r.width * rand(0.45, 0.7)) * (0.7 + 0.3 * strength);
      c.style.width = `${size.toFixed(0)}px`;
      c.style.height = `${(size * 0.55).toFixed(0)}px`;
      c.style.left = `${(r.left - hostBox.left + r.width * rand(0.15, 0.85) - size / 2).toFixed(1)}px`;
      c.style.top = `${(r.top - hostBox.top + r.height * 0.5 - size * 0.27).toFixed(1)}px`;
      frag.appendChild(c);
      jobs.push([c, [
        { transform: 'scale(.4)', opacity: 0 },
        { offset: 0.15, opacity: 0.55 * strength },
        { transform: `translate(${rand(-20, 20) * k}px, ${rand(6, 22) * k}px) scale(1.7)`, opacity: 0 },
      ], { duration: rand(1300, 1900), easing: 'cubic-bezier(.2,.6,.3,1)', fill: 'both' }]);
    }

    // specks: drift out a little, then sink slowly with a lazy sideways sway
    for (let i = 0; i < count; i++) {
      const s = document.createElement('span');
      s.className = 'dust';
      const size = rand(2, 5) * Math.max(0.8, k);
      s.style.width = s.style.height = `${size.toFixed(1)}px`;
      s.style.left = `${(r.left - hostBox.left + rand(0, r.width)).toFixed(1)}px`;
      s.style.top = `${(r.top - hostBox.top + rand(0, r.height)).toFixed(1)}px`;
      frag.appendChild(s);
      const dx = rand(-40, 40) * k, up = rand(-18, 6) * k, down = rand(40, 120) * k, sway = rand(6, 16) * k;
      jobs.push([s, [
        { transform: 'translate(0,0)', opacity: 0 },
        { offset: 0.1, transform: `translate(${(dx * 0.3).toFixed(1)}px, ${up.toFixed(1)}px)`, opacity: rand(0.45, 0.85) * strength },
        { offset: 0.5, transform: `translate(${(dx * 0.7 + sway).toFixed(1)}px, ${(up + down * 0.45).toFixed(1)}px)` },
        { transform: `translate(${(dx - sway).toFixed(1)}px, ${(up + down).toFixed(1)}px)`, opacity: 0 },
      ], { duration: rand(1500, 2600), delay: rand(0, 160), easing: 'ease-out', fill: 'both' }]);
    }

    host.appendChild(frag);
    for (const [el, kf, opt] of jobs) el.animate(kf, opt).finished.then(() => el.remove(), () => {});
  }

  // ---------- What falls out of a book ----------
  //
  // Every section has its own folder of sprites: SPRITE_DIR/<section id>/1.png, 2.png, ...
  // numbered without gaps. They are found by trying 1.png, 2.png, ... until one is missing,
  // so new sprites only need to be dropped into the folder. Yellowed paper slips are drawn
  // by CSS and fall out of every book.
  const sprites = new Map(); // section id -> [urls]
  function findSprites(id) {
    if (sprites.has(id)) return;
    const list = [];
    sprites.set(id, list);
    const tryNext = (n) => {
      if (n > 40) return;
      const im = new Image();
      im.onload = () => { list.push(im.src); tryNext(n + 1); };
      im.src = `${SPRITE_DIR}${encodeURIComponent(id)}/${n}.png`;
    };
    tryNext(1);
  }
  // start looking right away, so the pictures are ready by the first click
  SECTIONS.forEach((s) => findSprites(s.id));

  function spill(bookEl, from) {
    const k = scale;
    const hostBox = host.getBoundingClientRect();
    const box = (from || bookEl).getBoundingClientRect(); // the falling 3D book, if there is one
    const cx = box.left + box.width / 2 - hostBox.left;
    const cy = box.top + box.height / 2 - hostBox.top;
    const list = sprites.get(bookEl.dataset.section) || [];
    const frag = document.createDocumentFragment();
    const jobs = [];

    const pieces = [];
    const nSprites = list.length ? Math.round(rand(8, 13)) : 0;
    for (let i = 0; i < nSprites; i++) pieces.push({ src: list[Math.floor(Math.random() * list.length)] });
    // pages torn out of this very book: the spine length is the height of its pages
    const pageH = bookEl.offsetWidth * 0.9;
    const nSlips = Math.round(rand(1, 2)); // just a page or two: big flickering pages tire the eyes
    for (let i = 0; i < nSlips; i++) pieces.push({ slip: true });

    for (const p of pieces) {
      let el;
      if (p.slip) {
        el = document.createElement('span');
        el.className = 'spill spill--slip';
        const h = pageH * rand(0.96, 1.02), w = h * rand(0.64, 0.7);
        el.style.width = `${Math.round(w)}px`;
        el.style.height = `${Math.round(h)}px`;
        el.style.padding = `${Math.round(h * 0.08)}px ${Math.round(w * 0.11)}px ${Math.round(h * 0.1)}px`; // margins
        el.style.setProperty('--line', `${Math.max(4, h / 34).toFixed(1)}px`); // line spacing grows with the page
        el.style.setProperty('--tone', rand(0.75, 1).toFixed(2));
      } else {
        el = document.createElement('img');
        el.className = 'spill';
        el.src = p.src;
        el.alt = '';
        el.style.width = `${Math.round(rand(30, 54) * k)}px`;
      }
      // centred on a point of the book (pages are as big as the book, so this matters)
      const half = parseFloat(el.style.width) / 2;
      const halfH = p.slip ? parseFloat(el.style.height) / 2 : half;
      el.style.left = `${(cx + rand(-0.25, 0.25) * box.width - half).toFixed(1)}px`;
      el.style.top = `${(cy + rand(-8, 8) - halfH).toFixed(1)}px`;
      frag.appendChild(el);

      // burst out of the pages, then drift down swaying like a leaf
      // out of a book that is already falling: a short toss sideways and only a little up
      const bx = rand(-150, 150) * k, by = -rand(0, 50) * k;
      const fallY = rand(260, 420) * k;
      const sway = (p.slip ? rand(30, 60) : rand(14, 34)) * k * (Math.random() < 0.5 ? -1 : 1);
      // pages glide calmly, without turning over (turning big pages flickers); leaves turn slowly
      const r0 = rand(-40, 40);
      const spin = (p.slip ? rand(15, 45) : rand(90, 220)) * (Math.random() < 0.5 ? -1 : 1);
      const flipDeg = p.slip ? 0 : rand(120, 200);
      const flip = Math.random() < 0.5 ? 'rotateX' : 'rotateY';
      const at = (t, x, y, extra = {}) => ({
        offset: t,
        transform: `translate(${x.toFixed(1)}px, ${y.toFixed(1)}px) rotate(${(r0 + spin * t).toFixed(0)}deg) ${flip}(${(t * flipDeg).toFixed(0)}deg)`,
        ...extra,
      });
      jobs.push([el, [
        at(0, 0, 0, { opacity: 0 }),
        at(0.06, bx * 0.25, by * 0.4, { opacity: 1, easing: 'cubic-bezier(.2,.8,.4,1)' }),
        at(0.2, bx, by, { easing: 'ease-in-out' }),
        at(0.42, bx + sway, by + fallY * 0.3, { easing: 'ease-in-out' }),
        at(0.64, bx - sway * 0.8, by + fallY * 0.6, { easing: 'ease-in-out' }),
        at(0.84, bx + sway * 0.6, by + fallY * 0.85, { opacity: 1, easing: 'ease-in' }),
        at(1, bx, by + fallY, { opacity: 0 }),
      ], { duration: p.slip ? rand(2000, 2600) : rand(1400, 2000), delay: rand(60, 220), fill: 'both' }]);
    }

    host.appendChild(frag);
    for (const [el, kf, opt] of jobs) {
      el.animate(kf, opt).finished.then(() => el.remove(), () => {});
    }
  }

  // The falling book is a 3D box, so as it tumbles we see its cover, the block of
  // pages and the ends, not only the spine. Size: length = spine length, thickness =
  // spine height, depth (the width of the pages) about 2/3 of the length.
  //   front  – the spine (a copy of the book element itself)
  //   back   – the fore-edge, the block of pages
  //   top / bottom – the covers
  //   left / right – head and tail: pages between the two boards
  function bookBox(el) {
    const L = el.offsetWidth, T = el.offsetHeight, D = Math.round(L * 0.66);
    // Two layers: the wrap gives perspective and fades out; the box inside turns in 3D.
    // (Opacity on the 3D box itself would flatten it into a plane.)
    const wrap = document.createElement('div');
    wrap.className = 'book3d-wrap';
    wrap.style.left = el.style.left;
    wrap.style.bottom = el.style.bottom;
    wrap.style.width = `${L}px`;
    wrap.style.height = `${T}px`;
    const box = document.createElement('div');
    box.className = 'book3d';
    wrap.appendChild(box);
    for (const v of ['--c', '--bx', '--by', '--sat', '--br', '--wear', '--wx', '--wy', '--scratch']) box.style.setProperty(v, el.style.getPropertyValue(v));

    const spine = el.cloneNode(true);
    spine.className = el.className.replace(/\bis-falling\b/, '') + ' book3d__face';
    spine.removeAttribute('href');
    ['stack', 'level', 'drop', 'section', 'z'].forEach((k) => delete spine.dataset[k]); // not a book of the heap
    spine.querySelector('.book__top')?.remove(); // the 3D box has a real cover
    spine.removeAttribute('style');
    spine.style.cssText = el.style.cssText;
    Object.assign(spine.style, { left: 0, top: 0, bottom: 'auto', width: `${L}px`, height: `${T}px`, zIndex: 'auto', transition: 'none' });
    spine.style.transform = `translateZ(${D / 2}px)`;

    const face = (cls, w, h, transform) => {
      const f = document.createElement('div');
      f.className = `book3d__face ${cls}`;
      Object.assign(f.style, {
        width: `${w}px`, height: `${h}px`,
        left: `${(L - w) / 2}px`, top: `${(T - h) / 2}px`,
        transform,
      });
      return f;
    };
    box.append(
      spine,
      face('book3d__pages', L, T, `rotateY(180deg) translateZ(${D / 2}px)`),
      face('book3d__cover', L, D, `rotateX(90deg) translateZ(${T / 2}px)`),
      face('book3d__cover', L, D, `rotateX(-90deg) translateZ(${T / 2}px)`),
      face('book3d__pages book3d__pages--end', D, T, `rotateY(-90deg) translateZ(${L / 2}px)`),
      face('book3d__pages book3d__pages--end', D, T, `rotateY(90deg) translateZ(${L / 2}px)`),
    );
    return { wrap, box };
  }

  // The clicked book slides out of its stack, tips over the front edge of the desk
  // and falls, tumbling, so its cover and pages show. Only transform and opacity are animated.
  function fall(el) {
    el.classList.add('is-falling');
    el.dataset.z = el.style.zIndex;
    el.style.transition = 'none';
    const r = parseFloat(el.style.getPropertyValue('--r')) || 0;
    const dir = Math.random() < 0.5 ? -1 : 1;
    const drop = parseFloat(el.style.bottom) + 340; // past the desk and further down

    const { wrap, box } = bookBox(el);
    host.appendChild(wrap);
    el.style.visibility = 'hidden'; // the 3D copy takes its place

    // tips towards the viewer around its long axis (rotateX), turns a bit on the way
    // (rotateY) and keeps its tilt in the picture plane (rotateZ)
    const tipX = -rand(150, 230);
    const turnY = dir * rand(15, 45);
    const tiltZ = r + dir * rand(15, 40);
    const at = (x, y, rx, ry, rz) =>
      `translate3d(${x.toFixed(1)}px, ${y.toFixed(1)}px, 0) rotateZ(${rz.toFixed(1)}deg) rotateY(${ry.toFixed(1)}deg) rotateX(${rx.toFixed(1)}deg)`;
    const dur = FALL_MS + 150;
    box.animate([
      { transform: at(22, 0, 0, 0, r), easing: 'cubic-bezier(.2,.7,.4,1)' },
      // pulled out and lifted a little, the front edge starts to tip
      { offset: 0.18, transform: at(22 + dir * 14, -14, -18, turnY * 0.2, r - dir * 3), easing: 'cubic-bezier(.45,0,.85,.55)' },
      { transform: at(22 + dir * 70, drop, tipX, turnY, tiltZ) },
    ], { duration: dur, fill: 'forwards' });
    wrap.animate([{ opacity: 1 }, { offset: 0.8, opacity: 1 }, { opacity: 0 }], { duration: dur, fill: 'forwards' });
    // pages and leaves fall out once the book is already falling, from where it is then
    setTimeout(() => spill(el, box), SPILL_AT_MS);
    // dust: a puff where it slid out of the stack, then a trail behind it as it falls
    dust(el, 26, 1);
    setTimeout(() => dust(box, 14, 0.6), SPILL_AT_MS - 80);

    // the books that lay on it drop into the gap, with a small bump at the end
    const stack = el.dataset.stack, level = +el.dataset.level, gap = +el.dataset.drop;
    if (stack === undefined || !gap) return;
    host.querySelectorAll(`.book[data-stack="${stack}"]`).forEach((b) => {
      if (+b.dataset.level <= level) return;
      const br = parseFloat(b.style.getPropertyValue('--r')) || 0;
      const wob = (Math.random() - 0.5) * 3;
      b.animate([
        { transform: `rotate(${br}deg)` },
        { offset: 0.75, transform: `translateY(${gap}px) rotate(${br + wob}deg)`, easing: 'ease-out' },
        { offset: 0.88, transform: `translateY(${gap - 3}px) rotate(${br + wob * 0.6}deg)`, easing: 'ease-in' },
        { transform: `translateY(${gap}px) rotate(${br + wob * 0.5}deg)` },
      ], { duration: 360, delay: 110 + (+b.dataset.level - level) * 25, easing: 'cubic-bezier(.5,0,.9,.6)', fill: 'forwards' });
    });
  }

  // A section page: for now just its name and a note that it is empty
  function showSection(s) {
    document.title = `${s.label} · Thorntales`;
    const h1 = document.querySelector('.top h1');
    const sub = document.querySelector('.top .sub');
    if (h1) h1.textContent = s.label;
    if (sub) sub.textContent = t('texts.empty', 'Эта книга пока пустая.');
    if (window.Layout) Layout.setCrumbs([{ label: s.label }]);
    host.outerHTML = `<p class="texts-back"><a href="/texts.html">${escapeHtml(t('texts.back', '← к стопкам книг'))}</a></p>`;
  }

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }

  // after a click the page may come back from the back/forward cache: put the book back
  window.addEventListener('pageshow', (e) => {
    if (!e.persisted) return;
    host.querySelectorAll('.book3d-wrap, .dust, .spill').forEach((b) => b.remove());
    host.querySelectorAll('.book').forEach((b) => b.getAnimations().forEach((a) => a.cancel()));
    host.querySelectorAll('.is-falling').forEach((b) => {
      b.getAnimations().forEach((a) => a.cancel());
      b.classList.remove('is-falling');
      b.style.transition = '';
      b.style.zIndex = b.dataset.z;
      b.style.visibility = '';
    });
  });
})();
