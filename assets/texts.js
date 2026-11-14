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
  const TAKE_DELAY_MS = 1150; // the section opens while the leaves are still in the air
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

    stacks.forEach((stack, si) => {
      const cx = slot * (si + 0.5) + rand(-0.12, 0.12) * slot;
      let y = deskH; // from the bottom of the desk area
      let widest = 0;
      stack.forEach((b, level) => {
        const minLen = b.section ? textWidth(b.section.label, titlePx) + 64 * k : 0;
        const len = Math.min(slot * 0.92, Math.max(minLen, rand(170, 280) * k));
        const th = Math.round((b.section ? rand(44, 62) : rand(34, 64)) * k); // thick old volumes
        // now and then a book sticks far out of its stack
        const shift = (Math.random() < 0.18 ? rand(30, 48) * (Math.random() < 0.5 ? -1 : 1) : rand(-20, 20)) * k;
        const rot = rand(-4, 4) * (level === 0 ? 0.35 : 1);
        els.push(bookEl(b, {
          left: Math.max(2, Math.min(W - len - 2, cx - len / 2 + shift)), bottom: y, w: len, h: th, rot, z: 10 + level,
          color: pick(COLORS), titlePx, stack: si, level, drop: th - 1,
        }));
        y += th - 1; // books sit tight, a slight overlap hides gaps from rotation
        widest = Math.max(widest, len / 2 + Math.abs(shift));
      });
      maxTop = Math.max(maxTop, y);
      stack.cx = cx;
      stack.half = widest;
      stack.top = y;
    });

    // a nameless book thrown across the tops of two neighbouring stacks
    if (stacksN >= 2) {
      const i = Math.floor(rand(0, stacksN - 1));
      const a = stacks[i], b = stacks[i + 1];
      const x1 = a.cx + rand(-0.15, 0.15) * a.half, x2 = b.cx + rand(-0.15, 0.15) * b.half;
      const dx = x2 - x1, dy = b.top - a.top;
      const len = Math.hypot(dx, dy) + rand(40, 90) * k;
      const th = Math.round(rand(36, 50) * k);
      const ang = -Math.atan2(dy, dx) * 180 / Math.PI; // CSS rotation, y up
      els.push(bookEl({ section: null }, {
        left: Math.max(2, Math.min(W - len - 2, (x1 + x2) / 2 - len / 2)), bottom: (a.top + b.top) / 2 - 2, w: len, h: th,
        rot: ang + rand(-1.5, 1.5), z: 25, color: pick(COLORS), titlePx,
      }));
      maxTop = Math.max(maxTop, Math.max(a.top, b.top) + th);
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
        left: pivot, bottom: deskH, w: len, h: th, color: pick(COLORS), titlePx, z: 5, // behind the stack: never covers a title
        upright: true, rot: side > 0 ? -90 - lean : -90 + lean,
      }));
      maxTop = Math.max(maxTop, deskH + len);
    });

    const desk = document.createElement('div');
    desk.className = 'desk__top';
    desk.style.height = `${deskH}px`;
    host.appendChild(desk);
    els.forEach((el) => host.appendChild(el));
    host.style.height = `${Math.ceil(maxTop + 40 * k)}px`;
  }

  function bookEl(b, o) {
    const el = document.createElement(b.section ? 'a' : 'span');
    el.className = 'book' + (b.section ? ' book--section' : '') + (o.upright ? ' book--upright' : '');
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
    return el;
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

  function spill(bookEl) {
    const k = scale;
    const hostBox = host.getBoundingClientRect();
    const box = bookEl.getBoundingClientRect();
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
    const nSlips = Math.round(rand(4, 7));
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
      const bx = rand(-150, 150) * k, by = -rand(30, 120) * k;
      const fallY = rand(260, 420) * k;
      const sway = rand(14, 34) * k * (Math.random() < 0.5 ? -1 : 1);
      const r0 = rand(-40, 40), spin = rand(120, 320) * (Math.random() < 0.5 ? -1 : 1);
      const flip = Math.random() < 0.5 ? 'rotateX' : 'rotateY'; // paper and leaves turn over as they fall
      const at = (t, x, y, extra = {}) => ({
        offset: t,
        transform: `translate(${x.toFixed(1)}px, ${y.toFixed(1)}px) rotate(${(r0 + spin * t).toFixed(0)}deg) ${flip}(${(t * 540).toFixed(0)}deg)`,
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
      ], { duration: rand(1300, 1900), delay: rand(60, 220), fill: 'both' }]);
    }

    host.appendChild(frag);
    for (const [el, kf, opt] of jobs) {
      el.animate(kf, opt).finished.then(() => el.remove(), () => {});
    }
  }

  // The clicked book slides out of its stack, tips over the edge and falls off the
  // desk, tumbling, in front of everything. Only transform and opacity are animated.
  function fall(el) {
    el.classList.add('is-falling');
    el.dataset.z = el.style.zIndex;
    el.style.zIndex = 100;
    el.style.transition = 'none';
    const r = parseFloat(el.style.getPropertyValue('--r')) || 0;
    const dir = Math.random() < 0.5 ? -1 : 1;
    const drop = parseFloat(el.style.bottom) + 320; // past the desk and further down
    const out = 'translateX(22px)';                  // where hover left it
    el.animate([
      { transform: `rotate(${r}deg) ${out}`, opacity: 1, easing: 'cubic-bezier(.2,.7,.4,1)' },
      { offset: 0.2, transform: `translate(${dir * 16}px, -12px) rotate(${r - dir * 5}deg) ${out}`, opacity: 1, easing: 'cubic-bezier(.5,0,.9,.5)' },
      { offset: 0.82, opacity: 1 },
      { transform: `translate(${dir * 80}px, ${drop}px) rotate(${r + dir * 62}deg) ${out}`, opacity: 0 },
    ], { duration: FALL_MS, fill: 'forwards' });
    spill(el);

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
    host.querySelectorAll('.book').forEach((b) => b.getAnimations().forEach((a) => a.cancel()));
    host.querySelectorAll('.is-falling').forEach((b) => {
      b.getAnimations().forEach((a) => a.cancel());
      b.classList.remove('is-falling');
      b.style.transition = '';
      b.style.zIndex = b.dataset.z;
    });
  });
})();
