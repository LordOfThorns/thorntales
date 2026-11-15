/* search.js — search.html: find poems and texts by words and by tags.

   Sources (everything the reader may see, nothing else):
     poems   data/poetry/catalogs.json + <cat>.txt, with site edits (poem_edits) applied;
             poems hidden from this reader (PoetryData.hiddenPoems) are dropped, and
             poem_edits rows the reader may not read never arrive (RLS)
     texts   data/texts/<section>.json + text_edits (hidden texts are left out)
     polar   data/polar-night/texts.json (same format as a text section; may not exist yet)
   Tags come from the poems the reader sees, so a restricted tag shows only to those
   who may see its poems.

   URL: search.html?q=<words>&tag=<tag>[,<tag>]&in=all|poems|texts|polar
*/
// SiteFetch.get comes from layout.js; right after a deploy the browser may still run an
// older cached layout.js without it, so fall back to a plain fetch with the same answers.
window.siteGet = window.siteGet || ((url, o = {}) => (window.SiteFetch
  ? SiteFetch.get(url, o)
  : fetch(url, { cache: 'no-store' }).then((r) => {
    if (r.ok) return r;
    if (o.optional || (o.missingOk && r.status === 404)) return null;
    throw new Error(`${url}: ${r.status}`);
  }, (e) => { if (o.optional) return null; throw e; })));

(async function () {
  'use strict';

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const form = document.getElementById('searchForm');
  const input = document.getElementById('searchInput');
  const scopeSel = document.getElementById('searchScope');
  const tagsHost = document.getElementById('searchTags');
  const statusEl = document.getElementById('searchStatus');
  const results = document.getElementById('searchResults');
  const moreBtn = document.getElementById('searchMore');

  const SECTIONS = ['stories', 'tales', 'repatriate', 'dreams', 'chatter', 'other'];
  const PAGE = 40;
  const ENDINGS = ['ями', 'ами', 'ого', 'его', 'ому', 'ему', 'ыми', 'ими', 'ах', 'ях', 'ов', 'ев', 'ой', 'ей',
    'ий', 'ый', 'ая', 'яя', 'ое', 'ее', 'ом', 'ем', 'ам', 'ям', 'ую', 'юю', 'ть', 'а', 'я', 'ы', 'и', 'у', 'ю', 'е', 'о', 'ь'];
  const db = window.Auth ? Auth.client : null;

  // State from the URL
  const params = new URLSearchParams(location.search);
  input.value = params.get('q') || '';
  scopeSel.value = ['all', 'poems', 'texts', 'polar'].includes(params.get('in')) ? params.get('in') : 'all';
  const activeTags = new Set((params.get('tag') || '').split(',').map((s) => s.trim()).filter(Boolean));
  let shown = PAGE;

  let items = [];
  try {
    if (window.Auth) await Auth.ready; // what the reader sees depends on who is signed in
    items = await loadAll();
  } catch (e) {
    console.error(e);
    statusEl.textContent = t('poems.loadFailed') + ': ' + (e.message || e);
    return;
  }

  renderTags();
  run();

  let timer = 0;
  input.addEventListener('input', () => { clearTimeout(timer); timer = setTimeout(() => { shown = PAGE; run(); }, 180); });
  scopeSel.addEventListener('change', () => { shown = PAGE; renderTags(); run(); });
  form.addEventListener('submit', (e) => { e.preventDefault(); shown = PAGE; run(); });
  tagsHost.addEventListener('click', (e) => {
    const b = e.target.closest('[data-tag]');
    if (!b) return;
    const tag = b.dataset.tag;
    if (activeTags.has(tag)) activeTags.delete(tag); else activeTags.add(tag);
    shown = PAGE;
    renderTags();
    run();
  });
  moreBtn.addEventListener('click', () => { shown += PAGE; run(); });

  // ---------- Loading ----------

  async function loadAll() {
    const [poems, texts, polar] = await Promise.all([loadPoems(), loadTexts(), loadPolar()]);
    return [...poems, ...texts, ...polar];
  }

  async function loadPoems() {
    const [cats, hidden, rows] = await Promise.all([
      PoetryData.catalogs(),
      PoetryData.hiddenPoems(),
      selectAll('poem_edits'),
    ]);
    const editsByCat = new Map();
    for (const r of rows) {
      if (!editsByCat.has(r.catalog)) editsByCat.set(r.catalog, new Map());
      editsByCat.get(r.catalog).set(r.poem_id, r);
    }
    const out = [];
    await Promise.all(cats.map(async (c) => {
      const raw = await fetchText(`./data/poetry/${encodeURIComponent(c.id)}.txt`);
      const hiddenIds = hidden.get(c.id) || new Set();
      const edits = editsByCat.get(c.id) || new Map();
      const seen = new Set();
      const label = (window.I18n ? I18n.field(c, 'label') : c.label) || c.id;
      const push = (id, title, date, text, tags, tr) => out.push({
        kind: 'poems', id, title: title || '', date: date || '', tags: tags || [], where: label,
        href: `/poems.html?cat=${encodeURIComponent(c.id)}#${encodeURIComponent(id)}`,
        text, extra: tr || '',
      });
      for (const { meta, text } of PoetryData.parseTxt(raw || '')) {
        const id = meta.id || '';
        if (!id || hiddenIds.has(id)) continue;
        seen.add(id);
        const e = edits.get(id);
        if (e) push(id, e.title, e.date, plain(e.body), e.tags, trText(e.tr));
        else push(id, meta.title, meta.date, meta.format === 'html' ? plain(text) : text, splitTags(meta.tags));
      }
      for (const e of edits.values()) {
        if (!seen.has(e.poem_id) && !hiddenIds.has(e.poem_id)) push(e.poem_id, e.title, e.date, plain(e.body), e.tags, trText(e.tr));
      }
    }));
    return out;
  }

  async function loadTexts() {
    const rows = await selectAll('text_edits');
    const out = [];
    await Promise.all(SECTIONS.map(async (s) => {
      const file = await fetchJson(`./data/texts/${s}.json`);
      const byId = new Map((file?.texts || []).map((x) => [x.id, {
        id: x.id, title: x.title, date: x.date, html: x.html,
        tr: Object.values(x.i18n || {}).map((v) => `${v.title || ''} ${plain(v.html)}`).join(' '),
      }]));
      for (const r of rows.filter((r) => r.section === s)) {
        if (r.deleted) { byId.delete(r.text_id); continue; }
        byId.set(r.text_id, { id: r.text_id, title: r.title, date: r.date, html: r.body, tr: trText(r.tr) });
      }
      const where = t('texts.' + s);
      for (const x of byId.values()) {
        out.push({
          kind: 'texts', id: x.id, title: x.title || '', date: x.date || '', tags: [], where,
          href: `/texts.html?s=${encodeURIComponent(s)}&t=${encodeURIComponent(x.id)}`,
          text: plain(x.html), extra: x.tr || '',
        });
      }
    }));
    return out;
  }

  async function loadPolar() {
    const file = await fetchJson('./data/polar-night/texts.json');
    return (file?.texts || []).map((x) => ({
      kind: 'polar', id: x.id, title: x.title || '', date: x.date || '', tags: x.tags || [],
      where: t('nav.polarNight'), href: x.href || '/polar-night.html',
      text: plain(x.html), extra: '',
    }));
  }

  // ---------- Searching ----------

  function inScope(it) {
    const s = scopeSel.value;
    return s === 'all' || it.kind === s;
  }

  function run() {
    const words = norm(input.value).split(' ').filter(Boolean).map(stem);
    syncUrl();
    results.innerHTML = '';
    moreBtn.hidden = true;

    const scoped = items.filter(inScope);
    if (scopeSel.value === 'polar' && !scoped.length) {
      statusEl.textContent = t('search.polarEmpty');
      return;
    }
    if (!words.length && !activeTags.size) {
      statusEl.textContent = t('search.hint');
      return;
    }

    const hits = [];
    for (const it of scoped) {
      if (activeTags.size && ![...activeTags].every((tag) => it.tags.includes(tag))) continue;
      let score = 0;
      if (words.length) {
        const title = norm(it.title);
        const hay = it._hay || (it._hay = ' ' + norm(`${it.title} ${it.text} ${it.extra} ${it.tags.join(' ')}`));
        if (!words.every((w) => hay.includes(' ' + w))) continue;
        score = words.reduce((s, w) => s + (title.includes(w) ? 10 : 0), 0);
        if (title.includes(words.join(' '))) score += 20;
      }
      hits.push({ it, score });
    }
    // best title matches first, then newest
    hits.sort((a, b) => b.score - a.score || (b.it.date || '').localeCompare(a.it.date || ''));

    statusEl.textContent = hits.length ? t('search.found', { n: hits.length }) : t('search.none');
    results.innerHTML = hits.slice(0, shown).map(({ it }) => renderHit(it, words)).join('');
    moreBtn.hidden = hits.length <= shown;
  }

  function renderHit(it, words) {
    const kind = { poems: t('search.kindPoem'), texts: t('search.kindText'), polar: t('nav.polarNight') }[it.kind];
    const title = it.title && !/^[\s*.…]+$/.test(it.title) ? it.title : firstLine(it.text);
    const tags = it.tags.length
      ? `<span class="search-hit__tags">${it.tags.map((tag) => `<span class="poem-tag${activeTags.has(tag) ? ' poem-tag--active' : ''}">#${esc(tag)}</span>`).join('')}</span>`
      : '';
    return `
      <li class="search-hit">
        <a class="search-hit__link" href="${esc(it.href)}">
          <span class="search-hit__title">${mark(title, words)}</span>
          <span class="search-hit__meta">${esc(kind)} · ${esc(it.where)}${it.date ? ' · ' + esc(it.date) : ''}</span>
          <span class="search-hit__snippet">${snippet(it.text, words)}</span>
          ${tags}
        </a>
      </li>`;
  }

  // A piece of the text around the first found word, words highlighted
  function snippet(text, words) {
    const flat = String(text || '').replace(/\s+/g, ' ').trim();
    if (!flat) return '';
    let at = 0;
    if (words.length) {
      const low = norm(flat);
      const i = Math.min(...words.map((w) => { const k = wordAt(low, w, 0); return k < 0 ? Infinity : k; }));
      if (i !== Infinity) at = i;
    }
    // norm() keeps the length of the text (one char -> one char), so positions match
    const from = Math.max(0, at - 70);
    const to = Math.min(flat.length, at + 160);
    return (from > 0 ? '…' : '') + mark(flat.slice(from, to), words) + (to < flat.length ? '…' : '');
  }

  function mark(text, words) {
    const s = String(text || '');
    if (!words.length) return esc(s);
    const low = norm(s);
    const ranges = [];
    for (const w of words) {
      let i = wordAt(low, w, 0);
      while (i >= 0) { ranges.push([i, i + w.length]); i = wordAt(low, w, i + w.length); }
    }
    if (!ranges.length) return esc(s);
    ranges.sort((a, b) => a[0] - b[0]);
    let out = '', pos = 0;
    for (const [a, b] of ranges) {
      if (a < pos) continue;
      out += esc(s.slice(pos, a)) + '<mark>' + esc(s.slice(a, b)) + '</mark>';
      pos = b;
    }
    return out + esc(s.slice(pos));
  }

  // ---------- Tags ----------

  function renderTags() {
    const counts = new Map();
    for (const it of items) if (inScope(it)) for (const tag of it.tags) counts.set(tag, (counts.get(tag) || 0) + 1);
    for (const tag of activeTags) if (!counts.has(tag)) counts.set(tag, 0);
    if (!counts.size) { tagsHost.innerHTML = ''; return; }
    const list = [...counts.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], 'ru'));
    tagsHost.innerHTML = `<span class="search__tags-label">${esc(t('search.tags'))}</span>` + list.map(([tag, n]) => `
      <button class="poem-tag${activeTags.has(tag) ? ' poem-tag--active' : ''}" type="button" data-tag="${esc(tag)}"
              aria-pressed="${activeTags.has(tag)}">#${esc(tag)} <small>${n}</small></button>`).join('');
  }

  function syncUrl() {
    const p = new URLSearchParams();
    if (input.value.trim()) p.set('q', input.value.trim());
    if (activeTags.size) p.set('tag', [...activeTags].join(','));
    if (scopeSel.value !== 'all') p.set('in', scopeSel.value);
    const qs = p.toString();
    history.replaceState(null, '', location.pathname + (qs ? '?' + qs : ''));
  }

  // ---------- Helpers ----------

  // Lower case, ё -> е, every non letter/digit -> space. Keeps the length (char for char),
  // so positions found in the normalized text are valid in the original one.
  function norm(s) {
    return String(s || '').toLowerCase().replace(/ё/g, 'е').replace(/[^\p{L}\p{N}]/gu, ' ');
  }

  // Index of `w` at the start of a word in normalized text (after a space or at 0), from `from`
  function wordAt(low, w, from) {
    let i = low.indexOf(w, from);
    while (i > 0 && low[i - 1] !== ' ') i = low.indexOf(w, i + 1);
    return i;
  }

  // A very light Russian stemmer: drop one common ending, so "луна" also finds "луну",
  // "луны", "лунный". Short words are kept as typed. (ENDINGS is declared at the top.)
  function stem(w) {
    if (w.length < 4 || !/[а-я]$/.test(w)) return w;
    for (const e of ENDINGS) if (w.endsWith(e) && w.length - e.length >= 3) return w.slice(0, -e.length);
    return w;
  }

  function plain(html) {
    const s = String(html || '').replace(/<br\s*\/?>|<\/(p|div|li)>/gi, '\n');
    return window.RichEdit ? RichEdit.plain(s) : s.replace(/<[^>]*>/g, '');
  }

  function trText(tr) {
    return Object.values(tr || {}).map((v) => `${v?.title || ''} ${plain(v?.body)}`).join(' ');
  }

  function firstLine(text) {
    const l = String(text || '').split('\n').map((s) => s.trim()).find(Boolean) || '…';
    return l.length > 80 ? l.slice(0, 77) + '…' : l;
  }

  function splitTags(s) {
    return String(s || '').split(/[,;]/).map((x) => x.trim().replace(/^#+/, '')).filter(Boolean);
  }

  async function selectAll(table) {
    if (!db) return [];
    const rows = [];
    for (let from = 0; ; from += 1000) {
      const { data, error } = await db.from(table).select('*').range(from, from + 999);
      if (error) { console.warn(`${table} failed`, error); return rows; }
      rows.push(...(data || []));
      if (!data || data.length < 1000) return rows;
    }
  }

  // A file that still fails after the retries is skipped: search the rest
  async function fetchText(url) {
    const r = await siteGet(url, { optional: true });
    return r ? r.text() : null;
  }

  async function fetchJson(url) {
    try {
      const r = await siteGet(url, { optional: true });
      return r ? await r.json() : null;
    } catch (e) {
      return null;
    }
  }

  function esc(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;' }[ch]));
  }
})();
