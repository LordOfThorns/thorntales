/* polar-night.js — the Polar Night cycle (polar-night.html).

   polar-night.html                         all volumes + "Notes on the universe" (a stub for now)
   polar-night.html?v=<volume>              one volume: its chapters and author notes
   polar-night.html?v=<volume>&e=<entry>    one chapter (with comments) or note
   Admin only:  ?v=new · ?v=<v>&edit · ?v=<v>&e=new&kind=chapter|note · ?v=<v>&e=<e>&edit
                ?v=<v>&e=<e>&history  (every saved version, what changed, restore)

   Data: tables pn_volumes, pn_entries, pn_revisions (supabase/schema.sql). status = 'draft'
   (only the admin), 'members' (signed-in readers) or 'public', checked by RLS.
   If Supabase does not answer, the backup copy
   data/polar-night/book.json is shown ({ volumes: [...], entries: [...] }, public only).
   Original texts are Russian; tr = { en: { title, body }, de: { ... } } holds translations.
   Comments: thread "pn/<volume>/<entry>" (chapters only). Selecting text in a chapter offers
   "Quote" and puts it into the comment form (comments.js renders "> " lines as a quote).
*/
// SiteFetch.get comes from layout.js; an older cached layout.js may not have it yet
window.siteGet = window.siteGet || ((url, o = {}) => (window.SiteFetch
  ? SiteFetch.get(url, o)
  : fetch(url, { cache: 'no-store' }).then((r) => {
    if (r.ok) return r;
    if (o.optional || (o.missingOk && r.status === 404)) return null;
    throw new Error(`${url}: ${r.status}`);
  }, (e) => { if (o.optional) return null; throw e; })));

(function () {
  'use strict';

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const siteLang = window.I18n ? I18n.lang : 'ru';
  const db = window.Auth ? Auth.client : null;
  const LANGS = ['ru', 'en', 'de'];
  const ORIG = 'ru';
  const { sanitize: rawSanitize, plain, esc } = window.RichEdit;
  const sanitize = (html) => rawSanitize(html).trim();

  const root = document.getElementById('pn');
  if (!root) return;

  const params = new URLSearchParams(location.search);
  const P = {
    v: params.get('v'),
    e: params.get('e'),
    kind: params.get('kind') === 'note' ? 'note' : 'chapter',
    section: params.get('section') === 'old' ? 'old' : 'main', // a new volume: main or "old versions"
    edit: params.has('edit'),
    history: params.has('history'),
  };

  // ---------- URLs and small helpers ----------

  const url = (v, e, extra) => {
    const q = new URLSearchParams();
    if (v) q.set('v', v);
    if (e) q.set('e', e);
    for (const [k, val] of Object.entries(extra || {})) q.set(k, val);
    const s = q.toString().replace(/=(?=&|$)/g, ''); // "&edit" rather than "&edit="
    return `${location.pathname}${s ? `?${s}` : ''}`; // polar-night.html (or a test copy of it)
  };

  function roman(n) {
    if (!(n > 0 && n < 4000)) return String(n);
    const m = [[1000, 'M'], [900, 'CM'], [500, 'D'], [400, 'CD'], [100, 'C'], [90, 'XC'],
      [50, 'L'], [40, 'XL'], [10, 'X'], [9, 'IX'], [5, 'V'], [4, 'IV'], [1, 'I']];
    let out = '';
    for (const [v, s] of m) while (n >= v) { out += s; n -= v; }
    return out;
  }

  // "old versions of the books": early drafts, never public (see schema.sql)
  const isOld = (vol) => vol.section === 'old';
  const volLabel = (vol) => (isOld(vol) ? t('pn.oldLabel') : t('pn.volume', { n: roman(vol.num) }));
  const entryLabel = (x) => (x.kind === 'note' ? t('pn.note') : t('pn.chapterN', { n: x.num }));

  // Title / body in the reader's language when translated, else the Russian original
  function version(x, l = siteLang, bodyKey = 'body') {
    const tr = x.tr?.[l];
    if (l !== ORIG && tr && (tr.title || plain(tr.body || '').trim())) {
      return { title: tr.title || x.title, body: tr.body || '', lang: l };
    }
    return { title: x.title, body: x[bodyKey] || '', lang: ORIG };
  }
  function langsOf(x, bodyKey = 'body') {
    const out = [ORIG];
    for (const l of LANGS) if (l !== ORIG && plain(x.tr?.[l]?.body || '').trim()) out.push(l);
    if (!plain(x[bodyKey] || '').trim()) return out.length > 1 ? out : [ORIG];
    return out;
  }

  function header(title, sub) {
    const h1 = document.querySelector('.top h1');
    const p = document.querySelector('.top .sub');
    if (h1) h1.textContent = title;
    if (p) p.innerHTML = sub || '';
    document.title = `${title} · Thorntales`;
  }
  const crumbs = (list) => { if (window.Layout) Layout.setCrumbs(list); };
  // who sees it: drafts only the admin, 'members' only signed-in readers
  // who sees it, as a small label: draft (only the admin), reviewers, signed-in readers
  const STATUS_LABEL = { draft: 'pn.draft', reviewers: 'pn.reviewers', members: 'pn.members' };
  const statusMark = (x) => (STATUS_LABEL[x.status]
    ? ` <span class="pn-status pn-status--${x.status}">${esc(t(STATUS_LABEL[x.status]))}</span>` : '');

  // ---------- Data ----------

  let isAdmin = false;
  let isReviewer = false; // reviewers read the history of edits (the admin is one too)
  let vols = [];
  let entries = []; // list rows: no bodies (see LIST_COLS)

  // Lists need titles only: bodies are fetched for the one entry that is open
  const LIST_COLS = 'id,volume_id,kind,num,title,status,updated_at,tr_en:tr->en->>title,tr_de:tr->de->>title';

  // Supabase may hang when it is down: after a few seconds go on without it (the backup copy)
  const SOFT_MS = 6000;
  const soft = (p, what) => Promise.race([
    Promise.resolve(p),
    new Promise((resolve) => setTimeout(() => resolve({ data: null, error: { message: `${what}: no answer` } }), SOFT_MS)),
  ]);

  async function loadLists() {
    if (db) {
      const [v, e] = await Promise.all([
        soft(db.from('pn_volumes').select('*').order('num'), 'volumes'),
        soft(db.from('pn_entries').select(LIST_COLS).order('num'), 'entries'),
      ]);
      if (!v.error && !e.error) {
        vols = v.data;
        entries = e.data.map((r) => ({ ...r, tr: { en: { title: r.tr_en || '' }, de: { title: r.tr_de || '' } } }));
        return;
      }
      console.warn('polar night tables failed', v.error || e.error);
    }
    const book = await loadBackup();
    vols = book.volumes || [];
    entries = book.entries || [];
  }

  let backup = null;
  async function loadBackup() {
    if (!backup) {
      const r = await siteGet('./data/polar-night/book.json', { missingOk: true });
      backup = r ? await r.json() : { volumes: [], entries: [] };
    }
    return backup;
  }

  async function loadEntry(id) {
    if (db) {
      const { data, error } = await soft(db.from('pn_entries').select('*').eq('id', id).maybeSingle(), 'entry');
      if (!error) return data;
      console.warn('entry failed', error);
    }
    return (await loadBackup()).entries?.find((x) => x.id === id) || null;
  }

  const volById = (id) => vols.find((v) => v.id === id);
  const entriesOf = (volId, kind) => entries
    .filter((x) => x.volume_id === volId && x.kind === kind && (isAdmin || x.status !== 'draft'))
    .sort((a, b) => a.num - b.num || a.title.localeCompare(b.title));

  // ---------- Routing ----------

  async function start() {
    if (window.Auth) await soft(Auth.ready, 'auth');
    isAdmin = !!(window.Auth && Auth.profile?.is_admin);
    isReviewer = isAdmin || !!(window.Auth && Auth.profile?.is_reviewer);
    root.innerHTML = `<p class="text-note">${esc(t('account.loading'))}</p>`;
    try {
      await loadLists();
    } catch (err) {
      root.innerHTML = `<p class="text-note">${esc(t('pn.loadFailed'))}: ${esc(err.message || err)}</p>`;
      return;
    }
    await route();
  }

  // Not found (HTML). A guest may be looking at something only for signed-in readers.
  function notFound() {
    if (window.Auth && Auth.user) return esc(t(isReviewer ? 'pn.notFound' : 'pn.notFoundMember'));
    const back = encodeURIComponent(location.pathname + location.search);
    return `${esc(t('pn.notFoundGuest'))} <a href="/login.html?next=${back}">${esc(t('pn.signIn'))}</a>`;
  }

  async function route() {
    if (P.v === 'new') return isAdmin ? volumeEditor(null) : cycle();
    if (!P.v) return cycle();
    const vol = volById(P.v);
    if (!vol) return cycle(notFound());
    if (!P.e) return P.edit && isAdmin ? volumeEditor(vol) : volume(vol);
    if (P.e === 'new') return isAdmin ? entryEditor(vol, null, P.kind) : volume(vol);
    const x = await loadEntry(P.e);
    if (!x || x.volume_id !== vol.id) return volume(vol, notFound());
    if (P.edit && isAdmin) return entryEditor(vol, x, x.kind);
    if (P.history && isReviewer) return historyView(vol, x);
    return reader(vol, x);
  }

  // ---------- The cycle: all volumes ----------

  function cycle(note) {
    header(t('nav.polarNight'), esc(t('polar.kind')));
    crumbs([]);
    const shown = vols.filter((v) => isAdmin || v.status !== 'draft');
    const card = (v) => {
      const ver = version(v, siteLang, 'descr');
      const desc = plain(ver.body).trim();
      const n = entriesOf(v.id, 'chapter').length;
      return `
        <li><a class="pn-vol${isOld(v) ? ' pn-vol--old' : ''}" href="${esc(url(v.id))}">
          <span class="pn-vol__num">${esc(volLabel(v))}${statusMark(v)}</span>
          <span class="pn-vol__title">${esc(ver.title)}</span>
          ${desc ? `<span class="pn-vol__desc">${esc(desc.length > 220 ? `${desc.slice(0, 220)}…` : desc)}</span>` : ''}
          <span class="pn-vol__count">${esc(t('pn.chapterCount', { n }))}</span>
        </a></li>`;
    };
    const main = shown.filter((v) => !isOld(v));
    const old = shown.filter(isOld);
    // old versions: guests only learn that the section exists (RLS gives them no rows)
    const signedIn = !!(window.Auth && Auth.user);
    const back = encodeURIComponent(location.pathname + location.search);
    let oldBody;
    if (!signedIn) oldBody = `<p class="text-note">${esc(t('pn.oldGuest'))} <a href="/login.html?next=${back}">${esc(t('pn.signIn'))}</a></p>`;
    else if (old.length) oldBody = `<ol class="pn-vols pn-vols--old">${old.map(card).join('')}</ol>`;
    else oldBody = `<p class="text-note">${esc(t('pn.noOld'))}</p>`;
    root.innerHTML = `
      ${note ? `<p class="text-note">${note}</p>` : ''}
      ${isAdmin ? `<div class="text-admin"><a class="btn" href="${esc(url('new'))}">${esc(t('pn.newVolume'))}</a></div>` : ''}
      ${main.length ? `<ol class="pn-vols">${main.map(card).join('')}</ol>` : `<p class="text-note">${esc(t('pn.noVolumes'))}</p>`}
      <section class="pn-old">
        <div class="pn-part__head">
          <h2 class="pn-part__title">${esc(t('pn.old'))}</h2>
          ${isAdmin ? `<a class="btn" href="${esc(url('new', null, { section: 'old' }))}">${esc(t('pn.newOld'))}</a>` : ''}
        </div>
        <p class="pn-old__sub">${esc(t('pn.oldSub'))}</p>
        ${oldBody}
      </section>
      <section class="pn-universe">
        <h2 class="pn-universe__title">${esc(t('pn.universe'))}</h2>
        <p class="pn-universe__text">${esc(t('pn.universeStub'))}</p>
      </section>`;
  }

  // ---------- One volume: contents ----------

  function volume(vol, note) {
    const ver = version(vol, siteLang, 'descr');
    header(ver.title, `<span class="text-meta">${esc(volLabel(vol))}</span>${statusMark(vol)}
      ${isAdmin ? ` <a class="poem__edit-btn" href="${esc(url(vol.id, null, { edit: '' }))}">${esc(t('pn.editVolume'))}</a>` : ''}`);
    crumbs([{ label: ver.title }]);
    const list = (kind) => {
      const items = entriesOf(vol.id, kind);
      if (!items.length) return `<p class="text-note">${esc(t(kind === 'note' ? 'pn.noNotes' : 'pn.noChapters'))}</p>`;
      return `<ol class="pn-toc${kind === 'note' ? ' pn-toc--notes' : ''}">${items.map((x) => `
        <li class="pn-toc__item is-${esc(x.status || 'public')}"><a href="${esc(url(vol.id, x.id))}">
          ${kind === 'chapter' ? `<span class="pn-toc__num">${esc(String(x.num))}</span>` : ''}
          <span class="pn-toc__title">${esc(version(x).title)}</span>${statusMark(x)}
        </a></li>`).join('')}</ol>`;
    };
    const add = (kind, key) => (isAdmin ? `<a class="btn" href="${esc(url(vol.id, 'new', { kind }))}">${esc(t(key))}</a>` : '');
    root.innerHTML = `
      ${note ? `<p class="text-note">${note}</p>` : ''}
      ${plain(ver.body).trim() ? `<div class="text-body pn-vol-desc" lang="${ver.lang}">${sanitize(ver.body)}</div>` : ''}
      <section class="pn-part">
        <div class="pn-part__head"><h2 class="pn-part__title">${esc(t('pn.chapters'))}</h2>${add('chapter', 'pn.newChapter')}</div>
        ${list('chapter')}
      </section>
      <section class="pn-part">
        <div class="pn-part__head"><h2 class="pn-part__title">${esc(t('pn.notes'))}</h2>${add('note', 'pn.newNote')}</div>
        ${list('note')}
      </section>
      <p class="texts-back"><a href="${esc(url())}">${esc(t('pn.toCycle'))}</a></p>`;
  }

  // ---------- One chapter or note ----------

  function reader(vol, x) {
    let shown = langsOf(x).includes(siteLang) ? siteLang : ORIG;
    const siblings = entriesOf(vol.id, x.kind);
    const i = siblings.findIndex((y) => y.id === x.id);
    const prev = siblings[i - 1];
    const next = siblings[i + 1];
    const volTitle = version(vol, siteLang, 'descr').title;
    const thread = `pn/${vol.id}/${x.id}`;

    const draw = () => {
      const v = version(x, shown);
      const langs = langsOf(x);
      let langControl = '';
      if (langs.length > 1) {
        langControl = `<span class="text-langs" role="group">${langs.map((l) => `
          <button class="text-langs__btn" type="button" data-act="lang" data-lang="${l}" aria-pressed="${l === shown}">${esc(l.toUpperCase())}${l === ORIG ? ` <span class="text-langs__orig">· ${esc(t('texts.original'))}</span>` : ''}</button>`).join('')}</span>`;
      }
      if (!langs.includes(siteLang)) langControl += ` <span class="poem__untranslated">${esc(t('poems.notTranslated'))}</span>`;
      // the admin edits; reviewers (and the admin) read the history of edits
      const admin = `${isAdmin ? `
        <a class="poem__edit-btn" href="${esc(url(vol.id, x.id, { edit: '' }))}">${esc(t('pn.edit'))}</a>` : ''}${isReviewer ? `
        <a class="poem__edit-btn" href="${esc(url(vol.id, x.id, { history: '' }))}">${esc(t('pn.history'))}</a>` : ''}`;
      header(v.title, `<span class="text-meta">${esc(volLabel(vol))} · ${esc(entryLabel(x))}</span>${statusMark(x)} ${langControl} ${admin}`);
      crumbs([{ label: volTitle, href: url(vol.id) }, { label: v.title }]);
      root.innerHTML = `
        <article class="text-body pn-text" lang="${v.lang}">${sanitize(v.body)}</article>
        ${x.kind === 'chapter' ? '<div class="text-comments" id="pnComments"></div>' : ''}
        <nav class="catalog-pager">
          ${prev ? `<a class="catalog-pager__link catalog-pager__link--prev" href="${esc(url(vol.id, prev.id))}">
            <span class="catalog-pager__hint">${esc(t('pn.prev'))}</span>
            <span class="catalog-pager__label">${esc(version(prev).title)}</span></a>` : '<span></span>'}
          <a class="catalog-pager__link catalog-pager__link--all" href="${esc(url(vol.id))}">${esc(t('pn.toVolume'))}</a>
          ${next ? `<a class="catalog-pager__link catalog-pager__link--next" href="${esc(url(vol.id, next.id))}">
            <span class="catalog-pager__hint">${esc(t('pn.next'))}</span>
            <span class="catalog-pager__label">${esc(version(next).title)}</span></a>` : '<span></span>'}
        </nav>`;
      if (canQuote) Comments.mount([{ el: document.getElementById('pnComments'), thread }], thread);
      marks.apply(v.lang);
    };
    const canQuote = x.kind === 'chapter' && x.status !== 'draft' && !!(window.Comments && Comments.quote);
    const marks = makeMarks(x);
    draw();
    document.querySelector('.top').onclick = (e) => {
      const b = e.target.closest('[data-act="lang"]');
      if (b) { shown = b.dataset.lang; draw(); }
    };
    // what a selection in the text offers: quote it into a comment; for the admin, a margin mark
    const actions = [];
    if (canQuote) {
      actions.push({ label: t('pn.quote'), run: (range, text) => {
        getSelection().removeAllRanges();
        Comments.quote(thread, text.length > 1500 ? `${text.slice(0, 1500)}…` : text);
      } });
    }
    if (isAdmin && db) actions.push({ label: t('pn.markAdd'), run: (range) => marks.add(range) });
    selectionBar(actions);
  }

  // A small bar under a selection in the text with the given actions
  function selectionBar(actions) {
    if (!actions.length) return;
    const bar = document.createElement('div');
    bar.className = 'pn-selbar';
    bar.hidden = true;
    for (const a of actions) {
      const b = document.createElement('button');
      b.type = 'button';
      b.textContent = a.label;
      b.addEventListener('mousedown', (e) => e.preventDefault()); // keep the selection until the click
      b.addEventListener('click', () => {
        const sel = getSelection();
        if (!sel.rangeCount) return;
        const range = sel.getRangeAt(0).cloneRange();
        bar.hidden = true;
        a.run(range, sel.toString().replace(/\n{3,}/g, '\n\n').trim());
      });
      bar.append(b);
    }
    document.body.append(bar);
    const place = () => {
      const sel = getSelection();
      const art = root.querySelector('.pn-text');
      if (!art || !sel || sel.isCollapsed || !sel.rangeCount
        || !art.contains(sel.anchorNode) || !art.contains(sel.focusNode) || !sel.toString().trim()) { bar.hidden = true; return; }
      const r = sel.getRangeAt(0).getBoundingClientRect();
      bar.hidden = false;
      const top = r.bottom + 10 + bar.offsetHeight < innerHeight ? r.bottom + 10 : r.top - bar.offsetHeight - 10;
      bar.style.top = `${Math.max(8, top)}px`;
      bar.style.left = `${Math.min(innerWidth - bar.offsetWidth - 8, Math.max(8, r.left + r.width / 2 - bar.offsetWidth / 2))}px`;
    };
    let timer = 0;
    document.addEventListener('selectionchange', () => { clearTimeout(timer); timer = setTimeout(place, 120); });
    addEventListener('scroll', () => { bar.hidden = true; }, { passive: true });
  }

  // ---------- Margin marks: the author's notes on a fragment of the text ----------
  // Table pn_marks. A mark is found again by its fragment (quote) plus a little text
  // before and after it, in the language version it was made in; one whose fragment is
  // gone after an edit is listed under the text for the admin ("its place changed").
  // The icon is a CSS ::after, so marks never change the text the offsets are counted in.

  const CONTEXT = 48;

  function makeMarks(x) {
    let list = [];
    let lang = ORIG;
    let pop = null;       // the open popup
    let pinned = false;   // opened by a click (stays) or by hovering (closes on leave)
    let hoverTimer = 0;
    const load = async () => {
      if (!db) return;
      const { data, error } = await soft(db.from('pn_marks').select('*').eq('entry_id', x.id).order('created_at'), 'marks');
      if (error) { console.warn('marks failed', error); return; } // no table yet: no marks
      list = data || [];
    };
    let loaded = load();
    const art = () => root.querySelector('.pn-text');
    const byId = (id) => list.find((m) => String(m.id) === String(id));

    // offset of a point (node, offset) in the text of the article
    function indexOf(a, node, offset) {
      const r = document.createRange();
      r.setStart(a, 0);
      r.setEnd(node, offset);
      return r.toString().length;
    }

    // where the fragment is now: every occurrence, the best fitting context wins
    function locate(S, m) {
      let best = -1;
      let bestScore = -1;
      for (let at = S.indexOf(m.quote); at >= 0; at = S.indexOf(m.quote, at + 1)) {
        let score = 0;
        for (let k = 1; k <= m.prefix.length && S[at - k] === m.prefix[m.prefix.length - k]; k++) score++;
        const end = at + m.quote.length;
        for (let k = 0; k < m.suffix.length && S[end + k] === m.suffix[k]; k++) score++;
        if (score > bestScore) { best = at; bestScore = score; }
      }
      return best;
    }

    // wrap text [start, end) in highlighted spans (one per text node it crosses)
    function wrap(a, start, end, m) {
      const walker = document.createTreeWalker(a, NodeFilter.SHOW_TEXT);
      const parts = [];
      let pos = 0;
      for (let n = walker.nextNode(); n; n = walker.nextNode()) {
        const len = n.data.length;
        const s = Math.max(start, pos);
        const e = Math.min(end, pos + len);
        if (s < e) parts.push([n, s - pos, e - pos]);
        pos += len;
        if (pos >= end) break;
      }
      parts.forEach(([n, s, e], i) => {
        let mid = n;
        if (e < mid.data.length) mid.splitText(e);
        if (s > 0) mid = mid.splitText(s);
        const span = document.createElement('span');
        span.className = `pn-mark pn-mark--${m.status}${i === parts.length - 1 ? ' pn-mark--end' : ''}`;
        span.dataset.mark = m.id;
        if (i === parts.length - 1) {
          span.tabIndex = 0;
          span.setAttribute('role', 'button');
          span.setAttribute('aria-label', t('pn.markTitle'));
        }
        mid.replaceWith(span);
        span.append(mid);
      });
    }

    function render() {
      const a = art();
      if (!a) return;
      close();
      a.querySelectorAll('.pn-mark').forEach((s) => s.replaceWith(...s.childNodes));
      a.normalize();
      const S = a.textContent;
      const found = [];
      const lost = [];
      for (const m of list.filter((y) => y.lang === lang)) {
        const at = locate(S, m);
        if (at < 0) lost.push(m); else found.push([at, m]);
      }
      for (const [at, m] of found) wrap(a, at, at + m.quote.length, m);
      renderLost(lost);
    }

    function renderLost(lost) {
      root.querySelector('.pn-lost')?.remove();
      if (!isAdmin || !lost.length) return;
      const box = document.createElement('section');
      box.className = 'pn-lost';
      box.innerHTML = `<h3 class="pn-lost__title">${esc(t('pn.marksLost'))}</h3>
        <ul class="pn-lost__list">${lost.map((m) => `
          <li class="pn-lost__item" data-mark="${m.id}">
            <blockquote class="comment__quote">${esc(m.quote)}</blockquote>
            <p class="pn-lost__note">${esc(m.note)}</p>
          </li>`).join('')}</ul>`;
      art().after(box);
    }

    // ----- the popup -----
    function close() {
      clearTimeout(hoverTimer);
      pop?.remove();
      pop = null;
      pinned = false;
    }
    function placeAt(el, rect) {
      const w = Math.min(340, innerWidth - 16);
      el.style.width = `${w}px`;
      el.style.left = `${Math.max(8, Math.min(innerWidth - w - 8, rect.left)) + scrollX}px`;
      el.style.top = `${rect.bottom + 8 + scrollY}px`;
    }
    const stamp = (iso) => new Date(iso).toLocaleDateString(window.I18n ? I18n.locale : 'ru-RU', { day: 'numeric', month: 'long', year: 'numeric' });

    function show(m, anchor, byClick) {
      close();
      pinned = byClick;
      pop = document.createElement('div');
      pop.className = 'pn-pop';
      pop.setAttribute('role', 'dialog');
      pop.innerHTML = `
        <div class="pn-pop__head"><span class="pn-pop__title">${esc(t('pn.markTitle'))}</span>${isAdmin ? markWho(m.status) : ''}</div>
        <div class="pn-pop__note">${esc(m.note)}</div>
        <div class="pn-pop__foot">
          <time datetime="${esc(m.created_at)}">${esc(stamp(m.created_at))}</time>
          ${isAdmin ? `<span class="pn-pop__acts">
            <button class="link-btn" type="button" data-act="edit">${esc(t('pn.markEdit'))}</button>
          </span>` : ''}
        </div>`;
      document.body.append(pop);
      placeAt(pop, anchor.getBoundingClientRect());
      pop.addEventListener('mouseenter', () => clearTimeout(hoverTimer));
      pop.addEventListener('mouseleave', () => { if (!pinned) hoverTimer = setTimeout(close, 300); });
      pop.querySelector('[data-act="edit"]')?.addEventListener('click', () => form(m, anchor.getBoundingClientRect()));
    }
    const markWho = (s) => `<span class="pn-status pn-status--${STATUS_LABEL[s] ? s : 'members'}">${esc(t(s === 'draft' ? 'pn.markOnlyMe' : STATUS_LABEL[s] || 'pn.members'))}</span>`;

    // the form: a new mark (with its fragment) or an existing one
    function form(m, rect, fresh) {
      close();
      pinned = true;
      pop = document.createElement('form');
      pop.className = 'pn-pop pn-pop--form';
      const st = m?.status || 'draft';
      pop.innerHTML = `
        <div class="pn-pop__head"><span class="pn-pop__title">${esc(t('pn.markTitle'))}</span></div>
        <blockquote class="comment__quote pn-pop__quote">${esc((m || fresh).quote.length > 300 ? `${(m || fresh).quote.slice(0, 300)}…` : (m || fresh).quote)}</blockquote>
        <textarea class="input pn-pop__text" name="note" rows="4" maxlength="8000" required placeholder="${esc(t('pn.markPh'))}"></textarea>
        <div class="pn-pop__who">
          ${['draft', 'reviewers', 'members'].map((s) => `<label class="pn-edit__opt"><input type="radio" name="status" value="${s}"${s === st ? ' checked' : ''}>
            <span>${esc(t(s === 'draft' ? 'pn.markOnlyMe' : `pn.status.${s}`))}</span></label>`).join('')}
        </div>
        <p class="poem-edit__error" hidden></p>
        <div class="pn-pop__foot">
          <button class="btn btn--primary" type="submit">${esc(t('pn.save'))}</button>
          <button class="btn btn--ghost" type="button" data-act="cancel">${esc(t('pn.cancel'))}</button>
          ${m ? `<button class="link-btn pn-pop__del" type="button" data-act="delete">${esc(t('pn.markDelete'))}</button>` : ''}
        </div>`;
      document.body.append(pop);
      placeAt(pop, rect);
      const f = pop;
      f.elements.note.value = m?.note || '';
      f.elements.note.focus({ preventScroll: true });
      const fail = (msg) => { const p = f.querySelector('.poem-edit__error'); p.textContent = msg; p.hidden = false; };
      f.querySelector('[data-act="cancel"]').addEventListener('click', close);
      const del = f.querySelector('[data-act="delete"]');
      del?.addEventListener('click', async () => {
        // two steps instead of a browser confirm() (some browsers block those)
        if (!del.classList.contains('is-sure')) { del.classList.add('is-sure'); del.textContent = t('pn.markDeleteSure'); return; }
        const { error } = await db.from('pn_marks').delete().eq('id', m.id);
        if (error) return fail(failed(error));
        await reload();
      });
      f.addEventListener('submit', async (e) => {
        e.preventDefault();
        const note = f.elements.note.value.trim();
        if (!note) return;
        const status = f.elements.status.value;
        f.querySelectorAll('button, textarea, input').forEach((el) => { el.disabled = true; });
        const { error } = m
          ? await db.from('pn_marks').update({ note, status }).eq('id', m.id)
          : await db.from('pn_marks').insert({ entry_id: x.id, lang, ...fresh, note, status });
        if (error) {
          f.querySelectorAll('button, textarea, input').forEach((el) => { el.disabled = false; });
          return fail(failed(error));
        }
        await reload();
      });
    }

    async function reload() {
      close();
      loaded = load();
      await loaded;
      render();
    }

    // clicks and hovering on marks (the root is re-rendered, the listeners stay)
    root.addEventListener('click', (e) => {
      const s = e.target.closest('.pn-mark, .pn-lost__item');
      if (!s) return;
      const m = byId(s.dataset.mark);
      if (!m) return;
      if (s.classList.contains('pn-lost__item')) { form(m, s.getBoundingClientRect()); return; }
      const end = root.querySelector(`.pn-mark--end[data-mark="${m.id}"]`) || s;
      show(m, end, true);
    });
    root.addEventListener('keydown', (e) => {
      if ((e.key === 'Enter' || e.key === ' ') && e.target.classList?.contains('pn-mark--end')) {
        e.preventDefault();
        const m = byId(e.target.dataset.mark);
        if (m) show(m, e.target, true);
      }
    });
    root.addEventListener('mouseover', (e) => {
      const s = e.target.closest('.pn-mark');
      if (!s || (pop && pinned) || !matchMedia('(hover: hover)').matches) return;
      clearTimeout(hoverTimer);
      hoverTimer = setTimeout(() => {
        const m = byId(s.dataset.mark);
        const end = root.querySelector(`.pn-mark--end[data-mark="${s.dataset.mark}"]`) || s;
        if (m && !(pop && pinned)) show(m, end, false);
      }, 250);
    });
    root.addEventListener('mouseout', (e) => {
      if (!e.target.closest('.pn-mark') || pinned) return;
      clearTimeout(hoverTimer);
      hoverTimer = setTimeout(() => { if (!pinned) close(); }, 300);
    });
    document.addEventListener('mousedown', (e) => {
      if (pop && !pop.contains(e.target) && !e.target.closest('.pn-mark, .pn-selbar, .pn-lost__item')) close();
    });
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape') close(); });

    return {
      async apply(l) {
        lang = l;
        await loaded;
        render();
      },
      // a new mark from a selection: the fragment and a little text around it
      add(range) {
        const a = art();
        if (!a) return;
        const S = a.textContent;
        const start = indexOf(a, range.startContainer, range.startOffset);
        const end = indexOf(a, range.endContainer, range.endOffset);
        const quote = S.slice(start, end);
        if (!quote.trim()) return;
        const rect = range.getBoundingClientRect();
        getSelection().removeAllRanges();
        form(null, rect, {
          quote: quote.slice(0, 4000),
          prefix: S.slice(Math.max(0, start - CONTEXT), start),
          suffix: S.slice(end, end + CONTEXT),
        });
      },
    };
  }

  // ---------- Editors (admin) ----------

  // One form for a volume (its description) and for a chapter or note (its text):
  // number, who sees it, Russian original and translations in tabs.
  function editorForm({ title, crumbsList, item, bodyKey, bodyLabel, numDefault, cancelHref, onSave, onDelete,
    statuses = ['draft', 'members', 'public'] }) {
    const status = statuses.includes(item?.status) ? item.status : 'draft';
    const st = { v: {} };
    for (const l of LANGS) st.v[l] = { title: '', body: '' };
    if (item) {
      st.v[ORIG] = { title: item.title, body: item[bodyKey] || '' };
      for (const [l, v] of Object.entries(item.tr || {})) if (l !== ORIG) st.v[l] = { title: v.title || '', body: v.body || '' };
    }
    let tab = ORIG;
    let field = null;
    header(title, '');
    crumbs(crumbsList);
    root.innerHTML = `
      <form class="panel text-edit pn-edit">
        <div class="text-edit__row">
          <label class="field">
            <span class="field__label">${esc(t('pn.fieldNum'))}</span>
            <input class="input" name="num" type="number" min="1" max="9999" step="1" value="${esc(String(item?.num ?? numDefault))}">
          </label>
        </div>
        <fieldset class="pn-edit__who">
          <legend class="field__label">${esc(t('pn.whoSees'))}</legend>
          ${statuses.map((s) => `
            <label class="pn-edit__opt">
              <input type="radio" name="status" value="${s}"${s === status ? ' checked' : ''}>
              <span>${esc(t(`pn.status.${s}`))}</span>
            </label>`).join('')}
        </fieldset>
        <div class="text-edit__tabs" role="tablist"></div>
        <label class="field">
          <span class="field__label">${esc(t('pn.fieldTitle'))}</span>
          <input class="input" name="title" maxlength="300">
        </label>
        <div class="field">
          <span class="field__label">${esc(bodyLabel)}</span>
          <div class="text-edit__body"></div>
        </div>
        <p class="poem-edit__error" hidden></p>
        <div class="poem-edit__actions">
          <button class="btn btn--primary" type="submit">${esc(t('pn.save'))}</button>
          <a class="btn btn--ghost" href="${esc(cancelHref)}">${esc(t('pn.cancel'))}</a>
          <span class="poem-edit__spacer"></span>
          ${onDelete ? `<button class="btn btn--ghost" type="button" data-act="delete">${esc(t('pn.delete'))}</button>` : ''}
        </div>
      </form>`;
    const form = root.querySelector('form');
    const f = form.elements;
    const saveTab = () => { st.v[tab] = { title: f.title.value.trim(), body: field.get() }; };
    const drawTabs = () => {
      form.querySelector('.text-edit__tabs').innerHTML = LANGS.map((l) => {
        const kind = l === ORIG ? t('texts.tabOriginal') : t('texts.tabTranslation');
        const filled = st.v[l].title || plain(st.v[l].body).trim() ? ' is-filled' : '';
        return `<button type="button" role="tab" class="text-edit__tab${l === tab ? ' is-active' : ''}${filled}" data-tab="${l}" aria-selected="${l === tab}">${l.toUpperCase()} <small>${esc(kind)}</small></button>`;
      }).join('');
    };
    const openTab = (l) => {
      if (field) saveTab();
      tab = l;
      f.title.value = st.v[l].title;
      field = RichEdit.field(form.querySelector('.text-edit__body'), st.v[l].body);
      drawTabs();
    };
    openTab(ORIG);

    const busy = (on) => form.querySelectorAll('button, input, select, textarea').forEach((el) => { el.disabled = on; });
    const fail = (msg) => {
      busy(false);
      const box = form.querySelector('.poem-edit__error');
      box.textContent = msg;
      box.hidden = false;
    };

    form.addEventListener('click', async (e) => {
      const tb = e.target.closest('[data-tab]');
      if (tb) { openTab(tb.dataset.tab); return; }
      if (e.target.closest('[data-act="delete"]')) {
        busy(true);
        const err = await onDelete();
        if (err) fail(err);
        else busy(false);
      }
    });

    form.addEventListener('submit', async (e) => {
      e.preventDefault();
      saveTab();
      if (!st.v[ORIG].title) return fail(t('pn.needTitle'));
      const tr = {};
      for (const l of LANGS) {
        if (l === ORIG) continue;
        const v = st.v[l];
        if (v.title || plain(v.body).trim()) tr[l] = { title: v.title, body: v.body };
      }
      busy(true);
      const err = await onSave({
        num: Math.max(1, Math.min(9999, Math.round(Number(f.num.value) || 1))),
        title: st.v[ORIG].title,
        body: st.v[ORIG].body,
        tr,
        status: f.status.value,
      });
      if (err) fail(err);
    });
  }

  const nextNum = (list) => list.reduce((m, x) => Math.max(m, x.num || 0), 0) + 1;
  const newId = (prefix) => `${prefix}-${Date.now().toString(36)}${Math.random().toString(36).slice(2, 5)}`;
  const failed = (error) => t('pn.saveFailed') + (error.message || error);

  // old versions of the books are never public
  const statusesFor = (vol) => (isOld(vol) ? ['draft', 'reviewers', 'members'] : ['draft', 'reviewers', 'members', 'public']);

  function volumeEditor(vol) {
    const section = vol?.section || P.section;
    const newTitle = t(section === 'old' ? 'pn.newOldTitle' : 'pn.newVolumeTitle');
    editorForm({
      title: vol ? version(vol, ORIG, 'descr').title : newTitle,
      crumbsList: [{ label: vol ? vol.title : newTitle }],
      item: vol,
      bodyKey: 'descr',
      bodyLabel: t('pn.fieldDescr'),
      numDefault: nextNum(vols.filter((v) => (v.section || 'main') === section)),
      cancelHref: vol ? url(vol.id) : url(),
      statuses: statusesFor({ section }),
      async onSave(d) {
        let id = vol?.id;
        if (!id) {
          const taken = new Set(vols.map((v) => v.id));
          const base = `${section === 'old' ? 'old' : 'tom'}-${d.num}`;
          id = base;
          for (let k = 2; taken.has(id); k++) id = `${base}-${k}`;
        }
        const row = { id, section, num: d.num, title: d.title, descr: d.body, tr: d.tr, status: d.status };
        const { error } = vol
          ? await db.from('pn_volumes').update(row).eq('id', id)
          : await db.from('pn_volumes').insert(row);
        if (error) return failed(error);
        location.href = url(id);
        return null;
      },
      onDelete: vol ? async () => {
        if (entries.some((x) => x.volume_id === vol.id)) return t('pn.volumeNotEmpty');
        if (!confirm(t('pn.deleteVolumeAsk', { title: vol.title }))) return null;
        const { error } = await db.from('pn_volumes').delete().eq('id', vol.id);
        if (error) return failed(error);
        location.href = url();
        return null;
      } : null,
    });
  }

  function entryEditor(vol, x, kind) {
    const isNote = kind === 'note';
    const newTitle = t(isNote ? 'pn.newNoteTitle' : 'pn.newChapterTitle');
    editorForm({
      title: x ? x.title : newTitle,
      crumbsList: [{ label: vol.title, href: url(vol.id) }, { label: x ? x.title : newTitle }],
      item: x,
      bodyKey: 'body',
      bodyLabel: t('pn.fieldBody'),
      numDefault: nextNum(entries.filter((y) => y.volume_id === vol.id && y.kind === kind)),
      cancelHref: x ? url(vol.id, x.id) : url(vol.id),
      statuses: statusesFor(vol),
      async onSave(d) {
        const id = x?.id || newId(isNote ? 'z' : 'g');
        const row = { id, volume_id: vol.id, kind, num: d.num, title: d.title, body: d.body, tr: d.tr, status: d.status };
        const { error } = x
          ? await db.from('pn_entries').update(row).eq('id', id)
          : await db.from('pn_entries').insert(row);
        if (error) return failed(error);
        location.href = url(vol.id, id);
        return null;
      },
      onDelete: x ? async () => {
        if (!confirm(t('pn.deleteEntryAsk', { title: x.title }))) return null;
        const { error } = await db.from('pn_entries').delete().eq('id', x.id);
        if (error) return failed(error);
        location.href = url(vol.id);
        return null;
      } : null,
    });
  }

  // ---------- History of an entry (admin) ----------

  async function historyView(vol, x) {
    header(t('pn.historyTitle', { title: x.title }), `<span class="text-meta">${esc(volLabel(vol))} · ${esc(entryLabel(x))}</span>
      <a class="poem__edit-btn" href="${esc(url(vol.id, x.id))}">${esc(t('pn.backToEntry'))}</a>`);
    crumbs([{ label: vol.title, href: url(vol.id) }, { label: x.title, href: url(vol.id, x.id) }, { label: t('pn.history') }]);
    root.innerHTML = `<p class="text-note">${esc(t('account.loading'))}</p>`;
    const { data, error } = await db.from('pn_revisions')
      .select('id,saved_at,title,body,tr').eq('entry_id', x.id)
      .order('saved_at', { ascending: false }).limit(300);
    if (error) { root.innerHTML = `<p class="text-note">${esc(t('pn.loadFailed'))}: ${esc(error.message)}</p>`; return; }
    const revs = data;
    if (!revs.length) { root.innerHTML = `<p class="text-note">${esc(t('pn.historyEmpty'))}</p>`; return; }

    let sel = 0;
    let lang = ORIG;
    let whole = false;
    const textOf = (r, l) => (l === ORIG ? { title: r.title, body: r.body } : { title: r.tr?.[l]?.title || '', body: r.tr?.[l]?.body || '' });
    const stamp = (iso) => new Date(iso).toLocaleString(window.I18n ? I18n.locale : 'ru-RU', {
      day: 'numeric', month: 'short', year: 'numeric', hour: '2-digit', minute: '2-digit',
    });
    // how much changed in the original: words added / removed against the version before
    const counts = revs.map((r, i) => {
      const older = revs[i + 1];
      if (!older) return null;
      return Diff.count(paras(older.body), paras(r.body));
    });

    const draw = () => {
      const r = revs[sel];
      const older = revs[sel + 1];
      const now = textOf(r, lang);
      const was = older ? textOf(older, lang) : { title: '', body: '' };
      let diffHtml;
      if (!older) {
        diffHtml = `<p class="pn-diff__note">${esc(t('pn.firstVersion'))}</p>${Diff.render([], paras(now.body), true)}`;
      } else {
        const titleChange = was.title !== now.title
          ? `<p class="pn-diff__title"><span>${esc(t('pn.titleChanged'))}:</span> <del>${esc(was.title)}</del> → <ins>${esc(now.title)}</ins></p>` : '';
        const body = Diff.render(paras(was.body), paras(now.body), whole);
        diffHtml = titleChange + (body || (titleChange ? '' : `<p class="pn-diff__note">${esc(t('pn.noChanges'))}</p>`));
      }
      root.innerHTML = `
        <div class="pn-history">
          <ol class="pn-history__list">${revs.map((rv, i) => {
            const c = counts[i];
            return `<li><button type="button" class="pn-history__item${i === sel ? ' is-active' : ''}" data-rev="${i}">
              <span class="pn-history__date">${esc(stamp(rv.saved_at))}</span>
              <span class="pn-history__count">${c ? `<ins>+${c.added}</ins> <del>−${c.removed}</del>` : esc(t('pn.firstVersion'))}</span>
            </button></li>`;
          }).join('')}</ol>
          <div class="pn-history__view">
            <div class="pn-history__bar">
              <span class="text-langs" role="group">${LANGS.map((l) => `<button class="text-langs__btn" type="button" data-lang="${l}" aria-pressed="${l === lang}">${l.toUpperCase()}</button>`).join('')}</span>
              <label class="lang-only"><input type="checkbox" data-act="whole"${whole ? ' checked' : ''}><span>${esc(t('pn.showAll'))}</span></label>
              ${sel > 0 && isAdmin ? `<button class="btn" type="button" data-act="restore">${esc(t('pn.restore'))}</button>` : ''}
            </div>
            <div class="pn-diff text-body">${diffHtml}</div>
          </div>
        </div>`;
    };
    draw();
    root.onclick = async (e) => {
      const rv = e.target.closest('[data-rev]');
      if (rv) { sel = Number(rv.dataset.rev); draw(); return; }
      const lb = e.target.closest('[data-lang]');
      if (lb) { lang = lb.dataset.lang; draw(); return; }
      if (e.target.closest('[data-act="restore"]')) {
        const r = revs[sel];
        if (!confirm(t('pn.restoreAsk', { date: stamp(r.saved_at) }))) return;
        const { error: err } = await db.from('pn_entries').update({ title: r.title, body: r.body, tr: r.tr }).eq('id', x.id);
        if (err) { alert(failed(err)); return; }
        location.href = url(vol.id, x.id, { history: '' });
      }
    };
    root.onchange = (e) => {
      if (e.target.matches('[data-act="whole"]')) { whole = e.target.checked; draw(); }
    };
  }

  // Paragraphs of a text as plain strings (the unit the history compares first)
  function paras(html) {
    const box = new DOMParser().parseFromString(`<div>${html || ''}</div>`, 'text/html').body.firstChild;
    const out = [];
    let loose = '';
    for (const n of box.childNodes) {
      const s = (n.textContent || '').replace(/\s+/g, ' ').trim();
      if (n.nodeType === 1 && /^(P|H3|H4|BLOCKQUOTE|LI|UL|OL|HR|DIV)$/.test(n.tagName)) {
        if (loose.trim()) { out.push(loose.trim()); loose = ''; }
        if (s) out.push(s);
      } else {
        loose += ` ${s}`;
      }
    }
    if (loose.trim()) out.push(loose.trim());
    return out;
  }

  // ---------- Diff: paragraphs first, then words inside changed paragraphs ----------

  const Diff = (() => {
    // Longest common subsequence on two arrays -> ops [['=', a, b] | ['-', a] | ['+', b]]
    function ops(a, b) {
      let s = 0;
      while (s < a.length && s < b.length && a[s] === b[s]) s++;
      let ea = a.length;
      let eb = b.length;
      while (ea > s && eb > s && a[ea - 1] === b[eb - 1]) { ea--; eb--; }
      const out = [];
      for (let i = 0; i < s; i++) out.push(['=', a[i], b[i]]);
      const A = a.slice(s, ea);
      const B = b.slice(s, eb);
      const n = A.length;
      const m = B.length;
      if (n * m > 4e6) { // too big to compare in the browser: show as replaced
        for (const x of A) out.push(['-', x]);
        for (const y of B) out.push(['+', y]);
      } else {
        const L = new Uint32Array((n + 1) * (m + 1));
        for (let i = n - 1; i >= 0; i--) {
          for (let j = m - 1; j >= 0; j--) {
            L[i * (m + 1) + j] = A[i] === B[j] ? L[(i + 1) * (m + 1) + j + 1] + 1
              : Math.max(L[(i + 1) * (m + 1) + j], L[i * (m + 1) + j + 1]);
          }
        }
        let i = 0;
        let j = 0;
        while (i < n && j < m) {
          if (A[i] === B[j]) { out.push(['=', A[i], B[j]]); i++; j++; }
          else if (L[(i + 1) * (m + 1) + j] >= L[i * (m + 1) + j + 1]) out.push(['-', A[i++]]);
          else out.push(['+', B[j++]]);
        }
        while (i < n) out.push(['-', A[i++]]);
        while (j < m) out.push(['+', B[j++]]);
      }
      for (let i = ea; i < a.length; i++) out.push(['=', a[i], b[i - ea + eb]]);
      return out;
    }

    const words = (s) => s.match(/\s+|[^\s]+/g) || [];
    const wordCount = (s) => (s.match(/[^\s]+/g) || []).length;

    // Groups of changed paragraphs: removed ones are paired with added ones for a word diff
    function blocks(was, now) {
      const list = ops(was, now);
      const out = [];
      for (let i = 0; i < list.length;) {
        if (list[i][0] === '=') { out.push({ same: list[i][1] }); i++; continue; }
        const del = [];
        const ins = [];
        while (i < list.length && list[i][0] !== '=') {
          if (list[i][0] === '-') del.push(list[i][1]); else ins.push(list[i][1]);
          i++;
        }
        const k = Math.max(del.length, ins.length);
        for (let p = 0; p < k; p++) out.push({ del: del[p], ins: ins[p] });
      }
      return out;
    }

    function count(was, now) {
      let added = 0;
      let removed = 0;
      for (const b of blocks(was, now)) {
        if (b.same !== undefined) continue;
        if (b.del !== undefined && b.ins !== undefined) {
          for (const [op, x] of ops(words(b.del), words(b.ins))) {
            if (op === '-') removed += wordCount(x);
            if (op === '+') added += wordCount(x);
          }
        } else if (b.del !== undefined) removed += wordCount(b.del);
        else added += wordCount(b.ins);
      }
      return { added, removed };
    }

    function paraHtml(b) {
      if (b.same !== undefined) return `<p>${esc(b.same)}</p>`;
      if (b.del !== undefined && b.ins !== undefined) {
        // ops: ['=', word, word] | ['-', removed word] | ['+', added word]
        return `<p class="pn-diff__changed">${ops(words(b.del), words(b.ins)).map(([op, x]) => {
          if (op === '=') return esc(x);
          if (op === '-') return `<del>${esc(x)}</del>`;
          return `<ins>${esc(x)}</ins>`;
        }).join('')}</p>`;
      }
      if (b.del !== undefined) return `<p class="pn-diff__changed"><del>${esc(b.del)}</del></p>`;
      return `<p class="pn-diff__changed"><ins>${esc(b.ins)}</ins></p>`;
    }

    // whole = false: changed paragraphs with one paragraph of context, the rest folded
    function render(was, now, whole) {
      const list = blocks(was, now);
      if (!list.some((b) => b.same === undefined)) return whole ? list.map(paraHtml).join('') : '';
      if (whole) return list.map(paraHtml).join('');
      const keep = list.map((b, i) => b.same === undefined
        || (list[i - 1] && list[i - 1].same === undefined) || (list[i + 1] && list[i + 1].same === undefined));
      let out = '';
      let skipped = 0;
      list.forEach((b, i) => {
        if (keep[i]) {
          if (skipped) out += `<p class="pn-diff__fold">${esc(t('pn.unchanged', { n: skipped }))}</p>`;
          skipped = 0;
          out += paraHtml(b);
        } else skipped++;
      });
      if (skipped) out += `<p class="pn-diff__fold">${esc(t('pn.unchanged', { n: skipped }))}</p>`;
      return out;
    }

    return { count, render };
  })();

  // ---------- Start ----------

  start().then(() => {
    // signing in or out as the admin changes drafts and buttons: simply reload
    if (window.Auth) {
      Auth.onChange(({ profile }) => {
        const admin = !!profile?.is_admin;
        if (admin !== isAdmin || (admin || !!profile?.is_reviewer) !== isReviewer) location.reload();
      });
    }
  });
})();
