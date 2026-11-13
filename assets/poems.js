/* poems.js — one poetry catalog: poems.html?cat=<catalog id>

   Original:     data/poetry/<cat>.txt
   Translation:  data/poetry/<cat>.<lang>.txt  (optional, e.g. salvia.en.txt)
   Poems are matched by their `id:` line, see data/TRANSLATIONS.md.

   Site edits: table public.poem_edits (supabase/schema.sql). A row replaces the poem
   with the same id from the txt file, or adds a new poem if the file has no such id.
   Only the admin (profiles.is_admin) can write there; the editor shows only for them.
*/
(async function () {
  const toc = document.getElementById('toc');
  const poemsHost = document.getElementById('poems');
  const pageTitle = document.getElementById('pageTitle');
  const pageSub = document.getElementById('pageSub');

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const field = (obj, name) => (obj ? (window.I18n ? I18n.field(obj, name) : obj[name]) : undefined);
  const lang = window.I18n ? I18n.lang : 'ru';

  const cat = new URLSearchParams(location.search).get('cat') || 'lyrical';
  const dataDir = './data/poetry/';
  const catalogsUrl = dataDir + 'catalogs.json';
  const db = window.Auth ? Auth.client : null;

  // Page state
  let filePoems = [];          // poems from the txt file, file order
  const edits = new Map();     // poem id -> row of poem_edits
  const trById = new Map();    // poem id -> translated poem
  let poems = [];              // what is shown: file + edits, sorted
  let isAdmin = false;
  let editing = null;          // poem id being edited, '' = new poem form, null = none
  let activeTag = '';
  let showUntranslatedMark = false;

  try {
    const [raw, catsJson, editRows] = await Promise.all([
      fetchText(`${dataDir}${encodeURIComponent(cat)}.txt`, true),
      fetch(catalogsUrl, { cache: 'no-store' }).then(r => {
        if (!r.ok) throw new Error(`${catalogsUrl}: ${r.status}`);
        return r.json();
      }),
      loadEdits(),
    ]);

    const catalogs = catsJson.catalogs || [];
    const catIndex = catalogs.findIndex(x => x.id === cat);
    const catInfo = catalogs[catIndex];

    // Title and subtitle of the catalog (subtitle may contain <br> only)
    const pageLabel = field(catInfo, 'label') || humanizeCat(cat);
    pageTitle.textContent = pageLabel;
    document.title = `${pageLabel} · Thorntales`;
    if (window.Layout) Layout.setCrumbs([{ label: pageLabel }]);
    const desc = field(catInfo, 'desc');
    pageSub.innerHTML = desc ? allowOnlyBr(desc) : '';

    filePoems = parsePoemsTxt(raw);
    for (const row of editRows) edits.set(row.poem_id, row);

    // Translation for the reader's language, if this catalog is written in another one
    const sourceLang = catInfo?.lang || 'ru';
    if (lang !== sourceLang) {
      const trRaw = await fetchText(`${dataDir}${encodeURIComponent(cat)}.${lang}.txt`, false);
      if (trRaw) {
        for (const p of parsePoemsTxt(trRaw)) if (p.text.trim()) trById.set(p.id, p);
      }
      const translatedCount = filePoems.filter(p => trById.has(p.id)).length;
      if (translatedCount) {
        pageSub.insertAdjacentHTML('afterend',
          `<p class="poetry-translated">${escapeHtml(t('poems.translatedCount', { n: translatedCount, total: filePoems.length }))}</p>`);
      }
    }
    showUntranslatedMark = lang !== sourceLang;

    render();
    renderPager(catalogs, catIndex);

    // Content is loaded after the browser tried to jump to #hash, so jump again
    if (location.hash) {
      const target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
      if (target) target.scrollIntoView();
    }

    // Show the editor once we know the reader is the admin (and hide it after sign-out)
    if (window.Auth) {
      Auth.onChange(({ profile }) => {
        const admin = !!profile?.is_admin;
        if (admin === isAdmin) return;
        isAdmin = admin;
        if (!admin) editing = null;
        render();
      });
    }

  } catch (e) {
    poemsHost.innerHTML = `
      <article class="poem">
        <header class="poem__head">
          <h2 class="poem__title">${escapeHtml(t('poems.loadFailed'))}</h2>
          <div class="poem__date">${escapeHtml(String(e))}</div>
        </header>
      </article>
    `;
  }

  // ---------- Data ----------

  // required=false: a missing file (404) returns null instead of throwing
  async function fetchText(url, required) {
    const res = await fetch(url, { cache: 'no-store' });
    if (!res.ok) {
      if (!required && res.status === 404) return null;
      throw new Error(`${url}: ${res.status}`);
    }
    return res.text();
  }

  // Site edits never break the page: without them the file is shown as is
  async function loadEdits() {
    if (!db) return [];
    const { data, error } = await db
      .from('poem_edits')
      .select('poem_id, title, date, body, tags, updated_at')
      .eq('catalog', cat);
    if (error) {
      console.warn('poem edits failed', error);
      return [];
    }
    return data || [];
  }

  // File poems with edits applied, plus poems added on the site, sorted by date
  function buildPoems() {
    const fileIds = new Set(filePoems.map(p => p.id));
    const list = filePoems.map(p => {
      const e = edits.get(p.id);
      return e ? { ...fromRow(e), edited: true } : { ...p };
    });
    for (const e of edits.values()) {
      if (!fileIds.has(e.poem_id)) list.push({ ...fromRow(e), added: true });
    }
    // A poem without a date sorts with the date of the poem before it, so it stays
    // next to its neighbour (forum imports are in publishing order). Sort is stable.
    let prev = '';
    const key = new Map(list.map(p => { if (p.date) prev = p.date; return [p, p.date || prev]; }));
    list.sort((a, b) => key.get(a).localeCompare(key.get(b)));
    for (const p of list) {
      const tr = trById.get(p.id);
      if (tr) p.tr = tr;
    }
    return list;
  }

  function fromRow(e) {
    return { id: e.poem_id, title: e.title, date: e.date || '', text: e.body, tags: e.tags || [] };
  }

  // ---------- Rendering ----------

  function render() {
    poems = buildPoems();
    if (activeTag && !poems.some(p => p.tags.includes(activeTag))) activeTag = '';
    const shown = activeTag ? poems.filter(p => p.tags.includes(activeTag)) : poems;

    // TOC: sidebar on wide screens, collapsible block on narrow ones
    const wasOpen = toc.querySelector('details')?.open;
    toc.innerHTML = `
      <details class="poetry-toc__details">
        <summary class="poetry-toc__summary">${escapeHtml(t('poems.toc'))} <span class="poetry-toc__count">${shown.length}</span></summary>
        <div class="poetry-toc__inner">
          ${shown.map(p => `
            <a class="poetry-toc__link" href="#${escapeAttr(p.id)}">
              <span class="poetry-toc__date">${escapeHtml(p.date || '')}</span>
              <span class="poetry-toc__title">${escapeHtml((p.tr || p).title || '')}</span>
            </a>
          `).join('')}
        </div>
      </details>
    `;
    setupTocMode(toc.querySelector('details'), wasOpen);

    let html = '';
    if (isAdmin) {
      html += editing === ''
        ? renderForm(null)
        : `<div class="poem-admin"><button class="btn" type="button" data-act="new">${escapeHtml(t('poems.add'))}</button></div>`;
    }
    if (activeTag) {
      html += `<div class="poem-filter">${escapeHtml(t('poems.tagFilter'))}:
        <span class="poem-tag poem-tag--active">#${escapeHtml(activeTag)}</span>
        <button class="link-btn" type="button" data-act="clear-tag">${escapeHtml(t('poems.clearTag'))}</button></div>`;
    }
    html += shown.map(p => (isAdmin && editing === p.id) ? renderForm(p) : renderPoem(p)).join('');
    poemsHost.innerHTML = html;

    wireLangToggles();

    // Comment threads: key = poem/<catalog>/<poem id>, shared by all languages,
    // so poem ids must stay stable (editing a poem never changes its id)
    if (window.Comments) {
      const prefix = `poem/${cat}/`;
      Comments.mount(
        shown.filter(p => editing !== p.id)
          .map(p => ({ el: document.getElementById(p.id), thread: prefix + p.id }))
          .filter(x => x.el),
        prefix
      );
    }
  }

  function renderPoem(p) {
    const shown = p.tr || p;
    let langControl = '';
    if (p.tr) {
      langControl = `<button class="poem__lang-toggle" type="button" data-poem="${escapeAttr(p.id)}">${escapeHtml(t('poems.showOriginal'))}</button>`;
    } else if (showUntranslatedMark) {
      langControl = `<span class="poem__untranslated">${escapeHtml(t('poems.notTranslated'))}</span>`;
    }

    let adminControl = '';
    if (isAdmin) {
      const mark = p.added ? t('poems.newMark') : p.edited ? t('poems.editedMark') : '';
      adminControl = `
        ${mark ? `<span class="poem__site-mark">${escapeHtml(mark)}</span>` : ''}
        <button class="poem__edit-btn" type="button" data-act="edit" data-poem="${escapeAttr(p.id)}">${escapeHtml(t('poems.edit'))}</button>`;
    }

    const tags = p.tags.length
      ? `<div class="poem__tags">${p.tags.map(tag =>
          `<button class="poem-tag${tag === activeTag ? ' poem-tag--active' : ''}" type="button" data-act="tag" data-tag="${escapeHtml(tag)}">#${escapeHtml(tag)}</button>`
        ).join('')}</div>`
      : '';

    return `
      <article class="poem" id="${escapeAttr(p.id)}"${p.tr ? ` lang="${lang}"` : ''}>
        <header class="poem__head">
          <h2 class="poem__title">${escapeHtml(shown.title || '')}</h2>
          <div class="poem__meta">
            ${adminControl}
            ${langControl}
            <span class="poem__date">${escapeHtml(p.date || '')}</span>
          </div>
        </header>
        ${tags}
        <div class="poem__body">
          <pre class="poem__text">${escapeHtml(shown.text || '')}</pre>
        </div>
      </article>
    `;
  }

  // Editor for an existing poem (p) or a new one (p = null). Always edits the original.
  function renderForm(p) {
    const rows = Math.min(40, Math.max(8, (p?.text || '').split('\n').length + 2));
    let extra = '';
    if (p?.edited) extra = `<button class="btn btn--ghost" type="button" data-act="revert">${escapeHtml(t('poems.revert'))}</button>`;
    if (p?.added) extra = `<button class="btn btn--ghost" type="button" data-act="revert">${escapeHtml(t('poems.delete'))}</button>`;
    return `
      <article class="poem poem--editing"${p ? ` id="${escapeAttr(p.id)}"` : ''}>
        <form class="poem-edit" data-poem="${escapeAttr(p?.id || '')}">
          <div class="poem-edit__row">
            <label class="field poem-edit__title">
              <span class="field__label">${escapeHtml(t('poems.fieldTitle'))}</span>
              <input class="input" name="title" maxlength="300" value="${escapeHtml(p?.title || '')}">
            </label>
            <label class="field poem-edit__date">
              <span class="field__label">${escapeHtml(t('poems.fieldDate'))}</span>
              <input class="input" name="date" maxlength="10" placeholder="2026-11-13" value="${escapeHtml(p?.date || '')}" title="${escapeHtml(t('poems.dateHint'))}">
            </label>
          </div>
          <label class="field">
            <span class="field__label">${escapeHtml(t('poems.fieldTags'))} <span class="field__hint">· ${escapeHtml(t('poems.tagsHint'))}</span></span>
            <input class="input" name="tags" value="${escapeHtml((p?.tags || []).join(', '))}">
          </label>
          <label class="field">
            <span class="field__label">${escapeHtml(t('poems.fieldText'))}</span>
            <textarea class="input poem-edit__text" name="body" rows="${rows}">${escapeHtml(p?.text || '')}</textarea>
          </label>
          <p class="poem-edit__error" hidden></p>
          <div class="poem-edit__actions">
            <button class="btn btn--primary" type="submit">${escapeHtml(t('poems.save'))}</button>
            <button class="btn btn--ghost" type="button" data-act="cancel">${escapeHtml(t('poems.cancel'))}</button>
            ${extra ? `<span class="poem-edit__spacer"></span>${extra}` : ''}
          </div>
        </form>
      </article>
    `;
  }

  // ---------- Actions ----------

  poemsHost.addEventListener('click', async (e) => {
    const btn = e.target.closest('[data-act]');
    if (!btn) return;
    const act = btn.dataset.act;

    if (act === 'tag') {
      activeTag = btn.dataset.tag === activeTag ? '' : btn.dataset.tag;
      render();
      poemsHost.scrollIntoView({ block: 'start' });
    } else if (act === 'clear-tag') {
      activeTag = '';
      render();
    } else if (act === 'new') {
      editing = '';
      render();
      poemsHost.querySelector('.poem-edit input[name="title"]')?.focus();
    } else if (act === 'edit') {
      editing = btn.dataset.poem;
      render();
      const form = poemsHost.querySelector('.poem-edit');
      form?.closest('.poem').scrollIntoView({ block: 'start' });
      form?.querySelector('textarea').focus({ preventScroll: true });
    } else if (act === 'cancel') {
      const id = editing;
      editing = null;
      render();
      if (id) document.getElementById(id)?.scrollIntoView({ block: 'start' });
    } else if (act === 'revert') {
      const id = editing;
      const p = poems.find(x => x.id === id);
      if (!p || !confirm(t(p.added ? 'poems.deleteConfirm' : 'poems.revertConfirm'))) return;
      const form = btn.closest('form');
      setBusy(form, true);
      const { error } = await db.from('poem_edits').delete().eq('catalog', cat).eq('poem_id', id);
      if (error) return showError(form, t('poems.saveFailed') + error.message);
      edits.delete(id);
      editing = null;
      render();
      document.getElementById(id)?.scrollIntoView({ block: 'start' });
    }
  });

  poemsHost.addEventListener('submit', async (e) => {
    const form = e.target.closest('.poem-edit');
    if (!form) return;
    e.preventDefault();

    // form.elements, not form.title: a form has its own `title` property
    const f = form.elements;
    const title = f.title.value.trim();
    const date = f.date.value.trim();
    const body = f.body.value.replace(/\r\n/g, '\n').replace(/\s+$/, '');
    const tags = parseTags(f.tags.value);
    if (!title) return showError(form, t('poems.needTitle'));
    if (date && !/^\d{4}(-(0[1-9]|1[0-2])(-(0[1-9]|[12]\d|3[01]))?)?$/.test(date)) return showError(form, t('poems.badDate'));

    // A new poem gets its id once, from the first date and title; it never changes later
    const id = form.dataset.poem || uniqueId(makeId(date, title));

    setBusy(form, true);
    const { data, error } = await db
      .from('poem_edits')
      .upsert({ catalog: cat, poem_id: id, title, date, body, tags })
      .select('poem_id, title, date, body, tags, updated_at')
      .single();
    if (error) return showError(form, t('poems.saveFailed') + error.message);

    edits.set(id, data);
    editing = null;
    render();
    document.getElementById(id)?.scrollIntoView({ block: 'start' });
  });

  function setBusy(form, busy) {
    form.querySelectorAll('button, input, textarea').forEach(el => { el.disabled = busy; });
  }

  function showError(form, msg) {
    setBusy(form, false);
    const box = form.querySelector('.poem-edit__error');
    box.textContent = msg;
    box.hidden = false;
  }

  function parseTags(s) {
    const out = [];
    for (let tag of s.split(/[,;\n]/)) {
      tag = tag.trim().replace(/^#+/, '').replace(/\s+/g, ' ').slice(0, 40);
      if (tag && !out.includes(tag)) out.push(tag);
    }
    return out.slice(0, 30);
  }

  function uniqueId(base) {
    const taken = new Set([...filePoems.map(p => p.id), ...edits.keys()]);
    base = base || 'poem';
    let id = base, n = 2;
    while (taken.has(id)) id = `${base.slice(0, 76)}-${n++}`;
    return id;
  }

  // ---------- Page parts ----------

  function setupTocMode(details, wasOpen){
    if (!details) return;
    const wide = window.matchMedia('(min-width: 1000px)');
    const apply = () => {
      details.open = wide.matches;
      details.classList.toggle('is-sidebar', wide.matches);
    };
    apply();
    if (!wide.matches && wasOpen) details.open = true;
    wide.addEventListener('change', apply);

    // On narrow screens close the TOC after jumping to a poem
    details.addEventListener('click', (e) => {
      if (!wide.matches && e.target.closest('.poetry-toc__link')) details.open = false;
    });
  }

  function renderPager(catalogs, i){
    const pager = document.getElementById('pager');
    if (!pager || i < 0) return;
    const prev = catalogs[i - 1];
    const next = catalogs[i + 1];
    pager.innerHTML = `
      ${prev ? `<a class="catalog-pager__link catalog-pager__link--prev" href="${escapeAttr(prev.href)}">
        <span class="catalog-pager__hint">${escapeHtml(t('poems.prev'))}</span>
        <span class="catalog-pager__label">${escapeHtml(field(prev, 'label'))}</span></a>` : '<span></span>'}
      <a class="catalog-pager__link catalog-pager__link--all" href="/poetry.html">${escapeHtml(t('poems.all'))}</a>
      ${next ? `<a class="catalog-pager__link catalog-pager__link--next" href="${escapeAttr(next.href)}">
        <span class="catalog-pager__hint">${escapeHtml(t('poems.next'))}</span>
        <span class="catalog-pager__label">${escapeHtml(field(next, 'label'))}</span></a>` : '<span></span>'}
    `;
  }

  function allowOnlyBr(html){
    // Escape everything, then bring back only <br>
    const safe = escapeHtml(html);
    return safe.replace(/&lt;br\s*\/?&gt;/gi, '<br>');
  }

  // Poems in file order. Header lines: date, title, id, tags (comma separated)
  function parsePoemsTxt(raw) {
    const text = raw.replace(/\r\n/g, '\n');

    // Split into blocks by the marker; anything before the first marker is ignored
    let blocks = text.split('\n=== poem ===\n');
    if (text.startsWith('=== poem ===')) {
      blocks = text.split('=== poem ===\n').slice(1);
    }

    const list = [];

    for (const blockRaw of blocks) {
      const block = blockRaw.trim();
      if (!block) continue;

      // Header lines ("key: value") go until the first empty line
      const parts = block.split('\n\n');
      const headerPart = parts[0] || '';
      const bodyPart = parts.slice(1).join('\n\n') || '';

      const meta = {};
      for (const line of headerPart.split('\n')) {
        const m = line.match(/^([a-zA-Z0-9_-]+)\s*:\s*(.*)$/);
        if (!m) continue;
        meta[m[1].toLowerCase()] = m[2].trim();
      }

      const title = meta.title || t('poems.untitled');
      const date = meta.date || '';
      const id = meta.id || makeId(date, title);
      const tags = meta.tags ? parseTags(meta.tags) : [];

      list.push({ id, date, title, text: bodyPart, tags });
    }

    return list;
  }

  // "Original" / "Translation" switch on translated poems
  function wireLangToggles(){
    const byId = new Map(poems.map(p => [p.id, p]));
    poemsHost.querySelectorAll('.poem__lang-toggle').forEach(btn => {
      let showingOriginal = false;
      btn.addEventListener('click', () => {
        const p = byId.get(btn.closest('.poem').id);
        if (!p) return;
        showingOriginal = !showingOriginal;
        const shown = showingOriginal ? p : p.tr;
        const article = btn.closest('.poem');
        article.querySelector('.poem__title').textContent = shown.title || '';
        article.querySelector('.poem__text').textContent = shown.text || '';
        if (showingOriginal) article.removeAttribute('lang'); else article.setAttribute('lang', lang);
        btn.textContent = t(showingOriginal ? 'poems.showTranslation' : 'poems.showOriginal');
      });
    });
  }

  function makeId(date, title){
    const base = `${date}-${title}`.trim() || title || 'poem';
    return base
      .toLowerCase()
      .replace(/ё/g, 'е')
      .replace(/[^a-z0-9Ѐ-ӿ]+/gi, '-')
      .replace(/^-+|-+$/g, '')
      .slice(0, 80);
  }

  function humanizeCat(cat){
    return cat === 'lyrical' ? t('poems.lyrical') : t('poems.catalog', { cat });
  }

  function escapeHtml(s){
    s = String(s ?? '');
    return s.replace(/[&<>"']/g, ch => ({
      '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'
    }[ch]));
  }
  function escapeAttr(s){
    return escapeHtml(String(s ?? '')).replace(/\s+/g, '-');
  }
})();
