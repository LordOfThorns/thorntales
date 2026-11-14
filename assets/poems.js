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
  const db = window.Auth ? Auth.client : null;

  // Page state
  let filePoems = [];          // poems from the txt file, file order
  const edits = new Map();     // poem id -> row of poem_edits
  const LANGS = ['ru', 'en', 'de'];
  const fileTr = new Map();    // lang -> Map(poem id -> translated poem from <cat>.<lang>.txt)
  let sourceLang = 'ru';       // the language the collection is written in
  let poems = [];              // what is shown: file + edits, sorted
  // "only what is in my language": a per-reader wish, kept in this browser (shared with texts)
  const ONLY_KEY = 'thorntales.onlyMyLang';
  let onlyMine = false;
  try { onlyMine = localStorage.getItem(ONLY_KEY) === '1'; } catch (e) { /* storage blocked */ }
  let isAdmin = false;
  let editing = null;          // poem id being edited, '' = new poem form, null = none
  let activeTag = '';
  let showUntranslatedMark = false;
  let restricted = new Map();  // admin only: tag -> access mode, for tags that are not public
  let catInfo = null;
  let catalogs = [];
  let catIndex = -1;
  let editingCatalog = false;
  const catalogAdmin = document.createElement('div');
  catalogAdmin.className = 'catalog-admin';
  pageSub.after(catalogAdmin);

  try {
    const [raw, cats, editRows, hidden] = await Promise.all([
      fetchText(`${dataDir}${encodeURIComponent(cat)}.txt`, true),
      PoetryData.catalogs(),
      loadEdits(),
      PoetryData.hiddenPoems(),
    ]);

    catalogs = cats;
    catIndex = catalogs.findIndex(x => x.id === cat);
    catInfo = catalogs[catIndex] || null;
    renderHeader();

    // Poems hidden from this reader by tag access are dropped from the file too
    const hiddenIds = hidden.get(cat) || new Set();
    filePoems = parsePoemsTxt(raw).filter(p => !hiddenIds.has(p.id));
    for (const row of editRows) edits.set(row.poem_id, row);

    // Translations from files, for every language but the original (missing files are fine);
    // the editor needs all of them, the page shows the reader's one
    sourceLang = catInfo?.lang || 'ru';
    await Promise.all(LANGS.filter(l => l !== sourceLang).map(async (l) => {
      const trRaw = await fetchText(`${dataDir}${encodeURIComponent(cat)}.${l}.txt`, false);
      const m = new Map();
      if (trRaw) for (const p of parsePoemsTxt(trRaw)) if (p.text.trim()) m.set(p.id, p);
      fileTr.set(l, m);
    }));
    if (lang !== sourceLang) {
      const translatedCount = buildPoems().filter(p => p.tr).length;
      if (translatedCount) {
        pageSub.insertAdjacentHTML('afterend',
          `<p class="poetry-translated">${escapeHtml(t('poems.translatedCount', { n: translatedCount, total: filePoems.length }))}</p>`);
      }
    }
    showUntranslatedMark = lang !== sourceLang;

    render();
    renderPager();

    // Content is loaded after the browser tried to jump to #hash, so jump again
    if (location.hash) {
      const target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
      if (target) target.scrollIntoView();
    }

    // Show the editor once we know the reader is the admin (and hide it after sign-out)
    if (window.Auth) {
      Auth.onChange(async ({ profile }) => {
        const admin = !!profile?.is_admin;
        if (admin === isAdmin) return;
        isAdmin = admin;
        if (!admin) { editing = null; editingCatalog = false; }
        restricted = admin ? await PoetryData.restrictedTags() : new Map();
        renderHeader();
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
      .select('*') // '*' and not a list: works before and after the tr column is added
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
      // a file poem is plain text: as HTML it is the escaped text (line breaks stay)
      return e ? { ...fromRow(e), edited: true } : { ...p, html: escapeHtml(p.text) };
    });
    for (const e of edits.values()) {
      if (!fileIds.has(e.poem_id)) list.push({ ...fromRow(e), added: true });
    }
    // A poem without a date sorts with the date of the poem before it, so it stays
    // next to its neighbour (forum imports are in publishing order). Sort is stable.
    let prev = '';
    const key = new Map(list.map(p => { if (p.date) prev = p.date; return [p, p.date || prev]; }));
    list.sort((a, b) => key.get(a).localeCompare(key.get(b)));
    if (lang !== sourceLang) {
      for (const p of list) {
        const tr = translationOf(p.id, lang);
        if (tr) p.tr = tr;
      }
    }
    return list;
  }

  function fromRow(e) {
    // a poem saved on the site is HTML (bold, italic...) with its line breaks as typed
    return { id: e.poem_id, title: e.title, date: e.date || '', html: e.body, tags: e.tags || [] };
  }

  // A poem in another language: a translation saved on the site wins over the file one
  function translationOf(id, l) {
    const e = edits.get(id)?.tr?.[l];
    if (e && RichEdit.plain(e.body).trim()) return { title: e.title || '', html: e.body };
    const f = fileTr.get(l)?.get(id);
    return f && f.text.trim() ? { title: f.title, html: escapeHtml(f.text) } : null;
  }

  // ---------- Rendering ----------

  // Collection title, description and (for the admin) its editor
  function renderHeader() {
    const pageLabel = field(catInfo, 'label') || humanizeCat(cat);
    pageTitle.textContent = pageLabel;
    document.title = `${pageLabel} · Thorntales`;
    if (window.Layout) Layout.setCrumbs([{ label: pageLabel }]);
    const desc = field(catInfo, 'desc');
    pageSub.innerHTML = desc ? allowOnlyBr(desc) : '';
    pageSub.hidden = editingCatalog;
    pageTitle.hidden = editingCatalog;

    if (!isAdmin || !catInfo) {
      catalogAdmin.innerHTML = '';
      return;
    }
    if (!editingCatalog) {
      catalogAdmin.innerHTML = `
        <button class="poem__edit-btn" type="button" data-cat-act="edit">${escapeHtml(t('catalog.edit'))}</button>
        <a class="poem__edit-btn" href="/tags.html">${escapeHtml(t('catalog.tags'))}</a>`;
      return;
    }
    const img = catInfo.image || '';
    catalogAdmin.innerHTML = `
      <form class="panel poem-edit catalog-edit">
        <div class="catalog-edit__top">
          <div class="catalog-edit__img">
            ${img ? `<img src="${escapeHtml(img)}" alt="">` : `<span>${escapeHtml(t('catalog.noImage'))}</span>`}
          </div>
          <div class="catalog-edit__img-actions">
            <label class="btn">
              <span>${escapeHtml(t('catalog.uploadImage'))}</span>
              <input type="file" name="image" accept="image/jpeg,image/png,image/webp" hidden>
            </label>
            ${img ? `<button class="btn btn--ghost" type="button" data-cat-act="remove-image">${escapeHtml(t('catalog.removeImage'))}</button>` : ''}
            <p class="field__hint">${escapeHtml(t('catalog.imageHint'))}</p>
          </div>
        </div>
        <label class="field">
          <span class="field__label">${escapeHtml(t('catalog.fieldLabel'))}</span>
          <input class="input" name="label" maxlength="200" value="${escapeHtml(catInfo.label || '')}">
        </label>
        <label class="field">
          <span class="field__label">${escapeHtml(t('catalog.fieldDesc'))}</span>
          <textarea class="input" name="descr" rows="5">${escapeHtml(PoetryData.descToText(catInfo.desc))}</textarea>
        </label>
        <p class="poem-edit__error" hidden></p>
        <div class="poem-edit__actions">
          <button class="btn btn--primary" type="submit">${escapeHtml(t('poems.save'))}</button>
          <button class="btn btn--ghost" type="button" data-cat-act="cancel">${escapeHtml(t('poems.cancel'))}</button>
          ${catInfo.edited ? `<span class="poem-edit__spacer"></span>
            <button class="btn btn--ghost" type="button" data-cat-act="revert">${escapeHtml(t('poems.revert'))}</button>` : ''}
        </div>
      </form>`;
  }

  // Picture chosen in the form but not saved yet
  let pendingImage = null; // File, or '' = remove the picture

  catalogAdmin.addEventListener('change', (e) => {
    if (e.target.name !== 'image' || !e.target.files[0]) return;
    pendingImage = e.target.files[0];
    const box = catalogAdmin.querySelector('.catalog-edit__img');
    box.innerHTML = `<img src="${URL.createObjectURL(pendingImage)}" alt="">`;
  });

  catalogAdmin.addEventListener('click', async (e) => {
    const btn = e.target.closest('[data-cat-act]');
    if (!btn) return;
    const act = btn.dataset.catAct;
    if (act === 'edit') {
      editingCatalog = true;
      pendingImage = null;
      renderHeader();
      catalogAdmin.querySelector('input[name="label"]').focus();
    } else if (act === 'cancel') {
      editingCatalog = false;
      renderHeader();
    } else if (act === 'remove-image') {
      pendingImage = '';
      catalogAdmin.querySelector('.catalog-edit__img').innerHTML = `<span>${escapeHtml(t('catalog.noImage'))}</span>`;
      btn.remove();
    } else if (act === 'revert') {
      if (!confirm(t('catalog.revertConfirm'))) return;
      const form = btn.closest('form');
      setBusy(form, true);
      const { error } = await db.from('catalog_edits').delete().eq('catalog', cat);
      if (error) return showError(form, t('poems.saveFailed') + error.message);
      await PoetryData.removeCatalogImages(cat).catch(() => {});
      await reloadCatalog();
    }
  });

  catalogAdmin.addEventListener('submit', async (e) => {
    e.preventDefault();
    const form = e.target;
    const f = form.elements;
    const label = f.label.value.trim();
    const descr = f.descr.value.replace(/\r\n/g, '\n').trim();
    if (!label) return showError(form, t('poems.needTitle'));

    setBusy(form, true);
    try {
      let image_url = catInfo.image || null;
      if (pendingImage) image_url = await PoetryData.uploadCatalogImage(cat, pendingImage);
      else if (pendingImage === '') {
        image_url = null;
        await PoetryData.removeCatalogImages(cat);
      }
      const { error } = await db.from('catalog_edits').upsert({ catalog: cat, label, descr, image_url });
      if (error) throw error;
    } catch (err) {
      return showError(form, t('poems.saveFailed') + (err.message || err));
    }
    await reloadCatalog();
  });

  async function reloadCatalog() {
    catalogs = await PoetryData.catalogs();
    catInfo = catalogs[catIndex] || catalogs.find(x => x.id === cat) || null;
    editingCatalog = false;
    pendingImage = null;
    renderHeader();
    renderPager();
  }

  function render() {
    poems = buildPoems();
    if (activeTag && !poems.some(p => p.tags.includes(activeTag))) activeTag = '';
    let shown = activeTag ? poems.filter(p => p.tags.includes(activeTag)) : poems;
    // "only in my language": for a reader whose language is not the collection's one
    const canFilter = lang !== sourceLang;
    if (canFilter && onlyMine) shown = shown.filter(p => p.tr || (isAdmin && editing === p.id));

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
    if (canFilter) {
      html += `<label class="lang-only"><input type="checkbox" data-act="only-mine"${onlyMine ? ' checked' : ''}>
        <span>${escapeHtml(t('texts.onlyMyLang'))}</span></label>`;
      if (onlyMine && !shown.length) html += `<p class="text-note">${escapeHtml(t('texts.noneInMyLang'))}</p>`;
    }
    if (activeTag) {
      html += `<div class="poem-filter">${escapeHtml(t('poems.tagFilter'))}:
        <span class="poem-tag poem-tag--active">#${escapeHtml(activeTag)}</span>
        <button class="link-btn" type="button" data-act="clear-tag">${escapeHtml(t('poems.clearTag'))}</button></div>`;
    }
    html += shown.map(p => (isAdmin && editing === p.id) ? renderForm(p) : renderPoem(p)).join('');
    poemsHost.innerHTML = html;
    const openForm = poemsHost.querySelector('.poem-edit');
    if (openForm) initForm(openForm, editing === '' ? null : poems.find(x => x.id === editing));

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
      ? `<div class="poem__tags">${p.tags.map(tag => {
          // The admin sees which tags hide the poem from others
          const mode = restricted.get(tag);
          const cls = (tag === activeTag ? ' poem-tag--active' : '') + (mode ? ' poem-tag--restricted' : '');
          const title = mode ? ` title="${escapeHtml(t('tags.mode.' + mode))}"` : '';
          return `<button class="poem-tag${cls}" type="button" data-act="tag" data-tag="${escapeHtml(tag)}"${title}>#${escapeHtml(tag)}</button>`;
        }).join('')}</div>`
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
          <div class="poem__text">${RichEdit.sanitize(shown.html || '')}</div>
        </div>
      </article>
    `;
  }

  // Editor for an existing poem (p) or a new one (p = null): the original and its
  // translations on tabs (RU / EN / DE), each with a title and a text field with bold,
  // italic and an HTML mode. Date and tags are shared. The fields are set up by
  // initForm() once the form is in the page.
  function renderForm(p) {
    let extra = '';
    if (p?.edited) extra = `<button class="btn btn--ghost" type="button" data-act="revert">${escapeHtml(t('poems.revert'))}</button>`;
    if (p?.added) extra = `<button class="btn btn--ghost" type="button" data-act="revert">${escapeHtml(t('poems.delete'))}</button>`;
    return `
      <article class="poem poem--editing"${p ? ` id="${escapeAttr(p.id)}"` : ''}>
        <form class="poem-edit" data-poem="${escapeAttr(p?.id || '')}">
          <div class="poem-edit__row">
            <label class="field poem-edit__date">
              <span class="field__label">${escapeHtml(t('poems.fieldDate'))}</span>
              <input class="input" name="date" maxlength="10" placeholder="2026-11-13" value="${escapeHtml(p?.date || '')}" title="${escapeHtml(t('poems.dateHint'))}">
            </label>
            <label class="field">
              <span class="field__label">${escapeHtml(t('poems.fieldTags'))} <span class="field__hint">· ${escapeHtml(t('poems.tagsHint'))}</span></span>
              <input class="input" name="tags" value="${escapeHtml((p?.tags || []).join(', '))}">
            </label>
          </div>
          <div class="text-edit__tabs" role="tablist"></div>
          <label class="field">
            <span class="field__label">${escapeHtml(t('poems.fieldTitle'))}</span>
            <input class="input" name="title" maxlength="300">
          </label>
          <div class="field">
            <span class="field__label">${escapeHtml(t('poems.fieldText'))}</span>
            <div class="poem-edit__body"></div>
          </div>
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

  // State of the open poem form: the text in every language, and the open tab
  let formState = null;
  function initForm(form, p) {
    const v = {};
    for (const l of LANGS) v[l] = { title: '', html: '' };
    if (p) {
      v[sourceLang] = { title: p.title || '', html: p.html || '' };
      for (const l of LANGS) {
        if (l === sourceLang) continue;
        const tr = translationOf(p.id, l);
        if (tr) v[l] = { title: tr.title || '', html: tr.html || '' };
      }
    }
    formState = { v, tab: sourceLang, field: null, form };
    openTab(sourceLang);
  }
  function saveTab() {
    const s = formState;
    if (!s?.field) return;
    s.v[s.tab] = { title: s.form.elements.title.value.trim(), html: s.field.get() };
  }
  function openTab(l) {
    const s = formState;
    saveTab();
    s.tab = l;
    s.form.elements.title.value = s.v[l].title;
    s.field = RichEdit.field(s.form.querySelector('.poem-edit__body'), s.v[l].html, { lines: true });
    s.form.querySelector('.text-edit__tabs').innerHTML = LANGS.map((x) => {
      const kind = x === sourceLang ? t('texts.tabOriginal', 'оригинал') : t('texts.tabTranslation', 'перевод');
      const filled = s.v[x].title || RichEdit.plain(s.v[x].html).trim() ? ' is-filled' : '';
      return `<button type="button" role="tab" class="text-edit__tab${x === l ? ' is-active' : ''}${filled}" data-tab="${x}" aria-selected="${x === l}">${x.toUpperCase()} <small>${escapeHtml(kind)}</small></button>`;
    }).join('');
  }

  // ---------- Actions ----------

  poemsHost.addEventListener('click', async (e) => {
    const tabBtn = e.target.closest('[data-tab]');
    if (tabBtn && formState) { openTab(tabBtn.dataset.tab); return; }
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
      form?.querySelector('.rich__area')?.focus({ preventScroll: true });
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

  poemsHost.addEventListener('change', (e) => {
    if (!e.target.matches('[data-act="only-mine"]')) return;
    onlyMine = e.target.checked;
    try { localStorage.setItem(ONLY_KEY, onlyMine ? '1' : '0'); } catch (err) { /* storage blocked */ }
    render();
  });

  poemsHost.addEventListener('submit', async (e) => {
    const form = e.target.closest('.poem-edit');
    if (!form) return;
    e.preventDefault();

    // form.elements, not form.title: a form has its own `title` property
    const f = form.elements;
    saveTab();
    const orig = formState.v[sourceLang];
    const title = orig.title;
    const date = f.date.value.trim();
    const body = orig.html.replace(/\s+$/, '');
    const tags = parseTags(f.tags.value);
    // translations: only the languages that have something in them
    const tr = {};
    for (const l of LANGS) {
      if (l === sourceLang) continue;
      const x = formState.v[l];
      if (x.title || RichEdit.plain(x.html).trim()) tr[l] = { title: x.title, body: x.html.replace(/\s+$/, '') };
    }
    if (!title) {
      if (formState.tab !== sourceLang) openTab(sourceLang);
      return showError(form, t('poems.needTitle'));
    }
    if (date && !/^\d{4}(-(0[1-9]|1[0-2])(-(0[1-9]|[12]\d|3[01]))?)?$/.test(date)) return showError(form, t('poems.badDate'));

    // A new poem gets its id once, from the first date and title; it never changes later
    const id = form.dataset.poem || uniqueId(makeId(date, title));

    setBusy(form, true);
    const { data, error } = await db
      .from('poem_edits')
      .upsert({ catalog: cat, poem_id: id, title, date, body, tags, tr })
      .select('*')
      .single();
    if (error) return showError(form, t('poems.saveFailed') + error.message);

    edits.set(id, data);
    editing = null;
    formState = null;
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

  function renderPager(){
    const i = catIndex;
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
        article.querySelector('.poem__text').innerHTML = RichEdit.sanitize(shown.html || '');
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
