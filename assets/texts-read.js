/* texts-read.js — one book of texts: texts.html?s=<section> (the list) and
   texts.html?s=<section>&t=<text id> (one text with comments).

   Texts:  data/texts/<section>.json  { "texts": [ { id, title, date, lang, html,
             i18n: { en: { title, html }, ... } } ] }   (a missing file = an empty book)
   Edits:  table public.text_edits (supabase/schema.sql): a row replaces the text with
           the same id from the file, adds a new one, or hides one (deleted).
   Comments: thread "text/<section>/<text id>", shared by all languages.

   The admin (profiles.is_admin) sees "New text" and "Edit": an editor with bold,
   italic and an HTML mode, and tabs for the original and its translations.
   A translation is shown to readers whose site language it is.

   Used by texts.js: TextBooks.show(section, textId).
*/
(function () {
  'use strict';

  const t = (key, fallback, vars) => {
    if (window.I18n && I18n.has(key)) return I18n.t(key, vars);
    return vars ? fallback.replace(/\{(\w+)\}/g, (_, k) => vars[k] ?? '') : fallback;
  };
  const siteLang = window.I18n ? I18n.lang : 'ru';
  const db = window.Auth ? Auth.client : null;
  const LANGS = ['ru', 'en', 'de'];

  // HTML of a text: only a small, safe set of tags (rich-editor.js)
  const { sanitize: rawSanitize, plain, esc } = window.RichEdit;
  const sanitize = (html) => rawSanitize(html).trim();

  function makeId(date, title) {
    return `${date}-${title}`.toLowerCase().replace(/ё/g, 'е')
      .replace(/[^a-z0-9Ѐ-ӿ]+/gi, '-').replace(/^-+|-+$/g, '').slice(0, 80) || 'text';
  }

  // ---------- Data ----------

  async function load(section) {
    const [file, rows] = await Promise.all([
      fetch(`./data/texts/${encodeURIComponent(section)}.json`, { cache: 'no-store' })
        .then((r) => (r.ok ? r.json() : { texts: [] })).catch(() => ({ texts: [] })),
      db ? db.from('text_edits').select('*').eq('section', section)
        .then(({ data, error }) => { if (error) console.warn('text edits failed', error); return data || []; })
        : Promise.resolve([]),
    ]);
    const fromFile = (x) => ({
      id: x.id, title: x.title, date: x.date || '', lang: x.lang || 'ru', body: x.html || '',
      tr: Object.fromEntries(Object.entries(x.i18n || {}).map(([l, v]) => [l, { title: v.title || '', body: v.html || '' }])),
      inFile: true,
    });
    const fromRow = (r) => ({
      id: r.text_id, title: r.title, date: r.date || '', lang: r.lang || 'ru', body: r.body || '',
      tr: r.tr || {}, edited: true, deleted: r.deleted,
    });
    const byId = new Map();
    for (const x of file.texts || []) byId.set(x.id, fromFile(x));
    for (const r of rows) {
      const had = byId.get(r.text_id);
      byId.set(r.text_id, { ...fromRow(r), inFile: !!had?.inFile, added: !had });
    }
    const all = [...byId.values()];
    all.sort((a, b) => (a.date || '').localeCompare(b.date || '') || a.title.localeCompare(b.title));
    return all;
  }

  // Languages a text exists in: its original first, then its translations
  function langsOf(x) {
    const out = [x.lang];
    for (const l of LANGS) {
      if (l !== x.lang && x.tr?.[l] && plain(x.tr[l].body).trim()) out.push(l);
    }
    return out;
  }
  // The language a reader sees first: the site language if the text has it, else the original
  const defaultLang = (x) => (langsOf(x).includes(siteLang) ? siteLang : x.lang);

  // The text in one of its languages
  function version(x, l = defaultLang(x)) {
    if (l !== x.lang && x.tr?.[l]) return { title: x.tr[l].title || x.title, body: x.tr[l].body, lang: l };
    return { title: x.title, body: x.body, lang: x.lang };
  }

  // ---------- Page ----------

  let root, sec, texts = [], isAdmin = false;

  async function show(section, textId) {
    sec = section;
    root = document.getElementById('desk');
    root.className = 'text-book';
    root.removeAttribute('style');
    root.innerHTML = `<p class="sub">${esc(t('account.loading', 'Загрузка…'))}</p>`;
    document.title = `${section.label} · Thorntales`;
    texts = await load(section.id);
    isAdmin = !!(window.Auth && Auth.profile?.is_admin);
    route(textId);
    if (window.Auth) {
      Auth.onChange(({ profile }) => {
        const admin = !!profile?.is_admin;
        if (admin === isAdmin) return;
        isAdmin = admin;
        route(new URLSearchParams(location.search).get('t'));
      });
    }
  }

  function route(textId) {
    if (textId === 'new') return isAdmin ? editor(null) : list();
    if (textId) {
      const x = texts.find((y) => y.id === textId && !y.deleted);
      return x ? reader(x) : list(t('texts.notFound', 'Такого текста нет.'));
    }
    list();
  }

  const urlFor = (id) => `/texts.html?s=${encodeURIComponent(sec.id)}${id ? `&t=${encodeURIComponent(id)}` : ''}`;

  function header(title, sub) {
    const h1 = document.querySelector('.top h1');
    const p = document.querySelector('.top .sub');
    if (h1) h1.textContent = title;
    if (p) p.innerHTML = sub || '';
  }

  // ---------- The list of texts in a book ----------

  function list(note) {
    header(sec.label, '');
    if (window.Layout) Layout.setCrumbs([{ label: sec.label }]);
    document.title = `${sec.label} · Thorntales`;
    const all = texts.filter((x) => !x.deleted || isAdmin);
    // "only in my language": shown when some texts are not in the site language
    const someMissing = all.some((x) => !langsOf(x).includes(siteLang));
    const only = onlyMine();
    const visible = only ? all.filter((x) => langsOf(x).includes(siteLang)) : all;
    const items = visible.map((x) => {
      const v = version(x);
      // not in the site language: show every language it is there in
      // A quiet column of the same width in every row, filled only when the text is not
      // in the site language: the original first (warm, with "orig."), translations after.
      const langs = langsOf(x);
      const langInfo = langs.includes(siteLang) ? '' : langs.map((l) => (l === x.lang
        ? `<span class="text-list__lang text-list__lang--orig" title="${esc(t('texts.originalLang', 'язык оригинала'))}">${esc(l.toUpperCase())}<small>${esc(t('texts.origShort', 'ориг.'))}</small></span>`
        : `<span class="text-list__lang" title="${esc(t('texts.translationLang', 'перевод'))}">${esc(l.toUpperCase())}</span>`))
        .join('<span class="text-list__sep">·</span>');
      const mark = isAdmin && x.deleted ? `<span class="text-list__lang">${esc(t('texts.hiddenMark', 'скрыт'))}</span>` : '';
      return `
        <li class="text-list__item${x.deleted ? ' is-deleted' : ''}">
          <a href="${esc(urlFor(x.id))}">
            <span class="text-list__date">${esc(x.date)}</span>
            <span class="text-list__title">${esc(v.title)}</span>
            <span class="text-list__langs">${langInfo}${mark}</span>
          </a>
        </li>`;
    }).join('');
    root.innerHTML = `
      ${note ? `<p class="text-note">${esc(note)}</p>` : ''}
      ${isAdmin ? `<div class="text-admin"><a class="btn" href="${esc(urlFor('new'))}">${esc(t('texts.newText', '+ Новый текст'))}</a></div>` : ''}
      ${someMissing || only ? langOnlySwitch(only) : ''}
      ${visible.length ? `<ol class="text-list">${items}</ol>`
        : `<p class="text-note">${esc(only && all.length ? t('texts.noneInMyLang', 'На вашем языке здесь пока ничего нет.') : t('texts.empty', 'Эта книга пока пустая.'))}</p>`}
      <p class="texts-back"><a href="/texts.html">${esc(t('texts.back', '← к стопкам книг'))}</a></p>`;
  }

  // "Only what is in my language": a per-reader wish kept in this browser, the same
  // switch as on the poem pages (poems.js reads the same key)
  const ONLY_KEY = 'thorntales.onlyMyLang';
  function onlyMine() {
    try { return localStorage.getItem(ONLY_KEY) === '1'; } catch (e) { return false; }
  }
  function langOnlySwitch(on) {
    return `<label class="lang-only"><input type="checkbox" data-act="only-mine"${on ? ' checked' : ''}>
      <span>${esc(t('texts.onlyMyLang', 'Только на моём языке'))}</span></label>`;
  }
  document.addEventListener('change', (e) => {
    if (!e.target.matches('#desk [data-act="only-mine"]')) return;
    try { localStorage.setItem(ONLY_KEY, e.target.checked ? '1' : '0'); } catch (err) { /* storage blocked */ }
    list();
  });

  // ---------- One text ----------

  function reader(x) {
    let shown = defaultLang(x);
    const draw = () => {
      const v = version(x, shown);
      const i = texts.filter((y) => !y.deleted).indexOf(x);
      const vis = texts.filter((y) => !y.deleted);
      const prev = vis[i - 1], next = vis[i + 1];
      document.title = `${v.title} · Thorntales`;
      if (window.Layout) Layout.setCrumbs([{ label: sec.label, href: urlFor() }, { label: v.title }]);
      // one button per language the text exists in (the original is marked)
      const langs = langsOf(x);
      let langControl = '';
      if (langs.length > 1) {
        langControl = `<span class="text-langs" role="group">${langs.map((l) => `
          <button class="text-langs__btn" type="button" data-act="lang" data-lang="${l}" aria-pressed="${l === shown}"
            title="${esc(l === x.lang ? t('texts.originalLang', 'язык оригинала') : t('texts.translationLang', 'перевод'))}">${esc(l.toUpperCase())}${l === x.lang ? ` <span class="text-langs__orig">· ${esc(t('texts.original', 'оригинал'))}</span>` : ''}</button>`).join('')}</span>`;
      }
      if (!langs.includes(siteLang)) {
        langControl += ` <span class="poem__untranslated">${esc(t('poems.notTranslated', 'ещё не переведено'))}</span>`;
      }
      header(v.title, `<span class="text-meta">${esc(x.date)}</span> ${langControl}
        ${isAdmin ? `<button class="poem__edit-btn" type="button" data-act="edit">${esc(t('texts.edit', 'Править'))}</button>` : ''}`);
      root.innerHTML = `
        <article class="text-body" lang="${v.lang}">${sanitize(v.body)}</article>
        <div class="text-comments" id="textComments"></div>
        <nav class="catalog-pager">
          ${prev ? `<a class="catalog-pager__link catalog-pager__link--prev" href="${esc(urlFor(prev.id))}">
            <span class="catalog-pager__hint">${esc(t('texts.prev', '← предыдущий'))}</span>
            <span class="catalog-pager__label">${esc(version(prev).title)}</span></a>` : '<span></span>'}
          <a class="catalog-pager__link catalog-pager__link--all" href="${esc(urlFor())}">${esc(t('texts.toList', 'Все тексты книги'))}</a>
          ${next ? `<a class="catalog-pager__link catalog-pager__link--next" href="${esc(urlFor(next.id))}">
            <span class="catalog-pager__hint">${esc(t('texts.next', 'следующий →'))}</span>
            <span class="catalog-pager__label">${esc(version(next).title)}</span></a>` : '<span></span>'}
        </nav>`;
      if (window.Comments) {
        const prefix = `text/${sec.id}/`;
        Comments.mount([{ el: document.getElementById('textComments'), thread: prefix + x.id }], prefix);
      }
    };
    draw();
    const top = document.querySelector('.top');
    top.onclick = (e) => {
      const act = e.target.closest('[data-act]')?.dataset.act;
      if (act === 'lang') { shown = e.target.closest('[data-lang]').dataset.lang; draw(); }
      if (act === 'edit') editor(x);
    };
  }

  // ---------- Editor (admin) ----------

  const richField = (container, html) => RichEdit.field(container, html);

  function editor(x) {
    const isNew = !x;
    const st = {
      lang: x?.lang || 'ru',
      date: x?.date || new Date().toISOString().slice(0, 10),
      v: {}, // lang -> { title, body }
    };
    for (const l of LANGS) st.v[l] = { title: '', body: '' };
    if (x) {
      st.v[x.lang] = { title: x.title, body: x.body };
      for (const [l, v] of Object.entries(x.tr || {})) if (l !== x.lang) st.v[l] = { title: v.title || '', body: v.body || '' };
    }
    let tab = st.lang;
    let field = null;

    header(isNew ? t('texts.newText', '+ Новый текст').replace(/^\+\s*/, '') : x.title, '');
    if (window.Layout) Layout.setCrumbs([{ label: sec.label, href: urlFor() }, { label: isNew ? t('texts.newTextCrumb', 'Новый текст') : x.title }]);
    root.innerHTML = `
      <form class="panel text-edit">
        <div class="text-edit__row">
          <label class="field">
            <span class="field__label">${esc(t('texts.origLang', 'Язык оригинала'))}</span>
            <select class="input" name="lang">${LANGS.map((l) => `<option value="${l}"${l === st.lang ? ' selected' : ''}>${l.toUpperCase()}</option>`).join('')}</select>
          </label>
          <label class="field">
            <span class="field__label">${esc(t('poems.fieldDate', 'Дата'))}</span>
            <input class="input" name="date" maxlength="10" value="${esc(st.date)}" placeholder="2026-11-14">
          </label>
        </div>
        <div class="text-edit__tabs" role="tablist"></div>
        <label class="field">
          <span class="field__label">${esc(t('poems.fieldTitle', 'Название'))}</span>
          <input class="input" name="title" maxlength="300">
        </label>
        <div class="field">
          <span class="field__label">${esc(t('poems.fieldText', 'Текст'))}</span>
          <div class="text-edit__body"></div>
        </div>
        <p class="poem-edit__error" hidden></p>
        <div class="poem-edit__actions">
          <button class="btn btn--primary" type="submit">${esc(t('poems.save', 'Сохранить'))}</button>
          <a class="btn btn--ghost" href="${esc(isNew ? urlFor() : urlFor(x.id))}">${esc(t('poems.cancel', 'Отмена'))}</a>
          <span class="poem-edit__spacer"></span>
          ${x?.edited && x.inFile ? `<button class="btn btn--ghost" type="button" data-act="revert">${esc(t('poems.revert', 'Вернуть из файла'))}</button>` : ''}
          ${x ? `<button class="btn btn--ghost" type="button" data-act="delete">${esc(t('texts.delete', 'Удалить текст'))}</button>` : ''}
        </div>
      </form>`;
    const form = root.querySelector('form');
    const f = form.elements;

    const saveTab = () => {
      st.v[tab] = { title: f.title.value.trim(), body: field.get() };
    };
    const drawTabs = () => {
      form.querySelector('.text-edit__tabs').innerHTML = LANGS.map((l) => {
        const kind = l === st.lang ? t('texts.tabOriginal', 'оригинал') : t('texts.tabTranslation', 'перевод');
        const filled = st.v[l].title || plain(st.v[l].body).trim() ? ' is-filled' : '';
        return `<button type="button" role="tab" class="text-edit__tab${l === tab ? ' is-active' : ''}${filled}" data-tab="${l}" aria-selected="${l === tab}">${l.toUpperCase()} <small>${esc(kind)}</small></button>`;
      }).join('');
    };
    const openTab = (l) => {
      if (field) saveTab();
      tab = l;
      f.title.value = st.v[l].title;
      field = richField(form.querySelector('.text-edit__body'), st.v[l].body);
      drawTabs();
    };
    openTab(tab);

    form.addEventListener('click', async (e) => {
      const tb = e.target.closest('[data-tab]');
      if (tb) { openTab(tb.dataset.tab); return; }
      const act = e.target.closest('[data-act]')?.dataset.act;
      if (act === 'delete' || act === 'revert') {
        const msg = act === 'revert' ? t('poems.revertConfirm', 'Убрать правки с сайта и вернуть текст из файла?') : t('texts.deleteConfirm', 'Удалить этот текст?');
        if (!confirm(msg)) return;
        busy(true);
        let error;
        if (act === 'revert' || !x.inFile) {
          ({ error } = await db.from('text_edits').delete().eq('section', sec.id).eq('text_id', x.id));
        } else {
          // a text from the file can only be hidden: the row says "deleted"
          ({ error } = await db.from('text_edits').upsert(rowOf(x.id, true)));
        }
        if (error) return fail(t('poems.saveFailed', 'Не сохранилось: ') + error.message);
        texts = await load(sec.id);
        history.replaceState(null, '', urlFor(act === 'revert' ? x.id : null));
        route(act === 'revert' ? x.id : null);
      }
    });
    f.lang.addEventListener('change', () => { st.lang = f.lang.value; drawTabs(); });

    function rowOf(id, deleted = false) {
      const tr = {};
      for (const l of LANGS) {
        if (l === st.lang) continue;
        const v = st.v[l];
        if (v.title || plain(v.body).trim()) tr[l] = { title: v.title, body: v.body };
      }
      return {
        section: sec.id, text_id: id, title: st.v[st.lang].title || x?.title || id,
        date: st.date, lang: st.lang, body: st.v[st.lang].body, tr, deleted,
      };
    }
    function busy(on) { form.querySelectorAll('button, input, select, textarea').forEach((el) => { el.disabled = on; }); }
    function fail(msg) {
      busy(false);
      const box = form.querySelector('.poem-edit__error');
      box.textContent = msg;
      box.hidden = false;
    }

    form.addEventListener('submit', async (e) => {
      e.preventDefault();
      saveTab();
      st.date = f.date.value.trim();
      st.lang = f.lang.value;
      if (!st.v[st.lang].title) return fail(t('texts.needTitle', 'Нужно название на языке оригинала'));
      if (st.date && !/^\d{4}(-(0[1-9]|1[0-2])(-(0[1-9]|[12]\d|3[01]))?)?$/.test(st.date)) {
        return fail(t('poems.badDate', 'Дата должна быть вида 2024-05-17, 2024-05 или 2024'));
      }
      let id = x?.id;
      if (!id) {
        const taken = new Set(texts.map((y) => y.id));
        const base = makeId(st.date, st.v[st.lang].title);
        id = base;
        for (let n = 2; taken.has(id); n++) id = `${base.slice(0, 76)}-${n}`;
      }
      busy(true);
      const { error } = await db.from('text_edits').upsert(rowOf(id));
      if (error) return fail(t('poems.saveFailed', 'Не сохранилось: ') + error.message);
      texts = await load(sec.id);
      history.replaceState(null, '', urlFor(id));
      route(id);
      window.scrollTo(0, 0);
    });
  }

  window.TextBooks = { show, sanitize };
})();
