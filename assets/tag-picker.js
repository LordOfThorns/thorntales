/* tag-picker.js — clickable tags under a "tags, comma separated" field (admin editors of
   poems and texts).

   TagPicker.attach(input) puts a row of chips under the input: every tag used on the site
   (poems in the txt files and edited on the site, texts in the json files and edited on
   the site, tags with an access rule), the most used first. A click adds the tag to the
   input or takes it out; the input stays editable, so a new tag is still typed by hand.
   Chips of tags already in the field are lit; typing in the field updates them.
   A long list is folded to the first FOLD tags, with a filter by the word being typed.

   The list of tags is loaded once per page (TagPicker.loadAll), on the first editor.
*/
(function () {
  'use strict';

  const FOLD = 24;
  const SECTIONS = ['stories', 'tales', 'repatriate', 'dreams', 'chatter', 'other'];
  const t = (key, fallback) => (window.I18n && I18n.has(key) ? I18n.t(key) : fallback);

  const split = (s) => String(s || '').split(/[,;\n]/)
    .map((x) => x.trim().replace(/^#+/, '').replace(/\s+/g, ' ').slice(0, 40))
    .filter(Boolean);

  const get = (url) => (window.SiteFetch ? SiteFetch.get(url, { optional: true }) : fetch(url).then((r) => (r.ok ? r : null)));

  let loading = null;
  // tag -> how many poems and texts carry it
  function loadAll() {
    if (loading) return loading;
    loading = (async () => {
      const counts = new Map();
      const add = (tag) => counts.set(tag, (counts.get(tag) || 0) + 1);
      const jobs = [];
      // poems in the txt files of every collection
      jobs.push((async () => {
        const res = await get('/data/poetry/catalogs.json');
        const list = res ? (await res.json()).catalogs || [] : [];
        await Promise.all(list.map(async (c) => {
          const r = await get(`/data/poetry/${encodeURIComponent(c.id)}.txt`);
          if (!r) return;
          for (const m of (await r.text()).matchAll(/^tags:\s*(.+)$/gm)) split(m[1]).forEach(add);
        }));
      })());
      // texts in the json files of the sections
      jobs.push(...SECTIONS.map(async (s) => {
        const r = await get(`/data/texts/${s}.json`);
        if (!r) return;
        for (const x of (await r.json()).texts || []) (x.tags || []).forEach(add);
      }));
      // edits made on the site and tags with an access rule (the admin may read them)
      const db = window.Auth && Auth.client;
      if (db) {
        const rows = async (table, cols) => {
          try { const { data } = await db.from(table).select(cols); return data || []; } catch (e) { return []; }
        };
        jobs.push(rows('poem_edits', 'tags').then((d) => d.forEach((r) => (r.tags || []).forEach(add))));
        jobs.push(rows('text_edits', 'tags').then((d) => d.forEach((r) => (r.tags || []).forEach(add))));
        jobs.push(rows('tag_access', 'tag').then((d) => d.forEach((r) => { if (!counts.has(r.tag)) counts.set(r.tag, 0); })));
      }
      await Promise.all(jobs.map((j) => j.catch((e) => console.warn('tags', e))));
      return counts;
    })();
    return loading;
  }

  function attach(input) {
    if (!input || input.dataset.tagPicker) return;
    input.dataset.tagPicker = '1';
    const box = document.createElement('div');
    box.className = 'tag-picker';
    box.innerHTML = `<span class="tag-picker__loading">${esc(t('tags.pickLoading', 'Загружаю теги…'))}</span>`;
    input.insertAdjacentElement('afterend', box);

    let all = [];        // tags sorted by use
    let open = false;    // the whole list or only the first FOLD

    const current = () => split(input.value);
    function setTags(list) {
      input.value = list.join(', ');
      input.dispatchEvent(new Event('input', { bubbles: true }));
    }

    function render() {
      const have = current();
      // the word being typed after the last comma filters the list
      const typing = /,\s*$/.test(input.value) ? '' : (input.value.split(/[,;]/).pop() || '').trim().toLowerCase();
      // tags in the field come first (new ones too), then the rest by use
      let rest = all.filter((x) => !have.includes(x));
      // a word that is not a known tag yet filters the list, while some tag still matches
      // it (no match: it is a new tag, and the whole list stays)
      const matches = typing && !all.includes(typing) ? rest.filter((x) => x.toLowerCase().includes(typing)) : [];
      const filtering = matches.length > 0;
      if (filtering) rest = matches;
      const shown = open || filtering ? rest : rest.slice(0, FOLD);
      const chip = (tag, on) => `<button type="button" class="tag-chip${on ? ' is-on' : ''}" aria-pressed="${on}" data-tag="${esc(tag)}">${on ? '✓ ' : '+ '}${esc(tag)}</button>`;
      const more = !open && !filtering && rest.length > FOLD
        ? `<button type="button" class="tag-picker__more" data-more="1">${esc(t('tags.pickMore', 'ещё'))} ${rest.length - FOLD}</button>` : '';
      const hint = `<span class="tag-picker__hint">${esc(t('tags.pickHint', 'нажмите тег, чтобы добавить или убрать'))}</span>`;
      // the word still being typed is not a tag yet
      const lit = filtering ? have.filter((x) => x.toLowerCase() !== typing) : have;
      box.innerHTML = hint + lit.map((x) => chip(x, true)).join('') + shown.map((x) => chip(x, false)).join('') + more;
    }

    box.addEventListener('click', (e) => {
      const btn = e.target.closest('button');
      if (!btn) return;
      e.preventDefault();
      if (btn.dataset.more) { open = true; render(); return; }
      const tag = btn.dataset.tag;
      const have = current();
      // a word half typed is replaced by the tag picked for it (only when the tag goes on
      // from it: a whole new tag typed without a comma stays)
      const typing = /,\s*$/.test(input.value) ? '' : (input.value.split(/[,;]/).pop() || '').trim();
      const half = typing && !all.includes(typing) && tag.toLowerCase().includes(typing.toLowerCase()) ? typing : null;
      const list = have.includes(tag) ? have.filter((x) => x !== tag) : [...have.filter((x) => x !== half), tag];
      setTags(list);
      render();
    });
    input.addEventListener('input', render);

    loadAll().then((counts) => {
      all = [...counts.entries()]
        .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], 'ru'))
        .map(([tag]) => tag);
      render();
    });
  }

  function esc(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }

  window.TagPicker = { attach, loadAll };
})();
