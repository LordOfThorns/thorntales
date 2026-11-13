/* tags.js — admin page: who sees poems with a given tag.
   Tables: tag_access (tag, mode), tag_readers (tag, user_id), see supabase/schema.sql.
   Changes are saved right away. */
(async function () {
  'use strict';

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const sub = document.getElementById('tagsSub');
  const page = document.getElementById('tagsPage');
  const listHost = document.getElementById('tagList');
  const addForm = document.getElementById('tagAdd');
  const MODES = ['public', 'signed_in', 'listed'];

  if (!window.Auth) {
    sub.textContent = t('poems.loadFailed');
    return;
  }
  const db = Auth.client;
  await Auth.ready;
  if (!Auth.profile?.is_admin) {
    sub.textContent = t('tags.adminOnly');
    return;
  }

  let poems = [];             // { catalog, poem_id, title, tags }
  const access = new Map();   // tag -> mode
  const readers = new Map();  // tag -> Set of user ids
  let users = [];             // non-admin profiles
  const fileIds = new Map();  // catalog -> Set of poem ids in its txt file
  let catLabels = new Map();
  const extraTags = new Set(); // added on this page, no rule and no poems yet

  try {
    const [p, a, r, u, cats] = await Promise.all([
      must(db.from('poem_edits').select('catalog, poem_id, title, tags')),
      must(db.from('tag_access').select('tag, mode')),
      must(db.from('tag_readers').select('tag, user_id')),
      must(db.from('profiles').select('id, display_name, avatar_url, is_admin').order('display_name')),
      PoetryData.catalogs(),
    ]);
    poems = p.filter(x => x.tags && x.tags.length);
    for (const x of a) access.set(x.tag, x.mode);
    for (const x of r) {
      if (!readers.has(x.tag)) readers.set(x.tag, new Set());
      readers.get(x.tag).add(x.user_id);
    }
    users = u.filter(x => !x.is_admin);
    catLabels = new Map(cats.map(c => [c.id, c.label]));

    // Which tagged poems also sit in a public txt file
    await Promise.all([...new Set(poems.map(x => x.catalog))].map(async (cat) => {
      const res = await fetch(`./data/poetry/${encodeURIComponent(cat)}.txt`);
      const ids = new Set();
      if (res.ok) for (const m of (await res.text()).matchAll(/^id:\s*(.+)$/gm)) ids.add(m[1].trim());
      fileIds.set(cat, ids);
    }));
  } catch (e) {
    sub.textContent = t('poems.loadFailed') + ': ' + (e.message || e);
    return;
  }

  sub.textContent = t('tags.sub');
  page.hidden = false;
  render();

  async function must(query) {
    const { data, error } = await query;
    if (error) throw error;
    return data || [];
  }

  function allTags() {
    const set = new Set([...access.keys(), ...extraTags]);
    for (const p of poems) for (const tag of p.tags) set.add(tag);
    return [...set].sort((a, b) => a.localeCompare(b, 'ru'));
  }

  function render() {
    const tags = allTags();
    if (!tags.length) {
      listHost.innerHTML = `<p class="sub">${escapeHtml(t('tags.empty'))}</p>`;
      return;
    }
    listHost.innerHTML = tags.map(renderTag).join('');
  }

  function renderTag(tag) {
    const mode = access.get(tag) || 'public';
    const tagged = poems.filter(p => p.tags.includes(tag));
    const who = readers.get(tag) || new Set();
    const key = encodeURIComponent(tag);

    const modes = MODES.map(m => `
      <label class="tag-mode">
        <input type="radio" name="mode-${escapeHtml(key)}" value="${m}" data-tag="${escapeHtml(tag)}" ${m === mode ? 'checked' : ''}>
        <span>${escapeHtml(t('tags.mode.' + m))}</span>
      </label>`).join('');

    const userList = mode !== 'listed' ? '' : `
      <div class="tag-readers">
        ${users.length ? users.map(u => `
          <label class="tag-reader">
            <input type="checkbox" data-tag="${escapeHtml(tag)}" data-user="${escapeHtml(u.id)}" ${who.has(u.id) ? 'checked' : ''}>
            ${Auth.avatarHtml(u, 'tag-reader__avatar')}
            <span>${escapeHtml(u.display_name)}</span>
          </label>`).join('') : `<p class="field__hint">${escapeHtml(t('tags.noUsers'))}</p>`}
        <p class="field__hint">${escapeHtml(t('tags.adminSeesAll'))}</p>
      </div>`;

    const poemList = tagged.length ? `
      <details class="tag-poems">
        <summary>${escapeHtml(t('tags.poemCount', { n: tagged.length }))}</summary>
        <ul>
          ${tagged.map(p => {
            const inFile = fileIds.get(p.catalog)?.has(p.poem_id);
            return `<li>
              <a href="/poems.html?cat=${encodeURIComponent(p.catalog)}#${encodeURIComponent(p.poem_id)}">${escapeHtml(p.title)}</a>
              <span class="tag-poems__cat">${escapeHtml(catLabels.get(p.catalog) || p.catalog)}</span>
              ${inFile && mode !== 'public' ? `<span class="tag-poems__file" title="${escapeHtml(t('tags.inFileHint'))}">${escapeHtml(t('tags.inFile'))}</span>` : ''}
            </li>`;
          }).join('')}
        </ul>
      </details>` : `<p class="field__hint">${escapeHtml(t('tags.noPoems'))}</p>`;

    return `
      <section class="panel tag-card${mode !== 'public' ? ' tag-card--restricted' : ''}" data-tag="${escapeHtml(tag)}">
        <h2 class="tag-card__name">#${escapeHtml(tag)}</h2>
        <div class="tag-modes">${modes}</div>
        ${userList}
        ${poemList}
        <p class="form__msg" role="status"></p>
      </section>`;
  }

  // ---------- Saving ----------

  listHost.addEventListener('change', async (e) => {
    const el = e.target;
    const tag = el.dataset.tag;
    if (!tag) return;
    const card = el.closest('.tag-card');
    const msg = card.querySelector('.form__msg');
    card.querySelectorAll('input').forEach(i => { i.disabled = true; });
    msg.classList.remove('is-error');
    msg.textContent = t('tags.saving');

    try {
      if (el.type === 'radio') {
        await must(db.from('tag_access').upsert({ tag, mode: el.value }).select());
        access.set(tag, el.value);
      } else {
        // A reader can be added only to a tag that has a rule row
        if (!access.has(tag)) {
          await must(db.from('tag_access').upsert({ tag, mode: 'listed' }).select());
          access.set(tag, 'listed');
        }
        const userId = el.dataset.user;
        if (el.checked) {
          await must(db.from('tag_readers').upsert({ tag, user_id: userId }).select());
          if (!readers.has(tag)) readers.set(tag, new Set());
          readers.get(tag).add(userId);
        } else {
          await must(db.from('tag_readers').delete().eq('tag', tag).eq('user_id', userId).select());
          readers.get(tag)?.delete(userId);
        }
      }
    } catch (err) {
      render();
      const again = listHost.querySelector(`.tag-card[data-tag="${CSS.escape(tag)}"] .form__msg`);
      if (again) {
        again.textContent = t('poems.saveFailed') + (err.message || err);
        again.classList.add('is-error');
      }
      return;
    }
    render();
    const done = listHost.querySelector(`.tag-card[data-tag="${CSS.escape(tag)}"] .form__msg`);
    if (done) done.textContent = t('tags.saved');
  });

  addForm.addEventListener('submit', (e) => {
    e.preventDefault();
    const tag = addForm.elements.tag.value.trim().replace(/^#+/, '').replace(/\s+/g, ' ').slice(0, 40);
    if (!tag) return;
    extraTags.add(tag);
    addForm.reset();
    render();
    listHost.querySelector(`.tag-card[data-tag="${CSS.escape(tag)}"]`)?.scrollIntoView({ block: 'center' });
  });

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, ch => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }
})();
