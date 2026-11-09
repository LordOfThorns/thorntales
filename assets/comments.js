/* comments.js
   Collapsible comment threads under content blocks.

   Usage (after the blocks are in the DOM):
     Comments.mount([{ el: <container element>, thread: 'poem/salvia/<poem id>' }, ...], 'poem/salvia/');
   The second argument is a common thread prefix used to load all counters in one request.
*/
(function () {
  'use strict';

  const MAX_LEN = 4000;
  const t = (key) => (window.I18n ? I18n.t(key) : key);
  const items = new Map(); // thread -> { el, thread, count, open, list, loaded }

  function mount(list, prefix) {
    if (!window.Auth) return; // supabase-js failed to load: no comments, page still works

    for (const { el, thread } of list) {
      const box = document.createElement('section');
      box.className = 'comments';
      box.innerHTML = `
        <button class="comments__toggle" type="button" aria-expanded="false">
          <span class="comments__toggle-text">${t('comments.title')}</span>
          <span class="comments__count" hidden></span>
        </button>
        <div class="comments__panel" hidden></div>`;
      el.appendChild(box);

      const item = { el: box, thread, count: 0, open: false, list: [], loaded: false };
      items.set(thread, item);
      box.querySelector('.comments__toggle').addEventListener('click', () => toggle(item));
    }

    loadCounts(prefix);
    Auth.onChange(() => {
      for (const item of items.values()) if (item.open) renderPanel(item);
    });
  }

  // ---------- Data ----------

  async function loadCounts(prefix) {
    const { data, error } = await Auth.client
      .from('comments')
      .select('thread')
      .like('thread', prefix + '%');
    if (error) {
      console.warn('comment counts failed', error);
      return;
    }
    const counts = {};
    for (const row of data) counts[row.thread] = (counts[row.thread] || 0) + 1;
    for (const item of items.values()) {
      item.count = counts[item.thread] || 0;
      renderToggle(item);
    }
  }

  async function loadThread(item) {
    const { data, error } = await Auth.client
      .from('comments')
      .select('id, body, created_at, updated_at, user_id, author:profiles(display_name, avatar_url, is_admin)')
      .eq('thread', item.thread)
      .order('created_at', { ascending: true });
    if (error) throw error;
    item.list = data;
    item.count = data.length;
    item.loaded = true;
  }

  // ---------- UI ----------

  async function toggle(item) {
    item.open = !item.open;
    const panel = item.el.querySelector('.comments__panel');
    item.el.querySelector('.comments__toggle').setAttribute('aria-expanded', String(item.open));
    item.el.classList.toggle('is-open', item.open);
    panel.hidden = !item.open;
    if (!item.open) return;

    if (!item.loaded) {
      panel.innerHTML = `<p class="comments__note">${t('comments.loading')}</p>`;
      try {
        await loadThread(item);
      } catch (e) {
        panel.innerHTML = `<p class="comments__note is-error">${t('comments.loadFailed')}</p>`;
        return;
      }
    }
    renderToggle(item);
    renderPanel(item);
  }

  function renderToggle(item) {
    const c = item.el.querySelector('.comments__count');
    c.hidden = item.count === 0;
    c.textContent = item.count;
  }

  function renderPanel(item) {
    const panel = item.el.querySelector('.comments__panel');
    const me = Auth.user;
    const admin = !!Auth.profile?.is_admin;

    const listHtml = item.list.length
      ? `<ol class="comments__list">${item.list.map((c) => commentHtml(c, me, admin)).join('')}</ol>`
      : `<p class="comments__note">${t('comments.empty')}</p>`;

    const formHtml = me
      ? `<form class="comment-form">
           <textarea class="input comment-form__text" name="body" rows="3" maxlength="${MAX_LEN}"
             placeholder="${t('comments.placeholder')}" required></textarea>
           <div class="comment-form__row">
             <span class="comment-form__hint">${t('comments.hint')}</span>
             <button class="btn btn--primary" type="submit">${t('comments.send')}</button>
           </div>
           <p class="form__msg" role="status"></p>
         </form>`
      : `<p class="comments__note"><a href="/login.html?next=${encodeURIComponent(location.pathname + location.search + '#' + (item.el.closest('[id]')?.id || ''))}">${t('comments.signInLink')}</a>${t('comments.signInRest')}</p>`;

    // Keep a draft that was being typed when the panel re-renders
    const draft = panel.querySelector('.comment-form__text')?.value || '';
    panel.innerHTML = listHtml + formHtml;
    const ta = panel.querySelector('.comment-form__text');
    if (ta) ta.value = draft;

    wirePanel(item, panel);
  }

  function commentHtml(c, me, admin) {
    const author = c.author || {};
    const mine = me && c.user_id === me.id;
    const actions = [
      mine ? `<button class="link-btn" data-act="edit">${t('comments.edit')}</button>` : '',
      mine || admin ? `<button class="link-btn" data-act="delete">${t('comments.delete')}</button>` : '',
    ].join('');

    return `
      <li class="comment" data-id="${c.id}">
        ${Auth.avatarHtml(author, 'comment__avatar')}
        <div class="comment__main">
          <div class="comment__head">
            <span class="comment__name">${Auth.escapeHtml(author.display_name || t('comments.someone'))}</span>
            ${author.is_admin ? `<span class="comment__badge">${t('comments.author')}</span>` : ''}
            <time class="comment__date" datetime="${c.created_at}">${formatDate(c.created_at)}${c.updated_at ? ' · ' + t('comments.edited') : ''}</time>
          </div>
          <div class="comment__body">${Auth.escapeHtml(c.body)}</div>
          ${actions ? `<div class="comment__actions">${actions}</div>` : ''}
        </div>
      </li>`;
  }

  function wirePanel(item, panel) {
    const form = panel.querySelector('.comment-form');
    if (form) {
      const ta = form.querySelector('textarea');
      ta.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) form.requestSubmit();
      });
      form.addEventListener('submit', async (e) => {
        e.preventDefault();
        const body = ta.value.trim();
        const msg = form.querySelector('.form__msg');
        if (!body) return;
        const btn = form.querySelector('button');
        btn.disabled = true;
        const { error } = await Auth.client.from('comments').insert({ thread: item.thread, body });
        btn.disabled = false;
        if (error) {
          msg.textContent = t('comments.sendFailed') + error.message;
          msg.classList.add('is-error');
          return;
        }
        ta.value = '';
        await reload(item);
      });
    }

    panel.querySelectorAll('.comment [data-act]').forEach((btn) => {
      const li = btn.closest('.comment');
      const id = Number(li.dataset.id);
      if (btn.dataset.act === 'delete') {
        btn.addEventListener('click', async () => {
          if (!confirm(t('comments.confirmDelete'))) return;
          const { error } = await Auth.client.from('comments').delete().eq('id', id);
          if (error) return alert(t('comments.deleteFailed') + error.message);
          await reload(item);
        });
      } else if (btn.dataset.act === 'edit') {
        btn.addEventListener('click', () => startEdit(item, li, id));
      }
    });
  }

  function startEdit(item, li, id) {
    const c = item.list.find((x) => x.id === id);
    const bodyEl = li.querySelector('.comment__body');
    const actions = li.querySelector('.comment__actions');
    actions.hidden = true;
    bodyEl.innerHTML = `
      <form class="comment-form comment-form--edit">
        <textarea class="input comment-form__text" rows="3" maxlength="${MAX_LEN}" required></textarea>
        <div class="comment-form__row">
          <button class="link-btn" type="button" data-cancel>${t('comments.cancel')}</button>
          <button class="btn btn--primary" type="submit">${t('comments.save')}</button>
        </div>
      </form>`.trim(); // .comment__body keeps whitespace (pre-wrap), so no stray spaces around the form
    const form = bodyEl.querySelector('form');
    const ta = form.querySelector('textarea');
    ta.value = c.body;
    ta.focus();
    ta.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) form.requestSubmit();
      if (e.key === 'Escape') renderPanel(item);
    });
    form.querySelector('[data-cancel]').addEventListener('click', () => renderPanel(item));
    form.addEventListener('submit', async (e) => {
      e.preventDefault();
      const body = ta.value.trim();
      if (!body) return;
      const { error } = await Auth.client.from('comments').update({ body }).eq('id', id);
      if (error) return alert(t('comments.saveFailed') + error.message);
      await reload(item);
    });
  }

  async function reload(item) {
    try {
      await loadThread(item);
    } catch (e) {
      console.warn(e);
    }
    renderToggle(item);
    renderPanel(item);
  }

  function formatDate(iso) {
    const d = new Date(iso);
    const sameYear = d.getFullYear() === new Date().getFullYear();
    return d.toLocaleString(window.I18n ? I18n.locale : 'ru-RU', {
      day: 'numeric',
      month: 'long',
      ...(sameYear ? {} : { year: 'numeric' }),
      hour: '2-digit',
      minute: '2-digit',
    });
  }

  window.Comments = { mount };
})();
