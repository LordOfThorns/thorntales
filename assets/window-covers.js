/* window-covers.js — pictures of the stained-glass windows on the home page.

   A window shows its picture from index.html (--img) unless table public.site_images has
   a row 'window:<page>' (window:poetry, window:texts, window:polar-night, window:projects):
   then that picture (bucket site-images) is used. The last known picture is remembered in
   this browser, so a replaced cover shows at once on the next visit.

   The admin gets two small buttons at the top of each window, shown on hover (always on
   touch screens): "change cover" (a file is scaled down and saved as WebP) and "back to
   the original". They stop the click, so the glass does not break.
*/
(function () {
  'use strict';

  const t = (key) => (window.I18n ? I18n.t(key) : key);
  const db = window.Auth ? Auth.client : null;
  const CACHE_KEY = 'thorntales.windowCovers';
  const MAX_W = 800, MAX_H = 1334; // the window is 3:5, a 2x screen shows up to ~500x830

  const windows = [...document.querySelectorAll('.window')].map((el) => ({
    el,
    key: 'window:' + ((el.getAttribute('href') || '').replace(/^\/|\.html$/g, '') || 'home'),
    original: el.style.getPropertyValue('--img'),
    custom: '',
  }));
  if (!windows.length) return;

  // 1. Pictures remembered in this browser: no flash of the old picture
  let cache = {};
  try { cache = JSON.parse(localStorage.getItem(CACHE_KEY) || '{}') || {}; } catch (e) { cache = {}; }
  for (const w of windows) if (cache[w.key]) apply(w, cache[w.key]);

  // 2. The real state from the site
  load();

  async function load() {
    if (!db) return;
    const { data, error } = await db.from('site_images').select('key, url').like('key', 'window:%');
    if (error) { console.warn('site images failed', error); return; }
    const byKey = new Map((data || []).map((r) => [r.key, r.url]));
    for (const w of windows) apply(w, byKey.get(w.key) || '');
    remember();
  }

  function apply(w, url) {
    w.custom = url || '';
    if (url) w.el.style.setProperty('--img', `url("${url.replace(/"/g, '%22')}")`);
    else w.el.style.setProperty('--img', w.original);
    const reset = w.el.querySelector('[data-cover="reset"]');
    if (reset) reset.hidden = !w.custom;
  }

  function remember() {
    const out = {};
    for (const w of windows) if (w.custom) out[w.key] = w.custom;
    try { localStorage.setItem(CACHE_KEY, JSON.stringify(out)); } catch (e) { /* storage blocked */ }
  }

  // ---------- Admin buttons ----------

  if (window.Auth) {
    Auth.onChange(({ profile }) => {
      const admin = !!profile?.is_admin;
      for (const w of windows) {
        const box = w.el.querySelector('.cover-tools');
        if (admin && !box) addTools(w);
        if (!admin && box) box.remove();
      }
    });
  }

  function addTools(w) {
    const box = document.createElement('span');
    box.className = 'cover-tools';
    box.innerHTML = `
      <span class="cover-tools__btn" role="button" tabindex="0" data-cover="change"
            title="${esc(t('home.coverChange'))}" aria-label="${esc(t('home.coverChange'))}">
        <svg viewBox="0 0 16 16" aria-hidden="true"><path d="M2 12.5V14h1.5l8-8L10 4.5l-8 8zM12.6 3.4l-1-1a.9.9 0 0 0-1.3 0l-.8.8L11 4.7l.8-.8a.9.9 0 0 0 0-1.3z"/></svg>
      </span>
      <span class="cover-tools__btn" role="button" tabindex="0" data-cover="reset" ${w.custom ? '' : 'hidden'}
            title="${esc(t('home.coverReset'))}" aria-label="${esc(t('home.coverReset'))}">
        <svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 3a5 5 0 1 1-4.6 3H1.8A6.5 6.5 0 1 0 3.6 3.5L2 2v4h4L4.6 4.6A5 5 0 0 1 8 3z"/></svg>
      </span>
      <span class="cover-tools__msg" role="status"></span>
      <input type="file" accept="image/jpeg,image/png,image/webp" hidden>`;
    w.el.appendChild(box);

    const file = box.querySelector('input[type=file]');
    const msg = box.querySelector('.cover-tools__msg');

    // Clicks here must not break the glass nor follow the link
    const stop = (e) => { e.preventDefault(); e.stopPropagation(); };
    box.addEventListener('click', (e) => {
      if (e.target === file) { e.stopPropagation(); return; }
      stop(e);
      const act = e.target.closest('[data-cover]')?.dataset.cover;
      if (act === 'change') file.click();
      if (act === 'reset') reset(w, msg);
    });
    box.addEventListener('keydown', (e) => {
      if ((e.key === 'Enter' || e.key === ' ') && e.target.dataset.cover) { stop(e); e.target.click(); }
    });
    file.addEventListener('change', () => {
      const f = file.files[0];
      file.value = '';
      if (f) upload(w, f, msg);
    });
  }

  async function upload(w, f, msg) {
    if (!/^image\/(jpeg|png|webp)$/.test(f.type)) { say(msg, t('home.coverNeedImage'), true); return; }
    say(msg, t('home.coverSaving'));
    try {
      const blob = await fitResize(f, MAX_W, MAX_H);
      const folder = 'windows/' + w.key.slice('window:'.length);
      const path = `${folder}/${Date.now()}.webp`;
      const bucket = db.storage.from('site-images');
      const up = await bucket.upload(path, blob, { contentType: 'image/webp', upsert: false });
      if (up.error) throw up.error;
      const url = bucket.getPublicUrl(path).data.publicUrl;
      const { error } = await db.from('site_images').upsert({ key: w.key, url });
      if (error) throw error;
      await removeOld(folder, path);
      apply(w, url);
      remember();
      say(msg, '');
    } catch (e) {
      console.error(e);
      say(msg, t('home.coverFailed') + ' ' + (e.message || e), true);
    }
  }

  async function reset(w, msg) {
    if (!confirm(t('home.coverConfirmReset'))) return;
    say(msg, t('home.coverSaving'));
    try {
      const { error } = await db.from('site_images').delete().eq('key', w.key);
      if (error) throw error;
      await removeOld('windows/' + w.key.slice('window:'.length), null);
      apply(w, '');
      remember();
      say(msg, '');
    } catch (e) {
      console.error(e);
      say(msg, t('home.coverFailed') + ' ' + (e.message || e), true);
    }
  }

  async function removeOld(folder, keepPath) {
    const bucket = db.storage.from('site-images');
    const { data } = await bucket.list(folder, { limit: 100 });
    const old = (data || []).map((x) => `${folder}/${x.name}`).filter((p) => p !== keepPath);
    if (old.length) await bucket.remove(old);
  }

  function say(el, text, isError) {
    el.textContent = text;
    el.classList.toggle('is-error', !!isError);
  }

  // Scale down to fit maxW x maxH, keeping proportions; WebP
  async function fitResize(file, maxW, maxH) {
    const bitmap = await createImageBitmap(file);
    const k = Math.min(1, maxW / bitmap.width, maxH / bitmap.height);
    const canvas = document.createElement('canvas');
    canvas.width = Math.max(1, Math.round(bitmap.width * k));
    canvas.height = Math.max(1, Math.round(bitmap.height * k));
    const ctx = canvas.getContext('2d');
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(bitmap, 0, 0, canvas.width, canvas.height);
    bitmap.close();
    return new Promise((resolve, reject) => {
      canvas.toBlob((b) => (b ? resolve(b) : reject(new Error('canvas'))), 'image/webp', 0.9);
    });
  }

  function esc(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;' }[ch]));
  }
})();
