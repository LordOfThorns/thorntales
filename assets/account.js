/* account.js — personal page: display name, avatar, password, logout. */
(function () {
  'use strict';

  const AVATAR_PX = 256;

  const root = document.getElementById('account');
  const sub = document.getElementById('accountSub');
  const preview = document.getElementById('avatarPreview');
  const avatarInput = document.getElementById('avatarInput');
  const avatarRemove = document.getElementById('avatarRemove');
  const nameForm = document.getElementById('nameForm');
  const passForm = document.getElementById('passForm');
  const profileMsg = document.getElementById('profileMsg');
  const passMsg = document.getElementById('passMsg');

  const db = Auth.client;
  const t = I18n.t;

  Auth.onChange(({ user, profile }) => {
    if (!user) {
      location.replace('/login.html?next=' + encodeURIComponent('/account.html'));
      return;
    }
    root.hidden = false;
    sub.textContent = user.email + (profile?.is_admin ? ' · ' + t('account.moderator') : '');
    preview.innerHTML = Auth.avatarHtml(profile, 'avatar--xl');
    avatarRemove.hidden = !profile?.avatar_url;
    if (document.activeElement !== nameForm.name) {
      nameForm.name.value = profile?.display_name || '';
    }
  });

  // ---------- Name ----------

  nameForm.addEventListener('submit', async (e) => {
    e.preventDefault();
    const name = nameForm.name.value.trim().replace(/\s+/g, ' ');
    if (!name) return show(profileMsg, t('account.nameEmpty'), true);
    if (name.length > 40) return show(profileMsg, t('account.nameLong'), true);

    const { error } = await db.from('profiles').update({ display_name: name }).eq('id', Auth.user.id);
    if (error) return show(profileMsg, t('account.saveFailed') + error.message, true);
    show(profileMsg, t('account.saved'));
    Auth.refreshProfile();
  });

  // ---------- Avatar ----------

  avatarInput.addEventListener('change', async () => {
    const file = avatarInput.files[0];
    avatarInput.value = '';
    if (!file) return;
    if (!/^image\/(jpeg|png|webp)$/.test(file.type)) {
      return show(profileMsg, t('account.needImage'), true);
    }

    show(profileMsg, t('account.uploading'));
    try {
      const blob = await squareResize(file, AVATAR_PX);
      const uid = Auth.user.id;
      // New file name each time so browsers/CDN don't show the cached old avatar
      // Safari can't encode WebP and silently returns PNG
      const ext = blob.type === 'image/png' ? 'png' : 'webp';
      const path = `${uid}/avatar-${Date.now()}.${ext}`;

      const up = await db.storage.from('avatars').upload(path, blob, {
        contentType: blob.type || 'image/webp',
        upsert: false,
      });
      if (up.error) throw up.error;

      const { data } = db.storage.from('avatars').getPublicUrl(path);
      const { error } = await db.from('profiles').update({ avatar_url: data.publicUrl }).eq('id', uid);
      if (error) throw error;

      await removeOldAvatars(path);
      show(profileMsg, t('account.avatarUpdated'));
      Auth.refreshProfile();
    } catch (err) {
      show(profileMsg, t('account.uploadFailed') + (err.message || err), true);
    }
  });

  avatarRemove.addEventListener('click', async () => {
    const { error } = await db.from('profiles').update({ avatar_url: null }).eq('id', Auth.user.id);
    if (error) return show(profileMsg, t('account.failed') + error.message, true);
    await removeOldAvatars(null);
    show(profileMsg, t('account.avatarRemoved'));
    Auth.refreshProfile();
  });

  // Delete every file in the user's folder except `keepPath`
  async function removeOldAvatars(keepPath) {
    const uid = Auth.user.id;
    const { data, error } = await db.storage.from('avatars').list(uid, { limit: 100 });
    if (error || !data) return;
    const old = data.map((f) => `${uid}/${f.name}`).filter((p) => p !== keepPath);
    if (old.length) await db.storage.from('avatars').remove(old);
  }

  // Center-crop to a square and scale down, encode as WebP
  async function squareResize(file, size) {
    const bitmap = await createImageBitmap(file);
    const side = Math.min(bitmap.width, bitmap.height);
    const sx = (bitmap.width - side) / 2;
    const sy = (bitmap.height - side) / 2;
    const out = Math.min(size, side);

    const canvas = document.createElement('canvas');
    canvas.width = out;
    canvas.height = out;
    const ctx = canvas.getContext('2d');
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(bitmap, sx, sy, side, side, 0, 0, out, out);
    bitmap.close();

    return new Promise((resolve, reject) => {
      canvas.toBlob((b) => (b ? resolve(b) : reject(new Error('canvas'))), 'image/webp', 0.88);
    });
  }

  // ---------- Password ----------

  passForm.addEventListener('submit', async (e) => {
    e.preventDefault();
    const pass = passForm.pass.value;
    if (pass.length < 8) return show(passMsg, t('account.passShort'), true);
    if (pass !== passForm.pass2.value) return show(passMsg, t('account.passMismatch'), true);

    const { error } = await db.auth.updateUser({ password: pass });
    if (error) return show(passMsg, t('account.failed') + error.message, true);
    passForm.reset();
    show(passMsg, t('account.passChanged'));
  });

  // ---------- Logout ----------

  document.getElementById('logout').addEventListener('click', async () => {
    await Auth.signOut();
    location.href = '/';
  });

  function show(el, text, isError) {
    el.textContent = text;
    el.classList.toggle('is-error', !!isError);
  }
})();
