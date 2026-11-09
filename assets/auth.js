/* auth.js
   Supabase client + current user/profile, shared by all pages.
   Needs (in this order): layout.js, supabase-js UMD bundle, auth.js.

   API:
     Auth.client                 // supabase client
     Auth.ready                  // Promise resolved after the session/profile is known
     Auth.user, Auth.profile     // current user / row from public.profiles (or null)
     Auth.onChange(fn)           // fn({ user, profile }) now and on every change
     Auth.refreshProfile()
     Auth.signOut()
*/
(function () {
  'use strict';

  // Public by design: the publishable key only allows what the database policies allow.
  const SUPABASE_URL = 'https://hnggliqiezorttzowuaq.supabase.co';
  const SUPABASE_KEY = 'sb_publishable_nP0LO5i-tID9J7GMyI9czg_Db2vZDTE';

  if (!window.supabase || !window.supabase.createClient) {
    console.error('supabase-js is not loaded');
    return;
  }

  const client = window.supabase.createClient(SUPABASE_URL, SUPABASE_KEY, {
    auth: { persistSession: true, autoRefreshToken: true },
  });

  const listeners = new Set();
  const state = { user: null, profile: null };
  let isReady = false;
  let resolveReady;
  const ready = new Promise((r) => { resolveReady = () => { isReady = true; r(); }; });

  async function loadProfile(user) {
    if (!user) return null;
    const { data, error } = await client
      .from('profiles')
      .select('id, display_name, avatar_url, is_admin')
      .eq('id', user.id)
      .maybeSingle();
    if (error) console.warn('profile load failed', error);
    return data || null;
  }

  async function setUser(user) {
    state.user = user || null;
    state.profile = await loadProfile(state.user);
    renderNavUser();
    listeners.forEach((fn) => fn({ ...state }));
  }

  // ---------- Menu slot ----------

  const tr = (key, fallback) => (window.I18n ? I18n.t(key) : fallback);

  function renderNavUser() {
    const slot = window.Layout && Layout.userSlot();
    if (!slot) return;

    if (!state.user) {
      const back = encodeURIComponent(location.pathname + location.search + location.hash);
      slot.innerHTML = `<a class="nav-login" href="/login.html?next=${back}">${escapeHtml(tr('nav.signIn', 'Войти'))}</a>`;
      return;
    }

    const name = state.profile?.display_name || state.user.email || tr('nav.account', 'Кабинет');
    slot.innerHTML = `
      <a class="nav-user" href="/account.html" title="${escapeHtml(tr('nav.accountTitle', 'Личный кабинет'))}">
        ${avatarHtml(state.profile, 'nav-user__avatar')}
        <span class="nav-user__name">${escapeHtml(name)}</span>
      </a>`;
  }

  // Avatar image or a letter circle (shared with comments/account)
  function avatarHtml(profile, cls) {
    const name = profile?.display_name || '?';
    if (profile?.avatar_url) {
      return `<img class="avatar ${cls || ''}" src="${escapeHtml(profile.avatar_url)}" alt="" loading="lazy">`;
    }
    const letter = name.trim().charAt(0).toUpperCase() || '?';
    // Stable hue per name so letter avatars differ
    let h = 0;
    for (const ch of name) h = (h * 31 + ch.codePointAt(0)) % 360;
    return `<span class="avatar avatar--letter ${cls || ''}" style="--avatar-hue:${h}" aria-hidden="true">${escapeHtml(letter)}</span>`;
  }

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }

  // ---------- Start ----------

  client.auth.onAuthStateChange((event, session) => {
    // Supabase advises not to await other supabase calls inside this callback
    setTimeout(async () => {
      const prevId = state.user?.id;
      const nextUser = session?.user || null;
      if (event === 'INITIAL_SESSION' || prevId !== nextUser?.id || event === 'USER_UPDATED') {
        await setUser(nextUser);
      }
      resolveReady();
    }, 0);
  });

  window.Auth = {
    client,
    ready,
    get user() { return state.user; },
    get profile() { return state.profile; },
    onChange(fn) {
      listeners.add(fn);
      // Listeners added before the first session check get called by setUser()
      if (isReady) fn({ ...state });
      return () => listeners.delete(fn);
    },
    async refreshProfile() {
      state.profile = await loadProfile(state.user);
      renderNavUser();
      listeners.forEach((fn) => fn({ ...state }));
    },
    async signOut() {
      await client.auth.signOut();
    },
    avatarHtml,
    escapeHtml,
  };
})();
