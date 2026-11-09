/* layout.js
   Shared page chrome: top menu, breadcrumbs, back-to-top button.

   Usage on a page:
     <body data-page="poetry">            // id from SITE below
     ...
     <script src="/assets/layout.js"></script>   // before page scripts

   Pages with dynamic content can append crumbs:
     Layout.setCrumbs([{ label: 'Сальвии' }]);              // appended after the page itself
     Layout.setCrumbs([{ label: 'X', href: '/x.html' }, { label: 'Y' }]);
*/
(function () {
  'use strict';

  // Site tree. `parent` builds breadcrumbs (any depth), `menu: true` puts a page into the top menu.
  // To nest deeper later, just add a page with parent: '<id>'.
  const SITE = {
    home:       { label: 'Башня',          href: '/' },
    polarNight: { label: 'Полярная Ночь',  href: '/polar-night.html', parent: 'home', menu: true },
    texts:      { label: 'Тексты',         href: '/texts.html',       parent: 'home', menu: true },
    poetry:     { label: 'Стихи',          href: '/poetry.html',      parent: 'home', menu: true },
    projects:   { label: 'Проекты',        href: '/projects.html',    parent: 'home', menu: true },
    mapGen:     { label: 'Генератор карт', href: '/terrain.html',     parent: 'projects' },
    login:      { label: 'Вход',           href: '/login.html',       parent: 'home' },
    account:    { label: 'Кабинет',        href: '/account.html',     parent: 'home' },
  };

  const pageId = document.body.dataset.page || '';
  let extraCrumbs = [];

  // Labels come from i18n.js (key nav.<page id>); the Russian label above is the fallback
  const t = (key, fallback) =>
    window.I18n && I18n.has(key) ? I18n.t(key) : (fallback ?? (window.I18n ? I18n.t(key) : key));
  for (const [id, p] of Object.entries(SITE)) p.label = t('nav.' + id, p.label);

  function buildLangSwitch() {
    if (!window.I18n) return '';
    const buttons = I18n.langs.map((l) =>
      `<button type="button" data-lang="${l}" aria-pressed="${l === I18n.lang}">${l.toUpperCase()}</button>`
    ).join('');
    return `<div class="lang-switch" role="group" aria-label="${escapeHtml(t('nav.language', 'Язык'))}">${buttons}</div>`;
  }

  // ---------- Top menu ----------

  function chainOf(id) {
    const chain = [];
    let cur = id;
    while (cur && SITE[cur]) {
      chain.unshift(cur);
      cur = SITE[cur].parent;
    }
    return chain;
  }

  function buildNav() {
    const activeChain = chainOf(pageId);
    const items = Object.entries(SITE)
      .filter(([, p]) => p.menu)
      .map(([id, p]) => {
        const current = activeChain.includes(id) ? ' aria-current="page"' : '';
        return `<li><a href="${p.href}"${current}>${escapeHtml(p.label)}</a></li>`;
      })
      .join('');

    const nav = document.createElement('header');
    nav.className = 'site-nav' + (pageId === 'home' ? ' site-nav--home' : '');
    nav.innerHTML = `
      <div class="site-nav__inner">
        <a class="site-nav__brand" href="${SITE.home.href}">Thorntales</a>
        <nav class="site-nav__menu" id="site-nav-menu" aria-label="${escapeHtml(t('nav.mainMenu', 'Основное меню'))}">
          <ul>${items}</ul>
          <div class="site-nav__user" id="site-nav-user"></div>
        </nav>
        ${buildLangSwitch()}
        <button class="site-nav__toggle" type="button" aria-expanded="false" aria-controls="site-nav-menu" aria-label="${escapeHtml(t('nav.menu', 'Меню'))}">
          <span></span><span></span><span></span>
        </button>
      </div>`;

    document.body.prepend(nav);

    nav.querySelectorAll('.lang-switch [data-lang]').forEach((b) => {
      b.addEventListener('click', () => I18n.setLang(b.dataset.lang));
    });

    const toggle = nav.querySelector('.site-nav__toggle');
    toggle.addEventListener('click', () => {
      const open = nav.classList.toggle('is-open');
      toggle.setAttribute('aria-expanded', String(open));
    });
    // Close the mobile menu on Escape or outside click
    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') closeMenu();
    });
    document.addEventListener('click', (e) => {
      if (!nav.contains(e.target)) closeMenu();
    });
    function closeMenu() {
      nav.classList.remove('is-open');
      toggle.setAttribute('aria-expanded', 'false');
    }
  }

  // ---------- Breadcrumbs ----------

  function renderCrumbs() {
    const chain = chainOf(pageId);
    if (chain.length < 2 && extraCrumbs.length === 0) return; // nothing to show on home

    const crumbs = chain.map((id) => ({ label: SITE[id].label, href: SITE[id].href }));
    crumbs.push(...extraCrumbs);
    // Last crumb is the current location: no link
    const last = crumbs.length - 1;

    let host = document.querySelector('.crumbs');
    if (!host) {
      host = document.createElement('nav');
      host.className = 'crumbs';
      host.setAttribute('aria-label', t('nav.crumbs', 'Вы здесь'));
      const wrap = document.querySelector('main.wrap') || document.body;
      wrap.prepend(host);
    }
    host.innerHTML = '<ol>' + crumbs.map((c, i) =>
      i === last || !c.href
        ? `<li><span aria-current="page">${escapeHtml(c.label)}</span></li>`
        : `<li><a href="${c.href}">${escapeHtml(c.label)}</a></li>`
    ).join('') + '</ol>';
  }

  // ---------- Back to top ----------

  function buildBackToTop() {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'to-top';
    btn.setAttribute('aria-label', t('nav.toTop', 'Наверх'));
    btn.textContent = '↑';
    btn.addEventListener('click', () => {
      const reduce = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
      window.scrollTo({ top: 0, behavior: reduce ? 'auto' : 'smooth' });
    });
    document.body.appendChild(btn);

    const update = () => btn.classList.toggle('is-visible', window.scrollY > 700);
    window.addEventListener('scroll', update, { passive: true });
    update();
  }

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }

  buildNav();
  renderCrumbs();
  buildBackToTop();

  window.Layout = {
    site: SITE,
    setCrumbs(list) {
      extraCrumbs = Array.isArray(list) ? list : [];
      renderCrumbs();
    },
    // Slot for the future login / user name (Supabase)
    userSlot: () => document.getElementById('site-nav-user'),
  };
})();
