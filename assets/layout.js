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
        <a class="site-nav__brand" href="${SITE.home.href}"><span class="site-nav__brand-text">Thorntales</span></a>
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
    addThorns(nav.querySelector('.site-nav__brand'));

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

  // ---------- Thorn branches around the logo ----------

  // Two thin stems wind around "Thorntales": the back layer lies under the letters,
  // the front layer repeats the stems clipped to a few short pieces, so the branches
  // look like they pass in front of the word there. Thorns are placed along the
  // stems with getPointAtLength, so the paths must be in the DOM first.
  function addThorns(brand) {
    if (!brand) return;
    const NS = 'http://www.w3.org/2000/svg';
    const W = 132, H = 40;
    const stems = [
      'M3,30 C22,38 34,6 58,13 S92,34 129,9',
      'M6,11 C26,2 44,34 68,27 S104,6 127,31',
    ];
    const frontPieces = [[30, 13], [61, 9], [93, 12]]; // [x, width] in viewBox units

    const back = document.createElementNS(NS, 'svg');
    back.setAttribute('class', 'brand-thorns brand-thorns--back');
    back.setAttribute('viewBox', `0 0 ${W} ${H}`);
    back.setAttribute('aria-hidden', 'true');
    back.setAttribute('focusable', 'false');
    const g = document.createElementNS(NS, 'g');
    g.id = 'brand-thorns-branch';
    back.appendChild(g);
    brand.prepend(back);

    for (const d of stems) {
      const stem = document.createElementNS(NS, 'path');
      stem.setAttribute('d', d);
      stem.setAttribute('class', 'brand-thorns__stem');
      g.appendChild(stem);
      g.appendChild(thornsAlong(stem));
    }

    const front = document.createElementNS(NS, 'svg');
    front.setAttribute('class', 'brand-thorns brand-thorns--front');
    front.setAttribute('viewBox', `0 0 ${W} ${H}`);
    front.setAttribute('aria-hidden', 'true');
    front.setAttribute('focusable', 'false');
    front.innerHTML = `
      <clipPath id="brand-thorns-front-clip">
        ${frontPieces.map(([x, w]) => `<rect x="${x}" y="0" width="${w}" height="${H}"/>`).join('')}
      </clipPath>
      <use href="#brand-thorns-branch" clip-path="url(#brand-thorns-front-clip)"/>`;
    brand.appendChild(front);

    function thornsAlong(stem) {
      const out = document.createElementNS(NS, 'path');
      out.setAttribute('class', 'brand-thorns__thorn');
      const len = stem.getTotalLength();
      if (!len) return out; // not rendered (e.g. display:none), skip thorns
      // Fixed pseudo-random sequence: uneven like a real branch, but the same on every page
      let seed = Math.round(len * 7);
      const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
      const thorn = (x, y, tx, ty, nx, ny, h, b, lean) => {
        const tipX = x + nx * h + tx * h * lean, tipY = y + ny * h + ty * h * lean;
        return `M${(x - tx * b).toFixed(2)},${(y - ty * b).toFixed(2)}` +
               `L${tipX.toFixed(2)},${tipY.toFixed(2)}` +
               `L${(x + tx * b).toFixed(2)},${(y + ty * b).toFixed(2)}Z`;
      };
      let d = '';
      let side = rnd() < 0.5 ? 1 : -1;
      for (let s = 4 + rnd() * 6; s < len - 4; s += 6 + rnd() * 8) {
        const p = stem.getPointAtLength(s);
        const q = stem.getPointAtLength(s + 0.6);
        let tx = q.x - p.x, ty = q.y - p.y;
        const tl = Math.hypot(tx, ty) || 1;
        tx /= tl; ty /= tl;
        const nx = -ty * side, ny = tx * side;
        if (rnd() < 0.16) {
          // A short side twig ending in a sharp point, with a thorn of its own
          const h = 6 + rnd() * 3;
          d += thorn(p.x, p.y, tx, ty, nx, ny, h, 0.7, 0.8);
          const mx = p.x + (nx + tx * 0.8) * h * 0.5, my = p.y + (ny + ty * 0.8) * h * 0.5;
          d += thorn(mx, my, nx, ny, -tx, -ty, 2, 0.5, 0.7);
        } else {
          // thorns lean toward the tip of the stem
          d += thorn(p.x, p.y, tx, ty, nx, ny, 2.2 + rnd() * 1.8, 0.8 + rnd() * 0.3, 0.4 + rnd() * 0.4);
        }
        if (rnd() < 0.8) side = -side;
      }
      out.setAttribute('d', d);
      return out;
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
