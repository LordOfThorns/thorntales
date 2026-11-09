/* poems.js — one poetry catalog: poems.html?cat=<catalog id>

   Original:     data/poetry/<cat>.txt
   Translation:  data/poetry/<cat>.<lang>.txt  (optional, e.g. salvia.en.txt)
   Poems are matched by their `id:` line, see data/TRANSLATIONS.md.
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
  const catalogsUrl = dataDir + 'catalogs.json';

  try {
    const [raw, catsJson] = await Promise.all([
      fetchText(`${dataDir}${encodeURIComponent(cat)}.txt`, true),
      fetch(catalogsUrl, { cache: 'no-store' }).then(r => {
        if (!r.ok) throw new Error(`${catalogsUrl}: ${r.status}`);
        return r.json();
      }),
    ]);

    const catalogs = catsJson.catalogs || [];
    const catIndex = catalogs.findIndex(x => x.id === cat);
    const catInfo = catalogs[catIndex];

    // Title and subtitle of the catalog (subtitle may contain <br> only)
    const pageLabel = field(catInfo, 'label') || humanizeCat(cat);
    pageTitle.textContent = pageLabel;
    document.title = `${pageLabel} · Thorntales`;
    if (window.Layout) Layout.setCrumbs([{ label: pageLabel }]);
    const desc = field(catInfo, 'desc');
    pageSub.innerHTML = desc ? allowOnlyBr(desc) : '';

    const poems = parsePoemsTxt(raw);

    // Translation for the reader's language, if this catalog is written in another one
    const sourceLang = catInfo?.lang || 'ru';
    let translatedCount = 0;
    if (lang !== sourceLang) {
      const trRaw = await fetchText(`${dataDir}${encodeURIComponent(cat)}.${lang}.txt`, false);
      if (trRaw) {
        const byId = new Map(parsePoemsTxt(trRaw, { keepOrder: true }).map(p => [p.id, p]));
        for (const p of poems) {
          const tr = byId.get(p.id);
          if (tr && tr.text.trim()) {
            p.tr = tr;
            translatedCount++;
          }
        }
      }
      if (translatedCount) {
        pageSub.insertAdjacentHTML('afterend',
          `<p class="poetry-translated">${escapeHtml(t('poems.translatedCount', { n: translatedCount, total: poems.length }))}</p>`);
      }
    }
    const showUntranslatedMark = lang !== sourceLang;

    // TOC: sidebar on wide screens, collapsible block on narrow ones
    toc.innerHTML = `
      <details class="poetry-toc__details">
        <summary class="poetry-toc__summary">${escapeHtml(t('poems.toc'))} <span class="poetry-toc__count">${poems.length}</span></summary>
        <div class="poetry-toc__inner">
          ${poems.map(p => `
            <a class="poetry-toc__link" href="#${escapeAttr(p.id)}">
              <span class="poetry-toc__date">${escapeHtml(p.date || '')}</span>
              <span class="poetry-toc__title">${escapeHtml((p.tr || p).title || '')}</span>
            </a>
          `).join('')}
        </div>
      </details>
    `;
    setupTocMode(toc.querySelector('details'));

    poemsHost.innerHTML = poems.map(p => renderPoem(p, showUntranslatedMark)).join('');
    wireLangToggles(poems);

    // Comment threads: key = poem/<catalog>/<poem id>, shared by all languages,
    // so poem ids in the txt files must stay stable
    if (window.Comments) {
      const prefix = `poem/${cat}/`;
      Comments.mount(
        poems.map(p => ({ el: document.getElementById(p.id), thread: prefix + p.id })).filter(x => x.el),
        prefix
      );
    }

    // Neighbour catalogs (order as in catalogs.json)
    renderPager(catalogs, catIndex);

    // Content is loaded after the browser tried to jump to #hash, so jump again
    if (location.hash) {
      const target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
      if (target) target.scrollIntoView();
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

  // required=false: a missing file (404) returns null instead of throwing
  async function fetchText(url, required) {
    const res = await fetch(url, { cache: 'no-store' });
    if (!res.ok) {
      if (!required && res.status === 404) return null;
      throw new Error(`${url}: ${res.status}`);
    }
    return res.text();
  }

  function setupTocMode(details){
    if (!details) return;
    const wide = window.matchMedia('(min-width: 1000px)');
    const apply = () => {
      details.open = wide.matches;
      details.classList.toggle('is-sidebar', wide.matches);
    };
    apply();
    wide.addEventListener('change', apply);

    // On narrow screens close the TOC after jumping to a poem
    details.addEventListener('click', (e) => {
      if (!wide.matches && e.target.closest('.poetry-toc__link')) details.open = false;
    });
  }

  function renderPager(catalogs, i){
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

  function parsePoemsTxt(raw, opts = {}) {
    const text = raw.replace(/\r\n/g, '\n');

    // Split into blocks by the marker; anything before the first marker is ignored
    let blocks = text.split('\n=== poem ===\n');
    if (text.startsWith('=== poem ===')) {
      blocks = text.split('=== poem ===\n').slice(1);
    }

    const poems = [];

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

      poems.push({ id, date, title, text: bodyPart });
    }

    // Sort by date (ISO dates)
    if (!opts.keepOrder) poems.sort((a, b) => (a.date || '').localeCompare(b.date || ''));

    return poems;
  }

  function renderPoem(p, showUntranslatedMark){
    const shown = p.tr || p;
    let langControl = '';
    if (p.tr) {
      langControl = `<button class="poem__lang-toggle" type="button" data-poem="${escapeAttr(p.id)}">${escapeHtml(t('poems.showOriginal'))}</button>`;
    } else if (showUntranslatedMark) {
      langControl = `<span class="poem__untranslated">${escapeHtml(t('poems.notTranslated'))}</span>`;
    }

    return `
      <article class="poem" id="${escapeAttr(p.id)}"${p.tr ? ` lang="${lang}"` : ''}>
        <header class="poem__head">
          <h2 class="poem__title">${escapeHtml(shown.title || '')}</h2>
          <div class="poem__meta">
            ${langControl}
            <span class="poem__date">${escapeHtml(p.date || '')}</span>
          </div>
        </header>
        <div class="poem__body">
          <pre class="poem__text">${escapeHtml(shown.text || '')}</pre>
        </div>
      </article>
    `;
  }

  // "Original" / "Translation" switch on translated poems
  function wireLangToggles(poems){
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
        article.querySelector('.poem__text').textContent = shown.text || '';
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
