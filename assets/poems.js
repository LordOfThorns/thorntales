(async function () {
  const toc = document.getElementById('toc');
  const poemsHost = document.getElementById('poems');
  const pageTitle = document.getElementById('pageTitle');

  const cat = new URLSearchParams(location.search).get('cat') || 'lyrical';
  const url = `./data/poetry/${encodeURIComponent(cat)}.txt`;
  const pageSub = document.getElementById('pageSub');
  const catalogsUrl = './data/poetry/catalogs.json';

  try {
	  const res = await fetch(url, { cache: 'no-store' });
	  if (!res.ok) throw new Error(`${url}: ${res.status}`);
	  const raw = await res.text();

	  const catsRes = await fetch(catalogsUrl, { cache: 'no-store' });
	  if (!catsRes.ok) throw new Error(`${catalogsUrl}: ${catsRes.status}`);
	  const catsJson = await catsRes.json();

	  const catalogs = catsJson.catalogs || [];
	  const catIndex = catalogs.findIndex(x => x.id === cat);
	  const catInfo = catalogs[catIndex];

	  // Заголовки: label вместо "cat"
	  const pageLabel = catInfo?.label || humanizeCat(cat);
	  pageTitle.textContent = pageLabel;
	  document.title = `${pageLabel} · Thorntales`;
	  if (window.Layout) Layout.setCrumbs([{ label: pageLabel }]);

	  // Подзаголовок: desc (разрешаем только <br>)
	  pageSub.innerHTML = catInfo?.desc ? allowOnlyBr(catInfo.desc) : '';

	  const poems = parsePoemsTxt(raw);

	  // TOC: sidebar on wide screens, collapsible block on narrow ones
	  toc.innerHTML = `
		<details class="poetry-toc__details">
		  <summary class="poetry-toc__summary">Оглавление <span class="poetry-toc__count">${poems.length}</span></summary>
		  <div class="poetry-toc__inner">
		    ${poems.map(p => `
			  <a class="poetry-toc__link" href="#${escapeAttr(p.id)}">
			    <span class="poetry-toc__date">${escapeHtml(p.date || '')}</span>
			    <span class="poetry-toc__title">${escapeHtml(p.title || '')}</span>
			  </a>
		    `).join('')}
		  </div>
		</details>
	  `;
	  setupTocMode(toc.querySelector('details'));

	  // Poems
	  poemsHost.innerHTML = poems.map(p => renderPoem(p)).join('');

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
          <h2 class="poem__title">Не загрузилось</h2>
          <div class="poem__date">${escapeHtml(String(e))}</div>
        </header>
      </article>
    `;
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
        <span class="catalog-pager__hint">← предыдущий</span>
        <span class="catalog-pager__label">${escapeHtml(prev.label)}</span></a>` : '<span></span>'}
      <a class="catalog-pager__link catalog-pager__link--all" href="/poetry.html">Все каталоги</a>
      ${next ? `<a class="catalog-pager__link catalog-pager__link--next" href="${escapeAttr(next.href)}">
        <span class="catalog-pager__hint">следующий →</span>
        <span class="catalog-pager__label">${escapeHtml(next.label)}</span></a>` : '<span></span>'}
    `;
  }

  function allowOnlyBr(html){
  // 1) экранируем всё
  const safe = escapeHtml(html);
  // 2) возвращаем только <br> обратно (и <br/>)
  return safe.replace(/&lt;br\s*\/?&gt;/gi, '<br>');
}

  function parsePoemsTxt(raw) {
    // Нормализуем переводы строк
    const text = raw.replace(/\r\n/g, '\n');

    // Режем на блоки по маркеру
    const chunks = text.split('\n=== poem ===\n');

    // Если файл начинается с маркера, то первый chunk будет содержать "=== poem ===" без ведущего \n
    // Нормализуем: выкинем всё до первого маркера
    let blocks = chunks;
    if (!text.startsWith('=== poem ===')) {
      // если кто-то решил вставить комментарий сверху, то blocks уже ок
    } else {
      // startsWith: первый блок содержит "=== poem ===" в начале без split, поэтому поправим
      blocks = text.split('=== poem ===\n').slice(1);
    }

    const poems = [];

    for (const blockRaw of blocks) {
      const block = blockRaw.trim();
      if (!block) continue;

      // header до первой пустой строки
      const parts = block.split('\n\n');
      const headerPart = parts[0] || '';
      const bodyPart = parts.slice(1).join('\n\n') || '';

      const meta = {};
      for (const line of headerPart.split('\n')) {
        const m = line.match(/^([a-zA-Z0-9_-]+)\s*:\s*(.*)$/);
        if (!m) continue;
        meta[m[1].toLowerCase()] = m[2].trim();
      }

      const title = meta.title || 'Без названия';
      const date = meta.date || '';
      const id = meta.id || makeId(date, title);

      poems.push({
        id,
        date,
        title,
        text: bodyPart
      });
    }

    // Сортировка по дате (если даты ISO)
    poems.sort((a, b) => (a.date || '').localeCompare(b.date || ''));

    return poems;
  }

  function renderPoem(p){
    const lines = typeof p.text === 'string' ? p.text.split('\n') : [];
    return `
      <article class="poem" id="${escapeAttr(p.id)}">
        <header class="poem__head">
          <h2 class="poem__title">${escapeHtml(p.title || '')}</h2>
          <div class="poem__date">${escapeHtml(p.date || '')}</div>
        </header>
        <div class="poem__body">
          <pre class="poem__text">${escapeHtml(lines.join('\n'))}</pre>
        </div>
      </article>
    `;
  }

  function makeId(date, title){
    const base = `${date}-${title}`.trim() || title || 'poem';
    return base
      .toLowerCase()
      .replace(/ё/g, 'е')
      .replace(/[^a-z0-9\u0400-\u04ff]+/gi, '-')
      .replace(/^-+|-+$/g, '')
      .slice(0, 80);
  }

  function humanizeCat(cat){
    // можешь сделать карту названий, если хочешь
    return cat === 'lyrical' ? 'Лирика' : `Стихи: ${cat}`;
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
