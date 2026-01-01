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

	const catInfo = (catsJson.catalogs || []).find(x => x.id === cat);

	// Заголовки: label вместо "cat"
	const title = catInfo?.label || humanizeCat(cat);
	pageTitle.textContent = title;
	document.title = title;

	// Подзаголовок: desc (с разрешённым <br>)
	pageSub.innerHTML = catInfo?.desc
	  ? allowOnlyBr(catInfo.desc)
	  : '';

    const poems = parsePoemsTxt(raw);

    const title = humanizeCat(cat);
    pageTitle.textContent = title;
    document.title = title;

    // TOC
    toc.innerHTML = `
      <div class="poetry-toc__inner">
        ${poems.map(p => `
          <a class="poetry-toc__link" href="#${escapeAttr(p.id)}">
            <span class="poetry-toc__date">${escapeHtml(p.date || '')}</span>
            <span class="poetry-toc__title">${escapeHtml(p.title || '')}</span>
          </a>
        `).join('')}
      </div>
    `;

    // Poems
    poemsHost.innerHTML = poems.map(p => renderPoem(p)).join('');

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
