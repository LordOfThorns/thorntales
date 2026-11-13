(async function () {
  const host = document.getElementById('catalogs');
  if (!host) return;

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const field = (obj, name) => (window.I18n ? I18n.field(obj, name) : obj[name]);

  let catalogs;
  try {
    const res = await fetch('./data/poetry/catalogs.json', { cache: 'no-store' });
    if (!res.ok) throw new Error(`catalogs.json: ${res.status}`);
    catalogs = (await res.json()).catalogs;

    host.innerHTML = `
      <div class="poetry-grid">
        ${catalogs.map(c => {
          const desc = field(c, 'desc');
          return `
          <a class="poetry-card" href="${c.href}">
            <span class="poetry-card__count" data-cat="${escapeHtml(c.id)}"></span>
            <div class="poetry-card__title">${escapeHtml(field(c, 'label'))}</div>
            ${desc ? `<div class="poetry-card__desc">${desc}</div>` : ``}
          </a>`;
        }).join('')}
      </div>
    `;
  } catch (e) {
    host.innerHTML = `<p class="sub">${escapeHtml(t('poetry.catalogsFailed'))}: ${escapeHtml(String(e))}</p>`;
    return;
  }

  // Poem counts: ids from each txt file + poems added on the site (poem_edits rows
  // whose id is not in the file). Filled in after the cards are shown.
  const addedP = loadSiteIds();
  await Promise.all(catalogs.map(async (c) => {
    try {
      const res = await fetch(`./data/poetry/${encodeURIComponent(c.id)}.txt`);
      if (!res.ok) return;
      const ids = new Set();
      let blocks = 0;
      for (const block of (await res.text()).replace(/\r\n/g, '\n').split('=== poem ===\n').slice(1)) {
        blocks++;
        const m = block.match(/^id:\s*(.+)$/m);
        if (m) ids.add(m[1].trim());
      }
      const siteIds = (await addedP).get(c.id) || [];
      const n = blocks + siteIds.filter(id => !ids.has(id)).length;
      const el = host.querySelector(`.poetry-card__count[data-cat="${CSS.escape(c.id)}"]`);
      if (el) {
        el.textContent = n;
        el.title = t('poetry.count', { n });
      }
    } catch (e) {
      console.warn('poem count failed', c.id, e);
    }
  }));

  async function loadSiteIds() {
    const byCat = new Map();
    if (!window.Auth) return byCat;
    const { data, error } = await Auth.client.from('poem_edits').select('catalog, poem_id');
    if (error) { console.warn('poem edits failed', error); return byCat; }
    for (const r of data) {
      if (!byCat.has(r.catalog)) byCat.set(r.catalog, []);
      byCat.get(r.catalog).push(r.poem_id);
    }
    return byCat;
  }

  function escapeHtml(s){
    return String(s ?? '').replace(/[&<>"']/g, ch => ({
      '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'
    }[ch]));
  }
})();
