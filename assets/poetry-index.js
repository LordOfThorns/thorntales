(async function () {
  const host = document.getElementById('catalogs');
  if (!host) return;

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const field = (obj, name) => (window.I18n ? I18n.field(obj, name) : obj[name]);

  let catalogs;
  try {
    catalogs = await PoetryData.catalogs();
    host.innerHTML = `
      <div class="poetry-grid">
        ${catalogs.map(c => {
          const desc = field(c, 'desc');
          return `
          <a class="poetry-card" href="${escapeHtml(c.href)}">
            <span class="poetry-card__count" data-cat="${escapeHtml(c.id)}"></span>
            <div class="poetry-card__title">${escapeHtml(field(c, 'label'))}</div>
            ${desc ? `<div class="poetry-card__desc">${allowOnlyBr(desc)}</div>` : ``}
            ${c.image ? `<img class="poetry-card__img" src="${escapeHtml(c.image)}" alt="" loading="lazy">` : ''}
          </a>`;
        }).join('')}
      </div>
    `;
  } catch (e) {
    host.innerHTML = `<p class="sub">${escapeHtml(t('poetry.catalogsFailed'))}: ${escapeHtml(String(e))}</p>`;
    return;
  }

  // Link to the tag access page, for the admin only
  if (window.Auth) {
    Auth.onChange(({ profile }) => {
      let link = document.getElementById('tagsLink');
      if (profile?.is_admin && !link) {
        host.insertAdjacentHTML('beforebegin',
          `<p class="catalog-admin" id="tagsLink"><a class="poem__edit-btn" href="/tags.html">${escapeHtml(t('catalog.tags'))}</a></p>`);
      } else if (!profile?.is_admin && link) {
        link.remove();
      }
    });
  }

  // Poem counts: ids from each txt file, minus poems hidden from this reader by tag
  // access, plus poems added on the site (poem_edits rows whose id is not in the file;
  // rows the reader may not see are not returned at all). Filled in after the cards.
  const siteP = loadSiteIds();
  const hiddenP = PoetryData.hiddenPoems();
  await Promise.all(catalogs.map(async (c) => {
    try {
      const res = await fetch(`./data/poetry/${encodeURIComponent(c.id)}.txt`);
      if (!res.ok) return;
      const ids = [];
      for (const block of (await res.text()).replace(/\r\n/g, '\n').split('=== poem ===\n').slice(1)) {
        const m = block.match(/^id:\s*(.+)$/m);
        ids.push(m ? m[1].trim() : '');
      }
      const fileIds = new Set(ids);
      const hidden = (await hiddenP).get(c.id) || new Set();
      const siteIds = (await siteP).get(c.id) || [];
      const n = ids.filter(id => !hidden.has(id)).length + siteIds.filter(id => !fileIds.has(id)).length;
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
    await Auth.ready; // restricted rows are visible only to readers who may see them
    const { data, error } = await Auth.client.from('poem_edits').select('catalog, poem_id');
    if (error) { console.warn('poem edits failed', error); return byCat; }
    for (const r of data) {
      if (!byCat.has(r.catalog)) byCat.set(r.catalog, []);
      byCat.get(r.catalog).push(r.poem_id);
    }
    return byCat;
  }

  // Descriptions may contain <br> only
  function allowOnlyBr(html){
    return escapeHtml(html).replace(/&lt;br\s*\/?&gt;/gi, '<br>');
  }

  function escapeHtml(s){
    return String(s ?? '').replace(/[&<>"']/g, ch => ({
      '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'
    }[ch]));
  }
})();
