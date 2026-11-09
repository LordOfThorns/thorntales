(async function () {
  const host = document.getElementById('catalogs');
  if (!host) return;

  const t = (key, fallback) => (window.I18n ? I18n.t(key) : fallback);
  const field = (obj, name) => (window.I18n ? I18n.field(obj, name) : obj[name]);

  try {
    const res = await fetch('./data/poetry/catalogs.json', { cache: 'no-store' });
    if (!res.ok) throw new Error(`catalogs.json: ${res.status}`);
    const data = await res.json();

    host.innerHTML = `
      <div class="poetry-grid">
        ${data.catalogs.map(c => {
          const desc = field(c, 'desc');
          return `
          <a class="poetry-card" href="${c.href}">
            <div class="poetry-card__title">${escapeHtml(field(c, 'label'))}</div>
            ${desc ? `<div class="poetry-card__desc">${desc}</div>` : ``}
          </a>`;
        }).join('')}
      </div>
    `;
  } catch (e) {
    host.innerHTML = `<p class="sub">${escapeHtml(t('poetry.catalogsFailed', 'Не загрузились каталоги'))}: ${escapeHtml(String(e))}</p>`;
  }

  function escapeHtml(s){
    return String(s ?? '').replace(/[&<>"']/g, ch => ({
      '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'
    }[ch]));
  }
})();
