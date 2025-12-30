(async function () {
  const host = document.getElementById('catalogs');
  if (!host) return;

  try {
    const res = await fetch('./data/poetry/catalogs.json', { cache: 'no-store' });
    if (!res.ok) throw new Error(`catalogs.json: ${res.status}`);
    const data = await res.json();

    host.innerHTML = `
      <div class="poetry-grid">
        ${data.catalogs.map(c => `
          <a class="poetry-card" href="${c.href}">
            <div class="poetry-card__title">${escapeHtml(c.label)}</div>
            ${c.desc ? `<div class="poetry-card__desc">${escapeHtml(c.desc)}</div>` : ``}
          </a>
        `).join('')}
      </div>
    `;
  } catch (e) {
    host.innerHTML = `<p class="sub">Не загрузились каталоги: ${escapeHtml(String(e))}</p>`;
  }

  function escapeHtml(s){
    return s.replace(/[&<>"']/g, ch => ({
      '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'
    }[ch]));
  }
})();
