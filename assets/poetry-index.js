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
          <a class="poetry-card" href="${escapeHtml(c.href)}" style="${paperStyle()}">
            <span class="poetry-card__paper" aria-hidden="true"></span>
            <span class="poetry-card__fold" aria-hidden="true"></span>
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

  // ---------- Old paper ----------
  //
  // Every card is a sheet of old paper with its own past: a slight tilt, a shade of
  // yellow-grey, a few foxing spots and stains, edges that are not quite straight,
  // sometimes a torn edge or a folded corner. Returned as inline CSS variables.
  function paperStyle() {
    const rnd = (a, b) => a + Math.random() * (b - a);
    const f = (v) => v.toFixed(1);

    // edges: points along each side pushed in by a pixel or two; one side may be torn
    const torn = Math.random() < 0.4 ? Math.floor(rnd(0, 4)) : -1; // 0 top, 1 right, 2 bottom, 3 left
    const fold = Math.random() < 0.45 ? Math.round(rnd(14, 26)) : 0; // folded top-right corner
    const jag = (side) => (side === torn ? rnd(0, 5) : rnd(0, 1.6));
    const pts = [];
    const N = 9;
    for (let i = 0; i <= N; i++) { // top: left to right
      const x = (i / N) * 100;
      if (fold && i === N) { pts.push(`calc(100% - ${fold}px) ${f(jag(0))}px`); continue; }
      pts.push(`${f(x)}% ${f(jag(0))}px`);
    }
    if (fold) pts.push(`calc(100% - ${f(jag(1))}px) ${fold}px`);
    for (let i = 1; i <= N; i++) pts.push(`calc(100% - ${f(jag(1))}px) ${f((i / N) * 100)}%`);       // right
    for (let i = N - 1; i >= 0; i--) pts.push(`${f((i / N) * 100)}% calc(100% - ${f(jag(2))}px)`);  // bottom
    for (let i = N - 1; i >= 1; i--) pts.push(`${f(jag(3))}px ${f((i / N) * 100)}%`);               // left

    // stains and foxing spots: a few soft brown blots in random places
    const spots = [];
    const n = Math.round(rnd(1, 4));
    for (let i = 0; i < n; i++) {
      const s = rnd(4, 26), a = rnd(0.06, 0.16);
      spots.push(`radial-gradient(circle at ${f(rnd(5, 95))}% ${f(rnd(5, 95))}%, rgba(120,80,40,${a.toFixed(2)}) 0, rgba(120,80,40,0) ${f(s)}px)`);
    }
    if (Math.random() < 0.3) spots.push(coffeeRing(rnd, f));

    return [
      `--rot:${f(rnd(-1.8, 1.8))}deg`,
      `--tone:${rnd(0.82, 1).toFixed(2)}`,
      `--edge:polygon(${pts.join(',')})`,
      `--spots:${spots.join(',')}`,
      `--fold:${fold}px`,
      `--px:${Math.round(rnd(0, 240))}px`,
      `--py:${Math.round(rnd(0, 240))}px`,
    ].join(';');
  }

  // A coffee ring left by a cup, drawn as a small SVG (one background layer): a wobbly,
  // slightly oval outline, darkest where the coffee dried at the rim, a soft brown bleed
  // just inside it, a faint stain in the middle; now and then the rim fades out on one
  // side, as when the cup was tilted.
  function coffeeRing(rnd, f) {
    const r = rnd(18, 32), ry = r * rnd(0.88, 1.04);
    const pad = 4, w = 2 * r + 2 * pad, h = 2 * ry + 2 * pad, cx = w / 2, cy = h / 2;
    // wobble: a few slow waves round the circle with random phases
    const waves = [2, 3, 5].map((k) => [k, rnd(0, Math.PI * 2), rnd(0.02, 0.06) / (k / 2)]);
    const pts = [];
    for (let i = 0; i < 64; i++) {
      const t = (i / 64) * Math.PI * 2;
      let m = 1;
      for (const [k, ph, amp] of waves) m += Math.sin(t * k + ph) * amp;
      pts.push([cx + Math.cos(t) * r * m, cy + Math.sin(t) * ry * m]);
    }
    const d = `M${pts.map(([x, y]) => `${f(x)} ${f(y)}`).join('L')}Z`;
    const a = rnd(0.75, 1);
    const brown = (al) => `rgba(122,68,18,${(al * a).toFixed(2)})`;
    // a rim fading out on one side: a dash pattern that leaves a soft gap
    const len = Math.PI * (r + ry);
    const gap = Math.random() < 0.45 ? len * rnd(0.12, 0.3) : 0;
    const dash = gap ? ` stroke-dasharray='${f(len - gap)} ${f(gap)}' stroke-dashoffset='${f(rnd(0, len))}'` : '';
    const svg = `<svg xmlns='http://www.w3.org/2000/svg' width='${f(w)}' height='${f(h)}' viewBox='0 0 ${f(w)} ${f(h)}'>` +
      `<path d='${d}' fill='rgba(160,98,38,${(0.1 * a).toFixed(2)})'/>` +
      `<path d='${d}' fill='none' stroke='${brown(0.16)}' stroke-width='4.5' transform='translate(${f(cx * 0.03)} ${f(cy * 0.03)}) scale(.97)'/>` +
      `<path d='${d}' fill='none' stroke='${brown(0.55)}' stroke-width='1.3' stroke-linejoin='round'${dash}/>` +
      `</svg>`;
    // single-quoted url() with ' escaped: the value goes into a double-quoted style attribute
    const uri = encodeURIComponent(svg).replace(/'/g, '%27');
    return `url('data:image/svg+xml,${uri}') ${f(rnd(10, 90))}% ${f(rnd(15, 85))}% / ${f(w)}px ${f(h)}px no-repeat`;
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
