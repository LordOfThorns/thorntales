/* admin.js — the author's page: admin.html.

   "Download full backup" builds the site as it is now — files of the repository with all
   site edits from Supabase applied — and saves it as one zip:

     site/data/poetry/<cat>.txt, <cat>.<lang>.txt, catalogs.json, images/<cat>.<ext>
     site/data/texts/<section>.json
         what the public repository should contain: edits baked in, poems with a
         restricted tag (tag_access mode is not 'public') and hidden texts left out.
         tools/apply-backup.cmd copies this folder over data/ in the repository.
     full/data/...      the same with everything, hidden poems and texts included
     full/supabase/<table>.json   raw rows of every table (comments, profiles, ...)

   A poem saved on the site is HTML; in a txt file it gets "format: html" unless it is
   plain text (only line breaks). Poems with a restricted tag that exist only in a txt
   file are first saved to poem_edits, so they stay on the site once the file drops them.
*/
(async function () {
  'use strict';

  const t = (key, vars) => (window.I18n ? I18n.t(key, vars) : key);
  const sub = document.getElementById('adminSub');
  const page = document.getElementById('adminPage');
  const logHost = document.getElementById('backupLog');
  const button = document.getElementById('backupBtn');

  const SECTIONS = ['stories', 'tales', 'repatriate', 'dreams', 'chatter', 'other'];
  const LANGS = ['ru', 'en', 'de'];

  if (!window.Auth) { sub.textContent = t('poems.loadFailed'); return; }
  const db = Auth.client;
  await Auth.ready;
  if (!Auth.profile?.is_admin) { sub.textContent = t('tags.adminOnly'); return; }
  sub.textContent = t('admin.sub');
  page.hidden = false;

  button.addEventListener('click', async () => {
    button.disabled = true;
    logHost.textContent = '';
    try {
      const { zip, name } = await buildBackup(log);
      saveBlob(zip, name);
      log(t('admin.done', { name }));
    } catch (e) {
      console.error(e);
      log(t('admin.failed') + ' ' + (e.message || e));
    } finally {
      button.disabled = false;
    }
  });

  function log(line) {
    logHost.textContent += line + '\n';
    logHost.scrollTop = logHost.scrollHeight;
  }

  // ---------- Building ----------

  async function buildBackup(log) {
    const files = []; // { path, data: string | Uint8Array }
    const add = (path, data) => files.push({ path, data });

    log(t('admin.stepTables'));
    const tables = {};
    for (const name of ['poem_edits', 'catalog_edits', 'text_edits', 'tag_access', 'tag_readers', 'comments', 'profiles']) {
      tables[name] = await selectAll(name);
    }
    // newer tables: a missing one (SQL not run yet) does not stop the backup
    for (const name of ['site_images', 'pn_volumes', 'pn_entries', 'pn_revisions']) {
      try { tables[name] = await selectAll(name); } catch (e) { log(`${name}: ${e.message}`); }
    }
    const restricted = new Set(tables.tag_access.filter((r) => r.mode !== 'public').map((r) => r.tag));
    const isHidden = (tags) => (tags || []).some((tag) => restricted.has(tag));

    // ----- poetry -----
    log(t('admin.stepPoetry'));
    const catalogsJson = await fetchJson('./data/poetry/catalogs.json');
    const catEdits = new Map(tables.catalog_edits.map((r) => [r.catalog, r]));
    const editsByCat = new Map();
    for (const r of tables.poem_edits) {
      if (!editsByCat.has(r.catalog)) editsByCat.set(r.catalog, new Map());
      editsByCat.get(r.catalog).set(r.poem_id, r);
    }

    let saved = 0;
    for (const c of catalogsJson.catalogs || []) {
      const raw = (await fetchText(`./data/poetry/${encodeURIComponent(c.id)}.txt`)) || '';
      const edits = editsByCat.get(c.id) || new Map();
      const filePoems = PoetryData.parseTxt(raw).map(fromFile);

      // A restricted poem only in the file would disappear from the site: save it first
      for (const p of filePoems) {
        if (isHidden(p.tags) && !edits.has(p.id)) {
          const row = { catalog: c.id, poem_id: p.id, title: p.title || '…', date: p.date,
            body: p.format === 'html' ? p.text : escapeHtml(p.text), tags: p.tags };
          const { error } = await db.from('poem_edits').upsert(row);
          if (error) throw new Error(`${c.id}/${p.id}: ${error.message}`);
          edits.set(p.id, { ...row, tr: {} });
          tables.poem_edits.push(row);
          saved++;
        }
      }

      const poems = mergePoems(filePoems, edits, c.id);
      add(`full/data/poetry/${c.id}.txt`, serializePoems(poems));
      add(`site/data/poetry/${c.id}.txt`, serializePoems(poems.filter((p) => !isHidden(p.tags))));

      // Translations: <cat>.<lang>.txt from the file, site translations (tr) win
      for (const l of LANGS) {
        const trRaw = await fetchText(`./data/poetry/${encodeURIComponent(c.id)}.${l}.txt`);
        const byId = new Map(PoetryData.parseTxt(trRaw || '').map(fromFile).map((p) => [p.id, p]));
        for (const p of poems) {
          const e = edits.get(p.id)?.tr?.[l];
          if (e && (e.title || plain(e.body))) {
            byId.set(p.id, toPoem({ id: p.id, date: p.date, title: e.title || '', html: e.body, tags: [] }));
          }
        }
        if (!byId.size) continue;
        const order = new Map(poems.map((p, i) => [p.id, i]));
        const list = [...byId.values()].sort((a, b) => (order.get(a.id) ?? 1e9) - (order.get(b.id) ?? 1e9));
        const hiddenIds = new Set(poems.filter((p) => isHidden(p.tags)).map((p) => p.id));
        add(`full/data/poetry/${c.id}.${l}.txt`, serializePoems(list));
        add(`site/data/poetry/${c.id}.${l}.txt`, serializePoems(list.filter((p) => !hiddenIds.has(p.id))));
      }
    }
    if (saved) log(t('admin.savedHidden', { n: saved }));

    // Collection names, descriptions and pictures edited on the site
    const cats = [];
    for (const c of catalogsJson.catalogs || []) {
      const e = catEdits.get(c.id);
      if (!e) { cats.push(c); continue; }
      const out = { ...c, label: e.label, desc: e.descr.replace(/\r\n/g, '\n').split('\n').join('<br>') };
      if (e.image_url) {
        const img = await fetchBytes(e.image_url);
        if (img) {
          const ext = (e.image_url.match(/\.(webp|png|jpe?g)(\?|$)/i) || [0, 'webp'])[1].toLowerCase();
          const path = `data/poetry/images/${c.id}.${ext}`;
          add(`site/${path}`, img);
          add(`full/${path}`, img);
          out.image = `./${path}`;
        }
      } else {
        delete out.image;
      }
      cats.push(out);
    }
    const catsText = JSON.stringify({ ...catalogsJson, catalogs: cats }, null, 2) + '\n';
    add('site/data/poetry/catalogs.json', catsText);
    add('full/data/poetry/catalogs.json', catsText);

    // ----- texts -----
    log(t('admin.stepTexts'));
    for (const s of SECTIONS) {
      const file = (await fetchJson(`./data/texts/${s}.json`, true)) || { texts: [] };
      const byId = new Map((file.texts || []).map((x) => [x.id, { ...x }]));
      for (const r of tables.text_edits.filter((r) => r.section === s)) {
        const i18n = {};
        for (const [l, v] of Object.entries(r.tr || {})) {
          if (v && (v.title || plain(v.body))) i18n[l] = { title: v.title || '', html: v.body || '' };
        }
        const x = { id: r.text_id, title: r.title, date: r.date || '', lang: r.lang || 'ru', html: r.body || '' };
        if (Object.keys(i18n).length) x.i18n = i18n;
        if (r.deleted) x.hidden = true;
        byId.set(r.text_id, x);
      }
      const all = [...byId.values()].sort((a, b) => (a.date || '').localeCompare(b.date || ''));
      if (!all.length) continue;
      const pub = all.filter((x) => !x.hidden);
      add(`full/data/texts/${s}.json`, JSON.stringify({ texts: all }, null, 2) + '\n');
      if (pub.length) add(`site/data/texts/${s}.json`, JSON.stringify({ texts: pub }, null, 2) + '\n');
    }

    // ----- Polar Night: the copy the page shows when Supabase does not answer -----
    if (tables.pn_volumes && tables.pn_entries) {
      const vols = [...tables.pn_volumes].sort((a, b) => a.num - b.num);
      const ents = [...tables.pn_entries].sort((a, b) => a.volume_id.localeCompare(b.volume_id) || a.kind.localeCompare(b.kind) || a.num - b.num);
      const book = (v, e) => JSON.stringify({ volumes: v, entries: e }, null, 2) + '\n';
      add('full/data/polar-night/book.json', book(vols, ents));
      // the site copy is a public file: only what everyone may see
      const pubVols = vols.filter((v) => v.status === 'public');
      const pubIds = new Set(pubVols.map((v) => v.id));
      add('site/data/polar-night/book.json', book(pubVols, ents.filter((e) => e.status === 'public' && pubIds.has(e.volume_id))));
    }

    // ----- pictures replaced on the site (home page windows, ...) -----
    for (const r of tables.site_images || []) {
      const img = await fetchBytes(r.url);
      if (!img) continue;
      const ext = (r.url.match(/\.(webp|png|jpe?g)(\?|$)/i) || [0, 'webp'])[1].toLowerCase();
      add(`full/site-images/${r.key.replace(/[^a-z0-9_-]+/gi, '-')}.${ext}`, img);
    }

    // ----- raw tables -----
    for (const [name, rows] of Object.entries(tables)) {
      add(`full/supabase/${name}.json`, JSON.stringify(rows, null, 2) + '\n');
    }
    add('README.txt', readme());

    log(t('admin.stepZip', { n: files.length }));
    const stamp = new Date().toISOString().slice(0, 16).replace('T', '-').replace(':', '');
    return { zip: makeZip(files), name: `thorntales-backup-${stamp}.zip` };
  }

  // A poem from a txt block
  function fromFile({ meta, text }) {
    const tags = (meta.tags || '').split(/[,;]/).map((s) => s.trim().replace(/^#+/, '')).filter(Boolean);
    return { id: meta.id || '', date: meta.date || '', title: meta.title || '', tags,
      format: meta.format === 'html' ? 'html' : '', text, extra: extraMeta(meta) };
  }
  const KNOWN = new Set(['date', 'title', 'id', 'tags', 'format']);
  const extraMeta = (meta) => Object.entries(meta).filter(([k]) => !KNOWN.has(k));

  // A site row (HTML) as a txt poem: plain text when the HTML has nothing but line breaks
  function toPoem({ id, date, title, html, tags }) {
    const body = String(html || '').replace(/\r\n/g, '\n').replace(/\s+$/, '');
    const onlyBr = !/<(?!br\s*\/?>)[a-z\/!]/i.test(body);
    const text = onlyBr ? unescapeHtml(body.replace(/<br\s*\/?>/gi, '\n')) : body;
    return { id, date: date || '', title: title || '', tags: tags || [], format: onlyBr ? '' : 'html', text, extra: [] };
  }

  // File order with site edits applied, poems added on the site after the file ones;
  // the page sorts by date anyway, the file keeps its order for a readable diff
  function mergePoems(filePoems, edits, cat) {
    const seen = new Set();
    const list = filePoems.map((p) => {
      seen.add(p.id);
      const e = edits.get(p.id);
      return e ? { ...toPoem({ id: p.id, date: e.date, title: e.title, html: e.body, tags: e.tags }), extra: p.extra } : p;
    });
    const added = [...edits.values()].filter((e) => !seen.has(e.poem_id))
      .map((e) => toPoem({ id: e.poem_id, date: e.date, title: e.title, html: e.body, tags: e.tags }))
      .sort((a, b) => (a.date || '').localeCompare(b.date || ''));
    if (added.length) log(t('admin.added', { cat, n: added.length }));
    return list.concat(added);
  }

  function serializePoems(list) {
    return list.map((p) => {
      const head = [];
      if (p.date) head.push(`date: ${p.date}`);
      head.push(`title: ${p.title}`);
      if (p.id) head.push(`id: ${p.id}`);
      if (p.tags.length) head.push(`tags: ${p.tags.join(', ')}`);
      if (p.format) head.push(`format: ${p.format}`);
      for (const [k, v] of p.extra || []) head.push(`${k}: ${v}`);
      // an empty line inside the text would end the header early only before the text;
      // inside the text it is fine, the parser joins the rest back
      return `=== poem ===\n${head.join('\n')}\n\n${p.text.replace(/\s+$/, '')}\n`;
    }).join('\n');
  }

  function readme() {
    return [
      'Thorntales backup, ' + new Date().toISOString(),
      '',
      'site/  - the public data folder with all site edits applied, without hidden poems,',
      '         hidden texts and not public Polar Night. tools/apply-backup.cmd copies it over data/.',
      'full/  - everything: the same files with hidden poems and texts, pictures replaced on the',
      '         site (full/site-images) and every Supabase',
      '         table as JSON (full/supabase). Keep it private.',
      '',
    ].join('\n');
  }

  // ---------- Data helpers ----------

  async function selectAll(table) {
    const rows = [];
    for (let from = 0; ; from += 1000) {
      const { data, error } = await db.from(table).select('*').range(from, from + 999);
      if (error) throw new Error(`${table}: ${error.message}`);
      rows.push(...(data || []));
      if (!data || data.length < 1000) return rows;
    }
  }

  async function fetchText(url) {
    const r = await fetch(url, { cache: 'no-store' });
    return r.ok ? r.text() : null;
  }

  async function fetchJson(url, optional) {
    const r = await fetch(url, { cache: 'no-store' });
    if (!r.ok) { if (optional) return null; throw new Error(`${url}: ${r.status}`); }
    return r.json();
  }

  async function fetchBytes(url) {
    try {
      const r = await fetch(url, { cache: 'no-store' });
      return r.ok ? new Uint8Array(await r.arrayBuffer()) : null;
    } catch (e) {
      return null;
    }
  }

  const plain = (html) => (window.RichEdit ? RichEdit.plain(html || '') : String(html || '')).trim();

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;' }[ch]));
  }

  function unescapeHtml(s) {
    const el = document.createElement('textarea');
    el.innerHTML = s;
    return el.value;
  }

  // ---------- Zip (stored, no compression) ----------

  const CRC_TABLE = (() => {
    const tbl = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
      tbl[n] = c >>> 0;
    }
    return tbl;
  })();

  function crc32(bytes) {
    let c = 0xffffffff;
    for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
    return (c ^ 0xffffffff) >>> 0;
  }

  function makeZip(files) {
    const enc = new TextEncoder();
    const now = new Date();
    const dosTime = (now.getHours() << 11) | (now.getMinutes() << 5) | (now.getSeconds() >> 1);
    const dosDate = ((now.getFullYear() - 1980) << 9) | ((now.getMonth() + 1) << 5) | now.getDate();
    const chunks = [];
    const central = [];
    let offset = 0;
    for (const f of files) {
      const name = enc.encode(f.path);
      const data = typeof f.data === 'string' ? enc.encode(f.data) : f.data;
      const crc = crc32(data);
      const local = new DataView(new ArrayBuffer(30));
      local.setUint32(0, 0x04034b50, true);
      local.setUint16(4, 20, true);
      local.setUint16(6, 0x0800, true); // UTF-8 names
      local.setUint16(8, 0, true);      // stored
      local.setUint16(10, dosTime, true);
      local.setUint16(12, dosDate, true);
      local.setUint32(14, crc, true);
      local.setUint32(18, data.length, true);
      local.setUint32(22, data.length, true);
      local.setUint16(26, name.length, true);
      local.setUint16(28, 0, true);
      chunks.push(new Uint8Array(local.buffer), name, data);

      const cd = new DataView(new ArrayBuffer(46));
      cd.setUint32(0, 0x02014b50, true);
      cd.setUint16(4, 20, true);
      cd.setUint16(6, 20, true);
      cd.setUint16(8, 0x0800, true);
      cd.setUint16(10, 0, true);
      cd.setUint16(12, dosTime, true);
      cd.setUint16(14, dosDate, true);
      cd.setUint32(16, crc, true);
      cd.setUint32(20, data.length, true);
      cd.setUint32(24, data.length, true);
      cd.setUint16(28, name.length, true);
      cd.setUint32(42, offset, true);
      central.push(new Uint8Array(cd.buffer), name);
      offset += 30 + name.length + data.length;
    }
    const cdSize = central.reduce((s, c) => s + c.length, 0);
    const end = new DataView(new ArrayBuffer(22));
    end.setUint32(0, 0x06054b50, true);
    end.setUint16(8, files.length, true);
    end.setUint16(10, files.length, true);
    end.setUint32(12, cdSize, true);
    end.setUint32(16, offset, true);
    return new Blob([...chunks, ...central, new Uint8Array(end.buffer)], { type: 'application/zip' });
  }

  function saveBlob(blob, name) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
  }
})();
