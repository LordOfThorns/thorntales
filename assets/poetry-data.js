/* poetry-data.js — shared loading for poetry pages.
   Needs auth.js before it (optional: without Supabase the files are used as they are).

     PoetryData.catalogs()        // catalogs.json with site edits (catalog_edits) applied
     PoetryData.hiddenPoems()     // Map catalog -> Set of poem ids the reader may not see
     PoetryData.restrictedTags()  // admin only: Map tag -> mode, for tags that are not public
     PoetryData.descToText(desc)  // "a<br>b" -> "a\nb" (for the editor)
     PoetryData.uploadCatalogImage(cat, file) -> public url
*/
(function () {
  'use strict';

  const db = () => (window.Auth ? Auth.client : null);

  async function catalogs() {
    const [json, rows] = await Promise.all([
      fetch('./data/poetry/catalogs.json', { cache: 'no-store' }).then((r) => {
        if (!r.ok) throw new Error(`catalogs.json: ${r.status}`);
        return r.json();
      }),
      select('catalog_edits', 'catalog, label, descr, image_url'),
    ]);
    const edits = new Map(rows.map((r) => [r.catalog, r]));
    return (json.catalogs || []).map((c) => {
      const e = edits.get(c.id);
      if (!e) return { ...c };
      // The site edit is written in Russian: it replaces the base fields only,
      // translations from catalogs.json (c.i18n) still win in their languages
      return {
        ...c,
        label: e.label,
        desc: e.descr.replace(/\r\n/g, '\n').split('\n').join('<br>'),
        image: e.image_url || '',
        edited: true,
      };
    });
  }

  async function hiddenPoems() {
    const map = new Map();
    if (!db()) return map;
    if (window.Auth) await Auth.ready; // the answer depends on who is signed in
    const { data, error } = await db().rpc('hidden_poems');
    if (error) {
      console.warn('hidden poems failed', error);
      return map;
    }
    for (const r of data || []) {
      if (!map.has(r.catalog)) map.set(r.catalog, new Set());
      map.get(r.catalog).add(r.poem_id);
    }
    return map;
  }

  async function restrictedTags() {
    const rows = await select('tag_access', 'tag, mode');
    return new Map(rows.filter((r) => r.mode !== 'public').map((r) => [r.tag, r.mode]));
  }

  async function select(table, cols) {
    if (!db()) return [];
    const { data, error } = await db().from(table).select(cols);
    if (error) {
      console.warn(`${table} failed`, error);
      return [];
    }
    return data || [];
  }

  function descToText(desc) {
    return String(desc || '').replace(/<br\s*\/?>/gi, '\n');
  }

  // Picture for a collection card: keeps its proportions, scaled down to fit
  // 240x128 (twice the size shown on the card), WebP with transparency kept.
  // Old pictures of this collection are removed after the new one is stored.
  async function uploadCatalogImage(cat, file) {
    const blob = await fitResize(file, 240, 128);
    const path = `${cat}/${Date.now()}.webp`;
    const bucket = db().storage.from('catalog-images');
    const up = await bucket.upload(path, blob, { contentType: 'image/webp', upsert: false });
    if (up.error) throw up.error;
    await removeCatalogImages(cat, path);
    return bucket.getPublicUrl(path).data.publicUrl;
  }

  async function removeCatalogImages(cat, keepPath) {
    const bucket = db().storage.from('catalog-images');
    const { data } = await bucket.list(cat, { limit: 100 });
    const old = (data || []).map((f) => `${cat}/${f.name}`).filter((p) => p !== keepPath);
    if (old.length) await bucket.remove(old);
  }

  async function fitResize(file, maxW, maxH) {
    const bitmap = await createImageBitmap(file);
    const k = Math.min(1, maxW / bitmap.width, maxH / bitmap.height);
    const canvas = document.createElement('canvas');
    canvas.width = Math.max(1, Math.round(bitmap.width * k));
    canvas.height = Math.max(1, Math.round(bitmap.height * k));
    const ctx = canvas.getContext('2d');
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(bitmap, 0, 0, canvas.width, canvas.height);
    bitmap.close();
    return new Promise((resolve, reject) => {
      canvas.toBlob((b) => (b ? resolve(b) : reject(new Error('canvas'))), 'image/webp', 0.9);
    });
  }


  window.PoetryData = { catalogs, hiddenPoems, restrictedTags, descToText, uploadCatalogImage, removeCatalogImages };
})();
