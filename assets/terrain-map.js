/* terrain-map.js
   Procedural physical-ish terrain map with rivers (flow accumulation),
   depression filling, thermal erosion, and hillshade.
   No deps. Canvas2D.

   Public API:
     TerrainMap.render({ canvas, size, seed, seaLevel, riverThreshold })
*/

(function () {
  'use strict';

  // -------------------- Utilities --------------------

  function clamp(x, a, b) { return x < a ? a : (x > b ? b : x); }
  function lerp(a, b, t) { return a + (b - a) * t; }
  function smoothstep(t) { return t * t * (3 - 2 * t); }

  // Fast-ish seeded RNG (Mulberry32)
  function mulberry32(seed) {
    let t = seed >>> 0;
    return function () {
      t += 0x6D2B79F5;
      let r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Hash for integer coords -> [0,1)
  function hash2i(x, y, seed) {
    // Robert Jenkins-ish mix + seed
    let n = (x * 374761393 + y * 668265263) ^ (seed * 1442695041);
    n = (n ^ (n >>> 13)) * 1274126177;
    n = (n ^ (n >>> 16)) >>> 0;
    return n / 4294967296;
  }

  // Value noise (grid) with smooth interpolation
  function valueNoise2(x, y, seed) {
    const x0 = Math.floor(x), y0 = Math.floor(y);
    const x1 = x0 + 1, y1 = y0 + 1;
    const sx = smoothstep(x - x0);
    const sy = smoothstep(y - y0);

    const v00 = hash2i(x0, y0, seed);
    const v10 = hash2i(x1, y0, seed);
    const v01 = hash2i(x0, y1, seed);
    const v11 = hash2i(x1, y1, seed);

    const ix0 = lerp(v00, v10, sx);
    const ix1 = lerp(v01, v11, sx);
    return lerp(ix0, ix1, sy);
  }

  // Fractal Brownian Motion (fBm)
  function fbm(x, y, seed, octaves, lacunarity, gain) {
    let amp = 0.5;
    let freq = 1.0;
    let sum = 0.0;
    for (let i = 0; i < octaves; i++) {
      sum += amp * valueNoise2(x * freq, y * freq, seed + i * 1013);
      freq *= lacunarity;
      amp *= gain;
    }
    return sum; // ~[0,1]
  }

  // 2D index helpers
  function idx(x, y, n) { return y * n + x; }

  // Neighbors (8-connected)
  const N8 = [
    [-1, -1], [0, -1], [1, -1],
    [-1, 0],           [1, 0],
    [-1, 1],  [0, 1],  [1, 1]
  ];

  // -------------------- Terrain generation --------------------

function buildHeightmap(n, seed, continentsMin, continentsMax) {
  const h = new Float32Array(n * n);
  const rng = mulberry32(seed);

  // Сколько материков реально делаем (равномерно в диапазоне)
  const cmin = Math.max(1, continentsMin | 0);
  const cmax = Math.max(cmin, continentsMax | 0);
  const continents = cmin + ((rng() * (cmax - cmin + 1)) | 0);

  // Параметры "бугров" материков
  const centers = [];
  for (let i = 0; i < continents; i++) {
    // Не ставим центры слишком близко к краям, чтобы материки не были “обрезками”
    const cx = 0.15 + rng() * 0.70;
    const cy = 0.15 + rng() * 0.70;

    // Разный размер материка
    const sx = 0.10 + rng() * 0.18;
    const sy = 0.10 + rng() * 0.18;

    // Вес (влияет на “высоту” суши)
    const w = 0.8 + rng() * 0.9;

    centers.push({ cx, cy, sx, sy, w });
  }

  // Domain warp, чтобы убрать “слишком ровные” формы
  // Это дешёвый варп на базе шума.
  function warp(u, v) {
    const wx = fbm(u * 2.0, v * 2.0, seed + 900, 3, 2.0, 0.5) - 0.5;
    const wy = fbm(u * 2.0, v * 2.0, seed + 901, 3, 2.0, 0.5) - 0.5;
    return [u + wx * 0.08, v + wy * 0.08];
  }

  for (let y = 0; y < n; y++) {
    for (let x = 0; x < n; x++) {
      const u0 = x / (n - 1);
      const v0 = y / (n - 1);
      const [u, v] = warp(u0, v0);

      // --- Continent mask: сумма гауссиан ---
      let cont = 0;
      for (let k = 0; k < centers.length; k++) {
        const c = centers[k];
        const dx = (u - c.cx) / c.sx;
        const dy = (v - c.cy) / c.sy;
        const g = Math.exp(-(dx * dx + dy * dy)); // 0..1
        cont += g * c.w;
      }

      // Нормируем и делаем края материков более “острыми”
      cont = cont / (continents * 0.95);
      cont = clamp(cont, 0, 1);
      cont = Math.pow(cont, 1.35);

      // --- Shoreline noise: рваные берега ---
      const shore = fbm(u * 6.0, v * 6.0, seed + 111, 5, 2.0, 0.5);
      const shore2 = fbm(u * 18.0, v * 18.0, seed + 112, 4, 2.0, 0.5);
      const coast = clamp(cont + (shore - 0.5) * 0.22 + (shore2 - 0.5) * 0.10, 0, 1);

      // --- Mountains + detail ---
      const m = fbm(u * 8.0, v * 8.0, seed + 20, 6, 2.0, 0.5);
      const ridge = 1.0 - Math.abs(m * 2.0 - 1.0);
      const mountains = Math.pow(ridge, 2.2);

      const d = fbm(u * 24.0, v * 24.0, seed + 30, 5, 2.0, 0.5);

      // Высота: материк задаёт “сушу/воду”, горы добавляются на суше
      // Важно: больше НЕ делаем радиальный edge-falloff, иначе опять будет “один круглый материк”.
      let height = 0.78 * coast + 0.18 * d + 0.20 * mountains * coast;

      h[idx(x, y, n)] = height;
    }
  }

  // Normalize to [0,1]
  let min = Infinity, max = -Infinity;
  for (let i = 0; i < h.length; i++) { min = Math.min(min, h[i]); max = Math.max(max, h[i]); }
  const inv = 1 / (max - min + 1e-9);
  for (let i = 0; i < h.length; i++) h[i] = (h[i] - min) * inv;

  return h;
}

  // Thermal erosion: simple slope-limited relaxation
  function thermalErode(h, n, iterations, talus) {
    const tmp = new Float32Array(h.length);

    for (let it = 0; it < iterations; it++) {
      tmp.set(h);

      for (let y = 1; y < n - 1; y++) {
        for (let x = 1; x < n - 1; x++) {
          const i = idx(x, y, n);
          const hi = tmp[i];

          // Find steepest drop neighbor
          let maxDrop = 0;
          let tx = 0, ty = 0;

          for (let k = 0; k < 8; k++) {
            const nx = x + N8[k][0];
            const ny = y + N8[k][1];
            const j = idx(nx, ny, n);
            const drop = hi - tmp[j];
            if (drop > maxDrop) {
              maxDrop = drop;
              tx = nx; ty = ny;
            }
          }

          if (maxDrop > talus) {
            // Move a small fraction "downhill"
            const move = (maxDrop - talus) * 0.25;
            h[i] -= move;
            h[idx(tx, ty, n)] += move;
          }
        }
      }

      // Clamp
      for (let i = 0; i < h.length; i++) h[i] = clamp(h[i], 0, 1);
    }
  }

  // -------------------- Hydrology --------------------

  // Depression filling (priority-flood style). Simple and effective.
  function fillDepressions(h, n) {
    // We will raise interior pits so every cell can drain to boundary.
    // Algorithm: Priority queue seeded with boundary cells.
    // For simplicity: implement binary heap.

    const size = n * n;
    const visited = new Uint8Array(size);
    const out = new Float32Array(h); // filled heights

    // Min-heap of [height, index]
    const heapH = new Float32Array(size);
    const heapI = new Int32Array(size);
    let heapSize = 0;

    function heapPush(height, index) {
      let i = heapSize++;
      heapH[i] = height;
      heapI[i] = index;
      while (i > 0) {
        const p = (i - 1) >> 1;
        if (heapH[p] <= heapH[i]) break;
        // swap
        const th = heapH[p]; heapH[p] = heapH[i]; heapH[i] = th;
        const ti = heapI[p]; heapI[p] = heapI[i]; heapI[i] = ti;
        i = p;
      }
    }

    function heapPop() {
      const resH = heapH[0];
      const resI = heapI[0];
      heapSize--;
      if (heapSize > 0) {
        heapH[0] = heapH[heapSize];
        heapI[0] = heapI[heapSize];
        // sift down
        let i = 0;
        while (true) {
          const l = i * 2 + 1;
          const r = l + 1;
          if (l >= heapSize) break;
          let m = l;
          if (r < heapSize && heapH[r] < heapH[l]) m = r;
          if (heapH[i] <= heapH[m]) break;
          const th = heapH[i]; heapH[i] = heapH[m]; heapH[m] = th;
          const ti = heapI[i]; heapI[i] = heapI[m]; heapI[m] = ti;
          i = m;
        }
      }
      return [resH, resI];
    }

    function markAndPush(i) {
      visited[i] = 1;
      heapPush(out[i], i);
    }

    // Seed boundaries
    for (let x = 0; x < n; x++) {
      markAndPush(idx(x, 0, n));
      markAndPush(idx(x, n - 1, n));
    }
    for (let y = 1; y < n - 1; y++) {
      markAndPush(idx(0, y, n));
      markAndPush(idx(n - 1, y, n));
    }

    while (heapSize > 0) {
      const [ch, ci] = heapPop();
      const cx = ci % n;
      const cy = (ci / n) | 0;

      for (let k = 0; k < 8; k++) {
        const nx = cx + N8[k][0];
        const ny = cy + N8[k][1];
        if (nx < 0 || ny < 0 || nx >= n || ny >= n) continue;
        const ni = idx(nx, ny, n);
        if (visited[ni]) continue;
        visited[ni] = 1;

        // If neighbor is lower than current, raise it to current (spill point)
        if (out[ni] < ch) out[ni] = ch;
        heapPush(out[ni], ni);
      }
    }

    return out;
  }

  // Flow direction: for each cell pick the lowest neighbor (D8)
  function computeFlowDir(h, n) {
    // store downstream index, or -1 if none
    const down = new Int32Array(n * n);
    down.fill(-1);

    for (let y = 1; y < n - 1; y++) {
      for (let x = 1; x < n - 1; x++) {
        const i = idx(x, y, n);
        const hi = h[i];

        let best = -1;
        let bestH = hi;

        for (let k = 0; k < 8; k++) {
          const nx = x + N8[k][0];
          const ny = y + N8[k][1];
          const j = idx(nx, ny, n);
          const hj = h[j];
          if (hj < bestH) {
            bestH = hj;
            best = j;
          }
        }

        down[i] = best; // may be -1 on flat/peak
      }
    }

    return down;
  }

  // Flow accumulation using topological order by height (sort indices by height desc)
  function computeFlowAccum(h, down) {
    const nCells = h.length;
    const order = new Int32Array(nCells);
    for (let i = 0; i < nCells; i++) order[i] = i;

    // JS sort on typed array is annoying: convert to normal array for order indices.
    // This is O(N log N). For 512-1024, OK.
    const ord = Array.from(order);
    ord.sort((a, b) => h[b] - h[a]); // descending heights

    const acc = new Float32Array(nCells);
    // Each cell contributes 1 unit of "rain"
    for (let i = 0; i < nCells; i++) acc[i] = 1;

    for (let k = 0; k < ord.length; k++) {
      const i = ord[k];
      const d = down[i];
      if (d >= 0) acc[d] += acc[i];
    }

    return acc;
  }

  // -------------------- Rendering --------------------

  function hillshade(h, n, x, y) {
    // simple gradient-based shade
    const xm1 = clamp(x - 1, 0, n - 1);
    const xp1 = clamp(x + 1, 0, n - 1);
    const ym1 = clamp(y - 1, 0, n - 1);
    const yp1 = clamp(y + 1, 0, n - 1);

    const dzdx = h[idx(xp1, y, n)] - h[idx(xm1, y, n)];
    const dzdy = h[idx(x, yp1, n)] - h[idx(x, ym1, n)];

    // Light from NW-ish
    const lx = -0.6, ly = -0.6, lz = 0.7;
    // normal approx
    let nx = -dzdx, ny = -dzdy, nz = 1.0;
    const invLen = 1 / Math.sqrt(nx * nx + ny * ny + nz * nz);
    nx *= invLen; ny *= invLen; nz *= invLen;

    let dot = nx * lx + ny * ly + nz * lz;
    dot = clamp(dot, 0, 1);
    return dot;
  }

  function colorFor(height, seaLevel) {
    // height in [0,1]
    // returns [r,g,b] 0..255
    if (height <= seaLevel) {
      // Ocean gradient
      const t = height / (seaLevel + 1e-9);
      const deep = [10, 30, 70];
      const shallow = [40, 120, 180];
      return [
        (deep[0] + (shallow[0] - deep[0]) * t) | 0,
        (deep[1] + (shallow[1] - deep[1]) * t) | 0,
        (deep[2] + (shallow[2] - deep[2]) * t) | 0
      ];
    }

    const t = (height - seaLevel) / (1 - seaLevel + 1e-9);

    // Land bands
    if (t < 0.18) return [70, 120, 60];      // lowlands
    if (t < 0.40) return [60, 140, 70];      // plains
    if (t < 0.62) return [90, 130, 80];      // hills
    if (t < 0.80) return [120, 120, 105];    // mountains
    return [235, 235, 240];                  // snow
  }

  function drawMap(ctx, h, acc, n, seaLevel, riverThreshold) {
    const img = ctx.createImageData(n, n);
    const data = img.data;

    // Normalize acc for visuals
    let accMax = 1;
    for (let i = 0; i < acc.length; i++) if (acc[i] > accMax) accMax = acc[i];

    for (let y = 0; y < n; y++) {
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        const height = h[i];

        let [r, g, b] = colorFor(height, seaLevel);

        // Shade (makes terrain readable)
        const sh = hillshade(h, n, x, y);
        const shade = 0.55 + 0.45 * sh;
        r = (r * shade) | 0;
        g = (g * shade) | 0;
        b = (b * shade) | 0;

        // Rivers: only on land, and only where accumulation is big enough
        if (height > seaLevel && acc[i] > riverThreshold) {
          // thickness grows with log(acc)
          const t = Math.log(acc[i]) / Math.log(accMax + 1e-9);
          const w = clamp((t - 0.25) * 2.0, 0, 1); // 0..1
          // Blend towards river blue
          const rr = 30, rg = 90, rb = 160;
          r = lerp(r, rr, 0.55 + 0.35 * w) | 0;
          g = lerp(g, rg, 0.55 + 0.35 * w) | 0;
          b = lerp(b, rb, 0.55 + 0.35 * w) | 0;
        }

        const o = i * 4;
        data[o + 0] = r;
        data[o + 1] = g;
        data[o + 2] = b;
        data[o + 3] = 255;
      }
    }

    ctx.putImageData(img, 0, 0);
  }

  // -------------------- Public API --------------------

  function render(opts) {
    const canvas = opts.canvas;
    if (!canvas) throw new Error('TerrainMap.render: canvas is required');
    const n = opts.size || 768;
    const seed = (opts.seed ?? 1337) | 0;
    const seaLevel = clamp(opts.seaLevel ?? 0.48, 0.05, 0.95);
    const riverThreshold = opts.riverThreshold ?? 1800;

    // Prepare canvas for crisp pixels
    canvas.width = n;
    canvas.height = n;

    const ctx = canvas.getContext('2d', { alpha: false });

    // 1) Height
    let h = buildHeightmap(n, seed, opts.continentsMin ?? 2, opts.continentsMax ?? 6);

    // 2) Erosion (light)
    thermalErode(h, n, 6, 0.018);

    // 3) Fill depressions for drainage correctness
    const hf = fillDepressions(h, n);

    // 4) Flow + accumulation
    const down = computeFlowDir(hf, n);
    const acc = computeFlowAccum(hf, down);

    // 5) Draw
    drawMap(ctx, h, acc, n, seaLevel, riverThreshold);
  }

  window.TerrainMap = { render };
})();