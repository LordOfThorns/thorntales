/* terrain-map.js
   Procedural terrain + rivers + climate/biomes + pseudo-tectonics (Voronoi plates).
   No deps. Canvas2D.

   Public API:
     TerrainMap.render({
       canvas, size, seed, seaLevel, riverThreshold,
       continentsMin, continentsMax,
       plates, windDir
     })
*/

(function () {
  'use strict';

  // -------------------- Utilities --------------------

  function clamp(x, a, b) { return x < a ? a : (x > b ? b : x); }
  function lerp(a, b, t) { return a + (b - a) * t; }
  function smoothstep(t) { return t * t * (3 - 2 * t); }
  function fract(x) { return x - Math.floor(x); }
  function length2(x, y) { return Math.sqrt(x * x + y * y); }

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

  // -------------------- Continents / base height --------------------

  function buildHeightmap(n, seed, continentsMin, continentsMax) {
    const h = new Float32Array(n * n);
    const rng = mulberry32(seed);

    const cmin = Math.max(1, continentsMin | 0);
    const cmax = Math.max(cmin, continentsMax | 0);
    const continents = cmin + ((rng() * (cmax - cmin + 1)) | 0);

    const centers = [];
    for (let i = 0; i < continents; i++) {
      const cx = 0.12 + rng() * 0.76;
      const cy = 0.12 + rng() * 0.76;
      const sx = 0.10 + rng() * 0.20;
      const sy = 0.10 + rng() * 0.20;
      const w = 0.8 + rng() * 1.0;
      centers.push({ cx, cy, sx, sy, w });
    }

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

        let cont = 0;
        for (let k = 0; k < centers.length; k++) {
          const c = centers[k];
          const dx = (u - c.cx) / c.sx;
          const dy = (v - c.cy) / c.sy;
          const g = Math.exp(-(dx * dx + dy * dy));
          cont += g * c.w;
        }

        cont = cont / (continents * 0.95);
        cont = clamp(cont, 0, 1);
        cont = Math.pow(cont, 1.35);

        const shore = fbm(u * 6.0, v * 6.0, seed + 111, 5, 2.0, 0.5);
        const shore2 = fbm(u * 18.0, v * 18.0, seed + 112, 4, 2.0, 0.5);
        const coast = clamp(cont + (shore - 0.5) * 0.22 + (shore2 - 0.5) * 0.10, 0, 1);

        const d = fbm(u * 22.0, v * 22.0, seed + 30, 5, 2.0, 0.5);
        let height = 0.82 * coast + 0.18 * d;

        h[idx(x, y, n)] = height;
      }
    }

    normalize01(h);
    return h;
  }

  function normalize01(a) {
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < a.length; i++) { if (a[i] < min) min = a[i]; if (a[i] > max) max = a[i]; }
    const inv = 1 / (max - min + 1e-9);
    for (let i = 0; i < a.length; i++) a[i] = (a[i] - min) * inv;
  }

  // -------------------- Pseudo-tectonics (Voronoi plates) --------------------

  function generatePlates(seed, platesCount) {
    const rng = mulberry32(seed + 777);
    const p = [];
    const count = clamp(platesCount | 0, 4, 40);

    for (let i = 0; i < count; i++) {
      const x = rng();
      const y = rng();
      const ang = rng() * Math.PI * 2;
      const spd = 0.25 + rng() * 0.85; // relative speed
      const vx = Math.cos(ang) * spd;
      const vy = Math.sin(ang) * spd;

      // oceanic vs continental tendency (affects baseline uplift a bit)
      const type = rng() < 0.55 ? 0 : 1; // 0 oceanic, 1 continental
      p.push({ x, y, vx, vy, type });
    }
    return p;
  }

  // For each cell: nearest plate id and second nearest distance for boundary strength.
  function computePlateFields(n, plates) {
    const ids = new Int16Array(n * n);
    const d2min = new Float32Array(n * n);
    const d2second = new Float32Array(n * n);

    for (let y = 0; y < n; y++) {
      const v = y / (n - 1);
      for (let x = 0; x < n; x++) {
        const u = x / (n - 1);

        let bestId = 0;
        let best = Infinity;
        let second = Infinity;

        for (let k = 0; k < plates.length; k++) {
          const dx = u - plates[k].x;
          const dy = v - plates[k].y;
          const d2 = dx * dx + dy * dy;
          if (d2 < best) {
            second = best;
            best = d2;
            bestId = k;
          } else if (d2 < second) {
            second = d2;
          }
        }

        const i = idx(x, y, n);
        ids[i] = bestId;
        d2min[i] = best;
        d2second[i] = second;
      }
    }

    return { ids, d2min, d2second };
  }

  // Build uplift map: mountains on convergent boundaries, rifts on divergent.
  function buildUpliftMap(n, plates, plateFields, seed) {
    const { ids, d2min, d2second } = plateFields;
    const uplift = new Float32Array(n * n);

    // boundary strength based on how close 1st/2nd nearest are (Voronoi edge)
    // smaller gap -> stronger boundary
    for (let y = 1; y < n - 1; y++) {
      for (let x = 1; x < n - 1; x++) {
        const i = idx(x, y, n);
        const idA = ids[i];

        // detect if near boundary: neighbour has different plate id
        let idB = idA;
        for (let k = 0; k < 8; k++) {
          const j = idx(x + N8[k][0], y + N8[k][1], n);
          if (ids[j] !== idA) { idB = ids[j]; break; }
        }

        if (idB === idA) continue;

        const A = plates[idA];
        const B = plates[idB];

        // boundary weight: 0..1
		const f1 = Math.sqrt(d2min[i]);
		const f2 = Math.sqrt(d2second[i]);
		const df = Math.max(1e-6, f2 - f1); // чем меньше, тем ближе к границе

		// ширина границы в "долях карты" (подбирается)
		const width = 0.020; // попробуй 0.015..0.030
		let bw = clamp(1.0 - (df / width), 0, 1);

		// сделаем профиль “горного пояса” более естественным
		bw = bw * bw; // или bw = Math.pow(bw, 2.0);

        if (bw <= 0) continue;

        // boundary normal approx: from A to B
        const px = x / (n - 1);
        const py = y / (n - 1);
        // use vector between plate centers as proxy for boundary normal
        let nx = B.x - A.x;
        let ny = B.y - A.y;
        const invn = 1 / (Math.sqrt(nx * nx + ny * ny) + 1e-9);
        nx *= invn; ny *= invn;

        // relative velocity
        const rvx = B.vx - A.vx;
        const rvy = B.vy - A.vy;

        // convergence positive if moving towards each other along normal
        const conv = -(rvx * nx + rvy * ny); // >0 convergent, <0 divergent

        // shear magnitude (transform faults)
        const shear = Math.abs(rvx * (-ny) + rvy * nx);

        // convert to uplift
        let u = 0;

        if (conv > 0) {
          // mountains
          u += conv * 0.9;
          // continental-continental more mountainy
          if (A.type === 1 && B.type === 1) u *= 1.15;
        } else {
          // rifts / trenches (negative)
          u += conv * 0.45; // conv is negative
          // oceanic divergence tends to create mid-ocean ridges, keep small negative
          if (A.type === 0 && B.type === 0) u *= 0.6;
        }

        // shear adds slight roughness
        u += shear * 0.18;

        // add some noise so boundaries aren’t perfect lines
        const nu = fbm(px * 10.0, py * 10.0, seed + 500, 3, 2.0, 0.5) - 0.5;
        u *= (0.85 + nu * 0.3);

        uplift[i] = u * bw;
      }
    }

    // blur a bit to spread ridges
    boxBlur(uplift, n, 4);
	boxBlur(uplift, n, 2);

    // normalise uplift into roughly [-1,1] (preserve sign)
    let maxAbs = 1e-6;
    for (let i = 0; i < uplift.length; i++) {
      const a = Math.abs(uplift[i]);
      if (a > maxAbs) maxAbs = a;
    }
    const inv = 1 / maxAbs;
    for (let i = 0; i < uplift.length; i++) uplift[i] *= inv;

    return uplift;
  }

  function boxBlur(a, n, radius) {
    if (radius <= 0) return;
    const tmp = new Float32Array(a.length);
    const r = radius;

    // horizontal
    for (let y = 0; y < n; y++) {
      let sum = 0;
      for (let x = -r; x <= r; x++) {
        const xx = clamp(x, 0, n - 1);
        sum += a[idx(xx, y, n)];
      }
      for (let x = 0; x < n; x++) {
        tmp[idx(x, y, n)] = sum / (r * 2 + 1);
        const xOut = x - r;
        const xIn = x + r + 1;
        if (xOut >= 0) sum -= a[idx(xOut, y, n)];
        if (xIn < n) sum += a[idx(xIn, y, n)];
      }
    }

    // vertical
    for (let x = 0; x < n; x++) {
      let sum = 0;
      for (let y = -r; y <= r; y++) {
        const yy = clamp(y, 0, n - 1);
        sum += tmp[idx(x, yy, n)];
      }
      for (let y = 0; y < n; y++) {
        a[idx(x, y, n)] = sum / (r * 2 + 1);
        const yOut = y - r;
        const yIn = y + r + 1;
        if (yOut >= 0) sum -= tmp[idx(x, yOut, n)];
        if (yIn < n) sum += tmp[idx(x, yIn, n)];
      }
    }
  }

  // -------------------- Erosion --------------------

  function thermalErode(h, n, iterations, talus) {
    const tmp = new Float32Array(h.length);

    for (let it = 0; it < iterations; it++) {
      tmp.set(h);

      for (let y = 1; y < n - 1; y++) {
        for (let x = 1; x < n - 1; x++) {
          const i = idx(x, y, n);
          const hi = tmp[i];

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
            const move = (maxDrop - talus) * 0.25;
            h[i] -= move;
            h[idx(tx, ty, n)] += move;
          }
        }
      }

      for (let i = 0; i < h.length; i++) h[i] = clamp(h[i], 0, 1);
    }
  }

  // -------------------- Hydrology --------------------

  function fillDepressions(h, n) {
    const size = n * n;
    const visited = new Uint8Array(size);
    const out = new Float32Array(h);

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

        if (out[ni] < ch) out[ni] = ch;
        heapPush(out[ni], ni);
      }
    }

    return out;
  }

  function computeFlowDir(h, n) {
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

        down[i] = best;
      }
    }

    return down;
  }

  function computeFlowAccum(h, down) {
    const nCells = h.length;
    const order = new Int32Array(nCells);
    for (let i = 0; i < nCells; i++) order[i] = i;

    const ord = Array.from(order);
    ord.sort((a, b) => h[b] - h[a]);

    const acc = new Float32Array(nCells);
    for (let i = 0; i < nCells; i++) acc[i] = 1;

    for (let k = 0; k < ord.length; k++) {
      const i = ord[k];
      const d = down[i];
      if (d >= 0) acc[d] += acc[i];
    }

    return acc;
  }

  // -------------------- Climate --------------------

  function computeTemperature(h, n, seaLevel) {
    const t = new Float32Array(n * n);
    for (let y = 0; y < n; y++) {
      const lat = Math.abs((y / (n - 1)) - 0.5) * 2.0; // 0 equator -> 1 poles
      // base temperature curve (nonlinear: wider tropics)
      const base = 1.0 - Math.pow(lat, 1.2);

      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        const height = h[i];

        // altitude cooling (only above sea)
        const alt = Math.max(0, (height - seaLevel) / (1 - seaLevel + 1e-9));
        let temp = base - alt * 0.65;

        // slight coastal moderation: ocean areas less extreme
        if (height <= seaLevel) temp = base * 0.92 + 0.08;

        t[i] = clamp(temp, 0, 1);
      }
    }
    return t;
  }

  // windDir: 'W', 'E', 'N', 'S' (dominant wind direction)
  function computeRain(h, n, seaLevel, windDir, desertStrength) {
    const rain = new Float32Array(n * n);

    // helpers to scan lines
    function scanWestToEast() {
      for (let y = 0; y < n; y++) {
        let m = 0.0; // moisture carried
        let prevH = h[idx(0, y, n)];
        for (let x = 0; x < n; x++) {
          const i = idx(x, y, n);
          const hi = h[i];

          if (hi <= seaLevel) {
            m = 1.0;
            rain[i] = 1.0;
            prevH = hi;
            continue;
          }

          // decay as we move inland
          m *= 0.985;

          // orographic precipitation: if slope upward, dump moisture
          const up = Math.max(0, hi - prevH);
          const drop = up * 3.0; // tune
          const r = clamp(m * (0.22 + drop), 0, 1);
          rain[i] = r;

          // rain shadow: after dumping, moisture decreases
          m = clamp(m - r * 0.55, 0, 1);

          prevH = hi;
        }
      }
    }

    function scanEastToWest() {
      for (let y = 0; y < n; y++) {
        let m = 0.0;
        let prevH = h[idx(n - 1, y, n)];
        for (let x = n - 1; x >= 0; x--) {
          const i = idx(x, y, n);
          const hi = h[i];

          if (hi <= seaLevel) {
            m = 1.0;
            rain[i] = 1.0;
            prevH = hi;
            continue;
          }

          m *= 0.985;
          const up = Math.max(0, hi - prevH);
          const drop = up * 3.0;
          const r = clamp(m * (0.22 + drop), 0, 1);
          rain[i] = r;
          m = clamp(m - r * 0.55, 0, 1);
          prevH = hi;
        }
      }
    }

    function scanNorthToSouth() {
      for (let x = 0; x < n; x++) {
        let m = 0.0;
        let prevH = h[idx(x, 0, n)];
        for (let y = 0; y < n; y++) {
          const i = idx(x, y, n);
          const hi = h[i];

          if (hi <= seaLevel) {
            m = 1.0;
            rain[i] = 1.0;
            prevH = hi;
            continue;
          }

          m *= 0.985;
          const up = Math.max(0, hi - prevH);
          const drop = up * 3.0;
          const r = clamp(m * (0.22 + drop), 0, 1);
          rain[i] = r;
          m = clamp(m - r * 0.55, 0, 1);
          prevH = hi;
        }
      }
    }

    function scanSouthToNorth() {
      for (let x = 0; x < n; x++) {
        let m = 0.0;
        let prevH = h[idx(x, n - 1, n)];
        for (let y = n - 1; y >= 0; y--) {
          const i = idx(x, y, n);
          const hi = h[i];

          if (hi <= seaLevel) {
            m = 1.0;
            rain[i] = 1.0;
            prevH = hi;
            continue;
          }

          m *= 0.985;
          const up = Math.max(0, hi - prevH);
          const drop = up * 3.0;
          const r = clamp(m * (0.22 + drop), 0, 1);
          rain[i] = r;
          m = clamp(m - r * 0.55, 0, 1);
          prevH = hi;
        }
      }
    }

    // initialise with chosen wind scan
    if (windDir === 'E') scanEastToWest();
    else if (windDir === 'N') scanNorthToSouth();
    else if (windDir === 'S') scanSouthToNorth();
    else scanWestToEast(); // default W

    // add latitude-driven dryness bands (Hadley cells-ish):
    for (let y = 0; y < n; y++) {
      const lat = Math.abs((y / (n - 1)) - 0.5) * 2.0;
      // deserts around ~0.3 lat, wetter near equator and subpolar
      const hadleyDry = Math.exp(-Math.pow((lat - 0.35) / 0.14, 2));
      const equatorWet = Math.exp(-Math.pow(lat / 0.22, 2));
      const polarDry = Math.exp(-Math.pow((lat - 1.0) / 0.22, 2));

		const ds = clamp(desertStrength ?? 0.35, 0, 1);

		const mod = clamp(
		  1.0 + equatorWet * 0.25
			  - hadleyDry * ds
			  - polarDry * 0.10,
		  0.45,
		  1.15
		);


      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        // oceans stay wet
        if (h[i] <= seaLevel) continue;
        rain[i] = clamp(rain[i] * mod, 0, 1);
      }
    }
	
	// --- desertStrength should visibly matter ---
const ds2 = clamp(desertStrength ?? 0.35, 0, 1);

// dryness noise field (patchy deserts, not a clean band)
const drySeed = 12345;

for (let y = 0; y < n; y++) {
  const v = y / (n - 1);
  for (let x = 0; x < n; x++) {
    const u = x / (n - 1);
    const i = idx(x, y, n);

    if (h[i] <= seaLevel) continue; // ocean unaffected

    // patchiness: 0..1, makes deserts cluster instead of uniform wash
    const patch = fbm(u * 3.5, v * 3.5, drySeed, 4, 2.0, 0.5);

    // continentality proxy: less rain where current rain is already low (inland/shadow)
    const inland = clamp(1.0 - rain[i], 0, 1);

    // reduce rain: stronger in already-dry places + patchy field
    const dryness = clamp(0.35 + 0.65 * inland, 0, 1) * (0.55 + 0.45 * patch);

    // apply
    const factor = 1.0 - ds2 * 0.65 * dryness; // 0.65 = strength knob
    rain[i] = clamp(rain[i] * factor, 0, 1);
  }
}


    return rain;
  }

  // -------------------- Rendering helpers --------------------

  function hillshade(h, n, x, y) {
    const xm1 = clamp(x - 1, 0, n - 1);
    const xp1 = clamp(x + 1, 0, n - 1);
    const ym1 = clamp(y - 1, 0, n - 1);
    const yp1 = clamp(y + 1, 0, n - 1);

    const dzdx = h[idx(xp1, y, n)] - h[idx(xm1, y, n)];
    const dzdy = h[idx(x, yp1, n)] - h[idx(x, ym1, n)];

    const lx = -0.6, ly = -0.6, lz = 0.7;
    let nx = -dzdx, ny = -dzdy, nz = 1.0;
    const invLen = 1 / Math.sqrt(nx * nx + ny * ny + nz * nz);
    nx *= invLen; ny *= invLen; nz *= invLen;

    let dot = nx * lx + ny * ly + nz * lz;
    dot = clamp(dot, 0, 1);
    return dot;
  }
  
  function heightColor(height, seaLevel) {
  if (height <= seaLevel) return oceanColor(height, seaLevel);

  const t = (height - seaLevel) / (1 - seaLevel + 1e-9);

  // Hypsometric tint: green -> yellow -> brown -> white
  if (t < 0.20) return [70, 140, 80];     // lowlands
  if (t < 0.45) return [145, 170, 90];    // plains / savanna-ish
  if (t < 0.65) return [170, 165, 110];   // high plains
  if (t < 0.82) return [140, 125, 105];   // mountains
  return [235, 235, 240];                 // high peaks / snow
}

  function oceanColor(height, seaLevel) {
    const t = height / (seaLevel + 1e-9);
    const deep = [10, 30, 70];
    const shallow = [40, 120, 180];
    return [
      (deep[0] + (shallow[0] - deep[0]) * t) | 0,
      (deep[1] + (shallow[1] - deep[1]) * t) | 0,
      (deep[2] + (shallow[2] - deep[2]) * t) | 0
    ];
  }

  function biomeColor(height, seaLevel, temp, rain) {
    if (height <= seaLevel) return oceanColor(height, seaLevel);

    const alt = (height - seaLevel) / (1 - seaLevel + 1e-9);

    // snow by altitude or cold
    if (temp < 0.14 || alt > 0.88) return [235, 235, 240];

    // tundra / cold steppe
    if (temp < 0.28) {
      if (rain < 0.30) return [150, 155, 130]; // cold dry
      return [120, 145, 120];                  // taiga-ish
    }

    // hot deserts
    if (rain < 0.20) {
      if (temp > 0.55) return [200, 180, 110]; // sand
      return [175, 170, 135];                  // semi-arid / cold desert
    }

    // grassland / savanna
    if (rain < 0.38) {
      if (temp > 0.55) return [165, 175, 80];  // savanna
      return [120, 165, 95];                   // steppe
    }

    // forests
    if (rain < 0.65) {
      if (temp > 0.55) return [55, 140, 80];   // tropical seasonal forest
      return [45, 125, 70];                    // temperate forest
    }

    // rainforests / very wet
    if (temp > 0.55) return [35, 120, 75];     // rainforest
    return [55, 120, 95];                      // wet temperate
  }

  function drawMap(ctx, h, acc, temp, rain, n, seaLevel, riverThreshold, mode) {
    const img = ctx.createImageData(n, n);
    const data = img.data;

    let accMax = 1;
    for (let i = 0; i < acc.length; i++) if (acc[i] > accMax) accMax = acc[i];

    for (let y = 0; y < n; y++) {
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        const height = h[i];

        let rgb;
		if (mode === 'height') rgb = heightColor(height, seaLevel);
		else rgb = biomeColor(height, seaLevel, temp[i], rain[i]);
		let [r, g, b] = rgb;

        // Shade
        const sh = hillshade(h, n, x, y);
        const shade = 0.55 + 0.45 * sh;
        r = (r * shade) | 0;
        g = (g * shade) | 0;
        b = (b * shade) | 0;

        // Rivers
        if (height > seaLevel && acc[i] > riverThreshold) {
          const t = Math.log(acc[i]) / Math.log(accMax + 1e-9);
          const w = clamp((t - 0.25) * 2.0, 0, 1);
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
    const seaLevel = clamp(opts.seaLevel ?? 0.50, 0.05, 0.95);
    const riverThreshold = opts.riverThreshold ?? 1800;

    const continentsMin = opts.continentsMin ?? 2;
    const continentsMax = opts.continentsMax ?? 6;

    const platesCount = clamp((opts.plates ?? 16) | 0, 4, 40);
    const windDir = (opts.windDir ?? 'W'); // 'W','E','N','S'

    canvas.width = n;
    canvas.height = n;
    const ctx = canvas.getContext('2d', { alpha: false });

    // 1) Base height from continents
    let h = buildHeightmap(n, seed, continentsMin, continentsMax);

    // 2) Pseudo-tectonics uplift layer
    const plates = generatePlates(seed, platesCount);
    const plateFields = computePlateFields(n, plates);
    const uplift = buildUpliftMap(n, plates, plateFields, seed);

    // apply uplift (mountain chains / rifts)
	const mountainsStrength = clamp(opts.mountainsStrength ?? 0.18, 0, 0.6);

	for (let i = 0; i < h.length; i++) {
	  // landMask: 0 в океане, 1 на суше, плавный переход у берега
	  const m = clamp((h[i] - seaLevel + 0.06) / 0.18, 0, 1);
	  const landMask = m * m * (3 - 2 * m); // smoothstep

	  // немного оставим океанические хребты, но сильно слабее
	  const oceanFactor = 0.15;
	  const factor = oceanFactor + (1 - oceanFactor) * landMask;

	  h[i] = clamp(h[i] + uplift[i] * mountainsStrength * factor, 0, 1);
	}

    // 3) Erosion (light) after tectonics
    thermalErode(h, n, 7, 0.017);

    // 4) Fill depressions for drainage
    const hf = fillDepressions(h, n);

    // 5) Hydrology
    const down = computeFlowDir(hf, n);
    const acc = computeFlowAccum(hf, down);

    // 6) Climate
    const temp = computeTemperature(h, n, seaLevel);
	const desertStrength = clamp(opts.desertStrength ?? 0.35, 0, 1);
	const rain = computeRain(h, n, seaLevel, windDir, desertStrength);
    // 7) Draw
    const mode = (opts.mode === 'height') ? 'height' : 'biomes';
	drawMap(ctx, h, acc, temp, rain, n, seaLevel, riverThreshold, mode);
  }

  window.TerrainMap = { render };
})();
