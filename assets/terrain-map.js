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
  
  function buildWindField(g, seed) {
  const vx = new Float32Array(g * g);
  const vy = new Float32Array(g * g);
  const spd = new Float32Array(g * g);

  // stream function ψ(u,v)
  const psi = new Float32Array(g * g);

  const freq1 = 3.0;
  const freq2 = 7.0;

  for (let y = 0; y < g; y++) {
    const v = y / (g - 1);
    for (let x = 0; x < g; x++) {
      const u = x / (g - 1);

      // ψ: чем плавнее, тем крупнее вихри (freq ниже = крупнее)
      const p1 = fbm(u * freq1, v * freq1, seed + 2000, 4, 2.0, 0.5);
      const p2 = fbm(u * freq2, v * freq2, seed + 3000, 3, 2.0, 0.5);
      psi[idx(x, y, g)] = 0.75 * p1 + 0.25 * p2;
    }
  }

  // curl from ψ via finite differences
  const eps = 1 / (g - 1);

  for (let y = 0; y < g; y++) {
    const v = y / (g - 1);
    const lat = (v - 0.5) * 2.0;          // -1..1
    const latAbs = Math.abs(lat);         // 0..1

    // Широтные пояса: W->E в умеренных, E->W на экваторе/полюсах
    // -cos(2π*|lat|): -1 (экватор), +1 (средние), -1 (полюса)
    const baseX = -Math.cos(2 * Math.PI * latAbs);

    // Умеренные широты более “вихревые”
    const belt = Math.exp(-Math.pow((latAbs - 0.55) / 0.22, 2)); // пик ~0.55
    const swirlK = 0.55 + 0.85 * belt; // 0.55..1.40

    for (let x = 0; x < g; x++) {
      const xm1 = clamp(x - 1, 0, g - 1);
      const xp1 = clamp(x + 1, 0, g - 1);
      const ym1 = clamp(y - 1, 0, g - 1);
      const yp1 = clamp(y + 1, 0, g - 1);

      const dpsi_dx = (psi[idx(xp1, y, g)] - psi[idx(xm1, y, g)]) / (2 * eps + 1e-9);
      const dpsi_dy = (psi[idx(x, yp1, g)] - psi[idx(x, ym1, g)]) / (2 * eps + 1e-9);

      // curl: v = ( dψ/dy, -dψ/dx )
      let wx = dpsi_dy;
      let wy = -dpsi_dx;

      // Смешиваем с широтной “струёй”
      wx = wx * swirlK + baseX * 1.15;
      wy = wy * swirlK;

      // Нормализация
      const len = Math.sqrt(wx * wx + wy * wy) + 1e-9;
      wx /= len; wy /= len;

      const i = idx(x, y, g);
      vx[i] = wx;
      vy[i] = wy;

      // скорость: чуть больше в умеренных широтах + шум
      const sNoise = fbm((x / (g - 1)) * 2.0, (y / (g - 1)) * 2.0, seed + 4000, 3, 2.0, 0.5);
      spd[i] = clamp(0.55 + 0.55 * belt + (sNoise - 0.5) * 0.25, 0.25, 1.25);
    }
  }

  return { g, vx, vy, spd };
}

function sampleWindBilinear(wind, u, v) {
  const g = wind.g;
  const x = clamp(u * (g - 1), 0, g - 1);
  const y = clamp(v * (g - 1), 0, g - 1);

  const x0 = Math.floor(x), y0 = Math.floor(y);
  const x1 = Math.min(x0 + 1, g - 1);
  const y1 = Math.min(y0 + 1, g - 1);

  const tx = x - x0, ty = y - y0;

  const i00 = idx(x0, y0, g), i10 = idx(x1, y0, g);
  const i01 = idx(x0, y1, g), i11 = idx(x1, y1, g);

  const vx0 = lerp(wind.vx[i00], wind.vx[i10], tx);
  const vx1 = lerp(wind.vx[i01], wind.vx[i11], tx);
  const vy0 = lerp(wind.vy[i00], wind.vy[i10], tx);
  const vy1 = lerp(wind.vy[i01], wind.vy[i11], tx);
  const sp0 = lerp(wind.spd[i00], wind.spd[i10], tx);
  const sp1 = lerp(wind.spd[i01], wind.spd[i11], tx);

  return {
    vx: lerp(vx0, vx1, ty),
    vy: lerp(vy0, vy1, ty),
    spd: lerp(sp0, sp1, ty),
  };
}

function sampleHeightNearest(h, n, u, v) {
  const x = clamp((u * (n - 1) + 0.5) | 0, 0, n - 1);
  const y = clamp((v * (n - 1) + 0.5) | 0, 0, n - 1);
  return h[idx(x, y, n)];
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
function computeRainWind(h, n, seaLevel, wind, desertStrength) {
  const rain = new Float32Array(n * n);

  const ds = clamp(desertStrength ?? 0.35, 0, 1);
  const dsK = ds * ds;

  // параметры: можно крутить
  const steps = 26;                 // сколько “смотрим по ветру”
  const baseStep = 1.8 / wind.g;    // шаг в UV
  const decay = 0.985 - 0.10 * dsK; // чем больше deserts, тем быстрее сохнет
  const oroK = 2.8;                 // орография
  const bg = 0.015;                 // фон (чтобы не было 0 везде)

  for (let y = 0; y < n; y++) {
    const v0 = y / (n - 1);
    for (let x = 0; x < n; x++) {
      const u0 = x / (n - 1);
      const i = idx(x, y, n);
      const h0 = h[i];

      if (h0 <= seaLevel) {
        rain[i] = 1.0;
        continue;
      }

      let u = u0, v = v0;
      let m = 0.0;          // влага
      let accR = 0.0;       // осадки
      let prevH = h0;

      for (let s = 0; s < steps; s++) {
        const w = sampleWindBilinear(wind, u, v);

        // backtrace: идём ПРОТИВ ветра, ищем “откуда пришёл воздух”
        const step = baseStep * w.spd;
        u = clamp(u - w.vx * step, 0, 1);
        v = clamp(v - w.vy * step, 0, 1);

        const hh = sampleHeightNearest(h, n, u, v);

        if (hh <= seaLevel) {
          // над океаном воздух насыщается
          m = 1.0;
          prevH = hh;
          continue;
        }

        // воздух сушеет по мере движения
        m *= decay;

        // осадки на подъёме (орография)
        const up = Math.max(0, hh - prevH);
        const r = clamp(m * (0.06 + up * oroK), 0, 1);

        accR += r;
        m = Math.max(0, m - r * 0.65);

        prevH = hh;
      }

      rain[i] = clamp(accR + m * bg, 0, 1);
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
  
  function biomeColor(height, seaLevel, temp, rain, desertStrength) {
    if (height <= seaLevel) return oceanColor(height, seaLevel);

    const alt = (height - seaLevel) / (1 - seaLevel + 1e-9);
	
	const ds = clamp(desertStrength ?? 0.35, 0, 1);
	const k = ds * ds;

	// Сдвигаем пороги: при больших deserts пустыня появляется чаще.
	const desertCut = 0.18 + 0.18 * k;   // 0.18 .. 0.36
	const steppeCut = 0.36 + 0.12 * k;   // 0.36 .. 0.48

    // snow by altitude or cold
    if (temp < 0.14 || alt > 0.88) return [235, 235, 240];

    // tundra / cold steppe
    if (temp < 0.28) {
      if (rain < 0.30) return [150, 155, 130]; // cold dry
      return [120, 145, 120];                  // taiga-ish
    }

	// hot deserts
	if (rain < desertCut) {
	  if (temp > 0.55) return [200, 180, 110]; // sand
	  return [175, 170, 135];                  // semi-arid / cold desert
	}

	// grassland / savanna
	if (rain < steppeCut) {
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
function drawMap(ctx, h, acc, temp, rain, biomeIds, n, seaLevel, riverThreshold, mode, desertStrength) {
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
      else rgb = biomeIdToColor(biomeIds[i], height, seaLevel);

      let [r, g, b] = rgb;

      const sh = hillshade(h, n, x, y);
      const shade = 0.55 + 0.45 * sh;
      r = (r * shade) | 0;
      g = (g * shade) | 0;
      b = (b * shade) | 0;

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

  
  function debugBiomeStats(h, temp, rain, n, seaLevel, desertStrength) {
  const ds = clamp(desertStrength ?? 0.35, 0, 1);
  const k = ds * ds;
  const desertCut = 0.18 + 0.18 * k;
  const steppeCut = 0.36 + 0.12 * k;

  let land = 0, desert = 0, steppe = 0, forest = 0, snow = 0;
  let rMin =  1e9, rMax = -1e9, rSum = 0;

  for (let i = 0; i < h.length; i++) {
    if (h[i] <= seaLevel) continue;
    land++;

    const alt = (h[i] - seaLevel) / (1 - seaLevel + 1e-9);
    const t = temp[i];
    const r = rain[i];

    rMin = Math.min(rMin, r);
    rMax = Math.max(rMax, r);
    rSum += r;

    if (t < 0.14 || alt > 0.88) { snow++; continue; }
    if (t < 0.28) { /* tundra/taiga */ continue; }

    if (r < desertCut) desert++;
    else if (r < steppeCut) steppe++;
    else if (r < 0.65) forest++;
  }

  const rAvg = land ? (rSum / land) : 0;
  console.log('[BiomeStats] ds=', ds.toFixed(2),
    'desertCut=', desertCut.toFixed(3),
    'steppeCut=', steppeCut.toFixed(3),
    'rain(min/avg/max)=', rMin.toFixed(3), rAvg.toFixed(3), rMax.toFixed(3),
    'land=', land,
    'desert%=', land ? (100*desert/land).toFixed(1) : '0',
    'steppe%=', land ? (100*steppe/land).toFixed(1) : '0',
    'forest%=', land ? (100*forest/land).toFixed(1) : '0',
    'snow%=', land ? (100*snow/land).toFixed(1) : '0'
  );
}

function normalizeRainOverLand(rain, h, n, seaLevel) {
  let min = Infinity, max = -Infinity;

  for (let i = 0; i < rain.length; i++) {
    if (h[i] <= seaLevel) continue;
    const r = rain[i];
    if (r < min) min = r;
    if (r > max) max = r;
  }

  const inv = 1 / (max - min + 1e-9);

  for (let i = 0; i < rain.length; i++) {
    if (h[i] <= seaLevel) {
      rain[i] = 1.0;          // океан всегда “мокрый”
    } else {
      rain[i] = clamp((rain[i] - min) * inv, 0, 1);
    }
  }
}


function drawWindDebug(ctx, wind, n) {
  ctx.clearRect(0, 0, n, n);

  // лёгкая сетка линий
  ctx.fillStyle = '#0b1020';
  ctx.fillRect(0, 0, n, n);

  const stepPx = Math.max(14, (n / 48) | 0);
  ctx.strokeStyle = 'rgba(120, 255, 160, 0.65)';
  ctx.lineWidth = 1;

  for (let y = stepPx; y < n; y += stepPx) {
    for (let x = stepPx; x < n; x += stepPx) {
      const u = x / (n - 1);
      const v = y / (n - 1);
      const w = sampleWindBilinear(wind, u, v);

      const len = 10 + 18 * (w.spd - 0.25) / (1.25 - 0.25);
      const x2 = x + w.vx * len;
      const y2 = y + w.vy * len;

      ctx.beginPath();
      ctx.moveTo(x, y);
      ctx.lineTo(x2, y2);
      ctx.stroke();
    }
  }
}

function biomeId(height, seaLevel, temp, rain, desertStrength) {
  if (height <= seaLevel) return 0; // ocean

  const alt = (height - seaLevel) / (1 - seaLevel + 1e-9);
  const ds = clamp(desertStrength ?? 0.35, 0, 1);
  const k = ds * ds;

  const desertCut = 0.18 + 0.18 * k;
  const steppeCut = 0.36 + 0.12 * k;

  if (temp < 0.14 || alt > 0.88) return 6; // snow
  if (temp < 0.28) return (rain < 0.30) ? 5 : 4; // cold desert / taiga

  if (rain < desertCut) return (temp > 0.55) ? 1 : 2; // hot desert / cold semi-arid
  if (rain < steppeCut) return (temp > 0.55) ? 3 : 7; // savanna / steppe

  if (rain < 0.65) return (temp > 0.55) ? 8 : 9; // seasonal / temperate forest
  return (temp > 0.55) ? 10 : 11; // rainforest / wet temperate
}

function biomeIdToColor(id, height, seaLevel) {
  if (id === 0) return oceanColor(height, seaLevel);
  // подставь свои цвета (ниже примерные)
  const pal = {
    1:[200,180,110], 2:[175,170,135], 3:[165,175, 80], 7:[120,165, 95],
    4:[120,145,120], 5:[150,155,130], 6:[235,235,240],
    8:[ 55,140, 80], 9:[ 45,125, 70],10:[ 35,120, 75],11:[ 55,120, 95],
  };
  return pal[id] ?? [255, 0, 255];
}

function smoothBiomeIds(ids, h, n, seaLevel, iterations) {
  const out = new Int16Array(ids.length);
  const counts = new Int16Array(32); // хватит на небольшую палитру

  for (let it = 0; it < iterations; it++) {
    out.set(ids);

    for (let y = 1; y < n - 1; y++) {
      for (let x = 1; x < n - 1; x++) {
        const i = idx(x, y, n);
        if (h[i] <= seaLevel) continue; // океан не трогаем

        counts.fill(0);

        // текущий + 8 соседей
        for (let k = 0; k < 8; k++) {
          const j = idx(x + N8[k][0], y + N8[k][1], n);
          const id = ids[j];
          if (id > 0) counts[id]++;
        }
        const self = ids[i];
        if (self > 0) counts[self]++;

        // найти победителя
        let bestId = self;
        let best = -1;
        for (let id = 1; id < counts.length; id++) {
          if (counts[id] > best) { best = counts[id]; bestId = id; }
        }

        // меняем только если реально “большинство”, чтобы не смазывать границы
        if (best >= 6) out[i] = bestId; // 6 из 9
      }
    }

    ids.set(out);
  }
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
	const windGrid = buildWindField(Math.max(96, Math.min(160, (n / 6) | 0)), seed);

    const rain = computeRainWind(h, n, seaLevel, windGrid, desertStrength);
	normalizeRainOverLand(rain, h, n, seaLevel);
	
	debugBiomeStats(h, temp, rain, n, seaLevel, desertStrength);

    // 7) Draw	
	const mode = (opts.mode === 'height') ? 'height' : (opts.mode === 'wind' ? 'wind' : 'biomes');

	if (mode === 'wind') {
	  drawWindDebug(ctx, windGrid, n);
	  return;
	}
	
	const biomeIds = new Int16Array(n * n);
	for (let i = 0; i < biomeIds.length; i++) {
	  biomeIds[i] = biomeId(h[i], seaLevel, temp[i], rain[i], desertStrength);
	}
	smoothBiomeIds(biomeIds, h, n, seaLevel, 3);

drawMap(ctx, h, acc, temp, rain, biomeIds, n, seaLevel, riverThreshold, mode, desertStrength);
	
  }

  window.TerrainMap = { render };
})();
