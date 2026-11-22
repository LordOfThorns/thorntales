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

  // The world wraps around east-west like a real planet: the left edge of the map
  // continues at the right edge. Top and bottom are the poles and do not wrap.
  // So every x coordinate is periodic: u = x / n is in [0, 1) and u = 1 is u = 0 again,
  // noise repeats across the width, and neighbours of the first column are in the last.
  function wrapX(x, n) { return ((x % n) + n) % n; }
  // shortest signed east-west distance between two u in [0, 1)
  function wrapD(du) { return du - Math.round(du); }

  // Value noise (grid) with smooth interpolation.
  // px > 0: the lattice repeats every px cells along x (seamless east-west)
  function valueNoise2(x, y, seed, px) {
    let x0 = Math.floor(x);
    const y0 = Math.floor(y);
    let x1 = x0 + 1;
    const y1 = y0 + 1;
    const sx = smoothstep(x - x0);
    const sy = smoothstep(y - y0);
    if (px) { x0 = wrapX(x0, px); x1 = (x0 + 1) % px; }

    const v00 = hash2i(x0, y0, seed);
    const v10 = hash2i(x1, y0, seed);
    const v01 = hash2i(x0, y1, seed);
    const v11 = hash2i(x1, y1, seed);

    const ix0 = lerp(v00, v10, sx);
    const ix1 = lerp(v01, v11, sx);
    return lerp(ix0, ix1, sy);
  }

  // Fractal Brownian Motion (fBm). px > 0: seamless east-west, the lattice of the first
  // octave repeats every px cells (lacunarity must be a whole number for this)
  function fbm(x, y, seed, octaves, lacunarity, gain, px) {
    let amp = 0.5;
    let freq = 1.0;
    let sum = 0.0;
    for (let i = 0; i < octaves; i++) {
      sum += amp * valueNoise2(x * freq, y * freq, seed + i * 1013, px ? px * freq : 0);
      freq *= lacunarity;
      amp *= gain;
    }
    return sum; // ~[0,1]
  }

  // Seamless fBm over the map: u in [0, 1) wraps, s = features across the map (whole number)
  function wfbm(u, v, s, seed, octaves, gain = 0.5) {
    return fbm(u * s * ASPECT, v * s, seed, octaves, 2.0, gain, s * ASPECT);
  }
  
  // -------------------- Wind --------------------
  // Planetary wind belts plus a weak, large-scale meander.
  // lat 0 = equator (map middle), 1 = pole (top/bottom), same as computeTemperature.
  //   easterly trades 0..1/3 (30 deg), westerlies 1/3..2/3, polar easterlies 2/3..1 (as on Earth);
  //   meridional part: trades blow towards the equator, westerlies towards the pole.
  // The meander is the curl of a smooth noise field, scaled to at most WIND_MEANDER of the
  // belt wind, so it bends the flow but cannot turn it into closed eddies.
  // Speed is NOT renormalised: calm belts (~30 and 60 deg) stay calm and become dry,
  // like the subtropical highs on Earth.
  const WIND_MEANDER = 0.35;
  const WIND_MERIDIONAL = 0.35;

  // gh = grid height; the grid is ASPECT times wider, cells stay square
  function buildWindField(gh, seed) {
    const g = Math.round(gh * ASPECT); // grid width (the row length)
    const vx = new Float32Array(g * gh);
    const vy = new Float32Array(g * gh);
    const spd = new Float32Array(g * gh);

    // Low frequency, 2 octaves: broad meanders only (high octaves made small eddies)
    const psi = new Float32Array(g * gh);
    for (let y = 0; y < gh; y++) {
      for (let x = 0; x < g; x++) {
        psi[idx(x, y, g)] = wfbm(x / g, y / (gh - 1), 2, seed + 2000, 2, 0.45);
      }
    }

    // curl(psi) is divergence free; scale it by its maximum
    const cx = new Float32Array(g * gh);
    const cy = new Float32Array(g * gh);
    let cmax = 1e-9;
    for (let y = 0; y < gh; y++) {
      for (let x = 0; x < g; x++) {
        const xm = wrapX(x - 1, g), xp = wrapX(x + 1, g);
        const ym = Math.max(0, y - 1), yp = Math.min(gh - 1, y + 1);
        const dx = psi[idx(xp, y, g)] - psi[idx(xm, y, g)];
        const dy = psi[idx(x, yp, g)] - psi[idx(x, ym, g)];
        const i = idx(x, y, g);
        cx[i] = dy;
        cy[i] = -dx;
        cmax = Math.max(cmax, Math.hypot(dx, dy));
      }
    }

    for (let y = 0; y < gh; y++) {
      const v = y / (gh - 1);
      const north = v < 0.5;
      const lat = Math.abs(v - 0.5) * 2;

      // c cells per hemisphere (3 on Earth): the belt wind turns at the cell borders,
      // easterly in the tropical cell; one slow cell has weak easterlies everywhere
      const c = PLANET.cells;
      const pts = [0, ...cellBorders(c), 1];
      let j = 0;
      while (j < pts.length - 2 && lat > pts[j + 1]) j++;
      // belt wind: strongest inside a cell, calm at its borders (-cos over a stretched
      // latitude that puts the borders at half steps)
      const fl = j === 0 ? (c === 1 ? lat : 0.5 * lat / pts[1])
        : (j - 0.5) + (lat - pts[j]) / (pts[j + 1] - pts[j]) * (j === pts.length - 2 ? 0.5 : 1);
      const zonal = c === 1 ? -0.5 : -Math.cos(Math.PI * fl);
      // surface flow: towards the equator in the tropical cell, then alternating
      const poleward = (j % 2 === 0 ? -1 : 1) * Math.sin(Math.PI * (lat - pts[j]) / (pts[j + 1] - pts[j]));
      const bx = zonal;
      const by = WIND_MERIDIONAL * poleward * (north ? -1 : 1); // screen y points south

      for (let x = 0; x < g; x++) {
        const i = idx(x, y, g);
        let ax = bx, ay = by;
        if (PLANET.locked) {
          // surface wind from the night side towards the point under the star (map centre):
          // the way the sunlight grows; calm under the star and on the far side
          const lon = (x / g - 0.5) * 2 * Math.PI, la = (0.5 - v) * Math.PI;
          ax = -Math.sin(lon);
          ay = Math.sin(la) * Math.cos(lon); // screen y points south
        }
        const wx = ax + WIND_MEANDER * cx[i] / cmax;
        const wy = ay + WIND_MEANDER * cy[i] / cmax;
        const len = Math.hypot(wx, wy) + 1e-9;
        vx[i] = wx / len;
        vy[i] = wy / len;
        spd[i] = clamp(len, 0.05, 1.4);
      }
    }

    return { g, gh, vx, vy, spd };
  }

// u wraps east-west (u = x / width), v is clamped (v = y / (height - 1))
function sampleWindBilinear(wind, u, v) {
  const g = wind.g;
  const x = (u - Math.floor(u)) * g;
  const gh = wind.gh || g;
  const y = clamp(v * (gh - 1), 0, gh - 1);

  const x0 = Math.floor(x) % g, y0 = Math.floor(y);
  const x1 = (x0 + 1) % g;
  const y1 = Math.min(y0 + 1, gh - 1);

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
  const nh = h.length / n;
  const x = wrapX(Math.floor(u * n + 0.5), n);
  const y = clamp((v * (nh - 1) + 0.5) | 0, 0, nh - 1);
  return h[idx(x, y, n)];
}


  // 2D index helpers. Maps are W x H (the world is twice as wide as it is high, like a
  // map of a planet): in every function n is the WIDTH (the row length), the height is
  // passed as nh or taken from an array (nh = array.length / n).
  function idx(x, y, n) { return y * n + x; }

  // Width / height of the world. Noise and distances along x are scaled by it, so
  // features stay round although u (0..1 across) covers twice the distance of v.
  let ASPECT = 1;

  // Neighbors (8-connected)
  const N8 = [
    [-1, -1], [0, -1], [1, -1],
    [-1, 0],           [1, 0],
    [-1, 1],  [0, 1],  [1, 1]
  ];

  // -------------------- Continents / base height --------------------

  function buildHeightmap(n, nh, seed, continentsMin, continentsMax) {
    const h = new Float32Array(n * nh);
    const rng = mulberry32(seed);

    const cmin = Math.max(1, continentsMin | 0);
    const cmax = Math.max(cmin, continentsMax | 0);
    const continents = cmin + ((rng() * (cmax - cmin + 1)) | 0);

    // Each continent is a chain of rotated, elongated blobs along a bending spine:
    // bulky in the middle, tapering to the ends (peninsulas). Sometimes a side branch
    // (a large peninsula) and sometimes a negative blob that bites out a bay or inland sea.
    // Positions are in u (0..1 across the width) and v (0..1 down the height); sizes and
    // steps are in v units, so moves along x are divided by ASPECT
    const blobs = [];
    for (let c = 0; c < continents; c++) {
      // the world is ASPECT times larger than a square one: continents grow by sqrt(ASPECT)
      // so they cover the same share of it
      const size = (0.08 + rng() * 0.10) * Math.sqrt(ASPECT);
      const w = 0.8 + rng() * 1.0;
      const nodes = 3 + ((rng() * 4) | 0);
      const step = size * 0.85;
      let dir = rng() * Math.PI * 2;
      let x = rng() - Math.cos(dir) * step * (nodes - 1) / 2 / ASPECT; // anywhere east-west: the world wraps
      let y = 0.18 + rng() * 0.64 - Math.sin(dir) * step * (nodes - 1) / 2;
      const spine = [];

      for (let j = 0; j < nodes; j++) {
        const t = nodes > 1 ? j / (nodes - 1) : 0.5;
        const taper = (0.5 + 0.5 * Math.sin(Math.PI * t)) * (0.85 + 0.3 * rng());
        blobs.push({
          x, y, ang: dir,
          rx: size * taper * (1.2 + rng() * 0.8),   // along the spine
          ry: size * taper * (0.7 + rng() * 0.4),   // across
          w: w * 0.62,
        });
        spine.push({ x, y, dir, taper });
        dir += (rng() - 0.5) * 1.1;
        x += Math.cos(dir) * step / ASPECT;
        y += Math.sin(dir) * step;
      }

      // side branch: a large peninsula
      if (rng() < 0.6) {
        const s = spine[(rng() * spine.length) | 0];
        let bdir = s.dir + (rng() < 0.5 ? -1 : 1) * (0.8 + rng() * 0.8);
        let bx = s.x, by = s.y;
        const parts = 2 + ((rng() * 2) | 0);
        for (let j = 1; j <= parts; j++) {
          const k = 1 - j / (parts + 1);
          bx += Math.cos(bdir) * step * 0.7 / ASPECT;
          by += Math.sin(bdir) * step * 0.7;
          blobs.push({ x: bx, y: by, ang: bdir, rx: size * 0.55 * k + 0.01, ry: size * 0.32 * k + 0.008, w: w * 0.62 });
          bdir += (rng() - 0.5) * 0.8;
        }
      }

      // bay or inland sea
      if (rng() < 0.45) {
        const s = spine[(rng() * spine.length) | 0];
        const side = rng() < 0.5 ? -1 : 1;
        const off = size * (0.3 + rng() * 0.5);
        blobs.push({
          x: s.x + Math.cos(s.dir + side * Math.PI / 2) * off / ASPECT,
          y: s.y + Math.sin(s.dir + side * Math.PI / 2) * off,
          ang: s.dir, rx: size * (0.35 + rng() * 0.3), ry: size * (0.25 + rng() * 0.2), w: -w * 0.55,
        });
      }
    }
    for (const b of blobs) { b.cos = Math.cos(b.ang); b.sin = Math.sin(b.ang); }

    // Two scales of bending: the broad one moves whole continents, the middle one bends
    // their long sides (a stretched blob alone gave coasts straight as a ruler)
    function warp(u, v) {
      const wx = wfbm(u, v, 2, seed + 900, 3, 0.5) - 0.5;
      const wy = wfbm(u, v, 2, seed + 901, 3, 0.5) - 0.5;
      const mx = wfbm(u, v, 6, seed + 902, 2, 0.5) - 0.5;
      const my = wfbm(u, v, 6, seed + 903, 2, 0.5) - 0.5;
      return [u + (wx * 0.10 + mx * 0.07) / ASPECT, v + wy * 0.10 + my * 0.07];
    }

    for (let y = 0; y < nh; y++) {
      for (let x = 0; x < n; x++) {
        const u0 = x / n; // periodic east-west
        const v0 = y / (nh - 1);
        const [u, v] = warp(u0, v0);

        let cont = 0;
        for (let k = 0; k < blobs.length; k++) {
          const b = blobs[k];
          const dx = wrapD(u - b.x) * ASPECT, dy = v - b.y; // through the east-west seam, in v units
          const along = (dx * b.cos + dy * b.sin) / b.rx;
          const across = (-dx * b.sin + dy * b.cos) / b.ry;
          const q = along * along + across * across;
          if (q > 12) continue; // negligible, skip exp
          cont += Math.exp(-q) * b.w;
        }

        cont = cont / (continents * 0.95);
        cont = clamp(cont, 0, 1);
        cont = Math.pow(cont, 1.35);

        const shore = wfbm(u, v, 6, seed + 111, 5, 0.5);
        const shore2 = wfbm(u, v, 18, seed + 112, 4, 0.5);
        // fine detail: small bays, capes and inlets along every coast
        const shore3 = wfbm(u, v, 46, seed + 113, 2, 0.5);
        const coast = clamp(cont + (shore - 0.5) * 0.22 + (shore2 - 0.5) * 0.10 + (shore3 - 0.5) * 0.06, 0, 1);

        const d = wfbm(u, v, 22, seed + 30, 5, 0.5);
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

      // oceanic vs continental tendency (decides which side of a collision gets the range)
      const type = rng() < 0.55 ? 0 : 1; // 0 oceanic, 1 continental
      // Size radius with a heavy tail: most plates are small, a few are giant
      // (like Earth). Bigger r pushes the plate's edges outwards.
      const r = PLATE_R_MAX * Math.pow(rng(), 3);
      p.push({ x, y, vx, vy, type, r2: r * r });
    }
    return p;
  }

  const PLATE_R_MAX = 0.32;

  // Power diagram (weighted Voronoi with score = distance^2 - r^2): plates get very
  // different sizes, cells stay convex (no ring-shaped microplates) and the edge
  // between two plates is still a straight line, so the distance to it is exact.
  // For each cell: nearest plate A, second nearest B and the distance to the A|B edge.
  // The lookup position is domain-warped so edges wobble.
  const PLATE_WARP = 0.10;

  function computePlateFields(n, nh, plates, seed) {
    const ids = new Int16Array(n * nh);
    const ids2 = new Int16Array(n * nh);
    const edge = new Float32Array(n * nh);

    for (let y = 0; y < nh; y++) {
      const v0 = y / (nh - 1);
      for (let x = 0; x < n; x++) {
        const u0 = x / n; // periodic east-west
        const u = u0 + (wfbm(u0, v0, 3, seed + 610, 3, 0.5) - 0.5) * 2 * PLATE_WARP / ASPECT;
        const v = v0 + (wfbm(u0, v0, 3, seed + 611, 3, 0.5) - 0.5) * 2 * PLATE_WARP;

        let a = 0, b = 0, da = Infinity, db = Infinity;
        for (let k = 0; k < plates.length; k++) {
          const dx = wrapD(u - plates[k].x) * ASPECT, dy = v - plates[k].y; // in v units
          const s = dx * dx + dy * dy - plates[k].r2;
          if (s < da) { db = da; b = a; da = s; a = k; }
          else if (s < db) { db = s; b = k; }
        }

        const i = idx(x, y, n);
        ids[i] = a;
        ids2[i] = b;
        // distance from p to the A|B power line
        const A = plates[a], B = plates[b];
        const ab = Math.hypot((wrapD(u - A.x) - wrapD(u - B.x)) * ASPECT, B.y - A.y) + 1e-9; // A and B as seen from here
        edge[i] = (db - da) / (2 * ab);
      }
    }

    return { ids, ids2, edge };
  }

  // Ridged multifractal noise: sharp crests with valleys between them, 0..~1.
  // Seamless east-west: x = u * s, the lattice repeats every s cells (octaves double it)
  function ridged(x, y, seed, octaves, s) {
    let sum = 0, amp = 0.5, freq = 1, weight = 1, norm = 0;
    for (let o = 0; o < octaves; o++) {
      let r = 1 - Math.abs(2 * valueNoise2(x * freq, y * freq, seed + o * 131, s * freq) - 1);
      r *= r;
      r *= weight;
      weight = clamp(r * 1.6, 0, 1);
      sum += r * amp;
      norm += amp;
      amp *= 0.5;
      freq *= 2;
    }
    return sum / norm;
  }

  // Mountain/rift field along plate edges, roughly -0.5..1.2 (not normalised:
  // collision strength matters). Profiles depend on the kind of boundary:
  //   convergent, both continental:  wide range centred on the edge
  //   convergent, mixed:             range on the continental side, trench on the other
  //   divergent:                     rift valley (negative), slight shoulders
  //   transform (sliding):           low broken hills
  // Inside a belt, ridged noise makes separate crests, peaks and passes, and the belt
  // strength changes along the edge, so ranges have gaps and high sections.
  function buildUpliftMap(n, plates, plateFields, seed) {
    const nh = plateFields.ids.length / n; // map height (n is the width)
    const { ids, ids2, edge } = plateFields;
    const uplift = new Float32Array(n * nh);

    for (let y = 0; y < nh; y++) {
      const v = y / (nh - 1);
      for (let x = 0; x < n; x++) {
        const u = x / n; // periodic east-west
        const i = idx(x, y, n);
        const d = edge[i];
        if (d > 0.16) continue;

        const A = plates[ids[i]], B = plates[ids2[i]];
        let nx = wrapD(B.x - A.x) * ASPECT, ny = B.y - A.y;
        const inv = 1 / (Math.hypot(nx, ny) + 1e-9);
        nx *= inv; ny *= inv;
        const rvx = B.vx - A.vx, rvy = B.vy - A.vy;
        const conv = -(rvx * nx + rvy * ny);          // > 0 plates collide
        const shear = Math.abs(rvx * -ny + rvy * nx);

        // belt strength and width vary along the edge
        const along = wfbm(u, v, 5, seed + 520, 3, 0.5);       // ~0.2..0.8
        // varies along the edge, but not so much that a range breaks into separate lumps
        const strength = clamp(0.45 + (along - 0.5) * 1.6, 0.4, 1.25);
        const width = 0.035 + 0.03 * wfbm(u, v, 4, seed + 521, 2, 0.5);

        // crest detail; warped so ridges are not grid aligned
        const wu = u + (wfbm(u, v, 6, seed + 530, 2, 0.5) - 0.5) * 0.05 / ASPECT;
        const wv = v + (wfbm(u, v, 6, seed + 531, 2, 0.5) - 0.5) * 0.05;
        const crest = ridged(wu * 22 * ASPECT, wv * 22, seed + 540, 5, 22 * ASPECT);         // ~0..1, crests ~10-30 px apart

        let up = 0;
        if (conv > 0.05) {
          let center = 0;
          if (A.type === 1 && B.type === 0) center = width * 0.9;      // range inland on A
          const t = (d - center) / width;
          const env = Math.exp(-t * t);
          const both = A.type === 1 && B.type === 1 ? 1.35 : 1;        // continental collision
          up += conv * strength * both * env * (0.25 + 1.15 * crest);
          if (A.type === 0 && B.type === 1) {
            up -= conv * 0.3 * Math.exp(-Math.pow(d / (width * 0.9), 2));   // trench
          }
        } else if (conv < -0.05) {
          const t = d / (width * 0.6);
          up += conv * 0.18 * Math.exp(-t * t);                        // rift valley (conv < 0); sometimes a rift lake
          const s = (d - width * 0.9) / (width * 0.5);
          up += -conv * 0.12 * Math.exp(-s * s) * (0.5 + crest);       // rift shoulders
        }
        const tt = d / (width * 0.7);
        up += shear * 0.12 * Math.exp(-tt * tt) * crest;               // transform hills

        uplift[i] = up;
      }
    }

    const out = blurred(uplift, n, Math.max(1, Math.round(n / 384)));

    // Normalise by the top 1% of uplift, so the Mountains slider sets peak height
    // predictably (the raw value is a product of several random factors, ~0.2 typical)
    const pos = [];
    for (let i = 0; i < out.length; i += 3) if (out[i] > 0) pos.push(out[i]);
    pos.sort((a, b) => a - b);
    const p99 = pos.length ? pos[Math.floor(pos.length * 0.99)] : 1;
    const inv = 1 / Math.max(p99, 1e-6);
    // rifts and trenches stay shallow: lowlands are only ~0.07..0.17 above the sea
    for (let i = 0; i < out.length; i++) {
      const v = out[i] * inv;
      out[i] = v > 0 ? Math.min(v, 1.3) : Math.max(v * 0.35, -0.3);
    }
    return out;
  }

  // A plate counts as continental when enough of it is land in the base heightmap,
  // so collisions and coastal ranges happen where the continents actually are.
  function assignPlateTypes(plates, plateFields, h, seaLevel) {
    const land = new Float64Array(plates.length);
    const all = new Float64Array(plates.length);
    for (let i = 0; i < h.length; i += 3) {
      const id = plateFields.ids[i];
      all[id]++;
      if (h[i] > seaLevel) land[id]++;
    }
    for (let k = 0; k < plates.length; k++) {
      plates[k].type = all[k] > 0 && land[k] / all[k] > 0.3 ? 1 : 0;
    }
  }

  // Box blur into a new array: wraps east-west, clamped at the poles
  function blurred(src, n, radius) {
    const nh = src.length / n; // map height (n is the width)
    const r = Math.max(0, radius | 0);
    if (r === 0) return Float32Array.from(src);
    const tmp = new Float32Array(src.length);
    const out = new Float32Array(src.length);
    const inv = 1 / (2 * r + 1);

    for (let y = 0; y < nh; y++) {
      const row = y * n;
      let sum = 0;
      for (let k = -r; k <= r; k++) sum += src[row + wrapX(k, n)];
      for (let x = 0; x < n; x++) {
        tmp[row + x] = sum * inv;
        sum += src[row + wrapX(x + r + 1, n)] - src[row + wrapX(x - r, n)];
      }
    }
    for (let x = 0; x < n; x++) {
      let sum = 0;
      for (let k = -r; k <= r; k++) sum += tmp[clamp(k, 0, nh - 1) * n + x];
      for (let y = 0; y < nh; y++) {
        out[y * n + x] = sum * inv;
        sum += tmp[Math.min(y + r + 1, nh - 1) * n + x] - tmp[Math.max(y - r, 0) * n + x];
      }
    }
    return out;
  }

  // -------------------- Relief on land --------------------
  // The continent field defines the coastline only. Its inland part is compressed into
  // lowlands/plateaus (a shelf instead of a dome), then hills and mountains are added.
  const PLAINS_HEIGHT = 0.07;   // how high inland lowlands get above the sea
  const PLAINS_SCALE = 0.05;    // how fast land rises from the coast to that level
  const PLATEAU_HEIGHT = 0.10;  // broad uplands in some regions
  const HILLS_HEIGHT = 0.07;

  // -------------------- Coast types --------------------
  // Fjords on cold mountainous coasts; smooth coasts with sand bars and lagoons on warm
  // lowland coasts (only in some stretches, chosen by noise). Distances are in px at 768.
  const FJORD_REACH = 34;       // how far inland fjords cut
  const FJORD_DEPTH = 0.32;
  const SMOOTH_REACH = 7;       // band around the coastline that gets smoothed
  const BAR_FROM = 4, BAR_TO = 6; // sand bar distance from the shore
  const BAR_MAX_DEPTH = 0.08;   // only over shallow water (depth as share of sea level)

  // Breadth-first distance (4-neighbour, in cells) from every cell where src is 1
  function bfsDistance(src, n, maxDist) {
    const nh = src.length / n; // map height (n is the width)
    const dist = new Int16Array(n * nh).fill(-1);
    const queue = new Int32Array(n * nh);
    let qh = 0, qt = 0;
    for (let i = 0; i < src.length; i++) if (src[i]) { dist[i] = 0; queue[qt++] = i; }
    while (qh < qt) {
      const i = queue[qh++];
      const d = dist[i];
      if (d >= maxDist) continue;
      const x = i % n, y = (i / n) | 0;
      const l = x > 0 ? i - 1 : i + n - 1;   // east-west neighbours wrap
      const r = x < n - 1 ? i + 1 : i - n + 1;
      if (dist[l] < 0) { dist[l] = d + 1; queue[qt++] = l; }
      if (dist[r] < 0) { dist[r] = d + 1; queue[qt++] = r; }
      if (y > 0 && dist[i - n] < 0) { dist[i - n] = d + 1; queue[qt++] = i - n; }
      if (y < nh - 1 && dist[i + n] < 0) { dist[i + n] = d + 1; queue[qt++] = i + n; }
    }
    return dist;
  }

  function shapeCoasts(h, n, seed, seaLevel) {
    const nh = h.length / n; // map height (n is the width)
    const N = n * nh;
    const px = nh / 768;
    const isSea = new Uint8Array(N), isLand = new Uint8Array(N);
    for (let i = 0; i < N; i++) { if (h[i] <= seaLevel) isSea[i] = 1; else isLand[i] = 1; }
    const reach = Math.ceil(Math.max(FJORD_REACH, SMOOTH_REACH, BAR_TO + 1) * px);
    const toSea = bfsDistance(isSea, n, reach);   // on land: distance to the sea
    const toLand = bfsDistance(isLand, n, reach); // in the sea: distance to land

    // Per cell: how cold (same bent latitude as computeTemperature) and which coast style
    const cold = new Float32Array(N);
    const smoothType = new Float32Array(N);
    for (let y = 0; y < nh; y++) {
      const v = y / (nh - 1);
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        if (toSea[i] < 0 && toLand[i] < 0) continue; // far from any coast
        const u = x / n; // periodic east-west
        const warp = (wfbm(u, v, 3, seed + 800, 3, 0.5) - 0.5) * 2 * TEMP_LAT_WARP;
        const lat = clamp(Math.abs(v - 0.5) * 2.0 + warp, 0, 1);
        cold[i] = smoothstep(clamp((lat - 0.55) / 0.2, 0, 1));
        smoothType[i] = smoothstep(clamp((wfbm(u, v, 3, seed + 1250, 3, 0.5) - 0.48) * 5, 0, 1))
          * (1 - cold[i]);
      }
    }

    // Fjords: narrow branching channels (ridge lines of ridged noise) cut into cold
    // mountainous coasts, deepest at the shore, fading inland
    const fjordReach = FJORD_REACH * px;
    for (let y = 0; y < nh; y++) {
      const v = y / (nh - 1);
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        const d = toSea[i];
        if (d <= 0 || d > fjordReach || cold[i] <= 0) continue;
        const relief = clamp((h[i] - seaLevel - 0.015) / 0.05, 0, 1); // hilly or mountainous coasts
        if (relief <= 0) continue;
        const u = x / n; // periodic east-west
        const lines = ridged(u * 26 * ASPECT, v * 26, seed + 1300, 3, 26 * ASPECT);
        const channel = smoothstep(clamp((lines - 0.6) / 0.2, 0, 1));
        const fade = Math.pow(1 - d / fjordReach, 0.6);
        h[i] -= channel * fade * Math.sqrt(cold[i] * relief) * FJORD_DEPTH;
      }
    }

    // Smooth coasts: blend heights near the shoreline towards a blurred copy
    const soft = blurred(h, n, Math.max(2, Math.round(5 * px)));
    const smoothReach = SMOOTH_REACH * px;
    for (let i = 0; i < N; i++) {
      const s = smoothType[i];
      if (s <= 0) continue;
      const d = toSea[i] >= 0 && isLand[i] ? toSea[i] : (toLand[i] >= 0 ? toLand[i] : 999);
      if (d > smoothReach) continue;
      const k = s * (1 - d / (smoothReach + 1));
      // never raise mountains or deepen open sea much: blend only near sea level
      if (Math.abs(h[i] - seaLevel) > 0.06) continue;
      h[i] = lerp(h[i], soft[i], k);
    }

    // Sand bars over shallow water in front of smooth coasts; gaps keep lagoons open
    const barFrom = BAR_FROM * px, barTo = BAR_TO * px;
    for (let y = 0; y < nh; y++) {
      const v = y / (nh - 1);
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        if (h[i] > seaLevel || smoothType[i] < 0.35) continue;
        const d = toLand[i];
        if (d < barFrom || d > barTo) continue;
        const depth = (seaLevel - h[i]) / seaLevel;
        if (depth > BAR_MAX_DEPTH) continue;
        const u = x / n; // periodic east-west
        if (wfbm(u, v, 40, seed + 1270, 2, 0.5) < 0.38) continue; // inlets
        h[i] = seaLevel + 0.003;
      }
    }
  }

  const MOUNTAINS_SCALE =3.0;   // normalised uplift (top 1% = 1) x Mountains slider (0.18) x this
  const MOUNTAINS_AT_SEA = 0.25;

  function shapeLand(h, n, seed, seaLevel) {
    const nh = h.length / n; // map height (n is the width)
    for (let y = 0; y < nh; y++) {
      const v = y / (nh - 1);
      for (let x = 0; x < n; x++) {
        const u = x / n; // periodic east-west
        const i = idx(x, y, n);
        const e = h[i] - seaLevel;
        if (e <= 0) continue;
        const coastRise = 1 - Math.exp(-e / PLAINS_SCALE);             // 0 at coast -> 1 inland
        const plateau = smoothstep(clamp((wfbm(u, v, 3, seed + 700, 3, 0.5) - 0.45) * 4, 0, 1));
        const hillMask = smoothstep(clamp((wfbm(u, v, 4, seed + 710, 2, 0.5) - 0.4) * 3, 0, 1));
        const hills = ridged(u * 16 * ASPECT, v * 16, seed + 720, 4, 16 * ASPECT);
        h[i] = seaLevel + coastRise * (PLAINS_HEIGHT + PLATEAU_HEIGHT * plateau + HILLS_HEIGHT * hillMask * hills);
      }
    }
  }

  // -------------------- Erosion --------------------

  // Stream-power erosion: each pass routes water over the terrain and lowers every land
  // cell by k * sqrt(flow) * slope towards the cell it drains into. Big streams cut
  // valleys, the land left between them becomes sharp, branching ridges. On flat
  // lowlands the slope is ~0, so they barely change.
  const STREAM_ERODE_PASSES = 4;
  const STREAM_ERODE_K = 0.15;   // upper slopes: flow ~1..100 cells, slope ~0.005

  // Erosion is computed at half resolution (4x cheaper) and the removed material is
  // applied to the full map with bilinear interpolation, which also softens D8 steps.
  function streamErodeHalf(h, n, seaLevel, passes, k) {
    const nh = h.length / n; // map height (n is the width)
    if (nh < 384) { streamErode(h, n, seaLevel, passes, k); return; }
    const m = n >> 1, mh = nh >> 1; // half width, half height
    const small = new Float32Array(m * mh);
    for (let y = 0; y < mh; y++) {
      for (let x = 0; x < m; x++) {
        const i = idx(x * 2, y * 2, n);
        small[idx(x, y, m)] = (h[i] + h[i + 1] + h[i + n] + h[i + n + 1]) * 0.25;
      }
    }
    const before = Float32Array.from(small);
    // k is per cell; at half resolution cells are twice as long, flow counts 4x smaller
    streamErode(small, m, seaLevel, passes, k * 1.4);
    for (let y = 0; y < nh; y++) {
      const v = clamp((y + 0.5) / 2 - 0.5, 0, mh - 1);
      const y0 = Math.floor(v), y1 = Math.min(y0 + 1, mh - 1), ty = v - y0;
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        if (h[i] <= seaLevel) continue;
        const u = (x + 0.5) / 2 - 0.5;           // wraps east-west
        const fx = Math.floor(u), tx = u - fx;
        const x0 = wrapX(fx, m), x1 = (x0 + 1) % m;
        const d = (a, b) => small[idx(a, b, m)] - before[idx(a, b, m)];
        const delta = lerp(lerp(d(x0, y0), d(x1, y0), tx), lerp(d(x0, y1), d(x1, y1), tx), ty);
        h[i] = Math.max(seaLevel + 1e-4, h[i] + delta);
      }
    }
  }

  function streamErode(h, n, seaLevel, passes, k) {
    for (let p = 0; p < passes; p++) {
      const { filled, order } = fillDepressions(h, n, seaLevel);
      const down = computeFlowDir(filled, n);
      const acc = computeFlowAccum(down, order);
      for (let i = 0; i < h.length; i++) {
        const j = down[i];
        if (j < 0 || h[i] <= seaLevel) continue;
        const drop = h[i] - h[j];
        if (drop <= 0) continue;
        const dist = (i % n !== j % n) && (((i / n) | 0) !== ((j / n) | 0)) ? 1.414 : 1;
        const e = k * Math.sqrt(acc[i]) * (drop / dist);
        h[i] -= Math.min(e, drop * 0.9);
      }
    }
  }

  // Index of the neighbour (x + dx, y + dy): wraps east-west, -1 beyond a pole
  function nbr(x, y, dx, dy, n, nh) {
    const ny = y + dy;
    if (ny < 0 || ny >= nh) return -1;
    return ny * n + wrapX(x + dx, n);
  }

  function thermalErode(h, n, iterations, talus) {
    const nh = h.length / n; // map height (n is the width)
    const tmp = new Float32Array(h.length);

    for (let it = 0; it < iterations; it++) {
      tmp.set(h);

      for (let y = 1; y < nh - 1; y++) {
        for (let x = 0; x < n; x++) {          // every column: east-west wraps
          const i = idx(x, y, n);
          const hi = tmp[i];

          let maxDrop = 0;
          let to = -1;

          for (let k = 0; k < 8; k++) {
            const j = nbr(x, y, N8[k][0], N8[k][1], n, nh);
            const drop = hi - tmp[j];
            if (drop > maxDrop) {
              maxDrop = drop;
              to = j;
            }
          }

          if (maxDrop > talus) {
            const move = (maxDrop - talus) * 0.25;
            h[i] -= move;
            h[to] += move;
          }
        }
      }

      for (let i = 0; i < h.length; i++) h[i] = clamp(h[i], 0, 1);
    }
  }

  // -------------------- Hydrology --------------------

  // Priority-Flood + epsilon (Barnes et al. 2014): depressions are filled up to their
  // spill level, and every filled cell gets a tiny extra rise away from the outlet, so
  // water keeps flowing across filled areas instead of stopping on a perfect flat.
  // Returns the filled surface (Float64 so the tiny steps are representable) and
  // the order cells were processed in (non-decreasing height), used for accumulation.
  const FILL_EPSILON = 1e-7;

  // Outlets (where water leaves the map): the polar rows and, when seaLevel is given,
  // every cell of the sea. East and west are not edges any more: the world wraps.
  function fillDepressions(h, n, seaLevel) {
    const nh = h.length / n; // map height (n is the width)
    const size = n * nh;
    const visited = new Uint8Array(size);
    const out = Float64Array.from(h);
    const order = new Int32Array(size);
    let orderLen = 0;

    const heapH = new Float64Array(size);
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
      markAndPush(idx(x, nh - 1, n));
    }
    if (seaLevel !== undefined) {
      for (let i = n; i < size - n; i++) if (h[i] <= seaLevel) markAndPush(i);
    }

    while (heapSize > 0) {
      const [ch, ci] = heapPop();
      order[orderLen++] = ci;
      const cx = ci % n;
      const cy = (ci / n) | 0;

      for (let k = 0; k < 8; k++) {
        const ni = nbr(cx, cy, N8[k][0], N8[k][1], n, nh);
        if (ni < 0 || visited[ni]) continue;
        visited[ni] = 1;

        if (out[ni] <= ch) out[ni] = ch + FILL_EPSILON;
        heapPush(out[ni], ni);
      }
    }

    return { filled: out, order };
  }

  function computeFlowDir(h, n) {
    const nh = h.length / n; // map height (n is the width)
    const down = new Int32Array(n * nh);
    down.fill(-1);

    for (let y = 1; y < nh - 1; y++) {
      for (let x = 0; x < n; x++) {          // every column: east-west wraps
        const i = idx(x, y, n);
        const hi = h[i];

        let best = -1;
        let bestH = hi;

        for (let k = 0; k < 8; k++) {
          const j = nbr(x, y, N8[k][0], N8[k][1], n, nh);
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

  // Water a land cell adds to rivers: a little everywhere, much more where it rains.
  // Average land rain is ~0.4, so a cell gives ~0.65 on average and the Rivers
  // threshold keeps roughly its old meaning (drainage area in cells).
  const RUNOFF_BASE = 0.15;
  const RUNOFF_RAIN = 1.25;

  // order: cells by non-decreasing filled height (from fillDepressions), so walking it
  // backwards visits every cell before the lower cell it drains into.
  // runoff (optional): water each cell contributes; 1 per cell if omitted.
  function computeFlowAccum(down, order, runoff) {
    const nCells = down.length;
    const acc = new Float32Array(nCells);
    if (runoff) acc.set(runoff); else acc.fill(1);

    for (let k = order.length - 1; k >= 0; k--) {
      const i = order[k];
      const d = down[i];
      if (d >= 0) acc[d] += acc[i];
    }

    return acc;
  }

  // -------------------- Planet (prototype) --------------------
  // opts.planet = { tilt (degrees), flux (sunlight, Earth = 1), cells (wind cells per
  // hemisphere: 1 slow rotation, 3 like Earth, 5 fast) }. Everything is measured against
  // Earth: with tilt 23.44, flux 1 and 3 cells the maps are exactly as without it.
  //  - Temperature: annual mean sunlight by latitude for the tilt (the standard daily
  //    insolation formula averaged over a year), times flux, as an effective temperature
  //    (sunlight^1/4), turned into the map's 0..1 temperature through Earth's own curve.
  //    At tilt > ~54 degrees the poles get more sun in a year than the equator.
  //  - Slow rotation carries heat better: the equator-pole contrast is softened.
  //  - Wind cells: the belts (and the calm, dry belts between them) follow the cell count;
  //    rain is also raised where surface winds meet and lowered where they part, relative
  //    to Earth's pattern.
  //  - Atmosphere (opts.planet.atmo): greenhouse warming, how well heat is spread between
  //    the equator and the poles (and between day and night), how strong the seasons are.
  //  - Tidal locking (opts.planet.locked): one side always faces the star. Zones then lie
  //    in rings around the point under the star (the middle of the map): a hot day side,
  //    a frozen night side, a liveable ring near the edge of the light. Winds blow from the
  //    night side towards that point, where the air rises and it rains most.
  //  - Seasons (seasonLayer below): snow and sea ice for one moment of the year, dry seasons
  //    of savannas and shrubs; drawn on top of the map, so it can change without a new map.
  const EARTH = { tilt: 23.44, flux: 1, cells: 3, locked: false, atmo: 'earth', ecc: 0.0167 };
  const ATMO = {
    //        greenhouse  heat spread  season strength
    thin:  { green: -0.15, spread: 1.2,  season: 1.3 },
    earth: { green: 0,     spread: 1,    season: 1 },
    thick: { green: 0.14,  spread: 0.7,  season: 0.6 },
    venus: { green: 0.4,   spread: 0.35, season: 0.25 },
  };
  // Earth curve of the 0..1 temperature by latitude (0 equator, 1 pole). With the zone
  // thresholds of biomeId (and cooling with height) it puts hot zones up to ~35 deg, temperate
  // forests and steppes around 45 deg, taiga from ~55 deg, tundra from ~65 deg, as on Earth
  // (1.2 put rainforest at 40 deg, 0.85 taiga and tundra at 45 deg)
  const EARTH_TEMP_EXP = 1.0;
  let PLANET = EARTH;
  let LAT_TEMP = null;      // temperature at lat 0..1 (LAT_STEPS + 1 values), null = Earth
  const LAT_STEPS = 256;
  const BELT_RAIN = 0.7;    // how strongly rising / sinking air changes rain
  const HEAT_SPREAD = 0.8;  // share of the difference from Earth's curve left after heat transport

  // Annual mean daily sunlight at latitude phi for axial tilt eps (radians); 1/pi of the
  // solar constant units, circular orbit
  function annualInsolation(phi, eps) {
    const N = 120;
    const p = Math.min(phi, Math.PI / 2 - 1e-4);
    let s = 0;
    for (let k = 0; k < N; k++) {
      const dec = Math.asin(Math.sin(eps) * Math.sin(((k + 0.5) / N) * 2 * Math.PI));
      const c = -Math.tan(p) * Math.tan(dec);
      const h0 = c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c); // half day length
      s += h0 * Math.sin(p) * Math.sin(dec) + Math.cos(p) * Math.cos(dec) * Math.sin(h0);
    }
    return s / N / Math.PI;
  }

  function buildLatTemp(planet) {
    const earthTemp = (lat) => 1 - Math.pow(lat, EARTH_TEMP_EXP); // the generator's Earth curve
    const tauOf = (tilt, flux) => {
      const a = new Float64Array(LAT_STEPS + 1);
      for (let k = 0; k <= LAT_STEPS; k++) {
        a[k] = Math.pow(flux * annualInsolation((k / LAT_STEPS) * Math.PI / 2, tilt * Math.PI / 180), 0.25);
      }
      return a;
    };
    const te = tauOf(EARTH.tilt, 1);
    const norm = te[0];
    for (let k = 0; k <= LAT_STEPS; k++) te[k] /= norm; // Earth: 1 at the equator, falling
    const tp = tauOf(planet.tilt, planet.flux);
    const slope = 1 / (te[0] - te[LAT_STEPS]); // temperature per unit of tau beyond Earth's range
    const out = new Float32Array(LAT_STEPS + 1);
    for (let k = 0; k <= LAT_STEPS; k++) {
      const tau = tp[k] / norm;
      let temp;
      if (tau >= te[0]) temp = 1 + (tau - te[0]) * slope;
      else if (tau <= te[LAT_STEPS]) temp = (tau - te[LAT_STEPS]) * slope;
      else {
        // Earth's latitude with the same sunlight (te falls monotonically)
        let lo = 0, hi = LAT_STEPS;
        while (hi - lo > 1) { const m = (lo + hi) >> 1; if (te[m] >= tau) lo = m; else hi = m; }
        const f = (te[lo] - tau) / (te[lo] - te[hi] + 1e-12);
        temp = earthTemp((lo + f) / LAT_STEPS);
      }
      out[k] = clamp(temp, -0.8, 1.6);
    }
    // Winds and currents carry heat: Earth's curve already includes that, but the
    // difference from Earth is spread out as well (pulled towards its global mean)
    {
      let mean = 0, wsum = 0;
      for (let k = 0; k <= LAT_STEPS; k++) {
        const w = Math.cos((k / LAT_STEPS) * Math.PI / 2); // area of the latitude band
        mean += (out[k] - earthTemp(k / LAT_STEPS)) * w; wsum += w;
      }
      mean /= wsum;
      for (let k = 0; k <= LAT_STEPS; k++) {
        const e = earthTemp(k / LAT_STEPS);
        out[k] = e + mean + HEAT_SPREAD * (out[k] - e - mean);
      }
    }
    // slow rotation spreads heat towards the poles, fast rotation keeps it in the tropics;
    // a thick atmosphere spreads it too, a thin one hardly
    const atmo = ATMO[planet.atmo] || ATMO.earth;
    const spread = (planet.cells <= 1 ? 0.7 : planet.cells >= 5 ? 1.15 : 1) * atmo.spread;
    spreadAround(out, spread, (k) => Math.cos((k / LAT_STEPS) * Math.PI / 2));
    for (let k = 0; k <= LAT_STEPS; k++) out[k] += atmo.green;
    return out;
  }

  // Pull a table towards its area-weighted mean (heat transport); w(k) = area of step k
  function spreadAround(out, spread, w) {
    if (spread === 1) return;
    let mean = 0, wsum = 0;
    for (let k = 0; k <= LAT_STEPS; k++) { mean += out[k] * w(k); wsum += w(k); }
    mean /= wsum;
    for (let k = 0; k <= LAT_STEPS; k++) out[k] = mean + (out[k] - mean) * spread;
  }

  // Effective temperature (sunlight^1/4, Earth's equator = 1) turned into the map's 0..1
  // temperature through Earth's curve, extended in a straight line beyond Earth's range
  function tauToTemp(tau, te) {
    const slope = 1 / (te[0] - te[LAT_STEPS]);
    if (tau >= te[0]) return 1 + (tau - te[0]) * slope;
    if (tau <= te[LAT_STEPS]) return (tau - te[LAT_STEPS]) * slope;
    let lo = 0, hi = LAT_STEPS;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (te[m] >= tau) lo = m; else hi = m; }
    const f = (te[lo] - tau) / (te[lo] - te[hi] + 1e-12);
    return 1 - Math.pow((lo + f) / LAT_STEPS, EARTH_TEMP_EXP);
  }

  // Earth's annual sunlight by latitude as tau (equator = 1) and its scale
  let EARTH_TAU = null;
  function earthTau() {
    if (EARTH_TAU) return EARTH_TAU;
    const te = new Float64Array(LAT_STEPS + 1);
    for (let k = 0; k <= LAT_STEPS; k++) {
      te[k] = Math.pow(annualInsolation((k / LAT_STEPS) * Math.PI / 2, EARTH.tilt * Math.PI / 180), 0.25);
    }
    const norm = te[0];
    for (let k = 0; k <= LAT_STEPS; k++) te[k] /= norm;
    EARTH_TAU = { te, norm };
    return EARTH_TAU;
  }

  // Tidally locked planet: temperature by the angle from the point under the star,
  // s = angle / 180 deg (0 under the star, 0.5 the edge of the light, 1 the far side)
  function buildLockedTemp(planet) {
    const { te, norm } = earthTau();
    const atmo = ATMO[planet.atmo] || ATMO.earth;
    const out = new Float32Array(LAT_STEPS + 1);
    for (let k = 0; k <= LAT_STEPS; k++) {
      const sun = Math.max(0, Math.cos((k / LAT_STEPS) * Math.PI)); // the star always at this height
      const tau = Math.pow(planet.flux * sun, 0.25) / norm;
      out[k] = clamp(tauToTemp(tau, te), -0.8, 1.6);
    }
    // wind carries heat from the day side to the night side; much more with thick air
    const rot = planet.cells <= 1 ? 0.85 : planet.cells >= 5 ? 1.1 : 1;
    spreadAround(out, clamp(0.5 * atmo.spread * rot, 0.15, 0.9), (k) => Math.sin((k / LAT_STEPS) * Math.PI) + 1e-3);
    for (let k = 0; k <= LAT_STEPS; k++) out[k] += atmo.green;
    return out;
  }

  // Angle from the point under the star (map centre) over 180 deg, for a locked planet
  function lockedS(u, v) {
    const lon = (u - 0.5) * 2 * Math.PI, lat = (0.5 - v) * Math.PI;
    return Math.acos(clamp(Math.cos(lat) * Math.cos(lon), -1, 1)) / Math.PI;
  }

  function latTemp(lat) {
    if (!LAT_TEMP) return 1.0 - Math.pow(lat, EARTH_TEMP_EXP);
    const f = clamp(lat, 0, 1) * LAT_STEPS;
    const k = Math.min(LAT_STEPS - 1, Math.floor(f));
    return lerp(LAT_TEMP[k], LAT_TEMP[k + 1], f - k);
  }

  function isEarth(p) {
    return Math.abs(p.tilt - EARTH.tilt) < 0.01 && Math.abs(p.flux - 1) < 1e-4 && p.cells === 3
      && !p.locked && p.atmo === 'earth';
  }

  // Wind cell borders (lat 0..1) for c cells per hemisphere: where the belt wind turns.
  // c = 3 gives 0.25 and 0.75, the borders the generator always used.
  function cellBorders(c) {
    const b = [];
    for (let k = 0; k < c - 1; k++) b.push((k + 1) / c); // Earth (3): 30 and 60 degrees
    return b;
  }

  // +1 where surface winds meet and air rises (wet: the equator, polar fronts),
  // -1 where they part and air sinks (dry: subtropical highs, the poles)
  function convergence(c, lat) {
    const pts = [0, ...cellBorders(c), 1];
    let j = 0;
    while (j < pts.length - 2 && lat > pts[j + 1]) j++;
    const t = clamp((lat - pts[j]) / (pts[j + 1] - pts[j] + 1e-9), 0, 1);
    const a = j % 2 === 0 ? 1 : -1;
    return lerp(a, -a, (1 - Math.cos(Math.PI * t)) / 2);
  }

  // -------------------- Seasons layer --------------------
  // One moment of the year on top of the drawn map (no new map needed): seasonal snow on
  // land, sea ice, frozen lakes, and the dry season of savannas, shrubs and steppes.
  // lam = the planet's position in its year (solar longitude, 0 = spring equinox of the
  // north), or null for the year as a whole (then only ice that stays all year).
  // The swing of sunlight is damped: land warms and cools with a delay, the sea far less;
  // a thick atmosphere evens the seasons out, a thin one sharpens them.
  const SEASON_GAIN = 0.5;   // share of the sunlight swing that reaches the land temperature
  const SEASON_SEA = 0.25;   // the sea follows the seasons this much less (it stores heat)
  const SEASON_MAX = 0.45;   // the swing is capped: in the polar night sunlight drops to zero
  const SNOW_T = 0.11;       // land below this is under snow (the year-mean ice line is ~0.12)
  const ICE_T = 0.16;        // sea below this freezes (Earth: all year above ~82 deg, in winter to ~72 deg)
  const PERIHELION = 283 * Math.PI / 180; // Earth's: early January (northern winter)

  function dailyInsolation(phi, dec) {
    const p = clamp(phi, -Math.PI / 2 + 1e-4, Math.PI / 2 - 1e-4);
    const c = -Math.tan(p) * Math.tan(dec);
    const h0 = c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c);
    return (h0 * Math.sin(p) * Math.sin(dec) + Math.cos(p) * Math.cos(dec) * Math.sin(h0)) / Math.PI;
  }

  // Change of temperature from the year mean at solar longitude lam, by signed latitude
  // (index LAT_STEPS = equator, 0 = south pole, 2 * LAT_STEPS = north pole)
  function seasonDelta(planet, lam) {
    const { te, norm } = earthTau();
    const slope = 1 / (te[0] - te[LAT_STEPS]);
    const eps = planet.tilt * Math.PI / 180;
    const dec = Math.asin(Math.sin(eps) * Math.sin(lam));
    const e = planet.ecc || 0;
    // closer to the star = more light; relative to the year mean, which flux already holds
    const near = Math.pow((1 + e * Math.cos(lam - PERIHELION)) / (1 - e * e), 2) * Math.sqrt(1 - e * e);
    const atmo = ATMO[planet.atmo] || ATMO.earth;
    const out = new Float32Array(2 * LAT_STEPS + 1);
    for (let k = -LAT_STEPS; k <= LAT_STEPS; k++) {
      const phi = (k / LAT_STEPS) * Math.PI / 2;
      const now = Math.pow(planet.flux * near * dailyInsolation(phi, dec), 0.25);
      const year = Math.pow(planet.flux * annualInsolation(Math.abs(phi), eps), 0.25);
      out[k + LAT_STEPS] = clamp((now - year) / norm * slope * SEASON_GAIN, -SEASON_MAX, SEASON_MAX) * atmo.season;
    }
    return out;
  }

  function seasonLayer(lam) {
    const TM = window.TerrainMap;
    const temp = TM.lastTemp, ids = TM.lastBiomeIds;
    if (!temp || !ids) return null;
    const n = TM.lastSize, nh = TM.lastHeight, planet = TM.lastPlanet || EARTH, lake = TM.lastLake;
    const seasonal = lam != null && !planet.locked;
    const delta = seasonal ? seasonDelta(planet, lam) : null;
    const canvas = document.createElement('canvas');
    canvas.width = n; canvas.height = nh;
    const ctx = canvas.getContext('2d');
    const img = ctx.createImageData(n, nh);
    const d = img.data;
    const ice = new Uint8Array(n * nh); // 1 sea ice / frozen lake, 2 seasonal snow
    const tiltK = Math.min(1, planet.tilt / 15); // little tilt, little dry season
    const put = (i, c, a) => { const o = i * 4; d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = a * 255; };
    for (let y = 0; y < nh; y++) {
      const sl = 0.5 - y / (nh - 1);                       // signed latitude, + north
      const dRow = delta ? delta[Math.round((sl * 2 + 1) * LAT_STEPS)] : 0;
      // +1 summer, -1 winter in this hemisphere
      const summer = seasonal ? Math.sin(lam) * (sl >= 0 ? 1 : -1) * tiltK : 0;
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        const id = ids[i];
        if (id === 0 || (lake && lake[i])) {
          // sea (and lakes): ice where it is cold enough
          const t = temp[i] + dRow * (id === 0 ? SEASON_SEA : 0.7);
          const cover = 1 - smoothstep01(ICE_T - 0.025, ICE_T + 0.02, t);
          if (cover > 0.02) { put(i, SEA_ICE, cover * 0.92); if (cover > 0.5) ice[i] = 1; }
          continue;
        }
        if (id === BIOME.SNOW) continue; // ice all year already
        const t = temp[i] + dRow;
        const cover = 1 - smoothstep01(SNOW_T - 0.035, SNOW_T + 0.02, t);
        if (cover > 0.02) { put(i, ATLAS.snow, cover * 0.88); if (cover > 0.5) ice[i] = 2; continue; }
        if (!seasonal) continue;
        // dry seasons: savannas and seasonal forests dry out in winter (the rains follow the
        // sun), shrubs and steppes in summer (Mediterranean and continental droughts)
        let dry = 0;
        if (id === BIOME.SAVANNA || id === BIOME.SEASONAL_FOREST) dry = Math.max(0, -summer);
        else if (id === BIOME.SHRUBLAND) dry = Math.max(0, summer);
        else if (id === BIOME.STEPPE) dry = Math.max(0, summer) * 0.6;
        if (dry > 0.02) put(i, DRY_GRASS, dry * 0.4);
      }
    }
    ctx.putImageData(img, 0, 0);
    TM.lastIce = ice;
    return canvas;
  }
  const SEA_ICE = [224, 234, 242];
  const DRY_GRASS = [186, 156, 98];

  // -------------------- Climate --------------------

  // Latitude is bent by noise (warm and cold tongues, like ocean currents and air masses),
  // so climate zones are not straight horizontal bands.
  const TEMP_LAT_WARP = 0.16;
  const TEMP_LOCAL = 0.05;

  function computeTemperature(h, n, seaLevel, seed) {
    const nh = h.length / n; // map height (n is the width)
    const t = new Float32Array(n * nh);
    for (let y = 0; y < nh; y++) {
      const v = y / (nh - 1);
      for (let x = 0; x < n; x++) {
        const u = x / n; // periodic east-west
        const i = idx(x, y, n);
        const height = h[i];
        const warp = (wfbm(u, v, 3, seed + 800, 3, 0.5) - 0.5) * 2 * TEMP_LAT_WARP;
        // 0 equator -> 1 poles; on a locked planet: 0 under the star -> 1 the far side
        const lat = PLANET.locked
          ? clamp(lockedS(u, v) + warp * 0.5, 0, 1)
          : clamp(Math.abs(v - 0.5) * 2.0 + warp, 0, 1);
        // base temperature curve (nonlinear: wider tropics) + local variation
        const base = latTemp(lat) // Earth: 1 - lat^EARTH_TEMP_EXP; other planets: see Planet above
          + (wfbm(u, v, 8, seed + 810, 2, 0.5) - 0.5) * 2 * TEMP_LOCAL;

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

  // -------------------- Moisture / rain --------------------
  // Moisture is transported by the wind on a coarse grid (semi-Lagrangian advection),
  // iterated until it reaches a steady state:
  //   ocean      recharges the air (warm sea -> more moisture);
  //   land       rains out a share every step: base rate + orographic lift (air going
  //              uphill along the wind), much less on the lee side (rain shadow);
  //              part of the rain evaporates back (plants, lakes);
  //   cold air   holds less, the excess falls out.
  // Every cell gets its moisture from the whole upwind area, not from one traced line,
  // so the result follows coasts and mountains instead of wind streamlines.
  // Returns rain index 0..1 per pixel (ocean = 1).
  // Measured on 5 seeds, share of land pixels on a zone border: per-gully lift 14.7%,
  // no orography 5.7%, smoothed (6 cells, x2) 7.9%: ranges still cast rain shadows.
  const OROGRAPHY_SMOOTH = 6;  // coarse cells (1 cell = 4 px at 768)
  const OROGRAPHY_GAIN = 2;    // compensates the lower slopes of the smoothed height

  function computeRain(h, n, seaLevel, temp, wind, desertStrength) {
    const ds = clamp(desertStrength ?? 0.35, 0, 1);
    const nh = h.length / n;           // map height (n is the width)
    const cgh = Math.min(nh, 192);     // coarse grid: height, and square cells
    const scale = nh / cgh;
    const cg = Math.round(n / scale);  // coarse grid width (twice the height)
    const dt = 1.2;                    // coarse cells the air moves per step at speed 1

    const pBase = 0.003 + 0.016 * ds;  // more deserts -> land dries the air faster
    const kOro = 5;                    // lift is ~0.02 per step on the steepest 10% of slopes
    const kShadow = 45;                // descending air: base rate drops to ~30%
    const recycle = 0.3;
    const recharge = 0.25;
    const mix = 0.25;                  // light mixing with neighbours, removes thin filaments
    // desertStrength changes how far inland moisture gets, not how wet coasts are
    const pRef = 0.0035 * (pBase / (0.003 + 0.016 * 0.35));

    // Coarse height (sea level as floor) and temperature
    const H = new Float32Array(cg * cgh);
    const T = new Float32Array(cg * cgh);
    const land = new Uint8Array(cg * cgh);
    for (let y = 0; y < cgh; y++) {
      const y0 = Math.floor(y * scale), y1 = Math.max(y0 + 1, Math.floor((y + 1) * scale));
      for (let x = 0; x < cg; x++) {
        const x0 = Math.floor(x * scale), x1 = Math.max(x0 + 1, Math.floor((x + 1) * scale));
        let sh = 0, st = 0, c = 0;
        for (let yy = y0; yy < y1; yy++) {
          for (let xx = x0; xx < x1; xx++) {
            const i = idx(xx, yy, n);
            sh += Math.max(h[i], seaLevel);
            st += temp[i];
            c++;
          }
        }
        const j = idx(x, y, cg);
        H[j] = sh / c;
        T[j] = st / c;
        land[j] = H[j] > seaLevel + 1e-6 ? 1 : 0;
      }
    }

    // Wind displacement per step and orographic lift (wind . grad H)
    const WX = new Float32Array(cg * cgh);
    const WY = new Float32Array(cg * cgh);
    const LIFT = new Float32Array(cg * cgh);
    for (let y = 0; y < cgh; y++) {
      for (let x = 0; x < cg; x++) {
        const j = idx(x, y, cg);
        const w = sampleWindBilinear(wind, x / cg, y / (cgh - 1));
        WX[j] = w.vx * w.spd * dt;
        WY[j] = w.vy * w.spd * dt;
      }
    }
    // Lift is taken from a smoothed height: rain shadows belong to whole ranges,
    // not to every gully (per-gully lift made wet/dry stripes on every eroded slope)
    const oroR = OROGRAPHY_SMOOTH;
    const HS = oroR > 0 ? blurred(H, cg, oroR) : H;
    const oroK = oroR > 0 ? OROGRAPHY_GAIN : 1;
    for (let y = 0; y < cgh; y++) {
      for (let x = 0; x < cg; x++) {
        const xm = wrapX(x - 1, cg), xp = wrapX(x + 1, cg);
        const ym = Math.max(0, y - 1), yp = Math.min(cgh - 1, y + 1);
        const gx = (HS[idx(xp, y, cg)] - HS[idx(xm, y, cg)]) * 0.5 * oroK;
        const gy = (HS[idx(x, yp, cg)] - HS[idx(x, ym, cg)]) * 0.5 * oroK;
        const j = idx(x, y, cg);
        LIFT[j] = WX[j] * gx + WY[j] * gy;
      }
    }

    let Q = new Float32Array(cg * cgh);   // moisture in the air
    let Qn = new Float32Array(cg * cgh);
    const PR = new Float32Array(cg * cgh); // precipitation per step
    // rising or sinking air by latitude for the planet's wind cells, relative to Earth
    // (all 1 on Earth-like planets)
    // (a locked planet: air rises under the star and sinks on the far side)
    const belt = new Float32Array(cg * cgh).fill(1);
    if (PLANET.locked) {
      for (let y = 0; y < cgh; y++) {
        const lat = Math.abs(y / (cgh - 1) - 0.5) * 2;
        for (let x = 0; x < cg; x++) {
          const s = lockedS(x / cg, y / (cgh - 1));
          const conv = clamp(Math.cos(s * Math.PI) * 1.6 - 0.3, -1, 1);
          belt[idx(x, y, cg)] = Math.exp(BELT_RAIN * 1.3 * (conv - convergence(EARTH.cells, lat)));
        }
      }
    } else if (PLANET.cells !== EARTH.cells) {
      for (let y = 0; y < cgh; y++) {
        const lat = Math.abs(y / (cgh - 1) - 0.5) * 2;
        const b = Math.exp(BELT_RAIN * (convergence(PLANET.cells, lat) - convergence(EARTH.cells, lat)));
        for (let x = 0; x < cg; x++) belt[idx(x, y, cg)] = b;
      }
    }
    const evap = (j) => 0.35 + 0.65 * T[j];
    const capacity = (j) => 0.2 + 0.8 * T[j];
    for (let j = 0; j < Q.length; j++) Q[j] = land[j] ? 0 : evap(j);

    function sampleQ(u, v) {
      // air from beyond a pole is open-ocean air; east-west it comes round the planet
      if (v < 0 || v > cgh - 1) return 0.6;
      const fx = Math.floor(u), y0 = Math.floor(v);
      const tx = u - fx, ty = v - y0;
      const x0 = wrapX(fx, cg), x1 = (x0 + 1) % cg, y1 = Math.min(y0 + 1, cgh - 1);
      return lerp(
        lerp(Q[idx(x0, y0, cg)], Q[idx(x1, y0, cg)], tx),
        lerp(Q[idx(x0, y1, cg)], Q[idx(x1, y1, cg)], tx),
        ty
      );
    }

    // Enough steps for air to cross the whole map
    const iterations = Math.ceil((cgh * 1.3) / dt);
    for (let it = 0; it < iterations; it++) {
      for (let y = 0; y < cgh; y++) {
        for (let x = 0; x < cg; x++) {
          const j = idx(x, y, cg);
          let q = sampleQ(x - WX[j], y - WY[j]); // where this air was one step ago
          const nb = (
            Q[idx(wrapX(x - 1, cg), y, cg)] + Q[idx(wrapX(x + 1, cg), y, cg)] +
            Q[idx(x, Math.max(0, y - 1), cg)] + Q[idx(x, Math.min(cgh - 1, y + 1), cg)]
          ) * 0.25;
          q = q * (1 - mix) + nb * mix;

          if (!land[j]) {
            q += (evap(j) - q) * recharge;
            PR[j] = 0;
          } else {
            const lift = LIFT[j];
            let rate = pBase + (lift > 0 ? kOro * lift : 0);
            if (lift < 0) rate *= Math.exp(kShadow * lift);
            rate *= belt[j];
            let p = q * clamp(rate, 0, 0.9);
            const excess = q - p - capacity(j);
            if (excess > 0) p += excess;
            q -= p * (1 - recycle);
            PR[j] = p;
          }
          Qn[j] = q;
        }
      }
      const tmp = Q; Q = Qn; Qn = tmp;
    }

    // Back to full resolution (bilinear over land cells only) and to a 0..1 index
    const rain = new Float32Array(n * nh);
    for (let y = 0; y < nh; y++) {
      const v = clamp((y + 0.5) / scale - 0.5, 0, cgh - 1);
      const y0 = Math.floor(v), y1 = Math.min(y0 + 1, cgh - 1), ty = v - y0;
      for (let x = 0; x < n; x++) {
        const i = idx(x, y, n);
        if (h[i] <= seaLevel) { rain[i] = 1; continue; }
        const u = (x + 0.5) / scale - 0.5;       // wraps east-west
        const fx = Math.floor(u), tx = u - fx;
        const x0 = wrapX(fx, cg), x1 = (x0 + 1) % cg;
        let s = 0, wsum = 0;
        const corners = [
          [idx(x0, y0, cg), (1 - tx) * (1 - ty)], [idx(x1, y0, cg), tx * (1 - ty)],
          [idx(x0, y1, cg), (1 - tx) * ty],       [idx(x1, y1, cg), tx * ty],
        ];
        for (const [j, wgt] of corners) {
          if (land[j]) { s += PR[j] * wgt; wsum += wgt; }
        }
        const p = wsum > 0 ? s / wsum : 0;
        rain[i] = 1 - Math.exp(-p / pRef);
      }
    }
    return rain;
  }


  // -------------------- Rendering helpers --------------------

  function hillshade(h, n, x, y) {
    const nh = h.length / n; // map height (n is the width)
    const xm1 = wrapX(x - 1, n);
    const xp1 = wrapX(x + 1, n);
    const ym1 = clamp(y - 1, 0, nh - 1);
    const yp1 = clamp(y + 1, 0, nh - 1);

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
  // -------------------- Atlas style rendering --------------------
  // Hand-drawn atlas look: muted biome colours with soft borders, elevation as stacked
  // translucent "ink" layers (like layers in an image editor), very gentle hillshade,
  // saturated ocean with a light glow along the coast, periwinkle rivers.

  const ATLAS = {
    lowland: [158, 194, 132],      // height mode base (no biomes)
    ink: [62, 66, 44],             // elevation layers
    inkAlpha: 0.12,
    layers: [0.08, 0.18, 0.30, 0.44, 0.60, 0.78],
    snow: [238, 238, 244],
    // Sea depth steps, like a nautical chart: depth = (sea - h) / sea, band up to `to`
    depthBands: [
      { to: 0.03, c: [96, 178, 192] },   // shallows
      { to: 0.08, c: [70, 150, 194] },   // shelf
      { to: 0.15, c: [52, 120, 186] },
      // deep steps close in colour: open sea stays calm instead of blotchy
      { to: 0.25, c: [40, 92, 168] },
      { to: 0.40, c: [36, 82, 158] },
      { to: 1.01, c: [32, 72, 148] },    // deep sea
    ],
    shore: [226, 232, 190],        // pale edge on the land side
    river: [56, 118, 188],         // same blue family as the sea
    lake: [74, 152, 196],          // like the shelf band
  };

  // Seabed depth as a share of sea level. Measured near coasts: 1-3 px out the median
  // depth is ~0.03, 8-12 px out it ranges ~0.08..0.28, so the shelf width varies a lot.
  const SHORE_SHALLOW_DEPTH = 0.06; // pale shore on land only next to water this shallow

  const LAKE_MIN_DEPTH = 0.006; // filled depth (height units) to count as a lake
  // Lakes are rare on Earth (about 2% of land, mostly in glacial or rift country), so
  // only the deepest, largest basins are kept until this share of the land is water
  const LAKE_SHARE = 0.012;
  const LAKE_MIN_RAIN = 0.3;    // drier basins stay dry
  const LAKE_MIN_AREA = 30;     // pixels at size 768 (scaled with map area); smaller ones are noise

  // -------------------- Meandering rivers --------------------
  // Every river is traced from its source along the flow directions until it joins
  // another river, a lake or the sea. The cell path is smoothed (Chaikin corner
  // cutting), then bent sideways by a wave whose wavelength and amplitude grow with
  // the river's size; on steep ground the bends fade out. Ends stay fixed so
  // tributaries meet their river. Lines are drawn anti-aliased, wider downstream.
  const MEANDER_AMP = 3.8;        // px at 768, for a big lowland river
  const MEANDER_WAVE = 20;        // px at 768, for a big river
  const MEANDER_FLAT_SLOPE = 0.006; // above this slope (per px at 768) bends fade out

  function chaikin(pts, passes) {
    let p = pts;
    for (let k = 0; k < passes; k++) {
      if (p.length < 3) return p;
      const out = [p[0]];
      for (let i = 0; i < p.length - 1; i++) {
        const a = p[i], b = p[i + 1];
        out.push({ x: a.x * 0.75 + b.x * 0.25, y: a.y * 0.75 + b.y * 0.25, f: a.f * 0.75 + b.f * 0.25 });
        out.push({ x: a.x * 0.25 + b.x * 0.75, y: a.y * 0.25 + b.y * 0.75, f: a.f * 0.25 + b.f * 0.75 });
      }
      out.push(p[p.length - 1]);
      p = out;
    }
    return p;
  }

  function rasterizeRivers(h, acc, down, n, seaLevel, threshold, lake, openSea, seed) {
    const nh = h.length / n; // map height (n is the width)
    const N = n * nh;
    const px = nh / 768;
    const out = new Float32Array(N);
    if (!down) return out;

    const isRiver = new Uint8Array(N);
    let accMax = 1;
    for (let i = 0; i < N; i++) {
      if (openSea[i] || lake[i] || acc[i] <= threshold) continue;
      isRiver[i] = 1;
      if (acc[i] > accMax) accMax = acc[i];
    }
    const logThr = Math.log(Math.max(2, threshold));
    const logMax = Math.log(accMax + 1e-9);
    const flowT = (a) => clamp((Math.log(Math.max(a, 1)) - logThr) / (logMax - logThr + 1e-9), 0, 1);

    const hasDonor = new Uint8Array(N);
    for (let i = 0; i < N; i++) if (isRiver[i] && down[i] >= 0 && isRiver[down[i]]) hasDonor[down[i]] = 1;

    const canvas = document.createElement('canvas');
    canvas.width = n; canvas.height = nh;
    const g = canvas.getContext('2d');
    g.strokeStyle = '#fff';
    g.lineCap = 'round';
    g.lineJoin = 'round';

    const visited = new Uint8Array(N);
    const rng = mulberry32(seed + 4242);
    for (let s = 0; s < N; s++) {
      if (!isRiver[s] || hasDonor[s]) continue;
      // trace downstream
      const cells = [];
      let i = s;
      while (i >= 0 && isRiver[i] && !visited[i]) {
        visited[i] = 1;
        cells.push(i);
        i = down[i];
      }
      if (i >= 0) cells.push(i); // the junction, lake shore or sea cell we ran into
      if (cells.length < 2) continue;

      // a river may cross the east-west seam: unwrap x so the path stays continuous
      // (it is drawn three times below, shifted by the map width, and clipped by the canvas)
      let pts = [];
      for (const c of cells) {
        let x = (c % n) + 0.5;
        if (pts.length) {
          const prev = pts[pts.length - 1].x;
          x += Math.round((prev - x) / n) * n;
        }
        pts.push({ x, y: ((c / n) | 0) + 0.5, f: flowT(acc[c]) });
      }
      pts = chaikin(pts, 2);

      // length along the path, for the wave phase and the end taper
      const len = [0];
      for (let k = 1; k < pts.length; k++) {
        len.push(len[k - 1] + Math.hypot(pts[k].x - pts[k - 1].x, pts[k].y - pts[k - 1].y));
      }
      const total = len[len.length - 1];
      let phase = rng() * Math.PI * 2;
      const bent = pts.map((p, k) => {
        const a = pts[Math.max(0, k - 1)], b = pts[Math.min(pts.length - 1, k + 1)];
        let tx = b.x - a.x, ty = b.y - a.y;
        const tl = Math.hypot(tx, ty) || 1;
        tx /= tl; ty /= tl;
        const cx = wrapX(Math.floor(p.x), n), cy = clamp(p.y | 0, 0, nh - 1);
        const slope = Math.hypot(
          h[idx(wrapX(cx + 1, n), cy, n)] - h[idx(wrapX(cx - 1, n), cy, n)],
          h[idx(cx, Math.min(nh - 1, cy + 1), n)] - h[idx(cx, Math.max(0, cy - 1), n)]
        ) * 0.5 / px;
        const flat = 1 - smoothstep(clamp(slope / MEANDER_FLAT_SLOPE, 0, 1));
        const wave = (8 + MEANDER_WAVE * p.f) * px;
        if (k > 0) phase += ((len[k] - len[k - 1]) / wave) * Math.PI * 2 * (0.75 + 0.5 * valueNoise2(len[k] * 0.05, s * 0.01, seed + 4300));
        const taper = clamp(Math.min(len[k], total - len[k]) / (wave * 0.5), 0, 1);
        const amp = (0.6 + MEANDER_AMP * p.f) * px * flat * taper;
        const off = Math.sin(phase) * amp;
        return { x: p.x - ty * off, y: p.y + tx * off, f: p.f };
      });

      for (const shift of [-n, 0, n]) {
        for (let k = 1; k < bent.length; k++) {
          const a = bent[k - 1], b = bent[k];
          if (Math.max(a.x, b.x) + shift < -8 || Math.min(a.x, b.x) + shift > n + 8) continue;
          g.lineWidth = (0.9 + 2.1 * b.f) * px;
          g.beginPath();
          g.moveTo(a.x + shift, a.y);
          g.lineTo(b.x + shift, b.y);
          g.stroke();
        }
      }
    }

    const data = g.getImageData(0, 0, n, nh).data;
    for (let i = 0; i < N; i++) out[i] = data[i * 4 + 3] / 255;
    return out;
  }

  // 1 for open sea: connected water (4-neighbour, wrapping east-west) at or below sea level
  // that reaches a pole, or is large. Now that east and west wrap, an ocean may circle the
  // planet without touching any map edge; small enclosed basins stay lakes.
  const OPEN_SEA_SHARE = 0.02; // a water body this large (share of the map) is sea anyway

  function connectedToEdge(h, n, seaLevel) {
    const nh = h.length / n; // map height (n is the width)
    const N = n * nh;
    const sea = new Uint8Array(N);
    const seen = new Uint8Array(N);
    const queue = new Int32Array(N);
    const big = Math.max(64, N * OPEN_SEA_SHARE);
    for (let s = 0; s < N; s++) {
      if (seen[s] || h[s] > seaLevel) continue;
      let qh = 0, qt = 0, pole = false;
      seen[s] = 1;
      queue[qt++] = s;
      while (qh < qt) {
        const i = queue[qh++];
        const x = i % n, y = (i / n) | 0;
        if (y === 0 || y === nh - 1) pole = true;
        for (const [dx, dy] of [[-1, 0], [1, 0], [0, -1], [0, 1]]) {
          const j = nbr(x, y, dx, dy, n, nh);
          if (j >= 0 && !seen[j] && h[j] <= seaLevel) { seen[j] = 1; queue[qt++] = j; }
        }
      }
      if (pole || qt >= big) for (let k = 0; k < qt; k++) sea[queue[k]] = 1;
    }
    return sea;
  }

  // Remove connected regions (8-neighbour) smaller than minArea from a 0/1 mask
  // Keep the lakes that matter: regions are ranked by depth x size and kept until their
  // total area reaches `budget` cells; shallow puddles in every hollow are dropped
  function limitLakes(mask, hf, h, n, budget, seaLevel) {
    const nh = mask.length / n;
    const seen = new Uint8Array(mask.length);
    const stack = new Int32Array(mask.length);
    const regions = [];
    for (let s = 0; s < mask.length; s++) {
      if (!mask[s] || seen[s]) continue;
      const cells = [];
      let deepest = 0, top = 0;
      stack[top++] = s;
      seen[s] = 1;
      while (top > 0) {
        const i = stack[--top];
        cells.push(i);
        deepest = Math.max(deepest, hf[i] - h[i], seaLevel - h[i]);
        const x = i % n, y = (i / n) | 0;
        for (let k = 0; k < 8; k++) {
          const j = nbr(x, y, N8[k][0], N8[k][1], n, nh);
          if (j >= 0 && mask[j] && !seen[j]) { seen[j] = 1; stack[top++] = j; }
        }
      }
      regions.push({ cells, score: deepest * Math.sqrt(cells.length) });
    }
    regions.sort((a, b) => b.score - a.score);
    let total = 0;
    for (const r of regions) {
      // the first (deepest, largest) lake may go somewhat over the budget, a giant may not
      const limit = total === 0 ? budget * 2.5 : budget;
      if (total + r.cells.length > limit) { for (const i of r.cells) mask[i] = 0; continue; }
      total += r.cells.length;
    }
  }

  // Transition belts between zones. On Earth a desert does not meet a dense forest:
  // between them lie savanna, steppe or shrubs; and lowland ice gives way to tundra
  // before the taiga. Forest cells within `r` px of a desert become the open zone that
  // fits their climate; forest cells near lowland ice become tundra.
  function addTransitions(ids, h, n, seaLevel, r) {
    const N = ids.length;
    const desert = new Uint8Array(N), ice = new Uint8Array(N);
    for (let i = 0; i < N; i++) {
      const id = ids[i];
      if (id === BIOME.HOT_DESERT || id === BIOME.COLD_DESERT) desert[i] = 1;
      // ice low down = polar ice cap (ice on high mountains already has rock and meadows around it)
      else if (id === BIOME.SNOW && (h[i] - seaLevel) / (1 - seaLevel + 1e-9) < 0.45) ice[i] = 1;
    }
    const toDesert = bfsDistance(desert, n, r);
    const toIce = bfsDistance(ice, n, r);
    for (let i = 0; i < N; i++) {
      const id = ids[i];
      const nearIce = toIce[i] > 0;
      const nearDesert = toDesert[i] > 0;
      if (nearIce && (id === BIOME.TAIGA || id === BIOME.TEMPERATE_FOREST || id === BIOME.WET_TEMPERATE)) {
        ids[i] = BIOME.TUNDRA;
      } else if (nearDesert) {
        if (id === BIOME.RAINFOREST || id === BIOME.SEASONAL_FOREST) ids[i] = BIOME.SAVANNA;
        else if (id === BIOME.TEMPERATE_FOREST || id === BIOME.WET_TEMPERATE) ids[i] = BIOME.SHRUBLAND;
        else if (id === BIOME.TAIGA) ids[i] = BIOME.STEPPE;
      }
    }
  }

  // Volcanic island chains over hotspots (like Hawaii): a plate slides over a fixed plume,
  // so a line of islands forms, the youngest high and the older ones worn down, sinking
  // into seamounts. Placed only in open ocean, away from the coasts.
  function addIslandChains(h, n, seed, seaLevel) {
    const nh = h.length / n;
    const rng = mulberry32(seed + 4242);
    const chains = 3 + ((rng() * 4) | 0);
    for (let c = 0; c < chains; c++) {
      let u = rng(), v = 0.12 + rng() * 0.76;
      const dir = rng() * Math.PI * 2;
      const count = 4 + ((rng() * 6) | 0);
      const step = 0.018 + rng() * 0.016;           // v units between islands
      const r0 = 0.007 + rng() * 0.007;             // radius of the youngest island
      for (let k = 0; k < count; k++) {
        const age = k / Math.max(1, count - 1);     // 0 youngest .. 1 oldest
        const cx = wrapX(Math.round(u * n), n), cy = Math.round(v * (nh - 1));
        u += Math.cos(dir) * step / ASPECT + (rng() - 0.5) * 0.006;
        v += Math.sin(dir) * step + (rng() - 0.5) * 0.006;
        if (cy < 2 || cy > nh - 3) break;
        if (h[idx(cx, cy, n)] > seaLevel - 0.05) continue; // no islands on shelves or land
        const peak = seaLevel + (0.07 - 0.11 * age) * (0.7 + rng() * 0.6); // old ones drown
        const rad = r0 * (1 - 0.45 * age) * (0.8 + rng() * 0.5) * (nh - 1); // px
        const R = Math.ceil(rad * 2.2);
        for (let dy = -R; dy <= R; dy++) {
          const y = cy + dy;
          if (y < 0 || y >= nh) continue;
          for (let dx = -R; dx <= R; dx++) {
            const x = wrapX(cx + dx, n);
            const i = idx(x, y, n);
            // irregular outline: the distance is bent by noise
            const wob = 1 + 0.35 * (wfbm(x / n, y / (nh - 1), 40, seed + 4300 + c, 2, 0.5) - 0.5) * 2;
            const d = Math.hypot(dx, dy) / (rad * wob);
            const w = Math.exp(-d * d * 1.6);
            if (w < 0.01) continue;
            h[i] = Math.max(h[i], lerp(h[i], peak, w));
          }
        }
      }
    }
  }

  function dropSmallRegions(mask, n, minArea) {
    const nh = mask.length / n; // map height (n is the width)
    const seen = new Uint8Array(mask.length);
    const stack = new Int32Array(mask.length);
    const region = [];
    for (let s = 0; s < mask.length; s++) {
      if (!mask[s] || seen[s]) continue;
      let top = 0;
      stack[top++] = s;
      seen[s] = 1;
      region.length = 0;
      while (top > 0) {
        const i = stack[--top];
        region.push(i);
        const x = i % n, y = (i / n) | 0;
        for (let k = 0; k < 8; k++) {
          const j = nbr(x, y, N8[k][0], N8[k][1], n, nh);
          if (j >= 0 && mask[j] && !seen[j]) { seen[j] = 1; stack[top++] = j; }
        }
      }
      if (region.length < minArea) for (const i of region) mask[i] = 0;
    }
  }

  // Zone colours: told apart by lightness as well as hue, in families. Forests are the
  // darker greens (rainforest darkest, taiga bluish, temperate rainforest teal), open dry
  // land is light and warm (steppe yellow-green, savanna straw, deserts sand / grey-beige),
  // cold and high ground is greyish (tundra sage, rock brownish, meadows fresh light green).
  const ATLAS_BIOMES = {
    1: [218, 196, 144],  // hot desert: sand
    2: [192, 180, 160],  // cold desert / semi-arid: grey-beige
    3: [206, 186, 104],  // savanna: straw
    4: [50, 88, 82],     // taiga: dark blue-green
    5: [152, 166, 150],  // tundra: cool sage grey
    6: ATLAS.snow,       // snow and ice
    7: [178, 198, 122],  // steppe / grassland: light yellow-green
    8: [108, 152, 72],   // seasonal forest: yellow-green
    9: [72, 130, 74],    // temperate forest: mid green
    10: [26, 86, 56],    // rainforest: darkest green
    11: [42, 108, 98],   // temperate rainforest: teal
    12: [154, 156, 96],  // shrubland: olive
    13: [166, 190, 142], // alpine meadow: fresh light green
    14: [134, 122, 112], // bare rock: brownish grey
    15: [92, 140, 146],  // wetland: blue-teal
  };

  // Ragged borders: colours and heights are read from a point shifted by fine noise,
  // so every border becomes jagged and crisp (no blending, no smudging).
  // Amplitudes in px at size 768.
  const BORDER_JAG = 3;      // ~25 px features
  const BORDER_JAG_FINE = 1.6; // ~8 px features

  function smoothstep01(a, b, x) {
    const t = clamp((x - a) / (b - a), 0, 1);
    return t * t * (3 - 2 * t);
  }

function drawMap(ctx, h, acc, temp, rain, biomeIds, n, seaLevel, riverThreshold, mode, desertStrength, hf, seed, down) {
  const nh = h.length / n; // map height (n is the width)
  const img = ctx.createImageData(n, nh);
  const data = img.data;
  const N = n * nh;
  const px = nh / 768; // scale effects with map size

  // Sea = water below sea level that is connected to the map edge. Basins inside a
  // continent that dip below sea level are lakes, not sea (otherwise they got ocean
  // colours and a shelf glow inside lakes).
  const openSea = connectedToEdge(h, n, seaLevel);

  // Lakes: filled depressions that are deep enough and not in dry climates
  // (in deserts they stay dry basins), plus every inland basin below sea level
  const lake = new Float32Array(N);
  if (hf) {
    for (let i = 0; i < N; i++) {
      if (h[i] <= seaLevel) continue;
      const depth = hf[i] - h[i];
      if (depth > LAKE_MIN_DEPTH && rain[i] > LAKE_MIN_RAIN) lake[i] = 1;
    }
    dropSmallRegions(lake, n, Math.round(LAKE_MIN_AREA * px * px));
  }
  // inland water below sea level joins the other lake candidates (the budget decides)
  for (let i = 0; i < N; i++) if (h[i] <= seaLevel && !openSea[i]) lake[i] = 1;
  if (hf) {
    let landCells = 0;
    for (let i = 0; i < N; i++) if (h[i] > seaLevel) landCells++;
    limitLakes(lake, hf, h, n, Math.round(LAKE_SHARE * landCells), seaLevel);
  }
  const lakeSoft = blurred(lake, n, Math.max(1, Math.round(1 * px)));
  window.TerrainMap.lastLake = lake; // for the hover tooltip

  // Rivers: traced as smooth meandering lines (see rasterizeRivers), 0..1 coverage
  const river = rasterizeRivers(h, acc, down, n, seaLevel, riverThreshold, lake, openSea, seed);

  // Land mask and per-pixel base colour
  const land = new Float32Array(N);
  const R = new Float32Array(N), G = new Float32Array(N), B = new Float32Array(N);
  const snow = new Float32Array(N);
  for (let i = 0; i < N; i++) {
    if (openSea[i]) continue;
    land[i] = 1; // includes inland water, drawn as lake on top
    const c = mode === 'height' ? ATLAS.lowland : (ATLAS_BIOMES[biomeIds[i]] || ATLAS.lowland);
    R[i] = c[0]; G[i] = c[1]; B[i] = c[2];
    if (mode !== 'height' && biomeIds[i] === 6) snow[i] = 1;
  }

  const snowSoft = blurred(snow, n, Math.max(1, Math.round(2 * px)));

  // Softened zone colours (land only, so coasts don't pull in black): borders between
  // zones get a narrow blend instead of a hard cut, zones themselves stay flat
  const ZONE_BLEND = 0.55;
  const blendR = Math.max(1, Math.round(2 * px));
  const landSoft = blurred(land, n, blendR);
  const RS = blurred(R, n, blendR), GS = blurred(G, n, blendR), BS = blurred(B, n, blendR);
  const zones = mode !== 'height';

  // Heights for the elevation layers: almost raw, so contours stay crisp
  const hc = blurred(h, n, 1);

  // Jitter offset per pixel: colour and height are taken from a nearby point,
  // which turns every border into a crisp jagged line
  const jx = new Int32Array(N);
  const jy = new Int32Array(N);
  const jag = BORDER_JAG * px, jagF = BORDER_JAG_FINE * px;
  for (let y = 0; y < nh; y++) {
    const v = y / (nh - 1);
    for (let x = 0; x < n; x++) {
      const i = idx(x, y, n);
      const u = x / n; // periodic east-west
      const ox = (wfbm(u, v, 30, seed + 9110, 2, 0.5) - 0.5) * 2 * jag
               + (wfbm(u, v, 95, seed + 9111, 1, 0.5) - 0.5) * 2 * jagF;
      const oy = (wfbm(u, v, 30, seed + 9120, 2, 0.5) - 0.5) * 2 * jag
               + (wfbm(u, v, 95, seed + 9121, 1, 0.5) - 0.5) * 2 * jagF;
      const sx = wrapX(Math.round(x + ox), n);
      const sy = clamp(Math.round(y + oy), 0, nh - 1);
      const j = idx(sx, sy, n);
      // land never pulls a colour from the sea (water may pull from land: depth 0)
      jx[i] = land[i] && !land[j] ? x : sx;
      jy[i] = land[i] && !land[j] ? y : sy;
    }
  }

  // Sea depth (follows the real seabed, so shelves are wide where the bottom is shallow
  // far out and absent off steep coasts); a light blur keeps band edges from speckling
  const hSea = blurred(h, n, Math.max(1, Math.round(2 * px)));
  const seaDepth = (j) => (land[j] ? 0 : clamp((seaLevel - hSea[j]) / (seaLevel + 1e-9), 0, 1));
  const shelf = new Float32Array(N); // 1 where the water is shallow
  for (let i = 0; i < N; i++) {
    if (!land[i] && seaDepth(i) < SHORE_SHALLOW_DEPTH) shelf[i] = 1;
  }
  // On land: how shallow the nearby water is (beaches on shallow coasts only),
  // and how close the water is
  const nearShelf = blurred(shelf, n, Math.max(1, Math.round(3 * px)));
  const nearWater = blurred(land, n, Math.max(1, Math.round(3 * px)));


  for (let y = 0; y < nh; y++) {
    for (let x = 0; x < n; x++) {
      const i = idx(x, y, n);
      let r, g, b;

      if (!land[i]) {
        // Ocean: stepped depth bands with jagged edges, a slight shade inside each band
        const depth = seaDepth(idx(jx[i], jy[i], n));
        const bands = ATLAS.depthBands;
        let k = 0;
        while (k < bands.length - 1 && depth > bands[k].to) k++;
        const from = k > 0 ? bands[k - 1].to : 0;
        const inBand = clamp((depth - from) / (bands[k].to - from), 0, 1);
        const c = bands[k].c, cn = bands[Math.min(k + 1, bands.length - 1)].c;
        const a = inBand * 0.18;
        r = lerp(c[0], cn[0], a); g = lerp(c[1], cn[1], a); b = lerp(c[2], cn[2], a);
      } else {
        // Colour from the jittered point: crisp, jagged borders
        const j = idx(jx[i], jy[i], n);
        r = R[j]; g = G[j]; b = B[j];
        if (zones && landSoft[j] > 0.01) {
          r = lerp(r, RS[j] / landSoft[j], ZONE_BLEND);
          g = lerp(g, GS[j] / landSoft[j], ZONE_BLEND);
          b = lerp(b, BS[j] / landSoft[j], ZONE_BLEND);
        }

        // Elevation: stacked translucent ink layers with crisp jagged contours, no rim
        const alt = clamp((hc[j] - seaLevel) / (1 - seaLevel + 1e-9), 0, 1);
        const ew = 0.002;
        const inkK = 1 - 0.75 * snowSoft[i]; // keep snow fields light
        for (const level of ATLAS.layers) {
          const cover = smoothstep01(level - ew, level + ew, alt);
          if (cover <= 0) break;
          const a = ATLAS.inkAlpha * cover * inkK;
          r = lerp(r, ATLAS.ink[0], a); g = lerp(g, ATLAS.ink[1], a); b = lerp(b, ATLAS.ink[2], a);
        }

        // Relief shading (flat ground has hillshade ~0.7)
        const sh = hillshade(h, n, x, y);
        const shade = clamp(1 + (sh - 0.7) * 0.72, 0.76, 1.14); // a bit stronger: relief reads better
        r *= shade; g *= shade; b *= shade;

        // Pale edge along the shore
        const shoreK = (1 - smoothstep01(0.55, 0.95, nearWater[i])) * 0.55 * clamp(nearShelf[i] * 2.5, 0, 1);
        r = lerp(r, ATLAS.shore[0], shoreK); g = lerp(g, ATLAS.shore[1], shoreK); b = lerp(b, ATLAS.shore[2], shoreK);

        // Lakes
        if (lakeSoft[i] > 0) {
          const a = smoothstep01(0.15, 0.6, lakeSoft[i]) * 0.92;
          r = lerp(r, ATLAS.lake[0], a); g = lerp(g, ATLAS.lake[1], a); b = lerp(b, ATLAS.lake[2], a);
        }

        // Rivers: the core line, and a wider soft band for strong rivers
        const rv = river[i] * 0.92 * (1 - lakeSoft[i]);
        if (rv > 0) {
          r = lerp(r, ATLAS.river[0], rv); g = lerp(g, ATLAS.river[1], rv); b = lerp(b, ATLAS.river[2], rv);
        }
      }

      const o = i * 4;
      data[o + 0] = clamp(r, 0, 255) | 0;
      data[o + 1] = clamp(g, 0, 255) | 0;
      data[o + 2] = clamp(b, 0, 255) | 0;
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


function drawWindDebug(ctx, wind, n, nh) {
  ctx.clearRect(0, 0, n, nh);

  // лёгкая сетка линий
  ctx.fillStyle = '#0b1020';
  ctx.fillRect(0, 0, n, nh);

  const stepPx = Math.max(14, (nh / 48) | 0);
  ctx.strokeStyle = 'rgba(120, 255, 160, 0.65)';
  ctx.lineWidth = 1;

  for (let y = stepPx; y < nh; y += stepPx) {
    for (let x = stepPx; x < n; x += stepPx) {
      const u = x / n; // periodic east-west
      const v = y / (nh - 1);
      const w = sampleWindBilinear(wind, u, v);

      // speed 0.05..1.4: calm belts get short strokes
      const len = 4 + 22 * w.spd / 1.4;
      const x2 = x + w.vx * len;
      const y2 = y + w.vy * len;

      ctx.beginPath();
      ctx.moveTo(x, y);
      ctx.lineTo(x2, y2);
      ctx.stroke();
    }
  }
}

// Biome ids (colours in ATLAS_BIOMES, names in i18n.js as biome.<id>)
const BIOME = {
  OCEAN: 0, HOT_DESERT: 1, COLD_DESERT: 2, SAVANNA: 3, TAIGA: 4, TUNDRA: 5, SNOW: 6,
  STEPPE: 7, SEASONAL_FOREST: 8, TEMPERATE_FOREST: 9, RAINFOREST: 10, WET_TEMPERATE: 11,
  SHRUBLAND: 12, ALPINE: 13, ROCK: 14, WETLAND: 15,
};

// Altitudinal zonation: a mountain rising above its surroundings is wetter than the
// plain around it (more orographic rain, less evaporation in cool air), so its slopes
// get belts: desert -> steppe/shrubs -> forest -> alpine meadows -> rock -> snow.
// Driven by height above the regional average (not absolute height), so wide high
// plateaus can stay dry, like Tibet or the Altiplano.
const MOUNTAIN_WET = 0.55;  // extra effective moisture on a strongly raised slope
const RISE_FULL = 0.22;     // height above the surroundings (share of land range) for full effect

// slope: height difference to neighbours (height units per pixel at 768)
// wet: the cell is on a big, slow river or by a lake (flat, lots of water)
// rise: height above the regional average, as a share of the land height range
function biomeId(height, seaLevel, temp, rain, desertStrength, slope, wet, rise) {
  if (height <= seaLevel) return BIOME.OCEAN;

  const alt = (height - seaLevel) / (1 - seaLevel + 1e-9);
  const ds = clamp(desertStrength ?? 0.35, 0, 1);
  const k = ds * ds;

  const desertCut = 0.18 + 0.18 * k;
  const steppeCut = 0.36 + 0.12 * k;
  const shrubCut = steppeCut + 0.10;

  if (temp < 0.12 || alt > 0.86) return BIOME.SNOW;

  // effective moisture: rain plus the mountain bonus
  const t = clamp((rise || 0) / RISE_FULL, 0, 1);
  const moist = rain + MOUNTAIN_WET * t * t * (3 - 2 * t);

  // high mountains: bare rock on steep or dry slopes, meadows elsewhere
  if (alt > 0.68 || (alt > 0.45 && slope > 0.03)) return BIOME.ROCK;
  if (alt > 0.5) return moist < 0.3 ? BIOME.ROCK : BIOME.ALPINE;

  if (wet && temp > 0.25 && rain > 0.3) return BIOME.WETLAND;

  if (temp < 0.2) return BIOME.TUNDRA;
  if (temp < 0.34) return moist < 0.28 ? BIOME.TUNDRA : BIOME.TAIGA;

  const hot = temp > 0.56; // ~35 deg on lowland with the linear Earth curve
  if (moist < desertCut) return hot ? BIOME.HOT_DESERT : BIOME.COLD_DESERT;
  if (moist < steppeCut) return hot ? BIOME.SAVANNA : BIOME.STEPPE;
  if (moist < shrubCut) return temp > 0.45 ? BIOME.SHRUBLAND : BIOME.STEPPE;
  if (temp < 0.42) return BIOME.TAIGA;
  if (moist < 0.68) return hot ? BIOME.SEASONAL_FOREST : BIOME.TEMPERATE_FOREST;
  return hot ? BIOME.RAINFOREST : BIOME.WET_TEMPERATE;
}

function smoothBiomeIds(ids, h, n, seaLevel, iterations) {
  const nh = h.length / n; // map height (n is the width)
  const out = new Int16Array(ids.length);
  const counts = new Int16Array(32); // хватит на небольшую палитру

  for (let it = 0; it < iterations; it++) {
    out.set(ids);

    for (let y = 1; y < nh - 1; y++) {
      for (let x = 0; x < n; x++) {          // every column: east-west wraps
        const i = idx(x, y, n);
        if (h[i] <= seaLevel) continue; // океан не трогаем

        counts.fill(0);

        // текущий + 8 соседей
        for (let k = 0; k < 8; k++) {
          const j = nbr(x, y, N8[k][0], N8[k][1], n, nh);
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

    // The world is twice as wide as it is high (like a map of a planet): size = width
    const n = opts.size || 1280;
    const nh = n >> 1;
    ASPECT = n / nh;
    // the planet (prototype): Earth when not given
    const pl = opts.planet || {};
    PLANET = {
      tilt: clamp(Number(pl.tilt ?? EARTH.tilt), 0, 90),
      flux: clamp(Number(pl.flux ?? 1), 0.05, 20),
      cells: [1, 3, 5].includes(pl.cells) ? pl.cells : 3,
      locked: !!pl.locked,
      atmo: ATMO[pl.atmo] ? pl.atmo : 'earth',
      ecc: clamp(Number(pl.ecc ?? EARTH.ecc), 0, 0.9),
    };
    LAT_TEMP = isEarth(PLANET) ? null : PLANET.locked ? buildLockedTemp(PLANET) : buildLatTemp(PLANET);
    const seed = (opts.seed ?? 1337) | 0;
    const seaLevel = clamp(opts.seaLevel ?? 0.50, 0.05, 0.95);
    const riverThreshold = opts.riverThreshold ?? 600;

    const continentsMin = opts.continentsMin ?? 2;
    const continentsMax = opts.continentsMax ?? 6;

    const platesCount = clamp((opts.plates ?? 16) | 0, 4, 40);
    const windDir = (opts.windDir ?? 'W'); // 'W','E','N','S'

    canvas.width = n;
    canvas.height = nh;
    const ctx = canvas.getContext('2d', { alpha: false });

    // 1) Base height from continents
    let h = buildHeightmap(n, nh, seed, continentsMin, continentsMax);

    // 2) Land relief: lowlands, plateaus and hills instead of one dome per continent
    shapeLand(h, n, seed, seaLevel);
    // 2b) Volcanic island chains in the open ocean
    addIslandChains(h, n, seed, seaLevel);

    // 3) Mountains, rifts and island arcs along plate edges
    const plates = generatePlates(seed, platesCount);
    const plateFields = computePlateFields(n, nh, plates, seed);
    assignPlateTypes(plates, plateFields, h, seaLevel);
    const uplift = buildUpliftMap(n, plates, plateFields, seed);
	const mountainsStrength = clamp(opts.mountainsStrength ?? 0.18, 0, 0.6);

	for (let i = 0; i < h.length; i++) {
	  // full strength on land, weaker at sea (island arcs, mid-ocean ridges)
	  const m = clamp((h[i] - seaLevel + 0.04) / 0.08, 0, 1);
	  const factor = MOUNTAINS_AT_SEA + (1 - MOUNTAINS_AT_SEA) * m * m * (3 - 2 * m);
	  h[i] = clamp(h[i] + uplift[i] * mountainsStrength * MOUNTAINS_SCALE * factor, 0, 1);
	}

    // 4) Erosion (light) after tectonics
    streamErodeHalf(h, n, seaLevel, STREAM_ERODE_PASSES, STREAM_ERODE_K);
    thermalErode(h, n, 7, 0.017);
    shapeCoasts(h, n, seed, seaLevel);
    // 4b) Inland basins below sea level (not joined to the ocean) were all drawn as lakes:
    // up to 8% of the land under water. They become ordinary hollows just above the sea,
    // so only the deepest and largest ones fill with water (the lake budget in drawMap).
    {
      const open = connectedToEdge(h, n, seaLevel);
      for (let i = 0; i < h.length; i++) {
        if (h[i] <= seaLevel && !open[i]) {
          // the old shape squeezed to a gentle slope (exact flats make the water flow slow
          // to compute), high enough that rivers do not cut it below the sea
          h[i] = seaLevel + Math.max(0.006, 0.012 + (h[i] - seaLevel) * 0.25);
        }
      }
    }

    // 5) Climate (before hydrology: river flow depends on rain)
    const temp = computeTemperature(h, n, seaLevel, seed);
	const desertStrength = clamp(opts.desertStrength ?? 0.35, 0, 1);
	const windGrid = buildWindField(Math.max(64, Math.min(128, (nh / 6) | 0)), seed);

    // Absolute 0..1 index, no per-map min/max stretch (it used to amplify tiny differences)
    const rain = computeRain(h, n, seaLevel, temp, windGrid, desertStrength);

    // 6) Hydrology: fill depressions (with a tiny slope), D8 flow, rain-weighted accumulation
    const { filled: hf, order: flowOrder } = fillDepressions(h, n, seaLevel);
    const down = computeFlowDir(hf, n);
    const runoff = new Float32Array(n * nh);
    for (let i = 0; i < runoff.length; i++) {
      runoff[i] = h[i] > seaLevel ? RUNOFF_BASE + RUNOFF_RAIN * rain[i] : 0;
    }
    const acc = computeFlowAccum(down, flowOrder, runoff);

	debugBiomeStats(h, temp, rain, n, seaLevel, desertStrength);

    // 7) Draw	
	const mode = (opts.mode === 'height') ? 'height' : (opts.mode === 'wind' ? 'wind' : 'biomes');

	if (mode === 'wind') {
	  drawWindDebug(ctx, windGrid, n, nh);
	  return;
	}
	
	const biomeIds = new Int16Array(n * nh);
	const slopeScale = nh / 768;
	// wetlands: along rivers on flat ground (was flow > 2x threshold and slope < 0.004,
	// which almost never happened: 0.0% of land on 5 test seeds)
	const wetFlow = Math.max(150, riverThreshold);
	// a floodplain is a band a few px wide around the river (a 1 px strip was erased by
	// the zone smoothing below)
	const bigRiver = new Float32Array(n * nh);
	for (let i = 0; i < bigRiver.length; i++) if (h[i] > seaLevel && acc[i] > wetFlow) bigRiver[i] = 1;
	const nearRiver = blurred(bigRiver, n, Math.max(2, Math.round(3 * slopeScale)));
	// regional average height (sea counts as sea level) for the mountain moisture bonus
	const hFloor = new Float32Array(n * nh);
	for (let i = 0; i < hFloor.length; i++) hFloor[i] = Math.max(h[i], seaLevel);
	const hRegion = blurred(hFloor, n, Math.max(4, Math.round(32 * slopeScale)));
	const landRange = 1 - seaLevel + 1e-9;
	for (let y = 0; y < nh; y++) {
	  for (let x = 0; x < n; x++) {
	    const i = idx(x, y, n);
	    const xm = wrapX(x - 1, n), xp = wrapX(x + 1, n); // east-west wraps
	    const ym = Math.max(0, y - 1), yp = Math.min(nh - 1, y + 1);
	    const slope = Math.hypot(h[idx(xp, y, n)] - h[idx(xm, y, n)], h[idx(x, yp, n)] - h[idx(x, ym, n)]) * 0.5 * slopeScale;
	    const wet = nearRiver[i] > 0.05 && slope < 0.012;
	    const rise = (h[i] - hRegion[i]) / landRange;
	    biomeIds[i] = biomeId(h[i], seaLevel, temp[i], rain[i], desertStrength, slope, wet, rise);
	  }
	}
	// savanna / steppe / shrubs between deserts and forests, tundra between ice and forest
	addTransitions(biomeIds, h, n, seaLevel, Math.max(4, Math.round(11 * slopeScale)));
	smoothBiomeIds(biomeIds, h, n, seaLevel, 2);

	window.TerrainMap.lastMode = mode;
	window.TerrainMap.lastBiomeIds = biomeIds; // for the hover tooltip (0 = sea)
	window.TerrainMap.lastSize = n;
	window.TerrainMap.lastHeight = nh;
	window.TerrainMap.lastHeights = h;
	window.TerrainMap.lastSeaLevel = seaLevel;
	// for the seasons layer (seasonLayer): the year-mean temperature and the planet
	window.TerrainMap.lastTemp = temp;
	window.TerrainMap.lastPlanet = { ...PLANET };

drawMap(ctx, h, acc, temp, rain, biomeIds, n, seaLevel, riverThreshold, mode, desertStrength, hf, seed, down);
	
  }

  // Colours and the last mode, for the legend (map-legend.js)
  window.TerrainMap = {
    render, biomeColors: ATLAS_BIOMES, atlas: ATLAS,
    lastMode: null, lastBiomeIds: null, lastLake: null, lastSize: 0,
    lastHeights: null, lastSeaLevel: 0,
    // the planet module (planet-system.js) uses these for its explanations
    annualInsolation, latTempFor: (p) => buildLatTemp({ ...EARTH, ...p }),
    seasonLayer, lastIce: null, seaIceColor: SEA_ICE,
  };
})();
