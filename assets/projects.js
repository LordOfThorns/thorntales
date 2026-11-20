/* projects.js
   Floating soap-bubble field for projects.html.
   Real projects are links; small decorative bubbles can be popped for fun and respawn.
*/
(function () {
  'use strict';

  const tr = (key, fallback) => (window.I18n ? I18n.t(key) : fallback);

  // To add a project, add an entry here (texts: keys in i18n.js).
  const PROJECTS = [
    {
      label: tr('nav.mapGen', 'Генератор карт'),
      desc: tr('projects.mapGen.desc', 'странные карты несуществующих миров'),
      href: '/terrain.html',
    },
    {
      label: tr('nav.randomGen', 'Генератор всего'),
      desc: tr('projects.randomGen.desc', 'случайные вещи и существа'),
      href: '/random.html',
    },
  ];

  const DECOR_COUNT = 11;
  const POP_DELAY_MS = 420;

  const field = document.getElementById('bubbles');
  if (!field) return;

  const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');
  const bubbles = [];
  let W = 0;
  let H = 0;

  // Bubbles fly above the grass: their lowest point is a bit below the grass tips
  let floorY = 0;

  function measure() {
    W = field.clientWidth;
    H = field.clientHeight;
    grass.resize();
    floorY = H - grass.height * 0.55;
  }

  function projectSize() {
    return Math.max(150, Math.min(210, W * 0.28));
  }

  function rand(a, b) { return a + Math.random() * (b - a); }

  // ---------- Creation ----------

  function createProjectBubble(p, i) {
    const el = document.createElement('a');
    el.className = 'bubble bubble--project';
    el.href = p.href;
    el.innerHTML = `
      <span class="bubble__label">${escapeHtml(p.label)}</span>
      ${p.desc ? `<span class="bubble__desc">${escapeHtml(p.desc)}</span>` : ''}`;
    field.appendChild(el);

    const b = makeBody(el, projectSize(), { project: true });
    // Start projects roughly in the middle so they are visible right away
    b.x = W / 2 + (i - (PROJECTS.length - 1) / 2) * (b.r * 2.4);
    b.y = floorY / 2;
    bubbles.push(b);

    el.addEventListener('click', (e) => {
      if (e.button !== 0 || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) return;
      if (reducedMotion.matches) return;
      e.preventDefault();
      pop(b);
      setTimeout(() => { window.location.href = p.href; }, POP_DELAY_MS);
    });
  }

  function createDecorBubble() {
    const el = document.createElement('span');
    el.className = 'bubble bubble--decor';
    el.setAttribute('aria-hidden', 'true');
    field.appendChild(el);

    const b = makeBody(el, rand(18, 64), { project: false });
    placeFree(b);
    bubbles.push(b);

    el.addEventListener('pointerdown', () => {
      if (b.popping) return;
      pop(b);
      setTimeout(() => {
        // Respawn somewhere else
        el.classList.remove('is-popping');
        b.popping = false;
        b.r = rand(9, 32);
        applySize(b);
        placeFree(b);
        el.classList.add('is-appearing');
        setTimeout(() => el.classList.remove('is-appearing'), 600);
      }, POP_DELAY_MS + 900);
    });
  }

  function makeBody(el, size, opts) {
    const b = {
      el,
      r: size / 2,
      x: 0,
      y: 0,
      vx: rand(-1, 1),
      vy: rand(-1, 1),
      phase: rand(0, Math.PI * 2),   // for the wobble
      hue: rand(0, 360),             // iridescence offset
      project: opts.project,
      hovered: false,
      popping: false,
    };
    // Base speed: projects drift slower than small bubbles
    const speed = opts.project ? 0.22 : rand(0.25, 0.6);
    const len = Math.hypot(b.vx, b.vy) || 1;
    b.vx = (b.vx / len) * speed;
    b.vy = (b.vy / len) * speed;
    b.speed = speed;

    el.style.setProperty('--hue', `${b.hue.toFixed(0)}deg`);
    el.style.setProperty('--tw', `${(-rand(0, 3.4)).toFixed(2)}s`); // sparkles twinkle out of step
    applySize(b);

    el.addEventListener('pointerenter', () => { b.hovered = true; });
    el.addEventListener('pointerleave', () => { b.hovered = false; });
    el.addEventListener('focus', () => { b.hovered = true; });
    el.addEventListener('blur', () => { b.hovered = false; });
    return b;
  }

  function applySize(b) {
    b.el.style.width = `${b.r * 2}px`;
    b.el.style.height = `${b.r * 2}px`;
  }

  function placeFree(b) {
    // Try a few random spots that don't overlap other bubbles
    for (let attempt = 0; attempt < 40; attempt++) {
      b.x = rand(b.r, Math.max(b.r, W - b.r));
      b.y = rand(b.r, Math.max(b.r, floorY - b.r));
      const clash = bubbles.some((o) => o !== b && Math.hypot(o.x - b.x, o.y - b.y) < o.r + b.r + 6);
      if (!clash) return;
    }
  }

  function pop(b) {
    b.popping = true;
    b.el.classList.add('is-popping');
    burst(b);
  }

  // The film vanishes at once; what is left is a thin ring and droplets flying off its rim
  function burst(b) {
    const el = document.createElement('span');
    el.className = 'bubble-burst';
    el.setAttribute('aria-hidden', 'true');
    el.style.transform = `translate(${b.x.toFixed(1)}px, ${b.y.toFixed(1)}px)`;
    el.style.setProperty('--r', `${b.r.toFixed(1)}px`);
    const ring = document.createElement('span');
    ring.className = 'bubble-burst__ring';
    el.append(ring);
    const drops = Math.max(6, Math.min(16, Math.round(b.r / 3)));
    for (let i = 0; i < drops; i++) {
      const a = (i / drops) * Math.PI * 2 + rand(-0.3, 0.3);
      const far = b.r * rand(1.25, 1.85);
      const d = document.createElement('span');
      d.className = 'bubble-burst__drop';
      d.style.setProperty('--sx', `${(Math.cos(a) * b.r).toFixed(1)}px`);
      d.style.setProperty('--sy', `${(Math.sin(a) * b.r).toFixed(1)}px`);
      d.style.setProperty('--dx', `${(Math.cos(a) * far).toFixed(1)}px`);
      d.style.setProperty('--dy', `${(Math.sin(a) * far + rand(4, 14)).toFixed(1)}px`); // droplets fall a little
      d.style.setProperty('--dur', `${Math.round(rand(320, 520))}ms`);
      el.append(d);
    }
    field.append(el);
    setTimeout(() => el.remove(), 650);
  }

  // ---------- Night: sky, grass, fireflies ----------
  // Canvases inside the field, bottom to top: the sky with twinkling stars, the back grass
  // (behind the bubbles), then the front grass with fireflies (in front of them).
  // The wind is a slow wave plus gusts running across the field; blades under a low bubble
  // lean away from it, and the cursor ruffles them: they swing with its motion and spring back.

  const grass = (() => {
    const make = (cls) => {
      const cv = document.createElement('canvas');
      cv.className = `bubble-field__night bubble-field__night--${cls}`;
      cv.setAttribute('aria-hidden', 'true');
      field.append(cv);
      return cv;
    };
    const sky = { cv: make('sky') };
    const layers = [
      { cv: make('back'), step: 3.2, hMin: 0.45, hMax: 1, wMin: 2.2, wMax: 4,
        tops: ['rgba(70,125,150,.55)', 'rgba(100,155,175,.6)'], foot: 'rgba(5,10,20,.95)' },
      { cv: make('front'), step: 7, hMin: 0.22, hMax: 0.6, wMin: 2.6, wMax: 4.6,
        tops: ['rgba(110,165,185,.75)', 'rgba(140,195,210,.8)'], foot: 'rgba(3,6,12,1)' },
    ];
    const g = { height: 0 };
    let gw = 0;
    let gh = 0;      // field height
    let ch = 0;      // grass canvases are taller than the grass: fireflies fly above it
    let stars = [];
    let flies = [];

    // the cursor over the field (field coordinates); vx = its last horizontal move
    const ptr = { x: -1e4, y: -1e4, vx: 0, at: 0, inside: false };
    field.addEventListener('pointermove', (e) => {
      const r = field.getBoundingClientRect();
      const x = e.clientX - r.left;
      ptr.vx = ptr.inside ? x - ptr.x : 0;
      ptr.x = x;
      ptr.y = e.clientY - r.top;
      ptr.at = performance.now();
      ptr.inside = true;
    });
    field.addEventListener('pointerleave', () => { ptr.inside = false; ptr.x = -1e4; });

    function setup(cv, w, h) {
      const dpr = Math.min(2, window.devicePixelRatio || 1);
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
      cv.style.height = `${h}px`;
      const ctx = cv.getContext('2d');
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      return ctx;
    }

    g.resize = () => {
      gw = field.clientWidth;
      gh = field.clientHeight;
      g.height = Math.round(Math.max(70, Math.min(130, gh * 0.2)));
      ch = Math.round(g.height * 2.2);

      sky.ctx = setup(sky.cv, gw, gh);
      stars = [];
      const n = Math.round((gw * gh) / 3600);
      for (let i = 0; i < n; i++) {
        const big = Math.random() < 0.08;
        stars.push({
          x: rand(0, gw),
          y: rand(0, gh - g.height * 0.6),
          r: big ? rand(1.2, 1.8) : rand(0.5, 1.1),
          a: big ? rand(0.7, 1) : rand(0.35, 0.85),
          speed: rand(0.4, 1.6),
          phase: rand(0, Math.PI * 2),
          blue: Math.random() < 0.5,
          big,
        });
      }

      for (const L of layers) {
        L.ctx = setup(L.cv, gw, ch);
        L.fills = L.tops.map((top) => {
          const gr = L.ctx.createLinearGradient(0, ch, 0, ch - g.height);
          gr.addColorStop(0, L.foot);
          gr.addColorStop(1, top);
          return gr;
        });
        L.blades = [];
        for (let x = -6; x < gw + 6; x += L.step * rand(0.6, 1.4)) {
          L.blades.push({
            x,
            h: g.height * rand(L.hMin, L.hMax),
            w: rand(L.wMin, L.wMax),
            flex: rand(0.7, 1.3),       // some blades give in to the wind more
            lean: rand(-0.12, 0.12),    // own slant at rest
            phase: rand(0, Math.PI * 2),
            shade: Math.random() < 0.35 ? 1 : 0,
            ro: 0,                      // ruffle: extra bend from the cursor, springs back
            rv: 0,
          });
        }
      }

      // a few fireflies here and there over the grass
      const count = Math.max(3, Math.min(11, Math.round(gw / 130)));
      flies = Array.from({ length: count }, () => ({
        x: rand(0, gw),
        y: rand(ch - g.height * 2, ch - g.height * 0.6),
        dir: rand(0, Math.PI * 2),
        speed: rand(0.15, 0.35),
        phase: rand(0, Math.PI * 2),
        blink: rand(0.5, 1.1),
      }));
    };

    function wind(x, t) {
      const span = gw + 800;
      const gustX = ((t * 110) % span) - 400; // a gust runs left to right every ~15 s
      const gust = Math.exp(-(((x - gustX) / 170) ** 2)) * 0.16;
      return 0.07 + 0.05 * Math.sin(t * 0.7 + x * 0.006) + 0.03 * Math.sin(t * 1.7 + x * 0.021) + gust;
    }

    function drawSky(t) {
      const ctx = sky.ctx;
      ctx.clearRect(0, 0, gw, gh);
      for (const s of stars) {
        const tw = 0.55 + 0.45 * Math.sin(t * s.speed + s.phase);
        ctx.globalAlpha = s.a * tw;
        ctx.fillStyle = s.blue ? '#cfe4ff' : '#f2f0ff';
        ctx.beginPath();
        ctx.arc(s.x, s.y, s.r, 0, Math.PI * 2);
        ctx.fill();
        if (s.big) { // a faint cross of light on the brightest ones
          ctx.globalAlpha = s.a * tw * 0.35;
          ctx.fillRect(s.x - s.r * 4, s.y - 0.3, s.r * 8, 0.6);
          ctx.fillRect(s.x - 0.3, s.y - s.r * 4, 0.6, s.r * 8);
        }
      }
      ctx.globalAlpha = 1;
    }

    function drawGrass(L, t, low, ruffle) {
      const ctx = L.ctx;
      const base = ch + 1;
      ctx.clearRect(0, 0, gw, ch);
      for (let s = 0; s < 2; s++) {
        ctx.beginPath();
        for (const bl of L.blades) {
          if (s === 0) {
            // the cursor: a move across the blade tips swings them, standing in the grass parts it
            if (ptr.inside) {
              const dx = bl.x - ptr.x;
              const reach = 46;
              const depth = ptr.y - (gh - bl.h); // > 0: the cursor is below this blade's tip
              if (Math.abs(dx) < reach && depth > -22) {
                const k = (1 - Math.abs(dx) / reach) * Math.min(1, (depth + 22) / 30);
                if (ruffle) bl.rv += Math.max(-30, Math.min(30, ptr.vx)) * 0.0045 * k;
                bl.rv += Math.sign(dx || 1) * 0.012 * k;
              }
            }
            bl.rv += -bl.ro * 0.07;   // spring back
            bl.rv *= 0.88;            // and calm down
            bl.ro += bl.rv;
          }
          if (bl.shade !== s) continue;
          let bend = bl.lean + wind(bl.x, t) * bl.flex + Math.sin(t * 2.6 + bl.phase) * 0.02 + bl.ro;
          for (const b of low) {
            const dx = bl.x - b.x;
            const reach = b.r * 1.15;
            if (Math.abs(dx) > reach) continue;
            const dy = (b.y + b.r) - (gh - bl.h); // how deep the bubble sinks into this blade
            if (dy < -14) continue;
            const k = Math.min(1, (dy + 14) / 34) * (1 - Math.abs(dx) / reach);
            bend += Math.sign(dx || 1) * 0.45 * k;
          }
          bend = Math.max(-1.1, Math.min(1.1, bend));
          const h = bl.h * (1 - Math.min(0.4, Math.abs(bend) * 0.32)); // a bent blade is lower
          const tipX = bl.x + bend * bl.h;
          const tipY = base - h;
          const midX = bl.x + bend * bl.h * 0.35;
          const midY = base - h * 0.55;
          ctx.moveTo(bl.x - bl.w / 2, base);
          ctx.quadraticCurveTo(midX - bl.w * 0.2, midY, tipX, tipY);
          ctx.quadraticCurveTo(midX + bl.w * 0.2, midY, bl.x + bl.w / 2, base);
          ctx.closePath();
        }
        ctx.fillStyle = L.fills[s];
        ctx.fill();
      }
    }

    function drawFlies(ctx, t, moving) {
      ctx.globalCompositeOperation = 'lighter';
      for (const f of flies) {
        if (moving) {
          f.dir += Math.sin(t * 0.9 + f.phase) * 0.03 + rand(-0.05, 0.05);
          f.x += Math.cos(f.dir) * f.speed;
          f.y += Math.sin(f.dir) * f.speed * 0.6;
          // keep to the band just over the grass, and inside the field
          const top = ch - g.height * 2.05;
          const bottom = ch - g.height * 0.55;
          if (f.y < top) { f.y = top; f.dir = Math.abs(f.dir) % Math.PI; }
          if (f.y > bottom) { f.y = bottom; f.dir = -Math.abs(f.dir) % Math.PI; }
          if (f.x < -10) f.x = gw + 10;
          if (f.x > gw + 10) f.x = -10;
        }
        // a slow glow, then dark for a while
        const glow = Math.max(0, Math.sin(t * f.blink + f.phase));
        const a = glow * glow;
        if (a < 0.02) continue;
        const r = 9 + 5 * a;
        const gr = ctx.createRadialGradient(f.x, f.y, 0, f.x, f.y, r);
        gr.addColorStop(0, `rgba(235,255,190,${0.95 * a})`);
        gr.addColorStop(0.18, `rgba(200,250,150,${0.55 * a})`);
        gr.addColorStop(1, 'rgba(160,230,120,0)');
        ctx.fillStyle = gr;
        ctx.fillRect(f.x - r, f.y - r, r * 2, r * 2);
      }
      ctx.globalCompositeOperation = 'source-over';
    }

    g.draw = (t, moving) => {
      const ruffle = performance.now() - ptr.at < 80; // the cursor has just moved
      // only bubbles low enough to touch the grass can push it
      const low = bubbles.filter((b) => !b.popping && b.y + b.r > gh - g.height - 14);
      drawSky(t);
      drawGrass(layers[0], t, low, ruffle);
      drawGrass(layers[1], t, low, ruffle);
      drawFlies(layers[1].ctx, t, moving);
    };
    return g;
  })();

  // ---------- Simulation ----------

  function step(dt) {
    const t = performance.now() / 1000;

    for (const b of bubbles) {
      if (b.popping) continue;
      // Hovered bubbles almost stop so they are easy to click
      const k = b.hovered ? 0.08 : 1;
      // Gentle wobble so paths are not straight lines
      const wobX = Math.cos(t * 0.6 + b.phase) * 0.06;
      const wobY = Math.sin(t * 0.8 + b.phase) * 0.08;
      b.x += (b.vx + wobX) * k * dt;
      b.y += (b.vy + wobY) * k * dt;

      // Walls
      if (b.x < b.r) { b.x = b.r; b.vx = Math.abs(b.vx); }
      if (b.x > W - b.r) { b.x = W - b.r; b.vx = -Math.abs(b.vx); }
      if (b.y < b.r) { b.y = b.r; b.vy = Math.abs(b.vy); }
      if (b.y > floorY - b.r) { b.y = floorY - b.r; b.vy = -Math.abs(b.vy); }
    }

    // Soft collisions between bubbles
    for (let i = 0; i < bubbles.length; i++) {
      const a = bubbles[i];
      if (a.popping) continue;
      for (let j = i + 1; j < bubbles.length; j++) {
        const c = bubbles[j];
        if (c.popping) continue;
        const dx = c.x - a.x;
        const dy = c.y - a.y;
        const dist = Math.hypot(dx, dy) || 0.001;
        const overlap = a.r + c.r - dist;
        if (overlap <= 0) continue;

        const nx = dx / dist;
        const ny = dy / dist;
        // Heavier (bigger) bubbles move less
        const ma = a.r * a.r;
        const mc = c.r * c.r;
        const total = ma + mc;
        a.x -= nx * overlap * (mc / total);
        a.y -= ny * overlap * (mc / total);
        c.x += nx * overlap * (ma / total);
        c.y += ny * overlap * (ma / total);

        // Exchange velocity along the normal (elastic-ish)
        const rel = (c.vx - a.vx) * nx + (c.vy - a.vy) * ny;
        if (rel < 0) {
          const imp = (2 * rel) / total;
          a.vx += imp * mc * nx; a.vy += imp * mc * ny;
          c.vx -= imp * ma * nx; c.vy -= imp * ma * ny;
        }
      }
    }

    // Keep each bubble near its own cruising speed (no slow-down or runaway after collisions)
    for (const b of bubbles) {
      const v = Math.hypot(b.vx, b.vy) || 0.001;
      const target = b.speed;
      const f = 1 + (target / v - 1) * 0.02;
      b.vx *= f;
      b.vy *= f;
    }
  }

  function render() {
    for (const b of bubbles) {
      b.el.style.transform = `translate(${(b.x - b.r).toFixed(1)}px, ${(b.y - b.r).toFixed(1)}px)`;
    }
  }

  let last = performance.now();
  function frame(now) {
    // dt in "60fps frames", capped so a background tab doesn't teleport bubbles
    const dt = Math.min(3, (now - last) / 16.67);
    last = now;
    if (!reducedMotion.matches) step(dt);
    render();
    grass.draw(reducedMotion.matches ? 0 : now / 1000, !reducedMotion.matches);
    requestAnimationFrame(frame);
  }

  // ---------- Start ----------

  measure();
  PROJECTS.forEach(createProjectBubble);
  for (let i = 0; i < DECOR_COUNT; i++) createDecorBubble();
  render();
  field.classList.add('is-ready');
  requestAnimationFrame(frame);

  window.addEventListener('resize', () => {
    measure();
    for (const b of bubbles) {
      if (b.project) {
        b.r = projectSize() / 2;
        applySize(b);
      }
      b.x = Math.min(Math.max(b.x, b.r), W - b.r);
      b.y = Math.min(Math.max(b.y, b.r), floorY - b.r);
    }
  });

  // Back/forward cache: un-pop the project bubble when coming back
  window.addEventListener('pageshow', (e) => {
    if (!e.persisted) return;
    for (const b of bubbles) {
      if (b.project && b.popping) {
        b.popping = false;
        b.el.classList.remove('is-popping');
      }
    }
  });

  function escapeHtml(s) {
    return String(s ?? '').replace(/[&<>"']/g, (ch) => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
    }[ch]));
  }
})();
