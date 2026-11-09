/* projects.js
   Floating soap-bubble field for projects.html.
   Real projects are links; small decorative bubbles can be popped for fun and respawn.
*/
(function () {
  'use strict';

  // To add a project, add an entry here.
  const PROJECTS = [
    {
      label: 'Генератор карт',
      desc: 'странные карты несуществующих миров',
      href: '/terrain.html',
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

  function measure() {
    W = field.clientWidth;
    H = field.clientHeight;
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
    b.y = H / 2;
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
      b.y = rand(b.r, Math.max(b.r, H - b.r));
      const clash = bubbles.some((o) => o !== b && Math.hypot(o.x - b.x, o.y - b.y) < o.r + b.r + 6);
      if (!clash) return;
    }
  }

  function pop(b) {
    b.popping = true;
    b.el.classList.add('is-popping');
  }

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
      if (b.y > H - b.r) { b.y = H - b.r; b.vy = -Math.abs(b.vy); }
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
      b.y = Math.min(Math.max(b.y, b.r), H - b.r);
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
