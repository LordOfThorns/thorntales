/* main.js — stained-glass windows on the home page: click breaks the glass, then
   navigates. The break effect:
     - impact flash at the hit point;
     - shards cut from the window's own picture (each shows the piece it came from),
       tumbling in 3D, glinting as they turn, flying out and falling in an arc;
     - a few twinkling sparkles.

   Smoothness: only transform and opacity are animated (the GPU handles both without
   repainting), the effect lives in its own layer on top of the page (inside the window
   it would sit under the window's drop-shadow filter, which then has to be redrawn
   every frame), the glint is a white layer whose opacity changes (not a brightness
   filter), and the flight ends before the browser leaves the page. The target page is
   prefetched when the pointer comes over a window.
*/
(function () {
  "use strict";

  const BREAK_DELAY_MS = 680;  // navigation starts after this; flights end before it
  const SHARDS = 30;
  const SPARKLES = 12;
  const GLASS_INSET = 12;      // .glass is inset by 12px inside .window (style.css)
  const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");

  const rand = (a, b) => a + Math.random() * (b - a);
  const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);

  // ---------- Frame marks: tiny embossed curls and words in an unknown script ----------
  //
  // A few marks per window, placed at random spots of the frame band (the 12px between
  // the outer edge and the glass): along the straight sides and the bottom (at the apex
  // the band is too thin). Rusty gold, pressed into the frame (dark edge above, faint light below, as in a groove).
  // Seeded by the window index, so they stay in the same places on every visit.
  const FRAME = 12;                  // frame width in px, as .glass inset in style.css
  const BOX_W = 250, BOX_H = 250 * 5 / 3;

  // glyphs of the unknown script: rounded loops, hooks and waves in a 4 x 6 cell
  const GLYPHS = [
    "M0,5 C0,1 4,1 4,5",                       // arch
    "M0,1 C0,5 4,5 4,1",                       // cup
    "M0,3 C0,0 4,0 4,3 C4,6 0,6 0,3",          // loop
    "M0,6 C0,2 2,0 3,1 C4,2 2,4 0,3",          // curl
    "M0,3 C1,0 3,6 4,3",                       // wave
    "M3,0 C0,0 0,3 2,3 C4,3 4,6 1,6",          // s
    "M0,6 C1,3 1,1 3,0 C4.5,0 4,2 3,2",        // hook
    "M1,0 C-0.5,2 -0.5,4 1,6 M3,1 C4.3,2.5 4.3,3.5 3,5", // ( )
    "M0,4 C1,2 3,2 4,4 C3,6 1,6 0,4 M2,0.4 C2.4,0.4 2.4,1.2 2,1.2", // eye with a dot
    "M0,2 C2,-1 4,2 2,3 C0,4 2,7 4,5",         // squiggle
    "M0,5 C1,5 1,1 2,1 C3,1 3,5 4,5",          // hump
    "M4,1 C1,0 0,2 2,3 C3.5,3.8 3,6 0.5,5.5",  // reversed s
  ];
  // scrolls and spirals, about 10 x 8, centred on 0,0
  const CURLS = [
    "M-5,2 C-5,-3 1,-4 2,0 C2.6,2.4 -1,3 -1.2,0.6 C-1.3,-0.8 0.4,-1 0.6,0",
    "M-5,0 C-3,-4 2,4 5,0 M-5,0 C-6.4,-1 -6,1.6 -4.6,1.2 M5,0 C6.4,1 6,-1.6 4.6,-1.2",
    "M0,4 C-4,2 -4,-4 0,-4 C4,-4 4,1.5 0.5,1.5 C-1.5,1.5 -1.5,-1.5 0.4,-1.5",
    "M-5,3 C-1,3 -2,-4 1.5,-4 C4.5,-4 4.5,0 2,0 C0.8,0 0.8,-1.6 2,-1.6",
    "M-5,0 C-4,-3 -1,-3 0,0 C1,3 4,3 5,0 C4,-2 2.5,-1 3.2,0.4",
  ];

  function frameMarks(index) {
    let s = (index + 3) * 7919;
    const rnd = () => ((s = (s * 16807) % 2147483647) / 2147483647);
    const STEP = 5.4; // glyph cell plus spacing
    const word = () => {
      // mostly middle length, now and then a long inscription
      const n = 4 + Math.floor(Math.pow(rnd(), 1.4) * 13);
      let d = "";
      for (let i = 0; i < n; i++) {
        const g = GLYPHS[Math.floor(rnd() * GLYPHS.length)];
        // shift every glyph to its cell
        d += g.replace(/(-?\d+(?:\.\d+)?),(-?\d+(?:\.\d+)?)/g, (_, x, y) => `${(+x + i * STEP).toFixed(1)},${y}`) + " ";
      }
      const len = n * STEP - 1.4;
      // a flourish under some words: a slow wave ending in a little curl
      if (rnd() < 0.45) {
        const a = len / 3, b = (2 * len) / 3;
        d += `M-1,8 C${a.toFixed(1)},10 ${b.toFixed(1)},6 ${(len + 1).toFixed(1)},8 ` +
             `C${(len + 3).toFixed(1)},9 ${(len + 3).toFixed(1)},6.4 ${(len + 1.4).toFixed(1)},6.8`;
      }
      return { d, len };
    };

    const marks = [];
    const count = 3 + Math.floor(rnd() * 3);
    const spots = ["left", "right", "bottom", "left", "right"];
    for (let k = 0; k < count; k++) {
      const spot = spots[Math.floor(rnd() * spots.length)];
      const isCurl = rnd() < 0.3;
      const m = isCurl
        ? { d: CURLS[Math.floor(rnd() * CURLS.length)], len: 10, scale: 0.7 }
        : { ...word(), scale: 0.6 };
      const half = (m.len * m.scale) / 2;
      let x, y, rot;
      if (spot === "bottom") { x = BOX_W * (0.2 + rnd() * 0.6); y = BOX_H - FRAME / 2; rot = 0; }
      else {
        // the straight part of the sides, below the shoulders
        x = spot === "left" ? FRAME / 2 : BOX_W - FRAME / 2;
        const top = BOX_H * 0.36 + half, bottom = BOX_H - FRAME - half; // keep the whole mark on the side
        y = top + rnd() * Math.max(0, bottom - top);
        rot = isCurl ? 0 : (spot === "left" ? -90 : 90);
      }
      const tx = isCurl ? 0 : -half / m.scale, ty = isCurl ? 0 : -4; // centre a word on its spot
      marks.push(`<g transform="translate(${x.toFixed(1)} ${y.toFixed(1)}) rotate(${rot}) scale(${m.scale}) translate(${tx.toFixed(1)} ${ty})"><path d="${m.d}"/></g>`);
    }
    const body = marks.join("");
    return `
      <svg class="frame-marks" viewBox="0 0 ${BOX_W} ${BOX_H.toFixed(2)}" preserveAspectRatio="none" aria-hidden="true" focusable="false">
        <g class="frame-marks__shade" transform="translate(-.4 -.5)">${body}</g>
        <g class="frame-marks__light" transform="translate(.4 .5)">${body}</g>
        <g class="frame-marks__gold">${body}</g>
      </svg>`;
  }

  // One overlay per window, placed over it in page coordinates
  function effectLayer(winEl) {
    let layer = winEl._shardLayer;
    if (!layer) {
      layer = document.createElement("div");
      layer.className = "shards";
      document.body.appendChild(layer);
      winEl._shardLayer = layer;
    }
    const r = winEl.getBoundingClientRect();
    layer.style.left = `${r.left + window.scrollX}px`;
    layer.style.top = `${r.top + window.scrollY}px`;
    layer.style.width = `${r.width}px`;
    layer.style.height = `${r.height}px`;
    layer.innerHTML = "";
    return layer;
  }

  // Picture of a window and its natural size, so shards can show the right piece.
  // The CSS --img url is resolved by the browser on the .img element.
  const pictureCache = new Map();
  function pictureOf(winEl) {
    const imgEl = winEl.querySelector(".img");
    if (!imgEl) return null;
    const bg = getComputedStyle(imgEl).backgroundImage;
    const m = bg && bg.match(/url\(["']?([^"')]+)["']?\)/);
    if (!m) return null;
    const url = m[1];
    if (!pictureCache.has(url)) {
      const im = new Image();
      const entry = { url, w: 0, h: 0 };
      im.onload = () => { entry.w = im.naturalWidth; entry.h = im.naturalHeight; };
      im.src = url;
      pictureCache.set(url, entry);
    }
    return pictureCache.get(url);
  }

  // Random shard outline: a jagged triangle or quad
  function shardShape() {
    const pts = Math.random() < 0.55
      ? [[rand(0, 30), rand(0, 20)], [rand(70, 100), rand(0, 45)], [rand(30, 80), rand(80, 100)]]
      : [[rand(0, 25), rand(0, 25)], [rand(60, 100), rand(0, 20)], [rand(75, 100), rand(60, 100)], [rand(0, 35), rand(70, 100)]];
    return `polygon(${pts.map(([x, y]) => `${x.toFixed(0)}% ${y.toFixed(0)}%`).join(", ")})`;
  }

  function spawnBreak(winEl) {
    const layer = effectLayer(winEl);

    const W = winEl.clientWidth;
    const H = winEl.clientHeight;
    const cx = W / 2;
    const cy = H * 0.42; // impact a little above the centre looks better

    // background geometry of the picture as drawn in .glass (background-size: cover)
    const pic = pictureOf(winEl);
    const gw = W - GLASS_INSET * 2, gh = H - GLASS_INSET * 2;
    let bgW = gw, bgH = gh, bgX = 0, bgY = 0;
    if (pic && pic.w && pic.h) {
      const s = Math.max(gw / pic.w, gh / pic.h);
      bgW = pic.w * s; bgH = pic.h * s;
      bgX = (gw - bgW) / 2; bgY = (gh - bgH) / 2;
    }

    const frag = document.createDocumentFragment();
    const jobs = []; // [element, keyframes, options], started together after insertion

    // impact flash
    const flash = document.createElement("div");
    flash.className = "impact";
    flash.style.left = `${cx - 60}px`;
    flash.style.top = `${cy - 60}px`;
    frag.appendChild(flash);
    jobs.push([flash,
      [{ transform: "scale(0.2)", opacity: 0.95 }, { transform: "scale(1.6)", opacity: 0 }],
      { duration: 220, easing: "cubic-bezier(.2,.7,.3,1)", fill: "forwards" }]);

    // shards
    for (let i = 0; i < SHARDS; i++) {
      const shard = document.createElement("div");
      shard.className = "shard";

      // size: a few big pieces, lots of crumbs
      const size = 7 + Math.pow(Math.random(), 1.8) * 42;
      const w = size * rand(0.7, 1.3), h = size * rand(0.7, 1.3);
      const spread = 18 + Math.max(0, 48 - size) * 0.8;
      const x = clamp(cx + rand(-spread, spread) - w / 2, GLASS_INSET, W - GLASS_INSET - w);
      const y = clamp(cy + rand(-spread, spread) - h / 2, GLASS_INSET, H - GLASS_INSET - h);

      shard.style.left = `${x}px`;
      shard.style.top = `${y}px`;
      shard.style.width = `${w}px`;
      shard.style.height = `${h}px`;
      shard.style.clipPath = shardShape();
      if (pic) {
        // the piece of the picture that was under this shard, with a glassy sheen
        shard.style.backgroundImage =
          `linear-gradient(135deg, rgba(255,255,255,.6), rgba(255,255,255,.05) 35%, rgba(255,255,255,0) 55%, rgba(255,255,255,.25) 80%, rgba(255,255,255,.5)), url("${pic.url}")`;
        shard.style.backgroundSize = `100% 100%, ${bgW}px ${bgH}px`;
        shard.style.backgroundPosition = `0 0, ${bgX + GLASS_INSET - x}px ${bgY + GLASS_INSET - y}px`;
      }
      const glint = document.createElement("div");
      glint.className = "shard__glint";
      shard.appendChild(glint);
      frag.appendChild(shard);

      // flight: out from the impact, up, then falling; small pieces fly further
      const ang = Math.atan2(y + h / 2 - cy, x + w / 2 - cx) + rand(-0.5, 0.5);
      const power = rand(130, 260) * (1.35 - size / 60);
      const dx = Math.cos(ang) * power;
      const up = rand(50, 120) * (1.3 - size / 70);
      const fall = rand(120, 220);
      const ax = rand(-1, 1).toFixed(2), ay = rand(-1, 1).toFixed(2), az = rand(-0.3, 0.3).toFixed(2);
      const spin = rand(240, 620) * (Math.random() < 0.5 ? -1 : 1);
      const dur = rand(480, 640);
      const glintAt = rand(0.2, 0.55);

      jobs.push([shard,
        [
          { offset: 0, transform: "translate3d(0,0,0) rotate3d(0,0,1,0deg) scale(1)", opacity: 1, easing: "cubic-bezier(.15,.8,.4,1)" },
          { offset: 0.35, transform: `translate3d(${dx * 0.55}px, ${Math.sin(ang) * power * 0.45 - up}px, 0) rotate3d(${ax},${ay},${az},${spin * 0.35}deg) scale(1)`, opacity: 1, easing: "cubic-bezier(.5,0,.9,.55)" },
          { offset: 0.72, opacity: 1 },
          { offset: 1, transform: `translate3d(${dx}px, ${Math.sin(ang) * power * 0.6 - up + fall}px, 0) rotate3d(${ax},${ay},${az},${spin}deg) scale(0.85)`, opacity: 0 },
        ],
        { duration: dur, fill: "forwards" }]);
      // glint: the white layer flashes once while the shard tumbles
      jobs.push([glint,
        [
          { offset: 0, opacity: 0 },
          { offset: Math.max(0.05, glintAt - 0.08), opacity: 0 },
          { offset: glintAt, opacity: 0.85 },
          { offset: Math.min(0.95, glintAt + 0.1), opacity: 0 },
          { offset: 1, opacity: 0 },
        ],
        { duration: dur, fill: "forwards" }]);
    }

    // sparkles
    for (let i = 0; i < SPARKLES; i++) {
      const s = document.createElement("div");
      s.className = "sparkle";
      const size = rand(6, 13);
      s.style.width = s.style.height = `${size}px`;
      s.style.left = `${cx - size / 2 + rand(-14, 14)}px`;
      s.style.top = `${cy - size / 2 + rand(-14, 14)}px`;
      frag.appendChild(s);
      const ang = rand(0, Math.PI * 2);
      const dist = rand(40, 120);
      const p = (k) => `translate3d(${Math.cos(ang) * dist * k}px, ${Math.sin(ang) * dist * k}px, 0)`;
      jobs.push([s,
        [
          { transform: `${p(0)} scale(0) rotate(0deg)`, opacity: 0 },
          { offset: 0.25, transform: `${p(0.4)} scale(1.15) rotate(45deg)`, opacity: 1 },
          { offset: 0.5, transform: `${p(0.7)} scale(0.5) rotate(90deg)`, opacity: 0.6 },
          { offset: 0.7, transform: `${p(0.85)} scale(1) rotate(135deg)`, opacity: 1 },
          { transform: `translate3d(${Math.cos(ang) * dist}px, ${Math.sin(ang) * dist + 20}px, 0) scale(0) rotate(180deg)`, opacity: 0 },
        ],
        { duration: rand(400, 560), delay: rand(0, 120), easing: "ease-out", fill: "both" }]);
    }

    // insert everything at once, then start all animations in the same frame
    layer.appendChild(frag);
    for (const [el, kf, opt] of jobs) el.animate(kf, opt);

    // clean up after the animation so the DOM does not grow
    setTimeout(() => { layer.innerHTML = ""; }, 1000);
  }

  // Start downloading the target page early, so leaving does not stall the animation
  const prefetched = new Set();
  function prefetch(href) {
    if (!href || prefetched.has(href)) return;
    prefetched.add(href);
    const link = document.createElement("link");
    link.rel = "prefetch";
    link.href = href;
    document.head.appendChild(link);
  }

  function onReady(fn) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", fn);
    } else {
      fn();
    }
  }

  // When the page is restored from the back/forward cache, mend the broken glass
  window.addEventListener("pageshow", (e) => {
    if (!e.persisted) return;
    document.querySelectorAll(".window.is-broken").forEach((w) => {
      w.classList.remove("is-broken");
      if (w._shardLayer) w._shardLayer.innerHTML = "";
    });
  });

  onReady(() => {
    document.querySelectorAll(".window").forEach((w, i) => {
      w.querySelector(".lead")?.insertAdjacentHTML("afterend", frameMarks(i));
      pictureOf(w); // preload the picture size for the shards
      const href = w.getAttribute("href");
      w.addEventListener("pointerenter", () => prefetch(href));
      w.addEventListener("focus", () => prefetch(href));

      w.addEventListener("click", (e) => {
        // Let the browser handle "open in new tab/window" clicks and users who disabled animations
        if (e.button !== 0 || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) return;
        if (reducedMotion.matches) return;

        e.preventDefault();

        // already broken: go right away
        if (w.classList.contains("is-broken")) {
          if (href) window.location.href = href;
          return;
        }

        spawnBreak(w);
        w.classList.add("is-broken"); // hides the lead grid, zooms the picture inside

        if (href) {
          setTimeout(() => { window.location.href = href; }, BREAK_DELAY_MS);
        }
      });
    });
  });

  // for testing the effect without navigating: ThornWindows.demo(index)
  window.ThornWindows = {
    demo(i = 0) {
      const w = document.querySelectorAll(".window")[i];
      if (!w) return;
      w.classList.remove("is-broken");
      spawnBreak(w);
      w.classList.add("is-broken");
    },
  };
})();
