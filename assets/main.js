/* main.js — stained-glass windows on the home page: click breaks the glass, then
   navigates. The break effect:
     - impact flash at the hit point;
     - shards cut from the window's own picture (each shows the piece it came from),
       tumbling in 3D, glinting as they turn, flying out and falling in an arc;
     - a few twinkling sparkles.
   Only transforms, opacity and filter are animated (Web Animations API), ~40 elements.
*/
(function () {
  "use strict";

  const BREAK_DELAY_MS = 650;  // navigation starts after this
  const SHARDS = 30;
  const SPARKLES = 12;
  const GLASS_INSET = 12;      // .glass is inset by 12px inside .window (style.css)
  const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");

  const rand = (a, b) => a + Math.random() * (b - a);

  function ensureShardsLayer(winEl) {
    let layer = winEl.querySelector(".shards");
    if (!layer) {
      layer = document.createElement("div");
      layer.className = "shards";
      winEl.appendChild(layer);
    }
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
    const layer = ensureShardsLayer(winEl);
    layer.innerHTML = "";

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

    // impact flash
    const flash = document.createElement("div");
    flash.className = "impact";
    flash.style.left = `${cx}px`;
    flash.style.top = `${cy}px`;
    layer.appendChild(flash);
    flash.animate(
      [
        { transform: "translate(-50%, -50%) scale(0.2)", opacity: 0.95 },
        { transform: "translate(-50%, -50%) scale(1.6)", opacity: 0 },
      ],
      { duration: 220, easing: "cubic-bezier(.2,.7,.3,1)", fill: "forwards" }
    );

    // shards
    for (let i = 0; i < SHARDS; i++) {
      const shard = document.createElement("div");
      shard.className = "shard";

      // size: a few big pieces, lots of crumbs
      const size = 7 + Math.pow(Math.random(), 1.8) * 42;
      const w = size * rand(0.7, 1.3), h = size * rand(0.7, 1.3);
      // start: around the impact point, wider for small pieces
      const spread = 18 + Math.max(0, 48 - size) * 0.8;
      const x = clamp(cx + rand(-spread, spread) - w / 2, GLASS_INSET, W - GLASS_INSET - w);
      const y = clamp(cy + rand(-spread, spread) - h / 2, GLASS_INSET, H - GLASS_INSET - h);

      shard.style.left = `${x}px`;
      shard.style.top = `${y}px`;
      shard.style.width = `${w}px`;
      shard.style.height = `${h}px`;
      shard.style.clipPath = shardShape();
      if (pic) {
        // the piece of the picture that was under this shard
        shard.style.backgroundImage =
          `linear-gradient(135deg, rgba(255,255,255,.6), rgba(255,255,255,.05) 35%, rgba(255,255,255,0) 55%, rgba(255,255,255,.25) 80%, rgba(255,255,255,.5)), url("${pic.url}")`;
        shard.style.backgroundSize = `100% 100%, ${bgW}px ${bgH}px`;
        shard.style.backgroundPosition =
          `0 0, ${bgX + GLASS_INSET - x}px ${bgY + GLASS_INSET - y}px`;
      }
      layer.appendChild(shard);

      // flight: out from the impact, up, then falling (gravity); small pieces fly further
      const ang = Math.atan2(y + h / 2 - cy, x + w / 2 - cx) + rand(-0.5, 0.5);
      const power = rand(130, 260) * (1.35 - size / 60);
      const dx = Math.cos(ang) * power;
      const up = rand(50, 120) * (1.3 - size / 70);
      const fall = rand(120, 220);
      const ax = rand(-1, 1).toFixed(2), ay = rand(-1, 1).toFixed(2), az = rand(-0.3, 0.3).toFixed(2);
      const spin = rand(240, 620) * (Math.random() < 0.5 ? -1 : 1);
      const dur = rand(560, 760);
      const glintAt = rand(0.2, 0.55); // when the face catches the light

      // flight: fast burst decelerating to the top of the arc, then accelerating fall;
      // fully visible for most of the flight, fading only at the end
      shard.animate(
        [
          { offset: 0, transform: "translate(0,0) rotate3d(0,0,1,0deg) scale(1)", easing: "cubic-bezier(.15,.8,.4,1)" },
          { offset: 0.35, transform: `translate(${dx * 0.55}px, ${Math.sin(ang) * power * 0.45 - up}px) rotate3d(${ax},${ay},${az},${spin * 0.35}deg) scale(1)`, easing: "cubic-bezier(.5,0,.9,.55)" },
          { offset: 1, transform: `translate(${dx}px, ${Math.sin(ang) * power * 0.6 - up + fall}px) rotate3d(${ax},${ay},${az},${spin}deg) scale(0.85)` },
        ],
        { duration: dur, fill: "forwards" }
      );
      shard.animate(
        [{ opacity: 1 }, { offset: 0.72, opacity: 1 }, { opacity: 0 }],
        { duration: dur, fill: "forwards" }
      );
      // glint: a quick flash of brightness while tumbling; slightly brighter than the
      // window at rest so the shards read as glass from the first frame
      shard.animate(
        [
          { offset: 0, filter: "brightness(1.25)" },
          { offset: Math.max(0.05, glintAt - 0.08), filter: "brightness(1.15)" },
          { offset: glintAt, filter: "brightness(2.6) saturate(.6)" },
          { offset: Math.min(0.95, glintAt + 0.1), filter: "brightness(1.15)" },
          { offset: 1, filter: "brightness(1.1)" },
        ],
        { duration: dur, fill: "forwards" }
      );
    }

    // sparkles
    for (let i = 0; i < SPARKLES; i++) {
      const s = document.createElement("div");
      s.className = "sparkle";
      const size = rand(6, 13);
      s.style.width = s.style.height = `${size}px`;
      s.style.left = `${cx - size / 2 + rand(-14, 14)}px`;
      s.style.top = `${cy - size / 2 + rand(-14, 14)}px`;
      layer.appendChild(s);
      const ang = rand(0, Math.PI * 2);
      const dist = rand(40, 120);
      const delay = rand(0, 160);
      s.animate(
        [
          { transform: "translate(0,0) scale(0) rotate(0deg)", opacity: 0 },
          { offset: 0.25, transform: `translate(${Math.cos(ang) * dist * 0.4}px, ${Math.sin(ang) * dist * 0.4}px) scale(1.15) rotate(45deg)`, opacity: 1 },
          { offset: 0.5, transform: `translate(${Math.cos(ang) * dist * 0.7}px, ${Math.sin(ang) * dist * 0.7}px) scale(0.5) rotate(90deg)`, opacity: 0.6 },
          { offset: 0.7, transform: `translate(${Math.cos(ang) * dist * 0.85}px, ${Math.sin(ang) * dist * 0.85}px) scale(1) rotate(135deg)`, opacity: 1 },
          { transform: `translate(${Math.cos(ang) * dist}px, ${Math.sin(ang) * dist + 20}px) scale(0) rotate(180deg)`, opacity: 0 },
        ],
        { duration: rand(420, 620), delay, easing: "ease-out", fill: "both" }
      );
    }

    // clean up after the animation so the DOM does not grow
    setTimeout(() => { layer.innerHTML = ""; }, 1000);
  }

  function clamp(v, a, b) { return v < a ? a : v > b ? b : v; }

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
      const layer = w.querySelector(".shards");
      if (layer) layer.innerHTML = "";
    });
  });

  onReady(() => {
    document.querySelectorAll(".window").forEach((w) => {
      pictureOf(w); // preload the picture size for the shards

      w.addEventListener("click", (e) => {
        // Let the browser handle "open in new tab/window" clicks and users who disabled animations
        if (e.button !== 0 || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) return;
        if (reducedMotion.matches) return;

        e.preventDefault();
        const href = w.getAttribute("href");

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
