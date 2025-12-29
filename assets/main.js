(function () {
  function ensureShardsLayer(winEl) {
    let layer = winEl.querySelector(".shards");
    if (!layer) {
      layer = document.createElement("div");
      layer.className = "shards";
      winEl.appendChild(layer);
    }
    return layer;
  }

  function spawnShards(winEl, count = 16) {
    const layer = ensureShardsLayer(winEl);

    // чистим старые осколки, если кликали раньше
    layer.innerHTML = "";

    const rect = winEl.getBoundingClientRect();
    const cx = rect.width / 2;
    const cy = rect.height * 0.42; // "точка удара" чуть выше центра, выглядит лучше

    for (let i = 0; i < count; i++) {
      const shard = document.createElement("div");
      shard.className = "shard";

      // стартовая позиция вокруг точки удара
      const sx = cx + (Math.random() - 0.5) * 28;
      const sy = cy + (Math.random() - 0.5) * 28;

      // направление: наружу + немного вверх
      const angle = (Math.random() * Math.PI * 2);
      const speed = 90 + Math.random() * 160;

      const dx = Math.cos(angle) * speed;
      const dy = Math.sin(angle) * speed - (50 + Math.random() * 60);

      // размеры и вращение
      const w = 10 + Math.random() * 18;
      const h = 12 + Math.random() * 22;
      const r0 = (Math.random() * 40 - 20).toFixed(1);
      const r1 = (180 + Math.random() * 220).toFixed(1);

      shard.style.left = `${sx}px`;
      shard.style.top = `${sy}px`;

      shard.style.setProperty("--dx", `${dx.toFixed(1)}px`);
      shard.style.setProperty("--dy", `${dy.toFixed(1)}px`);
      shard.style.setProperty("--w", `${w.toFixed(1)}px`);
      shard.style.setProperty("--h", `${h.toFixed(1)}px`);
      shard.style.setProperty("--r0", `${r0}deg`);
      shard.style.setProperty("--r1", `${r1}deg`);

      layer.appendChild(shard);
    }

    // убрать осколки после анимации, чтобы DOM не разрастался
    setTimeout(() => {
      layer.innerHTML = "";
    }, 800);
  }

  function onReady(fn) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", fn);
    } else {
      fn();
    }
  }

  onReady(() => {
    document.querySelectorAll(".window").forEach((w) => {
      w.addEventListener("click", (e) => {
        e.preventDefault();

        // если уже разбито: сразу перейти
        if (w.classList.contains("is-broken")) {
          const href = w.getAttribute("href");
          if (href) window.location.href = href;
          return;
        }

        // 1) осколки
        spawnShards(w, 18);

        // 2) включаем состояние (прячет сетку, зумит картинку внутри)
        w.classList.add("is-broken");

        // 3) перейти по ссылке через паузу
        const href = w.getAttribute("href");
        if (href) {
          setTimeout(() => {
            window.location.href = href;
          }, 1000);
        }
      });
    });
  });
})();
