/* rich-editor.js — shared by the text and poem editors (texts-read.js, poems.js).

   RichEdit.sanitize(html)          only a small, safe set of tags is kept
   RichEdit.plain(html)             the text without tags
   RichEdit.field(container, html, { lines })
       A field with Bold / Italic buttons (also Ctrl+B / Ctrl+I) and an HTML mode that
       shows the same text as source. Returns { get() } -> sanitized HTML.
       lines: true for poems — the text keeps its line breaks as typed (white-space:
       pre-wrap), Enter adds a new line instead of a new paragraph.
*/
(function () {
  'use strict';

  const t = (key, fallback) => (window.I18n && I18n.has(key) ? I18n.t(key) : fallback);

  const TAGS = {
    B: 'b', STRONG: 'b', I: 'i', EM: 'i', U: 'u', S: 's', BR: 'br', P: 'p', DIV: 'p',
    BLOCKQUOTE: 'blockquote', HR: 'hr', A: 'a', UL: 'ul', OL: 'ol', LI: 'li', H3: 'h3', H4: 'h4',
  };
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (ch) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;',
  }[ch]));

  function sanitize(html) {
    const doc = new DOMParser().parseFromString(`<div>${html || ''}</div>`, 'text/html');
    const walk = (node) => {
      let out = '';
      for (const c of node.childNodes) {
        if (c.nodeType === 3) { out += esc(c.textContent); continue; }
        if (c.nodeType !== 1) continue;
        if (/^(SCRIPT|STYLE|TEMPLATE|IFRAME|OBJECT|EMBED|NOSCRIPT|SVG|MATH)$/i.test(c.tagName)) continue; // dropped with their content
        const tag = TAGS[c.tagName];
        const inner = walk(c);
        if (!tag) { out += inner; continue; }
        if (tag === 'br' || tag === 'hr') { out += `<${tag}>`; continue; }
        if (tag === 'a') {
          const href = c.getAttribute('href') || '';
          out += /^(https?:|mailto:)/i.test(href)
            ? `<a href="${esc(href)}" target="_blank" rel="noopener">${inner}</a>` : inner;
          continue;
        }
        if (tag === 'p') {
          // centred paragraphs (class "c" from the import, or text-align from pasting)
          const centred = c.classList.contains('c') || /text-align:\s*center/i.test(c.getAttribute('style') || '');
          if (!inner.trim()) continue;
          out += `<p${centred ? ' class="c"' : ''}>${inner}</p>`;
          continue;
        }
        out += `<${tag}>${inner}</${tag}>`;
      }
      return out;
    };
    return walk(doc.body.firstChild);
  }

  const plain = (html) => new DOMParser().parseFromString(html || '', 'text/html').body.textContent || '';

  function field(container, html, opts = {}) {
    const lines = !!opts.lines;
    container.innerHTML = `
      <div class="rich${lines ? ' rich--lines' : ''}">
        <div class="rich__bar" role="toolbar">
          <button type="button" class="rich__btn" data-cmd="bold" title="${esc(t('texts.bold', 'Жирный (Ctrl+B)'))}"><b>Ж</b></button>
          <button type="button" class="rich__btn" data-cmd="italic" title="${esc(t('texts.italic', 'Курсив (Ctrl+I)'))}"><i>К</i></button>
          <span class="rich__sep"></span>
          <button type="button" class="rich__btn rich__btn--wide" data-cmd="html" aria-pressed="false">HTML</button>
        </div>
        <div class="rich__area input" contenteditable="true"></div>
        <textarea class="rich__source input" spellcheck="false" hidden></textarea>
      </div>`;
    const area = container.querySelector('.rich__area');
    const src = container.querySelector('.rich__source');
    const htmlBtn = container.querySelector('[data-cmd="html"]');
    area.innerHTML = sanitize(html);
    let htmlMode = false;
    if (!lines) { try { document.execCommand('defaultParagraphSeparator', false, 'p'); } catch (e) { /* old browsers */ } }

    const bar = container.querySelector('.rich__bar');
    bar.addEventListener('mousedown', (e) => e.preventDefault()); // keep the selection
    bar.addEventListener('click', (e) => {
      const cmd = e.target.closest('[data-cmd]')?.dataset.cmd;
      if (!cmd) return;
      if (cmd === 'html') {
        htmlMode = !htmlMode;
        if (htmlMode) {
          const s = sanitize(area.innerHTML);
          src.value = lines ? s : s.replace(/(<\/(p|blockquote|h3|h4|li|ul|ol)>|<hr>)/g, '$1\n');
          src.style.height = `${Math.max(area.offsetHeight, 200)}px`;
        } else {
          area.innerHTML = sanitize(src.value);
        }
        area.hidden = htmlMode;
        src.hidden = !htmlMode;
        htmlBtn.setAttribute('aria-pressed', String(htmlMode));
        container.querySelectorAll('[data-cmd="bold"], [data-cmd="italic"]').forEach((b) => { b.disabled = htmlMode; });
        return;
      }
      if (!htmlMode) { area.focus(); document.execCommand(cmd); }
    });

    // pasting keeps only the allowed formatting (a poem keeps its line breaks)
    area.addEventListener('paste', (e) => {
      e.preventDefault();
      const cb = e.clipboardData;
      const pasted = cb.getData('text/html');
      if (pasted && !lines) document.execCommand('insertHTML', false, sanitize(pasted));
      else document.execCommand('insertText', false, cb.getData('text/plain'));
    });
    area.addEventListener('keydown', (e) => {
      if ((e.ctrlKey || e.metaKey) && !e.altKey && ['b', 'i'].includes(e.key.toLowerCase())) {
        e.preventDefault();
        document.execCommand(e.key.toLowerCase() === 'b' ? 'bold' : 'italic');
      } else if (lines && e.key === 'Enter') {
        // a poem: a plain line break in the text, not a new <div>/<p>
        e.preventDefault();
        document.execCommand('insertText', false, '\n');
      }
    });
    return { get: () => sanitize(htmlMode ? src.value : area.innerHTML) };
  }

  window.RichEdit = { sanitize, plain, field, esc };
})();
