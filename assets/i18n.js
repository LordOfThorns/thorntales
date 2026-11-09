/* i18n.js
   Interface language (RU / EN / DE) and interface strings.
   Load first, before layout.js.

   Language choice: ?lang=xx in the URL > saved choice > browser language > ru.

   In HTML:
     <h1 data-i18n="key">Русский текст</h1>           -> textContent
     <p data-i18n-html="key">текст<br>с br</p>       -> innerHTML (only for strings from this file)
     <input data-i18n-attr="placeholder:key">        -> attribute(s), comma separated
   The Russian text stays in the HTML as a fallback: if a key has no translation, nothing changes.

   In JS:
     I18n.t('key', { n: 3 })   // "{n}" placeholders
     I18n.lang                 // 'ru' | 'en' | 'de'
     I18n.pick(value)          // value = 'text' or { ru, en, de } -> best match

   Content translations (poems etc.) are NOT here, they live next to the originals,
   see data/TRANSLATIONS.md.
*/
(function () {
  'use strict';

  const LANGS = ['ru', 'en', 'de'];
  const DEFAULT = 'ru';
  const LOCALES = { ru: 'ru-RU', en: 'en-GB', de: 'de-DE' };
  const STORE_KEY = 'thorntales.lang';

  const DICT = {
    ru: {
      'nav.home': 'Башня',
      'nav.polarNight': 'Полярная Ночь',
      'nav.texts': 'Тексты',
      'nav.poetry': 'Стихи',
      'nav.projects': 'Проекты',
      'nav.mapGen': 'Генератор карт',
      'nav.login': 'Вход',
      'nav.account': 'Кабинет',
      'nav.menu': 'Меню',
      'nav.mainMenu': 'Основное меню',
      'nav.crumbs': 'Вы здесь',
      'nav.toTop': 'Наверх',
      'nav.signIn': 'Войти',
      'nav.accountTitle': 'Личный кабинет',
      'nav.language': 'Язык',

      'poetry.catalogsFailed': 'Не загрузились каталоги',
      'poems.toc': 'Оглавление',
      'poems.prev': '← предыдущий',
      'poems.next': 'следующий →',
      'poems.all': 'Все каталоги',
      'poems.loadFailed': 'Не загрузилось',
      'poems.untitled': 'Без названия',
      'poems.showOriginal': 'Оригинал',
      'poems.showTranslation': 'Перевод',
      'poems.notTranslated': 'ещё не переведено',
      'poems.translatedCount': 'Переведено {n} из {total}',
      'poems.lyrical': 'Лирика',
      'poems.catalog': 'Стихи: {cat}',

      'projects.mapGen.desc': 'странные карты несуществующих миров',

      'login.empty': 'Введи почту и пароль.',
      'login.progress': 'Входим…',
      'login.bad': 'Неверная почта или пароль.',
      'login.failed': 'Не получилось войти: ',

      'account.moderator': 'модератор',
      'account.nameEmpty': 'Имя не может быть пустым.',
      'account.nameLong': 'Не длиннее 40 символов.',
      'account.saved': 'Сохранено.',
      'account.saveFailed': 'Не сохранилось: ',
      'account.needImage': 'Нужна картинка JPG, PNG или WebP.',
      'account.uploading': 'Загружаем…',
      'account.avatarUpdated': 'Аватарка обновлена.',
      'account.uploadFailed': 'Не получилось загрузить: ',
      'account.avatarRemoved': 'Аватарка убрана.',
      'account.failed': 'Не получилось: ',
      'account.passShort': 'Минимум 8 символов.',
      'account.passMismatch': 'Пароли не совпадают.',
      'account.passChanged': 'Пароль изменён.',

      'comments.title': 'Комментарии',
      'comments.loading': 'Загружаем…',
      'comments.loadFailed': 'Комментарии не загрузились. Попробуй обновить страницу.',
      'comments.empty': 'Пока никто ничего не написал.',
      'comments.placeholder': 'Написать комментарий…',
      'comments.hint': 'Ctrl+Enter — отправить',
      'comments.send': 'Отправить',
      'comments.signInLink': 'Войди',
      'comments.signInRest': ', чтобы написать.',
      'comments.edit': 'Изменить',
      'comments.delete': 'Удалить',
      'comments.cancel': 'Отмена',
      'comments.save': 'Сохранить',
      'comments.someone': 'Кто-то',
      'comments.author': 'автор',
      'comments.edited': 'изменено',
      'comments.sendFailed': 'Не отправилось: ',
      'comments.saveFailed': 'Не сохранилось: ',
      'comments.confirmDelete': 'Удалить комментарий?',
      'comments.deleteFailed': 'Не удалилось: ',
    },

    en: {
      'nav.home': 'Tower',
      'nav.polarNight': 'Polar Night',
      'nav.texts': 'Texts',
      'nav.poetry': 'Poems',
      'nav.projects': 'Projects',
      'nav.mapGen': 'Map generator',
      'nav.login': 'Sign in',
      'nav.account': 'Account',
      'nav.menu': 'Menu',
      'nav.mainMenu': 'Main menu',
      'nav.crumbs': 'You are here',
      'nav.toTop': 'Back to top',
      'nav.signIn': 'Sign in',
      'nav.accountTitle': 'Your account',
      'nav.language': 'Language',

      'home.title': 'You are in the Tower of a Thousand Stained-Glass Windows',
      'home.sub': 'if you remember that station<br>where this land opens up<br>a window in the seven-sided stained-glass tower<br>tell me<br>i will break the glass',
      'home.windows': 'Main navigation',
      'home.window.polarNight': 'Fantasy cycle “Polar Night”',
      'stub.text': 'Oops. You fell through a false window. This section is not ready yet.',
      'polar.kind': 'A fantasy cycle.',

      'poetry.title': 'Poems',
      'poetry.sub': 'collections',
      'poetry.catalogsFailed': 'Collections failed to load',
      'poems.toc': 'Contents',
      'poems.prev': '← previous',
      'poems.next': 'next →',
      'poems.all': 'All collections',
      'poems.loadFailed': 'Failed to load',
      'poems.untitled': 'Untitled',
      'poems.showOriginal': 'Original',
      'poems.showTranslation': 'Translation',
      'poems.notTranslated': 'not translated yet',
      'poems.translatedCount': '{n} of {total} translated',
      'poems.lyrical': 'Lyrics',
      'poems.catalog': 'Poems: {cat}',

      'projects.title': 'Projects',
      'projects.mapGen.desc': 'strange maps of nonexistent worlds',
      'terrain.sub': 'A strange map of a nonexistent world. Change the seed to get a different random generation.',

      'login.title': 'Sign in',
      'login.sub': 'Accounts are handed out by the keeper of the tower. If you don’t have one, write to me.',
      'login.email': 'Email',
      'login.password': 'Password',
      'login.submit': 'Sign in',
      'login.empty': 'Enter your email and password.',
      'login.progress': 'Signing in…',
      'login.bad': 'Wrong email or password.',
      'login.failed': 'Could not sign in: ',

      'account.title': 'Account',
      'account.loading': 'Loading…',
      'account.profile': 'Profile',
      'account.uploadAvatar': 'Upload avatar',
      'account.removeAvatar': 'Remove',
      'account.avatarHint': 'JPG, PNG or WebP. The image will be cropped to a square.',
      'account.nameLabel': 'Name other people see',
      'account.save': 'Save',
      'account.password': 'Password',
      'account.newPass': 'New password',
      'account.repeatPass': 'Repeat',
      'account.changePass': 'Change password',
      'account.logout': 'Sign out',
      'account.moderator': 'moderator',
      'account.nameEmpty': 'The name can’t be empty.',
      'account.nameLong': '40 characters at most.',
      'account.saved': 'Saved.',
      'account.saveFailed': 'Not saved: ',
      'account.needImage': 'Please choose a JPG, PNG or WebP image.',
      'account.uploading': 'Uploading…',
      'account.avatarUpdated': 'Avatar updated.',
      'account.uploadFailed': 'Upload failed: ',
      'account.avatarRemoved': 'Avatar removed.',
      'account.failed': 'Failed: ',
      'account.passShort': 'At least 8 characters.',
      'account.passMismatch': 'Passwords don’t match.',
      'account.passChanged': 'Password changed.',

      'comments.title': 'Comments',
      'comments.loading': 'Loading…',
      'comments.loadFailed': 'Comments failed to load. Try reloading the page.',
      'comments.empty': 'No comments yet.',
      'comments.placeholder': 'Write a comment…',
      'comments.hint': 'Ctrl+Enter to send',
      'comments.send': 'Send',
      'comments.signInLink': 'Sign in',
      'comments.signInRest': ' to comment.',
      'comments.edit': 'Edit',
      'comments.delete': 'Delete',
      'comments.cancel': 'Cancel',
      'comments.save': 'Save',
      'comments.someone': 'Someone',
      'comments.author': 'author',
      'comments.edited': 'edited',
      'comments.sendFailed': 'Not sent: ',
      'comments.saveFailed': 'Not saved: ',
      'comments.confirmDelete': 'Delete this comment?',
      'comments.deleteFailed': 'Could not delete: ',
    },

    de: {
      'nav.home': 'Turm',
      'nav.polarNight': 'Polarnacht',
      'nav.texts': 'Texte',
      'nav.poetry': 'Gedichte',
      'nav.projects': 'Projekte',
      'nav.mapGen': 'Kartengenerator',
      'nav.login': 'Anmelden',
      'nav.account': 'Konto',
      'nav.menu': 'Menü',
      'nav.mainMenu': 'Hauptmenü',
      'nav.crumbs': 'Sie sind hier',
      'nav.toTop': 'Nach oben',
      'nav.signIn': 'Anmelden',
      'nav.accountTitle': 'Dein Konto',
      'nav.language': 'Sprache',

      'home.title': 'Sie sind im Turm der tausend Buntglasfenster',
      'home.sub': 'wenn du dich an jene Station erinnerst<br>wo sich dieses Land öffnet<br>ein Fenster im siebenseitigen Buntglasturm<br>sag es mir<br>ich zerschlage das Glas',
      'home.windows': 'Hauptnavigation',
      'home.window.polarNight': 'Fantasy-Zyklus „Polarnacht“',
      'stub.text': 'Hoppla. Sie sind durch ein falsches Fenster gefallen. Dieser Bereich ist noch nicht fertig.',
      'polar.kind': 'Ein Fantasy-Zyklus.',

      'poetry.title': 'Gedichte',
      'poetry.sub': 'Sammlungen',
      'poetry.catalogsFailed': 'Sammlungen konnten nicht geladen werden',
      'poems.toc': 'Inhalt',
      'poems.prev': '← vorherige',
      'poems.next': 'nächste →',
      'poems.all': 'Alle Sammlungen',
      'poems.loadFailed': 'Konnte nicht geladen werden',
      'poems.untitled': 'Ohne Titel',
      'poems.showOriginal': 'Original',
      'poems.showTranslation': 'Übersetzung',
      'poems.notTranslated': 'noch nicht übersetzt',
      'poems.translatedCount': '{n} von {total} übersetzt',
      'poems.lyrical': 'Lyrik',
      'poems.catalog': 'Gedichte: {cat}',

      'projects.title': 'Projekte',
      'projects.mapGen.desc': 'seltsame Karten nicht existierender Welten',
      'terrain.sub': 'Eine seltsame Karte einer nicht existierenden Welt. Ändern Sie den Seed für eine andere Zufallsgenerierung.',

      'login.title': 'Anmelden',
      'login.sub': 'Konten vergibt der Hüter des Turms. Wenn du keins hast, schreib mir.',
      'login.email': 'E-Mail',
      'login.password': 'Passwort',
      'login.submit': 'Anmelden',
      'login.empty': 'Gib E-Mail und Passwort ein.',
      'login.progress': 'Anmeldung…',
      'login.bad': 'Falsche E-Mail oder falsches Passwort.',
      'login.failed': 'Anmeldung fehlgeschlagen: ',

      'account.title': 'Konto',
      'account.loading': 'Wird geladen…',
      'account.profile': 'Profil',
      'account.uploadAvatar': 'Avatar hochladen',
      'account.removeAvatar': 'Entfernen',
      'account.avatarHint': 'JPG, PNG oder WebP. Das Bild wird quadratisch zugeschnitten.',
      'account.nameLabel': 'Name, den andere sehen',
      'account.save': 'Speichern',
      'account.password': 'Passwort',
      'account.newPass': 'Neues Passwort',
      'account.repeatPass': 'Wiederholen',
      'account.changePass': 'Passwort ändern',
      'account.logout': 'Abmelden',
      'account.moderator': 'Moderator',
      'account.nameEmpty': 'Der Name darf nicht leer sein.',
      'account.nameLong': 'Höchstens 40 Zeichen.',
      'account.saved': 'Gespeichert.',
      'account.saveFailed': 'Nicht gespeichert: ',
      'account.needImage': 'Bitte ein JPG-, PNG- oder WebP-Bild wählen.',
      'account.uploading': 'Wird hochgeladen…',
      'account.avatarUpdated': 'Avatar aktualisiert.',
      'account.uploadFailed': 'Hochladen fehlgeschlagen: ',
      'account.avatarRemoved': 'Avatar entfernt.',
      'account.failed': 'Fehlgeschlagen: ',
      'account.passShort': 'Mindestens 8 Zeichen.',
      'account.passMismatch': 'Die Passwörter stimmen nicht überein.',
      'account.passChanged': 'Passwort geändert.',

      'comments.title': 'Kommentare',
      'comments.loading': 'Wird geladen…',
      'comments.loadFailed': 'Kommentare konnten nicht geladen werden. Lade die Seite neu.',
      'comments.empty': 'Noch keine Kommentare.',
      'comments.placeholder': 'Kommentar schreiben…',
      'comments.hint': 'Strg+Enter zum Senden',
      'comments.send': 'Senden',
      'comments.signInLink': 'Melde dich an',
      'comments.signInRest': ', um zu kommentieren.',
      'comments.edit': 'Bearbeiten',
      'comments.delete': 'Löschen',
      'comments.cancel': 'Abbrechen',
      'comments.save': 'Speichern',
      'comments.someone': 'Jemand',
      'comments.author': 'Autor',
      'comments.edited': 'bearbeitet',
      'comments.sendFailed': 'Nicht gesendet: ',
      'comments.saveFailed': 'Nicht gespeichert: ',
      'comments.confirmDelete': 'Diesen Kommentar löschen?',
      'comments.deleteFailed': 'Löschen fehlgeschlagen: ',
    },
  };

  // ---------- Current language ----------

  function detect() {
    const fromUrl = new URLSearchParams(location.search).get('lang');
    if (LANGS.includes(fromUrl)) {
      save(fromUrl);
      return fromUrl;
    }
    const saved = read();
    if (LANGS.includes(saved)) return saved;
    for (const l of navigator.languages || [navigator.language || '']) {
      const short = String(l).slice(0, 2).toLowerCase();
      if (LANGS.includes(short)) return short;
    }
    return DEFAULT;
  }

  // localStorage can be blocked (private mode etc.): never let that break the page
  function read() {
    try { return localStorage.getItem(STORE_KEY); } catch (e) { return null; }
  }
  function save(l) {
    try { localStorage.setItem(STORE_KEY, l); } catch (e) { /* ignore */ }
  }

  const lang = detect();
  document.documentElement.lang = lang;

  // ---------- Lookup ----------

  function t(key, vars) {
    let s = DICT[lang][key];
    if (s == null) s = DICT[DEFAULT][key];
    if (s == null) return key;
    if (vars) s = s.replace(/\{(\w+)\}/g, (m, k) => (k in vars ? String(vars[k]) : m));
    return s;
  }

  function has(key) {
    return DICT[lang][key] != null;
  }

  // 'text' or { ru: '...', en: '...' } -> best string for the current language
  function pick(value) {
    if (value == null || typeof value !== 'object') return value;
    return value[lang] ?? value[DEFAULT] ?? Object.values(value)[0];
  }

  // Localized field of a data object:
  //   { label: 'Сальвии', i18n: { en: { label: 'Salvias' } } }  ->  field(obj, 'label')
  function field(obj, name) {
    if (!obj) return undefined;
    return obj.i18n?.[lang]?.[name] ?? obj[name];
  }

  function apply(root) {
    const scope = root || document;
    scope.querySelectorAll('[data-i18n]').forEach((el) => {
      const key = el.dataset.i18n;
      if (has(key)) el.textContent = t(key);
    });
    scope.querySelectorAll('[data-i18n-html]').forEach((el) => {
      const key = el.dataset.i18nHtml;
      if (has(key)) el.innerHTML = t(key);
    });
    scope.querySelectorAll('[data-i18n-attr]').forEach((el) => {
      for (const pair of el.dataset.i18nAttr.split(',')) {
        const [attr, key] = pair.split(':').map((s) => s.trim());
        if (attr && key && has(key)) el.setAttribute(attr, t(key));
      }
    });
    // <title data-i18n> lives in <head>
    const title = document.querySelector('title[data-i18n]');
    if (title && has(title.dataset.i18n)) {
      document.title = `${t(title.dataset.i18n)} · Thorntales`;
    }
  }

  function setLang(next) {
    if (!LANGS.includes(next) || next === lang) return;
    save(next);
    // Keep the URL shareable: update ?lang only if it was there
    const url = new URL(location.href);
    if (url.searchParams.has('lang')) url.searchParams.set('lang', next);
    location.replace(url.toString());
  }

  window.I18n = {
    langs: LANGS,
    defaultLang: DEFAULT,
    lang,
    locale: LOCALES[lang],
    t,
    has,
    pick,
    field,
    apply,
    setLang,
  };

  apply();
})();
