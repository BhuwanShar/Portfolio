/* ═══════════════════════════════════════════════════════════════════
   BHUWAN SHARMA — Portfolio scripts (minimal)
   ═══════════════════════════════════════════════════════════════════ */

(() => {
  'use strict';

  // ─── FONTS ───
  if (!document.querySelector('link[href*="Source+Serif"]')) {
    const link = document.createElement('link');
    link.rel = 'stylesheet';
    link.href = 'https://fonts.googleapis.com/css2?family=Source+Serif+Pro:wght@400;600;700&family=Source+Sans+Pro:wght@300;400;600;700&family=JetBrains+Mono:wght@400;600&display=swap';
    document.head.appendChild(link);
  }

  // ─── THEME ───
  const html = document.documentElement;
  const themeBtn = document.getElementById('themeBtn');
  const stored = localStorage.getItem('theme');
  if (stored) html.dataset.theme = stored;
  else if (matchMedia('(prefers-color-scheme: dark)').matches) html.dataset.theme = 'dark';

  if (themeBtn) {
    themeBtn.addEventListener('click', () => {
      const next = html.dataset.theme === 'dark' ? 'light' : 'dark';
      html.dataset.theme = next;
      localStorage.setItem('theme', next);
    });
  }

  // ─── MOBILE MENU ───
  const mobBtn = document.getElementById('mobBtn');
  const mobMenu = document.getElementById('mobMenu');
  if (mobBtn && mobMenu) {
    mobBtn.addEventListener('click', () => {
      mobBtn.classList.toggle('open');
      mobMenu.classList.toggle('open');
    });
    mobMenu.querySelectorAll('a').forEach(a => {
      a.addEventListener('click', () => {
        mobBtn.classList.remove('open');
        mobMenu.classList.remove('open');
      });
    });
  }

  // ─── ACTIVE NAV LINK ───
  const path = location.pathname.split('/').pop() || 'index.html';
  document.querySelectorAll('.nav-links a, .mob-menu a').forEach(a => {
    const href = a.getAttribute('href');
    if (href === path || (path === '' && href === 'index.html')) {
      a.classList.add('active');
    }
  });

  // ─── TOAST helper ───
  window.showToast = (msg) => {
    let toast = document.getElementById('toast');
    if (!toast) {
      toast = document.createElement('div');
      toast.id = 'toast';
      toast.className = 'toast';
      document.body.appendChild(toast);
    }
    toast.textContent = msg;
    toast.classList.add('show');
    setTimeout(() => toast.classList.remove('show'), 2200);
  };

  // ─── Copy email helper ───
  document.querySelectorAll('[data-copy]').forEach(el => {
    el.addEventListener('click', (e) => {
      e.preventDefault();
      const text = el.dataset.copy;
      navigator.clipboard.writeText(text).then(() => {
        window.showToast('Copied to clipboard');
      });
    });
  });

})();