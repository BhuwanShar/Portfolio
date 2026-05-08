/* ══════════════════════════════════════════════════
   BHUWAN SHARMA — PORTFOLIO SCRIPT
   Computational Biology 
   ══════════════════════════════════════════════════ */

(function () {
    'use strict';

    // ── THEME ─────────────────────────────────────
    function initTheme() {
        const html = document.documentElement;
        const toggle = document.getElementById('themeToggle');
        const saved = localStorage.getItem('theme') || 'dark';

        html.setAttribute('data-theme', saved);

        if (toggle) {
            toggle.addEventListener('click', () => {
                const current = html.getAttribute('data-theme');
                const next = current === 'dark' ? 'light' : 'dark';
                html.setAttribute('data-theme', next);
                localStorage.setItem('theme', next);
            });
        }
    }

    // ── NAV ───────────────────────────────────────
    function initNav() {
        const nav = document.getElementById('mainNav');
        const hamburger = document.getElementById('hamburger');
        const navLinks = document.getElementById('navLinks');

        // Scroll shadow
        let lastScroll = 0;
        window.addEventListener('scroll', () => {
            const scrolled = window.scrollY > 20;
            nav.classList.toggle('scrolled', scrolled);
            lastScroll = window.scrollY;
        }, { passive: true });

        // Mobile toggle
        if (hamburger && navLinks) {
            hamburger.addEventListener('click', () => {
                hamburger.classList.toggle('open');
                navLinks.classList.toggle('open');
            });

            navLinks.querySelectorAll('a').forEach(link => {
                link.addEventListener('click', () => {
                    hamburger.classList.remove('open');
                    navLinks.classList.remove('open');
                });
            });

            document.addEventListener('click', (e) => {
                if (!hamburger.contains(e.target) && !navLinks.contains(e.target)) {
                    hamburger.classList.remove('open');
                    navLinks.classList.remove('open');
                }
            });
        }
    }

    // ── SCROLL SPY ────────────────────────────────
    function initScrollSpy() {
        const sections = document.querySelectorAll('section[id]');
        const navLinks = document.querySelectorAll('.nav-links a[data-section]');

        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    const id = entry.target.id;
                    navLinks.forEach(link => {
                        link.classList.toggle('active', link.dataset.section === id);
                    });
                }
            });
        }, {
            rootMargin: '-25% 0px -65% 0px',
            threshold: 0
        });

        sections.forEach(s => observer.observe(s));
    }

    // ── SCROLL REVEAL ─────────────────────────────
    function initReveal() {
        const elements = document.querySelectorAll('.section');
        const pipelineSteps = document.querySelectorAll('.pipeline-step');

        const sectionObserver = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    entry.target.classList.add('visible');
                    sectionObserver.unobserve(entry.target);
                }
            });
        }, {
            threshold: 0.05,
            rootMargin: '0px 0px -60px 0px'
        });

        elements.forEach(el => sectionObserver.observe(el));

        // Pipeline steps staggered
        const pipeObserver = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    const step = parseInt(entry.target.dataset.step) || 0;
                    setTimeout(() => {
                        entry.target.classList.add('visible');
                    }, step * 100);
                    pipeObserver.unobserve(entry.target);
                }
            });
        }, { threshold: 0.2 });

        pipelineSteps.forEach(s => pipeObserver.observe(s));
    }

    // ── SMOOTH SCROLL ─────────────────────────────
    function initSmoothScroll() {
        document.querySelectorAll('a[href^="#"]').forEach(link => {
            link.addEventListener('click', (e) => {
                const id = link.getAttribute('href').slice(1);
                if (!id) return;
                const target = document.getElementById(id);
                if (target) {
                    e.preventDefault();
                    const navH = document.querySelector('.nav')?.offsetHeight || 64;
                    const top = target.getBoundingClientRect().top + window.scrollY - navH - 20;
                    window.scrollTo({ top, behavior: 'smooth' });
                }
            });
        });
    }

    // ── PHOTO FALLBACK ────────────────────────────
    function initPhoto() {
        const img = document.getElementById('heroPhoto');
        const fallback = document.getElementById('photoFallback');
        if (!img || !fallback) return;

        img.addEventListener('error', () => {
            img.style.display = 'none';
            fallback.style.display = 'flex';
        });
    }

    // ── CURSOR GLOW ───────────────────────────────
    function initCursorGlow() {
        const glow = document.getElementById('cursorGlow');
        if (!glow) return;

        let mouseX = 0, mouseY = 0;
        let glowX = 0, glowY = 0;

        document.addEventListener('mousemove', (e) => {
            mouseX = e.clientX;
            mouseY = e.clientY;
        }, { passive: true });

        function animate() {
            glowX += (mouseX - glowX) * 0.08;
            glowY += (mouseY - glowY) * 0.08;
            glow.style.left = glowX + 'px';
            glow.style.top = glowY + 'px';
            requestAnimationFrame(animate);
        }
        animate();
    }

    // ── STAT COUNTER ──────────────────────────────
    function initCounters() {
        const counters = document.querySelectorAll('[data-count]');
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    const el = entry.target;
                    const target = parseInt(el.dataset.count);
                    animateCounter(el, target);
                    observer.unobserve(el);
                }
            });
        }, { threshold: 0.5 });

        counters.forEach(c => observer.observe(c));
    }

    function animateCounter(el, target) {
        let current = 0;
        const duration = 1500;
        const step = target / (duration / 16);

        function update() {
            current += step;
            if (current >= target) {
                el.textContent = target + '+';
                return;
            }
            el.textContent = Math.floor(current);
            requestAnimationFrame(update);
        }
        update();
    }

    // ── PROJECT FILTERS ───────────────────────────
    function initFilters() {
        const buttons = document.querySelectorAll('.filter-btn');
        const cards = document.querySelectorAll('.project-card');

        buttons.forEach(btn => {
            btn.addEventListener('click', () => {
                buttons.forEach(b => b.classList.remove('active'));
                btn.classList.add('active');

                const filter = btn.dataset.filter;
                cards.forEach(card => {
                    if (filter === 'all' || card.dataset.tags === filter) {
                        card.classList.remove('hidden');
                        card.style.opacity = '0';
                        card.style.transform = 'translateY(10px)';
                        requestAnimationFrame(() => {
                            card.style.transition = 'opacity 0.4s ease, transform 0.4s ease';
                            card.style.opacity = '1';
                            card.style.transform = 'translateY(0)';
                        });
                    } else {
                        card.classList.add('hidden');
                    }
                });
            });
        });
    }

    // ── PROJECT EXPAND ────────────────────────────
    function initExpand() {
        document.querySelectorAll('.expand-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                const targetId = btn.dataset.target;
                const details = document.getElementById(targetId);
                if (!details) return;

                const isOpen = details.classList.contains('open');
                details.classList.toggle('open');
                btn.classList.toggle('open');
                btn.querySelector('span').textContent = isOpen ? 'Details' : 'Collapse';
            });
        });
    }

    // ── COPY EMAIL ────────────────────────────────
    function initCopy() {
        document.querySelectorAll('.copy-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.preventDefault();
                e.stopPropagation();
                const text = btn.dataset.copy;
                navigator.clipboard.writeText(text).then(() => {
                    showToast('Email copied to clipboard');
                }).catch(() => {
                    showToast('Failed to copy');
                });
            });
        });
    }

    // ── TOAST ─────────────────────────────────────
    function showToast(msg) {
        const toast = document.getElementById('toast');
        const toastMsg = document.getElementById('toastMsg');
        if (!toast || !toastMsg) return;
        toastMsg.textContent = msg;
        toast.classList.add('show');
        setTimeout(() => toast.classList.remove('show'), 2500);
    }

    // ── CONTACT FORM ──────────────────────────────
    function initForm() {
        const form = document.getElementById('contactForm');
        if (!form) return;

        form.addEventListener('submit', (e) => {
            e.preventDefault();
            const name = document.getElementById('formName').value;
            const email = document.getElementById('formEmail').value;
            const subject = document.getElementById('formSubject').value;
            const message = document.getElementById('formMessage').value;

            // Mailto fallback
            const mailtoLink = `mailto:bhuwangautam09@gmail.com?subject=${encodeURIComponent(subject + ' — from ' + name)}&body=${encodeURIComponent(message + '\n\nFrom: ' + name + ' (' + email + ')')}`;
            window.location.href = mailtoLink;

            showToast('Opening email client...');
            form.reset();
        });
    }

    // ── NETWORK BACKGROUND CANVAS ─────────────────
    function initNetworkBackground() {
        const canvas = document.getElementById('networkCanvas');
        if (!canvas) return;

        const ctx = canvas.getContext('2d');
        let width, height;
        let particles = [];
        let animationId;

        function resize() {
            width = canvas.width = window.innerWidth;
            height = canvas.height = window.innerHeight;
        }

        function createParticles() {
            const count = Math.min(Math.floor((width * height) / 25000), 60);
            particles = [];
            for (let i = 0; i < count; i++) {
                particles.push({
                    x: Math.random() * width,
                    y: Math.random() * height,
                    vx: (Math.random() - 0.5) * 0.3,
                    vy: (Math.random() - 0.5) * 0.3,
                    size: Math.random() * 2 + 1,
                    opacity: Math.random() * 0.4 + 0.1
                });
            }
        }

        function draw() {
            ctx.clearRect(0, 0, width, height);

            const theme = document.documentElement.getAttribute('data-theme');
            const nodeColor = theme === 'dark' ? '0, 229, 255' : '8, 145, 178';
            const lineColor = theme === 'dark' ? '0, 229, 255' : '8, 145, 178';

            // Draw connections
            for (let i = 0; i < particles.length; i++) {
                for (let j = i + 1; j < particles.length; j++) {
                    const dx = particles[i].x - particles[j].x;
                    const dy = particles[i].y - particles[j].y;
                    const dist = Math.sqrt(dx * dx + dy * dy);
                    if (dist < 150) {
                        const alpha = (1 - dist / 150) * 0.08;
                        ctx.strokeStyle = `rgba(${lineColor}, ${alpha})`;
                        ctx.lineWidth = 0.5;
                        ctx.beginPath();
                        ctx.moveTo(particles[i].x, particles[i].y);
                        ctx.lineTo(particles[j].x, particles[j].y);
                        ctx.stroke();
                    }
                }
            }

            // Draw particles
            particles.forEach(p => {
                ctx.fillStyle = `rgba(${nodeColor}, ${p.opacity})`;
                ctx.beginPath();
                ctx.arc(p.x, p.y, p.size, 0, Math.PI * 2);
                ctx.fill();

                // Move
                p.x += p.vx;
                p.y += p.vy;

                if (p.x < 0 || p.x > width) p.vx *= -1;
                if (p.y < 0 || p.y > height) p.vy *= -1;
            });

            animationId = requestAnimationFrame(draw);
        }

        resize();
        createParticles();
        draw();

        window.addEventListener('resize', () => {
            resize();
            createParticles();
        });
    }

    // ── INTERACTIVE NETWORK VISUALIZATION ─────────
    function initNetworkViz() {
        const canvas = document.getElementById('networkViz');
        if (!canvas) return;

        const ctx = canvas.getContext('2d');

        let rect = canvas.getBoundingClientRect();
        let mouseX = -1000;
        let mouseY = -1000;
        let animationFrame;

        // ── NETWORK CONFIG ─────────────────────────
        const nodeCount = 42;

        const nodes = [];
        const edges = [];

        const nodeTypes = [
            'hub',
            'regulatory',
            'biomarker',
            'pathway'
        ];

        const typeColors = {
            hub: null,
            regulatory: null,
            biomarker: null,
            pathway: null
        };

        const genePrefixes = [
            'BRCA',
            'TP53',
            'EGFR',
            'STAT',
            'MAPK',
            'AKT',
            'MYC',
            'CDK',
            'SMAD',
            'FOXO'
        ];

        // ── THEME COLORS ───────────────────────────
        function getColors() {
            const theme =
                document.documentElement.getAttribute('data-theme');

            if (theme === 'dark') {
                typeColors.hub = '#00e5ff';
                typeColors.regulatory = '#34d399';
                typeColors.biomarker = '#a78bfa';
                typeColors.pathway = '#f59e0b';
            } else {
                typeColors.hub = '#0891b2';
                typeColors.regulatory = '#059669';
                typeColors.biomarker = '#7c3aed';
                typeColors.pathway = '#d97706';
            }
        }

        // ── RESPONSIVE CANVAS ──────────────────────
        function resizeCanvas() {
            const container = canvas.parentElement;

            canvas.width = container.clientWidth - 48;
            canvas.height = 460;

            rect = canvas.getBoundingClientRect();
        }

        // ── GENERATE NETWORK ───────────────────────
        function createNetwork() {
            nodes.length = 0;
            edges.length = 0;

            for (let i = 0; i < nodeCount; i++) {

                const type =
                    nodeTypes[
                        Math.floor(Math.random() * nodeTypes.length)
                    ];

                const gene =
                    genePrefixes[
                        Math.floor(Math.random() * genePrefixes.length)
                    ];

                nodes.push({
                    x: 60 + Math.random() * (canvas.width - 120),
                    y: 60 + Math.random() * (canvas.height - 120),

                    vx: (Math.random() - 0.5) * 0.25,
                    vy: (Math.random() - 0.5) * 0.25,

                    type,

                    size:
                        type === 'hub'
                            ? 7 + Math.random() * 4
                            : 4 + Math.random() * 4,

                    label:
                        `${gene}${Math.floor(Math.random() * 9) + 1}`,

                    highlighted: false
                });
            }

            // ── CREATE EDGES ───────────────────────
            for (let i = 0; i < nodeCount; i++) {

                const connections =
                    2 + Math.floor(Math.random() * 4);

                for (let c = 0; c < connections; c++) {

                    const j =
                        Math.floor(Math.random() * nodeCount);

                    if (
                        j !== i &&
                        !edges.some(
                            e =>
                                (e.from === i && e.to === j) ||
                                (e.from === j && e.to === i)
                        )
                    ) {
                        edges.push({
                            from: i,
                            to: j
                        });
                    }
                }
            }
        }

        // ── BACKGROUND GRID ────────────────────────
        function drawGrid(theme) {

            const gridColor =
                theme === 'dark'
                    ? 'rgba(255,255,255,0.025)'
                    : 'rgba(0,0,0,0.03)';

            ctx.strokeStyle = gridColor;
            ctx.lineWidth = 1;

            const spacing = 40;

            for (let x = 0; x < canvas.width; x += spacing) {
                ctx.beginPath();
                ctx.moveTo(x, 0);
                ctx.lineTo(x, canvas.height);
                ctx.stroke();
            }

            for (let y = 0; y < canvas.height; y += spacing) {
                ctx.beginPath();
                ctx.moveTo(0, y);
                ctx.lineTo(canvas.width, y);
                ctx.stroke();
            }
        }

        // ── MAIN DRAW LOOP ─────────────────────────
        function drawNetwork() {

            ctx.clearRect(
                0,
                0,
                canvas.width,
                canvas.height
            );

            getColors();

            const theme =
                document.documentElement.getAttribute('data-theme');

            drawGrid(theme);

            const edgeColor =
                theme === 'dark'
                    ? 'rgba(255,255,255,0.055)'
                    : 'rgba(0,0,0,0.06)';

            const edgeHighlight =
                theme === 'dark'
                    ? 'rgba(0,229,255,0.38)'
                    : 'rgba(8,145,178,0.28)';

            // ── HOVER DETECTION ────────────────────
            let hoveredNode = -1;

            nodes.forEach((node, i) => {

                const dx = mouseX - node.x;
                const dy = mouseY - node.y;

                if (
                    Math.sqrt(dx * dx + dy * dy) <
                    node.size + 12
                ) {
                    hoveredNode = i;
                }

                node.highlighted = false;
            });

            if (hoveredNode >= 0) {

                nodes[hoveredNode].highlighted = true;

                edges.forEach(edge => {

                    if (edge.from === hoveredNode) {
                        nodes[edge.to].highlighted = true;
                    }

                    if (edge.to === hoveredNode) {
                        nodes[edge.from].highlighted = true;
                    }
                });
            }

            // ── DRAW EDGES ─────────────────────────
            edges.forEach(edge => {

                const from = nodes[edge.from];
                const to = nodes[edge.to];

                const highlighted =
                    from.highlighted && to.highlighted;

                ctx.strokeStyle =
                    highlighted
                        ? edgeHighlight
                        : edgeColor;

                ctx.lineWidth =
                    highlighted
                        ? 1.4
                        : 0.6;

                ctx.beginPath();
                ctx.moveTo(from.x, from.y);
                ctx.lineTo(to.x, to.y);
                ctx.stroke();
            });

            // ── DRAW NODES ─────────────────────────
            nodes.forEach(node => {

                const color = typeColors[node.type];

                const alpha =
                    node.highlighted
                        ? 1
                        : 0.72;

                const radius =
                    node.highlighted
                        ? node.size * 1.45
                        : node.size;

                // Glow
                if (node.highlighted) {
                    ctx.shadowColor = color;
                    ctx.shadowBlur = 22;
                }

                ctx.globalAlpha = alpha;

                // Outer halo
                ctx.beginPath();
                ctx.fillStyle = `${color}18`;
                ctx.arc(
                    node.x,
                    node.y,
                    radius * 2.2,
                    0,
                    Math.PI * 2
                );
                ctx.fill();

                // Core node
                ctx.beginPath();
                ctx.fillStyle = color;
                ctx.arc(
                    node.x,
                    node.y,
                    radius,
                    0,
                    Math.PI * 2
                );
                ctx.fill();

                ctx.shadowBlur = 0;
                ctx.globalAlpha = 1;

                // ── LABELS ─────────────────────────
                if (node.highlighted) {

                    ctx.fillStyle =
                        theme === 'dark'
                            ? '#e2e8f0'
                            : '#0f172a';

                    ctx.font =
                        '11px "JetBrains Mono", monospace';

                    ctx.textAlign = 'center';

                    ctx.fillText(
                        node.label,
                        node.x,
                        node.y - radius - 10
                    );
                }

                // ── MOVEMENT ───────────────────────
                node.x += node.vx;
                node.y += node.vy;

                if (
                    node.x < 24 ||
                    node.x > canvas.width - 24
                ) {
                    node.vx *= -1;
                }

                if (
                    node.y < 24 ||
                    node.y > canvas.height - 24
                ) {
                    node.vy *= -1;
                }
            });

            animationFrame =
                requestAnimationFrame(drawNetwork);
        }

        // ── MOUSE INTERACTION ─────────────────────
        canvas.addEventListener(
            'mousemove',
            e => {

                const r =
                    canvas.getBoundingClientRect();

                mouseX = e.clientX - r.left;
                mouseY = e.clientY - r.top;
            },
            { passive: true }
        );

        canvas.addEventListener(
            'mouseleave',
            () => {
                mouseX = -1000;
                mouseY = -1000;
            }
        );

        // ── INITIALIZE ────────────────────────────
        resizeCanvas();
        createNetwork();
        drawNetwork();

        // ── HANDLE RESIZE ─────────────────────────
        window.addEventListener('resize', () => {

            cancelAnimationFrame(animationFrame);

            resizeCanvas();
            createNetwork();
            drawNetwork();
        });
    }
    // ── SCROLL INDICATOR ──────────────────────────
    function initScrollIndicator() {
        const indicator = document.getElementById('scrollIndicator');
        if (!indicator) return;

        window.addEventListener('scroll', () => {
            if (window.scrollY > 100) {
                indicator.style.opacity = '0';
                indicator.style.pointerEvents = 'none';
            } else {
                indicator.style.opacity = '1';
                indicator.style.pointerEvents = 'auto';
            }
        }, { passive: true });
    }

    // ── DOWNLOAD CV ───────────────────────────────
    function initDownloadCV() {
        const btn = document.getElementById('downloadCV');
        if (!btn) return;

        btn.addEventListener('click', (e) => {
            e.preventDefault();
            showToast('CV download will be available soon');
        });
    }

    // ── INIT ──────────────────────────────────────
    document.addEventListener('DOMContentLoaded', () => {
        initTheme();
        initNav();
        initScrollSpy();
        initReveal();
        initSmoothScroll();
        initPhoto();
        initCursorGlow();
        initCounters();
        initFilters();
        initExpand();
        initCopy();
        initForm();
        initNetworkBackground();
        initNetworkViz();
        initScrollIndicator();
        initDownloadCV();
    });

})();
