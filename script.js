/* ════════════════════════════════════════════════════════════
   EIGENSPACE — Complete JavaScript v2
   Biological + Mathematical Visualizations
   ════════════════════════════════════════════════════════════ */

(function () {
    'use strict';

    const $ = (sel) => document.querySelector(sel);
    const $$ = (sel) => document.querySelectorAll(sel);
    const lerp = (a, b, t) => a + (b - a) * t;
    const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
    const TAU = Math.PI * 2;

    /* ── JACOBI EIGENVALUE ──────────────────────────────────── */
    function jacobiEigen(matrix) {
        const n = matrix.length;
        const A = matrix.map(r => [...r]);
        const V = Array.from({ length: n }, (_, i) =>
            Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))
        );

        for (let iter = 0; iter < 100; iter++) {
            let maxVal = 0, p = 0, q = 1;
            for (let i = 0; i < n; i++) {
                for (let j = i + 1; j < n; j++) {
                    if (Math.abs(A[i][j]) > maxVal) {
                        maxVal = Math.abs(A[i][j]);
                        p = i; q = j;
                    }
                }
            }
            if (maxVal < 1e-10) break;

            const theta = 0.5 * Math.atan2(2 * A[p][q], A[p][p] - A[q][q]);
            const c = Math.cos(theta), s = Math.sin(theta);

            const App = A[p][p], Aqq = A[q][q], Apq = A[p][q];
            A[p][p] = c * c * App + 2 * s * c * Apq + s * s * Aqq;
            A[q][q] = s * s * App - 2 * s * c * Apq + c * c * Aqq;
            A[p][q] = 0; A[q][p] = 0;

            for (let i = 0; i < n; i++) {
                if (i !== p && i !== q) {
                    const Aip = A[i][p], Aiq = A[i][q];
                    A[i][p] = c * Aip + s * Aiq;
                    A[p][i] = A[i][p];
                    A[i][q] = -s * Aip + c * Aiq;
                    A[q][i] = A[i][q];
                }
            }

            for (let i = 0; i < n; i++) {
                const Vip = V[i][p], Viq = V[i][q];
                V[i][p] = c * Vip + s * Viq;
                V[i][q] = -s * Vip + c * Viq;
            }
        }

        const eigenvalues = A.map((row, i) => row[i]);
        const indices = eigenvalues.map((_, i) => i);
        indices.sort((a, b) => eigenvalues[a] - eigenvalues[b]);

        return {
            values: indices.map(i => eigenvalues[i]),
            vectors: indices.map(i => V.map(row => row[i]))
        };
    }

    /* ── GRAPH LAPLACIAN ────────────────────────────────────── */
    function graphLaplacian(adj) {
        const n = adj.length;
        const L = Array.from({ length: n }, () => new Array(n).fill(0));
        for (let i = 0; i < n; i++) {
            let deg = 0;
            for (let j = 0; j < n; j++) {
                if (adj[i][j]) {
                    L[i][j] = -adj[i][j];
                    deg += adj[i][j];
                }
            }
            L[i][i] = deg;
        }
        return L;
    }

    /* ── LOADING SCREEN ─────────────────────────────────────── */
    function initLoader() {
        const loader = $('#loader');
        const bar = $('#loaderBar');
        const phase = $('#loaderPhase');
        const ids = ['loaderL', 'loaderD', 'loaderA'];
        const ops = document.querySelectorAll('.loader-op');

        const phases = [
            'initializing membrane potential',
            'constructing adjacency matrix',
            'computing graph Laplacian',
            'extracting eigenvectors',
            'embedding spectral coordinates',
            'system ready'
        ];

        const startTime = Date.now();

        // Animate loader cell
        const loaderCanvas = $('#loaderCanvas');
        if (loaderCanvas) {
            const lctx = loaderCanvas.getContext('2d');
            const dpr = window.devicePixelRatio || 1;
            loaderCanvas.width = 200 * dpr;
            loaderCanvas.height = 200 * dpr;
            lctx.scale(dpr, dpr);

            function drawLoaderCell(time) {
                if (loader.classList.contains('hidden')) return;
                lctx.clearRect(0, 0, 200, 200);
                const cx = 100, cy = 100;

                // Cell membrane
                lctx.beginPath();
                for (let i = 0; i <= 64; i++) {
                    const angle = (i / 64) * TAU;
                    const wobble = Math.sin(angle * 5 + time * 0.003) * 3
                        + Math.sin(angle * 3 - time * 0.002) * 2;
                    const r = 55 + wobble;
                    const x = cx + Math.cos(angle) * r;
                    const y = cy + Math.sin(angle) * r;
                    if (i === 0) lctx.moveTo(x, y);
                    else lctx.lineTo(x, y);
                }
                lctx.closePath();
                lctx.strokeStyle = 'rgba(0, 240, 255, 0.25)';
                lctx.lineWidth = 1.5;
                lctx.stroke();

                // Nucleus
                lctx.beginPath();
                for (let i = 0; i <= 32; i++) {
                    const angle = (i / 32) * TAU;
                    const wobble = Math.sin(angle * 3 + time * 0.004) * 2;
                    const r = 22 + wobble;
                    const x = cx + Math.cos(angle) * r;
                    const y = cy + Math.sin(angle) * r;
                    if (i === 0) lctx.moveTo(x, y);
                    else lctx.lineTo(x, y);
                }
                lctx.closePath();
                lctx.strokeStyle = 'rgba(255, 0, 170, 0.2)';
                lctx.lineWidth = 1;
                lctx.stroke();

                // DNA bits flowing from nucleus
                for (let i = 0; i < 8; i++) {
                    const t = ((time * 0.001 + i * 0.4) % 3) / 3;
                    if (t > 1) continue;
                    const angle = (i / 8) * TAU + time * 0.0005;
                    const r = 22 + t * 38;
                    const x = cx + Math.cos(angle) * r;
                    const y = cy + Math.sin(angle) * r;
                    const alpha = Math.sin(t * Math.PI) * 0.6;

                    // Binary bit
                    lctx.font = '7px "JetBrains Mono", monospace';
                    lctx.fillStyle = `rgba(0, 240, 255, ${alpha})`;
                    lctx.textAlign = 'center';
                    lctx.fillText(Math.random() > 0.5 ? '1' : '0', x, y + 3);
                }

                // Floating gene labels
                const genes = ['TP53', 'MYC', 'BRCA1', 'KRAS'];
                for (let i = 0; i < genes.length; i++) {
                    const t = ((time * 0.0008 + i * 0.7) % 4) / 4;
                    const angle = (i / genes.length) * TAU + time * 0.0003;
                    const r = 40 + Math.sin(t * Math.PI) * 20;
                    const x = cx + Math.cos(angle) * r;
                    const y = cy + Math.sin(angle) * r;
                    const alpha = Math.sin(t * Math.PI) * 0.3;

                    lctx.font = '6px "JetBrains Mono", monospace';
                    lctx.fillStyle = `rgba(0, 240, 255, ${alpha})`;
                    lctx.textAlign = 'center';
                    lctx.fillText(genes[i], x, y);
                }

                requestAnimationFrame(drawLoaderCell);
            }
            requestAnimationFrame(drawLoaderCell);
        }

        ids.forEach((id, i) => {
            setTimeout(() => {
                const el = document.getElementById(id);
                if (el) el.classList.add('visible');
            }, 300 + i * 350);
        });

        ops.forEach((op, i) => {
            setTimeout(() => op.classList.add('visible'), 400 + i * 350);
        });

        const interval = setInterval(() => {
            const elapsed = Date.now() - startTime;
            const progress = Math.min(100, (elapsed / 2800) * 100);
            bar.style.width = progress + '%';

            const phaseIdx = Math.min(phases.length - 1, Math.floor((progress / 100) * phases.length));
            if (phase) phase.textContent = phases[phaseIdx];

            if (progress >= 100) {
                clearInterval(interval);
                setTimeout(() => {
                    loader.classList.add('hidden');
                    const nav = $('.spec-nav');
                    const chromo = $('.chromo-scroll');
                    if (nav) nav.classList.add('visible');
                    if (chromo) chromo.classList.add('visible');
                    revealHero();
                    startGeneStream();
                }, 400);
            }
        }, 16);
    }

    function revealHero() {
        $$('.name-line').forEach((line, i) => {
            setTimeout(() => line.classList.add('revealed'), 200 + i * 300);
        });
    }

    /* ── GENE STREAM (floating bits) ────────────────────────── */
    function startGeneStream() {
        const container = $('#geneStream');
        if (!container) return;

        const genes = ['ATCG', '1010', 'TP53', '0110', 'BRCA', '1001', 'MYC·', '0011',
            'KRAS', '1100', 'PTEN', '0101', 'RB1·', '1110', 'AKT1'];

        function spawnBit() {
            const bit = document.createElement('div');
            bit.className = 'gene-bit';
            bit.textContent = genes[Math.floor(Math.random() * genes.length)];
            bit.style.left = (Math.random() * 100) + '%';
            bit.style.animationDuration = (10 + Math.random() * 8) + 's';
            bit.style.animationDelay = (Math.random() * 2) + 's';
            container.appendChild(bit);
            setTimeout(() => bit.remove(), 20000);
        }

        setInterval(spawnBit, 2500);
        // Initial batch
        for (let i = 0; i < 4; i++) setTimeout(spawnBit, i * 600);
    }

    /* ── BACKGROUND NETWORK ─────────────────────────────────── */
    function initNetwork() {
        const canvas = $('#networkCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        let width, height, nodes = [];
        const mouse = { x: -1000, y: -1000 };
        const NODE_COUNT = 55;
        const CONNECT_DIST = 140;

        function resize() {
            width = canvas.width = window.innerWidth;
            height = canvas.height = window.innerHeight;
        }

        function createNodes() {
            nodes = [];
            for (let i = 0; i < NODE_COUNT; i++) {
                nodes.push({
                    x: Math.random() * width,
                    y: Math.random() * height,
                    vx: (Math.random() - 0.5) * 0.2,
                    vy: (Math.random() - 0.5) * 0.2,
                    r: Math.random() * 1.5 + 0.5,
                    phase: Math.random() * TAU
                });
            }
        }

        function draw(time) {
            ctx.clearRect(0, 0, width, height);

            for (const node of nodes) {
                node.x += node.vx + Math.sin(time * 0.0003 + node.phase) * 0.06;
                node.y += node.vy + Math.cos(time * 0.0004 + node.phase) * 0.04;

                const dx = node.x - mouse.x;
                const dy = node.y - mouse.y;
                const dist = Math.sqrt(dx * dx + dy * dy);
                if (dist < 120 && dist > 0) {
                    node.x += (dx / dist) * (120 - dist) * 0.004;
                    node.y += (dy / dist) * (120 - dist) * 0.004;
                }

                if (node.x < -20) node.x = width + 20;
                if (node.x > width + 20) node.x = -20;
                if (node.y < -20) node.y = height + 20;
                if (node.y > height + 20) node.y = -20;
            }

            for (let i = 0; i < nodes.length; i++) {
                for (let j = i + 1; j < nodes.length; j++) {
                    const dx = nodes[i].x - nodes[j].x;
                    const dy = nodes[i].y - nodes[j].y;
                    const dist = Math.sqrt(dx * dx + dy * dy);
                    if (dist < CONNECT_DIST) {
                        ctx.beginPath();
                        ctx.moveTo(nodes[i].x, nodes[i].y);
                        ctx.lineTo(nodes[j].x, nodes[j].y);
                        ctx.strokeStyle = `rgba(0, 240, 255, ${(1 - dist / CONNECT_DIST) * 0.07})`;
                        ctx.lineWidth = 0.5;
                        ctx.stroke();
                    }
                }
            }

            for (const node of nodes) {
                const pulse = Math.sin(time * 0.002 + node.phase) * 0.3 + 0.7;
                ctx.beginPath();
                ctx.arc(node.x, node.y, node.r, 0, TAU);
                ctx.fillStyle = `rgba(0, 240, 255, ${0.2 * pulse})`;
                ctx.fill();
            }

            requestAnimationFrame(draw);
        }

        resize();
        createNodes();
        window.addEventListener('resize', () => { resize(); createNodes(); });
        window.addEventListener('mousemove', (e) => { mouse.x = e.clientX; mouse.y = e.clientY; });
        requestAnimationFrame(draw);
    }

    /* ── OVERLINE MINI HELIX ────────────────────────────────── */
    function initOverlineHelix() {
        const canvas = $('#overlineHelix');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        canvas.width = 24 * dpr;
        canvas.height = 24 * dpr;
        ctx.scale(dpr, dpr);

        function draw(time) {
            ctx.clearRect(0, 0, 24, 24);
            for (let i = 0; i < 12; i++) {
                const t = i / 12;
                const y = t * 24;
                const phase = time * 0.003 + t * Math.PI * 3;
                const x1 = 12 + Math.sin(phase) * 6;
                const x2 = 12 + Math.sin(phase + Math.PI) * 6;
                const d1 = Math.cos(phase);

                ctx.beginPath();
                ctx.arc(d1 > 0 ? x1 : x2, y, 1.2, 0, TAU);
                ctx.fillStyle = d1 > 0 ? 'rgba(0, 240, 255, 0.6)' : 'rgba(255, 0, 170, 0.6)';
                ctx.fill();
            }
            requestAnimationFrame(draw);
        }
        requestAnimationFrame(draw);
    }

    /* ── AVATAR RING ────────────────────────────────────────── */
    function initAvatarRing() {
        const canvas = $('#avatarRingCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        canvas.width = 96 * dpr;
        canvas.height = 96 * dpr;
        ctx.scale(dpr, dpr);

        function draw(time) {
            ctx.clearRect(0, 0, 96, 96);
            const cx = 48, cy = 48, r = 42;

            for (let i = 0; i < 64; i++) {
                const angle = (i / 64) * TAU;
                const pulse = Math.sin(time * 0.002 + angle * 3) * 0.3 + 0.7;
                const x = cx + Math.cos(angle) * r;
                const y = cy + Math.sin(angle) * r;

                const t = (angle / TAU);
                const cr = Math.round(lerp(0, 255, t));
                const cg = Math.round(lerp(240, 0, t));
                const cb = 255;

                ctx.beginPath();
                ctx.arc(x, y, 1.5, 0, TAU);
                ctx.fillStyle = `rgba(${cr}, ${cg}, ${cb}, ${0.3 * pulse})`;
                ctx.fill();
            }
            requestAnimationFrame(draw);
        }
        requestAnimationFrame(draw);
    }

    /* ── EIGENDECOMPOSITION VISUALIZATION ───────────────────── */
    function initEigenViz() {
        const canvas = $('#eigenCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;

        const W = 560, H = 520;
        canvas.width = W * dpr;
        canvas.height = H * dpr;
        ctx.scale(dpr, dpr);

        const N = 8;
        let adjacency = [
            [0, 1, 1, 0, 0, 0, 0, 1],
            [1, 0, 1, 1, 0, 0, 0, 0],
            [1, 1, 0, 1, 1, 0, 0, 0],
            [0, 1, 1, 0, 1, 0, 0, 0],
            [0, 0, 1, 1, 0, 1, 1, 0],
            [0, 0, 0, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 1, 1, 0, 1],
            [1, 0, 0, 0, 0, 1, 1, 0]
        ];

        const geneLabels = ['TP53', 'BRCA1', 'MYC', 'EGFR', 'KRAS', 'PTEN', 'RB1', 'AKT1'];

        const cx = 200, cy = 190, radius = 120;
        const nodePos = [];
        for (let i = 0; i < N; i++) {
            const angle = -Math.PI / 2 + (TAU / N) * i;
            nodePos.push({
                x: cx + Math.cos(angle) * radius,
                y: cy + Math.sin(angle) * radius
            });
        }

        const baseColors = [
            '#00f0ff', '#22d3ee', '#06b6d4', '#0891b2',
            '#ff6090', '#ec4899', '#a855f7', '#ffb000'
        ];

        let eigen = jacobiEigen(adjacency);
        let displayValues = eigen.values.map(() => 0);
        let targetValues = eigen.values.slice();
        let hoveredBar = -1;
        let animPhase = 0;
        const startTime = Date.now();

        setInterval(() => {
            let i = Math.floor(Math.random() * N);
            let j = Math.floor(Math.random() * N);
            if (i === j) return;

            let edgeCount = 0;
            for (let a = 0; a < N; a++)
                for (let b = a + 1; b < N; b++)
                    if (adjacency[a][b]) edgeCount++;

            if (adjacency[i][j] && edgeCount <= 8) return;

            adjacency[i][j] = adjacency[i][j] ? 0 : 1;
            adjacency[j][i] = adjacency[i][j];

            eigen = jacobiEigen(adjacency);
            targetValues = eigen.values.slice();
        }, 4500);

        canvas.addEventListener('mousemove', (e) => {
            const r = canvas.getBoundingClientRect();
            const mx = (e.clientX - r.left) * (W / r.width);
            const my = (e.clientY - r.top) * (H / r.height);

            hoveredBar = -1;
            const barSp = 55;
            const barStartX = (W - barSp * N) / 2 + 12;
            for (let i = 0; i < N; i++) {
                const bx = barStartX + i * barSp;
                if (mx >= bx && mx <= bx + 40 && my >= 300 && my <= 480) {
                    hoveredBar = i;
                    break;
                }
            }
        });

        canvas.addEventListener('mouseleave', () => { hoveredBar = -1; });

        let isVisible = true;
        const observer = new IntersectionObserver(([e]) => { isVisible = e.isIntersecting; }, { threshold: 0.05 });
        observer.observe(canvas);

        function draw() {
            if (!isVisible) { requestAnimationFrame(draw); return; }

            const elapsed = Date.now() - startTime;
            animPhase += 0.016;
            ctx.clearRect(0, 0, W, H);

            const nodesReady = elapsed > 600;
            const edgesReady = elapsed > 1100;
            const barsReady = elapsed > 1600;

            for (let i = 0; i < displayValues.length; i++) {
                displayValues[i] = lerp(displayValues[i], targetValues[i], 0.04);
            }

            // Edges
            if (edgesReady) {
                for (let i = 0; i < N; i++) {
                    for (let j = i + 1; j < N; j++) {
                        if (!adjacency[i][j]) continue;
                        let alpha = 0.18;
                        let lw = 1;
                        let color = '100, 140, 255';

                        if (hoveredBar >= 0) {
                            const vi = eigen.vectors[hoveredBar][i];
                            const vj = eigen.vectors[hoveredBar][j];
                            color = (vi * vj > 0) ? '0, 240, 255' : '255, 0, 170';
                            alpha = 0.5;
                            lw = 1.8;
                        }

                        const pulse = Math.sin(animPhase * 1.2 + i * 0.4) * 0.06;

                        // Draw animated signal along edge
                        ctx.beginPath();
                        ctx.moveTo(nodePos[i].x, nodePos[i].y);
                        ctx.lineTo(nodePos[j].x, nodePos[j].y);
                        ctx.strokeStyle = `rgba(${color}, ${alpha + pulse})`;
                        ctx.lineWidth = lw;
                        ctx.stroke();

                        // Signal particle moving along edge
                        if (hoveredBar < 0) {
                            const t = ((animPhase * 0.5 + i * 0.3 + j * 0.7) % 2) / 2;
                            if (t < 1) {
                                const sx = lerp(nodePos[i].x, nodePos[j].x, t);
                                const sy = lerp(nodePos[i].y, nodePos[j].y, t);
                                ctx.beginPath();
                                ctx.arc(sx, sy, 1.5, 0, TAU);
                                ctx.fillStyle = `rgba(0, 240, 255, ${Math.sin(t * Math.PI) * 0.4})`;
                                ctx.fill();
                            }
                        }
                    }
                }
            }

            // Nodes
            if (nodesReady) {
                for (let i = 0; i < N; i++) {
                    const pos = nodePos[i];
                    const revealAlpha = clamp((elapsed - 600 - i * 70) / 400, 0, 1);
                    const scale = revealAlpha;

                    let fillColor = baseColors[i];
                    let nodeR = 8;
                    let glowR = 20;

                    if (hoveredBar >= 0) {
                        const v = eigen.vectors[hoveredBar][i];
                        const absV = Math.abs(v);
                        nodeR = 5 + absV * 16;
                        glowR = 14 + absV * 24;
                        fillColor = v > 0.01 ? '#00f0ff' : v < -0.01 ? '#ff00aa' : '#555';
                    }

                    // Glow
                    const grd = ctx.createRadialGradient(pos.x, pos.y, 0, pos.x, pos.y, glowR * scale);
                    grd.addColorStop(0, fillColor + '33');
                    grd.addColorStop(1, fillColor + '00');
                    ctx.beginPath();
                    ctx.arc(pos.x, pos.y, glowR * scale, 0, TAU);
                    ctx.fillStyle = grd;
                    ctx.fill();

                    // Circle
                    ctx.beginPath();
                    ctx.arc(pos.x, pos.y, nodeR * scale, 0, TAU);
                    ctx.fillStyle = fillColor;
                    ctx.globalAlpha = revealAlpha;
                    ctx.fill();
                    ctx.globalAlpha = 1;

                    // Gene label
                    if (revealAlpha > 0.5) {
                        ctx.font = '9px "JetBrains Mono", monospace';
                        ctx.fillStyle = 'rgba(226, 232, 240, 0.65)';
                        ctx.textAlign = 'center';
                        ctx.fillText(geneLabels[i], pos.x, pos.y - nodeR * scale - 10);
                    }
                }
            }

            // Eigenvalue bars
            if (barsReady) {
                const barY = 400;
                const maxBarH = 90;
                const barW = 40;
                const barSp = 55;
                const barStartX = (W - barSp * N) / 2 + 12;
                const maxAbs = Math.max(...displayValues.map(Math.abs), 0.01);

                ctx.font = '9px "JetBrains Mono", monospace';
                ctx.fillStyle = 'rgba(100, 116, 139, 0.4)';
                ctx.textAlign = 'left';
                ctx.fillText('EIGENVALUE SPECTRUM', barStartX, barY - maxBarH - 18);

                // Zero line
                ctx.beginPath();
                ctx.moveTo(barStartX - 4, barY);
                ctx.lineTo(barStartX + barSp * N, barY);
                ctx.strokeStyle = 'rgba(100, 140, 255, 0.1)';
                ctx.lineWidth = 0.5;
                ctx.stroke();

                for (let i = 0; i < N; i++) {
                    const bx = barStartX + i * barSp;
                    const val = displayValues[i];
                    const barH = (val / maxAbs) * maxBarH;
                    const rt = clamp((elapsed - 1600 - i * 80) / 500, 0, 1);
                    const h = barH * rt;
                    const hovered = hoveredBar === i;

                    const t = (val / maxAbs + 1) / 2;
                    const r = Math.round(lerp(255, 0, t));
                    const g = Math.round(lerp(60, 240, t));
                    const b = Math.round(lerp(150, 255, t));

                    if (hovered) {
                        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, 0.08)`;
                        ctx.fillRect(bx - 4, Math.min(barY, barY - h) - 4, barW + 8, Math.abs(h) + 8);
                    }

                    // Bar fill with gradient
                    const barGrd = ctx.createLinearGradient(bx, barY, bx, barY - h);
                    barGrd.addColorStop(0, `rgba(${r}, ${g}, ${b}, ${hovered ? 0.8 : 0.45})`);
                    barGrd.addColorStop(1, `rgba(${r}, ${g}, ${b}, ${hovered ? 0.95 : 0.6})`);
                    ctx.fillStyle = barGrd;
                    ctx.fillRect(bx, barY - h, barW, h);

                    ctx.strokeStyle = hovered ? `rgb(${r}, ${g}, ${b})` : `rgba(${r}, ${g}, ${b}, 0.2)`;
                    ctx.lineWidth = hovered ? 1.5 : 0.5;
                    ctx.strokeRect(bx, barY - h, barW, h);

                    ctx.font = hovered ? 'bold 10px "JetBrains Mono", monospace' : '8px "JetBrains Mono", monospace';
                    ctx.fillStyle = hovered ? `rgb(${r}, ${g}, ${b})` : 'rgba(148, 163, 184, 0.5)';
                    ctx.textAlign = 'center';
                    ctx.fillText(val.toFixed(2), bx + barW / 2, barY + 16);

                    ctx.font = '10px "EB Garamond", serif';
                    ctx.fillText('λ' + String.fromCharCode(8321 + i), bx + barW / 2, barY + 30);
                }
            }

            // Title
            ctx.font = '9px "JetBrains Mono", monospace';
            ctx.fillStyle = 'rgba(100, 116, 139, 0.3)';
            ctx.textAlign = 'left';
            ctx.fillText('GENE INTERACTION NETWORK', 20, 22);

            // Hover hint
            if (hoveredBar < 0 && elapsed > 3000) {
                ctx.font = '8px "JetBrains Mono", monospace';
                ctx.fillStyle = `rgba(100, 116, 139, ${0.15 + Math.sin(animPhase) * 0.1})`;
                ctx.textAlign = 'center';
                ctx.fillText('hover eigenvalues to see eigenvector modes', W / 2, H - 16);
            }

            requestAnimationFrame(draw);
        }

        requestAnimationFrame(draw);
    }

    /* ── CELL SIGNALING VISUALIZATION ───────────────────────── */
    function initCellViz() {
        const canvas = $('#cellCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        const W = 280, H = 320;
        canvas.width = W * dpr;
        canvas.height = H * dpr;
        ctx.scale(dpr, dpr);

        let isVisible = false;
        const observer = new IntersectionObserver(([e]) => { isVisible = e.isIntersecting; }, { threshold: 0.1 });
        observer.observe(canvas);

        const particles = [];
        const genes = ['TP53', 'MYC', 'EGFR', 'KRAS', 'PTEN', 'AKT1', 'RB1', 'BRCA'];

        // Generate signaling particles
        for (let i = 0; i < 30; i++) {
            particles.push({
                x: W / 2 + (Math.random() - 0.5) * 60,
                y: H / 2 + (Math.random() - 0.5) * 60,
                vx: (Math.random() - 0.5) * 0.8,
                vy: (Math.random() - 0.5) * 0.8,
                type: Math.random() > 0.6 ? 'gene' : 'signal',
                label: genes[Math.floor(Math.random() * genes.length)],
                bits: Array.from({ length: 4 }, () => Math.random() > 0.5 ? '1' : '0').join(''),
                phase: Math.random() * TAU,
                r: Math.random() * 2 + 1
            });
        }

        function draw(time) {
            if (!isVisible) { requestAnimationFrame(draw); return; }
            ctx.clearRect(0, 0, W, H);

            const cx = W / 2, cy = H / 2;

            // Cell membrane (outer)
            ctx.beginPath();
            for (let i = 0; i <= 80; i++) {
                const angle = (i / 80) * TAU;
                const wobble = Math.sin(angle * 6 + time * 0.002) * 4
                    + Math.sin(angle * 4 - time * 0.003) * 3;
                const r = 110 + wobble;
                const x = cx + Math.cos(angle) * r;
                const y = cy + Math.sin(angle) * r;
                if (i === 0) ctx.moveTo(x, y);
                else ctx.lineTo(x, y);
            }
            ctx.closePath();
            const memGrd = ctx.createRadialGradient(cx, cy, 60, cx, cy, 120);
            memGrd.addColorStop(0, 'rgba(0, 240, 255, 0.02)');
            memGrd.addColorStop(1, 'rgba(0, 240, 255, 0.06)');
            ctx.fillStyle = memGrd;
            ctx.fill();
            ctx.strokeStyle = 'rgba(0, 240, 255, 0.2)';
            ctx.lineWidth = 1.2;
            ctx.stroke();

            // Phospholipid bilayer dots
            for (let i = 0; i < 48; i++) {
                const angle = (i / 48) * TAU;
                const wobble = Math.sin(angle * 6 + time * 0.002) * 4
                    + Math.sin(angle * 4 - time * 0.003) * 3;
                const r = 110 + wobble;
                const x = cx + Math.cos(angle) * r;
                const y = cy + Math.sin(angle) * r;
                ctx.beginPath();
                ctx.arc(x, y, 2, 0, TAU);
                ctx.fillStyle = 'rgba(0, 240, 255, 0.15)';
                ctx.fill();
            }

            // Nucleus
            ctx.beginPath();
            for (let i = 0; i <= 40; i++) {
                const angle = (i / 40) * TAU;
                const wobble = Math.sin(angle * 3 + time * 0.004) * 2;
                const r = 35 + wobble;
                const x = cx + Math.cos(angle) * r;
                const y = cy + Math.sin(angle) * r;
                if (i === 0) ctx.moveTo(x, y);
                else ctx.lineTo(x, y);
            }
            ctx.closePath();
            ctx.strokeStyle = 'rgba(255, 0, 170, 0.25)';
            ctx.lineWidth = 1;
            ctx.stroke();
            const nucGrd = ctx.createRadialGradient(cx, cy, 0, cx, cy, 35);
            nucGrd.addColorStop(0, 'rgba(255, 0, 170, 0.06)');
            nucGrd.addColorStop(1, 'rgba(255, 0, 170, 0.01)');
            ctx.fillStyle = nucGrd;
            ctx.fill();

            // DNA double helix inside nucleus
            for (let i = 0; i < 20; i++) {
                const t = i / 20;
                const y = cy - 25 + t * 50;
                const phase = time * 0.003 + t * Math.PI * 3;
                const x1 = cx + Math.sin(phase) * 12;
                const x2 = cx + Math.sin(phase + Math.PI) * 12;
                const d = Math.cos(phase);

                if (d > 0) {
                    ctx.beginPath();
                    ctx.arc(x1, y, 1.5, 0, TAU);
                    ctx.fillStyle = 'rgba(0, 240, 255, 0.4)';
                    ctx.fill();
                } else {
                    ctx.beginPath();
                    ctx.arc(x2, y, 1.5, 0, TAU);
                    ctx.fillStyle = 'rgba(255, 0, 170, 0.4)';
                    ctx.fill();
                }

                // Base pair connections
                if (i % 3 === 0) {
                    ctx.beginPath();
                    ctx.moveTo(x1, y);
                    ctx.lineTo(x2, y);
                    ctx.strokeStyle = 'rgba(255, 176, 0, 0.1)';
                    ctx.lineWidth = 0.5;
                    ctx.stroke();
                }
            }

            // Particles (genes + signals)
            for (const p of particles) {
                // Move
                p.x += p.vx;
                p.y += p.vy;

                // Keep inside cell
                const dx = p.x - cx;
                const dy = p.y - cy;
                const dist = Math.sqrt(dx * dx + dy * dy);
                if (dist > 95) {
                    p.vx -= (dx / dist) * 0.1;
                    p.vy -= (dy / dist) * 0.1;
                }

                // Brownian
                p.vx += (Math.random() - 0.5) * 0.1;
                p.vy += (Math.random() - 0.5) * 0.1;
                p.vx *= 0.98;
                p.vy *= 0.98;

                if (p.type === 'gene') {
                    // Gene with binary representation
                    const alpha = 0.5 + Math.sin(time * 0.002 + p.phase) * 0.2;
                    ctx.font = '7px "JetBrains Mono", monospace';
                    ctx.fillStyle = `rgba(0, 240, 255, ${alpha * 0.7})`;
                    ctx.textAlign = 'center';
                    ctx.fillText(p.label, p.x, p.y - 5);

                    ctx.font = '5px "JetBrains Mono", monospace';
                    ctx.fillStyle = `rgba(0, 240, 255, ${alpha * 0.35})`;
                    ctx.fillText(p.bits, p.x, p.y + 5);
                } else {
                    // Signal dot
                    const pulse = Math.sin(time * 0.003 + p.phase) * 0.3 + 0.6;
                    ctx.beginPath();
                    ctx.arc(p.x, p.y, p.r, 0, TAU);
                    ctx.fillStyle = `rgba(255, 176, 0, ${0.3 * pulse})`;
                    ctx.fill();
                }
            }

            // Receptor proteins on membrane
            for (let i = 0; i < 6; i++) {
                const angle = (i / 6) * TAU + time * 0.0003;
                const wobble = Math.sin(angle * 6 + time * 0.002) * 4;
                const r = 110 + wobble;
                const x = cx + Math.cos(angle) * r;
                const y = cy + Math.sin(angle) * r;

                ctx.beginPath();
                ctx.arc(x, y, 4, 0, TAU);
                ctx.strokeStyle = 'rgba(168, 85, 247, 0.35)';
                ctx.lineWidth = 1;
                ctx.stroke();

                // Signal arrow from receptor toward nucleus
                const progress = ((time * 0.001 + i) % 3) / 3;
                if (progress < 0.8) {
                    const sx = lerp(x, cx, progress);
                    const sy = lerp(y, cy, progress);
                    ctx.beginPath();
                    ctx.arc(sx, sy, 1.5, 0, TAU);
                    ctx.fillStyle = `rgba(168, 85, 247, ${Math.sin(progress * Math.PI) * 0.4})`;
                    ctx.fill();
                }
            }

            requestAnimationFrame(draw);
        }
        requestAnimationFrame(draw);
    }

    /* ── DNA HELIX ──────────────────────────────────────────── */
    function initHelix() {
        const canvas = $('#helixCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        const W = 80, H = 260;
        canvas.width = W * dpr;
        canvas.height = H * dpr;
        ctx.scale(dpr, dpr);

        let isVisible = false;
        const observer = new IntersectionObserver(([e]) => { isVisible = e.isIntersecting; }, { threshold: 0.1 });
        observer.observe(canvas);

        const basePairs = ['A-T', 'T-A', 'G-C', 'C-G', 'A-T', 'G-C', 'T-A', 'C-G'];

        function draw(time) {
            if (!isVisible) { requestAnimationFrame(draw); return; }
            ctx.clearRect(0, 0, W, H);

            const steps = 50;
            for (let i = 0; i < steps; i++) {
                const t = i / steps;
                const y = t * H;
                const phase = time * 0.001 + t * Math.PI * 5;
                const x1 = W / 2 + Math.sin(phase) * 22;
                const x2 = W / 2 + Math.sin(phase + Math.PI) * 22;
                const depth1 = Math.cos(phase);
                const depth2 = Math.cos(phase + Math.PI);

                // Base pair connections
                if (i % 5 === 0) {
                    ctx.beginPath();
                    ctx.moveTo(x1, y);
                    ctx.lineTo(x2, y);
                    ctx.strokeStyle = 'rgba(255, 176, 0, 0.12)';
                    ctx.lineWidth = 1;
                    ctx.stroke();

                    // Base pair label
                    const bpIdx = Math.floor(i / 5) % basePairs.length;
                    ctx.font = '5px "JetBrains Mono", monospace';
                    ctx.fillStyle = 'rgba(255, 176, 0, 0.15)';
                    ctx.textAlign = 'center';
                    ctx.fillText(basePairs[bpIdx], W / 2, y - 2);
                }

                // Strand 1 (behind)
                if (depth1 < 0) {
                    ctx.beginPath();
                    ctx.arc(x1, y, 2, 0, TAU);
                    ctx.fillStyle = `rgba(0, 240, 255, ${0.15 + depth1 * 0.1})`;
                    ctx.fill();
                }

                // Strand 2 (behind)
                if (depth2 < 0) {
                    ctx.beginPath();
                    ctx.arc(x2, y, 2, 0, TAU);
                    ctx.fillStyle = `rgba(255, 0, 170, ${0.15 + depth2 * 0.1})`;
                    ctx.fill();
                }

                // Strand 1 (front)
                if (depth1 >= 0) {
                    ctx.beginPath();
                    ctx.arc(x1, y, 2.8, 0, TAU);
                    ctx.fillStyle = `rgba(0, 240, 255, ${0.3 + depth1 * 0.35})`;
                    ctx.fill();
                }

                // Strand 2 (front)
                if (depth2 >= 0) {
                    ctx.beginPath();
                    ctx.arc(x2, y, 2.8, 0, TAU);
                    ctx.fillStyle = `rgba(255, 0, 170, ${0.3 + depth2 * 0.35})`;
                    ctx.fill();
                }
            }

            requestAnimationFrame(draw);
        }
        requestAnimationFrame(draw);
    }

    /* ── EQUATION ANIMATION ─────────────────────────────────── */
    function initEquation() {
        const canvas = $('#equationCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        const W = 400, H = 60;
        canvas.width = W * dpr;
        canvas.height = H * dpr;
        ctx.scale(dpr, dpr);

        let isVisible = false;
        const observer = new IntersectionObserver(([e]) => { isVisible = e.isIntersecting; }, { threshold: 0.1 });
        observer.observe(canvas);

        function draw(time) {
            if (!isVisible) { requestAnimationFrame(draw); return; }
            ctx.clearRect(0, 0, W, H);

            const pulse = Math.sin(time * 0.002) * 0.1 + 0.9;

            ctx.font = '16px "EB Garamond", serif';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';

            // A = UΣV^T  →  L = D - A  →  Lv = λv
            const eqs = [
                { text: 'A', x: 30, color: `rgba(226, 232, 240, ${0.7 * pulse})` },
                { text: '=', x: 50, color: `rgba(100, 116, 139, ${0.5 * pulse})` },
                { text: 'U', x: 70, color: `rgba(0, 240, 255, ${0.7 * pulse})` },
                { text: 'Σ', x: 90, color: `rgba(255, 0, 170, ${0.7 * pulse})` },
                { text: 'Vᵀ', x: 112, color: `rgba(255, 176, 0, ${0.7 * pulse})` },
                { text: '→', x: 150, color: `rgba(100, 116, 139, ${0.3 * pulse})` },
                { text: 'L', x: 185, color: `rgba(226, 232, 240, ${0.7 * pulse})` },
                { text: '=', x: 205, color: `rgba(100, 116, 139, ${0.5 * pulse})` },
                { text: 'D', x: 225, color: `rgba(0, 240, 255, ${0.7 * pulse})` },
                { text: '−', x: 245, color: `rgba(100, 116, 139, ${0.5 * pulse})` },
                { text: 'A', x: 265, color: `rgba(255, 0, 170, ${0.7 * pulse})` },
                { text: '→', x: 300, color: `rgba(100, 116, 139, ${0.3 * pulse})` },
                { text: 'Lv', x: 332, color: `rgba(226, 232, 240, ${0.7 * pulse})` },
                { text: '=', x: 355, color: `rgba(100, 116, 139, ${0.5 * pulse})` },
                { text: 'λv', x: 378, color: `rgba(0, 240, 255, ${0.8 * pulse})` },
            ];

            eqs.forEach((eq, i) => {
                const delay = i * 0.1;
                const alpha = clamp((time * 0.001 - delay) * 0.5, 0, 1);
                ctx.globalAlpha = alpha;
                ctx.fillStyle = eq.color;
                ctx.fillText(eq.text, eq.x, H / 2);
            });

            ctx.globalAlpha = 1;
            requestAnimationFrame(draw);
        }
        requestAnimationFrame(draw);
    }

    /* ── PROJECT MINI VISUALIZATIONS ────────────────────────── */
    function initProjectMiniViz() {
        const canvases = $$('.project-mini-viz');
        const dpr = window.devicePixelRatio || 1;

        canvases.forEach(canvas => {
            const ctx = canvas.getContext('2d');
            const type = canvas.dataset.type;
            canvas.width = 48 * dpr;
            canvas.height = 48 * dpr;
            ctx.scale(dpr, dpr);

            function draw(time) {
                ctx.clearRect(0, 0, 48, 48);

                if (type === 'network') {
                    const pts = [
                        { x: 12, y: 12 }, { x: 36, y: 10 }, { x: 24, y: 24 },
                        { x: 10, y: 38 }, { x: 38, y: 36 }
                    ];
                    const edges = [[0, 1], [0, 2], [1, 2], [2, 3], [2, 4], [3, 4]];

                    edges.forEach(([a, b]) => {
                        ctx.beginPath();
                        ctx.moveTo(pts[a].x, pts[a].y);
                        ctx.lineTo(pts[b].x, pts[b].y);
                        ctx.strokeStyle = 'rgba(0, 240, 255, 0.2)';
                        ctx.lineWidth = 0.8;
                        ctx.stroke();
                    });

                    pts.forEach((p, i) => {
                        const pulse = Math.sin(time * 0.003 + i) * 0.3 + 0.7;
                        ctx.beginPath();
                        ctx.arc(p.x, p.y, 3, 0, TAU);
                        ctx.fillStyle = `rgba(0, 240, 255, ${0.5 * pulse})`;
                        ctx.fill();
                    });

                } else if (type === 'eigen') {
                    const vals = [0.8, 0.5, 0.3, -0.2, -0.6];
                    vals.forEach((v, i) => {
                        const bx = 4 + i * 9;
                        const h = v * 16;
                        const pulse = Math.sin(time * 0.002 + i * 0.8) * 2;
                        ctx.fillStyle = v > 0 ? 'rgba(255, 0, 170, 0.5)' : 'rgba(0, 240, 255, 0.5)';
                        ctx.fillRect(bx, 24 - h + pulse, 6, Math.abs(h));
                    });

                    ctx.beginPath();
                    ctx.moveTo(2, 24);
                    ctx.lineTo(46, 24);
                    ctx.strokeStyle = 'rgba(100, 140, 255, 0.15)';
                    ctx.lineWidth = 0.5;
                    ctx.stroke();

                } else if (type === 'ontology') {
                    const drawBranch = (x, y, angle, len, depth) => {
                        if (depth <= 0 || len < 3) return;
                        const ex = x + Math.cos(angle) * len;
                        const ey = y + Math.sin(angle) * len;

                        ctx.beginPath();
                        ctx.moveTo(x, y);
                        ctx.lineTo(ex, ey);
                        ctx.strokeStyle = `rgba(255, 176, 0, ${0.15 + depth * 0.08})`;
                        ctx.lineWidth = depth * 0.4;
                        ctx.stroke();

                        ctx.beginPath();
                        ctx.arc(ex, ey, 1.5, 0, TAU);
                        ctx.fillStyle = `rgba(255, 176, 0, ${0.3 + depth * 0.1})`;
                        ctx.fill();

                        const sway = Math.sin(time * 0.001 + depth) * 0.1;
                        drawBranch(ex, ey, angle - 0.5 + sway, len * 0.65, depth - 1);
                        drawBranch(ex, ey, angle + 0.5 + sway, len * 0.65, depth - 1);
                    };

                    drawBranch(24, 44, -Math.PI / 2, 16, 4);
                }

                requestAnimationFrame(draw);
            }
            requestAnimationFrame(draw);
        });
    }

    /* ── STEP CANVAS VISUALIZATIONS ─────────────────────────── */
    function initStepViz() {
        const canvases = $$('.step-canvas');
        const dpr = window.devicePixelRatio || 1;

        canvases.forEach(canvas => {
            const ctx = canvas.getContext('2d');
            const viz = canvas.dataset.viz;
            canvas.width = 64 * dpr;
            canvas.height = 64 * dpr;
            ctx.scale(dpr, dpr);

            let isVisible = false;
            const observer = new IntersectionObserver(([e]) => { isVisible = e.isIntersecting; }, { threshold: 0.1 });
            observer.observe(canvas);

            function draw(time) {
                if (!isVisible) { requestAnimationFrame(draw); return; }
                ctx.clearRect(0, 0, 64, 64);

                if (viz === 'structure') {
                    // Matrix grid
                    for (let r = 0; r < 4; r++) {
                        for (let c = 0; c < 4; c++) {
                            const val = Math.sin(time * 0.002 + r * 0.5 + c * 0.7) * 0.5 + 0.5;
                            const x = 10 + c * 12;
                            const y = 10 + r * 12;
                            ctx.fillStyle = `rgba(0, 240, 255, ${val * 0.5})`;
                            ctx.fillRect(x, y, 10, 10);
                            ctx.strokeStyle = 'rgba(0, 240, 255, 0.1)';
                            ctx.lineWidth = 0.5;
                            ctx.strokeRect(x, y, 10, 10);
                        }
                    }
                } else if (viz === 'reduce') {
                    // High-dim to low-dim
                    const nPoints = 12;
                    for (let i = 0; i < nPoints; i++) {
                        const angle = (i / nPoints) * TAU + time * 0.001;
                        const r = 20 + Math.sin(time * 0.002 + i) * 5;
                        const x = 32 + Math.cos(angle) * r;
                        const y = 32 + Math.sin(angle) * r;
                        ctx.beginPath();
                        ctx.arc(x, y, 2, 0, TAU);
                        ctx.fillStyle = `rgba(255, 0, 170, ${0.3 + Math.sin(time * 0.003 + i) * 0.2})`;
                        ctx.fill();
                    }
                    // Arrow to center
                    ctx.beginPath();
                    ctx.arc(32, 32, 3, 0, TAU);
                    ctx.fillStyle = 'rgba(0, 240, 255, 0.6)';
                    ctx.fill();
                } else if (viz === 'biology') {
                    // Mini DNA
                    for (let i = 0; i < 16; i++) {
                        const t = i / 16;
                        const y = t * 64;
                        const phase = time * 0.003 + t * Math.PI * 3;
                        const x1 = 32 + Math.sin(phase) * 14;
                        const x2 = 32 + Math.sin(phase + Math.PI) * 14;
                        const d = Math.cos(phase);
                        ctx.beginPath();
                        ctx.arc(d > 0 ? x1 : x2, y, 1.5, 0, TAU);
                        ctx.fillStyle = d > 0 ? 'rgba(0, 240, 255, 0.5)' : 'rgba(255, 0, 170, 0.5)';
                        ctx.fill();
                    }
                } else if (viz === 'impact') {
                    // Pulse ring
                    const nRings = 3;
                    for (let i = 0; i < nRings; i++) {
                        const t = ((time * 0.001 + i * 1.2) % 3) / 3;
                        const r = t * 28;
                        const alpha = (1 - t) * 0.4;
                        ctx.beginPath();
                        ctx.arc(32, 32, r, 0, TAU);
                        ctx.strokeStyle = `rgba(0, 240, 255, ${alpha})`;
                        ctx.lineWidth = 1;
                        ctx.stroke();
                    }
                    ctx.beginPath();
                    ctx.arc(32, 32, 4, 0, TAU);
                    ctx.fillStyle = 'rgba(0, 240, 255, 0.5)';
                    ctx.fill();
                }

                requestAnimationFrame(draw);
            }
            requestAnimationFrame(draw);
        });
    }

    /* ── SPECTRAL EMBEDDING INTERACTIVE ─────────────────────── */
    function initSpectralEmbedding() {
        const canvas3D = $('#spectralCanvas');
        const canvasNet = $('#networkGraphCanvas');
        const canvasSpec = $('#spectrumCanvas');
        if (!canvas3D || !canvasNet || !canvasSpec) return;

        const ctx3D = canvas3D.getContext('2d');
        const ctxNet = canvasNet.getContext('2d');
        const ctxSpec = canvasSpec.getContext('2d');
        const dpr = window.devicePixelRatio || 1;

        // Setup canvases
        const W3D = 700, H3D = 500;
        canvas3D.width = W3D * dpr;
        canvas3D.height = H3D * dpr;
        ctx3D.scale(dpr, dpr);

        const WN = 360, HN = 360;
        canvasNet.width = WN * dpr;
        canvasNet.height = HN * dpr;
        ctxNet.scale(dpr, dpr);

        const WS = 360, HS = 140;
        canvasSpec.width = WS * dpr;
        canvasSpec.height = HS * dpr;
        ctxSpec.scale(dpr, dpr);

        // Generate a network with clear community structure
        const N = 24;
        const communitySize = [8, 8, 8];
        const communityColors = ['#00f0ff', '#ff4060', '#ffb000'];
        const communityAssignment = [];

        for (let c = 0; c < communitySize.length; c++) {
            for (let i = 0; i < communitySize[c]; i++) {
                communityAssignment.push(c);
            }
        }

        // Build adjacency
        const adj = Array.from({ length: N }, () => new Array(N).fill(0));
        const pIntra = 0.55;
        const pInter = 0.06;

        for (let i = 0; i < N; i++) {
            for (let j = i + 1; j < N; j++) {
                const same = communityAssignment[i] === communityAssignment[j];
                if (Math.random() < (same ? pIntra : pInter)) {
                    adj[i][j] = 1;
                    adj[j][i] = 1;
                }
            }
        }

        // Ensure connectivity: add edges for isolated nodes
        for (let i = 0; i < N; i++) {
            let deg = 0;
            for (let j = 0; j < N; j++) deg += adj[i][j];
            if (deg === 0) {
                // Connect to a random same-community node
                for (let j = 0; j < N; j++) {
                    if (j !== i && communityAssignment[j] === communityAssignment[i]) {
                        adj[i][j] = 1;
                        adj[j][i] = 1;
                        break;
                    }
                }
            }
        }

        // Identify bridge nodes (connected to multiple communities)
        const nodeTypes = communityAssignment.map((c, i) => {
            const neighborComms = new Set();
            for (let j = 0; j < N; j++) {
                if (adj[i][j]) neighborComms.add(communityAssignment[j]);
            }
            return neighborComms.size > 1 ? 'bridge' : 'member';
        });

        const nodeColors = communityAssignment.map((c, i) => {
            if (nodeTypes[i] === 'bridge') return '#a855f7';
            return communityColors[c];
        });

        // Compute Laplacian and eigendecomposition
        const L = graphLaplacian(adj);
        const eigen = jacobiEigen(L);

        // Use eigenvectors 1, 2, 3 for embedding (skip constant eigenvector 0)
        const embed3D = [];
        const embedScale = 120;
        for (let i = 0; i < N; i++) {
            embed3D.push({
                x: (eigen.vectors[1] ? eigen.vectors[1][i] : 0) * embedScale,
                y: (eigen.vectors[2] ? eigen.vectors[2][i] : 0) * embedScale,
                z: (eigen.vectors[3] ? eigen.vectors[3][i] : 0) * embedScale,
                color: nodeColors[i],
                community: communityAssignment[i],
                type: nodeTypes[i]
            });
        }

        // 2D network layout (force-directed simple)
        const netPos = [];
        for (let i = 0; i < N; i++) {
            const c = communityAssignment[i];
            const angle = Math.random() * TAU;
            const r = Math.random() * 40;
            const cx = WN / 2 + (c === 0 ? -70 : c === 1 ? 70 : 0);
            const cy = HN / 2 + (c === 2 ? 70 : c === 0 ? -30 : -30);
            netPos.push({ x: cx + Math.cos(angle) * r, y: cy + Math.sin(angle) * r });
        }

        // Simple force simulation
        for (let iter = 0; iter < 200; iter++) {
            for (let i = 0; i < N; i++) {
                let fx = 0, fy = 0;
                for (let j = 0; j < N; j++) {
                    if (i === j) continue;
                    const dx = netPos[i].x - netPos[j].x;
                    const dy = netPos[i].y - netPos[j].y;
                    const dist = Math.sqrt(dx * dx + dy * dy) + 0.1;

                    // Repulsion
                    fx += (dx / dist) * 800 / (dist * dist);
                    fy += (dy / dist) * 800 / (dist * dist);

                    // Attraction (edges)
                    if (adj[i][j]) {
                        fx -= (dx / dist) * dist * 0.01;
                        fy -= (dy / dist) * dist * 0.01;
                    }
                }

                // Center gravity
                fx -= (netPos[i].x - WN / 2) * 0.005;
                fy -= (netPos[i].y - HN / 2) * 0.005;

                netPos[i].x += fx * 0.3;
                netPos[i].y += fy * 0.3;
            }
        }

        // 3D rotation state
        let rotX = -0.2, rotY = 0.4;
        let isDragging = false;
        let lastMouse = { x: 0, y: 0 };
        let autoRotate = true;
        let hoveredNode = -1;

        canvas3D.addEventListener('mousedown', (e) => {
            isDragging = true;
            autoRotate = false;
            lastMouse = { x: e.clientX, y: e.clientY };
        });
        window.addEventListener('mouseup', () => { isDragging = false; });
        window.addEventListener('mousemove', (e) => {
            if (isDragging) {
                rotY += (e.clientX - lastMouse.x) * 0.005;
                rotX += (e.clientY - lastMouse.y) * 0.005;
                lastMouse = { x: e.clientX, y: e.clientY };
            }
        });

        canvas3D.addEventListener('touchstart', (e) => {
            if (e.touches.length === 1) {
                isDragging = true;
                autoRotate = false;
                lastMouse = { x: e.touches[0].clientX, y: e.touches[0].clientY };
                e.preventDefault();
            }
        }, { passive: false });
        canvas3D.addEventListener('touchmove', (e) => {
            if (!isDragging || e.touches.length !== 1) return;
            rotY += (e.touches[0].clientX - lastMouse.x) * 0.005;
            rotX += (e.touches[0].clientY - lastMouse.y) * 0.005;
            lastMouse = { x: e.touches[0].clientX, y: e.touches[0].clientY };
            e.preventDefault();
        }, { passive: false });
        canvas3D.addEventListener('touchend', () => { isDragging = false; });

        // Hover on network
        canvasNet.addEventListener('mousemove', (e) => {
            const r = canvasNet.getBoundingClientRect();
            const mx = (e.clientX - r.left) * (WN / r.width);
            const my = (e.clientY - r.top) * (HN / r.height);

            hoveredNode = -1;
            for (let i = 0; i < N; i++) {
                const dx = mx - netPos[i].x;
                const dy = my - netPos[i].y;
                if (dx * dx + dy * dy < 200) {
                    hoveredNode = i;
                    break;
                }
            }
        });

        canvasNet.addEventListener('mouseleave', () => { hoveredNode = -1; });

        function project(x, y, z) {
            let x1 = x * Math.cos(rotY) - z * Math.sin(rotY);
            let z1 = x * Math.sin(rotY) + z * Math.cos(rotY);
            let y1 = y * Math.cos(rotX) - z1 * Math.sin(rotX);
            let z2 = y * Math.sin(rotX) + z1 * Math.cos(rotX);
            const fov = 500;
            const scale = fov / (fov + z2 + 200);
            return { x: W3D / 2 + x1 * scale, y: H3D / 2 + y1 * scale, z: z2, scale };
        }

        let isVisible = false;
        const observer = new IntersectionObserver(([e]) => { isVisible = e.isIntersecting; }, { threshold: 0.05 });
        observer.observe(canvas3D);

        function draw(time) {
            if (!isVisible) { requestAnimationFrame(draw); return; }

            if (autoRotate) {
                rotY += 0.0015;
                rotX = -0.2 + Math.sin(time * 0.0003) * 0.06;
            }

            // ── 3D EMBEDDING ──
            ctx3D.clearRect(0, 0, W3D, H3D);

            const projected = embed3D.map((p, i) => ({
                ...project(p.x, p.y, p.z),
                idx: i,
                color: p.color,
                community: p.community,
                type: p.type
            }));

            // Draw edges in 3D
            for (let i = 0; i < N; i++) {
                for (let j = i + 1; j < N; j++) {
                    if (!adj[i][j]) continue;
                    const pi = projected.find(p => p.idx === i);
                    const pj = projected.find(p => p.idx === j);
                    const same = communityAssignment[i] === communityAssignment[j];
                    ctx3D.beginPath();
                    ctx3D.moveTo(pi.x, pi.y);
                    ctx3D.lineTo(pj.x, pj.y);
                    ctx3D.strokeStyle = same
                        ? `rgba(100, 140, 255, 0.08)`
                        : `rgba(168, 85, 247, 0.12)`;
                    ctx3D.lineWidth = hoveredNode === i || hoveredNode === j ? 1.5 : 0.5;
                    ctx3D.stroke();
                }
            }

            // Sort by z for proper rendering
            projected.sort((a, b) => a.z - b.z);

            for (const p of projected) {
                const r = Math.max(3, 7 * p.scale);
                const isHovered = hoveredNode === p.idx;
                const isNeighbor = hoveredNode >= 0 && adj[hoveredNode][p.idx];
                const dimmed = hoveredNode >= 0 && !isHovered && !isNeighbor;

                // Glow
                const grd = ctx3D.createRadialGradient(p.x, p.y, 0, p.x, p.y, r * 4);
                grd.addColorStop(0, p.color + (isHovered ? '60' : '30'));
                grd.addColorStop(1, p.color + '00');
                ctx3D.beginPath();
                ctx3D.arc(p.x, p.y, r * 4, 0, TAU);
                ctx3D.fillStyle = grd;
                ctx3D.fill();

                // Point
                ctx3D.beginPath();
                ctx3D.arc(p.x, p.y, isHovered ? r * 1.5 : r, 0, TAU);
                ctx3D.fillStyle = p.color;
                ctx3D.globalAlpha = dimmed ? 0.15 : (0.4 + p.scale * 0.6);
                ctx3D.fill();
                ctx3D.globalAlpha = 1;

                // Node label on hover
                if (isHovered) {
                    ctx3D.font = '10px "JetBrains Mono", monospace';
                    ctx3D.fillStyle = p.color;
                    ctx3D.textAlign = 'center';
                    ctx3D.fillText(`node ${p.idx}`, p.x, p.y - r * 2 - 6);
                    ctx3D.font = '8px "JetBrains Mono", monospace';
                    ctx3D.fillStyle = 'rgba(148, 163, 184, 0.6)';
                    ctx3D.fillText(p.type === 'bridge' ? 'bridge' : `community ${p.community + 1}`, p.x, p.y - r * 2 + 6);
                }
            }

            // Axes
            const axisLen = 140;
            const axes = [
                { label: 'v₂', end: project(axisLen, 0, 0) },
                { label: 'v₃', end: project(0, axisLen, 0) },
                { label: 'v₄', end: project(0, 0, axisLen) }
            ];
            const origin = project(0, 0, 0);

            axes.forEach(ax => {
                ctx3D.beginPath();
                ctx3D.moveTo(origin.x, origin.y);
                ctx3D.lineTo(ax.end.x, ax.end.y);
                ctx3D.strokeStyle = 'rgba(148, 163, 184, 0.15)';
                ctx3D.lineWidth = 0.8;
                ctx3D.setLineDash([3, 3]);
                ctx3D.stroke();
                ctx3D.setLineDash([]);

                ctx3D.font = '10px "JetBrains Mono", monospace';
                ctx3D.fillStyle = 'rgba(148, 163, 184, 0.35)';
                ctx3D.textAlign = 'center';
                ctx3D.fillText(ax.label, ax.end.x, ax.end.y - 8);
            });

            // ── NETWORK GRAPH ──
            ctxNet.clearRect(0, 0, WN, HN);

            // Edges
            for (let i = 0; i < N; i++) {
                for (let j = i + 1; j < N; j++) {
                    if (!adj[i][j]) continue;
                    const isActive = hoveredNode === i || hoveredNode === j;
                    const same = communityAssignment[i] === communityAssignment[j];

                    ctxNet.beginPath();
                    ctxNet.moveTo(netPos[i].x, netPos[i].y);
                    ctxNet.lineTo(netPos[j].x, netPos[j].y);
                    ctxNet.strokeStyle = isActive
                        ? (same ? 'rgba(0, 240, 255, 0.4)' : 'rgba(168, 85, 247, 0.5)')
                        : (same ? 'rgba(100, 140, 255, 0.1)' : 'rgba(168, 85, 247, 0.08)');
                    ctxNet.lineWidth = isActive ? 1.5 : 0.6;
                    ctxNet.stroke();
                }
            }

            // Nodes
            for (let i = 0; i < N; i++) {
                const isHovered = hoveredNode === i;
                const isNeighbor = hoveredNode >= 0 && adj[hoveredNode][i];
                const dimmed = hoveredNode >= 0 && !isHovered && !isNeighbor;
                const color = nodeColors[i];
                const r = isHovered ? 7 : 5;

                // Glow
                if (isHovered || isNeighbor) {
                    const grd = ctxNet.createRadialGradient(netPos[i].x, netPos[i].y, 0,
                        netPos[i].x, netPos[i].y, 18);
                    grd.addColorStop(0, color + '40');
                    grd.addColorStop(1, color + '00');
                    ctxNet.beginPath();
                    ctxNet.arc(netPos[i].x, netPos[i].y, 18, 0, TAU);
                    ctxNet.fillStyle = grd;
                    ctxNet.fill();
                }

                ctxNet.beginPath();
                ctxNet.arc(netPos[i].x, netPos[i].y, r, 0, TAU);
                ctxNet.fillStyle = color;
                ctxNet.globalAlpha = dimmed ? 0.15 : 0.7;
                ctxNet.fill();
                ctxNet.globalAlpha = 1;
            }

            // ── SPECTRUM BARS ──
            ctxSpec.clearRect(0, 0, WS, HS);

            const barCount = Math.min(N, 12);
            const barW = 20;
            const barSp = (WS - 32) / barCount;
            const barStartX = 16;
            const barBaseline = HS - 28;
            const maxEig = Math.max(...eigen.values.slice(0, barCount), 0.01);

            ctxSpec.font = '7px "JetBrains Mono", monospace';
            ctxSpec.fillStyle = 'rgba(100, 116, 139, 0.3)';
            ctxSpec.textAlign = 'left';
            ctxSpec.fillText('EIGENVALUES OF L', barStartX, 14);

            for (let i = 0; i < barCount; i++) {
                const bx = barStartX + i * barSp;
                const val = eigen.values[i];
                const h = (val / maxEig) * (HS - 50);

                const t = i / barCount;
                const cr = Math.round(lerp(0, 255, t));
                const cg = Math.round(lerp(240, 60, t));

                ctxSpec.fillStyle = `rgba(${cr}, ${cg}, 255, 0.4)`;
                ctxSpec.fillRect(bx, barBaseline - h, barW - 4, h);

                ctxSpec.strokeStyle = `rgba(${cr}, ${cg}, 255, 0.2)`;
                ctxSpec.lineWidth = 0.5;
                ctxSpec.strokeRect(bx, barBaseline - h, barW - 4, h);

                ctxSpec.font = '6px "JetBrains Mono", monospace';
                ctxSpec.fillStyle = 'rgba(148, 163, 184, 0.4)';
                ctxSpec.textAlign = 'center';
                ctxSpec.fillText(val.toFixed(1), bx + (barW - 4) / 2, barBaseline + 10);

                // Mark the zero eigenvalue
                if (Math.abs(val) < 0.01) {
                    ctxSpec.fillStyle = 'rgba(0, 240, 255, 0.5)';
                    ctxSpec.fillText('≈0', bx + (barW - 4) / 2, barBaseline + 20);
                }
            }

            // Algebraic connectivity marker
            if (eigen.values.length > 1) {
                const bx2 = barStartX + 1 * barSp;
                ctxSpec.font = '6px "JetBrains Mono", monospace';
                ctxSpec.fillStyle = 'rgba(0, 240, 255, 0.4)';
                ctxSpec.textAlign = 'center';
                ctxSpec.fillText('← Fiedler', bx2 + barSp + 8, barBaseline - (eigen.values[1] / maxEig) * (HS - 50) - 6);
            }

            requestAnimationFrame(draw);
        }

        requestAnimationFrame(draw);
    }

    /* ── CHROMOSOME SCROLL ──────────────────────────────────── */
    function initChromoScroll() {
        const canvas = $('#chromoCanvas');
        const marker = $('#chromoMarker');
        if (!canvas || !marker) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;

        function resize() {
            const h = Math.min(500, window.innerHeight * 0.7);
            canvas.width = 32 * dpr;
            canvas.height = h * dpr;
            canvas.style.width = '24px';
            canvas.style.height = h + 'px';
            ctx.setTransform(1, 0, 0, 1, 0, 0);
            ctx.scale(dpr, dpr);
        }
        resize();
        window.addEventListener('resize', resize);

        function draw(time) {
            const h = parseInt(canvas.style.height);
            ctx.clearRect(0, 0, 32, h);

            const scrollMax = document.documentElement.scrollHeight - window.innerHeight;
            const scrollT = scrollMax > 0 ? window.scrollY / scrollMax : 0;

            // Chromosome bands
            const bands = 22;
            const bandH = h / bands;
            for (let i = 0; i < bands; i++) {
                const y = i * bandH;
                const intensity = Math.sin(time * 0.001 + i * 0.5) * 0.1 + 0.15;
                const isActive = Math.abs(i / bands - scrollT) < 0.08;

                ctx.fillStyle = isActive
                    ? `rgba(0, 240, 255, ${intensity + 0.15})`
                    : `rgba(100, 140, 255, ${intensity})`;
                ctx.fillRect(4, y + 1, 16, bandH - 2);

                if (i % 4 === 0) {
                    ctx.fillStyle = 'rgba(100, 116, 139, 0.15)';
                    ctx.font = '5px "JetBrains Mono", monospace';
                    ctx.textAlign = 'right';
                    ctx.fillText(`${i + 1}`, 24, y + bandH / 2 + 2);
                }
            }

            // Centromere
            const centY = h * 0.4;
            ctx.beginPath();
            ctx.moveTo(4, centY - 3);
            ctx.lineTo(20, centY);
            ctx.lineTo(4, centY + 3);
            ctx.fillStyle = 'rgba(255, 0, 170, 0.15)';
            ctx.fill();

            // Marker position
            marker.style.top = (scrollT * h) + 'px';

            requestAnimationFrame(draw);
        }
        requestAnimationFrame(draw);
    }

    /* ── SCROLL REVEAL ──────────────────────────────────────── */
    function initScrollReveal() {
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(e => {
                if (e.isIntersecting) {
                    e.target.classList.add('in-view');

                    // Animate project bars
                    e.target.querySelectorAll('.project-bar-fill').forEach(bar => {
                        bar.style.setProperty('--bar-width', bar.style.width);
                    });
                }
            });
        }, { threshold: 0.1, rootMargin: '-40px' });

        $$('.section').forEach(sec => observer.observe(sec));
    }

    /* ── NAVIGATION SCROLLSPY ───────────────────────────────── */
    function initScrollSpy() {
        const navItems = $$('.nav-item');
        const sections = $$('.section');
        let currentSection = 'hero';

        function updateNav() {
            const scrollY = window.scrollY + window.innerHeight * 0.35;

            sections.forEach(sec => {
                const top = sec.offsetTop;
                const bottom = top + sec.offsetHeight;
                if (scrollY >= top && scrollY < bottom) {
                    currentSection = sec.id;
                }
            });

            navItems.forEach(item => {
                item.classList.toggle('active', item.dataset.section === currentSection);
            });
        }

        window.addEventListener('scroll', updateNav, { passive: true });
        updateNav();

        navItems.forEach(item => {
            item.addEventListener('click', (e) => {
                e.preventDefault();
                const target = document.getElementById(item.dataset.section);
                if (target) {
                    const y = target.offsetTop - 20;
                    window.scrollTo({ top: y, behavior: 'smooth' });
                }

                const nav = $('#specNav');
                const toggle = $('#navToggle');
                if (nav) nav.classList.remove('mobile-open');
                if (toggle) toggle.classList.remove('active');
            });
        });
    }

    /* ── MOBILE NAV ─────────────────────────────────────────── */
    function initMobileNav() {
        const toggle = $('#navToggle');
        const nav = $('#specNav');
        if (!toggle || !nav) return;

        toggle.addEventListener('click', () => {
            toggle.classList.toggle('active');
            nav.classList.toggle('mobile-open');
        });
    }

    /* ── INIT ───────────────────────────────────────────────── */
    document.addEventListener('DOMContentLoaded', () => {
        initLoader();
        initNetwork();
        initOverlineHelix();
        initAvatarRing();
        initEigenViz();
        initCellViz();
        initHelix();
        initEquation();
        initProjectMiniViz();
        initStepViz();
        initSpectralEmbedding();
        initChromoScroll();
        initScrollReveal();
        initScrollSpy();
        initMobileNav();
    });
})();