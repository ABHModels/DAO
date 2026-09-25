#!/usr/bin/env python3
"""
DAO — X-ray Reflection Spectroscopy Model
Web-based parameter configurator & launcher
Author: Yimin Huang
"""

import json
import os
import threading
import webbrowser
from flask import Flask, render_template_string, request, jsonify, abort

app = Flask(__name__, static_folder='image', static_url_path='/image')

# ─── Working directory (where maindaocl lives) ─────────────────
WORK_DIR = os.path.dirname(os.path.abspath(__file__))

# ─── Model definitions (mirrors params.cpp) ────────────────────
# "scale":"log" marks parameters whose sliders run in log10 space (the
# emitted command always carries the real linear value).
CORONA_MODELS = {
    "powerlaw": {
        "label": "Power Law",
        "desc": "N(E) ∝ E<sup>−Γ</sup> · exp(−E<sub>lo</sub>/E)  [photons cm⁻² s⁻¹ keV⁻¹]",
        "params": [
            {"name": "Gamma", "label": "Γ", "desc": "Photon index (photon-number spectrum)", "default": 2.0, "min": 1.0, "max": 4.0, "step": 0.01},
            {"name": "E_low_cut", "label": "E<sub>lo</sub>", "desc": "Low-energy exp cutoff [keV]", "default": 0.1, "min": 0.01, "max": 10, "step": 0.01},
        ],
    },
    "cutoffpl": {
        "label": "Cutoff Power Law",
        "desc": "N(E) ∝ E<sup>−Γ</sup> · exp(−E/E<sub>cut</sub>) · exp(−E<sub>lo</sub>/E)  [photons cm⁻² s⁻¹ keV⁻¹]",
        "params": [
            {"name": "Gamma", "label": "Γ", "desc": "Photon index (photon-number spectrum)", "default": 2.0, "min": 1.0, "max": 4.0, "step": 0.01},
            {"name": "Ecut",  "label": "E<sub>cut</sub>", "desc": "High-energy cutoff [keV]", "default": 300.0, "min": 10, "max": 2000, "step": 1, "scale": "log"},
            {"name": "E_low_cut", "label": "E<sub>lo</sub>", "desc": "Low-energy exp cutoff [keV]", "default": 0.1, "min": 0.01, "max": 10, "step": 0.01},
        ],
    },
    "nthcomp": {
        "label": "NTHComp",
        "desc": "Thermal Comptonisation (Zdziarski+ 1996)",
        "params": [
            {"name": "Gamma", "label": "Γ", "desc": "Photon index", "default": 2.0, "min": 1.0, "max": 4.0, "step": 0.01},
            {"name": "kT_e",  "label": "kT<sub>e</sub>", "desc": "Electron temperature [keV]", "default": 60.0, "min": 1, "max": 500, "step": 0.5, "scale": "log"},
            {"name": "kT_bb", "label": "kT<sub>bb</sub>", "desc": "Seed photon temperature [keV]", "default": 0.1, "min": 0.001, "max": 5, "step": 0.001, "scale": "log"},
        ],
    },
    "comptt": {
        "label": "CompTT",
        "desc": "Comptonisation model (Titarchuk 1994)",
        "params": [
            {"name": "kT_e",  "label": "kT<sub>e</sub>", "desc": "Plasma temperature [keV]", "default": 50.0, "min": 1, "max": 500, "step": 0.5, "scale": "log"},
            {"name": "kT_bb", "label": "kT<sub>bb</sub>", "desc": "Soft photon temperature [keV]", "default": 0.05, "min": 0.001, "max": 5, "step": 0.001, "scale": "log"},
            {"name": "taup",  "label": "τ<sub>p</sub>", "desc": "Plasma optical depth", "default": 1.0, "min": 0.01, "max": 20, "step": 0.01, "scale": "log"},
        ],
    },
    "blackbody": {
        "label": "Blackbody",
        "desc": "B(E) = (2E³/h²c²) / [exp(E/kT) − 1]",
        "params": [
            {"name": "kT_bb", "label": "kT<sub>bb</sub>", "desc": "Temperature [keV]", "default": 0.05, "min": 0.001, "max": 50, "step": 0.001, "scale": "log"},
        ],
    },
}

SLAB_PARAMS = [
    {"name": "nh",        "label": "log n<sub>H</sub>",  "desc": "Hydrogen density [cm⁻³]",      "default": 15.0,   "min": 10,  "max": 22,   "step": 0.1,    "scale": "log"},
    {"name": "zeta",      "label": "log ξ",              "desc": "Ionisation parameter",          "default": 3.0,    "min": 0,   "max": 6,    "step": 0.1,    "scale": "log"},
    {"name": "frac",      "label": "f<sub>cor</sub>",    "desc": "F<sub>corona</sub> / F<sub>disk</sub>; ≤ 0 → corona only (no disk)", "default": -1, "min": -1, "max": 100, "step": 0.01},
    {"name": "incidence", "label": "cos θ",              "desc": "Incidence angle cosine",        "default": 0.7071, "min": 0.01,"max": 1,    "step": 0.0001},
    {"name": "Afe",       "label": "A<sub>Fe</sub>",     "desc": "Iron abundance [solar]",        "default": 1.0,    "min": 0.1, "max": 10,   "step": 0.1},
]

DISK_PARAMS = [
    {"name": "kT_disk", "label": "kT<sub>disk</sub>", "desc": "Disk blackbody temperature [keV]", "default": 0.35, "min": 0.001, "max": 5, "step": 0.001, "scale": "log"},
]


# ═══════════════════════════════════════════════════════════════
#  DESIGN SYSTEM — single source of truth (injected into all pages)
# ═══════════════════════════════════════════════════════════════
# Shared instrument theme: warm light surfaces and readable data.

BASE_CSS = r"""
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  :root {
    color-scheme:light;
    --bg:#f5f2ed; --bg2:#f8f5f0; --card:#fffdfa; --card-hi:#fffdfa;
    --border:#ebe5dc; --border-hi:#d4c9bc;
    --accent:#915237; --accent2:#753b28; --accent-sft:#f5e9df;
    --cyan:#855335; --green:#57715a; --red:#ad4d42;
    --white:#26221f; --text:#4a443e; --dim:#6f675e; --faint:#7c746b;
    --r-sm:8px; --r-md:16px; --r-lg:18px;
    /* Typography */
    --mono:'SF Mono','JetBrains Mono','Fira Code','Cascadia Code',ui-monospace,monospace;
    --fs-xs:.8125rem; --fs-sm:.875rem; --fs-base:.92rem; --fs-md:1.05rem;
    --fs-h2:.80rem; --fs-h1:1.5rem; --fs-hero:2.6rem;
    /* Spacing (4px base) */
    --sp-1:4px; --sp-2:8px; --sp-3:12px; --sp-4:16px; --sp-5:24px; --sp-6:32px; --sp-7:48px;
    --shadow:0 8px 28px rgba(56,39,26,.045); --shadow-hi:0 12px 32px rgba(56,39,26,.065);
  }
  html { scroll-behavior:smooth; }
  body {
    font-family:'Inter','Helvetica Neue','Segoe UI',system-ui,sans-serif;
    background:var(--bg); color:var(--text); min-height:100vh;
    overflow-x:hidden; font-size:15px; line-height:1.6; letter-spacing:.02em;
    -webkit-font-smoothing:antialiased;
  }

  a { color:var(--cyan); text-decoration:none; }
  a:hover { color:var(--white); }
  code, .mono { font-family:var(--mono); }

  /* ── Accessibility primitives ───────────────────────── */
  :focus-visible { outline:2px solid var(--accent); outline-offset:2px; border-radius:2px; }
  .sr-only { position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap;border:0; }
  .skip-link { position:absolute;left:-999px;top:0;z-index:1000; }
  .skip-link:focus { left:8px;top:8px;background:var(--card);color:var(--white);padding:8px 14px;border-radius:var(--r-sm);border:1px solid var(--accent); }
  @media (prefers-reduced-motion: reduce){ *{animation:none!important; transition:none!important; scroll-behavior:auto!important;} }

  /* ── Sticky top nav (shared) ────────────────────────── */
  header.nav { position:sticky; top:0; z-index:50;
    background:var(--card); border-bottom:1px solid var(--border); }
  header.nav nav { max-width:1200px; margin:0 auto; padding:0 var(--sp-5);
    display:flex; align-items:center; gap:var(--sp-6); height:54px; }
  header.nav .brand { font-weight:600; font-size:1.05rem; color:var(--accent);
    letter-spacing:.14em; display:flex; align-items:center; gap:10px; flex:none; }
  header.nav .brand:hover { color:var(--accent); }
  header.nav .tabs { display:flex; gap:var(--sp-5); align-items:center;
    overflow-x:auto; white-space:nowrap; min-width:0; }
  header.nav .tabs a { font-size:var(--fs-sm); font-weight:500; color:var(--dim);
    padding:16px 2px; border-bottom:2px solid transparent; transition:color .2s,border-color .2s; }
  header.nav .tabs a:hover { color:var(--white); }
  header.nav .tabs a.active { color:var(--white); font-weight:600; border-bottom-color:var(--accent); }
  @media (max-width:420px){ header.nav .brand { font-size:.9rem; } header.nav nav { gap:var(--sp-4); padding:0 var(--sp-3); } }

  /* ── Page wrapper ───────────────────────────────────── */
  .wrapper { position:relative; z-index:1; max-width:1120px; margin:0 auto; padding:var(--sp-5) var(--sp-5) var(--sp-7); }
  .page-head { padding:var(--sp-5) 0 var(--sp-3); }
  .page-head h1 { font-size:var(--fs-h1); font-weight:700; letter-spacing:.06em; color:var(--white); }
  .page-head .sub { font-size:var(--fs-sm); color:var(--dim); margin-top:var(--sp-2); }

  /* ── Cards (shared physics) ─────────────────────────── */
  .card, .plot-card, .formula-box, .note-box, .params-box {
    background:var(--card); border:1px solid transparent; border-radius:var(--r-md);
    box-shadow:var(--shadow);
  }
  .card { padding:var(--sp-5); transition:border-color .2s,box-shadow .2s,background .2s; margin-bottom:var(--sp-5); }
  .card:hover { border-color:var(--border-hi); box-shadow:var(--shadow-hi); background:var(--card-hi); }
  .card.full { grid-column:1/-1; }
  .grid { display:grid; grid-template-columns:1fr 1fr; gap:var(--sp-5); }
  @media (max-width:760px){ .grid { grid-template-columns:1fr; } }

  .card h2, .plot-card h2 {
    font-size:var(--fs-h2); font-weight:600; text-transform:uppercase; letter-spacing:.14em;
    color:var(--accent); display:flex; align-items:center; gap:10px; margin-bottom:var(--sp-4);
  }
  .card h2 .tick, .plot-card h2 .tick { width:2px; height:14px; background:var(--accent); border-radius:1px; flex:none; }

  /* ── Footer (shared) ────────────────────────────────── */
  footer.site { text-align:center; color:var(--faint); font-size:var(--fs-xs);
    border-top:1px solid var(--border); padding:var(--sp-6) 0 var(--sp-4); margin-top:var(--sp-6); letter-spacing:.04em; }
  footer.site a { color:var(--faint); }
  footer.site a:hover { color:var(--dim); }

  /* ── Toast (shared) ─────────────────────────────────── */
  .toast { position:fixed; bottom:32px; left:50%; transform:translateX(-50%);
    background:var(--card); border:1px solid var(--border); padding:12px 24px;
    border-radius:var(--r-sm); font-size:var(--fs-sm); font-weight:500; opacity:0;
    transition:opacity .3s; pointer-events:none; z-index:99; box-shadow:var(--shadow-hi);
    display:flex; align-items:center; gap:10px; }
  .toast.show { opacity:1; }
  .toast.success { color:var(--green); border-color:rgba(34,197,94,.3); }
  .toast.error   { color:var(--red);   border-color:rgba(239,83,80,.3); }
  .toast.info    { color:var(--cyan);  border-color:rgba(145,82,55,.3); }

  /* ── Pills (shared, used on /plots) ─────────────────── */
  .pill { font-size:var(--fs-xs); padding:5px 14px; border-radius:999px; cursor:pointer;
    font-family:inherit; font-weight:500; letter-spacing:.03em; border:1px solid; transition:all .2s; background:transparent; }
  .pill-fe   { border-color:rgba(239,83,80,.35);  background:rgba(239,83,80,.10);  color:var(--red); }
  .pill-line { border-color:rgba(145,82,55,.35); background:rgba(145,82,55,.10); color:var(--cyan); }
  .pill.off  { border-color:var(--border); background:transparent; color:var(--faint); }

  /* ── Buttons (shared) ───────────────────────────────── */
  .btn { padding:10px 24px; border-radius:var(--r-sm); font-size:var(--fs-sm); font-weight:600;
    cursor:pointer; border:none; font-family:inherit; transition:all .2s; letter-spacing:.03em; }
  .btn-primary { background:var(--accent); color:white; }
  .btn-primary:hover { background:var(--accent2); }
  .btn-secondary { background:var(--bg2); color:var(--text); border:1px solid var(--border); }
  .btn-secondary:hover { border-color:var(--border-hi); background:var(--card-hi); color:var(--white); }

  /* ── Empty state (plots/convergence) ────────────────── */
  .empty-state { text-align:center; padding:var(--sp-7) var(--sp-5); }
  .empty-state .ico { font-size:2.4rem; color:var(--faint); margin-bottom:var(--sp-4); }
  .empty-state h3 { color:var(--white); font-size:1.1rem; margin-bottom:var(--sp-2); }
  .empty-state p { color:var(--dim); font-size:var(--fs-sm); max-width:420px; margin:0 auto var(--sp-5); }
  /* Shared light instrument theme: all routes use the same navigation and controls. */
  body { background:var(--bg); letter-spacing:0; font-size:14px; }
  body::before { display:none; }
  header.nav { background:#fffdfa; backdrop-filter:none; border-bottom-color:#eee9e2; }
  header.nav nav { max-width:1240px; height:60px; gap:40px; }
  header.nav .brand { font-size:1.05rem; letter-spacing:.1em; color:var(--white); }
  .brand-orbit { display:none; }
  header.nav .tabs { gap:28px; }
  header.nav .tabs a { font-size:14px; padding:19px 0; }
  header.nav .tabs a.active { color:var(--accent); }
  .wrapper { max-width:1200px; padding-top:20px; }
  .card { background:var(--card); padding:24px; margin-bottom:16px; transition:box-shadow .2s; }
  .card:hover { background:var(--card); border-color:transparent; box-shadow:var(--shadow-hi); }
  .card h2 { color:var(--white); text-transform:none; letter-spacing:0; font-size:1rem; margin-bottom:16px; }
  .btn { font-size:13px; padding:9px 14px; letter-spacing:0; border-radius:9px; }
  .btn-primary { background:var(--accent); color:white; box-shadow:none; }
  .btn-primary:hover { background:var(--accent2); box-shadow:none; transform:none; }
  .btn-secondary { background:white; }
  footer.site { margin-top:0; padding:18px 24px; font-size:12px; }
  .page-head { padding:8px 0 20px; margin-bottom:20px; border-bottom:1px solid var(--border); }
  .page-head h1 { font-size:1.7rem; font-weight:600; letter-spacing:-.035em; }
  .page-head .sub { font-size:var(--fs-sm); margin-top:4px; }
  .plot-card h2 { color:var(--white); text-transform:none; letter-spacing:0; font-size:1rem;
    margin-bottom:16px; flex-wrap:wrap; }
  .card h2 .tick, .plot-card h2 .tick { display:none; }
  .pill { border-radius:9px; font-size:12px; }
  .pill-fe { background:#fbf0ef; border-color:#e5c3c1; }
  .pill-line { background:var(--accent-sft); border-color:#dfc7b6; }
  .pill.off { background:transparent; color:var(--dim); border-color:var(--border); }
  @media(max-width:520px) {
    .wrapper { padding:16px 12px; }
    header.nav nav { padding:0 12px; gap:20px; }
    .page-head h1 { font-size:1.4rem; }
  }

  /* ── Run identity card (plots/convergence) ──────────── */
  .run-card { padding:14px 18px; margin-bottom:var(--sp-5); }
  .run-card .rc-top { display:flex; align-items:center; gap:10px; flex-wrap:wrap; margin-bottom:10px; }
  .run-card .rc-model { font-size:var(--fs-md); font-weight:700; color:var(--white); letter-spacing:.03em; }
  .rc-badge { font-size:var(--fs-xs); padding:3px 10px; border-radius:999px; border:1px solid var(--border);
    letter-spacing:.03em; white-space:nowrap; }
  .rc-badge.b-test   { color:var(--red);    border-color:rgba(173,77,66,.35);  background:rgba(173,77,66,.08); }
  .rc-badge.b-prod   { color:var(--cyan);   border-color:rgba(133,83,53,.3);  background:rgba(133,83,53,.06); }
  .rc-badge.b-kernel { color:var(--accent); border-color:rgba(145,82,55,.35); background:var(--accent-sft); }
  .rc-badge.b-illum  { color:var(--text);   border-color:var(--border-hi); }
  .rc-badge.b-ok     { color:var(--green);  border-color:rgba(87,113,90,.35);  background:rgba(87,113,90,.08); }
  .rc-badge.b-warn   { color:var(--accent); border-color:rgba(145,82,55,.45); background:var(--accent-sft); }
  .rc-badge.b-dim    { color:var(--dim); }
  .rc-badge.b-label  { color:var(--accent); border-color:var(--accent); background:var(--accent-sft); font-weight:600; }
  .run-card .rc-meta { margin-left:auto; font-size:var(--fs-xs); color:var(--dim);
    display:flex; gap:12px; align-items:center; }
  .run-card .rc-hash { font-family:var(--mono); color:var(--faint); }
  .run-card .rc-groups { display:flex; gap:var(--sp-6); flex-wrap:wrap; }
  .run-card .rg-title { font-size:var(--fs-xs); color:var(--accent); text-transform:uppercase;
    letter-spacing:.1em; font-weight:600; margin-bottom:2px; }
  .run-card .rg-items { font-size:var(--fs-sm); color:var(--text); line-height:1.9; }
  .run-card .rg-items .pk { color:var(--dim); }
  .run-card .rg-items .pv { color:var(--white); font-weight:600; font-family:var(--mono); }
"""


def TOP_NAV(active=""):
    """Sticky primary nav, injected as the first body child on every page."""
    tabs = [
        ("config", "/", "Configurator"),
        ("docs", "/docs", "Reference"),
        ("plots", "/plots", "Results"),
        ("compare", "/compare", "Compare"),
        ("conv", "/convergence", "Convergence"),
    ]
    items = []
    for key, href, label in tabs:
        cls = " active" if key == active else ""
        cur = ' aria-current="page"' if key == active else ""
        items.append(f'<a href="{href}" class="tab{cls}"{cur}>{label}</a>')
    return (
        '<a class="skip-link" href="#main">Skip to content</a>'
        '<header class="nav"><nav aria-label="Primary">'
        '<a class="brand" href="/">DAO</a>'
        '<div class="tabs">' + "".join(it.replace('class="tab', 'class="') for it in items) + "</div>"
        "</nav></header>"
    )


FOOTER = (
    '<footer class="site">DAO &nbsp;·&nbsp; Yimin Huang &nbsp;·&nbsp; '
    '<a href="mailto:huangym23@m.fudan.edu.cn">huangym23@m.fudan.edu.cn</a> '
    '&nbsp;·&nbsp; <a href="/docs#references">References</a> '
    '&nbsp;·&nbsp; <a href="/docs#license">MIT License</a></footer>'
)


# ─── HTML template (Configurator) ──────────────────────────────
HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<title>DAO — X-ray Reflection Model</title>
<style>
{{ base_css | safe }}
  /* Editorial opening: the artwork, title and scientific context form one composition. */
  .hero-scene { position:relative; isolation:isolate; overflow:hidden; min-height:334px;
    display:flex; align-items:center; margin:4px 0 24px; padding:38px 48px;
    border-radius:18px; background:#100d0b; box-shadow:0 18px 36px rgba(33,25,18,.12);
    color:#f4f1e9; }
  .hero-scene::before { content:''; position:absolute; inset:0; z-index:-1;
    background:linear-gradient(90deg,#100d0bf7 0,#100d0bea 30%,#100d0b8c 49%,#100d0b00 74%); }
  .hero-scene img { position:absolute; z-index:-2; width:110%; height:100%; right:-8%; top:0;
    object-fit:cover; object-position:center 50%; }
  .hero-copy { position:relative; z-index:1; max-width:475px; }
  .hero-index { display:block; margin-bottom:21px; color:#dbb08b; font:11px var(--mono);
    text-transform:uppercase; letter-spacing:.16em; }
  .hero-scene h1 { font-family:Georgia,'Times New Roman',serif; color:#fcfaf5;
    font-size:clamp(2.25rem,4vw,3.55rem); line-height:1.08; font-weight:400; letter-spacing:-.035em; }
  .hero-scene h1 em { font-style:italic; color:#e7b68c; }
  .hero-scene p { margin-top:18px; max-width:390px; color:#d5c9bc; font-size:.94rem; line-height:1.65; }
  .hero-meta { display:flex; align-items:center; gap:12px; margin-top:28px; color:#d0c1b0;
    font:10px var(--mono); text-transform:uppercase; letter-spacing:.1em; }
  .hero-meta i { width:18px; height:1px; background:#b98662; flex:none; }
  .hero-art-note { position:absolute; right:24px; bottom:18px; font:10px var(--mono);
    color:#e8d7c5; letter-spacing:.02em; text-shadow:0 1px 3px #000; }
  .hero-art-note a { color:inherit; text-decoration:underline; text-underline-offset:3px; }
  .hero-art-note a:hover { color:white; }
  .workspace-intro { display:flex; align-items:baseline; justify-content:space-between; gap:16px;
    margin:0 0 14px; }
  .workspace-intro h2 { color:var(--white); font-family:Georgia,'Times New Roman',serif;
    font-size:1.55rem; font-weight:400; letter-spacing:-.025em; }
  .workspace-intro span { color:var(--dim); font:11px var(--mono); letter-spacing:.04em; }
  .workspace { display:grid; grid-template-columns:minmax(0,1.65fr) minmax(320px,1fr); gap:20px; align-items:start; }
  .parameter-column, .context-column { min-width:0; }
  .context-column { position:sticky; top:78px; }
  .parameter-column .card h2 { margin-bottom:20px; }
  .section-number { font:11px var(--mono); color:#a46642; margin-right:4px; }
  .context-column > .card:first-child { background:#f7f3ed; border-color:transparent; }
  .context-column > .card:last-child { background:#f1ebe3; border-color:transparent; color:var(--text); }
  .context-column > .card:last-child h2 { color:var(--white); }
  .context-column > .card:last-child .cmd-cap,
  .context-column > .card:last-child .cmd-hint { color:var(--dim); }
  .context-column > .card:last-child .cmd-box { background:#fffdfa; border-color:transparent; color:var(--white); }
  .context-column > .card:last-child .cmd-box .prompt { color:var(--accent); }
  .context-column > .card:last-child .btn-primary { background:var(--accent); color:white; }
  .context-column > .card:last-child .btn-primary:hover { background:var(--accent2); }
  .context-column > .card:last-child .btn-secondary { background:#fffdfa; border-color:transparent; color:var(--text); }
  .context-column > .card:last-child .btn-secondary:hover { background:white; }
  .section-sub { font-size:.8125rem; color:var(--dim); margin:-4px 0 12px; }
  .model-description { font-size:.8125rem; color:var(--dim); margin:10px 0 16px; line-height:1.6; }
  .geometry { margin:0; }
  .geometry svg { display:block; width:100%; height:auto; margin:12px 0; }
  .geometry svg text { fill:var(--dim); font:11px var(--mono); }
  .geometry figcaption { color:var(--dim); font-size:.75rem; line-height:1.6; margin-top:10px; }
  .geometry-readout { display:flex; justify-content:space-between; gap:10px; padding:10px 0;
    border-top:1px solid var(--border); border-bottom:1px solid var(--border); font:12px var(--mono); color:var(--text); flex-wrap:wrap; }
  .derived-caption { color:var(--dim); font-size:.8125rem; margin:16px 0 4px; }
  @media(max-width:900px) {
    .workspace { grid-template-columns:1fr; }
    .context-column { position:static; }
    .geometry svg { max-width:440px; margin:12px auto; }
  }
  @media(max-width:760px) {
    .hero-scene { min-height:420px; padding:28px; align-items:flex-start; }
    .hero-copy { max-width:430px; }
    .hero-scene h1 { font-size:clamp(2.1rem,7vw,3rem); }
    .hero-scene img { width:100%; height:78%; top:auto; right:0; bottom:0; object-fit:cover; object-position:50% 54%; }
    .hero-scene::before { background:linear-gradient(180deg,#100d0bf7 0,#100d0be3 38%,#100d0b20 82%); }
    .hero-meta { margin-top:16px; }
    .hero-art-note { right:24px; bottom:18px; }
  }
  @media(max-width:520px) {
    .hero-scene { min-height:390px; padding:24px; }
    .hero-index { margin-bottom:14px; }
    .hero-scene h1 { font-size:2.2rem; }
    .hero-scene p { margin-top:10px; font-size:.875rem; }
    .hero-scene img { height:65%; bottom:0; }
    .hero-meta { display:none; }
    .workspace-intro h2 { font-size:1.35rem; }
    .workspace-intro span { font-size:10px; }
  }

  /* ── Model selector (radiogroup) ────────────── */
  .model-grid { display:flex; flex-wrap:wrap; gap:6px; }
  .model-btn { background:white; border:1px solid var(--border); border-radius:9px;
    padding:9px 12px; cursor:pointer; font-family:inherit; font-size:13px; font-weight:500;
    color:var(--dim); transition:border-color .15s,background .15s; }
  .model-btn:hover { border-color:var(--border-hi); color:var(--white); }
  .model-btn.active { border-color:var(--accent); background:var(--accent-sft); color:var(--accent); }
  .model-btn[disabled] { opacity:.4; cursor:not-allowed; }

  /* ── Parameter fields ───────────────────────── */
  .field { display:grid; grid-template-columns:84px minmax(0,1fr) 92px; align-items:center; column-gap:var(--sp-3); row-gap:4px; margin-bottom:var(--sp-1); }
  .field label { font-size:var(--fs-sm); font-weight:500; color:var(--text); text-align:left; white-space:nowrap; }
  .field input[type=range] { -webkit-appearance:none; appearance:none; width:100%; min-width:0; height:28px;
    background:transparent; cursor:pointer; border-radius:6px; outline:none; }
  .field input[type=range]::-webkit-slider-runnable-track { height:3px; background:#d9d0c5; border-radius:99px; }
  .field input[type=range]::-moz-range-track { height:3px; background:#d9d0c5; border-radius:99px; }
  .field input[type=range]::-webkit-slider-thumb { -webkit-appearance:none; appearance:none;
    width:13px; height:13px; margin-top:-5px; border:2px solid #a8785b; background:white; border-radius:50%;
    box-shadow:0 1px 3px #3d2d2118; transition:box-shadow .15s,border-color .15s; }
  .field input[type=range]::-moz-range-thumb { box-sizing:border-box; width:13px; height:13px;
    border:2px solid #a8785b; background:white; border-radius:50%;
    box-shadow:0 1px 3px #3d2d2118; transition:box-shadow .15s,border-color .15s; }
  .field input[type=range]:hover::-webkit-slider-thumb,
  .field input[type=range]:active::-webkit-slider-thumb { border-color:var(--accent); box-shadow:0 0 0 4px var(--accent-sft); }
  .field input[type=range]:hover::-moz-range-thumb,
  .field input[type=range]:active::-moz-range-thumb { border-color:var(--accent); box-shadow:0 0 0 4px var(--accent-sft); }
  .field input[type=range]:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
  .field .val { appearance:textfield; -moz-appearance:textfield; font-size:14px; color:var(--white); text-align:center; font-weight:500;
    background:#f4eee7; border:1px solid transparent; border-radius:10px; padding:8px;
    width:92px; height:36px; font-family:var(--mono); font-variant-numeric:tabular-nums;
    transition:background .15s,border-color .15s,box-shadow .15s; }
  .field .val::-webkit-inner-spin-button, .field .val::-webkit-outer-spin-button { -webkit-appearance:none; margin:0; }
  .field .val:hover { background:#eee4d9; border-color:#d9d0c5; }
  .field .val:focus { background:white; border-color:#a8785b; outline:none; box-shadow:0 0 0 3px var(--accent-sft); }
  .field .val.flash { border-color:var(--red); box-shadow:0 0 0 2px rgba(239,83,80,.2); }
  .field-desc { font-size:var(--fs-xs); color:var(--dim); grid-column:2/-1; margin:2px 0 14px; line-height:1.5; }
  #coronaParams > div:last-child .field-desc, #slabParams > div:last-child .field-desc,
  #diskParams > div:last-child .field-desc, #testParams > div:last-child .field-desc { margin-bottom:4px; }
  .field-ends { grid-column:2/3; display:flex; justify-content:space-between; font-size:var(--fs-xs);
    color:var(--faint); margin-top:-2px; }
  @media (max-width:520px){
    .field { grid-template-columns:1fr auto; }
    .field label { grid-column:1/-1; text-align:left; }
    .field-ends { grid-column:1/-1; }
    .wrapper { padding:var(--sp-4) var(--sp-3); }
  }

  /* ── Toggles ────────────────────────────────── */
  .toggle-row { display:flex; align-items:center; gap:var(--sp-4); margin-bottom:var(--sp-3); }
  .toggle-row > label:first-child { font-size:var(--fs-sm); font-weight:500; color:var(--text); min-width:120px; text-align:right; }
  .toggle { position:relative; width:44px; height:24px; cursor:pointer; flex:none; }
  .toggle input { opacity:0; width:0; height:0; }
  .toggle .slider { position:absolute; inset:0; background:var(--border); border-radius:12px; border:1px solid var(--border); transition:all .2s; }
  .toggle .slider::before { content:''; position:absolute; width:18px; height:18px; left:2px; top:2px; background:white; border-radius:50%; transition:all .2s; }
  .toggle input:checked+.slider { background:var(--accent-sft); border-color:var(--accent); }
  .toggle input:checked+.slider::before { transform:translateX(20px); background:var(--accent); }
  .toggle input:focus-visible+.slider { box-shadow:0 0 0 3px var(--accent-sft); }

  /* ── Derived panel ──────────────────────────── */
  .divider { height:1px; margin:var(--sp-4) 0; background:linear-gradient(90deg,transparent,var(--border),transparent); border:none; }
  .derived { display:grid; gap:0; }
  .derived .item { display:flex; justify-content:space-between; align-items:baseline; gap:12px; padding:9px 0; border-bottom:1px solid var(--border); }
  .derived .item:last-child { border-bottom:0; }
  .derived .item .dlabel { font-size:12px; color:var(--dim); }
  .derived .item .dval { font:13px var(--mono); color:var(--white); font-variant-numeric:tabular-nums; }

  /* ── select ─────────────────────────────────── */
  .field select { grid-column:2/-1; font-size:var(--fs-sm); color:var(--white); background:var(--bg2);
    border:1px solid var(--border); border-radius:var(--r-sm); padding:7px 12px; font-family:inherit; cursor:pointer; }

  /* ── Command box + sticky bar ───────────────── */
  .cmd-cap { font-size:var(--fs-xs); color:var(--dim); margin-bottom:var(--sp-2); }
  .cmd-box { background:var(--bg); border:1px solid var(--border); border-radius:var(--r-md);
    padding:12px 14px; font-size:.8125rem; color:var(--text); overflow-wrap:anywhere; line-height:1.7; font-family:var(--mono); }
  .cmd-box .prompt { color:var(--dim); user-select:none; }
  .cmd-hint { font-size:var(--fs-xs); color:var(--dim); margin-top:var(--sp-3); line-height:1.5; }
  .cmd-hint code { color:var(--cyan); font-size:.78rem; }
  .actions { display:flex; gap:10px; margin-top:var(--sp-4); flex-wrap:wrap; }

  /* ── details / accordion ────────────────────── */
  details.adv summary { cursor:pointer; font-size:var(--fs-sm); color:var(--dim); list-style:none;
    display:flex; align-items:center; gap:8px; font-weight:500; }
  details.adv summary::-webkit-details-marker { display:none; }
  details.adv summary::before { content:'▸'; color:var(--accent); transition:transform .2s; }
  details.adv[open] summary::before { transform:rotate(90deg); }
  details.adv .body { margin-top:var(--sp-4); }

  /* ── Queue table ────────────────────────────── */
  .queue-table { width:100%; border-collapse:collapse; font-size:var(--fs-sm); margin-top:var(--sp-3); }
  .queue-table th { text-align:left; color:var(--accent); font-weight:600; font-size:var(--fs-xs);
    text-transform:uppercase; letter-spacing:.08em; padding:10px 12px; border-bottom:1px solid var(--border); }
  .queue-table td { padding:10px 12px; border-bottom:1px solid var(--border); color:var(--text); vertical-align:top; }
  .queue-table tr:hover td { background:var(--accent-sft); }
  .queue-table .q-idx { color:var(--dim); font-weight:600; width:40px; text-align:center; }
  .queue-table .q-model { color:var(--cyan); font-weight:600; width:100px; }
  .queue-table .q-params { color:var(--text); font-size:.78rem; word-break:break-all; font-family:var(--mono); }
  .queue-table .q-actions { width:150px; text-align:left; white-space:nowrap; }
  .q-btn { background:var(--bg2); border:1px solid var(--border); cursor:pointer; font-size:var(--fs-xs);
    padding:4px 10px; border-radius:var(--r-sm); transition:all .15s; color:var(--text); margin-left:4px; }
  .q-btn:hover { border-color:var(--border-hi); color:var(--white); }
  .q-btn-del:hover { color:var(--red); border-color:rgba(239,83,80,.4); }
  .queue-empty { text-align:center; color:var(--dim); font-size:var(--fs-sm); padding:28px 0; }
  .batch-actions { display:flex; gap:10px; margin-top:var(--sp-4); flex-wrap:wrap; align-items:center; }
  .batch-count { font-size:var(--fs-sm); color:var(--dim); margin-left:auto; font-weight:500; }
  .grid > .card { min-width:0; }
  .toggle-row { flex-wrap:wrap; }
  .toggle-row > label:first-child { min-width:0; text-align:left; }
  #queueBody { overflow-x:auto; }
</style>
</head>
<body>
{{ nav | safe }}
<main id="main" class="wrapper">

  <section class="hero-scene" aria-labelledby="heroTitle">
    <img src="/image/nasa-black-hole-accretion-disk.png" width="1920" height="1080"
         alt="NASA visualization of a black hole distorting light from its thin, orange accretion disk">
    <div class="hero-copy">
      <span class="hero-index">DAO / reflection spectroscopy</span>
      <h1 id="heroTitle">Spectrum from an<br><em>irradiated disk.</em></h1>
      <p>Explore how coronal radiation interacts with a photoionised slab and emerges as an X-ray reflection spectrum.</p>
      <div class="hero-meta"><span>Compton transfer</span><i aria-hidden="true"></i><span>Photoionisation</span></div>
    </div>
    <span class="hero-art-note">Visualization: <a href="https://svs.gsfc.nasa.gov/14146" target="_blank" rel="noopener noreferrer">NASA's Goddard Space Flight Center / Jeremy Schnittman</a></span>
  </section>

  <div class="workspace-intro"><h2>Model configuration</h2><span>SET PARAMETERS / GENERATE COMMAND</span></div>
  <div class="workspace">
  <div class="parameter-column">
  <div class="card full">
    <h2><span class="section-number">01</span> Incident spectrum</h2>
    <div class="model-grid" id="modelGrid" role="radiogroup" aria-label="Corona spectrum model"></div>
    <div class="model-description" id="modelDescription"></div>
    <div id="coronaParams"></div>
  </div>

  <div class="card">
      <h2><span class="section-number">02</span> Slab parameters</h2>
      <div id="slabParams"></div>
      <div class="toggle-row" style="margin-top:6px">
        <label for="angsca">Angle-dependent scattering (<code>-angsca</code>)</label>
        <label class="toggle"><input type="checkbox" id="angsca" checked aria-label="Angle-dependent scattering" onchange="updateCmd();persist()"><span class="slider"></span></label>
      </div>
      <div class="field-desc">On = angle-dependent kernel; off = angle-averaged kernel.</div>

    </div>
    <div class="card">
      <h2><span class="section-number">03</span> Disk emission</h2>
      <div id="diskParams"></div>

    </div>

  <div class="card full">
    <h2><span class="section-number">04</span> Solver options</h2>
    <details class="adv">
      <summary>Advanced · compPS benchmark</summary>
      <div class="body">
        <div class="toggle-row">
          <label for="testRt">Enable test mode</label>
          <label class="toggle"><input type="checkbox" id="testRt" aria-label="Enable compPS test mode" onchange="toggleTest()"><span class="slider"></span></label>
        </div>
        <div class="field-desc">Isothermal pure-scattering slab vs Xspec compPS (Poutanen &amp; Svensson 1996),
          illuminated by a bottom blackbody seed. kT<sub>e</sub> is the slab temperature; the corona is set to
          <code>blackbody</code> automatically to supply the seed.</div>
        <div id="testParams" style="display:none"></div>
      </div>
    </details>
  </div>

  </div>
  <aside class="context-column" aria-label="Geometry and model configuration">
  <div class="card">
    <h2>Illumination geometry</h2>
    <p class="section-sub">Incident direction relative to the slab normal</p>
      <figure class="geometry">
        <svg viewBox="0 0 420 250" role="img" aria-labelledby="geometryTitle geometryDesc">
          <title id="geometryTitle">Incident radiation on a plane-parallel slab</title>
          <desc id="geometryDesc">The incident ray follows the selected angle relative to the surface normal. Brown rays represent illustrative outgoing radiation, not computed trajectories.</desc>
          <defs>
                    <marker id="incidentArrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M 0 0 L 10 5 L 0 10 Z" fill="#a57644"/></marker>
            <marker id="outgoingArrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M 0 0 L 10 5 L 0 10 Z" fill="#9a6b4c"/></marker>
            <pattern id="slabGrid" width="20" height="16" patternUnits="userSpaceOnUse"><path d="M20 0H0V16" fill="none" stroke="#a59685" stroke-opacity=".12"/></pattern>
          </defs>
          <path d="M30 172H390V230H30Z" fill="#eee6db"/>
          <path d="M30 172H390V230H30Z" fill="url(#slabGrid)"/>
          <path d="M30 172H390" stroke="#a59685"/>
          <path d="M210 172V22" stroke="#978979" stroke-dasharray="4 5"/>
          <text x="220" y="30">normal</text>
          <path id="incidentRay" d="M104 66L210 172" fill="none" stroke="#a57644" stroke-width="2" marker-end="url(#incidentArrow)"/>
          <circle id="coronaSource" cx="104" cy="66" r="3" fill="#a57644"/>
          <text x="32" y="22" fill="#a57644">incident photons</text>
          <path id="angleArc" d="M210 130A42 42 0 0 0 180 142" fill="none" stroke="#a57644"/>
          <text id="angleLabel" x="177" y="118">θ</text>
          <path d="M210 172L326 64 M210 172L370 123" stroke="#9a6b4c" stroke-opacity=".7" fill="none" marker-end="url(#outgoingArrow)"/>
          <text x="288" y="48">outgoing</text>
          <circle cx="210" cy="172" r="5" fill="#9a6b4c"/>
          <text x="46" y="207">illuminated slab</text>
          <text x="298" y="207">n<tspan baseline-shift="sub" font-size="8">H</tspan> · ξ · A<tspan baseline-shift="sub" font-size="8">Fe</tspan></text>
        </svg>
        <div class="geometry-readout"><span id="geometryAngle">θ = 45.0°</span><span id="geometryDisk"></span></div>
        <figcaption>Plane-parallel geometry · θ = arccos μ.<br>Outgoing directions are schematic. The solver snaps μ to its angular quadrature.</figcaption>
      </figure>
    <div class="derived-caption">Derived quantities</div>
    <div class="derived" id="derivedBlock"></div>
  </div>
  <div class="card full">
    <h2>Run configuration</h2>
    <label for="runLabel" style="display:block;font-size:var(--fs-sm);margin-bottom:6px;">Run label (optional)</label>
    <input type="text" id="runLabel" placeholder="e.g. Iron abundance comparison"
           style="width:100%;box-sizing:border-box;padding:10px;border:1px solid var(--border);border-radius:8px;background:white;color:var(--text);font:inherit;"
           aria-describedby="runLabelHint" oninput="updateCmd();persist()">
    <p class="cmd-hint" id="runLabelHint">Saved in RUN.txt and params.json. Labels do not change the result folder;
      rerunning the same physics replaces that folder’s label and summary.</p>
    <div class="cmd-cap">Run this in your terminal:</div>
    <div class="cmd-box" id="cmdBox"><span class="prompt">$ </span><span id="cmdText"></span></div>
    <div class="actions">
      <button class="btn btn-primary" onclick="copyCmd()">Copy Command</button>
      <button class="btn btn-secondary" onclick="addToQueue()">Add to Queue</button>
      <button class="btn btn-secondary" onclick="resetAll()">Reset</button>
    </div>
    <div class="cmd-hint">
      Make sure your runtime environment is set before running
      (e.g. <code>source $HEADAS/headas-init.sh</code>).
    </div>
    <span id="barCmd" hidden></span>
  </div>

  </aside>
  </div>

  <div class="card full">
    <h2>Run queue</h2>
    <div id="queueBody"></div>
    <div class="batch-actions">
      <button class="btn btn-primary" onclick="copyAllCmds()">Copy All</button>
      <button class="btn btn-secondary" onclick="exportScript()">Export .sh</button>
      <button class="btn btn-secondary" onclick="clearQueue()">Clear Queue</button>
      <span class="batch-count" id="queueCount"></span>
    </div>
  </div>

</main>
{{ footer | safe }}

<div class="toast" id="toast"></div>

<script>
// ── Data from server ──────────────────────
const MODELS = {{ models_json | safe }};
const SLAB   = {{ slab_json | safe }};
const DISK   = {{ disk_json | safe }};
// compps benchmark slab params: kT_e is the slab temperature (emitted as -kT_e),
// tau the vertical Thomson depth. Distinct field name 'kTe_slab' avoids a DOM-id
// clash with a corona that also exposes kT_e.
const TEST   = [
  {"name":"kTe_slab","flag":"kT_e","label":"kT<sub>e</sub> (slab)","desc":"Slab uniform temperature [keV]","default":60,"min":1,"max":500,"step":0.5,"scale":"log"},
  {"name":"tau","label":"τ<sub>slab</sub>","desc":"Slab vertical Thomson optical depth","default":0.5,"min":0.01,"max":5,"step":0.01},
];
const DEFAULTS = {nh:15,zeta:3,frac:-1,incidence:0.7071067811865476,Afe:1,kT_disk:0.35,kTe_slab:60,tau:0.5};
const TIPS = {
  Ecut:"≈ 2–3 kT_e",
  E_low_cut:"Suppresses the IR divergence of a bare power law",
  incidence:"Snaps to the nearest Gauss–Legendre quadrature node at runtime",
  Afe:"Scales the 6.4–6.97 keV Fe K complex",
  zeta:"ζ = log₁₀(ξ); ξ = 10^ζ erg cm s⁻¹",
  nh:"n_H = 10^nh cm⁻³",
};

const DEFAULT_CORONA = "{{ default_corona }}";
const SKEY = "dao.session.v1";
let corona = DEFAULT_CORONA;
let vals = {};

// ── helpers ───────────────────────────────
function paramByName(name) {
  for (const k in MODELS) for (const p of MODELS[k].params) if (p.name === name) return p;
  for (const p of SLAB) if (p.name === name) return p;
  for (const p of DISK) if (p.name === name) return p;
  for (const p of TEST) if (p.name === name) return p;
  return null;
}
function stripHtml(s){ return s.replace(/<[^>]+>/g,'').trim(); }
function fmtNum(v) {
  if (Math.abs(v) >= 1e4 || (Math.abs(v) < 0.01 && v !== 0)) return v.toExponential(2);
  return parseFloat(v.toPrecision(6)).toString();
}

// ── Build fields ──────────────────────────
function buildField(p, container) {
  const id = 'f_' + p.name;
  const isLog = p.scale === 'log';
  const aria = stripHtml(p.label);
  const tip = TIPS[p.name] ? ` title="${TIPS[p.name]}"` : '';
  let rmin = p.min, rmax = p.max, rstep = p.step, rval = p.default;
  if (isLog) {
    rmin = Math.log10(p.min); rmax = Math.log10(p.max);
    rstep = (rmax - rmin) / 200; rval = Math.log10(p.default);
  }
  const d = document.createElement('div');
  d.innerHTML = `
    <div class="field">
      <label for="${id}"${tip}>${p.label}</label>
      <input type="range" id="${id}" data-log="${isLog?1:0}" min="${rmin}" max="${rmax}" step="${rstep}" value="${rval}"
             aria-label="${aria}" oninput="syncVal('${p.name}','${id}')">
      <input type="number" class="val" id="${id}_v" min="${p.min}" max="${p.max}" step="${p.step}"
             value="${fmtNum(p.default)}" inputmode="decimal" aria-label="${aria} value"
             onchange="syncSlider('${p.name}','${id}')">
      <div class="field-ends"><span>${fmtNum(p.min)}</span><span>${fmtNum(p.max)}</span></div>
    </div>
    <div class="field-desc">${p.desc}</div>`;
  container.appendChild(d);
  vals[p.name] = p.default;
}

function setFieldValue(name, v) {
  const id = 'f_' + name;
  const r = document.getElementById(id);
  const box = document.getElementById(id + '_v');
  if (!r) return;
  const isLog = r.dataset.log === '1';
  r.value = isLog ? Math.log10(v) : v;
  if (box) box.value = fmtNum(v);
}

function syncVal(name, id) {
  const r = document.getElementById(id);
  const isLog = r.dataset.log === '1';
  let v = parseFloat(r.value);
  if (isLog) v = Math.pow(10, v);
  document.getElementById(id+'_v').value = fmtNum(v);
  vals[name] = v;
  updateCmd(); updateDerived(); persist();
}

function syncSlider(name, id) {
  const box = document.getElementById(id+'_v');
  const r = document.getElementById(id);
  const p = paramByName(name);
  let v = parseFloat(box.value);
  if (isNaN(v) || !p) return;
  const clamped = Math.min(p.max, Math.max(p.min, v));
  if (clamped !== v) {
    box.classList.add('flash');
    setTimeout(() => box.classList.remove('flash'), 600);
    v = clamped;
  }
  box.value = fmtNum(v);
  r.value = (r.dataset.log === '1') ? Math.log10(v) : v;
  vals[name] = v;
  updateCmd(); updateDerived(); persist();
}

// ── Corona model selector (radiogroup) ────
function buildModels() {
  const grid = document.getElementById('modelGrid');
  grid.innerHTML = '';
  document.getElementById('modelDescription').innerHTML = MODELS[corona].desc;
  const testOn = document.getElementById('testRt') && document.getElementById('testRt').checked;
  Object.entries(MODELS).forEach(([key, m]) => {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.setAttribute('role', 'radio');
    const isActive = key === corona;
    btn.setAttribute('aria-checked', isActive ? 'true' : 'false');
    btn.tabIndex = isActive ? 0 : -1;
    btn.className = 'model-btn' + (isActive ? ' active' : '');
    btn.dataset.key = key;
    if (testOn && key !== 'blackbody') btn.disabled = true;
    btn.textContent = m.label;
    btn.title = stripHtml(m.desc);
    btn.onclick = () => selectModel(key);
    btn.onkeydown = (e) => onModelKey(e, key);
    grid.appendChild(btn);
  });
}

function onModelKey(e, key) {
  const keys = Object.keys(MODELS);
  let i = keys.indexOf(key);
  if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); selectModel(key); }
  else if (e.key === 'ArrowRight' || e.key === 'ArrowDown') {
    e.preventDefault(); const nk = keys[(i + 1) % keys.length]; selectModel(nk); focusModel(nk);
  } else if (e.key === 'ArrowLeft' || e.key === 'ArrowUp') {
    e.preventDefault(); const nk = keys[(i - 1 + keys.length) % keys.length]; selectModel(nk); focusModel(nk);
  }
}
function focusModel(key) {
  const el = document.querySelector(`.model-btn[data-key="${key}"]`);
  if (el && !el.disabled) el.focus();
}

function selectModel(key) {
  corona = key;
  Object.keys(MODELS).forEach(k => MODELS[k].params.forEach(p => delete vals[p.name]));
  buildModels();
  buildCoronaParams();
  updateCmd(); updateDerived(); persist();
}

function buildCoronaParams() {
  const el = document.getElementById('coronaParams');
  el.innerHTML = '';
  MODELS[corona].params.forEach(p => buildField(p, el));
}

// ── Derived quantities ────────────────────
function updateDerived() {
  updateGeometry();
  const nh = vals.nh ?? 15, zeta = vals.zeta ?? 3;
  const xi = Math.pow(10, zeta);
  const nH = Math.pow(10, nh);
  const J = xi * nH / Math.pow(4 * Math.PI, 2);
  document.getElementById('derivedBlock').innerHTML = `
    <div class="item"><div class="dlabel">ξ [erg cm s⁻¹]</div><div class="dval">${xi.toExponential(2)}</div></div>
    <div class="item"><div class="dlabel">n<sub>H</sub> [cm⁻³]</div><div class="dval">${nH.toExponential(2)}</div></div>
    <div class="item" title="J = ξ·n_H/(4π)² — mean intensity that normalises the corona+disk illumination"><div class="dlabel">J [erg cm⁻² s⁻¹]</div><div class="dval">${J.toExponential(2)}</div></div>
  `;
}

// Diagram uses the requested incidence cosine; outgoing rays are illustrative.
function updateGeometry() {
  const mu = Math.min(1, Math.max(0.01, vals.incidence ?? DEFAULTS.incidence));
  const theta = Math.acos(mu), sinTheta = Math.sin(theta);
  const x = 210 - 150 * sinTheta, y = 172 - 150 * mu;
  const ray = `M${x} ${y}L210 172`;
  document.getElementById('incidentRay').setAttribute('d', ray);
  document.getElementById('coronaSource').setAttribute('cx', x);
  document.getElementById('coronaSource').setAttribute('cy', y);
  document.getElementById('angleArc').setAttribute('d',
    theta < 1e-8 ? '' : `M210 130A42 42 0 0 0 ${210 - 42 * sinTheta} ${172 - 42 * mu}`);
  document.getElementById('angleLabel').setAttribute('x', 204 - 57 * Math.sin(theta / 2));
  document.getElementById('angleLabel').setAttribute('y', 172 - 57 * Math.cos(theta / 2));
  document.getElementById('geometryAngle').textContent = `θ = ${(theta * 180 / Math.PI).toFixed(1)}° · μ = ${mu.toFixed(3)}`;
  const diskOn = (vals.frac ?? DEFAULTS.frac) > 0;
  document.getElementById('geometryDisk').textContent = diskOn
    ? `kTdisk = ${fmtNum(vals.kT_disk ?? DEFAULTS.kT_disk)} keV`
    : 'Disk illumination off';
}

// ── Test mode toggle ──────────────────────
function toggleTest() {
  const on = document.getElementById('testRt').checked;
  document.getElementById('testParams').style.display = on ? 'block' : 'none';
  if (on && corona !== 'blackbody') {
    selectModel('blackbody');
    showToast('Test mode uses a blackbody seed — corona set to blackbody.', 'info');
  }
  buildModels();
  updateCmd(); persist();
}

// ── Command builder ───────────────────────
function shellQuote(value) {
  return "'" + value.replaceAll("'", "'\"'\"'") + "'";
}
function escapeHtml(value) {
  return String(value).replaceAll('&', '&amp;').replaceAll('<', '&lt;')
    .replaceAll('>', '&gt;').replaceAll('"', '&quot;').replaceAll("'", '&#39;');
}
function updateCmd() {
  let parts = ['./maindaocl', '-corona ' + corona];
  MODELS[corona].params.forEach(p => {
    parts.push('-' + p.name + ' ' + fmtNum(vals[p.name] ?? p.default));
  });
  [SLAB, DISK].forEach(group => group.forEach(p => {
    const v = vals[p.name] ?? p.default;
    if (v !== DEFAULTS[p.name]) parts.push('-' + p.name + ' ' + fmtNum(v));
  }));
  if (!document.getElementById('angsca').checked) parts.push('-angsca 0');
  if (document.getElementById('testRt').checked) {
    parts.push('-test_rt compps');
    TEST.forEach(p => {
      const v = vals[p.name] ?? p.default;
      parts.push('-' + (p.flag || p.name) + ' ' + fmtNum(v));
    });
  }
  const label = document.getElementById('runLabel').value;
  if (label) parts.push('-label ' + shellQuote(label));
  const cmd = parts.join(' ');
  document.getElementById('cmdText').textContent = cmd;
  document.getElementById('barCmd').textContent = '$ ' + cmd;
}

// ── Clipboard / actions ───────────────────
function clipCopy(text) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(() => showToast('Copied to clipboard')).catch(() => fallbackCopy(text));
  } else { fallbackCopy(text); }
}
function fallbackCopy(text) {
  const ta = document.createElement('textarea');
  ta.value = text; ta.style.cssText = 'position:fixed;left:-9999px';
  document.body.appendChild(ta); ta.select();
  document.execCommand('copy'); document.body.removeChild(ta);
  showToast('Copied to clipboard');
}
function copyCmd() { clipCopy(document.getElementById('cmdText').textContent); }

function resetAll() {
  localStorage.removeItem(SKEY);
  document.getElementById('runLabel').value = '';
  corona = DEFAULT_CORONA;
  vals = {};
  document.getElementById('angsca').checked = true;
  document.getElementById('testRt').checked = false;
  buildModels();
  buildCoronaParams();
  ['slabParams','diskParams','testParams'].forEach(id => document.getElementById(id).innerHTML = '');
  SLAB.forEach(p => buildField(p, document.getElementById('slabParams')));
  DISK.forEach(p => buildField(p, document.getElementById('diskParams')));
  TEST.forEach(p => buildField(p, document.getElementById('testParams')));
  toggleTest();
  updateDerived(); updateCmd(); persist();
  showToast('Reset to defaults');
}

function showToast(msg, type='success') {
  const t = document.getElementById('toast');
  const glyph = type === 'error' ? '✕' : type === 'info' ? 'ℹ' : '✓';
  t.className = 'toast ' + type;
  t.textContent = glyph + '  ' + msg;
  t.classList.add('show');
  setTimeout(() => t.classList.remove('show'), 2400);
}

// ── localStorage persistence ──────────────
function persist() {
  try {
    localStorage.setItem(SKEY, JSON.stringify({
      corona, vals, label: document.getElementById('runLabel').value,
      angsca: document.getElementById('angsca').checked,
      testRt: document.getElementById('testRt').checked,
      queue,
    }));
  } catch(e){}
}

// ── Batch queue ──────────────────────────
let queue = [];
function getCmd() { return document.getElementById('cmdText').textContent; }
function getSnapshot() {
  return { corona, vals: {...vals}, label: document.getElementById('runLabel').value,
           testRt: document.getElementById('testRt').checked,
           angsca: document.getElementById('angsca').checked };
}
function cmdSummary(cmd) {
  const m = cmd.match(/-corona\s+(\S+)/);
  const model = m ? m[1] : '?';
  const rest = cmd.replace('./maindaocl', '').replace(/-corona\s+\S+/, '').trim();
  return { model, rest };
}
function addToQueue() {
  const cmd = getCmd();
  if (queue.some(q => q.cmd === cmd)) { showToast('Already in queue', 'info'); return; }
  queue.push({ cmd, snap: getSnapshot() });
  renderQueue(); persist();
  showToast(`Added to queue (${queue.length} total)`);
}
function removeFromQueue(idx) { queue.splice(idx, 1); renderQueue(); persist(); }
function loadFromQueue(idx) {
  const snap = queue[idx].snap;
  corona = snap.corona;
  vals = {...snap.vals};
  document.getElementById('runLabel').value = snap.label ?? '';
  document.getElementById('testRt').checked = snap.testRt;
  if (snap.angsca !== undefined) document.getElementById('angsca').checked = snap.angsca;
  document.getElementById('testParams').style.display = snap.testRt ? 'block' : 'none';
  buildModels();
  const cp = document.getElementById('coronaParams'); cp.innerHTML = '';
  MODELS[corona].params.forEach(p => { buildField(p, cp); setFieldValue(p.name, vals[p.name] ?? p.default); });
  ['slabParams','diskParams','testParams'].forEach(sec => {
    const el = document.getElementById(sec); el.innerHTML = '';
    const list = sec === 'slabParams' ? SLAB : sec === 'diskParams' ? DISK : TEST;
    list.forEach(p => { buildField(p, el); setFieldValue(p.name, vals[p.name] ?? p.default); });
  });
  updateDerived(); updateCmd(); persist();
  showToast(`Loaded run #${idx + 1}`, 'info');
}
function renderQueue() {
  const body = document.getElementById('queueBody');
  const count = document.getElementById('queueCount');
  count.textContent = queue.length > 0 ? `${queue.length} run${queue.length > 1 ? 's' : ''} queued` : '';
  if (queue.length === 0) {
    body.innerHTML = '<div class="queue-empty">Queue is empty — configure parameters above and click "Add to Queue"</div>';
    return;
  }
  let html = `<table class="queue-table"><thead><tr>
    <th>#</th><th>Model</th><th>Parameters</th><th></th>
  </tr></thead><tbody>`;
  queue.forEach((item, i) => {
    const s = cmdSummary(item.cmd);
    html += `<tr>
      <td class="q-idx">${i + 1}</td>
      <td class="q-model">${escapeHtml(s.model)}</td>
      <td class="q-params">${escapeHtml(s.rest)}</td>
      <td class="q-actions">
        <button class="q-btn q-btn-load" onclick="loadFromQueue(${i})" aria-label="Load run ${i+1} into editor">↩ Load</button>
        <button class="q-btn q-btn-del" onclick="removeFromQueue(${i})" aria-label="Remove run ${i+1}">✕ Remove</button>
      </td>
    </tr>`;
  });
  html += '</tbody></table>';
  body.innerHTML = html;
}
function copyAllCmds() {
  if (queue.length === 0) { showToast('Queue is empty', 'error'); return; }
  clipCopy(queue.map(q => q.cmd).join('\n'));
}
function exportScript() {
  if (queue.length === 0) { showToast('Queue is empty', 'error'); return; }
  let script = '#!/bin/bash\n# DAO batch run — ' + new Date().toISOString().slice(0,10) + '\nset -e\n\n';
  script += '# Set your runtime environment before running, e.g.:\n';
  script += '# export HEADAS=/path/to/heasoft/arch\n';
  script += '# source $HEADAS/headas-init.sh\n\n';
  queue.forEach((q, i) => { script += `echo "=== Run ${i+1}/${queue.length} ==="\n${q.cmd}\n\n`; });
  const blob = new Blob([script], {type: 'text/x-shellscript'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = 'dao_batch.sh'; a.click();
  URL.revokeObjectURL(a.href);
  showToast('Exported dao_batch.sh');
}
function clearQueue() { queue = []; renderQueue(); persist(); showToast('Queue cleared', 'info'); }

// ── Init ──────────────────────────────────
function rehydrate() {
  try {
    const raw = localStorage.getItem(SKEY);
    if (!raw) return false;
    const s = JSON.parse(raw);
    if (s.corona && MODELS[s.corona]) corona = s.corona;
    if (s.vals) vals = {...s.vals};
    document.getElementById('runLabel').value = s.label ?? '';
    if (Array.isArray(s.queue)) queue = s.queue;
    document.getElementById('angsca').checked = s.angsca !== false;
    document.getElementById('testRt').checked = !!s.testRt;
    return true;
  } catch(e){ return false; }
}

function init() {
  const restored = rehydrate();
  buildModels();
  const cp = document.getElementById('coronaParams'); cp.innerHTML = '';
  MODELS[corona].params.forEach(p => { buildField(p, cp); if (restored) setFieldValue(p.name, vals[p.name] ?? p.default); });
  ['slabParams','diskParams','testParams'].forEach(id => document.getElementById(id).innerHTML = '');
  SLAB.forEach(p => { buildField(p, document.getElementById('slabParams')); if (restored) setFieldValue(p.name, vals[p.name] ?? p.default); });
  DISK.forEach(p => { buildField(p, document.getElementById('diskParams')); if (restored) setFieldValue(p.name, vals[p.name] ?? p.default); });
  TEST.forEach(p => { buildField(p, document.getElementById('testParams')); if (restored) setFieldValue(p.name, vals[p.name] ?? p.default); });
  document.getElementById('testParams').style.display = document.getElementById('testRt').checked ? 'block' : 'none';
  updateDerived(); updateCmd(); renderQueue();
}

init();
</script>
</body>
</html>
"""


# ─── Flask routes ───────────────────────────────────────────────
@app.route("/")
def index():
    return render_template_string(
        HTML,
        base_css=BASE_CSS,
        nav=TOP_NAV("config"),
        footer=FOOTER,
        models_json=json.dumps(CORONA_MODELS),
        slab_json=json.dumps(SLAB_PARAMS),
        disk_json=json.dumps(DISK_PARAMS),
        default_corona="nthcomp",
    )


# ─── Docs page ──────────────────────────────────────────────────
DOCS_HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<title>DAO — Parameter Reference</title>
<style>
{{ base_css | safe }}
  /* docs reading column — deliberately narrower (documented exception) */
  main.wrapper { max-width:820px; }
  .section { margin-bottom:var(--sp-6); }
  .section h2 {
    font-size:1rem; text-transform:none; letter-spacing:0;
    color:var(--white); margin-bottom:var(--sp-4); padding-bottom:10px;
    border-bottom:1px solid var(--border); display:flex; align-items:center; gap:10px;
  }
  .section h2 .tick { display:none; }
  table { width:100%; border-collapse:collapse; margin-bottom:8px; }
  th { text-align:left; font-size:var(--fs-xs); text-transform:uppercase; letter-spacing:.08em;
    color:var(--dim); padding:8px 12px; border-bottom:1px solid var(--border); }
  td { padding:10px 12px; border-bottom:1px solid var(--border); font-size:var(--fs-base); vertical-align:top; }
  tr:hover td { background:var(--accent-sft); }
  .p-flag { color:var(--cyan); font-weight:600; white-space:nowrap; font-family:var(--mono); }
  .p-type { color:var(--dim); font-size:var(--fs-xs); }
  .p-default { color:var(--white); font-weight:500; }
  .p-desc { color:var(--text); line-height:1.5; }
  .p-note { color:var(--dim); font-size:var(--fs-xs); font-style:italic; }
  .formula-box { padding:14px 18px; margin:10px 0 16px; font-size:.9rem; color:var(--white); line-height:1.8; }
  .formula-box .flabel { color:var(--cyan); font-weight:600; margin-right:8px; font-family:var(--mono); }
  .note-box { background:var(--accent-sft); border-left:3px solid var(--accent);
    padding:12px 16px; border-radius:0 var(--r-md) var(--r-md) 0; margin:14px 0;
    font-size:var(--fs-sm); color:var(--text); line-height:1.6; box-shadow:none; border-top:0; border-right:0; border-bottom:0; }
  .note-box strong { color:var(--accent); }
  .tag { display:inline-block; font-size:.75rem; padding:2px 6px; border-radius:4px;
    font-weight:600; letter-spacing:.04em; vertical-align:middle; margin-left:4px; }
  .tag-required { background:rgba(239,83,80,.15); color:var(--red); }
  .tag-optional { background:rgba(145,82,55,.1); color:var(--cyan); }
  .tag-model { background:var(--accent-sft); color:var(--accent); }
  p.lead { font-size:var(--fs-base); color:var(--text); margin-bottom:14px; line-height:1.6; }
</style>
</head>
<body>
{{ nav | safe }}
<main id="main" class="wrapper">

  <div class="page-head">
    <h1>Parameter Reference</h1>
    <div class="sub">DAO — X-ray Reflection Spectroscopy Model &nbsp;·&nbsp; Yimin Huang &nbsp;·&nbsp; huangym23@m.fudan.edu.cn</div>
  </div>

  <!-- ── Corona Models ────────────────────── -->
  <div class="section" id="corona">
    <h2><span class="tick" aria-hidden="true"></span>Corona Spectrum Models</h2>
    <p class="lead">
      The corona model defines the illuminating X-ray continuum incident on the slab.
      Select with <span class="p-flag">-corona &lt;model&gt;</span>.
      Each model requires specific parameters listed below.
    </p>
    <div class="formula-box">
      <div><span class="flabel">powerlaw</span>
        N(E) ∝ E<sup>−Γ</sup> · exp(−E<sub>lo</sub> / E)  [photons cm<sup>−2</sup> s<sup>−1</sup> keV<sup>−1</sup>]</div>
      <div><span class="flabel">cutoffpl</span>
        N(E) ∝ E<sup>−Γ</sup> · exp(−E / E<sub>cut</sub>) · exp(−E<sub>lo</sub> / E)  [photons cm<sup>−2</sup> s<sup>−1</sup> keV<sup>−1</sup>]</div>
      <div><span class="flabel">nthcomp</span>
        Thermal Comptonisation — Zdziarski, Johnson &amp; Magdziarz (1996).
        Seed photons at kT<sub>bb</sub> Comptonised by electrons at kT<sub>e</sub>.</div>
      <div><span class="flabel">comptt</span>
        Comptonisation model — Titarchuk (1994).
        Wien-law seed photons, plasma temperature kT<sub>e</sub>, optical depth τ<sub>p</sub>.</div>
      <div><span class="flabel">blackbody</span>
        B(E) = (2E<sup>3</sup> / h<sup>2</sup>c<sup>2</sup>) · 1 / [exp(E / kT) − 1]</div>
    </div>
    <table>
      <tr><th>Flag</th><th>Type</th><th>Default</th><th>Models</th><th>Description</th></tr>
      <tr>
        <td class="p-flag">-corona</td><td class="p-type">string</td><td class="p-default">—</td>
        <td><span class="tag tag-required">REQUIRED</span></td>
        <td class="p-desc">Corona model name. One of: <code>powerlaw</code>, <code>cutoffpl</code>,
          <code>nthcomp</code>, <code>comptt</code>, <code>blackbody</code>.</td>
      </tr>
      <tr>
        <td class="p-flag">-Gamma</td><td class="p-type">float</td><td class="p-default">2.0</td>
        <td><span class="tag tag-model">powerlaw</span><span class="tag tag-model">cutoffpl</span><span class="tag tag-model">nthcomp</span></td>
        <td class="p-desc">Photon index Γ (photon-number spectrum, N(E) ∝ E<sup>−Γ</sup>). Typical range: 1.4–3.0 for AGN.</td>
      </tr>
      <tr>
        <td class="p-flag">-Ecut</td><td class="p-type">float</td><td class="p-default">300</td>
        <td><span class="tag tag-model">cutoffpl</span></td>
        <td class="p-desc">High-energy exponential cutoff E<sub>cut</sub> [keV].
          Related to the coronal electron temperature: E<sub>cut</sub> ≈ 2–3 kT<sub>e</sub>.</td>
      </tr>
      <tr>
        <td class="p-flag">-E_low_cut</td><td class="p-type">float</td><td class="p-default">0.1</td>
        <td><span class="tag tag-model">powerlaw</span><span class="tag tag-model">cutoffpl</span></td>
        <td class="p-desc">Low-energy exponential cutoff E<sub>lo</sub> [keV].
          Suppresses the spectrum below this energy as exp(−E<sub>lo</sub>/E) to
          avoid the unphysical infrared divergence of a bare power law.</td>
      </tr>
      <tr>
        <td class="p-flag">-kT_e</td><td class="p-type">float</td><td class="p-default">60 (nthcomp) / 50 (comptt)</td>
        <td><span class="tag tag-model">nthcomp</span><span class="tag tag-model">comptt</span></td>
        <td class="p-desc">Electron / plasma temperature [keV] of the Comptonising corona.</td>
      </tr>
      <tr>
        <td class="p-flag">-kT_bb</td><td class="p-type">float</td><td class="p-default">0.1 / 0.05</td>
        <td><span class="tag tag-model">nthcomp</span><span class="tag tag-model">comptt</span><span class="tag tag-model">blackbody</span></td>
        <td class="p-desc">Seed photon / blackbody temperature [keV]. For nthcomp and comptt this is the
          soft photon field entering the corona; for blackbody it is the emission temperature.</td>
      </tr>
      <tr>
        <td class="p-flag">-taup</td><td class="p-type">float</td><td class="p-default">1.0</td>
        <td><span class="tag tag-model">comptt</span></td>
        <td class="p-desc">Plasma optical depth τ<sub>p</sub> of the Comptonising medium.</td>
      </tr>
    </table>
  </div>

  <!-- ── Slab Physics ─────────────────────── -->
  <div class="section" id="slab">
    <h2><span class="tick" aria-hidden="true"></span>Slab Physics</h2>
    <p class="lead">
      These parameters define the physical conditions of the reflecting slab
      (the accretion disk atmosphere).
    </p>
    <table>
      <tr><th>Flag</th><th>Type</th><th>Default</th><th>Description</th></tr>
      <tr>
        <td class="p-flag">-nh</td><td class="p-type">float</td><td class="p-default">15</td>
        <td class="p-desc">Logarithmic hydrogen number density: n<sub>H</sub> = 10<sup>nh</sup> cm<sup>−3</sup>.
          Typical values: 15 (standard AGN disk) to 22 (very dense inner disk).
          <div class="p-note">Cloudy is called at each depth point with this density.</div></td>
      </tr>
      <tr>
        <td class="p-flag">-zeta</td><td class="p-type">float</td><td class="p-default">3.0</td>
        <td class="p-desc">Ionisation parameter exponent. ζ = log<sub>10</sub>(ξ); ξ = 10<sup>ζ</sup> erg cm s<sup>−1</sup>.
          Controls the ionisation state of the slab. Low ξ → neutral (cold reflection);
          high ξ → highly ionised (Compton-dominated).
          <div class="p-note">ξ = (4π)² J / n<sub>H</sub> (Tarter+ 1969 ionisation parameter).</div></td>
      </tr>
      <tr>
        <td class="p-flag">-frac</td><td class="p-type">float</td><td class="p-default">-1</td>
        <td class="p-desc">frac = F<sub>corona</sub> / F<sub>disk</sub> sets the corona-to-disk illumination ratio.
          frac ≤ 0 (the default, −1) illuminates with the corona only — no thermal disk component.</td>
      </tr>
      <tr>
        <td class="p-flag">-incidence</td><td class="p-type">float</td><td class="p-default">0.7071</td>
        <td class="p-desc">Cosine of the incidence angle cos θ (snapped to nearest Gauss-Legendre
          quadrature node at runtime). Default corresponds to θ = 45°.
          <div class="p-note">Negative μ = downward direction; the code uses |cos θ| and snaps.</div></td>
      </tr>
      <tr>
        <td class="p-flag">-Afe</td><td class="p-type">float</td><td class="p-default">1.0</td>
        <td class="p-desc">Iron abundance relative to solar: A<sub>Fe</sub> = Fe / Fe<sub>⊙</sub>.
          Values &gt; 1 enhance the Fe K complex (6.4–6.97 keV); values &lt; 1 suppress it.</td>
      </tr>
    </table>
    <div class="note-box">
      <strong>Derived quantities:</strong><br>
      ξ = 10<sup>ζ</sup>, &nbsp; n<sub>H</sub> = 10<sup>nh</sup>, &nbsp;
      J = ξ · n<sub>H</sub> / (4π)² &nbsp;
      (the mean intensity that normalises the corona + disk spectra).
    </div>
  </div>

  <!-- ── Accretion Disk ───────────────────── -->
  <div class="section" id="disk">
    <h2><span class="tick" aria-hidden="true"></span>Accretion Disk</h2>
    <table>
      <tr><th>Flag</th><th>Type</th><th>Default</th><th>Description</th></tr>
      <tr>
        <td class="p-flag">-kT_disk</td><td class="p-type">float</td><td class="p-default">0.35</td>
        <td class="p-desc">Effective temperature of the disk blackbody that illuminates the slab from below [keV].</td>
      </tr>
    </table>
  </div>

  <!-- ── Solver / Kernel ──────────────────── -->
  <div class="section" id="solver">
    <h2><span class="tick" aria-hidden="true"></span>Solver / Kernel</h2>
    <table>
      <tr><th>Flag</th><th>Type</th><th>Default</th><th>Description</th></tr>
      <tr>
        <td class="p-flag">-angsca</td><td class="p-type">bool</td><td class="p-default">true</td>
        <td class="p-desc">Compton scattering kernel: <code>true</code>/<code>1</code>/<code>yes</code> →
          angle-dependent kernel (<code>KernelCache</code>); <code>false</code>/<code>0</code>/<code>no</code> →
          angle-averaged kernel (<code>avgKernelCache</code>).</td>
      </tr>
    </table>
  </div>

  <!-- ── Test Mode ────────────────────────── -->
  <div class="section" id="test">
    <h2><span class="tick" aria-hidden="true"></span>Test Mode</h2>
    <p class="lead">
      Test mode bypasses Cloudy and uses a synthetic slab with analytic opacities.
      Useful for validating the RT solver and Compton kernel. The only mode token is
      <code>compps</code> — an isothermal pure-scattering slab illuminated by a bottom
      blackbody seed, benchmarked against Xspec <code>compPS</code> (Poutanen &amp; Svensson 1996).
    </p>
    <table>
      <tr><th>Flag</th><th>Type</th><th>Default</th><th>Description</th></tr>
      <tr>
        <td class="p-flag">-test_rt &lt;mode&gt;</td><td class="p-type">flag + string</td><td class="p-default">off</td>
        <td class="p-desc">Enable test mode (no Cloudy calls). The only mode token is
          <code>compps</code> (compPS benchmark slab).
          Pair with a <code>blackbody</code> corona for the seed.</td>
      </tr>
      <tr>
        <td class="p-flag">-kT_e</td><td class="p-type">float</td><td class="p-default">60</td>
        <td class="p-desc">Slab uniform temperature [keV] in <code>compps</code> test mode (re-uses the
          corona <code>-kT_e</code> flag). Sets the mean fractional energy gain per scattering,
          ≈ 4kT/m<sub>e</sub>c<sup>2</sup> (thermal Doppler width ∝ √(2kT/m<sub>e</sub>c<sup>2</sup>)).</td>
      </tr>
      <tr>
        <td class="p-flag">-tau</td><td class="p-type">float</td><td class="p-default">0.5</td>
        <td class="p-desc">Slab vertical Thomson optical depth in <code>compps</code> test mode
          (matches the <code>compPS</code> <code>tau</code> parameter).</td>
      </tr>
    </table>
  </div>

  <!-- ── Execution flow ───────────────────── -->
  <div class="section">
    <h2><span class="tick" aria-hidden="true"></span>Execution Flow</h2>
    <div class="note-box" style="background:rgba(145,82,55,.06);border-left-color:var(--cyan);">
      <strong style="color:var(--cyan);">Production mode:</strong><br>
      1. Parse CLI → <code>ModelParams</code><br>
      2. Initialise grids (angle, depth, energy)<br>
      3. Compute corona + disk illumination spectra<br>
      4. Precompute Compton kernel + scattering cross-sections (cached on disk)<br>
      5. Outer loop: Cloudy depth sweep → extract j<sub>ν</sub>, κ<sub>abs</sub>, κ<sub>sct</sub> → RT solve → update J, ξ, T → check convergence<br>
      6. Save results to <code>results/&lt;hash&gt;/</code> each iteration
    </div>
    <div class="note-box" style="background:rgba(145,82,55,.06);border-left-color:var(--cyan);">
      <strong style="color:var(--cyan);">Test mode (compps):</strong><br>
      1. Parse CLI → <code>ModelParams</code><br>
      2. Initialise grids with synthetic 1000-bin log-spaced energy grid<br>
      3. Compute bottom blackbody seed illumination<br>
      4. Single RT solve on the isothermal pure-scattering slab (T = kT<sub>e</sub>)<br>
      5. Save results to <code>results/&lt;hash&gt;/</code>
    </div>
  </div>

  <div class="section" id="run-labels">
    <h2><span class="tick" aria-hidden="true"></span>Run labels and summaries</h2>
    <p class="lead">Use the optional run label on the configuration page, or add
      <code>-label "Iron abundance comparison"</code> to a command.</p>
    <p class="lead">Each invocation writes <code>results/&lt;hash&gt;/RUN.txt</code> with its label,
      timestamp, model summary and parameters. The label is also saved in <code>params.json</code>
      and shown in the results selectors. RUN.txt describes the requested run; its presence does
      not mean the calculation has finished.</p>
    <div class="note-box">Labels are metadata and do not change the physics hash.
      Repeating the same physical parameters reuses the same folder and replaces its label and summary.
      The queue and exported commands retain each label.</div>
  </div>

  <!-- ── Example Commands ─────────────────── -->
  <div class="section">
    <h2><span class="tick" aria-hidden="true"></span>Example Commands</h2>
    <div class="formula-box" style="font-size:.78rem;line-height:2;">
      <div><span class="p-flag">./maindaocl</span> -corona cutoffpl -Gamma 2.0 -Ecut 300 -nh 16 -zeta 4 -frac 0.5</div>
      <div><span class="p-flag">./maindaocl</span> -corona nthcomp -Gamma 2.0 -kT_e 100 -kT_bb 0.05 -nh 15 -zeta 3</div>
      <div><span class="p-flag">./maindaocl</span> -corona comptt -kT_e 50 -kT_bb 0.05 -taup 1.0 -Afe 3.0</div>
      <div><span class="p-flag">./maindaocl</span> -test_rt compps -corona blackbody -kT_e 60 -kT_bb 0.1 -tau 0.5</div>
    </div>
  </div>

  <!-- ── References & Citation ────────────── -->
  <div class="section" id="references">
    <h2><span class="tick" aria-hidden="true"></span>References &amp; Citation</h2>
    <p style="font-size:.82rem;color:var(--text);margin-bottom:14px;line-height:1.7;">
      If you use DAO in published work, please cite the DAO paper together with the methods and
      back-ends it is built on:
    </p>
    <table>
      <tr><th>Component</th><th>Reference</th></tr>
      <tr><td class="p-desc">Exact Compton redistribution kernel</td>
          <td class="p-desc">Madej, Różańska, Majczyna &amp; Należyta 2017, MNRAS, 469, 2032 ·
            <a href="https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.2032M" target="_blank" rel="noopener">ADS</a></td></tr>
      <tr><td class="p-desc">Compton cross section &amp; <code>compPS</code> benchmark</td>
          <td class="p-desc">Poutanen &amp; Svensson 1996, ApJ, 470, 249 ·
            <a href="https://ui.adsabs.harvard.edu/abs/1996ApJ...470..249P" target="_blank" rel="noopener">ADS</a></td></tr>
      <tr><td class="p-desc">Atomic-physics back-end (Cloudy C25)</td>
          <td class="p-desc">Gunasekera et al. 2025, arXiv:2508.01102 ·
            <a href="https://ui.adsabs.harvard.edu/abs/2025arXiv250801102G" target="_blank" rel="noopener">ADS</a></td></tr>
      <tr><td class="p-desc">XSPEC models &amp; HEASoft</td>
          <td class="p-desc">Arnaud 1996, ASPC, 101, 17
            (<a href="https://ui.adsabs.harvard.edu/abs/1996ASPC..101...17A" target="_blank" rel="noopener">ADS</a>) ·
            HEASARC 2014, HEASoft, ascl:1408.004
            (<a href="https://ui.adsabs.harvard.edu/abs/2014ascl.soft08004N" target="_blank" rel="noopener">ADS</a>)</td></tr>
      <tr><td class="p-desc">Corona models — nthcomp / comptt</td>
          <td class="p-desc">Zdziarski, Johnson &amp; Magdziarz 1996, MNRAS, 283, 193
            (<a href="https://ui.adsabs.harvard.edu/abs/1996MNRAS.283..193Z" target="_blank" rel="noopener">ADS</a>) ·
            Titarchuk 1994, ApJ, 434, 570
            (<a href="https://ui.adsabs.harvard.edu/abs/1994ApJ...434..570T" target="_blank" rel="noopener">ADS</a>)</td></tr>
      <tr><td class="p-desc">Ionisation parameter ξ</td>
          <td class="p-desc">Tarter, Tucker &amp; Salpeter 1969, ApJ, 156, 943
            (<a href="https://ui.adsabs.harvard.edu/abs/1969ApJ...156..943T" target="_blank" rel="noopener">ADS</a>)</td></tr>
      <tr><td class="p-desc">RT solver — Bézier short characteristics</td>
          <td class="p-desc">Auer 2003
            (<a href="https://ui.adsabs.harvard.edu/abs/2003ASPC..288....3A" target="_blank" rel="noopener">ADS</a>) ·
            de la Cruz Rodríguez &amp; Piskunov 2013
            (<a href="https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D" target="_blank" rel="noopener">ADS</a>) ·
            Suleimanov, Poutanen &amp; Werner 2012
            (<a href="https://ui.adsabs.harvard.edu/abs/2012A%26A...545A.120S" target="_blank" rel="noopener">ADS</a>) ·
            Hubeny &amp; Mihalas 2015 §12.4
            (<a href="https://ui.adsabs.harvard.edu/abs/2014tsa..book.....H" target="_blank" rel="noopener">ADS</a>)</td></tr>
    </table>
  </div>

  <!-- ── License ──────────────────────────── -->
  <div class="section" id="license">
    <h2><span class="tick" aria-hidden="true"></span>License</h2>
    <p style="font-size:.82rem;color:var(--text);line-height:1.7;">
      DAO source code is released under the <strong>MIT License</strong> (see the repository
      <code>LICENSE</code> file). It <em>depends on, but does not bundle,</em> third-party software
      that carries its own separate license, obtained independently:
      <strong>Cloudy</strong> (<a href="https://opensource.org/licenses/Zlib" target="_blank" rel="noopener">zlib license</a>),
      <strong>HEASoft / XSPEC</strong> (NASA HEASARC), and a C++ port of Jerzy Madej's publicly
      available Compton-kernel Fortran code. See <code>README</code> for the full attribution.
    </p>
  </div>

</main>
{{ footer | safe }}
</body>
</html>
"""


@app.route("/docs")
def docs():
    return render_template_string(DOCS_HTML, base_css=BASE_CSS, nav=TOP_NAV("docs"), footer=FOOTER)


# ─── Plots page ─────────────────────────────────────────────────
import datetime
import math
import glob as globmod
import re

# ─── Run identity (server-built labels for /plots & /convergence) ──
# params.json keys per corona model. NOTE the keys differ from the CLI
# flags (E_cut / E_lo_cut / tau_slab); -1 is the "unused" sentinel.
CORONA_ID_PARAMS = {
    "powerlaw":  [("Gamma", "Γ"), ("E_lo_cut", "E_lo")],
    "cutoffpl":  [("Gamma", "Γ"), ("E_cut", "E_cut"), ("E_lo_cut", "E_lo")],
    "nthcomp":   [("Gamma", "Γ"), ("kT_e", "kT_e"), ("kT_bb", "kT_bb")],
    "comptt":    [("kT_e", "kT_e"), ("kT_bb", "kT_bb"), ("taup", "τ_p")],
    "blackbody": [("kT_bb", "kT_bb")],
}
SLAB_ID_PARAMS = [("nh", "log n_H"), ("zeta", "log ξ"),
                  ("incidence", "cos θ"), ("Afe", "A_Fe")]
KEV_KEYS = {"E_cut", "E_lo_cut", "kT_e", "kT_bb", "kT_disk"}
CONV_TOL = 1e-3  # max |ΔT|/T between the last two iterations → "converged"


def _fmt_id(v):
    if isinstance(v, bool):
        return "on" if v else "off"
    if isinstance(v, (int, float)):
        return f"{float(v):.4g}"
    return str(v)


def _time_str(meta, mtime):
    t = meta.get("time")
    if t:
        return str(t)[:16]
    if mtime:
        return datetime.datetime.fromtimestamp(mtime).strftime("%Y-%m-%d %H:%M")
    return ""


def _run_identity(name, meta, mtime, run_id=None):
    """Human-readable identity for one run: dropdown label, plot-title stamp,
    status badges and grouped parameters. Built server-side so both pages
    render the same thing. The plot-title stamp ends with run_id (the
    campaign-qualified relative path) so same-hash runs in different
    campaigns are distinguishable; the dropdown label keeps the basename
    since the <optgroup> already carries the campaign."""
    rid = run_id or name
    tstr = _time_str(meta, mtime)
    if not meta.get("corona"):
        return {"label": f"{name} · {tstr} · (no params.json)".strip(" ·"),
                "title": rid, "time_str": tstr, "groups": [],
                "badges": [{"text": "unknown parameters", "kind": "dim"}]}

    is_test = bool(meta.get("test_rt"))
    frac = meta.get("frac", -1)
    if is_test:
        head = f"{meta.get('test_mode', 'test')} test"
        toks = [f"{lab}={_fmt_id(meta[key])}"
                for key, lab in [("kT_e", "kT_e"), ("tau_slab", "τ")]
                if meta.get(key) not in (None, -1)]
    else:
        head = meta["corona"]
        toks = [f"{lab}={_fmt_id(meta[key])}"
                for key, lab in CORONA_ID_PARAMS.get(meta["corona"], [])
                if meta.get(key) not in (None, -1)]
        toks.append(f"logξ={_fmt_id(meta.get('zeta', '?'))}")
        toks.append(f"logn={_fmt_id(meta.get('nh', '?'))}")
    core = f"{head} · {' '.join(toks)}"

    user_label = str(meta.get("label", "") or "").strip()

    badges = []
    if user_label:
        badges.append({"text": user_label, "kind": "label"})
    if is_test:
        badges.append({"text": f"TEST · {meta.get('test_mode', '?')}", "kind": "test"})
    else:
        badges.append({"text": "production", "kind": "prod"})
    badges.append({"text": "angle-dependent kernel" if meta.get("angsca")
                   else "angle-averaged kernel", "kind": "kernel"})
    if not is_test:
        badges.append({"text": "corona only" if frac is None or frac <= 0
                       else f"F_cor/F_disk = {_fmt_id(frac)}", "kind": "illum"})

    # RT energy range — the quantity campaigns typically vary
    e_lo, e_hi = meta.get("E_rt_lo"), meta.get("E_rt_hi")
    ert_items = ([{"label": "E_RT", "value": f"{_fmt_id(e_lo)} eV – {_fmt_id(e_hi)} eV"}]
                 if e_lo not in (None, -1) and e_hi not in (None, -1) else [])

    groups = []
    if is_test:
        groups.append({"title": "compPS slab", "items": [
            {"label": lab, "value": _fmt_id(meta[key])}
            for key, lab in [("kT_e", "kT_e [keV]"), ("kT_bb", "kT_bb [keV]"),
                             ("tau_slab", "τ_slab")]
            if meta.get(key) not in (None, -1)] + ert_items})
    else:
        groups.append({"title": f"corona · {meta['corona']}", "items": [
            {"label": lab + (" [keV]" if key in KEV_KEYS else ""),
             "value": _fmt_id(meta[key])}
            for key, lab in CORONA_ID_PARAMS.get(meta["corona"], [])
            if meta.get(key) not in (None, -1)]})
        groups.append({"title": "slab", "items": [
            {"label": lab, "value": _fmt_id(meta[key])}
            for key, lab in SLAB_ID_PARAMS if meta.get(key) is not None] + ert_items})
        if frac is not None and frac > 0:
            groups.append({"title": "disk", "items": [
                {"label": "kT_disk [keV]", "value": _fmt_id(meta.get("kT_disk", "?"))}]})

    lbl_prefix = f"“{user_label}” · " if user_label else ""
    return {"label": f"{lbl_prefix}{core} · {tstr} · {name}",
            "title": f"{lbl_prefix}{core} · {rid}",
            "time_str": tstr, "badges": badges, "groups": groups}


def _run_convergence(run_path):
    """Max relative temperature change between the last two profile files."""
    try:
        files = sorted((f for f in os.listdir(run_path)
                        if re.fullmatch(r"profile_iter\d+\.dat", f)),
                       key=lambda f: int(re.search(r"iter(\d+)", f).group(1)))
    except OSError:
        files = []
    if not files:
        return {"status": "none", "dT": None, "n_iter": 0}
    if len(files) == 1:
        return {"status": "single", "dT": None, "n_iter": 1}

    def read_T(fn):
        T = []
        try:
            with open(os.path.join(run_path, fn)) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            T.append((int(parts[0]), float(parts[1]), float(parts[2])))
                        except ValueError:
                            pass
        except OSError:
            pass
        return T

    T_prev, T_last = read_T(files[-2]), read_T(files[-1])
    # A partially written profile must not be reported as converged.
    valid = (T_prev and len(T_prev) == len(T_last)
             and all(a[:2] == b[:2] and math.isfinite(a[2]) and math.isfinite(b[2])
                     and a[2] > 0 and b[2] > 0 for a, b in zip(T_prev, T_last)))
    if not valid:
        return {"status": "none", "dT": None, "n_iter": len(files)}
    dT = max(abs(b[2] - a[2]) / a[2] for a, b in zip(T_prev, T_last))
    return {"status": "converged" if dT <= CONV_TOL else "evolving",
            "dT": dT, "n_iter": len(files)}


def _scan_run_dirs(results_dir):
    """Yield (run_id, group, name, run_path, mtime) for every run directory.

    A run is a directory holding params.json, either directly under results/
    (group = "") or exactly one level down inside a user-named campaign folder
    (group = campaign dir name, run_id = "campaign/name"). Top-level runs come
    first; runs are newest-first within each group; campaign groups are ordered
    by their newest run.
    """
    def mt_of(p):
        try:
            return os.path.getmtime(p)
        except OSError:
            return 0

    root = os.path.realpath(results_dir) + os.sep

    def is_run_dir(path):
        return os.path.isdir(path) and os.path.realpath(path).startswith(root)

    top, campaigns = [], []
    if os.path.isdir(results_dir):
        for name in os.listdir(results_dir):
            path = os.path.join(results_dir, name)
            if not is_run_dir(path):
                continue
            if os.path.exists(os.path.join(path, "params.json")):
                top.append((mt_of(path), name, "", name, path))
            else:
                subs = []
                for sub in os.listdir(path):
                    sp = os.path.join(path, sub)
                    if is_run_dir(sp) and os.path.exists(os.path.join(sp, "params.json")):
                        subs.append((mt_of(sp), f"{name}/{sub}", name, sub, sp))
                if subs:
                    subs.sort(key=lambda e: e[0], reverse=True)
                    campaigns.append(subs)
    top.sort(key=lambda e: e[0], reverse=True)
    campaigns.sort(key=lambda g: g[0][0], reverse=True)
    for group in [top] + campaigns:
        for mt, run_id, grp, name, path in group:
            yield run_id, grp, name, path, mt


def _run_dir(run_id):
    """Resolve a run id ("hash" or "campaign/hash") to its directory under
    results/. 404s on any escape: '..' traversal, absolute paths, symlinks
    pointing outside the tree."""
    results_dir = os.path.realpath(os.path.join(WORK_DIR, "results"))
    p = os.path.realpath(os.path.join(results_dir, run_id))
    if not p.startswith(results_dir + os.sep):
        abort(404)
    return p


@app.route("/api/runs")
def api_runs():
    """Scan results/ (plus one campaign level) and read params.json from each.

    Each run carries "id" (relative path used in data-API URLs), "group"
    (campaign folder, "" for top-level), "hash" (basename, for display) and a
    server-built identity (label, title, badges, param groups) plus a
    convergence status so the client only renders.
    """
    results_dir = os.path.join(WORK_DIR, "results")
    runs = []
    for run_id, group, name, run_path, mt in _scan_run_dirs(results_dir):
        try:
            with open(os.path.join(run_path, "params.json")) as f:
                meta = json.load(f)
        except Exception:
            meta = {}
        # find available iterations
        iters = sorted(set(
            int(mm.group(1))
            for fn in os.listdir(run_path)
            for mm in [re.search(r"iter(\d+)", fn)]
            if mm
        ))
        runs.append({"id": run_id, "group": group, "hash": name, "iters": iters,
                     "meta": meta, "mtime": mt,
                     "convergence": _run_convergence(run_path),
                     **_run_identity(name, meta, mt, run_id=run_id)})
    return jsonify(runs=runs)


@app.route("/api/data/<path:run_id>/<int:iteration>")
def api_data(run_id, iteration):
    """Parse emergent, moments, profile data for a run/iteration."""
    rd = _run_dir(run_id)
    result = {}

    # 1. Emergent intensity
    ef = os.path.join(rd, f"emergent_iter{iteration:03d}.dat")
    if os.path.exists(ef):
        mu_vals = []
        rows = []
        with open(ef) as f:
            for line in f:
                if line.startswith("#"):
                    mm = re.search(r"I\(mu=([-\d.]+)\)", line)
                    if mm:
                        mu_vals.append(float(mm.group(1)))
                else:
                    rows.append([float(x) for x in line.split()])
        if rows:
            E = [r[0] for r in rows]
            I_corona = [r[1] for r in rows]
            I_disk = [r[2] for r in rows]
            angles = {}
            for j, mu in enumerate(mu_vals):
                angles[f"{mu:.4f}"] = [r[3 + j] for r in rows]
            result["emergent"] = {
                "E": E, "I_corona": I_corona, "I_disk": I_disk,
                "mu_vals": mu_vals, "angles": angles
            }

    # 2. Moments (J0 at depth=0 for mean intensity)
    mf = os.path.join(rd, f"moments_iter{iteration:03d}.dat")
    if os.path.exists(mf):
        # E, depth_id, tau_mid, T_K, J0, J2, J3
        j0_surface = {}  # E -> J0 at depth 0
        with open(mf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 5 and int(parts[1]) == 0:
                    j0_surface[float(parts[0])] = float(parts[4])
        if j0_surface:
            E_s = sorted(j0_surface.keys())
            result["mean_intensity"] = {
                "E": E_s, "J0": [j0_surface[e] for e in E_s]
            }

    # 3. Profile
    pf = os.path.join(rd, f"profile_iter{iteration:03d}.dat")
    if os.path.exists(pf):
        rows = []
        with open(pf) as f:
            for line in f:
                if not line.startswith("#"):
                    rows.append([float(x) for x in line.split()])
        if rows:
            result["profile"] = {
                "depth_id": [int(r[0]) for r in rows],
                "tau_mid": [r[1] for r in rows],
                "T_K": [r[2] for r in rows],
                "n_e": [r[3] for r in rows],
                "log_xi": [r[6] for r in rows],
            }

    # 4. Fe K lines — parse fe_lines.dat, find the block for this iteration
    fe_file = os.path.join(rd, "fe_lines.dat")
    if os.path.exists(fe_file):
        fe_lines = []
        in_block = False
        with open(fe_file) as f:
            for line in f:
                if line.startswith("# Fe lines") and re.search(rf"\biter={iteration}\b", line):
                    in_block = True
                    continue
                elif line.startswith("# Fe lines") and in_block:
                    break  # next iteration block
                elif line.startswith("#"):
                    continue
                elif in_block:
                    parts = line.split()
                    if len(parts) >= 5:
                        label = parts[0]
                        e_keV = float(parts[2])
                        relint = float(parts[3])
                        if relint > 0:
                            fe_lines.append({"label": label, "E_keV": e_keV, "relint": relint})
        if fe_lines:
            result["fe_lines"] = fe_lines

    return jsonify(result)


@app.route("/api/profiles/<path:run_id>")
def api_all_profiles(run_id):
    """Return temperature profiles for iterations of a run (every `step`-th)."""
    step = int(request.args.get("step", 2))
    rd = _run_dir(run_id)
    all_iters = []
    raw = {}
    if os.path.isdir(rd):
        for fn in sorted(os.listdir(rd)):
            mm = re.search(r"profile_iter(\d+)\.dat", fn)
            if not mm:
                continue
            it = int(mm.group(1))
            all_iters.append(it)
            rows = []
            with open(os.path.join(rd, fn)) as f:
                for line in f:
                    if not line.startswith("#"):
                        rows.append([float(x) for x in line.split()])
            if rows:
                raw[it] = rows
    # Keep every step-th iteration, always include first and last
    all_iters.sort()
    keep = set(all_iters[::step])
    if all_iters:
        keep.add(all_iters[0])
        keep.add(all_iters[-1])
    profiles = {}
    for it in keep:
        if it in raw:
            profiles[str(it)] = {
                "tau_mid": [r[1] for r in raw[it]],
                "T_K": [r[2] for r in raw[it]],
            }
    return jsonify(profiles=profiles)


@app.route("/api/mean_intensity/<path:run_id>")
def api_all_mean_intensity(run_id):
    """Return surface (depth=0) mean intensity J0(E) for iterations (every `step`-th)."""
    step = int(request.args.get("step", 2))
    rd = _run_dir(run_id)
    # First pass: collect available iteration numbers
    all_iters = []
    if os.path.isdir(rd):
        for fn in os.listdir(rd):
            mm = re.search(r"moments_iter(\d+)\.dat", fn)
            if mm:
                all_iters.append(int(mm.group(1)))
    all_iters.sort()
    keep = set(all_iters[::step])
    if all_iters:
        keep.add(all_iters[0])
        keep.add(all_iters[-1])
    # Second pass: only parse kept iterations
    moments = {}
    for it in sorted(keep):
        fn = os.path.join(rd, f"moments_iter{it:03d}.dat")
        if not os.path.exists(fn):
            continue
        E_list = []
        J0_list = []
        with open(fn) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 5 and int(parts[1]) == 0:
                    E_list.append(float(parts[0]))
                    J0_list.append(float(parts[4]))
        if E_list:
            moments[str(it)] = {"E": E_list, "J0": J0_list}
    return jsonify(moments=moments)


@app.route("/api/line_labels/<path:run_id>")
def api_line_labels(run_id):
    """Return line labels from the first available line_labels file."""
    import re as _re
    rd = _run_dir(run_id)
    if not os.path.isdir(rd):
        return jsonify(lines=[])
    # Find first line_labels file
    ll_files = sorted(f for f in os.listdir(rd) if f.startswith("line_labels_iter"))
    if not ll_files:
        return jsonify(lines=[])
    lines = []
    with open(os.path.join(rd, ll_files[0])) as f:
        for row in f:
            if row.startswith('#'):
                continue
            m = _re.match(r'"(.+?)"\s+(.+)', row)
            if not m:
                continue
            lab = m.group(1).strip()
            fields = m.group(2).split()
            if len(fields) < 6:
                continue
            try:
                eev = float(fields[1])
                emis = float(fields[4])
            except ValueError:
                continue
            ltype = fields[5]
            if emis > 0:
                lines.append({"label": lab, "eV": eev, "emiss": emis, "type": ltype})
    # Aggregate lines at the same energy bin: sum emissivities, keep
    # the label of the strongest contributor.  Then select top lines
    # PER ENERGY DECADE so every part of the spectrum gets representation.
    # Without per-decade selection, UV lines (high emissivity) dominate
    # and X-ray lines (Fe K-alpha, etc.) are excluded.
    import math

    # Group by energy (round to 4 significant figures to merge same-bin lines)
    bins = {}
    for ln in lines:
        if ln["eV"] <= 0:
            continue
        # Round energy to 4 sig figs to group lines in same Cloudy bin
        key = round(ln["eV"], -int(math.floor(math.log10(ln["eV"]))) + 3)
        if key not in bins:
            bins[key] = {"label": ln["label"], "eV": ln["eV"],
                         "emiss": ln["emiss"], "type": ln["type"]}
        else:
            bins[key]["emiss"] += ln["emiss"]
            # Keep label of strongest contributor
            if ln["emiss"] > bins[key]["emiss"] - ln["emiss"]:
                bins[key]["label"] = ln["label"]
                bins[key]["type"] = ln["type"]

    merged = list(bins.values())

    # Select top per decade
    decades = {}
    for ln in merged:
        dec = int(math.floor(math.log10(ln["eV"])))
        decades.setdefault(dec, []).append(ln)

    PER_DECADE = 80
    selected = []
    for dec in sorted(decades.keys()):
        bucket = decades[dec]
        bucket.sort(key=lambda x: -x["emiss"])
        selected.extend(bucket[:PER_DECADE])

    selected.sort(key=lambda x: -x["emiss"])
    return jsonify(lines=selected)


PLOTS_HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<title>DAO — Results Viewer</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
{{ base_css | safe }}
  main.wrapper { max-width:1200px; }
  .controls { display:flex; gap:var(--sp-4); align-items:center; flex-wrap:wrap; margin-bottom:var(--sp-5); }
  .controls label { font-size:var(--fs-sm); color:var(--dim); }
  .controls select { background:var(--bg2); color:var(--white); border:1px solid var(--border);
    border-radius:7px; padding:9px 12px; font-family:inherit; font-size:var(--fs-sm); min-width:0; max-width:100%; width:260px; }
  .controls select:hover { border-color:var(--border-hi); }
  .controls select:focus-visible { border-color:var(--accent); outline:2px solid var(--accent-sft); outline-offset:2px; }
  .params-box { padding:12px 16px; margin-bottom:var(--sp-5); font-size:var(--fs-xs); color:var(--text); line-height:1.8; }
  .params-box .pk { color:var(--cyan); }
  .params-box .pv { color:var(--white); font-weight:600; }
  .plot-card { padding:20px; margin-bottom:16px; }
  .plot-caption { font-size:var(--fs-xs); color:var(--dim); margin-top:8px; line-height:1.5; }
  .plot-area { width:100%; height:480px; }
  @media (max-width:520px){ .plot-area { height:min(60vh,480px); } }
  .status { text-align:center; padding:40px; color:var(--dim); font-size:var(--fs-sm); }
</style>
</head>
<body>
{{ nav | safe }}
<main id="main" class="wrapper">

  <div class="page-head"><h1>Results Viewer</h1><div class="sub">Emergent spectra &amp; temperature profiles per run/iteration</div></div>

  <div class="controls" id="controls">
    <label for="runSelect">Run:</label>
    <select id="runSelect" onchange="onRunChange()"></select>
    <label for="iterSelect">Iteration:</label>
    <select id="iterSelect" onchange="loadData()"></select>
  </div>

  <div class="params-box" id="paramsBox" style="display:none"></div>

  <div id="plotsContainer" style="display:none">
    <div class="plot-card">
      <h2 style="justify-content:flex-start;gap:12px;">
        <span class="tick" aria-hidden="true"></span>1. Emergent intensity at surface (outgoing μ &gt; 0)
        <button id="feToggle" class="pill pill-fe" onclick="toggleFeLines()">Fe lines: ON</button>
      </h2>
      <div class="plot-area" id="plotEmergent"></div>
      <div class="plot-caption">Incident corona shown as 2 I<sub>cor</sub>/μ<sub>inc</sub> (flux→intensity over two hemispheres);
        disk seed shown as I<sub>disk</sub>/2, so seed and emergent beams are comparable.</div>
    </div>
    <div class="plot-card">
      <h2 style="justify-content:flex-start;gap:12px;">
        <span class="tick" aria-hidden="true"></span>2. Angle-averaged outgoing intensity (mean over μ &gt; 0)
        <button id="btnLineLabels" class="pill pill-line off" onclick="toggleLineLabels()">Line IDs: OFF</button>
        <span id="lineLabelStatus" style="font-size:var(--fs-xs);color:var(--dim);"></span>
      </h2>
      <div class="plot-area" id="plotMean"></div>
    </div>
    <div class="plot-card">
      <h2><span class="tick" aria-hidden="true"></span>3. Temperature profile</h2>
      <div class="plot-area" id="plotProfile"></div>
    </div>
  </div>

  <div class="status" id="statusMsg">Loading…</div>
  <div id="emptyState" style="display:none"></div>

</main>
{{ footer | safe }}

<script>
const THEME = getComputedStyle(document.documentElement);
const themeColor = name => THEME.getPropertyValue(name).trim();
const PLOT_BG = themeColor('--card');
const GRID_COLOR = themeColor('--border');
const FONT_COLOR = themeColor('--text');
const PAPER_BG = themeColor('--card');
const DIM_HEX = themeColor('--dim');
const FAINT_HEX = themeColor('--faint');
const FE_HEX = themeColor('--red');
const LINE_HEX = themeColor('--accent');
const COLORS = ['#005A9C','#B34700','#20704A','#693C91','#A3264D','#755400','#006B73','#3F569B','#8D4038','#386447'];

let allRuns = [];
let runStamp = '';   // run identity stamped into every plot title

function esc(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

const BADGE_CLS = { label:'b-label', test:'b-test', prod:'b-prod', kernel:'b-kernel', illum:'b-illum', dim:'b-dim' };

function convBadge(run) {
  const c = run.convergence || {};
  const n = c.n_iter || (run.iters ? run.iters.length : 0);
  if (c.status === 'converged')
    return `<span class="rc-badge b-ok">temperature converged · max ΔT/T = ${esc(c.dT.toExponential(1))} · ${n} iter</span>`;
  if (c.status === 'evolving')
    return `<span class="rc-badge b-warn">temperature evolving · max ΔT/T = ${esc(c.dT.toExponential(1))} · ${n} iter</span>`;
  if (c.status === 'single') return '<span class="rc-badge b-dim">single iteration</span>';
  return '<span class="rc-badge b-dim">no profile data</span>';
}

function runSetKey(runs) {
  return JSON.stringify(runs);
}

function buildRunSelect(runs, keepId) {
  // top-level runs first (plain options), then one <optgroup> per campaign;
  // server order already groups them contiguously
  const sel = document.getElementById('runSelect');
  sel.innerHTML = '';
  let parent = sel;
  let curGroup = '';
  runs.forEach(run => {
    if ((run.group || '') !== curGroup) {
      curGroup = run.group || '';
      if (curGroup) {
        parent = document.createElement('optgroup');
        parent.label = curGroup;
        sel.appendChild(parent);
      } else parent = sel;
    }
    const o = document.createElement('option');
    o.value = run.id;
    o.textContent = run.label || run.id;
    o.title = run.label || run.id;
    parent.appendChild(o);
  });
  if (keepId && runs.some(r => r.id === keepId)) sel.value = keepId;
  else if (runs.length) sel.value = runs[0].id;
}

async function init() {
  const r = await fetch('/api/runs');
  const d = await r.json();
  allRuns = d.runs;  // server order: top-level newest-first, then campaigns
  if (!allRuns.length) { showEmpty(); return; }
  buildRunSelect(allRuns, null);
  onRunChange();
}

function showEmpty() {
  document.getElementById('runSelect').innerHTML = '';
  document.getElementById('plotsContainer').style.display = 'none';
  document.getElementById('paramsBox').style.display = 'none';
  document.getElementById('controls').style.display = 'none';
  document.getElementById('statusMsg').style.display = 'none';
  const e = document.getElementById('emptyState');
  e.style.display = 'block';
  e.innerHTML = `<div class="empty-state plot-card">
    <div class="ico" aria-hidden="true">◎</div>
    <h3>No results yet</h3>
    <p>Configure a model and run ./maindaocl — finished runs appear here automatically.</p>
    <a class="btn btn-primary" href="/">Open Configurator</a>
  </div>`;
}

function onRunChange() {
  const rid = document.getElementById('runSelect').value;
  const run = allRuns.find(r => r.id === rid);
  const iterSel = document.getElementById('iterSelect');
  iterSel.innerHTML = '';
  if (!run) return;
  run.iters.forEach(it => {
    const o = document.createElement('option');
    o.value = it; o.textContent = `Iteration ${it}`;
    iterSel.appendChild(o);
  });
  if (run.iters.length > 0) iterSel.value = run.iters[run.iters.length - 1];
  showRunCard(run);
  loadData();
}

function showRunCard(run) {
  const box = document.getElementById('paramsBox');
  if (!run) { box.style.display = 'none'; return; }
  const badges = (run.badges || []).map(b =>
    `<span class="rc-badge ${BADGE_CLS[b.kind] || 'b-dim'}">${esc(b.text)}</span>`).join('');
  const groups = (run.groups || []).map(g =>
    `<div class="rc-group"><div class="rg-title">${esc(g.title)}</div><div class="rg-items">` +
    g.items.map(it => `<span class="pk">${esc(it.label)}</span> = <span class="pv">${esc(it.value)}</span>`).join(' &nbsp;·&nbsp; ') +
    '</div></div>').join('');
  const model = run.title ? run.title.split(' · ')[0] : run.hash;
  const metaHash = run.meta && run.meta.hash;
  const mono = metaHash && metaHash !== run.hash ? `${run.id} (${metaHash})` : run.id;
  box.innerHTML = `
    <div class="rc-top">
      <span class="rc-model">${esc(model)}</span>
      ${badges}${convBadge(run)}
      <span class="rc-meta">${esc(run.time_str || '')} <span class="rc-hash">${esc(mono)}</span></span>
    </div>
    <div class="rc-groups">${groups}</div>`;
  box.classList.add('run-card');
  box.style.display = 'block';
}

async function loadData() {
  lineLabelsOn = false;
  document.getElementById('btnLineLabels').textContent = 'Line IDs: OFF';
  document.getElementById('btnLineLabels').className = 'pill pill-line off';
  document.getElementById('lineLabelStatus').textContent = '';
  const rid = document.getElementById('runSelect').value;
  const iter = document.getElementById('iterSelect').value;
  if (!rid || !iter) return;
  document.getElementById('statusMsg').style.display = 'block';
  document.getElementById('statusMsg').textContent = 'Loading data...';
  document.getElementById('plotsContainer').style.display = 'none';

  const r = await fetch(`/api/data/${encodeURIComponent(rid)}/${iter}`);
  const d = await r.json();

  if (!d.emergent && !d.profile) {
    document.getElementById('statusMsg').textContent = 'No data files found for this run/iteration.';
    return;
  }

  document.getElementById('statusMsg').style.display = 'none';
  document.getElementById('plotsContainer').style.display = 'block';

  const run = allRuns.find(r => r.id === rid);
  runStamp = ((run && run.title) || rid) + ' · iter ' + iter;
  const mu_inc = Math.abs((run && run.meta && run.meta.incidence) || 0.7071);

  plotEmergent(d.emergent, mu_inc, d.fe_lines || []);
  plotMeanOutgoing(d.emergent);
  plotProfile(d.profile);
}

let showFeLines = true;
let lastEmData = null;
let lastMuInc = 0.7071;
let lastFeLines = [];

function toggleFeLines() {
  showFeLines = !showFeLines;
  const btn = document.getElementById('feToggle');
  btn.textContent = showFeLines ? 'Fe lines: ON' : 'Fe lines: OFF';
  btn.className = 'pill ' + (showFeLines ? 'pill-fe' : 'off');
  plotEmergent(lastEmData, lastMuInc, lastFeLines);
}

function autoLogRange(traces) {
  let peak = -Infinity;
  traces.forEach(t => { t.y.forEach(v => { if (v > 0 && v > peak) peak = v; }); });
  if (!isFinite(peak) || peak <= 0) return undefined;
  return [Math.log10(peak * 1e-6), Math.log10(peak * 5)];
}

function plotEmergent(em, mu_inc, fe_lines) {
  lastEmData = em; lastMuInc = mu_inc; lastFeLines = fe_lines;
  const div = document.getElementById('plotEmergent');
  if (!em) { div.innerHTML = '<div class="status">No emergent data.</div>'; return; }

  const E_eV = em.E;
  const E_keV = E_eV.map(e => e / 1e3);
  const traces = [];

  traces.push({
    x: E_keV, y: E_eV.map((e, i) => e * 2.0 * em.I_corona[i] / mu_inc),
    name: 'E × 2J_cor/μ_inc', mode: 'lines',
    line: { color: DIM_HEX, width: 1.5, dash: 'dash' }
  });
  traces.push({
    x: E_keV, y: E_eV.map((e, i) => e * em.I_disk[i] / 2.0),
    name: 'E × I_disk/2', mode: 'lines',
    line: { color: FAINT_HEX, width: 1.5, dash: 'dot' }
  });

  let ci = 0;
  em.mu_vals.forEach(mu => {
    if (mu > 0) {
      const vals = em.angles[mu.toFixed(4)];
      traces.push({
        x: E_keV, y: E_eV.map((e, i) => e * vals[i]),
        name: `μ = ${mu.toFixed(4)}`, mode: 'lines',
        line: { color: COLORS[ci % COLORS.length], width: 1.8, dash: ['solid','dash','dot','dashdot'][Math.floor(ci / COLORS.length) % 4] }
      });
      ci++;
    }
  });

  const yr = autoLogRange(traces);
  const shapes = [];
  const annotations = [];
  if (showFeLines && fe_lines && fe_lines.length > 0 && yr) {
    const top6 = [...fe_lines].sort((a, b) => b.relint - a.relint).slice(0, 6);
    top6.forEach((fl, idx) => {
      shapes.push({
        type: 'line', xref: 'x', yref: 'paper',
        x0: fl.E_keV, x1: fl.E_keV, y0: 0, y1: 1,
        line: { color: 'rgba(173,77,66,0.5)', width: 1, dash: 'dot' }
      });
      annotations.push({
        x: Math.log10(fl.E_keV), xref: 'x', yref: 'paper',
        y: 1.0 - idx * 0.06, text: fl.label,
        showarrow: false, font: { color: FE_HEX, size: 11 },
        xanchor: 'left', xshift: 4
      });
    });
  }

  Plotly.react(div, traces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    title: { text: esc(runStamp), font: { size: 11, color: '#6f675e' }, x: 0.01, xanchor: 'left' },
    margin: { l:70, r:30, t:52, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12} },
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'E I_E  [erg cm⁻² s⁻¹ sr⁻¹]', range:yr, gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
    shapes: shapes, annotations: annotations,
  }, {responsive:true});
}

function plotMeanOutgoing(em) {
  const div = document.getElementById('plotMean');
  if (!em) { div.innerHTML = '<div class="status">No emergent data.</div>'; return; }

  const E_eV = em.E;
  const E_keV = E_eV.map(e => e / 1e3);
  const outgoing_keys = em.mu_vals.filter(mu => mu > 0).map(mu => mu.toFixed(4));
  const n_out = outgoing_keys.length;
  const mean_I = E_eV.map((_, ie) => {
    let sum = 0; outgoing_keys.forEach(k => { sum += em.angles[k][ie]; });
    return n_out > 0 ? sum / n_out : 0;
  });

  const traces = [{
    x: E_keV, y: E_eV.map((e, i) => e * mean_I[i]),
    name: 'E × Mean outgoing I', mode: 'lines',
    line: { color: '#915237', width: 2 }
  }];

  const yr = autoLogRange(traces);
  Plotly.react(div, traces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    title: { text: esc(runStamp), font: { size: 11, color: '#6f675e' }, x: 0.01, xanchor: 'left' },
    margin: { l:70, r:30, t:48, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12} },
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'E I_E  [erg cm⁻² s⁻¹ sr⁻¹]', range:yr, gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

function plotProfile(prof) {
  const div = document.getElementById('plotProfile');
  if (!prof) { div.innerHTML = '<div class="status">No profile data.</div>'; return; }
  const traces = [{
    x: prof.tau_mid, y: prof.T_K,
    name: 'T', mode: 'lines+markers',
    line: { color: '#753b28', width: 2 },
    marker: { size: 4, color: '#915237' }
  }];
  Plotly.react(div, traces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    title: { text: esc(runStamp), font: { size: 11, color: '#6f675e' }, x: 0.01, xanchor: 'left' },
    margin: { l:70, r:30, t:48, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12} },
    xaxis: { type:'log', title:'τ (Thomson)', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'T [K]', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

// ── Line label overlay ─────────────────────────────────────────
let lineLabelsOn = false;
let lineLabelsData = null;
const ROMAN = {1:'I',2:'II',3:'III',4:'IV',5:'V',6:'VI',7:'VII',8:'VIII',
               9:'IX',10:'X',11:'XI',12:'XII',13:'XIII',14:'XIV',15:'XV',
               16:'XVI',17:'XVII',18:'XVIII',19:'XIX',20:'XX',21:'XXI',
               22:'XXII',23:'XXIII',24:'XXIV',25:'XXV',26:'XXVI',27:'XXVII'};
function fmtSpecies(lab) {
  const m = lab.match(/^([A-Z][a-z]?)\s*(\d+)$/);
  if (m) { const n = parseInt(m[2]); if (ROMAN[n]) return m[1]+' '+ROMAN[n]; }
  return lab;
}
async function toggleLineLabels() {
  const btn = document.getElementById('btnLineLabels');
  const status = document.getElementById('lineLabelStatus');
  const div = document.getElementById('plotMean');
  const rid = document.getElementById('runSelect').value;
  if (lineLabelsOn) {
    lineLabelsOn = false;
    btn.textContent = 'Line IDs: OFF';
    btn.className = 'pill pill-line off';
    status.textContent = '';
    Plotly.relayout(div, { annotations: [], shapes: [] });
    return;
  }
  if (!lineLabelsData || lineLabelsData.hash !== rid) {
    status.textContent = 'Loading...';
    const r = await fetch(`/api/line_labels/${encodeURIComponent(rid)}`);
    const d = await r.json();
    lineLabelsData = { hash: rid, lines: d.lines || [] };
  }
  if (lineLabelsData.lines.length === 0) { status.textContent = 'No line labels found.'; return; }
  const annotations = [];
  const shapes = [];
  const MAX_LABELS = 60;
  const MIN_SEP = 0.015;
  const placedLogE = [];
  const topLines = lineLabelsData.lines.slice(0, 100);
  for (const ln of topLines) {
    if (annotations.length >= MAX_LABELS) break;
    const keV = ln.eV / 1e3;
    const logE = Math.log10(keV);
    if (placedLogE.some(p => Math.abs(logE - p) < MIN_SEP)) continue;
    placedLogE.push(logE);
    const isFluor = ln.type === 'F';
    const name = fmtSpecies(ln.label) + (isFluor ? ' (fl)' : '');
    const col = isFluor ? 'rgba(173,77,66,0.85)' : 'rgba(133,83,53,0.85)';
    const lcol = isFluor ? 'rgba(173,77,66,0.2)' : 'rgba(133,83,53,0.15)';
    shapes.push({ type:'line', xref:'x', yref:'paper', x0:keV, x1:keV, y0:0, y1:1, line:{color:lcol, width:0.8} });
    annotations.push({ x:Math.log10(keV), y:1, xref:'x', yref:'paper', text:name, showarrow:false,
      font:{size:9, color:col, family:'Inter, sans-serif'}, textangle:-90, xanchor:'left', yanchor:'top', yshift:-4 });
  }
  lineLabelsOn = true;
  btn.textContent = 'Line IDs: ON';
  btn.className = 'pill pill-line';
  status.textContent = `${annotations.length} lines`;
  Plotly.relayout(div, { annotations, shapes });
}

// ── Live refresh — runs keep growing while the page is open ────
let refreshBusy = false;

async function refreshRuns() {
  if (refreshBusy || document.hidden) return;
  refreshBusy = true;
  try {
    const r = await fetch('/api/runs');
    const d = await r.json();
    const runs = d.runs;
    if (!runs.length) { allRuns = []; showEmpty(); return; }
    document.getElementById('emptyState').style.display = 'none';
    document.getElementById('controls').style.display = '';
    const sel = document.getElementById('runSelect');
    const prevId = sel.value;
    const setChanged = runSetKey(runs) !== runSetKey(allRuns);
    allRuns = runs;
    if (setChanged) {
      buildRunSelect(runs, prevId);  // preserves selection when possible
      if (sel.value !== prevId) { onRunChange(); return; }
    }
    const run = allRuns.find(x => x.id === sel.value);
    if (!run) return;
    showRunCard(run);  // convergence badge / iter count may have advanced
    // append new iterations; follow the newest if the user was on the last one
    const iterSel = document.getElementById('iterSelect');
    const have = new Set(Array.from(iterSel.options).map(o => Number(o.value)));
    const wasLast = iterSel.selectedIndex < 0 ||
                    iterSel.selectedIndex === iterSel.options.length - 1;
    let added = false;
    run.iters.forEach(it => {
      if (!have.has(it)) {
        const o = document.createElement('option');
        o.value = it; o.textContent = `Iteration ${it}`;
        iterSel.appendChild(o);
        added = true;
      }
    });
    if (added && wasLast && run.iters.length) {
      iterSel.value = run.iters[run.iters.length - 1];
      loadData();
    }
  } catch (e) { /* transient poll failure — retry next tick */ }
  finally { refreshBusy = false; }
}

setInterval(refreshRuns, 5000);
document.addEventListener('visibilitychange',
  () => { if (!document.hidden) refreshRuns(); });

init();
</script>
</body>
</html>
"""

@app.route("/plots")
def plots():
    return render_template_string(PLOTS_HTML, base_css=BASE_CSS, nav=TOP_NAV("plots"), footer=FOOTER)


# ─── Compare page ────────────────────────────────────────────────
# Overlay emergent spectra (or temperature profiles) from any number of
# (run, iteration) selections. Reuses /api/runs + /api/data; all overlay
# maths run client-side.
COMPARE_HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<title>DAO — Compare Spectra</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
{{ base_css | safe }}
  main.wrapper { max-width:1200px; }
  .controls { display:flex; gap:var(--sp-3); align-items:center; flex-wrap:wrap; margin-bottom:var(--sp-4); }
  .controls label { font-size:var(--fs-sm); color:var(--dim); }
  .controls select { background:var(--bg2); color:var(--white); border:1px solid var(--border);
    border-radius:var(--r-sm); padding:7px 12px; font-family:inherit; font-size:var(--fs-sm); outline:none; min-width:200px; }
  .controls select:focus { border-color:var(--accent); }
  .controls select#iterSelect { min-width:120px; }
  .mode-row { display:flex; gap:var(--sp-2); align-items:center; flex-wrap:wrap; margin-bottom:var(--sp-5); }
  .seg { display:inline-flex; border:1px solid var(--border); border-radius:var(--r-sm); overflow:hidden; }
  .seg button { background:var(--bg2); color:var(--dim); border:none; font-family:inherit;
    font-size:var(--fs-sm); padding:7px 14px; cursor:pointer; transition:all .15s; }
  .seg button.active { background:var(--accent-sft); color:var(--accent); font-weight:600; }
  .seg button + button { border-left:1px solid var(--border); }

  /* ── Compare list ─────────────────────────────────── */
  .cmp-list { display:flex; flex-direction:column; gap:6px; margin-bottom:var(--sp-5); }
  .cmp-row { display:flex; align-items:center; gap:var(--sp-3); background:var(--card);
    border:1px solid var(--border); border-radius:var(--r-sm); padding:8px 12px; }
  .cmp-row .swatch { width:14px; height:14px; border-radius:3px; flex:none; box-shadow:0 0 6px rgba(0,0,0,.4); }
  .cmp-row .cmp-id { font-family:var(--mono); font-size:var(--fs-sm); color:var(--white); font-weight:600; }
  .cmp-row .cmp-sub { font-size:var(--fs-xs); color:var(--dim); }
  .cmp-row .ref-tag { font-size:var(--fs-xs); color:var(--cyan); border:1px solid rgba(133,83,53,.35);
    border-radius:999px; padding:1px 8px; }
  .cmp-row .cmp-spacer { margin-left:auto; }
  .cmp-row .cmp-mkref { background:var(--bg2); border:1px solid var(--border); color:var(--dim);
    font-size:var(--fs-xs); padding:4px 10px; border-radius:var(--r-sm); cursor:pointer; font-family:inherit; }
  .cmp-row .cmp-mkref:hover { color:var(--cyan); border-color:rgba(133,83,53,.4); }
  .cmp-row .cmp-del { background:var(--bg2); border:1px solid var(--border); color:var(--text);
    font-size:var(--fs-xs); padding:4px 10px; border-radius:var(--r-sm); cursor:pointer; font-family:inherit; }
  .cmp-row .cmp-del:hover { color:var(--red); border-color:rgba(173,77,66,.4); }
  @media (max-width:600px) {
    .cmp-row { flex-wrap:wrap; }
    .cmp-row > div { flex:1; min-width:190px; overflow-wrap:anywhere; }
    .cmp-row .cmp-spacer { display:none; }
    .controls select { max-width:100%; min-width:0; }
  }
  .cmp-empty { color:var(--dim); font-size:var(--fs-sm); padding:14px 4px; }
  .plot-card { padding:16px; margin-bottom:var(--sp-5); }
  .plot-area { width:100%; height:520px; }
  .plot-area.ratio { height:360px; }
  @media (max-width:520px){ .plot-area { height:min(60vh,480px); } }
  .plot-caption { font-size:var(--fs-xs); color:var(--dim); margin-top:8px; line-height:1.5; }
  .status { text-align:center; padding:40px; color:var(--dim); font-size:var(--fs-sm); }
</style>
</head>
<body>
{{ nav | safe }}
<main id="main" class="wrapper">

  <div class="page-head"><h1>Compare Spectra</h1>
    <div class="sub">Overlay emergent spectra or temperature profiles from any set of runs &amp; iterations</div></div>

  <div class="card full" id="picker">
    <h2><span class="tick" aria-hidden="true"></span>Add a spectrum</h2>
    <div class="controls">
      <label for="runSelect">Run:</label>
      <select id="runSelect" onchange="onRunChange()"></select>
      <label for="iterSelect">Iteration:</label>
      <select id="iterSelect"></select>
      <button class="btn btn-primary" onclick="addSpectrum()">+ Add to comparison</button>
    </div>
    <div class="cmp-list" id="cmpList"></div>
    <div class="controls" style="margin-bottom:0">
      <button class="btn btn-secondary" onclick="clearAll()">Clear all</button>
    </div>
  </div>

  <div class="mode-row">
    <span style="font-size:var(--fs-sm);color:var(--dim);">Quantity:</span>
    <div class="seg" id="quantSeg">
      <button data-q="emergent" class="active" onclick="setQuantity('emergent')">Emergent spectrum</button>
      <button data-q="temperature" onclick="setQuantity('temperature')">Temperature profile</button>
    </div>
    <span class="pill pill-line off" id="normBtn" onclick="toggleNorm()"
      style="margin-left:8px;cursor:pointer;">Normalize by integrated intensity: OFF</span>
    <span class="pill pill-line off" id="ratioBtn" onclick="toggleRatio()"
      style="cursor:pointer;">Ratio to reference: OFF</span>
  </div>

  <div id="plotsContainer">
    <div class="plot-card">
      <h2 style="justify-content:flex-start;gap:12px;">
        <span class="tick" aria-hidden="true"></span><span id="mainTitle">Emergent spectra — angle-averaged outgoing</span>
      </h2>
      <div class="plot-area" id="plotMain"></div>
      <div class="plot-caption" id="mainCaption">Mean over outgoing angles (μ &gt; 0) of each selected run/iteration, plotted as E·I<sub>E</sub>.</div>
    </div>
    <div class="plot-card" id="ratioCard" style="display:none">
      <h2 style="justify-content:flex-start;gap:12px;">
        <span class="tick" aria-hidden="true"></span>Ratio to reference
      </h2>
      <div class="plot-area ratio" id="plotRatio"></div>
      <div class="plot-caption" id="ratioCaption"></div>
    </div>
  </div>

  <div class="status" id="statusMsg">Loading…</div>
  <div id="emptyState" style="display:none"></div>

</main>
{{ footer | safe }}

<script>
const THEME = getComputedStyle(document.documentElement);
const themeColor = name => THEME.getPropertyValue(name).trim();
const PLOT_BG = themeColor('--card');
const GRID_COLOR = themeColor('--border');
const FONT_COLOR = themeColor('--text');
const PAPER_BG = themeColor('--card');
const COLORS = ['#005A9C','#B34700','#20704A','#693C91','#A3264D','#755400','#006B73','#3F569B','#8D4038','#386447'];
const CKEY = 'dao.compare.v1';

let allRuns = [];
let runById = {};
let cmp = [];                 // [{id, iter, color}]
const dataCache = {};         // `${id}/${iter}` -> /api/data json (or {__empty:true})
let quantity = 'emergent';
let ratioOn = false;
let normOn = false;
let lastRunSetKey = '';

function esc(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}
function keyOf(e) { return e.id + '/' + e.iter; }
function pickColor() {
  const used = new Set(cmp.map(e => e.color));
  for (const c of COLORS) if (!used.has(c)) return c;
  return COLORS[cmp.length % COLORS.length];
}
function paramPart(run) {
  // run.title = "<core> · <id>"; drop the trailing id token for a compact sublabel
  if (!run || !run.title) return '';
  const parts = run.title.split(' · ');
  return parts.length > 1 ? parts.slice(0, -1).join(' · ') : run.title;
}
function legendLabel(e) { return e.id + ' · it' + e.iter; }

// ── run set key (for live refresh of the add-picker) ─────────────
function runSetKey(runs) { return JSON.stringify(runs); }

function buildRunSelect(runs, keepId) {
  const sel = document.getElementById('runSelect');
  sel.innerHTML = '';
  let parent = sel, curGroup = '';
  runs.forEach(run => {
    if ((run.group || '') !== curGroup) {
      curGroup = run.group || '';
      if (curGroup) { parent = document.createElement('optgroup'); parent.label = curGroup; sel.appendChild(parent); }
      else parent = sel;
    }
    const o = document.createElement('option');
    o.value = run.id; o.textContent = run.label || run.id; o.title = run.label || run.id;
    parent.appendChild(o);
  });
  if (keepId && runs.some(r => r.id === keepId)) sel.value = keepId;
  else if (runs.length) sel.value = runs[0].id;
}

function onRunChange() {
  const rid = document.getElementById('runSelect').value;
  const run = runById[rid];
  const iterSel = document.getElementById('iterSelect');
  iterSel.innerHTML = '';
  if (!run) return;
  run.iters.forEach(it => {
    const o = document.createElement('option');
    o.value = it; o.textContent = 'Iteration ' + it;
    iterSel.appendChild(o);
  });
  if (run.iters.length) iterSel.value = run.iters[run.iters.length - 1];
}

async function init() {
  const r = await fetch('/api/runs');
  const d = await r.json();
  allRuns = d.runs || [];
  runById = {}; allRuns.forEach(run => runById[run.id] = run);
  lastRunSetKey = runSetKey(allRuns);
  document.getElementById('statusMsg').style.display = 'none';
  if (!allRuns.length) { showEmpty(); return; }
  buildRunSelect(allRuns, null);
  onRunChange();
  rehydrate();
  renderList();
  redraw();
}

function showEmpty() {
  document.getElementById('picker').style.display = 'none';
  document.querySelector('.mode-row').style.display = 'none';
  document.getElementById('plotsContainer').style.display = 'none';
  const e = document.getElementById('emptyState');
  e.style.display = 'block';
  e.innerHTML = `<div class="empty-state plot-card">
    <div class="ico" aria-hidden="true">◎</div>
    <h3>No results yet</h3>
    <p>Configure a model and run ./maindaocl — finished runs appear here, ready to compare.</p>
    <a class="btn btn-primary" href="/">Open Configurator</a>
  </div>`;
}

let refreshBusy = false;
async function refreshRuns() {
  if (refreshBusy || document.hidden) return;
  refreshBusy = true;
  try {
    const r = await fetch('/api/runs');
    const d = await r.json();
    const runs = d.runs || [];
    const key = runSetKey(runs);
    if (key !== lastRunSetKey) {
      lastRunSetKey = key;
      allRuns = runs;
      runById = {}; allRuns.forEach(run => runById[run.id] = run);
      const keep = document.getElementById('runSelect').value;
      const keepIter = document.getElementById('iterSelect').value;
      const wasEmpty = document.getElementById('picker').style.display === 'none';
      buildRunSelect(allRuns, keep);
      onRunChange();
      if (runById[keep]?.iters.includes(Number(keepIter))) {
        document.getElementById('iterSelect').value = keepIter;
      }
      if (!runs.length) { showEmpty(); return; }
      document.getElementById('picker').style.display = '';
      document.querySelector('.mode-row').style.display = '';
      document.getElementById('plotsContainer').style.display = '';
      document.getElementById('emptyState').style.display = 'none';
      if (wasEmpty) rehydrate();
      renderList();
    }
    // A selected iteration may still be writing; retry it on the next refresh.
    Object.keys(dataCache).forEach(k => delete dataCache[k]);
    if (cmp.length) await redraw();
  } catch (e) { /* transient; retry next tick */ }
  finally { refreshBusy = false; }
}

// ── compare set ──────────────────────────────────────────────────
function addSpectrum() {
  const id = document.getElementById('runSelect').value;
  const iterStr = document.getElementById('iterSelect').value;
  if (!id || iterStr === '') return;
  const iter = parseInt(iterStr, 10);
  const e = { id, iter, color: null };
  if (cmp.some(x => keyOf(x) === keyOf(e))) return;   // no duplicates
  e.color = pickColor();
  cmp.push(e);
  persist(); renderList(); redraw();
}
function removeSpectrum(idx) { cmp.splice(idx, 1); persist(); renderList(); redraw(); }
function makeRef(idx) {
  if (idx <= 0 || idx >= cmp.length) return;
  const [e] = cmp.splice(idx, 1);
  cmp.unshift(e);
  persist(); renderList(); redraw();
}
function clearAll() { cmp = []; persist(); renderList(); redraw(); }

function renderList() {
  const box = document.getElementById('cmpList');
  if (!cmp.length) {
    box.innerHTML = '<div class="cmp-empty">No spectra selected — pick a run and iteration above, then “Add to comparison”. Add two or more to overlay.</div>';
    return;
  }
  box.innerHTML = cmp.map((e, i) => {
    const run = runById[e.id];
    const sub = paramPart(run);
    const ref = i === 0 ? '<span class="ref-tag">reference</span>' : '';
    const mkref = i === 0 ? '' : `<button class="cmp-mkref" onclick="makeRef(${i})">Set as reference</button>`;
    return `<div class="cmp-row">
      <span class="swatch" style="background:${e.color}"></span>
      <div>
        <div><span class="cmp-id">${esc(e.id)} · iter ${esc(e.iter)}</span> ${ref}</div>
        <div class="cmp-sub">${esc(sub)}</div>
      </div>
      <span class="cmp-spacer"></span>
      ${mkref}
      <button class="cmp-del" onclick="removeSpectrum(${i})">✕ Remove</button>
    </div>`;
  }).join('');
}

// ── persistence ──────────────────────────────────────────────────
function persist() {
  try { localStorage.setItem(CKEY, JSON.stringify(cmp.map(e => ({id:e.id, iter:e.iter, color:e.color})))); } catch (x) {}
}
function rehydrate() {
  try {
    const raw = localStorage.getItem(CKEY);
    if (!raw) return;
    const saved = JSON.parse(raw);
    if (!Array.isArray(saved)) return;
    cmp = [];
    saved.forEach(s => {
      const run = runById[s.id];
      if (run && run.iters.includes(s.iter) && !cmp.some(x => keyOf(x) === s.id + '/' + s.iter)) {
        cmp.push({ id: s.id, iter: s.iter, color: COLORS.includes(s.color) ? s.color : pickColor() });
      }
    });
  } catch (x) { cmp = []; }
}

// ── quantity / ratio controls ────────────────────────────────────
function setQuantity(q) {
  quantity = q;
  document.querySelectorAll('#quantSeg button').forEach(b =>
    b.classList.toggle('active', b.dataset.q === q));
  const isEm = q === 'emergent';
  document.getElementById('mainTitle').textContent = isEm
    ? 'Emergent spectra — angle-averaged outgoing' : 'Temperature profiles';
  document.getElementById('mainCaption').innerHTML = isEm
    ? 'Mean over outgoing angles (μ &gt; 0) of each selected run/iteration, plotted as E·I<sub>E</sub>.'
    : 'Gas temperature vs Thomson depth τ for each selected run/iteration.';
  document.getElementById('ratioBtn').style.display = isEm ? '' : 'none';
  document.getElementById('normBtn').style.display = isEm ? '' : 'none';
  redraw();
}
function toggleRatio() {
  ratioOn = !ratioOn;
  const b = document.getElementById('ratioBtn');
  b.textContent = 'Ratio to reference: ' + (ratioOn ? 'ON' : 'OFF');
  b.className = 'pill ' + (ratioOn ? 'pill-line' : 'pill-line off');
  redraw();
}
function toggleNorm() {
  normOn = !normOn;
  const b = document.getElementById('normBtn');
  b.textContent = 'Normalize by integrated intensity: ' + (normOn ? 'ON' : 'OFF');
  b.className = 'pill ' + (normOn ? 'pill-line' : 'pill-line off');
  redraw();
}

// ── data fetch (cached) ──────────────────────────────────────────
async function fetchEntry(e) {
  const k = keyOf(e);
  if (dataCache[k]) return dataCache[k];
  try {
    const r = await fetch(`/api/data/${encodeURIComponent(e.id)}/${e.iter}`);
    if (!r.ok) throw new Error('Data not ready');
    const d = await r.json();
    dataCache[k] = d;
    return d;
  } catch (x) { return { __empty: true }; }
}

function meanOutgoing(em) {
  const keys = em.mu_vals.filter(mu => mu > 0).map(mu => mu.toFixed(4));
  const n = keys.length;
  return em.E.map((_, ie) => { let s = 0; keys.forEach(k => s += em.angles[k][ie]); return n ? s / n : 0; });
}
// total flux F = ∫ I_E dE (trapezoid, E ascending in eV)
function totalFlux(E, I) {
  let s = 0;
  for (let k = 1; k < E.length; k++) s += 0.5 * (I[k] + I[k - 1]) * (E[k] - E[k - 1]);
  return s;
}
// mean outgoing intensity, divided by its total flux when normalization is on
function seriesI(em) {
  const meanI = meanOutgoing(em);
  if (!normOn) return meanI;
  const F = totalFlux(em.E, meanI);
  return F > 0 ? meanI.map(v => v / F) : meanI;
}
function autoLogRange(traces) {
  let peak = -Infinity;
  traces.forEach(t => t.y.forEach(v => { if (v > 0 && v > peak) peak = v; }));
  if (!isFinite(peak) || peak <= 0) return undefined;
  return [Math.log10(peak * 1e-6), Math.log10(peak * 5)];
}
// linear interpolation in log-log space (falls back to linear near non-positive values)
function interpAt(xs, ys, x) {
  const N = xs.length;
  if (N === 0) return NaN;
  if (x < xs[0] || x > xs[N - 1]) return NaN;
  if (x === xs[0]) return ys[0];
  if (x === xs[N - 1]) return ys[N - 1];
  let lo = 0, hi = N - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  const x0 = xs[lo], x1 = xs[hi], y0 = ys[lo], y1 = ys[hi];
  if (x0 > 0 && x1 > 0 && y0 > 0 && y1 > 0 && x > 0) {
    const t = (Math.log10(x) - Math.log10(x0)) / (Math.log10(x1) - Math.log10(x0));
    return Math.pow(10, Math.log10(y0) + t * (Math.log10(y1) - Math.log10(y0)));
  }
  const t = (x - x0) / (x1 - x0);
  return y0 + t * (y1 - y0);
}

let drawVersion = 0;
async function redraw() {
  const version = ++drawVersion;
  const status = document.getElementById('statusMsg');
  if (!cmp.length) {
    status.style.display = 'none';
    document.getElementById('plotMain').innerHTML =
      '<div class="status">Add two or more spectra above to overlay them here.</div>';
    document.getElementById('ratioCard').style.display = 'none';
    return;
  }
  status.style.display = 'block'; status.textContent = 'Loading data…';
  const datas = await Promise.all(cmp.map(fetchEntry));
  if (version !== drawVersion) return;
  status.style.display = 'none';
  if (quantity === 'emergent') plotEmergentOverlay(datas);
  else plotTempOverlay(datas);
}

function plotEmergentOverlay(datas) {
  const traces = [];
  const missing = [];
  cmp.forEach((e, i) => {
    const d = datas[i];
    if (!d || !d.emergent) { missing.push(legendLabel(e)); return; }
    const em = d.emergent;
    const meanI = seriesI(em);
    traces.push({
      x: em.E.map(v => v / 1e3), y: em.E.map((v, k) => v * meanI[k]),
      name: legendLabel(e), mode: 'lines', line: { color: e.color, width: 2 },
    });
  });
  const yr = autoLogRange(traces);
  Plotly.react('plotMain', traces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:20, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12} },
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log',
      title: normOn ? 'E·I_E / F_tot  (normalized)' : 'E·I_E  [erg cm⁻² s⁻¹ sr⁻¹]',
      range:yr, gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
  document.getElementById('mainCaption').innerHTML =
    'Mean over outgoing angles (μ &gt; 0) of each selected run/iteration, plotted as E·I<sub>E</sub>.' +
    (normOn ? ' Each spectrum is divided by its own integrated intensity F<sub>tot</sub> = ∫I<sub>E</sub>dE '
            + '(over that run\'s energy grid), so curves compare shape independent of luminosity.' : '');
  noteMissing(missing);

  // ── ratio panel ──
  const card = document.getElementById('ratioCard');
  if (!ratioOn || cmp.length < 2) { card.style.display = 'none'; return; }
  const refD = datas[0];
  if (!refD || !refD.emergent) { card.style.display = 'none'; return; }
  const refEm = refD.emergent;
  const refMean = seriesI(refEm);
  const refE = refEm.E;                // eV, ascending
  const rtraces = [];
  cmp.forEach((e, i) => {
    if (i === 0) return;
    const d = datas[i];
    if (!d || !d.emergent) return;
    const em = d.emergent;
    const meanI = seriesI(em);
    const ratio = refE.map((Eev, k) => {
      const num = interpAt(em.E, meanI, Eev);
      const den = refMean[k];
      return (den > 0 && num > 0) ? num / den : null;
    });
    rtraces.push({
      x: refE.map(v => v / 1e3), y: ratio, name: legendLabel(e) + ' / ref',
      mode: 'lines', line: { color: e.color, width: 2 }, connectgaps: false,
    });
  });
  card.style.display = 'block';
  document.getElementById('ratioCaption').innerHTML =
    (normOn ? 'Intensity-normalized ' : 'Angle-averaged ') +
    'emergent intensity of each spectrum divided by the reference (<b>' +
    esc(legendLabel(cmp[0])) + '</b>), interpolated onto the reference energy grid within the overlapping energy range.';
  Plotly.react('plotRatio', rtraces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:20, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12} },
    shapes: [{ type:'line', xref:'paper', x0:0, x1:1, yref:'y', y0:1, y1:1,
               line:{ color:'#7c746b', width:1, dash:'dot' } }],
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'I / I_ref', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

function plotTempOverlay(datas) {
  document.getElementById('ratioCard').style.display = 'none';
  const traces = [];
  const missing = [];
  cmp.forEach((e, i) => {
    const d = datas[i];
    if (!d || !d.profile) { missing.push(legendLabel(e)); return; }
    traces.push({
      x: d.profile.tau_mid, y: d.profile.T_K, name: legendLabel(e),
      mode: 'lines+markers', line: { color: e.color, width: 2 }, marker: { size: 3, color: e.color },
    });
  });
  Plotly.react('plotMain', traces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:20, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12} },
    xaxis: { type:'log', title:'τ (Thomson)', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'T [K]', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
  noteMissing(missing);
}

function noteMissing(missing) {
  const cap = document.getElementById('mainCaption');
  cap.querySelectorAll('.miss-note').forEach(n => n.remove());
  if (missing.length) {
    const s = document.createElement('div');
    s.className = 'miss-note';
    s.style.cssText = 'color:#ad4d42;margin-top:4px;';
    s.textContent = 'No ' + (quantity === 'emergent' ? 'emergent' : 'profile') +
      ' data for: ' + missing.join(', ');
    cap.appendChild(s);
  }
}

setInterval(refreshRuns, 5000);
document.addEventListener('visibilitychange', () => { if (!document.hidden) refreshRuns(); });
init();
</script>
</body>
</html>
"""

@app.route("/compare")
def compare():
    return render_template_string(COMPARE_HTML, base_css=BASE_CSS, nav=TOP_NAV("compare"), footer=FOOTER)


# ─── Convergence page ────────────────────────────────────────────
CONVERGENCE_HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<title>DAO — Temperature Convergence</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
{{ base_css | safe }}
  main.wrapper { max-width:1200px; }
  .controls { display:flex; gap:var(--sp-4); align-items:center; flex-wrap:wrap; margin-bottom:var(--sp-5); }
  .controls label { font-size:var(--fs-sm); color:var(--dim); }
  .controls select { background:var(--bg2); color:var(--white); border:1px solid var(--border);
    border-radius:7px; padding:9px 12px; font-family:inherit; font-size:var(--fs-sm); min-width:0; max-width:100%; width:260px; }
  .controls select:hover { border-color:var(--border-hi); }
  .controls select:focus-visible { border-color:var(--accent); outline:2px solid var(--accent-sft); outline-offset:2px; }
  .params-box { padding:12px 16px; margin-bottom:var(--sp-5); font-size:var(--fs-xs); color:var(--text); line-height:1.8; }
  .params-box .pk { color:var(--cyan); }
  .params-box .pv { color:var(--white); font-weight:600; }
  .plot-card { padding:20px; margin-bottom:16px; }
  .plot-area { width:100%; height:540px; }
  @media (max-width:520px){ .plot-area { height:min(60vh,480px); } }
  .status { text-align:center; padding:40px; color:var(--dim); font-size:var(--fs-sm); }
</style>
</head>
<body>
{{ nav | safe }}
<main id="main" class="wrapper">

  <div class="page-head"><h1>Temperature Convergence</h1><div class="sub">Iteration-by-iteration temperature &amp; surface mean intensity</div></div>

  <div class="controls" id="controls">
    <label for="runSelect">Run:</label>
    <select id="runSelect" onchange="onRunChange()"></select>
  </div>

  <div class="params-box" id="paramsBox" style="display:none"></div>

  <div id="plotsContainer" style="display:none">
    <div class="plot-card">
      <h2><span class="tick" aria-hidden="true"></span>1. Temperature profile — all iterations</h2>
      <div class="plot-area" id="plotTempIter"></div>
    </div>
    <div class="plot-card">
      <h2><span class="tick" aria-hidden="true"></span>2. Surface mean intensity J<sub>0</sub>(E) — all iterations</h2>
      <div class="plot-area" id="plotJ0Iter"></div>
    </div>
  </div>

  <div class="status" id="statusMsg">Loading…</div>
  <div id="emptyState" style="display:none"></div>

</main>
{{ footer | safe }}

<script>
const THEME = getComputedStyle(document.documentElement);
const themeColor = name => THEME.getPropertyValue(name).trim();
const PLOT_BG = themeColor('--card');
const GRID_COLOR = themeColor('--border');
const FONT_COLOR = themeColor('--text');
const PAPER_BG = themeColor('--card');
const COLORS_POOL = ['#005A9C','#B34700','#20704A','#693C91','#A3264D','#755400','#006B73','#3F569B','#8D4038','#386447'];

let allRuns = [];
let runStamp = '';   // run identity stamped into every plot title

function esc(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

const BADGE_CLS = { label:'b-label', test:'b-test', prod:'b-prod', kernel:'b-kernel', illum:'b-illum', dim:'b-dim' };

function convBadge(run) {
  const c = run.convergence || {};
  const n = c.n_iter || (run.iters ? run.iters.length : 0);
  if (c.status === 'converged')
    return `<span class="rc-badge b-ok">temperature converged · max ΔT/T = ${esc(c.dT.toExponential(1))} · ${n} iter</span>`;
  if (c.status === 'evolving')
    return `<span class="rc-badge b-warn">temperature evolving · max ΔT/T = ${esc(c.dT.toExponential(1))} · ${n} iter</span>`;
  if (c.status === 'single') return '<span class="rc-badge b-dim">single iteration</span>';
  return '<span class="rc-badge b-dim">no profile data</span>';
}

function runSetKey(runs) {
  return JSON.stringify(runs);
}

function buildRunSelect(runs, keepId) {
  // top-level runs first (plain options), then one <optgroup> per campaign;
  // server order already groups them contiguously
  const sel = document.getElementById('runSelect');
  sel.innerHTML = '';
  let parent = sel;
  let curGroup = '';
  runs.forEach(run => {
    if ((run.group || '') !== curGroup) {
      curGroup = run.group || '';
      if (curGroup) {
        parent = document.createElement('optgroup');
        parent.label = curGroup;
        sel.appendChild(parent);
      } else parent = sel;
    }
    const o = document.createElement('option');
    o.value = run.id;
    o.textContent = run.label || run.id;
    o.title = run.label || run.id;
    parent.appendChild(o);
  });
  if (keepId && runs.some(r => r.id === keepId)) sel.value = keepId;
  else if (runs.length) sel.value = runs[0].id;
}

async function init() {
  const r = await fetch('/api/runs');
  const d = await r.json();
  allRuns = d.runs;  // server order: top-level newest-first, then campaigns
  if (!allRuns.length) { showEmpty(); return; }
  buildRunSelect(allRuns, null);
  onRunChange();
}

function showEmpty() {
  document.getElementById('runSelect').innerHTML = '';
  document.getElementById('plotsContainer').style.display = 'none';
  document.getElementById('paramsBox').style.display = 'none';
  document.getElementById('controls').style.display = 'none';
  document.getElementById('statusMsg').style.display = 'none';
  const e = document.getElementById('emptyState');
  e.style.display = 'block';
  e.innerHTML = `<div class="empty-state plot-card">
    <div class="ico" aria-hidden="true">◎</div>
    <h3>No results yet</h3>
    <p>Configure a model and run ./maindaocl — finished runs appear here automatically.</p>
    <a class="btn btn-primary" href="/">Open Configurator</a>
  </div>`;
}

function onRunChange() {
  const rid = document.getElementById('runSelect').value;
  const run = allRuns.find(r => r.id === rid);
  if (!run) return;
  runStamp = run.title || rid;
  shownIterCount = run.iters.length;
  showRunCard(run);
  loadProfiles(rid);
}

function showRunCard(run) {
  const box = document.getElementById('paramsBox');
  if (!run) { box.style.display = 'none'; return; }
  const badges = (run.badges || []).map(b =>
    `<span class="rc-badge ${BADGE_CLS[b.kind] || 'b-dim'}">${esc(b.text)}</span>`).join('');
  const groups = (run.groups || []).map(g =>
    `<div class="rc-group"><div class="rg-title">${esc(g.title)}</div><div class="rg-items">` +
    g.items.map(it => `<span class="pk">${esc(it.label)}</span> = <span class="pv">${esc(it.value)}</span>`).join(' &nbsp;·&nbsp; ') +
    '</div></div>').join('');
  const model = run.title ? run.title.split(' · ')[0] : run.hash;
  const metaHash = run.meta && run.meta.hash;
  const mono = metaHash && metaHash !== run.hash ? `${run.id} (${metaHash})` : run.id;
  box.innerHTML = `
    <div class="rc-top">
      <span class="rc-model">${esc(model)}</span>
      ${badges}${convBadge(run)}
      <span class="rc-meta">${esc(run.time_str || '')} <span class="rc-hash">${esc(mono)}</span></span>
    </div>
    <div class="rc-groups">${groups}</div>`;
  box.classList.add('run-card');
  box.style.display = 'block';
}

async function loadProfiles(rid) {
  document.getElementById('statusMsg').style.display = 'block';
  document.getElementById('statusMsg').textContent = 'Loading profiles...';
  document.getElementById('plotsContainer').style.display = 'none';

  const [rProf, rJ0] = await Promise.all([
    fetch(`/api/profiles/${encodeURIComponent(rid)}?step=2`),
    fetch(`/api/mean_intensity/${encodeURIComponent(rid)}?step=2`),
  ]);
  const dProf = await rProf.json();
  const dJ0 = await rJ0.json();

  const hasProf = dProf.profiles && Object.keys(dProf.profiles).length > 0;
  const hasJ0 = dJ0.moments && Object.keys(dJ0.moments).length > 0;

  if (!hasProf && !hasJ0) {
    document.getElementById('statusMsg').textContent = 'No profile data found for this run.';
    return;
  }

  document.getElementById('statusMsg').style.display = 'none';
  document.getElementById('plotsContainer').style.display = 'block';
  if (hasProf) plotTempIter(dProf.profiles);
  else document.getElementById('plotTempIter').innerHTML = '<div class="status">No profile data.</div>';
  plotJ0Iter(dJ0.moments);
}

function plotTempIter(profiles) {
  const div = document.getElementById('plotTempIter');
  const iters = Object.keys(profiles).map(Number).sort((a, b) => a - b);
  const n = iters.length;
  const traces = [];
  iters.forEach((it, idx) => {
    const color = iterColor(idx, n);
    const p = profiles[String(it)];
    traces.push({
      x: p.tau_mid, y: p.T_K, name: `iter ${it}`, mode: 'lines',
      line: { color: color, width: idx === n - 1 ? 2.5 : 1.2 },
      opacity: 0.4 + 0.6 * (n > 1 ? idx / (n - 1) : 1),
    });
  });
  Plotly.react(div, traces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    title: { text: esc(runStamp), font: { size: 11, color: '#6f675e' }, x: 0.01, xanchor: 'left' },
    margin: { l:70, r:30, t:48, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12}, orientation:'h', y:-0.15 },
    xaxis: { type:'log', title:'τ (Thomson)', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'T [K]', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

function iterColor(idx) {
  return COLORS_POOL[idx % COLORS_POOL.length];
}

function plotJ0Iter(moments) {
  const div = document.getElementById('plotJ0Iter');
  if (!moments || Object.keys(moments).length === 0) {
    div.innerHTML = '<div class="status">No moments data across iterations.</div>';
    return;
  }
  const iters = Object.keys(moments).map(Number).sort((a, b) => a - b);
  const n = iters.length;
  const traces = [];
  iters.forEach((it, idx) => {
    const color = iterColor(idx, n);
    const m = moments[String(it)];
    const E_keV = m.E.map(e => e / 1e3);
    const EJ0 = m.E.map((e, i) => e * m.J0[i]);
    traces.push({
      x: E_keV, y: EJ0, name: `iter ${it}`, mode: 'lines',
      line: { color: color, width: idx === n - 1 ? 2.5 : 1.2 },
      opacity: 0.4 + 0.6 * (n > 1 ? idx / (n - 1) : 1),
    });
  });
  let peak = -Infinity;
  traces.forEach(t => t.y.forEach(v => { if (v > 0 && v > peak) peak = v; }));
  const yr = isFinite(peak) && peak > 0 ? [Math.log10(peak * 1e-6), Math.log10(peak * 5)] : undefined;
  Plotly.react(div, traces, {
    paper_bgcolor: PAPER_BG, plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    title: { text: esc(runStamp), font: { size: 11, color: '#6f675e' }, x: 0.01, xanchor: 'left' },
    margin: { l:70, r:30, t:48, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:12}, orientation:'h', y:-0.15 },
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'E × J₀ [erg cm⁻² s⁻¹ sr⁻¹]', range:yr, gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

// ── Live refresh — runs keep growing while the page is open ────
let refreshBusy = false;
let shownIterCount = 0;   // iterations currently plotted for the selected run

async function refreshRuns() {
  if (refreshBusy || document.hidden) return;
  refreshBusy = true;
  try {
    const r = await fetch('/api/runs');
    const d = await r.json();
    const runs = d.runs;
    if (!runs.length) { allRuns = []; showEmpty(); return; }
    document.getElementById('emptyState').style.display = 'none';
    document.getElementById('controls').style.display = '';
    const sel = document.getElementById('runSelect');
    const prevId = sel.value;
    const setChanged = runSetKey(runs) !== runSetKey(allRuns);
    allRuns = runs;
    if (setChanged) {
      buildRunSelect(runs, prevId);  // preserves selection when possible
      if (sel.value !== prevId) { onRunChange(); return; }
    }
    const run = allRuns.find(x => x.id === sel.value);
    if (!run) return;
    showRunCard(run);  // convergence badge / iter count may have advanced
    if (run.iters.length !== shownIterCount) {
      shownIterCount = run.iters.length;
      loadProfiles(run.id);  // new iterations → re-fetch and re-plot
    }
  } catch (e) { /* transient poll failure — retry next tick */ }
  finally { refreshBusy = false; }
}

setInterval(refreshRuns, 5000);
document.addEventListener('visibilitychange',
  () => { if (!document.hidden) refreshRuns(); });

init();
</script>
</body>
</html>
"""

@app.route("/convergence")
def convergence():
    return render_template_string(CONVERGENCE_HTML, base_css=BASE_CSS, nav=TOP_NAV("conv"), footer=FOOTER)


# ─── Main ───────────────────────────────────────────────────────
def main():
    port = 5200
    url = f"http://127.0.0.1:{port}"
    print(f"\n  DAO — X-ray Reflection Spectroscopy Model")
    print(f"  Author: Yimin Huang\n")
    print(f"  Opening browser → {url}\n")
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    app.run(host="127.0.0.1", port=port, debug=False, threaded=True)


if __name__ == "__main__":
    main()
