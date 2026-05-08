#!/usr/bin/env python3
"""
DAOv2.0 — X-ray Reflection Spectroscopy Model
Web-based parameter configurator & launcher
Author: Yimin Huang
"""

import json
import os
import threading
import webbrowser
from flask import Flask, render_template_string, request, jsonify

app = Flask(__name__, static_folder='image', static_url_path='/image')

# ─── Working directory (where maindaocl lives) ─────────────────
WORK_DIR = os.path.dirname(os.path.abspath(__file__))

# ─── Model definitions (mirrors params.cpp) ────────────────────
CORONA_MODELS = {
    "powerlaw": {
        "label": "Power Law",
        "desc": "F(E) ∝ E<sup>−Γ</sup> · exp(−E<sub>lo</sub>/E)",
        "params": [
            {"name": "Gamma", "label": "Γ", "desc": "Photon index", "default": 2.0, "min": 1.0, "max": 4.0, "step": 0.01},
            {"name": "E_low_cut", "label": "E<sub>lo</sub>", "desc": "Low-energy exp cutoff [keV]", "default": 0.1, "min": 0.01, "max": 10, "step": 0.01},
        ],
    },
    "cutoffpl": {
        "label": "Cutoff Power Law",
        "desc": "F(E) ∝ E<sup>−Γ</sup> · exp(−E/E<sub>cut</sub>) · exp(−E<sub>lo</sub>/E)",
        "params": [
            {"name": "Gamma", "label": "Γ", "desc": "Photon index", "default": 2.0, "min": 1.0, "max": 4.0, "step": 0.01},
            {"name": "Ecut",  "label": "E<sub>cut</sub>", "desc": "High-energy cutoff [keV]", "default": 300.0, "min": 10, "max": 2000, "step": 1},
            {"name": "E_low_cut", "label": "E<sub>lo</sub>", "desc": "Low-energy exp cutoff [keV]", "default": 0.1, "min": 0.01, "max": 10, "step": 0.01},
        ],
    },
    "nthcomp": {
        "label": "NTHComp",
        "desc": "Thermal Comptonisation (Zdziarski+ 1996)",
        "params": [
            {"name": "Gamma", "label": "Γ", "desc": "Photon index", "default": 2.0, "min": 1.0, "max": 4.0, "step": 0.01},
            {"name": "kT_e",  "label": "kT<sub>e</sub>", "desc": "Electron temperature [keV]", "default": 60.0, "min": 1, "max": 500, "step": 0.5},
            {"name": "kT_bb", "label": "kT<sub>bb</sub>", "desc": "Seed photon temperature [keV]", "default": 0.1, "min": 0.001, "max": 5, "step": 0.001},
        ],
    },
    "comptt": {
        "label": "CompTT",
        "desc": "Comptonisation model (Titarchuk 1994)",
        "params": [
            {"name": "kT_e",  "label": "kT<sub>e</sub>", "desc": "Plasma temperature [keV]", "default": 50.0, "min": 1, "max": 500, "step": 0.5},
            {"name": "kT_bb", "label": "kT<sub>bb</sub>", "desc": "Soft photon temperature [keV]", "default": 0.05, "min": 0.001, "max": 5, "step": 0.001},
            {"name": "taup",  "label": "τ<sub>p</sub>", "desc": "Plasma optical depth", "default": 1.0, "min": 0.01, "max": 20, "step": 0.01},
        ],
    },
    "blackbody": {
        "label": "Blackbody",
        "desc": "B(E) = (2E³/h²c²) / [exp(E/kT) − 1]",
        "params": [
            {"name": "kT_bb", "label": "kT<sub>bb</sub>", "desc": "Temperature [keV]", "default": 0.05, "min": 0.001, "max": 50, "step": 0.001},
        ],
    },
}

SLAB_PARAMS = [
    {"name": "nh",        "label": "log n<sub>H</sub>",  "desc": "Hydrogen density [cm⁻³]",      "default": 15.0,   "min": 10,  "max": 22,   "step": 0.1},
    {"name": "zeta",      "label": "log ξ",              "desc": "Ionisation parameter",          "default": 3.0,    "min": 0,   "max": 6,    "step": 0.1},
    {"name": "frac",      "label": "f<sub>cor</sub>",    "desc": "F<sub>corona</sub> / F<sub>disk</sub>", "default": 100, "min": 0, "max": 100, "step": 0.01},
    {"name": "incidence", "label": "cos θ",              "desc": "Incidence angle cosine",        "default": 0.7071, "min": 0.01,"max": 1,    "step": 0.0001},
    {"name": "Afe",       "label": "A<sub>Fe</sub>",     "desc": "Iron abundance [solar]",        "default": 1.0,    "min": 0.1, "max": 10,   "step": 0.1},
]

DISK_PARAMS = [
    {"name": "kT_disk", "label": "kT<sub>disk</sub>", "desc": "Disk blackbody temperature [eV]", "default": 0.35, "min": 0.01, "max": 5, "step": 0.01},
]


# ─── HTML template ──────────────────────────────────────────────
HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<title>DAOv2.0 — X-ray Reflection Model</title>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  :root {
    --bg:#080a0f; --bg2:#0f1218; --card:#131720; --card-hi:#181d28;
    --border:#1e2433; --border-hi:#2a3244;
    --accent:#e8993a; --accent2:#d4622a; --cyan:#5cc8f0;
    --text:#b8bfcc; --dim:#5f6878; --white:#e8ecf2;
    --red:#e85050; --green:#3cc868; --radius:12px;
    --shadow:0 2px 12px rgba(0,0,0,.4);
  }
  body {
    font-family:'Inter','Helvetica Neue','Segoe UI',system-ui,sans-serif;
    background:var(--bg); color:var(--text); min-height:100vh;
    overflow-x:hidden; font-size:14px; line-height:1.6;
    -webkit-font-smoothing:antialiased;
  }
  body::before {
    content:''; position:fixed; inset:0; z-index:0;
    background:url('/image/bg_blackhole.png') center top / cover no-repeat fixed;
    opacity:.28; pointer-events:none;
  }
  body::after {
    content:''; position:fixed; inset:0; z-index:0;
    background:linear-gradient(180deg,transparent 0%,var(--bg) 55%);
    pointer-events:none;
  }
  .wrapper { position:relative; z-index:1; max-width:1100px; margin:0 auto; padding:24px 28px; }

  .hero { text-align:center; padding:72px 0 24px; }
  .hero h1 {
    font-size:2.4rem; font-weight:700; letter-spacing:.18em;
    background:linear-gradient(90deg,#e8993a,#d4622a,#e8993a);
    background-size:200% 100%;
    -webkit-background-clip:text; -webkit-text-fill-color:transparent;
    animation:shimmer 5s ease-in-out infinite;
  }
  @keyframes shimmer { 0%,100%{background-position:0% 50%} 50%{background-position:100% 50%} }
  .hero .sub { font-size:.9rem; color:var(--dim); margin-top:8px; letter-spacing:.04em; }
  .hero .author { font-size:.82rem; color:#484f5c; margin-top:6px; }
  .hero .author a { color:var(--accent); text-decoration:none; border-bottom:1px dotted rgba(232,153,58,.4); transition:border-color .2s; }
  .hero .author a:hover { border-color:var(--accent); }
  .hero .author a.link-cyan { color:var(--cyan); border-color:rgba(92,200,240,.4); }
  .hero .author a.link-cyan:hover { border-color:var(--cyan); }

  .divider { height:1px; margin:24px auto; max-width:520px; background:linear-gradient(90deg,transparent 0%,#2a3244 20%,var(--accent) 50%,#2a3244 80%,transparent 100%); border:none; }

  .grid { display:grid; grid-template-columns:1fr 1fr; gap:18px; margin-bottom:18px; }
  @media (max-width:760px) { .grid { grid-template-columns:1fr; } }

  .card {
    background:var(--card); border:1px solid var(--border); border-radius:var(--radius);
    padding:24px; box-shadow:var(--shadow); transition:border-color .25s,box-shadow .25s;
  }
  .card:hover { border-color:var(--border-hi); box-shadow:0 4px 20px rgba(0,0,0,.5); }
  .card.full { grid-column:1/-1; }
  .card h2 {
    font-size:.82rem; font-weight:600; text-transform:uppercase; letter-spacing:.13em;
    color:var(--accent); margin-bottom:18px; display:flex; align-items:center; gap:10px;
  }
  .card h2 .icon { font-size:1.05rem; opacity:.7; }

  .model-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(175px,1fr)); gap:10px; margin-bottom:20px; }
  .model-btn {
    background:var(--bg2); border:1.5px solid var(--border); border-radius:10px;
    padding:14px 12px; cursor:pointer; text-align:center; transition:all .25s;
  }
  .model-btn:hover { border-color:var(--border-hi); background:var(--card-hi); }
  .model-btn.active { border-color:var(--accent); background:rgba(232,153,58,.06); box-shadow:0 0 20px rgba(232,153,58,.08),inset 0 0 0 1px rgba(232,153,58,.1); }
  .model-btn .name { font-size:.88rem; font-weight:600; color:var(--white); }
  .model-btn .mdesc { font-size:.72rem; color:var(--dim); margin-top:4px; line-height:1.4; }

  .field { display:grid; grid-template-columns:120px 1fr 85px; align-items:center; gap:12px; margin-bottom:12px; }
  .field label { font-size:.82rem; font-weight:500; color:var(--cyan); text-align:right; white-space:nowrap; }
  .field input[type=range] { -webkit-appearance:none; width:100%; height:5px; background:linear-gradient(90deg,#1a2030,#2a3248); border-radius:3px; outline:none; }
  .field input[type=range]::-webkit-slider-thumb { -webkit-appearance:none; width:16px; height:16px; background:var(--accent); border-radius:50%; cursor:pointer; box-shadow:0 0 8px rgba(232,153,58,.35); transition:box-shadow .2s; }
  .field input[type=range]::-webkit-slider-thumb:hover { box-shadow:0 0 14px rgba(232,153,58,.5); }
  .field .val { font-size:.82rem; color:var(--white); text-align:right; font-weight:500; background:var(--bg2); border:1px solid var(--border); border-radius:8px; padding:6px 10px; width:85px; font-family:inherit; transition:border-color .2s; }
  .field .val:focus { border-color:var(--accent); outline:none; box-shadow:0 0 0 2px rgba(232,153,58,.15); }
  .field-desc { font-size:.72rem; color:var(--dim); grid-column:2/-1; margin-top:-6px; margin-bottom:4px; }

  .toggle-row { display:flex; align-items:center; gap:14px; margin-bottom:12px; }
  .toggle-row label { font-size:.82rem; font-weight:500; color:var(--cyan); min-width:120px; text-align:right; }
  .toggle { position:relative; width:44px; height:24px; cursor:pointer; }
  .toggle input { opacity:0; width:0; height:0; }
  .toggle .slider { position:absolute; inset:0; background:#1a2030; border-radius:12px; border:1px solid var(--border); transition:all .25s; }
  .toggle .slider::before { content:''; position:absolute; width:18px; height:18px; left:2px; top:2px; background:#4a5060; border-radius:50%; transition:all .25s; }
  .toggle input:checked+.slider { background:rgba(232,153,58,.2); border-color:rgba(232,153,58,.4); }
  .toggle input:checked+.slider::before { transform:translateX(20px); background:var(--accent); }

  .derived { display:flex; gap:28px; flex-wrap:wrap; margin-top:12px; padding-top:12px; }
  .derived .item { text-align:center; }
  .derived .item .dlabel { font-size:.72rem; color:var(--dim); font-weight:500; text-transform:uppercase; letter-spacing:.06em; }
  .derived .item .dval { font-size:1.15rem; font-weight:700; color:var(--white); margin-top:2px; }

  .cmd-box { background:var(--bg); border:1px solid var(--border); border-radius:10px; padding:16px 20px; font-size:.85rem; color:var(--green); word-break:break-all; line-height:1.7; font-family:'SF Mono','Fira Code','Cascadia Code',monospace; }
  .cmd-box .prompt { color:var(--dim); user-select:none; }
  .cmd-hint { font-size:.75rem; color:var(--dim); margin-top:10px; line-height:1.5; }
  .cmd-hint code { color:var(--cyan); font-family:'SF Mono','Fira Code',monospace; font-size:.78rem; }

  .actions { display:flex; gap:10px; margin-top:16px; flex-wrap:wrap; }
  .btn { padding:10px 26px; border-radius:9px; font-size:.84rem; font-weight:600; cursor:pointer; border:none; font-family:inherit; transition:all .25s; letter-spacing:.03em; }
  .btn-primary { background:linear-gradient(135deg,var(--accent),var(--accent2)); color:#000; box-shadow:0 2px 10px rgba(232,153,58,.2); }
  .btn-primary:hover { box-shadow:0 4px 24px rgba(232,153,58,.35); transform:translateY(-1px); }
  .btn-secondary { background:var(--bg2); color:var(--text); border:1px solid var(--border); }
  .btn-secondary:hover { border-color:var(--border-hi); background:var(--card-hi); }

  .toast { position:fixed; bottom:32px; left:50%; transform:translateX(-50%); background:var(--card); color:var(--green); border:1px solid rgba(60,200,104,.2); padding:12px 28px; border-radius:10px; font-size:.84rem; font-weight:500; opacity:0; transition:opacity .3s; pointer-events:none; z-index:99; box-shadow:0 4px 20px rgba(0,0,0,.5); }
  .toast.show { opacity:1; }

  .queue-table { width:100%; border-collapse:collapse; font-size:.82rem; margin-top:12px; }
  .queue-table th { text-align:left; color:var(--accent); font-weight:600; font-size:.75rem; text-transform:uppercase; letter-spacing:.08em; padding:10px 12px; border-bottom:1px solid var(--border); }
  .queue-table td { padding:10px 12px; border-bottom:1px solid #151a24; color:var(--text); vertical-align:top; }
  .queue-table tr:hover td { background:rgba(232,153,58,.03); }
  .queue-table .q-idx { color:var(--dim); font-weight:600; width:40px; text-align:center; }
  .queue-table .q-model { color:var(--cyan); font-weight:600; width:100px; }
  .queue-table .q-params { color:var(--text); font-size:.78rem; word-break:break-all; font-family:'SF Mono','Fira Code',monospace; }
  .queue-table .q-actions { width:80px; text-align:center; white-space:nowrap; }
  .queue-table .q-btn { background:none; border:none; cursor:pointer; font-size:.92rem; padding:4px 8px; border-radius:6px; transition:background .15s; }
  .queue-table .q-btn:hover { background:rgba(255,255,255,.06); }
  .queue-table .q-btn-del { color:var(--red); }
  .queue-table .q-btn-load { color:var(--cyan); }
  .queue-empty { text-align:center; color:var(--dim); font-size:.85rem; padding:28px 0; }
  .batch-actions { display:flex; gap:10px; margin-top:16px; flex-wrap:wrap; align-items:center; }
  .batch-count { font-size:.8rem; color:var(--dim); margin-left:auto; font-weight:500; }

  .footer { text-align:center; padding:32px 0 16px; font-size:.75rem; color:#363d4a; letter-spacing:.02em; }
</style>
</head>
<body>
<div class="wrapper">

  <div class="hero">
    <h1>DAOv2.0</h1>
    <div class="sub">X-ray Reflection Spectroscopy Model &nbsp;·&nbsp; Compton RT &nbsp;·&nbsp; Cloudy &nbsp;·&nbsp; HEASoft</div>
    <div class="author">Yimin Huang &nbsp;·&nbsp; <a href="/docs">Parameter Reference</a> &nbsp;·&nbsp; <a href="/plots" class="link-cyan">Results Viewer</a> &nbsp;·&nbsp; <a href="/convergence" class="link-cyan">Convergence</a></div>
    <hr class="divider">
  </div>

  <div class="card full">
    <h2><span class="icon">◉</span> Corona Spectrum</h2>
    <div class="model-grid" id="modelGrid"></div>
    <div id="coronaParams"></div>
  </div>

  <div class="grid">
    <div class="card">
      <h2><span class="icon">◈</span> Slab Physics</h2>
      <div id="slabParams"></div>
      <hr class="divider" style="max-width:100%">
      <div class="derived" id="derivedBlock"></div>
    </div>
    <div class="card">
      <h2><span class="icon">◎</span> Accretion Disk</h2>
      <div id="diskParams"></div>
    </div>
  </div>

  <div class="card full">
    <h2><span class="icon">⚙</span> Test Mode</h2>
    <div class="toggle-row">
      <label>Enable</label>
      <label class="toggle"><input type="checkbox" id="testRt" onchange="toggleTest()"><span class="slider"></span></label>
    </div>
    <div id="testParams" style="display:none"></div>
  </div>

  <div class="card full">
    <h2><span class="icon">▶</span> Current Configuration</h2>
    <div class="cmd-box" id="cmdBox"><span class="prompt">$ </span><span id="cmdText"></span></div>
    <div class="actions">
      <button class="btn btn-primary" onclick="addToQueue()">Add to Queue</button>
      <button class="btn btn-secondary" onclick="copyCmd()">Copy Command</button>
      <button class="btn btn-secondary" onclick="resetAll()">Reset</button>
    </div>
    <div class="cmd-hint">
      Make sure your runtime environment is set before running
      (e.g. <code>source $HEADAS/headas-init.sh</code>).
    </div>
  </div>

  <div class="card full">
    <h2><span class="icon">☰</span> Run Queue</h2>
    <div id="queueBody"></div>
    <div class="batch-actions">
      <button class="btn btn-primary" onclick="copyAllCmds()">Copy All</button>
      <button class="btn btn-primary" onclick="exportScript()">Export .sh</button>
      <button class="btn btn-secondary" onclick="clearQueue()">Clear Queue</button>
      <span class="batch-count" id="queueCount"></span>
    </div>
  </div>

  <div class="footer">DAOv2.0 &nbsp;·&nbsp; Yimin Huang &nbsp;·&nbsp; Compton scattering radiative transfer</div>
</div>

<div class="toast" id="toast"></div>

<script>
// ── Data from server ──────────────────────
const MODELS = {{ models_json | safe }};
const SLAB   = {{ slab_json | safe }};
const DISK   = {{ disk_json | safe }};
const TEST   = [{"name":"T_test","label":"T<sub>test</sub>","desc":"Slab temperature [K]","default":1e8,"min":1e4,"max":1e10,"step":1e4}];
const DEFAULTS = {nh:15,zeta:3,frac:1,incidence:0.7071067811865476,Afe:1,kT_disk:0.35,T_test:1e8};

let corona = "{{ default_corona }}";
let vals = {};

// ── Build fields ──────────────────────────
function buildField(p, container) {
  const id = 'f_' + p.name;
  const d = document.createElement('div');
  d.innerHTML = `
    <div class="field">
      <label for="${id}">${p.label}</label>
      <input type="range" id="${id}" min="${p.min}" max="${p.max}" step="${p.step}" value="${p.default}"
             oninput="syncVal('${p.name}','${id}')">
      <input type="text" class="val" id="${id}_v" value="${fmtNum(p.default)}"
             onchange="syncSlider('${p.name}','${id}')">
    </div>
    <div class="field-desc">${p.desc}</div>`;
  container.appendChild(d);
  vals[p.name] = p.default;
}

function fmtNum(v) {
  if (Math.abs(v) >= 1e4 || (Math.abs(v) < 0.01 && v !== 0)) return v.toExponential(2);
  return parseFloat(v.toPrecision(6)).toString();
}

function syncVal(name, id) {
  const v = parseFloat(document.getElementById(id).value);
  document.getElementById(id+'_v').value = fmtNum(v);
  vals[name] = v;
  updateCmd();
  updateDerived();
}

function syncSlider(name, id) {
  let v = parseFloat(document.getElementById(id+'_v').value);
  if (isNaN(v)) return;
  document.getElementById(id).value = v;
  vals[name] = v;
  updateCmd();
  updateDerived();
}

// ── Corona model selector ─────────────────
function buildModels() {
  const grid = document.getElementById('modelGrid');
  grid.innerHTML = '';
  Object.entries(MODELS).forEach(([key, m]) => {
    const btn = document.createElement('div');
    btn.className = 'model-btn' + (key === corona ? ' active' : '');
    btn.innerHTML = `<div class="name">${m.label}</div><div class="mdesc">${m.desc}</div>`;
    btn.onclick = () => selectModel(key);
    grid.appendChild(btn);
  });
}

function selectModel(key) {
  corona = key;
  // Remove old corona params from vals
  Object.keys(MODELS).forEach(k => MODELS[k].params.forEach(p => delete vals[p.name]));
  buildModels();
  buildCoronaParams();
  updateCmd();
}

function buildCoronaParams() {
  const el = document.getElementById('coronaParams');
  el.innerHTML = '';
  MODELS[corona].params.forEach(p => buildField(p, el));
}

// ── Derived quantities ────────────────────
function updateDerived() {
  const nh = vals.nh || 15, zeta = vals.zeta || 3;
  const xi = Math.pow(10, zeta);
  const nH = Math.pow(10, nh);
  const Fx = xi * nH / (4 * Math.PI);
  document.getElementById('derivedBlock').innerHTML = `
    <div class="item"><div class="dlabel">ξ</div><div class="dval">${xi.toExponential(2)}</div></div>
    <div class="item"><div class="dlabel">n<sub>H</sub></div><div class="dval">${nH.toExponential(2)}</div></div>
    <div class="item"><div class="dlabel">F<sub>x</sub></div><div class="dval">${Fx.toExponential(2)}</div></div>
  `;
}

// ── Test mode toggle ──────────────────────
function toggleTest() {
  const on = document.getElementById('testRt').checked;
  document.getElementById('testParams').style.display = on ? 'block' : 'none';
  updateCmd();
}

// ── Command builder ───────────────────────
function updateCmd() {
  let parts = ['./maindaocl', '-corona ' + corona];
  MODELS[corona].params.forEach(p => {
    parts.push('-' + p.name + ' ' + fmtNum(vals[p.name] ?? p.default));
  });
  [SLAB, DISK].forEach(group => group.forEach(p => {
    const v = vals[p.name] ?? p.default;
    if (v !== DEFAULTS[p.name]) parts.push('-' + p.name + ' ' + fmtNum(v));
  }));
  if (document.getElementById('testRt').checked) {
    parts.push('-test_rt');
    TEST.forEach(p => {
      const v = vals[p.name] ?? p.default;
      if (v !== DEFAULTS[p.name]) parts.push('-' + p.name + ' ' + fmtNum(v));
    });
  }
  document.getElementById('cmdText').textContent = parts.join(' ');
}

// ── Actions ───────────────────────────────
function clipCopy(text) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(() => showToast('Copied to clipboard')).catch(() => fallbackCopy(text));
  } else {
    fallbackCopy(text);
  }
}
function fallbackCopy(text) {
  const ta = document.createElement('textarea');
  ta.value = text; ta.style.cssText = 'position:fixed;left:-9999px';
  document.body.appendChild(ta); ta.select();
  document.execCommand('copy');
  document.body.removeChild(ta);
  showToast('Copied to clipboard');
}

function copyCmd() {
  clipCopy(document.getElementById('cmdText').textContent);
}

function resetAll() {
  corona = "{{ default_corona }}";
  vals = {};
  init();
  showToast('Reset to defaults');
}

function showToast(msg) {
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.classList.add('show');
  setTimeout(() => t.classList.remove('show'), 2000);
}

// ── Batch queue ──────────────────────────
let queue = [];

function getCmd() {
  return document.getElementById('cmdText').textContent;
}

function getSnapshot() {
  return { corona, vals: {...vals}, testRt: document.getElementById('testRt').checked };
}

function cmdSummary(cmd) {
  // Extract key params for display
  const m = cmd.match(/-corona\s+(\S+)/);
  const model = m ? m[1] : '?';
  const rest = cmd.replace('./maindaocl', '').replace(/-corona\s+\S+/, '').trim();
  return { model, rest };
}

function addToQueue() {
  const cmd = getCmd();
  const snap = getSnapshot();
  queue.push({ cmd, snap });
  renderQueue();
  showToast(`Added to queue (${queue.length} total)`);
}

function removeFromQueue(idx) {
  queue.splice(idx, 1);
  renderQueue();
}

function loadFromQueue(idx) {
  const snap = queue[idx].snap;
  corona = snap.corona;
  vals = {...snap.vals};
  document.getElementById('testRt').checked = snap.testRt;
  toggleTest();
  // Rebuild UI with snapshot values
  buildModels();
  const cp = document.getElementById('coronaParams');
  cp.innerHTML = '';
  MODELS[corona].params.forEach(p => {
    buildField(p, cp);
    const id = 'f_' + p.name;
    const v = vals[p.name] ?? p.default;
    document.getElementById(id).value = v;
    document.getElementById(id + '_v').value = fmtNum(v);
  });
  ['slabParams','diskParams','testParams'].forEach(sec => {
    const el = document.getElementById(sec);
    el.innerHTML = '';
    const list = sec === 'slabParams' ? SLAB : sec === 'diskParams' ? DISK : TEST;
    list.forEach(p => {
      buildField(p, el);
      const id = 'f_' + p.name;
      const v = vals[p.name] ?? p.default;
      document.getElementById(id).value = v;
      document.getElementById(id + '_v').value = fmtNum(v);
    });
  });
  updateDerived();
  updateCmd();
  showToast(`Loaded run #${idx + 1}`);
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
      <td class="q-model">${s.model}</td>
      <td class="q-params">${s.rest}</td>
      <td class="q-actions">
        <button class="q-btn q-btn-load" onclick="loadFromQueue(${i})" title="Load into editor">↩</button>
        <button class="q-btn q-btn-del" onclick="removeFromQueue(${i})" title="Remove">✕</button>
      </td>
    </tr>`;
  });
  html += '</tbody></table>';
  body.innerHTML = html;
}

function copyAllCmds() {
  if (queue.length === 0) { showToast('Queue is empty'); return; }
  clipCopy(queue.map(q => q.cmd).join('\n'));
}

function exportScript() {
  if (queue.length === 0) { showToast('Queue is empty'); return; }
  let script = '#!/bin/bash\n# DAOv2.0 batch run — ' + new Date().toISOString().slice(0,10) + '\nset -e\n\n';
  script += '# Set your runtime environment before running, e.g.:\n';
  script += '# export HEADAS=/path/to/heasoft/arch\n';
  script += '# source $HEADAS/headas-init.sh\n\n';
  queue.forEach((q, i) => {
    script += `echo "=== Run ${i+1}/${queue.length} ==="\n${q.cmd}\n\n`;
  });
  const blob = new Blob([script], {type: 'text/x-shellscript'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'dao_batch.sh';
  a.click();
  URL.revokeObjectURL(a.href);
  showToast('Exported dao_batch.sh');
}

function clearQueue() {
  queue = [];
  renderQueue();
  showToast('Queue cleared');
}

// ── Init ──────────────────────────────────
function init() {
  buildModels();
  buildCoronaParams();
  ['slabParams','diskParams','testParams'].forEach(id => document.getElementById(id).innerHTML = '');
  SLAB.forEach(p => buildField(p, document.getElementById('slabParams')));
  DISK.forEach(p => buildField(p, document.getElementById('diskParams')));
  TEST.forEach(p => buildField(p, document.getElementById('testParams')));
  updateDerived();
  updateCmd();
  renderQueue();
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
<title>DAOv2.0 — Parameter Reference</title>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  :root {
    --bg: #0a0c10; --bg2: #12151c; --card: #161a24; --border: #232838;
    --accent: #f0a030; --accent2: #e05828; --cyan: #38bdf8;
    --text: #c8cdd8; --dim: #6b7280; --white: #eef0f4;
  }
  body {
    font-family: 'Inter', 'Helvetica Neue', 'Segoe UI', system-ui, sans-serif;
    background: var(--bg); color: var(--text); min-height: 100vh; font-size: 15px;
  }
  body::before {
    content:''; position:fixed; inset:0; z-index:0;
    background: url('/image/bg_blackhole.png') center top / cover no-repeat fixed;
    opacity: 0.15; pointer-events: none;
  }
  .wrap { position:relative; z-index:1; max-width:900px; margin:0 auto; padding:30px 20px; }

  .back {
    display:inline-block; color:var(--accent); text-decoration:none; font-size:.82rem;
    margin-bottom:20px; padding:6px 14px; border:1px solid var(--border); border-radius:6px;
    transition: border-color .2s;
  }
  .back:hover { border-color: var(--accent); }

  h1 {
    font-size:1.6rem; font-weight:700; letter-spacing:.12em;
    background: linear-gradient(90deg, #f0a030, #e05828);
    -webkit-background-clip:text; -webkit-text-fill-color:transparent;
    margin-bottom:6px;
  }
  .subtitle { font-size:.82rem; color:var(--dim); margin-bottom:30px; }

  .section { margin-bottom: 36px; }
  .section h2 {
    font-size:.9rem; text-transform:uppercase; letter-spacing:.1em;
    color:var(--accent); margin-bottom:14px; padding-bottom:6px;
    border-bottom:1px solid var(--border);
  }

  table { width:100%; border-collapse:collapse; margin-bottom:8px; }
  th {
    text-align:left; font-size:.7rem; text-transform:uppercase; letter-spacing:.08em;
    color:var(--dim); padding:8px 12px; border-bottom:1px solid var(--border);
  }
  td { padding:10px 12px; border-bottom:1px solid #1a1e28; font-size:.78rem; vertical-align:top; }
  tr:hover td { background: rgba(240,160,48,.03); }
  .p-flag { color:var(--cyan); font-weight:600; white-space:nowrap; }
  .p-type { color:var(--dim); font-size:.7rem; }
  .p-default { color:var(--white); font-weight:500; }
  .p-desc { color:var(--text); line-height:1.5; }
  .p-note { color:var(--dim); font-size:.7rem; font-style:italic; }

  .formula-box {
    background:var(--bg2); border:1px solid var(--border); border-radius:8px;
    padding:14px 18px; margin:10px 0 16px; font-size:.85rem; color:var(--white);
    line-height:1.8;
  }
  .formula-box .flabel { color:var(--cyan); font-weight:600; margin-right:8px; }

  .note-box {
    background: rgba(240,160,48,.06); border-left:3px solid var(--accent);
    padding:12px 16px; border-radius:0 8px 8px 0; margin:14px 0;
    font-size:.78rem; color:var(--text); line-height:1.6;
  }
  .note-box strong { color:var(--accent); }

  .tag {
    display:inline-block; font-size:.6rem; padding:2px 6px; border-radius:4px;
    font-weight:600; letter-spacing:.04em; vertical-align:middle; margin-left:4px;
  }
  .tag-required { background:rgba(239,68,68,.15); color:#ef4444; }
  .tag-optional { background:rgba(56,189,248,.1); color:var(--cyan); }
  .tag-model { background:rgba(240,160,48,.12); color:var(--accent); }

  .footer { text-align:center; padding:30px 0 12px; font-size:.65rem; color:#333; }
</style>
</head>
<body>
<div class="wrap">

  <a class="back" href="/">← Back to Configurator</a>

  <h1>PARAMETER REFERENCE</h1>
  <div class="subtitle">DAOv2.0 — X-ray Reflection Spectroscopy Model &nbsp;·&nbsp; Yimin Huang</div>

  <!-- ── Corona Models ────────────────────── -->
  <div class="section">
    <h2>Corona Spectrum Models</h2>

    <p style="font-size:.78rem;color:var(--text);margin-bottom:14px;line-height:1.6;">
      The corona model defines the illuminating X-ray continuum incident on the slab.
      Select with <span class="p-flag">-corona &lt;model&gt;</span>.
      Each model requires specific parameters listed below.
    </p>

    <div class="formula-box">
      <div><span class="flabel">powerlaw</span>
        F(E) ∝ E<sup>−Γ</sup> · exp(−E<sub>lo</sub> / E)</div>
      <div><span class="flabel">cutoffpl</span>
        F(E) ∝ E<sup>−Γ</sup> · exp(−E / E<sub>cut</sub>) · exp(−E<sub>lo</sub> / E)</div>
      <div><span class="flabel">nthcomp</span>
        Thermal Comptonisation — Zdziarski, Johnson & Magdziarz (1996).
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
        <td class="p-desc">Photon index Γ of the power-law continuum. Typical range: 1.4–3.0 for AGN.</td>
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
        <td class="p-flag">-kT_e</td><td class="p-type">float</td><td class="p-default">60</td>
        <td><span class="tag tag-model">nthcomp</span><span class="tag tag-model">comptt</span></td>
        <td class="p-desc">Electron / plasma temperature [keV] of the Comptonising corona.</td>
      </tr>
      <tr>
        <td class="p-flag">-kT_bb</td><td class="p-type">float</td><td class="p-default">0.1</td>
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
  <div class="section">
    <h2>Slab Physics</h2>

    <p style="font-size:.78rem;color:var(--text);margin-bottom:14px;line-height:1.6;">
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
        <td class="p-desc">Ionisation parameter exponent: ξ = 10<sup>ζ</sup> erg cm s<sup>−1</sup>.
          Controls the ionisation state of the slab. Low ξ → neutral (cold reflection);
          high ξ → highly ionised (Compton-dominated).
          <div class="p-note">ξ = 4π F<sub>x</sub> / n<sub>H</sub> (Tarter+ 1969 definition).</div></td>
      </tr>
      <tr>
        <td class="p-flag">-frac</td><td class="p-type">float</td><td class="p-default">100</td>
        <td class="p-desc">Flux ratio f = F<sub>corona</sub> / F<sub>disk</sub>. Controls relative
          normalisation of the illuminating corona vs the thermal disk emission.
          If frac ≤ 0, no disk illumination (corona only, there are some problems now, please keep frac>0).</td>
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
      ξ = 10<sup>ζ</sup>, &nbsp;
      n<sub>H</sub> = 10<sup>nh</sup>, &nbsp;
      F<sub>x</sub> = ξ · n<sub>H</sub> / 4π &nbsp;
      (the illuminating flux that normalises the corona + disk spectra).
    </div>
  </div>

  <!-- ── Accretion Disk ───────────────────── -->
  <div class="section">
    <h2>Accretion Disk</h2>

    <table>
      <tr><th>Flag</th><th>Type</th><th>Default</th><th>Description</th></tr>
      <tr>
        <td class="p-flag">-kT_disk</td><td class="p-type">float</td><td class="p-default">0.35</td>
        <td class="p-desc">Disk blackbody temperature [eV]. This sets the peak of the thermal
          multi-colour disk emission that illuminates the slab from below.
          Typical: 0.1–1.0 eV for AGN; higher for X-ray binaries.</td>
      </tr>
    </table>
  </div>

  <!-- ── Test Mode ────────────────────────── -->
  <div class="section">
    <h2>Test Mode</h2>

    <p style="font-size:.78rem;color:var(--text);margin-bottom:14px;line-height:1.6;">
      Test mode bypasses Cloudy and uses a synthetic uniform-temperature slab with
      analytic opacities. Useful for validating the RT solver and Compton kernel.
    </p>

    <table>
      <tr><th>Flag</th><th>Type</th><th>Default</th><th>Description</th></tr>
      <tr>
        <td class="p-flag">-test_rt</td><td class="p-type">flag</td><td class="p-default">off</td>
        <td class="p-desc">Enable test mode. No Cloudy calls — uses a uniform slab with
          temperature T<sub>test</sub> and Thomson opacity.</td>
      </tr>
      <tr>
        <td class="p-flag">-T_test</td><td class="p-type">float</td><td class="p-default">1e8</td>
        <td class="p-desc">Uniform slab temperature [K] in test mode. Controls the Compton
          scattering redistribution width (Δε/ε ∼ 4kT/m<sub>e</sub>c<sup>2</sup>).</td>
      </tr>
    </table>
  </div>

  <!-- ── Execution flow ───────────────────── -->
  <div class="section">
    <h2>Execution Flow</h2>

    <div class="note-box" style="background:rgba(56,189,248,.05);border-left-color:var(--cyan);">
      <strong style="color:var(--cyan);">Production mode:</strong><br>
      1. Parse CLI → <code>ModelParams</code><br>
      2. Initialise grids (angle, depth, energy)<br>
      3. Compute corona + disk illumination spectra<br>
      4. Precompute Compton kernel + scattering cross-sections (cached on disk)<br>
      5. Outer loop: Cloudy depth sweep → extract j<sub>ν</sub>, κ<sub>abs</sub>, κ<sub>sct</sub> → RT solve → update J, ξ, T → check convergence<br>
      6. Save final spectra to <code>data/</code>
    </div>

    <div class="note-box" style="background:rgba(56,189,248,.05);border-left-color:var(--cyan);">
      <strong style="color:var(--cyan);">Test mode:</strong><br>
      1. Parse CLI → <code>ModelParams</code><br>
      2. Initialise grids with synthetic 1000-bin log-spaced energy grid<br>
      3. Compute illumination<br>
      4. Single RT solve with uniform temperature + Thomson opacity<br>
      5. Save results
    </div>
  </div>

  <!-- ── Example Commands ─────────────────── -->
  <div class="section">
    <h2>Example Commands</h2>

    <div class="formula-box" style="font-size:.78rem;line-height:2;">
      <div><span class="p-flag">./maindaocl</span> -corona cutoffpl -Gamma 2.0 -Ecut 300 -nh 16 -zeta 4 -frac 0.5</div>
      <div><span class="p-flag">./maindaocl</span> -corona nthcomp -Gamma 2.0 -kT_e 100 -kT_bb 0.05 -nh 15 -zeta 3</div>
      <div><span class="p-flag">./maindaocl</span> -corona comptt -kT_e 50 -kT_bb 0.05 -taup 1.0 -Afe 3.0</div>
      <div><span class="p-flag">./maindaocl</span> -test_rt -corona nthcomp -Gamma 2.0 -kT_e 100 -kT_bb 0.05</div>
    </div>
  </div>

  <div class="footer">DAOv2.0 &nbsp;·&nbsp; Yimin Huang &nbsp;·&nbsp; X-ray reflection model</div>
</div>
</body>
</html>
"""


@app.route("/docs")
def docs():
    return render_template_string(DOCS_HTML)


# ─── Plots page ─────────────────────────────────────────────────
import glob as globmod
import re

@app.route("/api/runs")
def api_runs():
    """Scan results/<hash>/ directories, read params.json from each."""
    results_dir = os.path.join(WORK_DIR, "results")
    runs = []
    if os.path.isdir(results_dir):
        for name in sorted(os.listdir(results_dir)):
            run_path = os.path.join(results_dir, name)
            pj = os.path.join(run_path, "params.json")
            if os.path.isdir(run_path) and os.path.exists(pj):
                try:
                    with open(pj) as f:
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
                runs.append({"hash": name, "iters": iters, "meta": meta})
    return jsonify(runs=runs)


@app.route("/api/data/<run_hash>/<int:iteration>")
def api_data(run_hash, iteration):
    """Parse emergent, moments, profile data for a run/iteration."""
    rd = os.path.join(WORK_DIR, "results", run_hash)
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
                if line.startswith("# Fe lines") and f"iter={iteration}" in line:
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


@app.route("/api/profiles/<run_hash>")
def api_all_profiles(run_hash):
    """Return temperature profiles for iterations of a run (every `step`-th)."""
    step = int(request.args.get("step", 2))
    rd = os.path.join(WORK_DIR, "results", run_hash)
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


@app.route("/api/mean_intensity/<run_hash>")
def api_all_mean_intensity(run_hash):
    """Return surface (depth=0) mean intensity J0(E) for iterations (every `step`-th)."""
    step = int(request.args.get("step", 2))
    rd = os.path.join(WORK_DIR, "results", run_hash)
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


@app.route("/api/line_labels/<run_hash>")
def api_line_labels(run_hash):
    """Return line labels from the first available line_labels file."""
    import re as _re
    rd = os.path.join(WORK_DIR, "results", run_hash)
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
<title>DAOv2.0 — Results Viewer</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  :root {
    --bg: #0a0c10; --bg2: #12151c; --card: #161a24; --border: #232838;
    --accent: #f0a030; --accent2: #e05828; --cyan: #38bdf8;
    --text: #c8cdd8; --dim: #6b7280; --white: #eef0f4;
    --green: #22c55e;
  }
  body {
    font-family: 'Inter', 'Helvetica Neue', 'Segoe UI', system-ui, sans-serif;
    background: var(--bg); color: var(--text); min-height: 100vh; font-size: 15px;
  }
  body::before {
    content:''; position:fixed; inset:0; z-index:0;
    background: url('/image/bg_blackhole.png') center top / cover no-repeat fixed;
    opacity: 0.1; pointer-events: none;
  }
  .wrap { position:relative; z-index:1; max-width:1200px; margin:0 auto; padding:20px; }

  .top-bar {
    display:flex; align-items:center; gap:16px; margin-bottom:20px; flex-wrap:wrap;
  }
  .back {
    color:var(--accent); text-decoration:none; font-size:.82rem;
    padding:6px 14px; border:1px solid var(--border); border-radius:6px;
  }
  .back:hover { border-color: var(--accent); }
  .top-bar h1 {
    font-size:1.3rem; font-weight:700; letter-spacing:.1em;
    background: linear-gradient(90deg, #f0a030, #e05828);
    -webkit-background-clip:text; -webkit-text-fill-color:transparent;
  }

  .controls {
    display:flex; gap:12px; align-items:center; flex-wrap:wrap;
    margin-bottom:20px;
  }
  .controls label { font-size:.78rem; color:var(--dim); }
  .controls select {
    background:var(--bg2); color:var(--white); border:1px solid var(--border);
    border-radius:6px; padding:6px 12px; font-family:inherit; font-size:.8rem;
    outline:none; min-width:180px;
  }
  .controls select:focus { border-color:var(--accent); }

  .params-box {
    background:var(--card); border:1px solid var(--border); border-radius:8px;
    padding:12px 16px; margin-bottom:20px; font-size:.75rem; color:var(--text);
    line-height:1.8;
  }
  .params-box .pk { color:var(--cyan); }
  .params-box .pv { color:var(--white); font-weight:600; }

  .plot-card {
    background:var(--card); border:1px solid var(--border); border-radius:10px;
    padding:16px; margin-bottom:20px;
  }
  .plot-card h2 {
    font-size:.82rem; text-transform:uppercase; letter-spacing:.1em;
    color:var(--accent); margin-bottom:10px;
  }
  .plot-area { width:100%; height:480px; }

  .status {
    text-align:center; padding:40px; color:var(--dim); font-size:.85rem;
  }

  .footer { text-align:center; padding:20px 0 10px; font-size:.65rem; color:#333; }
</style>
</head>
<body>
<div class="wrap">

  <div class="top-bar">
    <a class="back" href="/">← Configurator</a>
    <h1>RESULTS VIEWER</h1>
  </div>

  <div class="controls">
    <label>Run:</label>
    <select id="runSelect" onchange="onRunChange()">
      <option value="">— select a run —</option>
    </select>
    <label>Iteration:</label>
    <select id="iterSelect" onchange="loadData()">
      <option value="">—</option>
    </select>
  </div>

  <div class="params-box" id="paramsBox" style="display:none"></div>

  <div id="plotsContainer" style="display:none">
    <div class="plot-card">
      <h2 style="display:flex;align-items:center;gap:12px;">
        1. Emergent Intensity at Surface (outgoing &mu; &gt; 0)
        <button id="feToggle" onclick="toggleFeLines()" style="font-size:.72rem;padding:5px 14px;border-radius:20px;border:1px solid rgba(255,100,100,.3);background:rgba(255,100,100,.08);color:#ff8888;cursor:pointer;font-family:inherit;font-weight:500;letter-spacing:.03em;transition:all .2s;">Fe lines: ON</button>
      </h2>
      <div class="plot-area" id="plotEmergent"></div>
    </div>
    <div class="plot-card">
      <h2 style="display:flex;align-items:center;gap:12px;">
        2. Mean Outgoing Intensity
        <button id="btnLineLabels" onclick="toggleLineLabels()"
                style="font-size:.72rem;padding:5px 14px;border-radius:20px;
                       border:1px solid rgba(100,180,255,.3);background:rgba(100,180,255,.08);
                       color:#64b4ff;cursor:pointer;font-family:inherit;font-weight:500;
                       letter-spacing:.03em;transition:all .2s;">
          Line IDs: OFF
        </button>
        <span id="lineLabelStatus" style="font-size:.68rem;color:var(--dim);"></span>
      </h2>
      <div class="plot-area" id="plotMean"></div>
    </div>
    <div class="plot-card">
      <h2>3. Temperature Profile</h2>
      <div class="plot-area" id="plotProfile"></div>
    </div>
  </div>

  <div class="status" id="statusMsg">Select a run to view results.</div>

  <div class="footer">DAOv2.0 &nbsp;·&nbsp; Yimin Huang</div>
</div>

<script>
const PLOT_BG = '#161a24';
const GRID_COLOR = '#232838';
const FONT_COLOR = '#c8cdd8';
const LAYOUT_BASE = {
  paper_bgcolor: '#0a0c10', plot_bgcolor: PLOT_BG,
  font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
  margin: { l:70, r:30, t:30, b:60 },
  legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:10} },
  xaxis: { gridcolor: GRID_COLOR, zerolinecolor: GRID_COLOR },
  yaxis: { gridcolor: GRID_COLOR, zerolinecolor: GRID_COLOR },
};

const COLORS = ['#f0a030','#e05828','#38bdf8','#22c55e','#a78bfa','#f472b6','#facc15','#67e8f9'];

let allRuns = [];

async function init() {
  const r = await fetch('/api/runs');
  const d = await r.json();
  allRuns = d.runs;
  const sel = document.getElementById('runSelect');
  allRuns.forEach(run => {
    const o = document.createElement('option');
    o.value = run.hash;
    const m = run.meta;
    let label = `${run.hash}`;
    if (m.corona) label += ` — ${m.corona}`;
    if (m.nh !== undefined) label += `  nh=${m.nh}`;
    if (m.zeta !== undefined) label += `  z=${m.zeta}`;
    if (m.Gamma !== undefined && m.Gamma > 0) label += `  G=${m.Gamma}`;
    o.textContent = label;
    sel.appendChild(o);
  });
}

function onRunChange() {
  const rid = document.getElementById('runSelect').value;
  const run = allRuns.find(r => r.hash === rid);
  const iterSel = document.getElementById('iterSelect');
  iterSel.innerHTML = '';
  if (!run) {
    document.getElementById('plotsContainer').style.display = 'none';
    document.getElementById('paramsBox').style.display = 'none';
    document.getElementById('statusMsg').style.display = 'block';
    document.getElementById('statusMsg').textContent = 'Select a run to view results.';
    return;
  }
  run.iters.forEach(it => {
    const o = document.createElement('option');
    o.value = it;
    o.textContent = `Iteration ${it}`;
    iterSel.appendChild(o);
  });
  // Select last iteration
  if (run.iters.length > 0) iterSel.value = run.iters[run.iters.length - 1];
  // Show params
  showParams(run.meta);
  loadData();
}

function showParams(meta) {
  const box = document.getElementById('paramsBox');
  if (!meta || !meta.corona) { box.style.display = 'none'; return; }
  const skip = new Set(['hash']);
  let html = '';
  for (const [k,v] of Object.entries(meta)) {
    if (skip.has(k)) continue;
    const fv = typeof v === 'number' ? (Math.abs(v)>=1e4||(Math.abs(v)<0.01&&v!==0) ? v.toExponential(2) : v) : v;
    html += `<span class="pk">${k}</span>=<span class="pv">${fv}</span> &nbsp; `;
  }
  box.innerHTML = html;
  box.style.display = 'block';
}

async function loadData() {
  const rid = document.getElementById('runSelect').value;
  const iter = document.getElementById('iterSelect').value;
  if (!rid || !iter) return;
  document.getElementById('statusMsg').style.display = 'block';
  document.getElementById('statusMsg').textContent = 'Loading data...';
  document.getElementById('plotsContainer').style.display = 'none';

  const r = await fetch(`/api/data/${rid}/${iter}`);
  const d = await r.json();

  if (!d.emergent && !d.profile) {
    document.getElementById('statusMsg').textContent = 'No data files found for this run/iteration.';
    return;
  }

  document.getElementById('statusMsg').style.display = 'none';
  document.getElementById('plotsContainer').style.display = 'block';

  // Get incidence angle from run metadata
  const run = allRuns.find(r => r.hash === rid);
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
  btn.style.color = showFeLines ? '#ff8888' : '#555';
  btn.style.borderColor = showFeLines ? 'rgba(255,100,100,.3)' : 'rgba(100,100,100,.2)';
  btn.style.background = showFeLines ? 'rgba(255,100,100,.08)' : 'rgba(100,100,100,.05)';
  plotEmergent(lastEmData, lastMuInc, lastFeLines);
}

function autoLogRange(traces) {
  // Find peak across all traces, set range to [peak*1e-6, peak*5] (like plot_results.py)
  let peak = -Infinity;
  traces.forEach(t => {
    t.y.forEach(v => { if (v > 0 && v > peak) peak = v; });
  });
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

  // I_corona scaled: 2 * I_corona / mu_inc
  traces.push({
    x: E_keV, y: E_eV.map((e, i) => e * 2.0 * em.I_corona[i] / mu_inc),
    name: 'E × 2I_cor/μ_inc', mode: 'lines',
    line: { color: '#888', width: 1.5, dash: 'dash' }
  });
  // I_disk scaled: I_disk / 2
  traces.push({
    x: E_keV, y: E_eV.map((e, i) => e * em.I_disk[i] / 2.0),
    name: 'E × I_disk/2', mode: 'lines',
    line: { color: '#555', width: 1.5, dash: 'dot' }
  });

  // Outgoing angles (mu > 0): E[eV] * I
  let ci = 0;
  em.mu_vals.forEach(mu => {
    if (mu > 0) {
      const vals = em.angles[mu.toFixed(4)];
      traces.push({
        x: E_keV, y: E_eV.map((e, i) => e * vals[i]),
        name: `μ = ${mu.toFixed(4)}`, mode: 'lines',
        line: { color: COLORS[ci % COLORS.length], width: 1.8 }
      });
      ci++;
    }
  });

  const yr = autoLogRange(traces);

  // Fe K line vertical markers + labels
  const shapes = [];
  const annotations = [];
  if (showFeLines && fe_lines && fe_lines.length > 0 && yr) {
    // Show only the 6 strongest lines by relint
    const top6 = [...fe_lines].sort((a, b) => b.relint - a.relint).slice(0, 6);
    top6.forEach((fl, idx) => {
      shapes.push({
        type: 'line', xref: 'x', yref: 'paper',
        x0: fl.E_keV, x1: fl.E_keV, y0: 0, y1: 1,
        line: { color: 'rgba(255,100,100,0.5)', width: 1, dash: 'dot' }
      });
      annotations.push({
        x: Math.log10(fl.E_keV), xref: 'x', yref: 'paper',
        y: 1.0 - idx * 0.06, text: fl.label,
        showarrow: false, font: { color: '#ff8888', size: 11 },
        xanchor: 'left', xshift: 4
      });
    });
  }

  Plotly.react(div, traces, {
    paper_bgcolor: '#0a0c10', plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:40, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:10} },
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'E × I [erg² cm⁻² s⁻¹ erg⁻¹ sr⁻¹]', range:yr, gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
    shapes: shapes,
    annotations: annotations,
  }, {responsive:true});
}

function plotMeanOutgoing(em) {
  const div = document.getElementById('plotMean');
  if (!em) { div.innerHTML = '<div class="status">No emergent data.</div>'; return; }

  const E_eV = em.E;
  const E_keV = E_eV.map(e => e / 1e3);

  // Compute mean over outgoing angles (mu > 0)
  const outgoing_keys = em.mu_vals.filter(mu => mu > 0).map(mu => mu.toFixed(4));
  const n_out = outgoing_keys.length;
  const mean_I = E_eV.map((_, ie) => {
    let sum = 0;
    outgoing_keys.forEach(k => { sum += em.angles[k][ie]; });
    return n_out > 0 ? sum / n_out : 0;
  });

  const traces = [
    {
      x: E_keV, y: E_eV.map((e, i) => e * mean_I[i]),
      name: 'E × Mean outgoing I', mode: 'lines',
      line: { color: '#f0a030', width: 2 }
    },
  ];

  const yr = autoLogRange(traces);
  Plotly.react(div, traces, {
    paper_bgcolor: '#0a0c10', plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:30, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:10} },
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'E × I [erg² cm⁻² s⁻¹ erg⁻¹ sr⁻¹]', range:yr, gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

function plotProfile(prof) {
  const div = document.getElementById('plotProfile');
  if (!prof) { div.innerHTML = '<div class="status">No profile data.</div>'; return; }

  const traces = [{
    x: prof.tau_mid, y: prof.T_K,
    name: 'T', mode: 'lines+markers',
    line: { color: '#e05828', width: 2 },
    marker: { size: 4, color: '#f0a030' }
  }];

  Plotly.react(div, traces, {
    paper_bgcolor: '#0a0c10', plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:30, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:10} },
    xaxis: { type:'log', title:'τ (Thomson)', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'T [K]', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

// ── Line label overlay on Mean Outgoing Intensity plot ─────────
let lineLabelsOn = false;
let lineLabelsData = null;
const ROMAN = {1:'I',2:'II',3:'III',4:'IV',5:'V',6:'VI',7:'VII',8:'VIII',
               9:'IX',10:'X',11:'XI',12:'XII',13:'XIII',14:'XIV',15:'XV',
               16:'XVI',17:'XVII',18:'XVIII',19:'XIX',20:'XX',21:'XXI',
               22:'XXII',23:'XXIII',24:'XXIV',25:'XXV',26:'XXVI',27:'XXVII'};

function fmtSpecies(lab) {
  const m = lab.match(/^([A-Z][a-z]?)\s*(\d+)$/);
  if (m) { const n = parseInt(m[2]); if (ROMAN[n]) return m[1]+'\u2009'+ROMAN[n]; }
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
    btn.style.borderColor = 'rgba(100,180,255,.3)';
    btn.style.background = 'rgba(100,180,255,.08)';
    status.textContent = '';
    Plotly.relayout(div, { annotations: [], shapes: [] });
    return;
  }

  if (!lineLabelsData || lineLabelsData.hash !== rid) {
    status.textContent = 'Loading...';
    const r = await fetch(`/api/line_labels/${rid}`);
    const d = await r.json();
    lineLabelsData = { hash: rid, lines: d.lines || [] };
  }

  if (lineLabelsData.lines.length === 0) {
    status.textContent = 'No line labels found.';
    return;
  }

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
    const col = isFluor ? 'rgba(240,80,80,0.85)' : 'rgba(100,180,255,0.85)';
    const lcol = isFluor ? 'rgba(240,80,80,0.2)' : 'rgba(100,180,255,0.15)';

    shapes.push({
      type:'line', xref:'x', yref:'paper',
      x0:keV, x1:keV, y0:0, y1:1,
      line:{color:lcol, width:0.8},
    });
    annotations.push({
      x:Math.log10(keV), y:1, xref:'x', yref:'paper',
      text:name, showarrow:false,
      font:{size:9, color:col, family:'Inter, sans-serif'},
      textangle:-90, xanchor:'left', yanchor:'top', yshift:-4,
    });
  }

  lineLabelsOn = true;
  btn.textContent = 'Line IDs: ON';
  btn.style.borderColor = 'rgba(100,180,255,.6)';
  btn.style.background = 'rgba(100,180,255,.18)';
  status.textContent = `${annotations.length} lines`;

  Plotly.relayout(div, { annotations, shapes });
}

init();
</script>
</body>
</html>
"""

@app.route("/plots")
def plots():
    return render_template_string(PLOTS_HTML)


# ─── Convergence page ────────────────────────────────────────────
CONVERGENCE_HTML = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<title>DAOv2.0 — Temperature Convergence</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  :root {
    --bg: #0a0c10; --bg2: #12151c; --card: #161a24; --border: #232838;
    --accent: #f0a030; --accent2: #e05828; --cyan: #38bdf8;
    --text: #c8cdd8; --dim: #6b7280; --white: #eef0f4;
  }
  body {
    font-family: 'Inter', 'Helvetica Neue', 'Segoe UI', system-ui, sans-serif;
    background: var(--bg); color: var(--text); min-height: 100vh; font-size: 15px;
  }
  body::before {
    content:''; position:fixed; inset:0; z-index:0;
    background: url('/image/bg_blackhole.png') center top / cover no-repeat fixed;
    opacity: 0.1; pointer-events: none;
  }
  .wrap { position:relative; z-index:1; max-width:1200px; margin:0 auto; padding:20px; }

  .top-bar {
    display:flex; align-items:center; gap:16px; margin-bottom:20px; flex-wrap:wrap;
  }
  .back {
    color:var(--accent); text-decoration:none; font-size:.82rem;
    padding:6px 14px; border:1px solid var(--border); border-radius:6px;
  }
  .back:hover { border-color: var(--accent); }
  .top-bar h1 {
    font-size:1.3rem; font-weight:700; letter-spacing:.1em;
    background: linear-gradient(90deg, #f0a030, #e05828);
    -webkit-background-clip:text; -webkit-text-fill-color:transparent;
  }

  .controls {
    display:flex; gap:12px; align-items:center; flex-wrap:wrap;
    margin-bottom:20px;
  }
  .controls label { font-size:.78rem; color:var(--dim); }
  .controls select {
    background:var(--bg2); color:var(--white); border:1px solid var(--border);
    border-radius:6px; padding:6px 12px; font-family:inherit; font-size:.8rem;
    outline:none; min-width:180px;
  }
  .controls select:focus { border-color:var(--accent); }

  .params-box {
    background:var(--card); border:1px solid var(--border); border-radius:8px;
    padding:12px 16px; margin-bottom:20px; font-size:.75rem; color:var(--text);
    line-height:1.8;
  }
  .params-box .pk { color:var(--cyan); }
  .params-box .pv { color:var(--white); font-weight:600; }

  .plot-card {
    background:var(--card); border:1px solid var(--border); border-radius:10px;
    padding:16px; margin-bottom:20px;
  }
  .plot-card h2 {
    font-size:.82rem; text-transform:uppercase; letter-spacing:.1em;
    color:var(--accent); margin-bottom:10px;
  }
  .plot-area { width:100%; height:540px; }

  .status {
    text-align:center; padding:40px; color:var(--dim); font-size:.85rem;
  }

  .footer { text-align:center; padding:20px 0 10px; font-size:.65rem; color:#333; }
</style>
</head>
<body>
<div class="wrap">

  <div class="top-bar">
    <a class="back" href="/">← Configurator</a>
    <a class="back" href="/plots">← Results Viewer</a>
    <h1>TEMPERATURE CONVERGENCE</h1>
  </div>

  <div class="controls">
    <label>Run:</label>
    <select id="runSelect" onchange="onRunChange()">
      <option value="">— select a run —</option>
    </select>
  </div>

  <div class="params-box" id="paramsBox" style="display:none"></div>

  <div id="plotsContainer" style="display:none">
    <div class="plot-card">
      <h2>1. Temperature Profile — All Iterations</h2>
      <div class="plot-area" id="plotTempIter"></div>
    </div>
    <div class="plot-card">
      <h2>2. Surface Mean Intensity J<sub>0</sub>(E) — All Iterations</h2>
      <div class="plot-area" id="plotJ0Iter"></div>
    </div>
  </div>

  <div class="status" id="statusMsg">Select a run to view temperature convergence.</div>

  <div class="footer">DAOv2.0 &nbsp;·&nbsp; Yimin Huang</div>
</div>

<script>
const PLOT_BG = '#161a24';
const GRID_COLOR = '#232838';
const FONT_COLOR = '#c8cdd8';
const COLORS_POOL = ['#f0a030','#e05828','#38bdf8','#22c55e','#a78bfa','#f472b6','#facc15','#67e8f9'];

let allRuns = [];

async function init() {
  const r = await fetch('/api/runs');
  const d = await r.json();
  allRuns = d.runs;
  const sel = document.getElementById('runSelect');
  allRuns.forEach(run => {
    const o = document.createElement('option');
    o.value = run.hash;
    const m = run.meta;
    let label = `${run.hash}`;
    if (m.corona) label += ` — ${m.corona}`;
    if (m.nh !== undefined) label += `  nh=${m.nh}`;
    if (m.zeta !== undefined) label += `  z=${m.zeta}`;
    if (m.Gamma !== undefined && m.Gamma > 0) label += `  G=${m.Gamma}`;
    o.textContent = label;
    sel.appendChild(o);
  });
}

function onRunChange() {
  const rid = document.getElementById('runSelect').value;
  const run = allRuns.find(r => r.hash === rid);
  if (!run) {
    document.getElementById('plotsContainer').style.display = 'none';
    document.getElementById('paramsBox').style.display = 'none';
    document.getElementById('statusMsg').style.display = 'block';
    document.getElementById('statusMsg').textContent = 'Select a run to view temperature convergence.';
    return;
  }
  showParams(run.meta);
  loadProfiles(rid);
}

function showParams(meta) {
  const box = document.getElementById('paramsBox');
  if (!meta || !meta.corona) { box.style.display = 'none'; return; }
  const skip = new Set(['hash']);
  let html = '';
  for (const [k,v] of Object.entries(meta)) {
    if (skip.has(k)) continue;
    const fv = typeof v === 'number' ? (Math.abs(v)>=1e4||(Math.abs(v)<0.01&&v!==0) ? v.toExponential(2) : v) : v;
    html += `<span class="pk">${k}</span>=<span class="pv">${fv}</span> &nbsp; `;
  }
  box.innerHTML = html;
  box.style.display = 'block';
}

async function loadProfiles(rid) {
  document.getElementById('statusMsg').style.display = 'block';
  document.getElementById('statusMsg').textContent = 'Loading profiles...';
  document.getElementById('plotsContainer').style.display = 'none';

  const [rProf, rJ0] = await Promise.all([
    fetch(`/api/profiles/${rid}?step=2`),
    fetch(`/api/mean_intensity/${rid}?step=2`),
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
  if (hasJ0) plotJ0Iter(dJ0.moments);
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
      x: p.tau_mid, y: p.T_K,
      name: `iter ${it}`,
      mode: 'lines',
      line: { color: color, width: idx === n - 1 ? 2.5 : 1.2 },
      opacity: 0.4 + 0.6 * (n > 1 ? idx / (n - 1) : 1),
    });
  });

  Plotly.react(div, traces, {
    paper_bgcolor: '#0a0c10', plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:30, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:9}, orientation:'h', y:-0.15 },
    xaxis: { type:'log', title:'\u03c4 (Thomson)', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'T [K]', gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

function iterColor(idx, n) {
  const frac = n > 1 ? idx / (n - 1) : 1;
  const r = Math.round(0x38 + frac * (0xf0 - 0x38));
  const g = Math.round(0xbd + frac * (0xa0 - 0xbd));
  const b = Math.round(0xf8 + frac * (0x30 - 0xf8));
  return `rgb(${r},${g},${b})`;
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
      x: E_keV, y: EJ0,
      name: `iter ${it}`,
      mode: 'lines',
      line: { color: color, width: idx === n - 1 ? 2.5 : 1.2 },
      opacity: 0.4 + 0.6 * (n > 1 ? idx / (n - 1) : 1),
    });
  });

  // Auto y-range
  let peak = -Infinity;
  traces.forEach(t => t.y.forEach(v => { if (v > 0 && v > peak) peak = v; }));
  const yr = isFinite(peak) && peak > 0 ? [Math.log10(peak * 1e-6), Math.log10(peak * 5)] : undefined;

  Plotly.react(div, traces, {
    paper_bgcolor: '#0a0c10', plot_bgcolor: PLOT_BG,
    font: { family: 'Inter, Helvetica Neue, sans-serif', color: FONT_COLOR, size: 13 },
    margin: { l:70, r:30, t:30, b:60 },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: {size:9}, orientation:'h', y:-0.15 },
    xaxis: { type:'log', title:'E [keV]', range:[Math.log10(1e-3), Math.log10(1000)], gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR },
    yaxis: { type:'log', title:'E \u00d7 J\u2080 [erg cm\u207b\u00b2 s\u207b\u00b9 sr\u207b\u00b9]', range:yr, gridcolor:GRID_COLOR, zerolinecolor:GRID_COLOR, exponentformat:'e' },
  }, {responsive:true});
}

init();
</script>
</body>
</html>
"""

@app.route("/convergence")
def convergence():
    return render_template_string(CONVERGENCE_HTML)


# ─── Main ───────────────────────────────────────────────────────
def main():
    port = 5200
    url = f"http://127.0.0.1:{port}"
    print(f"\n  DAOv2.0 — X-ray Reflection Spectroscopy Model")
    print(f"  Author: Yimin Huang\n")
    print(f"  Opening browser → {url}\n")
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    app.run(host="127.0.0.1", port=port, debug=False, threaded=True)


if __name__ == "__main__":
    main()