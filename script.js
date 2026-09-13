/* Standalone browser build: React, Recharts, XLSX and Babel are loaded as
   global scripts by index.html (see the <script> tags there) instead of
   being imported as ES modules. */
const { useMemo, useRef, useState } = React;
const {
  ComposedChart, Line, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,
  BarChart, Bar, Cell,
} = Recharts;

/* ---------------------------------------------------------------------- */
/* Numeric helpers                                                        */
/* ---------------------------------------------------------------------- */

const toNum = (v) => {
  if (v === "" || v === null || v === undefined) return NaN;
  const n = parseFloat(v);
  return Number.isFinite(n) ? n : NaN;
};

const meanValid = (arr) => {
  const v = arr.map(toNum).filter((n) => !Number.isNaN(n));
  if (v.length === 0) return NaN;
  return v.reduce((a, b) => a + b, 0) / v.length;
};

const cvValid = (arr) => {
  const v = arr.map(toNum).filter((n) => !Number.isNaN(n));
  if (v.length < 2) return NaN;
  const m = v.reduce((a, b) => a + b, 0) / v.length;
  if (m === 0) return NaN;
  const sd = Math.sqrt(v.reduce((a, b) => a + (b - m) ** 2, 0) / (v.length - 1));
  return (sd / m) * 100;
};

const fmt = (n, d = 4) => (Number.isFinite(n) ? n.toFixed(d) : "—");
const fmtSig = (n, sig = 3) => {
  if (!Number.isFinite(n)) return "—";
  if (n === 0) return "0";
  const mag = Math.floor(Math.log10(Math.abs(n)));
  const factor = 10 ** (sig - 1 - mag);
  return (Math.round(n * factor) / factor).toString();
};

/* ---------------------------------------------------------------------- */
/* 4-parameter logistic fit: y = D + (A - D) / (1 + 10^((C - x) * B))     */
/* x = log10(concentration), y = blank-subtracted OD                      */
/* Fit via Nelder-Mead simplex (derivative-free, robust for 4 params)     */
/* ---------------------------------------------------------------------- */

function model4PL([A, B, C, D], x) {
  return D + (A - D) / (1 + Math.pow(10, (C - x) * B));
}

function sse(params, xs, ys) {
  let s = 0;
  for (let i = 0; i < xs.length; i++) {
    const e = model4PL(params, xs[i]) - ys[i];
    s += e * e;
  }
  return s;
}

function nelderMead(f, x0, opts = {}) {
  const n = x0.length;
  const alpha = 1, gamma = 2, rho = 0.5, sigma = 0.5;
  const maxIter = opts.maxIter || 4000;
  const tol = opts.tol || 1e-12;

  let simplex = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const p = x0.slice();
    p[i] += p[i] !== 0 ? p[i] * 0.1 : 0.1;
    simplex.push(p);
  }
  let fvals = simplex.map(f);

  for (let iter = 0; iter < maxIter; iter++) {
    const order = fvals.map((v, i) => i).sort((a, b) => fvals[a] - fvals[b]);
    simplex = order.map((i) => simplex[i]);
    fvals = order.map((i) => fvals[i]);

    if (Math.abs(fvals[n] - fvals[0]) < tol) break;

    const centroid = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) centroid[j] += simplex[i][j] / n;
    }

    const reflect = centroid.map((c, j) => c + alpha * (c - simplex[n][j]));
    const fReflect = f(reflect);

    if (fReflect < fvals[0]) {
      const expand = centroid.map((c, j) => c + gamma * (reflect[j] - c));
      const fExpand = f(expand);
      if (fExpand < fReflect) {
        simplex[n] = expand; fvals[n] = fExpand;
      } else {
        simplex[n] = reflect; fvals[n] = fReflect;
      }
    } else if (fReflect < fvals[n - 1]) {
      simplex[n] = reflect; fvals[n] = fReflect;
    } else {
      const contract = centroid.map((c, j) => c + rho * (simplex[n][j] - c));
      const fContract = f(contract);
      if (fContract < fvals[n]) {
        simplex[n] = contract; fvals[n] = fContract;
      } else {
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[0].map((c, j) => c + sigma * (simplex[i][j] - c));
          fvals[i] = f(simplex[i]);
        }
      }
    }
  }
  const order = fvals.map((v, i) => i).sort((a, b) => fvals[a] - fvals[b]);
  return { params: simplex[order[0]], value: fvals[order[0]] };
}

function fit4PL(xs, ys) {
  if (xs.length < 4) return null;
  const yMin = Math.min(...ys), yMax = Math.max(...ys);
  const xMed = xs.slice().sort((a, b) => a - b)[Math.floor(xs.length / 2)];

  const objective = (p) => sse(p, xs, ys);
  const seeds = [
    [yMax * 1.2 + 0.01, 1, xMed, Math.max(yMin - 0.01, 0)],
    [yMax * 1.5 + 0.01, 2, xMed, 0],
    [yMax * 1.2 + 0.01, -1, xMed, Math.max(yMin - 0.01, 0)],
  ];
  let best = null;
  for (const seed of seeds) {
    const r = nelderMead(objective, seed, { maxIter: 6000 });
    if (!best || r.value < best.value) best = r;
  }
  const [A, B, C, D] = best.params;

  const yPred = xs.map((x) => model4PL(best.params, x));
  const yBar = ys.reduce((a, b) => a + b, 0) / ys.length;
  const ssRes = ys.reduce((a, y, i) => a + (y - yPred[i]) ** 2, 0);
  const ssTot = ys.reduce((a, y) => a + (y - yBar) ** 2, 0);
  const r2 = ssTot === 0 ? 1 : 1 - ssRes / ssTot;

  return { A, B, C, D, r2, params: [A, B, C, D] };
}

function inverse4PL(fitRes, y) {
  const { A, B, C, D } = fitRes;
  const lo = Math.min(A, D), hi = Math.max(A, D);
  const eps = (hi - lo) * 1e-6;
  if (y <= lo + eps) return { status: "below", conc: null };
  if (y >= hi - eps) return { status: "above", conc: null };
  const ratio = (A - D) / (y - D) - 1;
  if (ratio <= 0) return { status: "below", conc: null };
  const x = C - Math.log10(ratio) / B;
  return { status: "ok", conc: Math.pow(10, x), x };
}

/* ---------------------------------------------------------------------- */
/* Default kit parameters — Vazyme EasyAna dsRNA (Modified) ELISA (2.0)   */
/* Manual DD3509EN, v24.2                                                  */
/* ---------------------------------------------------------------------- */

const DEFAULT_STD_CONCS = [1.5, 0.75, 0.375, 0.1875, 0.09375, 0.046875, 0.0234375];

const defaultKitParams = () => ({
  lodPg: 10,
  loqPg: 47,
  r2Threshold: 0.99,
  passPct: 1,
});

const emptyReps = () => ["", "", ""];

const defaultStandards = () =>
  DEFAULT_STD_CONCS.map((c) => ({ conc: c, od: emptyReps() }));

const newSample = (n) => ({
  id: `s${Date.now()}_${n}`,
  name: `Sample ${n}`,
  od: emptyReps(),
  dilution: 100,
  totalRna: "",
  expected: "",
});

/* ---------------------------------------------------------------------- */
/* Demo data — reproduces the worked example on p.6 of the kit insert     */
/* ---------------------------------------------------------------------- */

const DEMO_BLANK = ["0.050", "0.048", ""];
const DEMO_DIFFS = [3.3171, 1.6846, 0.8267, 0.4325, 0.2155, 0.1089, 0.0565];

function loadDemo(setBlankOD, setStandards, setSamples) {
  setBlankOD(DEMO_BLANK);
  const blankMean = meanValid(DEMO_BLANK);
  setStandards(
    DEFAULT_STD_CONCS.map((c, i) => ({
      conc: c,
      od: [(blankMean + DEMO_DIFFS[i]).toFixed(4), (blankMean + DEMO_DIFFS[i]).toFixed(4), ""],
    }))
  );
  setSamples([
    { id: "demo1", name: "IVT batch 12 (neat)", od: [(blankMean + 0.31).toFixed(4), (blankMean + 0.29).toFixed(4), ""], dilution: 1, totalRna: "480", expected: "" },
    { id: "demo2", name: "IVT batch 12 (1:50)", od: [(blankMean + 1.05).toFixed(4), (blankMean + 1.02).toFixed(4), ""], dilution: 50, totalRna: "480", expected: "" },
  ]);
}

/* ---------------------------------------------------------------------- */
/* Excel-style paste support                                              */
/* Lets a user copy a column (or block) of cells straight out of Excel/   */
/* Sheets and paste it into any grid cell here; the block is distributed  */
/* across rows/columns starting at the cell that received the paste,      */
/* growing the table if the pasted block runs past the last row.          */
/* ---------------------------------------------------------------------- */

function parsePastedGrid(text) {
  return text
    .replace(/\r/g, "")
    .split("\n")
    .filter((line, idx, arr) => !(idx === arr.length - 1 && line === ""))
    .map((line) => line.split("\t").map((c) => c.trim()));
}

// True only when the clipboard content looks like more than a single cell
// (i.e. it has row or column structure). A lone value falls through to the
// browser's normal single-cell paste behaviour.
function isMultiCellPaste(text) {
  return text.includes("\t") || text.includes("\n");
}

/* ---------------------------------------------------------------------- */
/* Chart-to-image capture (for PNG download and embedding into slides)    */
/* Recharts renders plain SVG, so we clone the live node, resolve the     */
/* CSS custom-property colors (var(--accent) etc. don't resolve once the  */
/* SVG is detached from this DOM), and rasterize via an offscreen canvas. */
/* ---------------------------------------------------------------------- */

const CSS_VAR_COLORS = {
  "var(--accent)": "#12645B",
  "var(--fail)": "#B3261E",
  "var(--warn)": "#9C6B12",
  "var(--ink)": "#15191A",
  "var(--ink-soft)": "#5B655F",
  "var(--border)": "#DBDFD7",
};

function resolveCssVarColors(svgString) {
  let s = svgString;
  Object.entries(CSS_VAR_COLORS).forEach(([k, v]) => { s = s.split(k).join(v); });
  return s;
}

function svgStringToPngDataUrl(svgString, width, height, scale = 2) {
  return new Promise((resolve, reject) => {
    const svg64 = btoa(unescape(encodeURIComponent(svgString)));
    const img = new Image();
    img.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = width * scale;
      canvas.height = height * scale;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.scale(scale, scale);
      ctx.drawImage(img, 0, 0, width, height);
      resolve(canvas.toDataURL("image/png"));
    };
    img.onerror = () => reject(new Error("Could not rasterize chart"));
    img.src = "data:image/svg+xml;base64," + svg64;
  });
}

async function captureChartPng(wrapperEl) {
  const svgEl = wrapperEl && wrapperEl.querySelector("svg");
  if (!svgEl) throw new Error("Chart isn't rendered yet — add sample data first.");
  const bbox = svgEl.getBoundingClientRect();
  const width = Math.max(Math.round(bbox.width), 1) || 700;
  const height = Math.max(Math.round(bbox.height), 1) || 300;
  const clone = svgEl.cloneNode(true);
  clone.setAttribute("width", width);
  clone.setAttribute("height", height);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const svgString = resolveCssVarColors(new XMLSerializer().serializeToString(clone));
  const dataUrl = await svgStringToPngDataUrl(svgString, width, height, 2);
  return { dataUrl, width, height };
}

function downloadDataUrl(dataUrl, filename) {
  const a = document.createElement("a");
  a.href = dataUrl; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}

/* ---------------------------------------------------------------------- */
/* Dependency-free .pptx writer                                           */
/* No CDN, no third-party library — builds a valid OOXML PowerPoint       */
/* package (a ZIP using the uncompressed STORE method) entirely by hand.  */
/* Each content slide is a title plus one full-bleed image (a chart, or a */
/* table rendered to PNG the same way charts are), so no OOXML table      */
/* markup is required.                                                    */
/* ---------------------------------------------------------------------- */

function crc32(bytes) {
  if (!crc32.table) {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
      t[n] = c;
    }
    crc32.table = t;
  }
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) crc = crc32.table[(crc ^ bytes[i]) & 0xFF] ^ (crc >>> 8);
  return (crc ^ 0xFFFFFFFF) >>> 0;
}

// Minimal ZIP writer (STORE method — no compression needed, so no zlib
// dependency). DOS timestamp is fixed at 1980-01-01; PowerPoint ignores it.
function zipStore(entries) {
  const enc = new TextEncoder();
  const chunks = [];
  const centrals = [];
  let offset = 0;
  for (const { name, data } of entries) {
    const nameBytes = enc.encode(name);
    const crc = crc32(data);
    const size = data.length;

    const local = new Uint8Array(30 + nameBytes.length);
    const lv = new DataView(local.buffer);
    lv.setUint32(0, 0x04034b50, true);
    lv.setUint16(4, 20, true);
    lv.setUint16(12, 0x21, true);
    lv.setUint32(14, crc, true);
    lv.setUint32(18, size, true);
    lv.setUint32(22, size, true);
    lv.setUint16(26, nameBytes.length, true);
    local.set(nameBytes, 30);
    chunks.push(local, data);

    const central = new Uint8Array(46 + nameBytes.length);
    const cv = new DataView(central.buffer);
    cv.setUint32(0, 0x02014b50, true);
    cv.setUint16(4, 20, true);
    cv.setUint16(6, 20, true);
    cv.setUint16(14, 0x21, true);
    cv.setUint32(16, crc, true);
    cv.setUint32(20, size, true);
    cv.setUint32(24, size, true);
    cv.setUint16(28, nameBytes.length, true);
    cv.setUint32(42, offset, true);
    central.set(nameBytes, 46);
    centrals.push(central);

    offset += local.length + data.length;
  }
  const centralStart = offset;
  let centralSize = 0;
  for (const c of centrals) { chunks.push(c); centralSize += c.length; }

  const eocd = new Uint8Array(22);
  const ev = new DataView(eocd.buffer);
  ev.setUint32(0, 0x06054b50, true);
  ev.setUint16(8, entries.length, true);
  ev.setUint16(10, entries.length, true);
  ev.setUint32(12, centralSize, true);
  ev.setUint32(16, centralStart, true);
  chunks.push(eocd);

  return new Blob(chunks, { type: "application/vnd.openxmlformats-officedocument.presentationml.presentation" });
}

function xmlEscape(s) {
  return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;").replace(/'/g, "&apos;");
}

function dataUrlToBytes(dataUrl) {
  const binary = atob(dataUrl.split(",")[1]);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return bytes;
}

// Renders a title + data table to a PNG via the same SVG/canvas pipeline
// used for chart export, so the summary table needs no extra dependency.
async function tableToPng(title, headers, rows, width = 1180) {
  const rowH = 30, headH = 40, pad = 56;
  const height = pad + headH + rows.length * rowH + 24;
  const th = headers.map((h) => `<th style="text-align:left;padding:6px 10px;font-size:13px;color:#5B655F;border-bottom:1px solid #DBDFD7;font-family:-apple-system,Arial,sans-serif;">${xmlEscape(h)}</th>`).join("");
  const trs = rows.map((r, i) => `<tr style="background:${i % 2 ? "#F5F6F3" : "#FFFFFF"}">${r.map((c) => `<td style="padding:6px 10px;font-size:13px;font-family:ui-monospace,Menlo,monospace;color:#15191A;border-bottom:1px solid #EEF0EC;white-space:nowrap;">${xmlEscape(c)}</td>`).join("")}</tr>`).join("");
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}"><foreignObject width="100%" height="100%"><div xmlns="http://www.w3.org/1999/xhtml" style="font-family:-apple-system,Arial,sans-serif;background:#fff;">
    <div style="font-size:20px;font-weight:700;color:#15191A;padding:10px 4px 16px 4px;">${xmlEscape(title)}</div>
    <table style="border-collapse:collapse;width:100%;"><thead><tr>${th}</tr></thead><tbody>${trs}</tbody></table>
  </div></foreignObject></svg>`;
  const dataUrl = await svgStringToPngDataUrl(svg, width, height, 2);
  return { dataUrl, width, height };
}

const EMU_PER_PX = 9525; // at 96 DPI
const SLIDE_W = Math.round(13.333 * 914400);
const SLIDE_H = Math.round(7.5 * 914400);

function fitImageBox(imgWpx, imgHpx, top, bottomMargin = 400000, sideMargin = 609600) {
  const imgW = imgWpx * EMU_PER_PX, imgH = imgHpx * EMU_PER_PX;
  const maxW = SLIDE_W - 2 * sideMargin;
  const maxH = SLIDE_H - top - bottomMargin;
  const scale = Math.min(maxW / imgW, maxH / imgH, 1);
  const w = Math.round(imgW * scale), h = Math.round(imgH * scale);
  return { x: Math.round(sideMargin + (maxW - w) / 2), y: Math.round(top + (maxH - h) / 2), cx: w, cy: h };
}

const OOXML_NS = 'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"';

function titleSlideXml(title, subtitle) {
  return `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sld ${OOXML_NS}><p:cSld><p:spTree>
<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/>
<p:sp><p:nvSpPr><p:cNvPr id="2" name="Title"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>
<p:spPr><a:xfrm><a:off x="609600" y="2286000"/><a:ext cx="10972800" cy="1097280"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr>
<p:txBody><a:bodyPr/><a:lstStyle/><a:p><a:r><a:rPr lang="en-US" sz="3600" b="1"><a:solidFill><a:srgbClr val="15191A"/></a:solidFill></a:rPr><a:t>${xmlEscape(title)}</a:t></a:r></a:p></p:txBody>
</p:sp>
<p:sp><p:nvSpPr><p:cNvPr id="3" name="Subtitle"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>
<p:spPr><a:xfrm><a:off x="609600" y="3474720"/><a:ext cx="10972800" cy="609600"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr>
<p:txBody><a:bodyPr/><a:lstStyle/><a:p><a:r><a:rPr lang="en-US" sz="1400"><a:solidFill><a:srgbClr val="5B655F"/></a:solidFill></a:rPr><a:t>${xmlEscape(subtitle)}</a:t></a:r></a:p></p:txBody>
</p:sp>
</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>`;
}

function imageSlideXml(title, box) {
  return `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sld ${OOXML_NS}><p:cSld><p:spTree>
<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/>
<p:sp><p:nvSpPr><p:cNvPr id="2" name="Title"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>
<p:spPr><a:xfrm><a:off x="457200" y="228600"/><a:ext cx="11277600" cy="640080"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr>
<p:txBody><a:bodyPr/><a:lstStyle/><a:p><a:r><a:rPr lang="en-US" sz="2000" b="1"><a:solidFill><a:srgbClr val="15191A"/></a:solidFill></a:rPr><a:t>${xmlEscape(title)}</a:t></a:r></a:p></p:txBody>
</p:sp>
<p:pic><p:nvPicPr><p:cNvPr id="4" name="Picture"/><p:cNvPicPr/><p:nvPr/></p:nvPicPr>
<p:blipFill><a:blip r:embed="rIdImg"/><a:stretch><a:fillRect/></a:stretch></p:blipFill>
<p:spPr><a:xfrm><a:off x="${box.x}" y="${box.y}"/><a:ext cx="${box.cx}" cy="${box.cy}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr>
</p:pic>
</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>`;
}

const SLIDE_RELS_NO_IMG = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rIdLayout" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>
</Relationships>`;

function slideRelsWithImg(imgFile) {
  return `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rIdLayout" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>
<Relationship Id="rIdImg" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/${imgFile}"/>
</Relationships>`;
}

const THEME_XML = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="dsRNA">
<a:themeElements>
<a:clrScheme name="dsRNA">
<a:dk1><a:sysClr val="windowText" lastClr="000000"/></a:dk1>
<a:lt1><a:sysClr val="window" lastClr="FFFFFF"/></a:lt1>
<a:dk2><a:srgbClr val="15191A"/></a:dk2>
<a:lt2><a:srgbClr val="EEF0EC"/></a:lt2>
<a:accent1><a:srgbClr val="12645B"/></a:accent1>
<a:accent2><a:srgbClr val="B3261E"/></a:accent2>
<a:accent3><a:srgbClr val="9C6B12"/></a:accent3>
<a:accent4><a:srgbClr val="5B655F"/></a:accent4>
<a:accent5><a:srgbClr val="DBDFD7"/></a:accent5>
<a:accent6><a:srgbClr val="E4EFEC"/></a:accent6>
<a:hlink><a:srgbClr val="0563C1"/></a:hlink>
<a:folHlink><a:srgbClr val="954F72"/></a:folHlink>
</a:clrScheme>
<a:fontScheme name="dsRNA">
<a:majorFont><a:latin typeface="Calibri"/><a:ea typeface=""/><a:cs typeface=""/></a:majorFont>
<a:minorFont><a:latin typeface="Calibri"/><a:ea typeface=""/><a:cs typeface=""/></a:minorFont>
</a:fontScheme>
<a:fmtScheme name="dsRNA">
<a:fillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:fillStyleLst>
<a:lnStyleLst><a:ln w="6350"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln><a:ln w="12700"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln><a:ln w="19050"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln></a:lnStyleLst>
<a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst>
<a:bgFillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:bgFillStyleLst>
</a:fmtScheme>
</a:themeElements>
</a:theme>`;

const SLIDE_MASTER_XML = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldMaster ${OOXML_NS}><p:cSld><p:bg><p:bgRef idx="1001"><a:schemeClr val="bg1"/></p:bgRef></p:bg><p:spTree>
<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/>
</p:spTree></p:cSld>
<p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/>
<p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst>
</p:sldMaster>`;

const SLIDE_MASTER_RELS_XML = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme" Target="../theme/theme1.xml"/>
</Relationships>`;

const SLIDE_LAYOUT_XML = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldLayout ${OOXML_NS} type="blank" preserve="1"><p:cSld><p:spTree>
<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/>
</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>`;

const SLIDE_LAYOUT_RELS_XML = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="../slideMasters/slideMaster1.xml"/>
</Relationships>`;

const CORE_XML = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<dc:title>dsRNA ELISA Results</dc:title><dc:creator>dsRNA ELISA Analyzer</dc:creator>
</cp:coreProperties>`;

const APP_XML = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"><Application>dsRNA ELISA Analyzer</Application></Properties>`;

// slides: [{ title, subtitle }] for the title slide, or
//         [{ title, imageDataUrl, imgW, imgH, top? }] for an image slide.
async function buildPptx(slides) {
  const enc = new TextEncoder();
  const entries = [];
  const add = (name, text) => entries.push({ name, data: enc.encode(text) });

  const slideParts = [];
  slides.forEach((s, i) => {
    if (s.imageDataUrl) {
      const box = fitImageBox(s.imgW, s.imgH, s.top || 1100000);
      const imgFile = `image${i + 1}.png`;
      slideParts.push({
        xml: imageSlideXml(s.title || "", box),
        rels: slideRelsWithImg(imgFile),
        imgEntry: { name: `ppt/media/${imgFile}`, data: dataUrlToBytes(s.imageDataUrl) },
      });
    } else {
      slideParts.push({ xml: titleSlideXml(s.title, s.subtitle || ""), rels: SLIDE_RELS_NO_IMG });
    }
  });

  add("[Content_Types].xml", `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Default Extension="png" ContentType="image/png"/>
<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>
<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>
<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>
<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>
${slideParts.map((_, i) => `<Override PartName="/ppt/slides/slide${i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>`).join("\n")}
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>`);

  add("_rels/.rels", `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="ppt/presentation.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>`);

  add("docProps/core.xml", CORE_XML);
  add("docProps/app.xml", APP_XML);

  add("ppt/presentation.xml", `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:presentation ${OOXML_NS}>
<p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rIdM"/></p:sldMasterIdLst>
<p:sldIdLst>${slideParts.map((_, i) => `<p:sldId id="${256 + i}" r:id="rIdS${i + 1}"/>`).join("")}</p:sldIdLst>
<p:sldSz cx="${SLIDE_W}" cy="${SLIDE_H}"/>
<p:notesSz cx="${SLIDE_H}" cy="${SLIDE_W}"/>
</p:presentation>`);

  add("ppt/_rels/presentation.xml.rels", `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rIdM" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="slideMasters/slideMaster1.xml"/>
${slideParts.map((_, i) => `<Relationship Id="rIdS${i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide" Target="slides/slide${i + 1}.xml"/>`).join("\n")}
</Relationships>`);

  add("ppt/slideMasters/slideMaster1.xml", SLIDE_MASTER_XML);
  add("ppt/slideMasters/_rels/slideMaster1.xml.rels", SLIDE_MASTER_RELS_XML);
  add("ppt/slideLayouts/slideLayout1.xml", SLIDE_LAYOUT_XML);
  add("ppt/slideLayouts/_rels/slideLayout1.xml.rels", SLIDE_LAYOUT_RELS_XML);
  add("ppt/theme/theme1.xml", THEME_XML);

  slideParts.forEach((part, i) => {
    add(`ppt/slides/slide${i + 1}.xml`, part.xml);
    add(`ppt/slides/_rels/slide${i + 1}.xml.rels`, part.rels);
    if (part.imgEntry) entries.push(part.imgEntry);
  });

  return zipStore(entries);
}

function downloadCsv(rows, filename) {
  const csv = rows.map((r) => r.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

/* ---------------------------------------------------------------------- */
/* UI atoms                                                                */
/* ---------------------------------------------------------------------- */

const Step = ({ n, title, sub, children }) => (
  <section style={{ display: "flex", gap: 20, marginBottom: 36 }}>
    <div style={{ flexShrink: 0, width: 30, textAlign: "right" }}>
      <div style={{
        fontFamily: "var(--mono)", fontSize: 13, color: "var(--accent)",
        borderRight: "2px solid var(--border)", paddingRight: 14, paddingTop: 3,
      }}>{n}</div>
    </div>
    <div style={{ flex: 1, minWidth: 0 }}>
      <h2 style={{ margin: "0 0 2px 0", fontSize: 17, fontWeight: 600, color: "var(--ink)" }}>{title}</h2>
      {sub && <p style={{ margin: "0 0 14px 0", fontSize: 13.5, color: "var(--ink-soft)", maxWidth: 640 }}>{sub}</p>}
      {children}
    </div>
  </section>
);

const Field = ({ label, children, width }) => (
  <label style={{ display: "block", fontSize: 12, color: "var(--ink-soft)", width }}>
    <span style={{ display: "block", marginBottom: 4 }}>{label}</span>
    {children}
  </label>
);

const inputStyle = {
  width: "100%", padding: "6px 8px", fontFamily: "var(--mono)", fontSize: 13,
  border: "1px solid var(--border)", borderRadius: 4, background: "#fff", color: "var(--ink)",
  outline: "none",
};

const NumIn = (props) => (
  <input
    type="text"
    inputMode="decimal"
    {...props}
    style={{ ...inputStyle, ...(props.style || {}) }}
    onFocus={(e) => (e.target.style.borderColor = "var(--accent)")}
    onBlur={(e) => { e.target.style.borderColor = "var(--border)"; props.onBlur && props.onBlur(e); }}
  />
);

const Th = ({ children, align = "left" }) => (
  <th style={{
    textAlign: align, fontSize: 11, fontWeight: 600, color: "var(--ink-soft)",
    textTransform: "none", padding: "0 8px 8px 8px", borderBottom: "1px solid var(--border)",
    whiteSpace: "nowrap",
  }}>{children}</th>
);

const Td = ({ children, align = "left", mono = true, style }) => (
  <td style={{
    padding: "6px 8px", fontSize: 13, textAlign: align, verticalAlign: "middle",
    fontFamily: mono ? "var(--mono)" : "inherit", color: "var(--ink)", ...style,
  }}>{children}</td>
);

const Badge = ({ tone, children }) => {
  const tones = {
    pass: { bg: "var(--accent-soft)", fg: "var(--accent)" },
    warn: { bg: "#FBF0DD", fg: "var(--warn)" },
    fail: { bg: "#F7E1DE", fg: "var(--fail)" },
    neutral: { bg: "#EEF0EC", fg: "var(--ink-soft)" },
  };
  const t = tones[tone] || tones.neutral;
  return (
    <span style={{
      display: "inline-block", padding: "2px 8px", borderRadius: 999, fontSize: 11.5,
      fontFamily: "var(--mono)", background: t.bg, color: t.fg, fontWeight: 600, whiteSpace: "nowrap",
    }}>{children}</span>
  );
};

const Btn = ({ children, onClick, variant = "default", ...rest }) => {
  const styles = {
    default: { background: "#fff", color: "var(--ink)", border: "1px solid var(--border)" },
    primary: { background: "var(--accent)", color: "#fff", border: "1px solid var(--accent)" },
    ghost: { background: "transparent", color: "var(--ink-soft)", border: "1px solid transparent" },
  };
  return (
    <button
      onClick={onClick}
      {...rest}
      style={{
        ...styles[variant], padding: "7px 14px", borderRadius: 6, fontSize: 13, fontWeight: 500,
        cursor: "pointer", fontFamily: "var(--sans)",
      }}
    >{children}</button>
  );
};

/* ---------------------------------------------------------------------- */
/* Main component                                                         */
/* ---------------------------------------------------------------------- */

function DsRnaElisaAnalyzer() {
  const [kit, setKit] = useState(defaultKitParams());
  const [showKit, setShowKit] = useState(false);
  const [blankOD, setBlankOD] = useState(["", "", ""]);
  const [standards, setStandards] = useState(defaultStandards());
  const [samples, setSamples] = useState([newSample(1)]);
  const concChartRef = useRef(null);
  const pctChartRef = useRef(null);
  const curveChartRef = useRef(null);

  const blankMean = useMemo(() => meanValid(blankOD), [blankOD]);

  const stdRows = useMemo(() => {
    return standards.map((s) => {
      const rawMean = meanValid(s.od);
      const sub = Number.isFinite(rawMean) && Number.isFinite(blankMean) ? rawMean - blankMean : NaN;
      const cv = cvValid(s.od);
      return { ...s, rawMean, sub, cv, logConc: s.conc > 0 ? Math.log10(s.conc) : NaN };
    });
  }, [standards, blankMean]);

  const fit = useMemo(() => {
    const pts = stdRows.filter((r) => Number.isFinite(r.logConc) && Number.isFinite(r.sub));
    if (pts.length < 4) return null;
    return fit4PL(pts.map((p) => p.logConc), pts.map((p) => p.sub));
  }, [stdRows]);

  const r2Pass = fit && fit.r2 >= kit.r2Threshold;

  const curveData = useMemo(() => {
    if (!fit) return [];
    const xs = stdRows.map((r) => r.logConc).filter(Number.isFinite);
    if (xs.length === 0) return [];
    const lo = Math.min(...xs), hi = Math.max(...xs);
    const pad = (hi - lo) * 0.08;
    const pts = [];
    const N = 80;
    for (let i = 0; i <= N; i++) {
      const x = lo - pad + ((hi - lo + 2 * pad) * i) / N;
      pts.push({ x, y: model4PL(fit.params, x) });
    }
    return pts;
  }, [fit, stdRows]);

  const standardPoints = useMemo(
    () => stdRows.filter((r) => Number.isFinite(r.logConc) && Number.isFinite(r.sub)).map((r) => ({ x: r.logConc, y: r.sub })),
    [stdRows]
  );

  const sampleResults = useMemo(() => {
    return samples.map((s) => {
      const rawMean = meanValid(s.od);
      const cv = cvValid(s.od);
      const sub = Number.isFinite(rawMean) && Number.isFinite(blankMean) ? rawMean - blankMean : NaN;
      const dilution = toNum(s.dilution) || 100;

      if (!fit || !Number.isFinite(sub)) {
        return { ...s, rawMean, sub, cv, dilution, status: "no-fit", finalConc: null, pct: null, recovery: null };
      }

      const inv = inverse4PL(fit, sub);
      let finalConc = null, belowLoq = false, pct = null, recovery = null;

      if (inv.status === "ok") {
        const wellConc = inv.conc; // ng/mL, as read in the well (post any pre-dilution)
        const loqNgml = kit.loqPg / 1000;
        belowLoq = wellConc < loqNgml;
        finalConc = wellConc * dilution; // ng/mL in original sample

        const totalRna = toNum(s.totalRna); // ng/uL
        if (Number.isFinite(totalRna) && totalRna > 0) {
          const dsRnaNgUl = finalConc / 1000;
          pct = (dsRnaNgUl / totalRna) * 100;
        }
        const expected = toNum(s.expected);
        if (Number.isFinite(expected) && expected > 0) {
          recovery = (finalConc / expected) * 100;
        }
      }

      return {
        ...s, rawMean, sub, cv, dilution,
        status: inv.status, x: inv.x, finalConc, belowLoq, pct, recovery,
      };
    });
  }, [samples, fit, blankMean, kit]);

  const passThreshold = toNum(kit.passPct);

  const passLabel = (r) => {
    if (r.pct === null || !Number.isFinite(r.pct)) return "n/a";
    if (!Number.isFinite(passThreshold)) return "n/a";
    return r.pct <= passThreshold ? "Pass" : "Fail";
  };

  const passBadge = (r) => {
    const label = passLabel(r);
    if (label === "n/a") return <Badge tone="neutral">n/a</Badge>;
    return <Badge tone={label === "Pass" ? "pass" : "fail"}>{label}</Badge>;
  };

  const comparisonConcData = useMemo(
    () =>
      sampleResults
        .filter((r) => r.status === "ok" && Number.isFinite(r.finalConc))
        .map((r) => ({ name: r.name, conc: r.finalConc })),
    [sampleResults]
  );

  const comparisonPctData = useMemo(
    () =>
      sampleResults
        .filter((r) => r.pct !== null && Number.isFinite(r.pct))
        .map((r) => ({
          name: r.name,
          pct: r.pct,
          pass: Number.isFinite(passThreshold) ? r.pct <= passThreshold : null,
        })),
    [sampleResults, passThreshold]
  );

  const samplePoints = useMemo(
    () => sampleResults.filter((r) => r.status === "ok" && Number.isFinite(r.x) && Number.isFinite(r.sub)).map((r) => ({ x: r.x, y: r.sub })),
    [sampleResults]
  );

  const updateStdOD = (i, repIdx, val) => {
    setStandards((prev) => {
      const next = prev.slice();
      const od = next[i].od.slice();
      od[repIdx] = val;
      next[i] = { ...next[i], od };
      return next;
    });
  };
  const updateStdConc = (i, val) => {
    setStandards((prev) => {
      const next = prev.slice();
      next[i] = { ...next[i], conc: toNum(val) };
      return next;
    });
  };
  const updateSample = (id, patch) => {
    setSamples((prev) => prev.map((s) => (s.id === id ? { ...s, ...patch } : s)));
  };
  const updateSampleOD = (id, repIdx, val) => {
    setSamples((prev) =>
      prev.map((s) => {
        if (s.id !== id) return s;
        const od = s.od.slice(); od[repIdx] = val;
        return { ...s, od };
      })
    );
  };
  const addSample = () => setSamples((prev) => [...prev, newSample(prev.length + 1)]);
  const removeSample = (id) => setSamples((prev) => (prev.length > 1 ? prev.filter((s) => s.id !== id) : prev));
  const clearAllSamples = () => setSamples([newSample(1)]);
  const clearBlank = () => setBlankOD(["", "", ""]);
  const clearStandards = () => setStandards(defaultStandards());

  // Blank-well paste: single conceptual column (Rep1/2/3), but accepts a
  // pasted row or column of up to 3 values from Excel.
  const handleBlankPaste = (e, colIndex) => {
    const text = e.clipboardData.getData("text");
    if (!isMultiCellPaste(text)) return;
    e.preventDefault();
    const flat = parsePastedGrid(text).flat();
    setBlankOD((prev) => {
      const next = prev.slice();
      flat.forEach((val, i) => {
        const idx = colIndex + i;
        while (next.length <= idx) next.push("");
        next[idx] = val;
      });
      return next.slice(0, Math.max(next.length, 3));
    });
  };

  // Standards grid paste. Columns: 0 = Conc, 1-3 = Rep1-3 OD.
  const handleStdPaste = (e, rowIndex, colIndex) => {
    const text = e.clipboardData.getData("text");
    if (!isMultiCellPaste(text)) return;
    e.preventDefault();
    const grid = parsePastedGrid(text);
    setStandards((prev) => {
      const next = prev.map((s) => ({ ...s, od: s.od.slice() }));
      grid.forEach((row, r) => {
        const targetRow = rowIndex + r;
        while (next.length <= targetRow) next.push({ conc: "", od: emptyReps() });
        row.forEach((val, c) => {
          const targetCol = colIndex + c;
          if (targetCol === 0) next[targetRow].conc = toNum(val);
          else if (targetCol >= 1 && targetCol <= 3) next[targetRow].od[targetCol - 1] = val;
        });
      });
      return next;
    });
  };

  // Samples grid paste. Columns: 0 Name, 1-3 Rep1-3 OD, 4 Dilution, 5 Total RNA, 6 Expected.
  const handleSamplePaste = (e, rowIndex, colIndex) => {
    const text = e.clipboardData.getData("text");
    if (!isMultiCellPaste(text)) return;
    e.preventDefault();
    const grid = parsePastedGrid(text);
    setSamples((prev) => {
      const next = prev.map((s) => ({ ...s, od: s.od.slice() }));
      grid.forEach((row, r) => {
        const targetRow = rowIndex + r;
        while (next.length <= targetRow) next.push(newSample(next.length + 1));
        row.forEach((val, c) => {
          const targetCol = colIndex + c;
          switch (targetCol) {
            case 0: next[targetRow].name = val; break;
            case 1: next[targetRow].od[0] = val; break;
            case 2: next[targetRow].od[1] = val; break;
            case 3: next[targetRow].od[2] = val; break;
            case 4: next[targetRow].dilution = val; break;
            case 5: next[targetRow].totalRna = val; break;
            case 6: next[targetRow].expected = val; break;
            default: break;
          }
        });
      });
      return next;
    });
  };

  const resultsTableHeader = [
    "Sample", "Mean raw OD", "Blank-subtracted OD", "Replicate CV %", "Dilution factor", "Status",
    "dsRNA concentration (ng/mL)", "Below LOQ", "Total RNA (ng/uL)", "% dsRNA (w/w)",
    `Pass threshold (% dsRNA ≤ ${fmt(passThreshold, 2)})`, "Pass/Fail", "Recovery %",
  ];
  const resultsTableRows = () =>
    sampleResults.map((r) => [
      r.name, fmt(r.rawMean, 4), fmt(r.sub, 4), fmt(r.cv, 1), r.dilution, r.status,
      r.finalConc !== null ? fmtSig(r.finalConc, 4) : "", r.belowLoq ? "yes" : "no",
      r.totalRna || "", r.pct !== null ? fmt(r.pct, 3) : "",
      fmt(passThreshold, 2), passLabel(r),
      r.recovery !== null ? fmt(r.recovery, 0) : "",
    ]);

  const exportCsv = () => {
    downloadCsv([resultsTableHeader, ...resultsTableRows()], "dsrna_elisa_results.csv");
  };

  const exportExcel = () => {
    const ws = XLSX.utils.aoa_to_sheet([resultsTableHeader, ...resultsTableRows()]);
    ws["!cols"] = resultsTableHeader.map(() => ({ wch: 20 }));
    const wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, ws, "Results");
    XLSX.writeFile(wb, "dsrna_elisa_results.xlsx");
  };

  // Downloads the Excel workbook plus a PNG of each comparison chart, so you
  // get the data and the figures together as a small set of files. (A true
  // embedded chart-in-cell requires a paid Excel library — this is the
  // closest equivalent without one.)
  const exportExcelWithCharts = async () => {
    exportExcel();
    try {
      if (fit && curveChartRef.current) {
        const { dataUrl } = await captureChartPng(curveChartRef.current);
        setTimeout(() => downloadDataUrl(dataUrl, "standard_curve.png"), 300);
      }
      if (comparisonConcData.length > 0 && concChartRef.current) {
        const { dataUrl } = await captureChartPng(concChartRef.current);
        setTimeout(() => downloadDataUrl(dataUrl, "dsrna_concentration_by_sample.png"), 700);
      }
      if (comparisonPctData.length > 0 && pctChartRef.current) {
        const { dataUrl } = await captureChartPng(pctChartRef.current);
        setTimeout(() => downloadDataUrl(dataUrl, "pct_dsrna_by_sample.png"), 1100);
      }
    } catch (err) {
      // Excel already downloaded; charts just aren't ready yet — not fatal.
      console.warn(err);
    }
  };

  const downloadChartPng = async (ref, filename) => {
    try {
      const { dataUrl } = await captureChartPng(ref.current);
      downloadDataUrl(dataUrl, filename);
    } catch (err) {
      alert(err.message);
    }
  };

  const exportSlides = async () => {
    try {
      const slides = [
        { title: "dsRNA ELISA — Results", subtitle: new Date().toLocaleDateString() },
      ];

      if (fit && curveChartRef.current) {
        const { dataUrl, width, height } = await captureChartPng(curveChartRef.current);
        slides.push({ title: `Standard curve (4PL, R² = ${fmt(fit.r2, 4)})`, imageDataUrl: dataUrl, imgW: width, imgH: height });
      }
      if (comparisonConcData.length > 0 && concChartRef.current) {
        const { dataUrl, width, height } = await captureChartPng(concChartRef.current);
        slides.push({ title: "dsRNA concentration by sample", imageDataUrl: dataUrl, imgW: width, imgH: height });
      }
      if (comparisonPctData.length > 0 && pctChartRef.current) {
        const { dataUrl, width, height } = await captureChartPng(pctChartRef.current);
        slides.push({ title: `% dsRNA (w/w) by sample — pass threshold ≤ ${fmt(passThreshold, 2)}%`, imageDataUrl: dataUrl, imgW: width, imgH: height });
      }

      const header = ["Sample", "dsRNA conc.", "% dsRNA", "Pass/Fail", "Recovery"];
      const rows = sampleResults.map((r) => [
        r.name,
        r.status === "ok" ? (r.belowLoq ? "< LOQ" : `${fmtSig(r.finalConc, 4)} ng/mL`) : "—",
        r.pct !== null ? `${fmt(r.pct, 3)}%` : "—",
        passLabel(r),
        r.recovery !== null ? `${fmt(r.recovery, 0)}%` : "—",
      ]);
      const table = await tableToPng("Summary", header, rows);
      slides.push({ imageDataUrl: table.dataUrl, imgW: table.width, imgH: table.height, top: 400000 });

      const blob = await buildPptx(slides);
      const url = URL.createObjectURL(blob);
      downloadDataUrl(url, "dsrna_elisa_summary.pptx");
      setTimeout(() => URL.revokeObjectURL(url), 2000);
    } catch (err) {
      alert("Could not export slides: " + err.message);
    }
  };

  const statusBadge = (r) => {
    if (r.status === "no-fit") return <Badge tone="neutral">enter data</Badge>;
    if (r.status === "above") return <Badge tone="fail">above curve — dilute &amp; retest</Badge>;
    if (r.status === "below") return <Badge tone="neutral">below curve</Badge>;
    if (r.belowLoq) return <Badge tone="warn">below LOQ</Badge>;
    return <Badge tone="pass">in range</Badge>;
  };

  return (
    <div style={{
      "--bg": "#F5F6F3", "--panel": "#FFFFFF", "--border": "#DBDFD7", "--ink": "#15191A",
      "--ink-soft": "#5B655F", "--accent": "#12645B", "--accent-soft": "#E4EFEC",
      "--warn": "#9C6B12", "--fail": "#B3261E",
      "--mono": "ui-monospace, SFMono-Regular, 'IBM Plex Mono', Menlo, Consolas, monospace",
      "--sans": "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif",
      background: "var(--bg)", color: "var(--ink)", fontFamily: "var(--sans)",
      minHeight: "100%", padding: "28px 24px 60px 24px",
    }}>
      <div style={{ maxWidth: 900, margin: "0 auto" }}>

        {/* Header */}
        <div style={{ marginBottom: 30, borderBottom: "1px solid var(--border)", paddingBottom: 18 }}>
          <div style={{ display: "flex", justifyContent: "space-between", alignItems: "flex-start", gap: 16, flexWrap: "wrap" }}>
            <div>
              <h1 style={{ margin: 0, fontSize: 22, fontWeight: 700, letterSpacing: "-0.01em" }}>dsRNA ELISA Analyzer</h1>
              <p style={{ margin: "4px 0 0 0", fontSize: 13.5, color: "var(--ink-soft)" }}>
                Sandwich ELISA quantitation of residual dsRNA in IVT mRNA — calibrated to the Vazyme EasyAna
                dsRNA (Modified) kit insert, DD3509EN v24.2.
              </p>
            </div>
            <div style={{ display: "flex", gap: 8 }}>
              <Btn onClick={() => loadDemo(setBlankOD, setStandards, setSamples)}>Load worked example</Btn>
              <Btn variant="ghost" onClick={() => setShowKit((v) => !v)}>{showKit ? "Hide" : "Kit parameters"}</Btn>
            </div>
          </div>

          {showKit && (
            <div style={{
              marginTop: 16, padding: 14, background: "var(--panel)", border: "1px solid var(--border)",
              borderRadius: 8, display: "flex", gap: 24, flexWrap: "wrap",
            }}>
              <Field label="LOD (pg/mL)" width={130}>
                <NumIn value={kit.lodPg} onChange={(e) => setKit((k) => ({ ...k, lodPg: e.target.value }))} />
              </Field>
              <Field label="LOQ (pg/mL)" width={130}>
                <NumIn value={kit.loqPg} onChange={(e) => setKit((k) => ({ ...k, loqPg: e.target.value }))} />
              </Field>
              <Field label="Min. acceptable R²" width={150}>
                <NumIn value={kit.r2Threshold} onChange={(e) => setKit((k) => ({ ...k, r2Threshold: e.target.value }))} />
              </Field>
              <Field label="Pass threshold, % dsRNA (w/w)" width={190}>
                <NumIn value={kit.passPct} onChange={(e) => setKit((k) => ({ ...k, passPct: e.target.value }))} />
              </Field>
              <div style={{ fontSize: 12, color: "var(--ink-soft)", maxWidth: 280, paddingTop: 18 }}>
                Defaults are the kit's stated sensitivity and QC spec (§09, §11-1). The pass threshold is not
                part of the kit spec — set it from your own process capability and safety/tox qualification data.
                Edit if a different lot, kit version, or lab-validated criterion applies.
              </div>
            </div>
          )}
        </div>

        {/* Step 1: Blank */}
        <Step n="01" title="Control well (blank)"
          sub="One diluent-only well per plate. Its OD is subtracted from every standard and sample well before curve fitting. Tip: copy a column of values from Excel and paste directly into any field below — it will fill across the remaining cells automatically.">
          <div style={{ display: "flex", gap: 10 }}>
            {blankOD.map((v, i) => (
              <Field key={i} label={`Rep ${i + 1}`} width={100}>
                <NumIn placeholder="OD450" value={v} onChange={(e) => {
                  const arr = blankOD.slice(); arr[i] = e.target.value; setBlankOD(arr);
                }} onPaste={(e) => handleBlankPaste(e, i)} />
              </Field>
            ))}
            <div style={{ paddingTop: 18, fontSize: 13, color: "var(--ink-soft)" }}>
              Mean blank OD: <span style={{ fontFamily: "var(--mono)", color: "var(--ink)" }}>{fmt(blankMean)}</span>
            </div>
            <div style={{ paddingTop: 14, marginLeft: "auto" }}>
              <Btn variant="ghost" onClick={clearBlank}>Clear</Btn>
            </div>
          </div>
        </Step>

        {/* Step 2: Standards */}
        <Step n="02" title="Standard curve"
          sub="Enter the OD450 read for each standard well. Concentrations default to the kit's serial-dilution ladder — edit if you prepared a different series. You can paste a whole column (or block) copied from Excel straight into any cell.">
          <div style={{ display: "flex", justifyContent: "flex-end", marginBottom: 8 }}>
            <Btn variant="ghost" onClick={clearStandards}>Clear</Btn>
          </div>
          <div style={{ overflowX: "auto" }}>
            <table style={{ borderCollapse: "collapse", width: "100%", minWidth: 560 }}>
              <thead>
                <tr>
                  <Th>Conc. (ng/mL)</Th><Th align="right">Rep 1</Th><Th align="right">Rep 2</Th><Th align="right">Rep 3</Th>
                  <Th align="right">Mean OD</Th><Th align="right">CV%</Th><Th align="right">Blank-sub OD</Th>
                </tr>
              </thead>
              <tbody>
                {stdRows.map((r, i) => (
                  <tr key={i} style={{ borderBottom: "1px solid #EEF0EC" }}>
                    <Td>
                      <NumIn value={r.conc} onChange={(e) => updateStdConc(i, e.target.value)}
                        onPaste={(e) => handleStdPaste(e, i, 0)} style={{ width: 90 }} />
                    </Td>
                    {[0, 1, 2].map((repIdx) => (
                      <Td key={repIdx} align="right">
                        <NumIn value={standards[i].od[repIdx]} onChange={(e) => updateStdOD(i, repIdx, e.target.value)}
                          onPaste={(e) => handleStdPaste(e, i, repIdx + 1)} style={{ width: 78, textAlign: "right" }} />
                      </Td>
                    ))}
                    <Td align="right">{fmt(r.rawMean)}</Td>
                    <Td align="right">
                      {Number.isFinite(r.cv) ? (r.cv > 20 ? <Badge tone="warn">{fmt(r.cv, 1)}</Badge> : fmt(r.cv, 1)) : "—"}
                    </Td>
                    <Td align="right">{fmt(r.sub)}</Td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <div style={{ marginTop: 14, display: "flex", gap: 10, alignItems: "center" }}>
            {fit ? (
              <>
                <Badge tone={r2Pass ? "pass" : "fail"}>R² = {fmt(fit.r2, 4)}</Badge>
                <span style={{ fontSize: 12.5, color: "var(--ink-soft)" }}>
                  {r2Pass
                    ? `meets the kit's acceptance criterion (≥ ${kit.r2Threshold}) — curve is valid.`
                    : `below the kit's acceptance criterion (≥ ${kit.r2Threshold}) — per §09 this run should be considered invalid.`}
                </span>
              </>
            ) : (
              <span style={{ fontSize: 12.5, color: "var(--ink-soft)" }}>Enter at least 4 standard OD values to fit the curve.</span>
            )}
          </div>

          {fit && (
            <div ref={curveChartRef} style={{ marginTop: 16, background: "var(--panel)", border: "1px solid var(--border)", borderRadius: 8, padding: "12px 16px 4px 4px" }}>
              <ResponsiveContainer width="100%" height={280}>
                <ComposedChart margin={{ top: 10, right: 20, bottom: 20, left: 34 }}>
                  <CartesianGrid stroke="#EEF0EC" />
                  <XAxis dataKey="x" type="number" domain={["auto", "auto"]}
                    tick={{ fontSize: 11, fontFamily: "var(--mono)" }}
                    label={{ value: "log₁₀ [dsRNA], ng/mL", position: "insideBottom", offset: -12, fontSize: 12, fill: "var(--ink-soft)" }} />
                  <YAxis dataKey="y" type="number"
                    tick={{ fontSize: 11, fontFamily: "var(--mono)" }}
                    label={{ value: "OD450 (blank-subtracted)", angle: -90, position: "insideLeft", offset: 8, style: { textAnchor: "middle" }, fontSize: 12, fill: "var(--ink-soft)" }} />
                  <Tooltip
                    formatter={(v, name) => [fmt(v, 4), name]}
                    labelFormatter={(v) => `log₁₀ x = ${fmt(v, 3)}`}
                    contentStyle={{ fontFamily: "var(--mono)", fontSize: 12, borderRadius: 6, border: "1px solid var(--border)" }} />
                  <Line data={curveData} dataKey="y" stroke="var(--accent)" strokeWidth={2} dot={false} isAnimationActive={false} name="4PL fit" />
                  <Scatter data={standardPoints} dataKey="y" fill="var(--ink)" name="Standards" />
                  <Scatter data={samplePoints} dataKey="y" fill="var(--fail)" shape="diamond" name="Samples" />
                </ComposedChart>
              </ResponsiveContainer>
            </div>
          )}
          {fit && (
            <div style={{ marginTop: 8 }}>
              <Btn variant="ghost" onClick={() => downloadChartPng(curveChartRef, "standard_curve.png")}>Download chart (PNG)</Btn>
            </div>
          )}
        </Step>

        {/* Step 3: Samples */}
        <Step n="03" title="Unknown samples"
          sub="Provide the dilution factor exactly as it was applied before loading. Total RNA concentration (optional) converts the result to % dsRNA (w/w); expected concentration (optional) computes % recovery for spiked QC samples. Paste a whole block copied from Excel (e.g. name + 3 OD columns for every sample at once) directly into any cell — extra rows are added automatically.">
          <div style={{ overflowX: "auto" }}>
            <table style={{ borderCollapse: "collapse", width: "100%", minWidth: 860 }}>
              <thead>
                <tr>
                  <Th>Sample name</Th><Th align="right">Rep 1</Th><Th align="right">Rep 2</Th><Th align="right">Rep 3</Th>
                  <Th align="right">Dilution</Th><Th align="right">Total RNA (ng/µL)</Th><Th align="right">Expected (ng/mL)</Th>
                  <Th>Status</Th><Th></Th>
                </tr>
              </thead>
              <tbody>
                {sampleResults.map((r, i) => (
                  <tr key={r.id} style={{ borderBottom: "1px solid #EEF0EC" }}>
                    <Td>
                      <NumIn value={r.name} inputMode="text" onChange={(e) => updateSample(r.id, { name: e.target.value })}
                        onPaste={(e) => handleSamplePaste(e, i, 0)}
                        style={{ width: 180, fontFamily: "var(--sans)" }} />
                    </Td>
                    {[0, 1, 2].map((repIdx) => (
                      <Td key={repIdx} align="right">
                        <NumIn value={r.od[repIdx]} onChange={(e) => updateSampleOD(r.id, repIdx, e.target.value)}
                          onPaste={(e) => handleSamplePaste(e, i, repIdx + 1)}
                          style={{ width: 76, textAlign: "right" }} />
                      </Td>
                    ))}
                    <Td align="right">
                      <NumIn value={r.dilution} onChange={(e) => updateSample(r.id, { dilution: e.target.value })}
                        onPaste={(e) => handleSamplePaste(e, i, 4)} style={{ width: 70, textAlign: "right" }} />
                    </Td>
                    <Td align="right">
                      <NumIn value={r.totalRna} onChange={(e) => updateSample(r.id, { totalRna: e.target.value })}
                        onPaste={(e) => handleSamplePaste(e, i, 5)} style={{ width: 90, textAlign: "right" }} />
                    </Td>
                    <Td align="right">
                      <NumIn value={r.expected} onChange={(e) => updateSample(r.id, { expected: e.target.value })}
                        onPaste={(e) => handleSamplePaste(e, i, 6)} style={{ width: 90, textAlign: "right" }} />
                    </Td>
                    <Td>{statusBadge(r)}</Td>
                    <Td>
                      {sampleResults.length > 1 && (
                        <Btn variant="ghost" onClick={() => removeSample(r.id)}>Remove</Btn>
                      )}
                    </Td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div style={{ marginTop: 14, display: "flex", gap: 10 }}>
            <Btn onClick={addSample}>+ Add sample</Btn>
            <Btn variant="ghost" onClick={clearAllSamples}>Clear all samples</Btn>
          </div>
        </Step>

        {/* Step 4: Sample comparison */}
        <Step n="04" title="Sample comparison"
          sub={`Compare dsRNA concentration across samples at a glance, and % dsRNA (w/w) wherever a total RNA concentration was provided. Bars in the % chart below are colored green (Pass) or red (Fail) against your pass threshold of ≤ ${fmt(passThreshold, 2)}% (set in Kit parameters).`}>
          {comparisonConcData.length === 0 ? (
            <span style={{ fontSize: 12.5, color: "var(--ink-soft)" }}>
              Enter OD values for at least one in-range sample above to see a comparison chart.
            </span>
          ) : (
            <>
              <div ref={concChartRef} style={{ background: "var(--panel)", border: "1px solid var(--border)", borderRadius: 8, padding: "12px 16px 4px 4px", marginBottom: 8 }}>
                <div style={{ fontSize: 12.5, fontWeight: 600, color: "var(--ink-soft)", padding: "0 0 4px 12px" }}>dsRNA concentration (ng/mL)</div>
                <ResponsiveContainer width="100%" height={Math.max(180, comparisonConcData.length * 42)}>
                  <BarChart data={comparisonConcData} layout="vertical" margin={{ top: 6, right: 24, bottom: 6, left: 10 }}>
                    <CartesianGrid stroke="#EEF0EC" horizontal={false} />
                    <XAxis type="number" tick={{ fontSize: 11, fontFamily: "var(--mono)" }} />
                    <YAxis dataKey="name" type="category" width={160} tick={{ fontSize: 11.5, fontFamily: "var(--sans)" }} />
                    <Tooltip formatter={(v) => [fmtSig(v, 4) + " ng/mL", "dsRNA conc."]} contentStyle={{ fontFamily: "var(--mono)", fontSize: 12, borderRadius: 6, border: "1px solid var(--border)" }} />
                    <Bar dataKey="conc" fill="var(--accent)" radius={[0, 4, 4, 0]} />
                  </BarChart>
                </ResponsiveContainer>
              </div>
              <div style={{ marginBottom: 16 }}>
                <Btn variant="ghost" onClick={() => downloadChartPng(concChartRef, "dsrna_concentration_by_sample.png")}>Download chart (PNG)</Btn>
              </div>

              {comparisonPctData.length > 0 && (
                <>
                  <div ref={pctChartRef} style={{ background: "var(--panel)", border: "1px solid var(--border)", borderRadius: 8, padding: "12px 16px 4px 4px" }}>
                    <div style={{ fontSize: 12.5, fontWeight: 600, color: "var(--ink-soft)", padding: "0 0 4px 12px" }}>% dsRNA (w/w) of total RNA — pass threshold ≤ {fmt(passThreshold, 2)}%</div>
                    <ResponsiveContainer width="100%" height={Math.max(180, comparisonPctData.length * 42)}>
                      <BarChart data={comparisonPctData} layout="vertical" margin={{ top: 6, right: 24, bottom: 6, left: 10 }}>
                        <CartesianGrid stroke="#EEF0EC" horizontal={false} />
                        <XAxis type="number" tick={{ fontSize: 11, fontFamily: "var(--mono)" }} unit="%" />
                        <YAxis dataKey="name" type="category" width={160} tick={{ fontSize: 11.5, fontFamily: "var(--sans)" }} />
                        <Tooltip formatter={(v) => [fmt(v, 4) + "%", "% dsRNA (w/w)"]} contentStyle={{ fontFamily: "var(--mono)", fontSize: 12, borderRadius: 6, border: "1px solid var(--border)" }} />
                        <Bar dataKey="pct" radius={[0, 4, 4, 0]}>
                          {comparisonPctData.map((d, i) => (
                            <Cell key={i} fill={d.pass === false ? "var(--fail)" : "var(--accent)"} />
                          ))}
                        </Bar>
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                  <div style={{ marginTop: 8 }}>
                    <Btn variant="ghost" onClick={() => downloadChartPng(pctChartRef, "pct_dsrna_by_sample.png")}>Download chart (PNG)</Btn>
                  </div>
                </>
              )}
              {comparisonPctData.length === 0 && (
                <div style={{ marginTop: 10, fontSize: 12.5, color: "var(--ink-soft)" }}>
                  Add a total RNA concentration (ng/µL) to any sample above to see % dsRNA (w/w) compared here.
                </div>
              )}
            </>
          )}
        </Step>

        {/* Step 5: Results */}
        <Step n="05" title="Results summary"
          sub={`Full table for records, plus export. Pass/Fail is judged against the pass threshold set in Kit parameters (currently ≤ ${fmt(passThreshold, 2)}% dsRNA w/w).`}>
          <div style={{ overflowX: "auto" }}>
            <table style={{ borderCollapse: "collapse", width: "100%", minWidth: 800 }}>
              <thead>
                <tr>
                  <Th>Sample</Th><Th align="right">Blank-sub OD</Th><Th align="right">CV%</Th>
                  <Th align="right">Dilution</Th><Th>Status</Th><Th align="right">dsRNA conc.</Th>
                  <Th align="right">% dsRNA</Th><Th>Pass/Fail</Th><Th align="right">Recovery</Th>
                </tr>
              </thead>
              <tbody>
                {sampleResults.map((r) => (
                  <tr key={r.id} style={{ borderBottom: "1px solid #EEF0EC" }}>
                    <Td mono={false}>{r.name}</Td>
                    <Td align="right">{fmt(r.sub)}</Td>
                    <Td align="right">{fmt(r.cv, 1)}</Td>
                    <Td align="right">{r.dilution}×</Td>
                    <Td>{statusBadge(r)}</Td>
                    <Td align="right">
                      {r.status === "ok"
                        ? (r.belowLoq ? `< ${fmtSig((kit.loqPg / 1000) * r.dilution * 1000, 3)} pg/mL` : `${fmtSig(r.finalConc, 4)} ng/mL`)
                        : "—"}
                    </Td>
                    <Td align="right">{r.pct !== null ? `${fmt(r.pct, 4)}%` : "—"}</Td>
                    <Td>{passBadge(r)}</Td>
                    <Td align="right">{r.recovery !== null ? `${fmt(r.recovery, 0)}%` : "—"}</Td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div style={{ marginTop: 14, display: "flex", gap: 10, flexWrap: "wrap" }}>
            <Btn variant="primary" onClick={exportCsv}>Export CSV</Btn>
            <Btn onClick={exportExcel}>Export Excel (.xlsx)</Btn>
            <Btn onClick={exportExcelWithCharts}>Export Excel + chart PNGs</Btn>
            <Btn onClick={exportSlides}>Export slides (.pptx)</Btn>
          </div>
          <div style={{ marginTop: 8, fontSize: 12, color: "var(--ink-soft)" }}>
            "Export slides" builds the .pptx file directly in your browser — no external library or network access needed.
          </div>
        </Step>

        <div style={{ borderTop: "1px solid var(--border)", paddingTop: 16, fontSize: 12, color: "var(--ink-soft)", lineHeight: 1.6 }}>
          Curve fit: four-parameter logistic on log₁₀(concentration) vs. blank-subtracted OD450, solved by
          Nelder–Mead least squares, matching §10 of the kit insert. "Above curve" samples should be diluted
          further and re-read; the tool then multiplies the re-read result by the new dilution factor when you
          update the row. This tool performs no regulatory determination — dsRNA acceptance criteria are
          process- and program-specific and should be set from your own characterization and toxicology data.
        </div>
      </div>
    </div>
  );
}

/* ---------------------------------------------------------------------- */
/* Mount                                                                   */
/* ---------------------------------------------------------------------- */

const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(React.createElement(DsRnaElisaAnalyzer));
