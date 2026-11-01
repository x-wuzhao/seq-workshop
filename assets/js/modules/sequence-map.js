var seqMapState = { zoom: 1, circular: false, rot: 0, sel: null, tracks: { seq: true, orf: true, enzyme: true, primer: true, mutation: true, feature: true, gc: true } };

function mapNiceStep(len, maxTicks) {
  var cands = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 20000, 50000];
  for (var i = 0; i < cands.length; i++) if (len / cands[i] <= maxTicks) return cands[i];
  return 100000;
}
function mapBaseColor(b) { return { A: '#2E8B4F', T: '#D64545', G: '#D98A16', C: '#2F7DC1' }[b] || '#7B8794'; }

function seqMapOrfs(seq) {
  var orfs = [];
  for (var f = 0; f < 3; f++) {
    findORFs(seq, f).forEach(function (o) { orfs.push({ s: o.start, e: o.end, aa: o.prot.length, strand: '+', frame: f + 1 }); });
    findORFs(reverseComplement(seq), f).forEach(function (o) {
      orfs.push({ s: seq.length - o.end, e: seq.length - o.start, aa: o.prot.length, strand: '-', frame: f + 1 });
    });
  }
  orfs = orfs.filter(function (o) { return o.aa >= 20; });
  if (orfs.length > 24) { orfs.sort(function (a, b) { return b.aa - a.aa; }); orfs = orfs.slice(0, 24); }
  orfs.sort(function (a, b) { return a.s - b.s; });
  return orfs;
}

function seqMapEnzymeTicks() {
  var ticks = [];
  (lastEnz || []).forEach(function (r) {
    if (!r.hits.length || r.enz[1].length < 6 || r.hits.length > 6) return;
    r.hits.forEach(function (h) { ticks.push({ pos: h.cut, name: r.enz[0] }); });
  });
  return ticks.slice(0, 80);
}

/* 环状视图：把同一批数据画到圆周上 */
var mapFeatures = [];
function mapFeatureList() { return mapFeatures; }

function parseGbLocation(loc) {
  if(!/^(?:\d+(?:\.\.\d+)?|complement\(\d+(?:\.\.\d+)?\))$/.test(loc))return null;
  var strand = /complement/i.test(loc) ? '-' : '+';
  var m = loc.match(/(\d+)\s*\.\.\s*(\d+)/);
  if (m) return { start: parseInt(m[1], 10), end: parseInt(m[2], 10), strand: strand };
  m = loc.match(/(\d+)/);
  if (m) return { start: parseInt(m[1], 10), end: parseInt(m[1], 10), strand: strand };
  return null;
}

function parseGenBank(text) {
  if((String(text).match(/^LOCUS\b/gm)||[]).length>1)return {error:'一次只能导入一个 GenBank 记录。',seq:'',features:[]};
  var lines = String(text).split(/\r?\n/);
  var seq = '', feats = [], inFeat = false, inOrigin = false, cur = null;
  lines.forEach(function (line) {
    if (/^LOCUS/.test(line)) { inFeat = false; return; }
    if (/^FEATURES/.test(line)) { inFeat = true; return; }
    if (/^ORIGIN/.test(line)) { inFeat = false; inOrigin = true; return; }
    if(/^\/\//.test(line)){inOrigin=false;return;}
    if (inOrigin) {
      var m = line.match(/[a-zA-Z]+/g);
      if (m) seq += m.join('');
      return;
    }
    if (!inFeat) return;
    var fk = line.match(/^ {5}(\S+)\s+(.+)$/);
    if (fk) {
      if (cur) feats.push(cur);
      cur = { key: fk[1], loc: fk[2].trim(), label: '' };
      return;
    }
    if (!cur) return;
    var q = line.match(/^\s+\/(\w+)=(.*)$/);
    if (q) {
      var val = q[2].replace(/^"|"$/g, '');
      if (q[1] === 'label' || q[1] === 'gene') cur.label = val;
      else if (q[1] === 'note' || q[1] === 'product') cur.note = (cur.note ? cur.note + ' ' : '') + val;
      return;
    }
    var cont = line.trim();
    if (cont && /^[<>a-zA-Z0-9(),.]+$/.test(cont) && !/^\//.test(cont) && cont.indexOf('"') < 0) {
      cur.loc += cont;
    }
  });
  if (cur) feats.push(cur);
  var mapped = [], skipped = 0;
  feats.forEach(function (f) {
    if (f.key === 'source') return;
    var loc = parseGbLocation(f.loc);
    if (!loc || loc.start<1 || loc.end>seq.length) {skipped++;return;}
    mapped.push({ start: Math.min(loc.start, loc.end), end: Math.max(loc.start, loc.end), strand: loc.strand, label: f.label || f.key, note: f.note || '', key: f.key });
  });
  var topology = 'linear';
  for (var i = 0; i < lines.length; i++) {
    if (/^LOCUS/.test(lines[i])) {
      if (/circular/i.test(lines[i])) topology = 'circular';
      break;
    }
  }
  return { seq: seq.toUpperCase(), features: mapped, topology: topology, skipped: skipped };
}

function renderMapFeatures() {
  var box = document.getElementById('mapFeatures');
  if (!box) return;
  if (!mapFeatures.length) {
    box.innerHTML = '<span class="mf-empty">还没有注释：框选一段后点「＋ 注释」即可添加；也可以导入 .gb 文件。</span>';
    return;
  }
  box.innerHTML = mapFeatures.map(function (f, i) {
    return '<span class="mf-chip" data-mf-sel="' + i + '"><b>' + SeqCore.escapeHTML(f.label) + '</b><em>' + f.start + '–' + f.end + '</em>'
      + (f.strand === '-' ? '<em class="mf-strand">反链</em>' : '')
      + '<button type="button" data-mf-del="' + i + '" aria-label="删除注释 ' + SeqCore.escapeHTML(f.label) + '">×</button></span>';
  }).join('');
}

function exportFeatureCSV() {
  if (!mapFeatures.length) {
    if (typeof showToast === 'function') showToast('还没有注释，先框选一段并用「＋ 注释」添加');
    return;
  }
  var rows = [['序号', '起点', '终点', '长度', '方向', '名称', '备注']];
  mapFeatures.forEach(function (f, i) {
    var s0 = Math.min(f.start, f.end), e0 = Math.max(f.start, f.end);
    rows.push([i + 1, s0, e0, e0 - s0 + 1, f.strand === '-' ? '反链' : '正链', f.label, f.note || '']);
  });
  exportCSV('序列工坊_注释.csv', rows);
  if (typeof showToast === 'function') showToast('已导出注释 CSV（' + mapFeatures.length + ' 条）');
}

function addMapFeature() {
  var sel = seqMapState.sel;
  if (!sel || Math.abs(sel.b - sel.a) < 1) {
    if (typeof showToast === 'function') showToast('请先在地图上拖动框选一段，再点「＋ 注释」');
    return;
  }
  var a = Math.min(sel.a, sel.b) + 1, b = Math.max(sel.a, sel.b);
  var label = window.prompt('给这段注释起个名字（例如：启动子 / CDS / 引物结合位点）', '注释 ' + (mapFeatures.length + 1));
  if (label === null) return;
  mapFeatures.push({ start: a, end: b, strand: '+', label: (label || '').trim() || ('注释 ' + (mapFeatures.length + 1)) });
  seqMapState.sel = null;
  renderSeqMap();
  renderMapFeatures();
  if (typeof scheduleWorkspaceSave === 'function') scheduleWorkspaceSave();
  if (typeof showToast === 'function') showToast('已添加注释：' + a + '–' + b);
}

function importGbText(text) {
  var r = parseGenBank(text);
  if(r.error){setMsg('err',r.error);return;}
  var gbParsed=enhancedParseInput(r.seq);
  if(gbParsed.error){setMsg('err',gbParsed.error);return;}
  if (!r.seq || r.seq.length < 3) {
    if (typeof showToast === 'function') showToast('没有从文件里读到有效序列');
    return;
  }
  window._pendingFeatures = r.features;
  document.getElementById('seq-input').value = r.seq;
  analyze();
  if(r.skipped)setMsg('info','已载入序列；'+r.skipped+' 条复杂或越界注释未导入。当前仅支持连续区间及其 complement 位置。');
  if (r.topology === 'circular') {
    seqMapState.circular = true;
    var lin = document.getElementById('mapLinear');
    var cir = document.getElementById('mapCircular');
    if (lin) lin.classList.remove('active');
    if (cir) cir.classList.add('active');
    var topo = document.getElementById('digestTopo');
    if (topo) {
      topo.value = 'circular';
      if (lastSeq && lastEnz) renderDigest(lastSeq, lastEnz);
    }
    renderSeqMap();
  }
  if (typeof showToast === 'function') showToast('已导入 ' + r.seq.length + ' bp · ' + r.features.length + ' 条注释' + (r.topology === 'circular' ? ' · 已切换环状视图' : ''));
}

function importGbFile(file) {
  var fr = new FileReader();
  fr.onload = function () { importGbText(fr.result); };
  fr.readAsText(file);
}

function applyPrimerDrag(dir, delta, base) {
  var p = (window._lastPrimers || [])[0];
  var seq = window._currentSeq || '';
  if (!p || !seq) return;
  var N = seq.length;
  if (dir === 'f') {
    var ns = Math.max(0, Math.min(p.rEnd - p.fLen - 1, base.fStart - 1 + delta));
    p.fStart = ns + 1;
    p.fSeq = seq.substr(ns, p.fLen);
    p.fTm = tmNN(p.fSeq);
    p.fGC = gcOf(p.fSeq);
    p.product = p.rEnd - ns;
  } else {
    var ne = Math.max(p.fStart + p.fLen + 1, Math.min(N, base.rEnd + delta));
    p.rEnd = ne;
    p.rStart = ne - p.rLen + 1;
    p.rSeq = reverseComplement(seq.substr(ne - p.rLen, p.rLen));
    p.rTm = tmNN(p.rSeq);
    p.rGC = gcOf(p.rSeq);
    p.product = ne - (p.fStart - 1);
  }
  p.manual = true;
}

function renderSeqMapCircular(seq) {
  var svg = document.getElementById('seqMap');
  var host = document.getElementById('mapScroll');
  if (!svg || !host) return;
  var T = seqMapState.tracks;
  var side = Math.min(700, Math.max(430, (host.clientWidth || 700) - 24));
  var W = side, H = side;
  var cx = W / 2, cy = H / 2, R = side * 0.28;
  var N = seq.length;
  var rot = seqMapState.rot || 0;
  function A(p) { return -Math.PI / 2 + rot + 2 * Math.PI * (p / N); }
  function P(p, r) { return [cx + r * Math.cos(A(p)), cy + r * Math.sin(A(p))]; }
  function f(pt) { return pt[0].toFixed(1) + ' ' + pt[1].toFixed(1); }
  function band(s, e, rIn, rOut, fill, op) {
    var a1 = A(s), a2 = A(e);
    var large = ((a2 - a1 + 2 * Math.PI) % (2 * Math.PI)) > Math.PI ? 1 : 0;
    return '<path d="M' + f(P(s, rOut)) + ' A' + rOut + ' ' + rOut + ' 0 ' + large + ' 1 ' + f(P(e, rOut))
      + ' L' + f(P(e, rIn)) + ' A' + rIn + ' ' + rIn + ' 0 ' + large + ' 0 ' + f(P(s, rIn)) + ' Z" fill="' + fill + '" opacity="' + (op || 1) + '"/>';
  }
  function wedge(p, r, color, reverse) {
    var a = reverse ? A(p) - 0.06 : A(p) + 0.06;
    var base = reverse ? A(p) + 0.03 : A(p) - 0.03;
    var t1 = [cx + (r - 5.5) * Math.cos(base), cy + (r - 5.5) * Math.sin(base)];
    var t2 = [cx + (r + 5.5) * Math.cos(base), cy + (r + 5.5) * Math.sin(base)];
    var tip = [cx + r * Math.cos(a), cy + r * Math.sin(a)];
    return '<path d="M' + f(t1) + ' L' + f(tip) + ' L' + f(t2) + ' Z" fill="' + color + '"/>';
  }
  var out = [];
  out.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + R + '" fill="none" stroke="var(--line)" stroke-width="1.6"/>');
  out.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + (R - 20) + '" fill="none" stroke="var(--line)" stroke-dasharray="2 5" opacity=".5"/>');

  /* 刻度与数字 */
  var step = mapNiceStep(N, 12);
  for (var p = 0; p < N; p += step) {
    var t1 = P(p, R), t2 = P(p, R + 5);
    out.push('<line x1="' + t1[0].toFixed(1) + '" y1="' + t1[1].toFixed(1) + '" x2="' + t2[0].toFixed(1) + '" y2="' + t2[1].toFixed(1) + '" stroke="var(--line)"/>');
    var lp = P(p, R + 13);
    out.push('<text x="' + lp[0].toFixed(1) + '" y="' + (lp[1] + 3).toFixed(1) + '" font-size="9" fill="var(--ink-soft)" text-anchor="middle">' + p + '</text>');
  }

  /* ORF */
  if (T.orf) {
    seqMapOrfs(seq).forEach(function (o) {
      var col = o.strand === '+' ? '#0E7C7B' : '#B9781B';
      var rIn = o.strand === '+' ? R + 17 : R - 29;
      var rOut = o.strand === '+' ? R + 27 : R - 19;
      out.push('<g data-jump="#nav-frames" style="cursor:pointer"><title>ORF ' + o.strand + o.frame + ' · ' + (o.s + 1) + '–' + o.e + ' · ' + o.aa + ' aa</title>'
        + band(o.s, o.e, rIn, rOut, col, .9) + wedge(o.e, (rIn + rOut) / 2, col, o.strand === '-') + '</g>');
    });
  }

  /* 酶切位点 + 名称（数量不多时直接标注） */
  var ticks = seqMapEnzymeTicks();
  if (T.enzyme) {
    ticks.forEach(function (t) {
      out.push('<g data-jump="#nav-enzyme" style="cursor:pointer"><title>' + t.name + ' · 切点 ' + t.pos + ' bp</title>'
        + '<line x1="' + P(t.pos, R - 6)[0].toFixed(1) + '" y1="' + P(t.pos, R - 6)[1].toFixed(1) + '" x2="' + P(t.pos, R + 6)[0].toFixed(1) + '" y2="' + P(t.pos, R + 6)[1].toFixed(1) + '" stroke="#C43E3E" stroke-width="2"/></g>');
    });
    if (ticks.length <= 10) {
      ticks.forEach(function (t) {
        var lp = P(t.pos, R - 50);
        var deg = (A(t.pos) * 180 / Math.PI + 90) % 360;
        if (deg > 90 && deg < 270) deg += 180;
        out.push('<text x="' + lp[0].toFixed(1) + '" y="' + (lp[1] + 3).toFixed(1) + '" font-size="9" fill="#8A3A3A" text-anchor="middle" transform="rotate(' + deg.toFixed(1) + ' ' + lp[0].toFixed(1) + ' ' + lp[1].toFixed(1) + ')" pointer-events="none">' + t.name + '</text>');
      });
    }
  }

  /* 用户注释 */
  if (T.feature && mapFeatures.length) {
    mapFeatures.forEach(function (f) {
      var s0 = Math.min(f.start, f.end) - 1, e0 = Math.max(f.start, f.end);
      var col = f.strand === '-' ? '#B9781B' : '#3E7CB1';
      out.push('<g style="cursor:pointer"><title>' + SeqCore.escapeHTML(f.label) + ' · ' + f.start + '–' + f.end + '</title>'
        + band(s0, e0, R - 39, R - 33, col, .9) + '</g>');
    });
  }

  /* 突变标记 */
  if (T.mutation && window._seqMapMutation) {
    var mm = window._seqMapMutation;
    var mp = P(mm.pos, R);
    out.push('<g data-jump="#nav-mutation" style="cursor:pointer"><title>' + mm.label + ' · ' + mm.detail + '</title>'
      + '<circle cx="' + mp[0].toFixed(1) + '" cy="' + mp[1].toFixed(1) + '" r="4.5" fill="#D64545" stroke="#fff" stroke-width="1.2"/></g>');
  }

  /* 引物 */
  if (T.primer) {
    var pr = (window._lastPrimers || [])[0];
    if (pr) {
      out.push('<g data-jump="#nav-primer" style="cursor:pointer"><title>正向引物 ' + pr.fStart + '–' + (pr.fStart + pr.fLen - 1) + '</title>'
        + band(pr.fStart - 1, pr.fStart - 1 + pr.fLen, R + 31, R + 38, '#7A4FA3', .95) + wedge(pr.fStart - 1 + pr.fLen, R + 34.5, '#7A4FA3', false) + '</g>');
      out.push('<g data-jump="#nav-primer" style="cursor:pointer"><title>反向引物 ' + pr.rStart + '–' + pr.rEnd + '</title>'
        + band(pr.rStart - 1, pr.rEnd, R + 31, R + 38, '#7A4FA3', .95) + wedge(pr.rStart - 1, R + 34.5, '#7A4FA3', true) + '</g>');
    }
  }

  /* GC 环 */
  if (T.gc) {
    var gw = Math.max(12, Math.round(N / 70));
    var gstep = Math.max(1, Math.round(gw / 4));
    var gpts = [];
    for (var g = 0; g + gw <= N; g += gstep) {
      var gc = 0;
      for (var k = g; k < g + gw; k++) { var c = seq.charAt(k); if (c === 'G' || c === 'C') gc++; }
      gpts.push({ p: g + gw / 2, v: gc / gw });
    }
    out.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + (R + 46) + '" fill="none" stroke="var(--line)" stroke-dasharray="3 4" opacity=".7"/>');
    if (gpts.length > 2) {
      var d = 'M' + f(P(gpts[0].p, R + 46));
      gpts.forEach(function (pt) { d += ' L' + f(P(pt.p, R + 46 + (pt.v - 0.5) * 26)); });
      for (var gi = gpts.length - 1; gi >= 0; gi--) { d += ' L' + f(P(gpts[gi].p, R + 46)); }
      d += ' Z';
      out.push('<path d="' + d + '" fill="rgba(127,179,213,.3)" stroke="#5B9BD5" stroke-width="1.2"/>');
    }
  }

  out.push('<text x="' + cx + '" y="' + (cy - 4) + '" text-anchor="middle" font-size="22" font-weight="700" fill="var(--ink)">' + N + '<tspan font-size="11" font-weight="400" fill="var(--ink-soft)"> bp</tspan></text>');
  out.push('<text x="' + cx + '" y="' + (cy + 16) + '" text-anchor="middle" font-size="11" fill="var(--ink-soft)">环状质粒视图</text>');

  svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
  svg.setAttribute('width', W);
  svg.setAttribute('height', H);
  svg.innerHTML = out.join('');

  var hint = document.getElementById('mapHint');
  if (hint) hint.textContent = '环状视图：拖动可旋转；外圈为 ORF 正向、内圈为反向，红线为酶切位点，外环折线为 GC 含量；点击特征可跳转，切回“线性”可逐碱基查看。';
}

function renderSeqMap(seq) {
  var svg = document.getElementById('seqMap');
  var host = document.getElementById('mapScroll');
  if (!svg || !host) return;
  seq = seq || window._currentSeq || '';
  if (!seq) { svg.innerHTML = ''; return; }
  var T = seqMapState.tracks;
  var ztools = document.getElementById('mapZoomTools');
  if (ztools) ztools.hidden = !!seqMapState.circular;
  if (seqMapState.circular) { renderSeqMapCircular(seq); return; }

  var avail = Math.max(560, (host.clientWidth || 900) - 22);
  var zoom = seqMapState.zoom;
  var px = (avail / seq.length) * zoom;
  var padL = 58, padR = 16;
  var W = Math.max(avail, Math.round(seq.length * px) + padL + padR);
  var N = seq.length;
  function X(i) { return padL + i * px; }
  var out = [];
  var T2 = seqMapState.tracks;

  /* --- 标尺（固定在最上方） --- */
  var step = mapNiceStep(N, 8);
  out.push('<line x1="' + padL + '" y1="20" x2="' + (W - padR) + '" y2="20" stroke="var(--line)"/>');
  for (var p = 0; p <= N; p += step) {
    var tx = X(p);
    out.push('<line x1="' + tx + '" y1="14" x2="' + tx + '" y2="20" stroke="var(--line)"/>');
    var anc = p === 0 ? 'start' : (p + step > N ? 'end' : 'middle');
    out.push('<text x="' + tx + '" y="10" font-size="9.5" fill="var(--ink-soft)" text-anchor="' + anc + '">' + p + '</text>');
  }
  out.push('<text x="10" y="14" font-size="10.5" fill="var(--ink-soft)">标尺</text>');

  /* --- 动态轨道布局：关闭的轨道不再占行高 --- */
  var specs = [];
  if (T2.seq) specs.push({ k: 'seq', label: '碱基', h: 18 });
  if (T2.orf) specs.push({ k: 'orf', label: 'ORF', h: 24 });
  if (T2.enzyme) specs.push({ k: 'enzyme', label: '酶切', h: 16 });
  if (T2.primer) specs.push({ k: 'primer', label: '引物', h: 16 });
  if (T2.mutation) specs.push({ k: 'mutation', label: '突变', h: 16 });
  if (T2.feature && mapFeatures.length) specs.push({ k: 'feature', label: '注释', h: 18 });
  if (T2.gc) specs.push({ k: 'gc', label: 'GC', h: 46 });

  var cur = 30, ys = {};
  specs.forEach(function (sp) { ys[sp.k] = cur; cur += sp.h + 8; });
  var H = cur + 4;
  seqMapState.metrics = { padL: padL, px: px, W: W, N: N, H: H };

  specs.forEach(function (sp) {
    out.push('<text x="10" y="' + (ys[sp.k] + sp.h / 2 + 4) + '" font-size="10.5" fill="var(--ink-soft)">' + sp.label + '</text>');
  });

  /* --- 碱基 --- */
  if (T2.seq) {
    var y0 = ys.seq;
    if (px >= 6.2 && N <= 3000) {
      for (var i = 0; i < N; i++) {
        out.push('<text x="' + (X(i) + px / 2) + '" y="' + (y0 + 13) + '" font-size="10.5" text-anchor="middle" fill="' + mapBaseColor(seq.charAt(i)) + '">' + seq.charAt(i) + '</text>');
      }
    } else if (N <= 3000) {
      for (var j = 0; j < N; j++) {
        out.push('<rect x="' + X(j) + '" y="' + y0 + '" width="' + (Math.max(px, 0.6) + 0.25) + '" height="16" fill="' + mapBaseColor(seq.charAt(j)) + '" opacity=".92"/>');
      }
    } else {
      out.push('<rect x="' + padL + '" y="' + y0 + '" width="' + (W - padL - padR) + '" height="16" fill="#0E7C7B" opacity=".55"/>');
      out.push('<text x="' + (padL + 10) + '" y="' + (y0 + 12) + '" font-size="10.5" fill="#fff">序列过长，已切换为整体色带（共 ' + N + ' bp）</text>');
    }
  }

  /* --- ORF --- */
  if (T2.orf) {
    var by0 = ys.orf;
    var orfs = seqMapOrfs(seq);
    var drawOrfs = function (list, baseY, color, strand) {
      var laneEnd = [0, 0];
      list.forEach(function (o) {
        var x1 = X(o.s), w = Math.max(2.5, X(o.e) - x1);
        var lane = 0;
        while (lane < laneEnd.length && laneEnd[lane] > x1 - 2) lane++;
        if (lane >= laneEnd.length) return;
        laneEnd[lane] = x1 + w;
        var y = baseY + lane * 11;
        out.push('<g data-jump="#nav-frames" style="cursor:pointer"><title>ORF ' + strand + o.frame + ' · ' + (o.s + 1) + '–' + o.e + ' · ' + o.aa + ' aa</title>'
          + '<rect x="' + x1 + '" y="' + y + '" width="' + w + '" height="9" rx="2.5" fill="' + color + '" opacity=".92"/></g>');
        if (w > 58) out.push('<text x="' + (x1 + 5) + '" y="' + (y + 7) + '" font-size="8.5" fill="#fff" pointer-events="none">' + strand + o.frame + ' · ' + o.aa + 'aa</text>');
      });
    };
    drawOrfs(orfs.filter(function (o) { return o.strand === '+'; }), by0, '#0E7C7B', '+');
    drawOrfs(orfs.filter(function (o) { return o.strand === '-'; }), by0 + 11, '#B9781B', '-');
  }

  /* --- 酶切 --- */
  if (T2.enzyme) {
    var ey = ys.enzyme;
    var ticks = seqMapEnzymeTicks();
    ticks.forEach(function (t) {
      out.push('<g data-jump="#nav-enzyme" style="cursor:pointer"><title>' + t.name + ' · 切点 ' + t.pos + ' bp</title>'
        + '<line x1="' + X(t.pos) + '" y1="' + ey + '" x2="' + X(t.pos) + '" y2="' + (ey + 14) + '" stroke="#C43E3E" stroke-width="2"/></g>');
    });
    if (!ticks.length) out.push('<text x="' + padL + '" y="' + (ey + 11) + '" font-size="10.5" fill="var(--ink-soft)">（未收录酶在该序列上有切点）</text>');
  }

  /* --- 引物 --- */
  if (T2.primer) {
    var py = ys.primer;
    var pr = (window._lastPrimers || [])[0];
    if (pr) {
      var arrow = function (x1, x2, dir, title) {
        var w = Math.max(4, x2 - x1), y = py, h = 10;
        var d = dir === 'f'
          ? 'M' + x1 + ',' + y + ' H' + (x1 + w - 6) + ' L' + (x1 + w) + ',' + (y + h / 2) + ' L' + (x1 + w - 6) + ',' + (y + h) + ' H' + x1 + ' Z'
          : 'M' + x2 + ',' + y + ' H' + (x1 + 6) + ' L' + x1 + ',' + (y + h / 2) + ' L' + (x1 + 6) + ',' + (y + h) + ' H' + x2 + ' Z';
        out.push('<g data-jump="#nav-primer" data-primer="' + dir + '" style="cursor:grab"><title>' + title + '</title>'
          + '<rect x="' + (x1 - 3) + '" y="' + (y - 4) + '" width="' + (w + 6) + '" height="18" fill="transparent"/>'
          + '<path d="' + d + '" fill="#7A4FA3"/></g>');
      };
      arrow(X(pr.fStart - 1), X(pr.fStart - 1 + pr.fLen), 'f', '正向引物 ' + pr.fStart + '–' + (pr.fStart + pr.fLen - 1) + '（可拖动）');
      arrow(X(pr.rStart - 1), X(pr.rEnd), 'r', '反向引物 ' + pr.rStart + '–' + pr.rEnd + '（可拖动）');
    } else {
      out.push('<text x="' + padL + '" y="' + (py + 11) + '" font-size="10.5" fill="var(--ink-soft)">（在“引物设计”里生成候选后，这里会显示引物位置）</text>');
    }
  }

  /* --- 突变 --- */
  if (T2.mutation) {
    var my = ys.mutation;
    var mut = window._seqMapMutation;
    if (mut) {
      var mx = X(mut.pos);
      out.push('<g data-jump="#nav-mutation" style="cursor:pointer"><title>' + mut.label + ' · ' + mut.detail + '</title>'
        + '<path d="M' + mx + ',' + (my + 2) + ' l-5,-9 h10 Z" fill="#D64545"/>'
        + '<line x1="' + mx + '" y1="' + (my + 2) + '" x2="' + mx + '" y2="' + (my + 15) + '" stroke="#D64545" stroke-width="1.5"/></g>');
      if (px > 1.2) out.push('<text x="' + (mx + 7) + '" y="' + (my + 12) + '" font-size="9" fill="#A33" pointer-events="none">' + mut.label + '</text>');
    } else {
      out.push('<text x="' + padL + '" y="' + (my + 11) + '" font-size="10.5" fill="var(--ink-soft)">（在“突变模拟”里应用一次突变后，这里会出现标记）</text>');
    }
  }

  /* --- 用户注释 --- */
  if (T2.feature && mapFeatures.length) {
    var fy = ys.feature;
    mapFeatures.forEach(function (f, idx) {
      var s0 = Math.min(f.start, f.end) - 1, e0 = Math.max(f.start, f.end);
      var x1 = X(s0), w = Math.max(3, X(e0) - x1);
      var col = f.strand === '-' ? '#B9781B' : '#3E7CB1';
      out.push('<g data-feature="' + idx + '" style="cursor:move"><title>' + SeqCore.escapeHTML(f.label) + ' · ' + f.start + '–' + f.end + '（拖动移动，拖两端改边界，双击改名）</title>'
        + '<rect x="' + x1 + '" y="' + fy + '" width="' + w + '" height="10" rx="2" fill="' + col + '" opacity=".85"/></g>');
      out.push('<rect data-feature-handle="' + idx + '" data-edge="s" x="' + (x1 - 3) + '" y="' + (fy - 2) + '" width="7" height="14" fill="transparent" style="cursor:ew-resize"/>');
      out.push('<rect data-feature-handle="' + idx + '" data-edge="e" x="' + (x1 + w - 4) + '" y="' + (fy - 2) + '" width="7" height="14" fill="transparent" style="cursor:ew-resize"/>');
      if (w > 52) out.push('<text x="' + (x1 + 4) + '" y="' + (fy + 8) + '" font-size="8.5" fill="#fff" pointer-events="none">' + SeqCore.escapeHTML(f.label) + '</text>');
    });
  }

  /* --- GC --- */
  if (T2.gc) {
    var gy0 = ys.gc + 40, gy1 = ys.gc + 6;
    var gw = Math.max(12, Math.round(N / 80));
    var gstep = Math.max(1, Math.round(gw / 4));
    var pts = [];
    for (var g = 0; g + gw <= N; g += gstep) {
      var gc = 0;
      for (var k = g; k < g + gw; k++) { var c = seq.charAt(k); if (c === 'G' || c === 'C') gc++; }
      pts.push({ x: X(g + gw / 2), v: gc / gw });
    }
    if (pts.length > 1) {
      var midY = gy0 - (gy0 - gy1) * 0.5;
      out.push('<line x1="' + padL + '" y1="' + midY + '" x2="' + (W - padR) + '" y2="' + midY + '" stroke="var(--line)" stroke-dasharray="3 3"/>');
      var d = 'M' + pts[0].x + ',' + gy0;
      pts.forEach(function (pt) { d += ' L' + pt.x + ',' + (gy0 - (gy0 - gy1) * pt.v); });
      d += ' L' + pts[pts.length - 1].x + ',' + gy0 + ' Z';
      out.push('<path d="' + d + '" fill="rgba(127,179,213,.35)"/>');
      var dl = 'M' + pts[0].x + ',' + (gy0 - (gy0 - gy1) * pts[0].v);
      pts.forEach(function (pt) { dl += ' L' + pt.x + ',' + (gy0 - (gy0 - gy1) * pt.v); });
      out.push('<path d="' + dl + '" fill="none" stroke="#5B9BD5" stroke-width="1.4"/>');
    }
    out.push('<text x="' + padL + '" y="' + (gy0 + 8) + '" font-size="9" fill="var(--ink-soft)">0%</text>');
    out.push('<text x="' + padL + '" y="' + (gy1 - 2) + '" font-size="9" fill="var(--ink-soft)">100%</text>');
  }

  /* --- 选区高亮 --- */
  var sel = seqMapState.sel;
  if (sel) {
    var a = Math.min(sel.a, sel.b), b = Math.max(sel.a, sel.b);
    if (b - a >= 1) {
      out.push('<rect x="' + X(a) + '" y="24" width="' + Math.max(1, X(b) - X(a)) + '" height="' + (H - 30) + '" fill="rgba(14,124,123,.10)" pointer-events="none"/>');
      out.push('<line x1="' + X(a) + '" y1="24" x2="' + X(a) + '" y2="' + (H - 6) + '" stroke="#0E7C7B" stroke-dasharray="3 3" pointer-events="none"/>');
      out.push('<line x1="' + X(b) + '" y1="24" x2="' + X(b) + '" y2="' + (H - 6) + '" stroke="#0E7C7B" stroke-dasharray="3 3" pointer-events="none"/>');
    }
  }

  svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
  svg.setAttribute('width', W);
  svg.setAttribute('height', H);
  svg.innerHTML = out.join('');

  var zl = document.getElementById('mapZoomLabel');
  if (zl) zl.textContent = zoom.toFixed(1) + '×';
  var hint = document.getElementById('mapHint');
  if (hint) {
    if (sel && Math.abs(sel.b - sel.a) >= 1) {
      var s1 = Math.min(sel.a, sel.b), s2 = Math.max(sel.a, sel.b);
      hint.textContent = '已选中 ' + (s1 + 1) + '–' + s2 + '（' + (s2 - s1) + ' bp）· 点击空白处取消；继续拖动可调整范围。';
    } else if (px >= 6.2) {
      hint.textContent = '已放大到单碱基视图：拖动可框选区段，点击 ORF / 酶切 / 引物 / 突变标记可跳到对应面板。';
    } else {
      hint.textContent = '提示：拖动框选区段；点 + 放大到约 ' + (6.2 / px).toFixed(1) + '× 可显示单个碱基；点击特征可跳转。';
    }
  }
}

function initSeqMap() {
  var svg = document.getElementById('seqMap');
  if (!svg || svg.getAttribute('data-bound') === '1') return;
  svg.setAttribute('data-bound', '1');

  /* 悬停提示 */
  var tip = document.getElementById('mapTip');
  if (!tip) {
    tip = document.createElement('div');
    tip.id = 'mapTip';
    tip.className = 'map-tip';
    tip.setAttribute('aria-hidden', 'true');
    document.body.appendChild(tip);
  }
  function showTip(text, x, y) {
    if (!text) { tip.classList.remove('show'); return; }
    tip.textContent = text;
    tip.style.left = Math.min(window.innerWidth - 300, x + 14) + 'px';
    tip.style.top = (y + 16) + 'px';
    tip.classList.add('show');
  }
  svg.addEventListener('mousemove', function (e) {
    var g = e.target && e.target.closest ? e.target.closest('[data-jump]') : null;
    if (!g) { showTip(''); return; }
    var t = g.querySelector('title');
    showTip(t ? t.textContent : '', e.clientX, e.clientY);
  });
  svg.addEventListener('mouseleave', function () { showTip(''); });

  /* 点击特征跳转（拖动后不触发） */
  svg.addEventListener('click', function (e) {
    if (Date.now() - (seqMapState._dragEnd || 0) < 250) return;
    var g = e.target && e.target.closest ? e.target.closest('[data-jump]') : null;
    if (!g) return;
    var target = document.querySelector(g.getAttribute('data-jump'));
    if (target && window.SeqModuleRouter) window.SeqModuleRouter.navigate(g.getAttribute('data-jump'));
    else if(target)target.scrollIntoView({behavior:'smooth',block:'start'});
  });

  /* 拖动：线性=框选，环状=旋转 */
  var drag = null;
  function toPos(e) {
    var m = seqMapState.metrics;
    if (!m) return 0;
    var r = svg.getBoundingClientRect();
    var vbW = (svg.viewBox && svg.viewBox.baseVal && svg.viewBox.baseVal.width) || r.width || 1;
    var sx = (e.clientX - r.left) * (vbW / Math.max(1, r.width));
    var p = Math.round((sx - m.padL) / (m.px || 1));
    return Math.max(0, Math.min(m.N, p));
  }
  svg.addEventListener('pointerdown', function (e) {
    if (e.button && e.button !== 0) return;
    if (e.target.closest && e.target.closest('[data-jump],[data-feature],[data-feature-handle]')) return;
    drag = { sx: e.clientX, sy: e.clientY, x: e.clientX, y: e.clientY, moved: 0, a: null, rot0: seqMapState.rot || 0 };
    if (!seqMapState.circular) {
      var p = toPos(e);
      drag.a = p;
      seqMapState.sel = { a: p, b: p };
      renderSeqMap();
    }
    try { svg.setPointerCapture(e.pointerId); } catch (err) {}
  });
  svg.addEventListener('pointermove', function (e) {
    if (!drag) return;
    drag.moved += Math.abs(e.clientX - drag.x) + Math.abs(e.clientY - drag.y);
    drag.x = e.clientX; drag.y = e.clientY;
    if (seqMapState.circular) {
      seqMapState.rot = (drag.rot0 || 0) + (e.clientX - drag.sx) * 0.006;
      renderSeqMap();
    } else {
      seqMapState.sel = { a: drag.a, b: toPos(e) };
      renderSeqMap();
    }
  });
  function endDrag(e) {
    if (!drag) return;
    var moved = drag.moved;
    drag = null;
    seqMapState._dragEnd = moved > 6 ? Date.now() : 0;
    if (moved < 6 && !seqMapState.circular) { seqMapState.sel = null; renderSeqMap(); }
    try { if (e && e.pointerId !== undefined) svg.releasePointerCapture(e.pointerId); } catch (err) {}
  }
  svg.addEventListener('pointerup', endDrag);
  svg.addEventListener('pointercancel', endDrag);

  /* 拖动引物：只在线性视图生效 */
  var pdrag = null;
  svg.addEventListener('pointerdown', function (e) {
    var pg = e.target && e.target.closest ? e.target.closest('[data-primer]') : null;
    if (!pg || seqMapState.circular) return;
    var p = (window._lastPrimers || [])[0];
    if (!p) return;
    pdrag = { dir: pg.getAttribute('data-primer'), sx: e.clientX, moved: 0, base: { fStart: p.fStart, rEnd: p.rEnd } };
    e.preventDefault();
    try { svg.setPointerCapture(e.pointerId); } catch (err) {}
  });
  svg.addEventListener('pointermove', function (e) {
    if (!pdrag) return;
    var m = seqMapState.metrics;
    if (!m) return;
    var r = svg.getBoundingClientRect();
    var vbW = (svg.viewBox && svg.viewBox.baseVal && svg.viewBox.baseVal.width) || r.width || 1;
    var deltaBp = ((e.clientX - pdrag.sx) * (vbW / Math.max(1, r.width))) / (m.px || 1);
    pdrag.moved = Math.abs(e.clientX - pdrag.sx);
    applyPrimerDrag(pdrag.dir, Math.round(deltaBp), pdrag.base);
    renderSeqMap();
  });
  function endPrimerDrag() {
    if (!pdrag) return;
    var moved = pdrag.moved;
    pdrag = null;
    if (moved > 6) seqMapState._dragEnd = Date.now();
    var p = (window._lastPrimers || [])[0];
    if (p && typeof showToast === 'function') {
      showToast('引物已移动：Tm ' + p.fTm.toFixed(1) + ' / ' + p.rTm.toFixed(1) + ' °C，产物 ' + p.product + ' bp');
    }
  }
  svg.addEventListener('pointerup', endPrimerDrag);
  svg.addEventListener('pointercancel', endPrimerDrag);

  /* 拖动注释：整体移动 / 改边界；双击改名 */
  var fdrag = null;
  svg.addEventListener('pointerdown', function (e) {
    var h = e.target.closest ? e.target.closest('[data-feature-handle]') : null;
    var bar = e.target.closest ? e.target.closest('[data-feature]') : null;
    if (!h && !bar) return;
    if (seqMapState.circular) return;
    var idx = parseInt((h || bar).getAttribute(h ? 'data-feature-handle' : 'data-feature'), 10);
    var f = mapFeatures[idx];
    if (!f) return;
    fdrag = {
      idx: idx,
      mode: h ? ('edge-' + h.getAttribute('data-edge')) : 'move',
      sx: e.clientX,
      start: Math.min(f.start, f.end),
      end: Math.max(f.start, f.end),
      moved: 0
    };
    e.preventDefault();
    try { svg.setPointerCapture(e.pointerId); } catch (err) {}
  });
  svg.addEventListener('pointermove', function (e) {
    if (!fdrag) return;
    var m = seqMapState.metrics;
    if (!m) return;
    var r = svg.getBoundingClientRect();
    var vbW = (svg.viewBox && svg.viewBox.baseVal && svg.viewBox.baseVal.width) || r.width || 1;
    var delta = Math.round(((e.clientX - fdrag.sx) * (vbW / Math.max(1, r.width))) / (m.px || 1));
    fdrag.moved = Math.abs(e.clientX - fdrag.sx);
    var f = mapFeatures[fdrag.idx];
    if (!f) return;
    if (fdrag.mode === 'move') {
      f.start = fdrag.start + delta;
      f.end = fdrag.end + delta;
    } else if (fdrag.mode === 'edge-s') {
      f.start = Math.min(fdrag.end - 1, fdrag.start + delta);
    } else {
      f.end = Math.max(fdrag.start + 1, fdrag.end + delta);
    }
    f.start = Math.max(1, Math.min(m.N, f.start));
    f.end = Math.max(1, Math.min(m.N, f.end));
    if (f.start > f.end) { var t = f.start; f.start = f.end; f.end = t; }
    renderSeqMap();
  });
  function endFeatureDrag() {
    if (!fdrag) return;
    var moved = fdrag.moved;
    fdrag = null;
    if (moved > 4) {
      seqMapState._dragEnd = Date.now();
      renderMapFeatures();
      if (typeof scheduleWorkspaceSave === 'function') scheduleWorkspaceSave();
    }
  }
  svg.addEventListener('pointerup', endFeatureDrag);
  svg.addEventListener('pointercancel', endFeatureDrag);
  svg.addEventListener('dblclick', function (e) {
    var bar = e.target.closest ? e.target.closest('[data-feature]') : null;
    if (!bar) return;
    var idx = parseInt(bar.getAttribute('data-feature'), 10);
    var f = mapFeatures[idx];
    if (!f) return;
    var label = window.prompt('修改注释名称', f.label);
    if (label === null) return;
    f.label = (label || '').trim() || f.label;
    renderSeqMap();
    renderMapFeatures();
  });

  /* 缩放 */
  function zoomBy(factor) {
    var host = document.getElementById('mapScroll');
    var before = host ? (host.scrollLeft + host.clientWidth / 2) / Math.max(1, host.scrollWidth) : 0.5;
    seqMapState.zoom = Math.max(1, Math.min(16, seqMapState.zoom * factor));
    renderSeqMap();
    if (host) host.scrollLeft = Math.max(0, before * host.scrollWidth - host.clientWidth / 2);
  }
  var zin = document.getElementById('mapZoomIn');
  var zout = document.getElementById('mapZoomOut');
  var fit = document.getElementById('mapFit');
  if (zin) zin.addEventListener('click', function () { zoomBy(1.6); });
  if (zout) zout.addEventListener('click', function () { zoomBy(1 / 1.6); });
  if (fit) fit.addEventListener('click', function () { seqMapState.zoom = 1; renderSeqMap(); });

  /* 线性 / 环状 */
  var lin = document.getElementById('mapLinear');
  var cir = document.getElementById('mapCircular');
  function setMapMode(isCircular) {
    seqMapState.circular = !!isCircular;
    seqMapState.sel = null;
    if (lin) lin.classList.toggle('active', !isCircular);
    if (cir) cir.classList.toggle('active', !!isCircular);
    renderSeqMap();
  }
  if (lin) lin.addEventListener('click', function () { setMapMode(false); });
  if (cir) cir.addEventListener('click', function () { setMapMode(true); });

  /* 注释：添加 / 导入 / 列表 */
  var addBtn = document.getElementById('mapAddFeature');
  if (addBtn) addBtn.addEventListener('click', addMapFeature);
  var expBtn = document.getElementById('mapExportFeatures');
  if (expBtn) expBtn.addEventListener('click', exportFeatureCSV);
  var impBtn = document.getElementById('mapImportGB');
  var fileInput = document.getElementById('mapGBFile');
  if (impBtn && fileInput) {
    impBtn.addEventListener('click', function () { fileInput.click(); });
    fileInput.addEventListener('change', function () {
      if (this.files && this.files[0]) importGbFile(this.files[0]);
      this.value = '';
    });
  }
  var fbox = document.getElementById('mapFeatures');
  if (fbox) fbox.addEventListener('click', function (e) {
    var del = e.target.closest ? e.target.closest('[data-mf-del]') : null;
    if (del) {
      var i = parseInt(del.getAttribute('data-mf-del'), 10);
      if (!isNaN(i)) {
        mapFeatures.splice(i, 1);
        renderSeqMap();
        renderMapFeatures();
        if (typeof scheduleWorkspaceSave === 'function') scheduleWorkspaceSave();
      }
      return;
    }
    var chip = e.target.closest ? e.target.closest('[data-mf-sel]') : null;
    if (chip) {
      var j = parseInt(chip.getAttribute('data-mf-sel'), 10);
      var f = mapFeatures[j];
      if (f) {
        seqMapState.sel = { a: Math.min(f.start, f.end) - 1, b: Math.max(f.start, f.end) };
        renderSeqMap();
      }
    }
  });

  /* 轨道开关 */
  var chips = document.querySelectorAll('.map-chip');
  for (var ci = 0; ci < chips.length; ci++) {
    (function (chip) {
      chip.addEventListener('click', function () {
        var k = chip.getAttribute('data-track');
        if (!k) return;
        seqMapState.tracks[k] = !seqMapState.tracks[k];
        chip.setAttribute('aria-pressed', seqMapState.tracks[k] ? 'true' : 'false');
        renderSeqMap();
      });
    })(chips[ci]);
  }

  var rt = 0;
  window.addEventListener('resize', function () {
    clearTimeout(rt);
    rt = setTimeout(function () { renderSeqMap(); }, 180);
  });
}


/* =========================================================
   十五、工作区：本地保存/恢复 + 文件上传
   ========================================================= */
