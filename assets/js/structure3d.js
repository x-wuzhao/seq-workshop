/* =====================================================================
   序列工坊 · 三维结构查看器（structure3d.js）
   ---------------------------------------------------------------------
   目标：把「DNA 空间结构」从平面 SVG 投影升级为真正的 WebGL 分子查看器：
     1) 用户序列 → 按 B-DNA 标准几何参数实时生成双螺旋（10.5 bp/圈、
        螺距 3.38 Å、直径约 20 Å），可见大沟/小沟、5′→3′ 方向、碱基对；
     2) 内置真实晶体/核磁结构（RCSB PDB 官方坐标）：B-DNA 十二聚体、
        tRNA、Trp-cage、GFP、血红蛋白，支持卡通/球棍、按链/按二级结构着色；
     3) 序列与结构联动：点击碱基对或残基回传坐标，面板显示注释。

   依赖：three.js（MIT，随包提供）。不联网、不上传任何序列。
   ===================================================================== */
(function () {
  'use strict';

  var THREE = window.THREE;
  if (!THREE) { window.Struct3D = null; return; }

  var BASE_COLORS = { A: 0x2E8B4F, T: 0xD64545, G: 0xD98A16, C: 0x2F7DC1, U: 0x7A4FA3, N: 0x7B8794 };
  var ELEM_COLORS = { C: 0xB9C6C4, N: 0x2F6FD0, O: 0xD64545, P: 0xD98A16, S: 0xD9A44A, FE: 0xB4531A, MG: 0x2FA36B, ZN: 0x8E7CC3 };
  var SS_COLORS = { H: 0x0E7C7B, E: 0xB9781B, C: 0x7C8B89 };

  var RISE = 3.38;                    // B-DNA 每碱基对上升（Å）
  var TWIST = Math.PI * 2 / 10.5;     // 每碱基对旋转角（10.5 bp/圈）
  var R_BACKBONE = 9.4;               // 磷酸骨架半径（Å）

  /* ------------------------------------------------------------------
     基础工具
     ------------------------------------------------------------------ */
  function hex(n) { return n; }

  function makeLabel(text, color, scale) {
    var c = document.createElement('canvas');
    c.width = 256; c.height = 128;
    var ctx = c.getContext('2d');
    ctx.clearRect(0, 0, 256, 128);
    ctx.font = '700 72px "JetBrains Mono", Consolas, monospace';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillStyle = color || '#FFFFFF';
    ctx.fillText(text, 128, 66);
    var tex = new THREE.CanvasTexture(c);
    tex.minFilter = THREE.LinearFilter;
    var mat = new THREE.SpriteMaterial({ map: tex, transparent: true, depthWrite: false });
    var sp = new THREE.Sprite(mat);
    var s = scale || 5;
    sp.scale.set(s, s / 2, 1);
    return sp;
  }

  /* 用 Catmull-Rom 曲线生成带状（卡通）几何体：截面为扁矩形，比圆管更像
     专业查看器里的 ribbon 表示。 */
  function ribbonGeometry(points, width, thick) {
    var curve = new THREE.CatmullRomCurve3(points, false, 'catmullrom', 0.5);
    var seg = Math.max(16, Math.min(420, points.length * 8));
    var frames = curve.computeFrenetFrames(seg, false);
    var pos = [], idx = [];
    var corners = [[-width / 2, -thick / 2], [width / 2, -thick / 2], [width / 2, thick / 2], [-width / 2, thick / 2]];
    for (var i = 0; i <= seg; i++) {
      var t = i / seg;
      var p = curve.getPointAt(t);
      var n = frames.normals[i], b = frames.binormals[i];
      for (var k = 0; k < 4; k++) {
        var a = corners[k][0], c = corners[k][1];
        pos.push(
          p.x + n.x * a + b.x * c,
          p.y + n.y * a + b.y * c,
          p.z + n.z * a + b.z * c
        );
      }
    }
    for (var j = 0; j < seg; j++) {
      var o0 = j * 4, o1 = (j + 1) * 4;
      for (var f = 0; f < 4; f++) {
        var g = (f + 1) % 4;
        idx.push(o0 + f, o1 + f, o1 + g, o0 + f, o1 + g, o0 + g);
      }
    }
    var geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    geo.setIndex(idx);
    geo.computeVertexNormals();
    return geo;
  }

  /* 键长判定：共价半径之和 + 容差 */
  var COV = { C: 0.76, N: 0.71, O: 0.66, P: 1.07, S: 1.05, FE: 1.32, MG: 1.41, ZN: 1.22, H: 0.31 };
  function covR(e) { return COV[e] || 0.8; }

  /* 空间网格加速近邻搜索，避免 O(n²) */
  function buildBonds(atoms, positions) {
    var cell = 3.2, grid = {}, bonds = [];
    var n = atoms.length;
    for (var i = 0; i < n; i++) {
      var p = positions[i];
      var key = Math.floor(p[0] / cell) + ',' + Math.floor(p[1] / cell) + ',' + Math.floor(p[2] / cell);
      (grid[key] || (grid[key] = [])).push(i);
    }
    var seen = {};
    for (var a = 0; a < n; a++) {
      var pa = positions[a];
      var gx = Math.floor(pa[0] / cell), gy = Math.floor(pa[1] / cell), gz = Math.floor(pa[2] / cell);
      for (var dx = -1; dx <= 1; dx++) for (var dy = -1; dy <= 1; dy++) for (var dz = -1; dz <= 1; dz++) {
        var arr = grid[(gx + dx) + ',' + (gy + dy) + ',' + (gz + dz)];
        if (!arr) continue;
        for (var t = 0; t < arr.length; t++) {
          var b = arr[t];
          if (b <= a) continue;
          var pb = positions[b];
          var ddx = pa[0] - pb[0], ddy = pa[1] - pb[1], ddz = pa[2] - pb[2];
          var d2 = ddx * ddx + ddy * ddy + ddz * ddz;
          if (d2 < 0.16 || d2 > 6.25) continue;
          var lim = covR(atoms[a].el) + covR(atoms[b].el) + 0.45;
          if (d2 <= lim * lim) {
            var kk = a + '_' + b;
            if (!seen[kk]) { seen[kk] = 1; bonds.push([a, b]); }
          }
        }
      }
    }
    return bonds;
  }

  /* InstancedMesh 的包围盒按实例矩阵计算（three 的 Box3 不展开实例），
     因此这里统一用坐标点直接求包围盒，保证相机取景正确。 */
  function boundsOf(points) {
    var minx = Infinity, miny = Infinity, minz = Infinity;
    var maxx = -Infinity, maxy = -Infinity, maxz = -Infinity;
    for (var i = 0; i < points.length; i++) {
      var p = points[i];
      if (p[0] < minx) minx = p[0];
      if (p[1] < miny) miny = p[1];
      if (p[2] < minz) minz = p[2];
      if (p[0] > maxx) maxx = p[0];
      if (p[1] > maxy) maxy = p[1];
      if (p[2] > maxz) maxz = p[2];
    }
    return {
      center: new THREE.Vector3((minx + maxx) / 2, (miny + maxy) / 2, (minz + maxz) / 2),
      size: new THREE.Vector3(maxx - minx, maxy - miny, maxz - minz)
    };
  }

  /* ------------------------------------------------------------------
     查看器
     ------------------------------------------------------------------ */
  function Viewer(host, opts) {
    this.host = host;
    this.opts = opts || {};
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(this.opts.dark ? 0x0A1312 : 0xF7FBFA);
    this.camera = new THREE.PerspectiveCamera(42, 1, 0.1, 4000);
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'high-performance' });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.domElement.className = 's3d-canvas';
    host.appendChild(this.renderer.domElement);

    var amb = new THREE.HemisphereLight(0xffffff, 0x8899aa, 0.85);
    this.scene.add(amb);
    var d1 = new THREE.DirectionalLight(0xffffff, 1.05); d1.position.set(60, 120, 140); this.scene.add(d1);
    var d2 = new THREE.DirectionalLight(0xbfe0dd, 0.4); d2.position.set(-120, -60, -80); this.scene.add(d2);

    this.root = new THREE.Group();
    this.scene.add(this.root);
    this.modelGroup = new THREE.Group();
    this.root.add(this.modelGroup);
    this._baseRotX = 0;

    this.fitDistance = 100;
    this.distance = 100;
    this.autoRotate = false;
    this.dirty = true;
    this.pickables = [];
    this.onPick = null;
    this.modelKind = null;
    this.currentId = null;

    this._bindInput();
    this._resizeObserver();
    this._tick();
  }

  Viewer.prototype._resizeObserver = function () {
    var self = this;
    function size() {
      var w = self.host.clientWidth || 600;
      var h = self.host.clientHeight || 420;
      self.renderer.setSize(w, h, true);
      self.camera.aspect = w / Math.max(1, h);
      self.camera.updateProjectionMatrix();
      self.dirty = true;
    }
    size();
    if (window.ResizeObserver) {
      this._ro = new ResizeObserver(size);
      this._ro.observe(this.host);
    } else {
      window.addEventListener('resize', size);
    }
  };

  Viewer.prototype._bindInput = function () {
    var self = this;
    var el = this.renderer.domElement;
    var drag = null;
    el.style.touchAction = 'none';
    el.addEventListener('pointerdown', function (e) {
      drag = { x: e.clientX, y: e.clientY, moved: 0, t: Date.now() };
      el.setPointerCapture(e.pointerId);
    });
    el.addEventListener('pointermove', function (e) {
      if (!drag) return;
      var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      drag.moved += Math.abs(dx) + Math.abs(dy);
      drag.x = e.clientX; drag.y = e.clientY;
      self.root.rotation.y += dx * 0.008;
      self.root.rotation.x += dy * 0.008;
      self.root.rotation.x = Math.max(-1.45, Math.min(1.45, self.root.rotation.x));
      self.dirty = true;
    });
    function up(e) {
      if (!drag) return;
      var wasClick = drag.moved < 6 && Date.now() - drag.t < 500;
      drag = null;
      if (wasClick) self._pick(e);
      try { el.releasePointerCapture(e.pointerId); } catch (err) { }
    }
    el.addEventListener('pointerup', up);
    el.addEventListener('pointercancel', function () { drag = null; });
    el.addEventListener('wheel', function (e) {
      e.preventDefault();
      self.zoom(e.deltaY > 0 ? 1.1 : 0.9);
    }, { passive: false });

    var pinch = null;
    el.addEventListener('touchstart', function (e) {
      if (e.touches.length === 2) {
        pinch = { d: dist(e.touches[0], e.touches[1]) };
      }
    }, { passive: true });
    el.addEventListener('touchmove', function (e) {
      if (e.touches.length === 2 && pinch) {
        var d = dist(e.touches[0], e.touches[1]);
        self.zoom(pinch.d / d);
        pinch.d = d;
      }
    }, { passive: true });
    el.addEventListener('touchend', function () { pinch = null; });

    function dist(a, b) {
      var dx = a.clientX - b.clientX, dy = a.clientY - b.clientY;
      return Math.sqrt(dx * dx + dy * dy);
    }
  };

  Viewer.prototype.zoom = function (factor) {
    this.distance = Math.max(this.fitDistance * 0.35, Math.min(this.fitDistance * 3.2, this.distance * factor));
    this.dirty = true;
  };

  Viewer.prototype.setAutoRotate = function (on) {
    this.autoRotate = !!on;
    this.dirty = true;
  };

  Viewer.prototype.resetView = function () {
    this.root.rotation.set(0, 0, 0);
    this.modelGroup.rotation.set(this._baseRotX || 0, 0, 0);
    this.distance = this.fitDistance;
    this.dirty = true;
  };

  Viewer.prototype.setTheme = function (dark) {
    this.scene.background = new THREE.Color(dark ? 0x0A1312 : 0xF7FBFA);
    this.dirty = true;
  };

  Viewer.prototype.clear = function () {
    var self = this;
    this.root.traverse(function (o) {
      if (o.geometry) o.geometry.dispose();
      if (o.material) {
        if (o.material.map) o.material.map.dispose();
        o.material.dispose();
      }
    });
    while (this.modelGroup.children.length) this.modelGroup.remove(this.modelGroup.children[0]);
    this.modelGroup.position.set(0, 0, 0);
    this.modelGroup.rotation.set(0, 0, 0);
    this._baseRotX = 0;
    this.pickables = [];
    this.modelKind = null;
    this.currentId = null;
    this._bounds = null;
    this.dirty = true;
  };

  Viewer.prototype._fit = function () {
    this.root.position.set(0, 0, 0);
    this.modelGroup.position.set(0, 0, 0);
    this.modelGroup.rotation.set(0, 0, 0);
    this.root.updateMatrixWorld(true);

    var center, size;
    if (this._bounds) {
      center = this._bounds.center.clone();
      size = this._bounds.size.clone();
    } else {
      var box = new THREE.Box3().setFromObject(this.modelGroup);
      if (box.isEmpty()) return;
      center = box.getCenter(new THREE.Vector3());
      size = box.getSize(new THREE.Vector3());
    }

    // 长轴若沿 Z 轴（如多数 PDB 核酸结构），把它转到竖直方向，默认视角更直观
    var rotX = 0;
    if (size.z > Math.max(size.x, size.y) * 1.15) {
      rotX = -Math.PI / 2;
    } else if (size.x > Math.max(size.y, size.z) * 1.15) {
      rotX = 0;
      this.modelGroup.rotation.set(0, Math.PI / 2, 0);
    }
    this._baseRotX = rotX;
    this.modelGroup.rotation.set(rotX, this.modelGroup.rotation.y, 0);
    // 旋转后重新居中：先转中心，再取反作为平移
    var cWorld = center.clone().applyEuler(new THREE.Euler(rotX, this.modelGroup.rotation.y, 0));
    this.modelGroup.position.copy(cWorld).multiplyScalar(-1);

    var ey = rotX ? size.z : size.y;
    var ez = rotX ? size.y : size.z;
    var ex = size.x;

    // 把包围盒 8 个角投影到视锥里，求出「刚好装下」的相机距离（留 12% 边距）
    var vFov = THREE.MathUtils.degToRad(this.camera.fov);
    var aspect = Math.max(0.2, this.camera.aspect || 1);
    var tanV = Math.tan(vFov / 2), tanH = tanV * aspect;
    var hx = Math.max(ex, 6) / 2, hy = Math.max(ey, 6) / 2, hz = Math.max(ez, 6) / 2;
    var limit = 0.88, need = 0;
    for (var sx = -1; sx <= 1; sx += 2) {
      for (var sy = -1; sy <= 1; sy += 2) {
        for (var sz = -1; sz <= 1; sz += 2) {
          need = Math.max(need, sz * hz + Math.abs(sy * hy) / (limit * tanV));
          need = Math.max(need, sz * hz + Math.abs(sx * hx) / (limit * tanH));
        }
      }
    }
    this.fitDistance = need;
    this.distance = this.fitDistance;
    this.dirty = true;
  };

  Viewer.prototype._tick = function () {
    var self = this;
    function loop() {
      if (self.autoRotate) {
        self.root.rotation.y += 0.006;
        self.dirty = true;
      }
      if (self.dirty) {
        var d = self.distance;
        self.camera.position.set(0, d * 0.07, d);
        self.camera.lookAt(0, 0, 0);
        self.renderer.render(self.scene, self.camera);
        self.dirty = false;
      }
      self._raf = requestAnimationFrame(loop);
    }
    loop();
  };

  Viewer.prototype._pick = function (e) {
    if (!this.pickables.length) return;
    var rect = this.renderer.domElement.getBoundingClientRect();
    var ndc = new THREE.Vector2(
      ((e.clientX - rect.left) / rect.width) * 2 - 1,
      -((e.clientY - rect.top) / rect.height) * 2 + 1
    );
    var ray = new THREE.Raycaster();
    ray.setFromCamera(ndc, this.camera);
    var hits = ray.intersectObjects(this.pickables, false);
    if (hits.length && this.onPick) this.onPick(hits[0].object.userData, hits[0]);
  };

  /* ------------------------------------------------------------------
     B-DNA：按标准几何参数实时生成
     ------------------------------------------------------------------ */
  Viewer.prototype.showDNA = function (seq, center, winSize, bare) {
    this.clear();
    this.modelKind = 'dna';
    var n = seq.length;
    if (!n) return;
    var win = Math.min(n, winSize || 24);
    var start = Math.max(0, Math.min(n - win, (center || Math.ceil(n / 2)) - 1 - Math.floor(win / 2)));
    var sub = seq.substring(start, start + win);

    var strand1 = [], strand2 = [];
    var plates = [];
    var i, th, y, p1, p2;

    for (i = 0; i < sub.length; i++) {
      th = i * TWIST;
      y = (i - (sub.length - 1) / 2) * RISE;
      p1 = new THREE.Vector3(Math.cos(th) * R_BACKBONE, y, Math.sin(th) * R_BACKBONE);
      p2 = new THREE.Vector3(Math.cos(th + Math.PI) * R_BACKBONE, y, Math.sin(th + Math.PI) * R_BACKBONE);
      strand1.push(p1); strand2.push(p2);
    }

    var matBack1 = new THREE.MeshStandardMaterial({ color: 0x0E7C7B, roughness: 0.42, metalness: 0.05 });
    var matBack2 = new THREE.MeshStandardMaterial({ color: 0xB9781B, roughness: 0.42, metalness: 0.05 });
    if (sub.length >= 2) {
      var c1 = new THREE.CatmullRomCurve3(strand1);
      var c2 = new THREE.CatmullRomCurve3(strand2);
      var tube1 = new THREE.Mesh(new THREE.TubeGeometry(c1, Math.max(24, sub.length * 10), 1.35, 12, false), matBack1);
      var tube2 = new THREE.Mesh(new THREE.TubeGeometry(c2, Math.max(24, sub.length * 10), 1.35, 12, false), matBack2);
      this.modelGroup.add(tube1, tube2);

      // 磷酸基团小球，让骨架更像真实的糖-磷酸链
      var sphGeo = new THREE.SphereGeometry(0.78, 10, 8);
      var count = strand1.length * 2;
      var phos = new THREE.InstancedMesh(sphGeo, new THREE.MeshStandardMaterial({ color: 0xffffff, roughness: 0.5 }), count);
      var pm = new THREE.Matrix4(), pc = new THREE.Color();
      var k = 0;
      for (var s = 0; s < strand1.length; s++) {
        pm.makeTranslation(strand1[s].x, strand1[s].y, strand1[s].z);
        phos.setMatrixAt(k, pm); pc.setHex(0x0E7C7B); phos.setColorAt(k, pc); k++;
        pm.makeTranslation(strand2[s].x, strand2[s].y, strand2[s].z);
        phos.setMatrixAt(k, pm); pc.setHex(0xB9781B); phos.setColorAt(k, pc); k++;
      }
      phos.instanceMatrix.needsUpdate = true;
      if (phos.instanceColor) phos.instanceColor.needsUpdate = true;
      this.modelGroup.add(phos);
    }

    // 碱基对：两块扁平碱基 + 氢键；中间留 ~2.8 Å 空隙画氢键
    var geomCache = {};
    var focusIdx = (center || 1) - 1 - start;
    for (i = 0; i < sub.length && !bare; i++) {
      var base = sub.charAt(i);
      var comp = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' }[base] || 'N';
      th = i * TWIST;
      y = strand1[i].y;
      var ux = Math.cos(th), uz = Math.sin(th);

      for (var side = 0; side < 2; side++) {
        var bb = side === 0 ? base : comp;
        var dir = side === 0 ? 1 : -1;
        var geo = geomCache[bb] || (geomCache[bb] = new THREE.BoxGeometry(3.1, 0.55, 7.0));
        var mat = new THREE.MeshStandardMaterial({
          color: BASE_COLORS[bb] || BASE_COLORS.N,
          roughness: 0.55, metalness: 0.0,
          emissive: (i === focusIdx ? 0x3a2a08 : 0x000000)
        });
        var m = new THREE.Mesh(geo, mat);
        m.position.set(ux * dir * 5.0, y, uz * dir * 5.0);
        m.lookAt(new THREE.Vector3(0, y, 0));
        m.userData = { kind: 'bp', index: start + i, base: base, comp: comp };
        this.modelGroup.add(m);
        plates.push(m);
        this.pickables.push(m);
      }

      // 氢键（A–T 2 条，G–C 3 条），横跨中间空隙
      var nH = (base === 'G' || base === 'C') ? 3 : 2;
      for (var h = 0; h < nH; h++) {
        var off = (h - (nH - 1) / 2) * 1.15;
        var hb = new THREE.Mesh(
          new THREE.CylinderGeometry(0.14, 0.14, 2.9, 6),
          new THREE.MeshStandardMaterial({ color: 0xEAF3F2, roughness: 0.6, transparent: true, opacity: 0.8 })
        );
        var tx = -uz, tz = ux;
        hb.position.set(tx * off, y, tz * off);
        hb.rotation.z = Math.PI / 2;
        hb.rotation.y = -th;
        this.modelGroup.add(hb);
      }
    }

    // 聚焦碱基对的高亮圆环
    if (!bare && focusIdx >= 0 && focusIdx < sub.length) {
      var ring = new THREE.Mesh(
        new THREE.TorusGeometry(R_BACKBONE + 2.0, 0.35, 8, 64),
        new THREE.MeshBasicMaterial({ color: 0xB9781B, transparent: true, opacity: 0.75 })
      );
      ring.rotation.x = Math.PI / 2;
      ring.position.y = strand1[focusIdx].y;
      this.modelGroup.add(ring);
    }

    // 5′/3′ 方向标注
    if (strand1.length) {
      var last = strand1.length - 1;
      var l1 = makeLabel('5′', '#0E7C7B', 11); l1.position.copy(strand1[0]).multiplyScalar(1.3); l1.position.y += 2.2;
      var l2 = makeLabel('3′', '#0E7C7B', 11); l2.position.copy(strand1[last]).multiplyScalar(1.3); l2.position.y -= 2.2;
      var l3 = makeLabel('3′', '#B9781B', 11); l3.position.copy(strand2[0]).multiplyScalar(1.3); l3.position.y += 2.2;
      var l4 = makeLabel('5′', '#B9781B', 11); l4.position.copy(strand2[last]).multiplyScalar(1.3); l4.position.y -= 2.2;
      this.modelGroup.add(l1, l2, l3, l4);
    }

    this._fit();
    this.currentId = null;
    this.dnaRange = { start: start + 1, end: start + sub.length };
  };

  /* ------------------------------------------------------------------
     PDB 结构：真实坐标
     ------------------------------------------------------------------ */
  Viewer.prototype.showStructure = function (id, style, colorMode) {
    var data = window.STRUCTURES && window.STRUCTURES[id];
    if (!data) return false;
    this.clear();
    this.modelKind = data.kind;
    this.currentId = id;
    this._data = data;

    var kind = data.kind;
    if (data.atoms) this._buildAtoms(data, style || 'stick', colorMode);
    if (data.ca) this._buildCartoon(data, style, colorMode);
    if (data.lig) this._buildLigands(data);

    this._fit();
    return true;
  };

  Viewer.prototype._atomColor = function (data, a, colorMode) {
    var res = data.res[a[3]];
    var resName = res[0], aa = res[1];
    var atomName = data.atab[a[4]];
    var elem = data.etab[a[5]];
    if (data.kind === 'dna' || data.kind === 'rna') {
      var backbone = /^(P|OP1|OP2|O5'|C5'|C4'|O4'|C3'|O3'|C2'|C1'|O2')$/;
      if (!backbone.test(atomName)) {
        return BASE_COLORS[aa] !== undefined ? BASE_COLORS[aa] : BASE_COLORS.N;
      }
      return elem === 'P' ? 0xD98A16 : 0x93A7A4;
    }
    if (colorMode === 'chain') {
      var idx = data.chains.indexOf(res[2]);
      return [0x0E7C7B, 0xB9781B, 0x2F7DC1, 0x2E8B4F, 0xD64545][idx % 5];
    }
    return ELEM_COLORS[elem] !== undefined ? ELEM_COLORS[elem] : 0x9AA6A4;
  };

  Viewer.prototype._buildAtoms = function (data, style, colorMode) {
    var atoms = data.atoms;
    var pos = new Array(atoms.length), i;
    var xs = 0, ys = 0, zs = 0;
    for (i = 0; i < atoms.length; i++) {
      pos[i] = [atoms[i][0], atoms[i][1], atoms[i][2]];
      xs += atoms[i][0]; ys += atoms[i][1]; zs += atoms[i][2];
    }
    xs /= atoms.length; ys /= atoms.length; zs /= atoms.length;
    for (i = 0; i < atoms.length; i++) {
      pos[i][0] -= xs; pos[i][1] -= ys; pos[i][2] -= zs;
    }
    this._bounds = boundsOf(pos);

    var showSticks = style !== 'sphere';
    var showSpheres = style !== 'stick';

    if (showSticks) {
      var bonds = buildBonds(atoms.map(function (a) {
        return { el: data.etab[a[5]] };
      }), pos);
      var cyl = new THREE.InstancedMesh(
        new THREE.CylinderGeometry(1, 1, 1, 6, 1, true),
        new THREE.MeshStandardMaterial({ color: 0xffffff, roughness: 0.45, metalness: 0.05 }),
        bonds.length * 2
      );
      var m4 = new THREE.Matrix4(), q = new THREE.Quaternion(), up = new THREE.Vector3(0, 1, 0);
      var pA = new THREE.Vector3(), pB = new THREE.Vector3(), mid = new THREE.Vector3(), dir = new THREE.Vector3();
      var col = new THREE.Color();
      var k = 0;
      for (i = 0; i < bonds.length; i++) {
        var a = bonds[i][0], b = bonds[i][1];
        pA.fromArray(pos[a]); pB.fromArray(pos[b]);
        dir.subVectors(pB, pA);
        var len = dir.length();
        q.setFromUnitVectors(up, dir.clone().normalize());
        mid.addVectors(pA, pB).multiplyScalar(0.5);
        m4.compose(mid, q, new THREE.Vector3(1, 1, 1));
        var radius = 0.17;
        // 半段各取一端颜色，模拟专业查看器的双色键
        [a, b].forEach(function (atomIdx, half) {
          var m2 = m4.clone();
          var p = half === 0 ? pA : pB;
          var c = half === 0 ? mid.clone().add(pA).multiplyScalar(0.5) : mid.clone().add(pB).multiplyScalar(0.5);
          m2.makeTranslation(0, 0, 0);
          q.setFromUnitVectors(up, dir.clone().normalize());
          var sc = new THREE.Vector3(radius, len / 2, radius);
          m2.compose(c, q, sc);
          cyl.setMatrixAt(k, m2);
          col.setHex(this._atomColor(data, atoms[atomIdx], colorMode));
          cyl.setColorAt(k, col);
          k++;
        }, this);
      }
      cyl.count = k;
      cyl.instanceMatrix.needsUpdate = true;
      if (cyl.instanceColor) cyl.instanceColor.needsUpdate = true;
      this.modelGroup.add(cyl);
    }

    if (showSpheres) {
      var sph = new THREE.InstancedMesh(
        new THREE.SphereGeometry(1, 12, 10),
        new THREE.MeshStandardMaterial({ color: 0xffffff, roughness: 0.4, metalness: 0.05 }),
        atoms.length
      );
      var mm = new THREE.Matrix4(), cc = new THREE.Color();
      for (i = 0; i < atoms.length; i++) {
        var elem = data.etab[atoms[i][5]];
        var rad = elem === 'P' ? 0.62 : (elem === 'C' ? 0.32 : 0.38);
        mm.makeTranslation(pos[i][0], pos[i][1], pos[i][2]);
        mm.scale(new THREE.Vector3(rad, rad, rad));
        sph.setMatrixAt(i, mm);
        cc.setHex(this._atomColor(data, atoms[i], colorMode));
        sph.setColorAt(i, cc);
      }
      sph.instanceMatrix.needsUpdate = true;
      if (sph.instanceColor) sph.instanceColor.needsUpdate = true;
      this.modelGroup.add(sph);
    }
  };

  Viewer.prototype._buildLigands = function (data) {
    var atoms = data.lig;
    var pos = new Array(atoms.length), i;
    for (i = 0; i < atoms.length; i++) pos[i] = [atoms[i][0], atoms[i][1], atoms[i][2]];
    var bonds = buildBonds(atoms.map(function (a) { return { el: data.etab[a[5]] }; }), pos);
    var cyl = new THREE.InstancedMesh(
      new THREE.CylinderGeometry(1, 1, 1, 8),
      new THREE.MeshStandardMaterial({ color: 0xC0392B, roughness: 0.35, metalness: 0.2 }),
      bonds.length
    );
    var m4 = new THREE.Matrix4(), q = new THREE.Quaternion(), up = new THREE.Vector3(0, 1, 0);
    var pA = new THREE.Vector3(), pB = new THREE.Vector3(), mid = new THREE.Vector3(), dir = new THREE.Vector3();
    for (i = 0; i < bonds.length; i++) {
      pA.fromArray(pos[bonds[i][0]]); pB.fromArray(pos[bonds[i][1]]);
      dir.subVectors(pB, pA);
      var len = dir.length();
      q.setFromUnitVectors(up, dir.clone().normalize());
      mid.addVectors(pA, pB).multiplyScalar(0.5);
      m4.compose(mid, q, new THREE.Vector3(0.19, len, 0.19));
      cyl.setMatrixAt(i, m4);
    }
    cyl.instanceMatrix.needsUpdate = true;
    this.modelGroup.add(cyl);
  };

  Viewer.prototype._buildCartoon = function (data, style, colorMode) {
    var ca = data.ca;
    if (!ca.length) return;
    // 按链、按残基序号排序
    var chains = {};
    for (var i = 0; i < ca.length; i++) {
      var res = data.res[ca[i][3]];
      var ch = res[2];
      (chains[ch] || (chains[ch] = [])).push({ p: new THREE.Vector3(ca[i][0], ca[i][1], ca[i][2]), seq: res[3], name: res[0], ri: ca[i][3] });
    }
    var chainIds = Object.keys(chains).sort();
    var self = this;
    var caForPick = [];

    if (style === 'stick') {
      // 卡通模式下也允许看主链：画 Cα 连线
      chainIds.forEach(function (ch) {
        var pts = chains[ch].map(function (r) { return r.p; });
        if (pts.length < 2) return;
        var curve = new THREE.CatmullRomCurve3(pts);
        var geo = new THREE.TubeGeometry(curve, pts.length * 6, 0.6, 8, false);
        var color = colorMode === 'chain'
          ? [0x0E7C7B, 0xB9781B, 0x2F7DC1, 0x2E8B4F, 0xD64545][chainIds.indexOf(ch) % 5]
          : 0x7C8B89;
        self.modelGroup.add(new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ color: color, roughness: 0.4 })));
      });
    } else {
      chainIds.forEach(function (ch) {
        var rows = chains[ch];
        var pts = rows.map(function (r) { return r.p; });
        if (pts.length < 3) return;
        // 按二级结构切成连续片段，分别着色
        var runs = [], cur = null;
        for (var j = 0; j < rows.length; j++) {
          var ss = data.ssmap[ch + ':' + rows[j].seq] || 'C';
          if (!cur || cur.ss !== ss) { cur = { ss: ss, pts: [] }; runs.push(cur); }
          cur.pts.push(pts[j]);
        }
        runs.forEach(function (run, ri) {
          var seg = run.pts;
          if (seg.length < 2) {
            seg = [seg[0], seg[0].clone().add(new THREE.Vector3(0.01, 0.01, 0))];
          }
          // 片段之间做一点过渡，避免断裂
          if (ri > 0) seg = [runs[ri - 1].pts[runs[ri - 1].pts.length - 1]].concat(seg);
          var color = colorMode === 'chain'
            ? [0x0E7C7B, 0xB9781B, 0x2F7DC1, 0x2E8B4F, 0xD64545][chainIds.indexOf(ch) % 5]
            : (SS_COLORS[run.ss] || SS_COLORS.C);
          var mat = new THREE.MeshStandardMaterial({ color: color, roughness: 0.38, metalness: 0.04, side: THREE.DoubleSide });
          if (run.ss === 'C') {
            // 无规卷曲用细圆管，避免与 β 折叠带子糊成一团
            var curveC = new THREE.CatmullRomCurve3(seg);
            var geoC = new THREE.TubeGeometry(curveC, Math.max(8, seg.length * 6), 0.45, 8, false);
            self.modelGroup.add(new THREE.Mesh(geoC, mat));
          } else {
            var w = run.ss === 'H' ? 2.4 : 1.9;
            var t = run.ss === 'E' ? 0.26 : 0.5;
            self.modelGroup.add(new THREE.Mesh(ribbonGeometry(seg, w, t), mat));
          }
        });
      });
    }

    // 透明拾取球：点击残基时返回注释
    var pick = new THREE.InstancedMesh(
      new THREE.SphereGeometry(0.9, 6, 5),
      new THREE.MeshBasicMaterial({ transparent: true, opacity: 0, depthWrite: false }),
      ca.length
    );
    var m = new THREE.Matrix4();
    for (var k = 0; k < ca.length; k++) {
      var r = data.res[ca[k][3]];
      m.makeTranslation(ca[k][0], ca[k][1], ca[k][2]);
      pick.setMatrixAt(k, m);
      caForPick.push(r);
    }
    pick.instanceMatrix.needsUpdate = true;
    pick.userData = { kind: 'residue', residues: caForPick };
    if (!this._bounds) {
      this._bounds = boundsOf(ca.map(function (a) { return [a[0], a[1], a[2]]; }));
    }
    this.modelGroup.add(pick);
    this.pickables.push(pick);
    this._pickInstanced = pick;
  };

  Viewer.prototype.showModel = function (mode, arg, style, colorMode) {
    if (mode === 'dna') this.showDNA(arg.seq, arg.center, arg.win);
    else this.showStructure(arg, style, colorMode);
  };

  window.Struct3D = {
    Viewer: Viewer,
    create: function (host, opts) {
      try {
        return new Viewer(host, opts);
      } catch (e) {
        console.warn('WebGL 初始化失败：', e);
        return null;
      }
    },
    baseColors: BASE_COLORS
  };
})();
