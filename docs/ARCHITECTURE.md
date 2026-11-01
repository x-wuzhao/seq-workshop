# 序列工坊系统结构

## 目标

序列工坊是一个面向移动端的 DNA 序列分析 Web/PWA。用户在浏览器中输入序列，所有分析在本地完成，再以可视化结果、教学解释和可接续的格式（FASTA / CSV / GenBank）导出；三维结构使用随包提供的真实 PDB 坐标渲染。

## 数据流

```text
输入文本 / FASTA
       |
       v
enhancedParseInput
       |
       v
核心算法层
  ├─ reverseComplement / transcribe
  ├─ translateFrame / findORFs（可切换遗传密码表）
  ├─ analyzeEnzymes / digestFragments（线性/环状、单/双酶切）
  ├─ nwAlign（全局）/ swAlign（局部）
  ├─ tmNN / pairDG / hairpinDG / designPrimersV2
  ├─ gcWindows / calcProtProps / optimizeCodons
  └─ buildGenBank
       |
       v
结果渲染层
  ├─ 序列与 ORF 展示
  ├─ ECharts GC / 氨基酸图表
  ├─ 酶切图谱、片段表、虚拟凝胶、点阵图
  ├─ 引物候选表、教学动画、突变模拟
  └─ 三维结构（assets/js/structure3d.js + Three.js）
       |
       v
导出层：FASTA / 分析 CSV / 蛋白 CSV / GenBank
```

## 三维结构的实现分层

```text
assets/structures/*.pdb          RCSB 官方坐标（CC0）
        |
        v  tools/pdb2json.py（去水、去氢、CA 或重原子、坐标两位小数）
assets/structures/structures.js  紧凑 JSON（window.STRUCTURES，约 125 KB）
        |
        v
assets/js/structure3d.js
  ├─ Viewer：场景/相机/灯光、自研轨道控制（拖动、滚轮、双指、自动旋转）
  ├─ showDNA：按 B 型标准参数（10.5 bp/圈、3.38 Å、骨架半径 9.4 Å）实时生成双螺旋
  ├─ showStructure：解析 PDB，生成卡通带（HELIX/SHEET）或球棍（原子 + 键）
  ├─ _fit：按包围盒 8 角点投影精确取景；长轴自动转竖直
  └─ 拾取：射线检测碱基对 / Cα 实例，回传坐标与残基注释
```

## 移动与网络应用实现

- `manifest.json` 提供应用名称、主题色、启动方式、快捷入口和安装图标。
- `service-worker.js` 预缓存应用壳、本地 ECharts、Three.js、结构数据与图标；同源静态资源采用 cache-first，页面导航采用 network-first，断网时回退到已缓存的 `index.html`。
- 应用不需要账号或后端 API，序列只在当前浏览器内存中参与计算。
- GitHub Pages 与 localhost 支持 Service Worker；`file://` 仅支持核心页面运行，不支持 PWA 注册。

## 安全与依赖边界

页面不上传序列，不写入序列历史，也不加载运行时 CDN。第三方库随包提供：ECharts 5.5.1（Apache-2.0）、Three.js 0.152.2（MIT）；结构坐标为 RCSB PDB 公共领域数据（CC0）。所有生信算法（反向互补、翻译、ORF、酶切、比对、Tm/ΔG、引物评分、GenBank 生成）与三维几何生成均为本项目自行实现。许可证与来源见 `THIRD-PARTY-NOTICES.txt`。
