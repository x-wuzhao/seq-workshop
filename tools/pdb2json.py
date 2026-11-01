# -*- coding: utf-8 -*-
"""把 RCSB 下载的 PDB 文本压缩成前端可直接引用的紧凑 JSON（写入 structures.js）。

设计取舍：
- 小分子量结构（DNA/RNA/小蛋白）保留全部重原子，供球棍模型使用；
- 大蛋白只保留 Cα 主链 + 二级结构区间，供卡通模型使用，体积从数百 KB 降到数十 KB；
- 坐标保留两位小数（0.01 Å），对教学可视化足够，且显著压缩体积。
"""
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "assets", "structures")
OUT = os.path.join(ROOT, "assets", "structures", "structures.js")

META = {
    "1BNA": {
        "title": "B-DNA 双螺旋（Drew–Dickerson 十二聚体）",
        "subtitle": "DNA 双螺旋的经典晶体结构",
        "method": "X 射线衍射",
        "resolution": 1.9,
        "kind": "dna",
        "mode": "full",
        "note": "序列 CGCGAATTCGCG，2 条链反向平行；可见 B 型构象的右手螺旋、大沟与小沟。",
        "url": "https://www.rcsb.org/structure/1BNA",
    },
    "1EHZ": {
        "title": "转运 RNA（tRNA-Phe）",
        "subtitle": "RNA 三级结构：三叶草折叠成 L 形",
        "method": "X 射线衍射",
        "resolution": 1.93,
        "kind": "rna",
        "mode": "full",
        "note": "76 个核苷酸；反密码子环位于一端，氨基酸臂位于另一端，是“翻译适配器”的结构证据。",
        "url": "https://www.rcsb.org/structure/1EHZ",
    },
    "1L2Y": {
        "title": "Trp-cage 迷你蛋白",
        "subtitle": "20 个氨基酸的微型蛋白（核磁结构）",
        "method": "溶液核磁共振",
        "resolution": None,
        "kind": "protein",
        "mode": "full",
        "note": "仅 20 个残基却形成稳定的疏水核心，常用来演示“序列决定构象”。",
        "url": "https://www.rcsb.org/structure/1L2Y",
    },
    "1GFL": {
        "title": "绿色荧光蛋白 GFP",
        "subtitle": "11 股 β 桶 + 中央 α 螺旋发包团",
        "method": "X 射线衍射",
        "resolution": 2.1,
        "kind": "protein",
        "mode": "ca",
        "note": "238 个氨基酸折叠成 β 桶，发包团被包在桶中央——这是它发光却不被淬灭的结构原因。",
        "url": "https://www.rcsb.org/structure/1GFL",
    },
    "2HHB": {
        "title": "血红蛋白（脱氧型）",
        "subtitle": "α2β2 四聚体 + 4 个血红素",
        "method": "X 射线衍射",
        "resolution": 1.74,
        "kind": "protein",
        "mode": "ca",
        "note": "四个亚基各结合一个血红素（红色球棍），铁原子位于卟啉环中央，负责结合 O₂。",
        "url": "https://www.rcsb.org/structure/2HHB",
    },
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
NA1 = {"DA": "A", "DT": "T", "DG": "G", "DC": "C", "A": "A", "U": "U", "G": "G", "C": "C"}


def parse(path):
    atoms = []       # dicts
    helices = []
    sheets = []
    header = ""
    for line in open(path, encoding="utf-8", errors="replace"):
        rec = line[:6]
        if rec == "HEADER":
            header = line[10:50].strip()
        elif rec == "HELIX ":
            try:
                ch = line[19]
                s = int(line[21:25])
                e = int(line[33:37])
                helices.append((ch, s, e))
            except ValueError:
                pass
        elif rec == "SHEET ":
            try:
                ch = line[21]
                s = int(line[22:26])
                e = int(line[33:37])
                sheets.append((ch, s, e))
            except ValueError:
                pass
        elif rec in ("ATOM  ", "HETATM"):
            name = line[12:16].strip()
            alt = line[16]
            resn = line[17:20].strip()
            ch = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            elem = line[76:78].strip().upper() or re.sub(r"[^A-Za-z]", "", name)[:1].upper()
            if alt not in (" ", "A", ""):
                continue
            if elem in ("H", "D") or resn in ("HOH", "DOD", "WAT"):
                continue
            atoms.append({
                "name": name, "resn": resn, "ch": ch, "resseq": resseq,
                "icode": icode, "x": x, "y": y, "z": z, "elem": elem,
                "het": rec.startswith("HETATM"),
            })
    return atoms, helices, sheets, header


def build(sid):
    m = dict(META[sid])
    path = os.path.join(SRC, sid + ".pdb")
    atoms, helices, sheets, header = parse(path)
    chains = []
    for a in atoms:
        if a["ch"] not in chains:
            chains.append(a["ch"])

    # residue index
    res_index = {}
    order = []
    for a in atoms:
        key = (a["ch"], a["resseq"], a["icode"], a["resn"])
        if key not in res_index:
            res_index[key] = len(order)
            order.append(key)

    atab, rtab, etab = [], [], []
    def tab(lst, v):
        if v not in lst:
            lst.append(v)
        return lst.index(v)

    residues = []
    for key in order:
        ch, resseq, icode, resn = key
        residues.append({
            "name": resn,
            "ch": ch,
            "seq": resseq,
            "aa": AA3.get(resn) or NA1.get(resn) or (resn[:1] if resn not in ("HEM", "SO4") else "X"),
            "het": resn in ("HEM",),
        })

    def pack_atom(a):
        return [
            round(a["x"], 2), round(a["y"], 2), round(a["z"], 2),
            res_index[(a["ch"], a["resseq"], a["icode"], a["resn"])],
            tab(atab, a["name"]),
            tab(etab, a["elem"]),
        ]

    ss = {}
    for ch, s, e in helices:
        for r in range(s, e + 1):
            ss[(ch, r)] = "H"
    for ch, s, e in sheets:
        for r in range(s, e + 1):
            ss[(ch, r)] = "E"

    out = dict(m)
    out["chains"] = chains
    out["atab"] = atab
    out["rtab"] = list(dict.fromkeys([r["name"] for r in residues]))
    out["etab"] = etab
    out["res"] = [[r["name"], r["aa"], r["ch"], r["seq"]] for r in residues]
    out["ssmap"] = {f"{k[0]}:{k[1]}": v for k, v in ss.items()}
    if m["mode"] == "full":
        out["atoms"] = [pack_atom(a) for a in atoms]
        out["kind_atoms"] = "all"
    else:
        ca = []
        for a in atoms:
            if a["name"] != "CA":
                continue
            ri = res_index[(a["ch"], a["resseq"], a["icode"], a["resn"])]
            ca.append([round(a["x"], 2), round(a["y"], 2), round(a["z"], 2), ri])
        out["ca"] = ca
        # 配体（血红素等）保留全部重原子
        lig = [pack_atom(a) for a in atoms if a["het"]]
        if lig:
            out["lig"] = lig
    out["n_atoms"] = len(atoms)
    out["n_res"] = len(residues)
    out["header"] = header
    out["pdb"] = None
    return out


def main():
    data = {}
    for sid in META:
        d = build(sid)
        data[sid] = d
        print(sid, "res", d["n_res"], "atoms", d["n_atoms"],
              "packed", len(d.get("atoms", d.get("ca", []))), "lig", len(d.get("lig", [])))
    js = ("/* 由 tools/pdb2json.py 从 RCSB PDB 官方坐标生成；仅保留重原子/主链与二级结构。\n"
          "   数据来源：https://www.rcsb.org/ （CC0 公共领域）。 */\n"
          "window.STRUCTURES=" + json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n")
    open(OUT, "w", encoding="utf-8").write(js)
    print("wrote", OUT, len(js.encode("utf-8")), "bytes")


if __name__ == "__main__":
    main()
