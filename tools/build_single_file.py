# -*- coding: utf-8 -*-
"""把 index.html 与它引用的本地脚本内联成一个可直接双击打开的离线单文件。

用法： python tools/build_single_file.py
输出： 序列工坊_离线单文件版.html
"""
import io
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "index.html")
OUT = os.path.join(ROOT, "序列工坊_离线单文件版.html")


def inline_scripts(html):
    pattern = re.compile(r'<script src="([^"]+?)(?:\?[^"]*)?"\s*></script>')

    def repl(m):
        rel = m.group(1)
        path = os.path.join(ROOT, rel.replace("/", os.sep))
        with io.open(path, encoding="utf-8") as f:
            code = f.read()
        # 防止脚本内容里的 </script> 提前结束标签
        code = code.replace("</script", "<\\/script")
        return "<script>\n/* inlined: %s */\n%s\n</script>" % (rel, code)

    return pattern.sub(repl, html)


def main():
    html = io.open(SRC, encoding="utf-8").read()
    out = inline_scripts(html)
    banner = ("<!-- 序列工坊 · 离线单文件版（由 tools/build_single_file.py 生成）\n"
              "     图表、三维引擎与 PDB 结构数据均已内联，双击即可运行；PWA 安装需 HTTPS/localhost。 -->\n")
    out = out.replace("<!DOCTYPE html>", "<!DOCTYPE html>\n" + banner, 1)
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(out)
    size = os.path.getsize(OUT)
    print("wrote", OUT, "%.2f MB" % (size / 1024 / 1024))
    print("remaining external scripts:", len(re.findall(r'<script src=', out)))


if __name__ == "__main__":
    main()
