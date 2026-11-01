# v4 验证记录 · 2026-11-01

## 独立计算对照

参考版本 Biopython 1.86，仅装在开发临时目录，不是运行依赖。生成脚本 tools/reference_fixtures.py；静态参考数据 tests/fixtures/biopython-1.86.json。

- 表1和表2共128个密码子翻译一致。
- 14组 Tm（固定示例、随机种子20261101、自互补、15–30bp）：与 MeltingTemp.Tm_NN 的 DNA_NN3、saltcorr=5、Na=50mM 对照，误差<0.01°C。
- 总寡核苷酸浓度250nM：普通配对 dnac1=125/dnac2=125，自互补 dnac1=250/dnac2=0，并应用相应对称修正。

对照发现并修正旧实现的 GC/CG 焓熵参数、混合初始化参数及无条件自互补修正。这不等于全部引物设计算法已经与 Primer3 一致。

## 自动化回归

npm test 覆盖语法、核心算法冒烟、v3生信、11组产品回归、Service Worker 离线回退与缓存隔离、静态资源/单文件、应用壳内容指纹。

| 输入/操作 | 核验结果 |
|---|---|
| ATGNNTAA | TTANNCAT；GC=12.5%、AT=62.5%、N=25%；未知Tm/分子量 |
| 多条FASTA | 报错，不拼接 |
| 390bp示例，第87位T→A，+1 | TAT→TAA，无义突变，同一ORF129aa→28aa |
| 相同突变的反链 -1 表示 | 同样判为无义突变 |
| ATGTTATAA，ecoli | ATGCTGTAA，改变1/3，排名分数0.707→1.000 |
| 重新分析 | 旧引物、突变、比对、优化结果与地图衍生标记清空 |
| 越界或未知格式备份 | 拒绝，保留原工作区 |
| 注释含HTML字符 | 按文字转义显示 |
| 48,000bp无完整ORF纯计算 | 返回无ORF；不代表完整UI压测通过 |

## 浏览器实际运行

经用户授权，独立 Edge 156 / Windows / localhost。桌面1440×1000，移动模拟320、390、768px。

- 首访引导自动聚焦，Tab约束在弹窗。
- 390bp示例分析、第87位无义突变和28aa产物。
- 引物20bp、目标55°C、产物300bp：5组候选及重复匹配风险提示。
- 真实输入验证N、多FASTA及旧结果清空。
- 1条注释、环状拓扑工作区完成备份与计算记录下载；JSON有效；错误schema被拒绝。
- 本轮检查的320/390/768px页面无整页横向溢出，抽屉可打开。
- 浏览器网络设为offline后刷新，页面及核心模块仍加载；9bp与1条注释工作区可以恢复并继续分析。此项为实际断网刷新验证。

## 未验证

实体手机、Safari/Firefox、真实PWA安装/卸载升级、多设备用户研究、完整读屏与对比度验收、全部复杂GenBank、完整50kb交互压力、基因组范围引物特异性、实验表达效果及科研级验证。

## 来源

- [NCBI遗传密码表](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- [Biopython 1.86 MeltingTemp](https://github.com/biopython/biopython/blob/biopython-186/Bio/SeqUtils/MeltingTemp.py)
- [SantaLucia 1998](https://doi.org/10.1073/pnas.95.4.1460)

修改后应重跑受影响的回归并更新缓存内容指纹，不将有限检查描述为全功能或科研准确性认证。
