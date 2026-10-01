# 交接记录

更新日期：2026 年 10 月 1 日。

接下来继续和用户逐页看 `ICAMSbeamer/talk.pdf`，讲清楚图和代码。今天主要看了正文第 1 页，也就是封面后面放着三张图的那一页。还没有开始逐页核对后面的内容。

用户明确要求：以后都用白话，简单直接。解释代码时，先说它拿了什么、做了什么，再说细节。要分清哪些是 Kanapy 自带的，哪些是制作报告时另外加的。不能把没有依据的选择说成教授的要求。

## 文件在哪里

| 文件 | 用途 |
| --- | --- |
| [ICAMSbeamer/talk.pdf](ICAMSbeamer/talk.pdf) | 当前报告，共 21 页：封面、14 页正文、6 页补充材料 |
| [ICAMSbeamer/talk.tex](ICAMSbeamer/talk.tex) | 报告正文和排版的可编辑源文件 |
| [ICAMSbeamer/speaker_notes_zh.md](ICAMSbeamer/speaker_notes_zh.md) | 中文逐页讲解 |
| [ICAMSbeamer/generate_figures.py](ICAMSbeamer/generate_figures.py) | 生成配图的脚本 |
| [ICAMSbeamer/README.md](ICAMSbeamer/README.md) | 报告说明和复现方法 |
| [ICAMSbeamer/figure_provenance.json](ICAMSbeamer/figure_provenance.json) | 图片来源、参数和计算记录 |
| [ICAMSbeamer/context_descriptor.json](ICAMSbeamer/context_descriptor.json) | 第 1 页例子的输入设置 |
| [ICAMSbeamer/context_generation.json](ICAMSbeamer/context_generation.json) | 该例子的源码版本和随机种子 |
| [ICAMSbeamer/context_apd.npz](ICAMSbeamer/context_apd.npz) | 保存的 APD 参数和体素归属 |
| [ICAMSbeamer/context_rve_low_fill.pkl](ICAMSbeamer/context_rve_low_fill.pkl) | 已生成 RVE 的缓存，重画时避免重复计算 |

`ICAMSbeamer` 还包含所有配图、主题文件和 `logos`。这些文件一起保留，笔记本上就能直接看报告和重新编译。

## 今天已经核对清楚的事

第 1 页的三张图来自 Kanapy 6.5.5 实际生成的同一个例子：23 个晶粒，盒子边长 24 微米，非周期边界，随机种子 `20261001`。原始几何生成时的源码提交是 `19713e5612fdaaea2acc4378b2c897be3f54dd2d`。体素为 `20 × 20 × 20`，表面共有 7554 个三角形。记录显示，摆放后的重叠接触数为零。

相关代码在 `generate_figures.py` 的 `contextual_rve()`，目前从第 401 行开始。

| 图 | 数据从哪里来 | 怎样画出来 |
| --- | --- | --- |
| `pipeline_ellipsoids.png`，左图 | `ms.pack()` 之后的椭球 | 第 457 行调用 Kanapy 自带的 `plot_ellipsoids_3D()` |
| `pipeline_voxels.png`，中图 | `ms.voxelize()` 得到的体素归属 | 第 460 行调用 Kanapy 自带的 `plot_voxels_3D()` |
| `pipeline_surface.png`，右图 | `ms.generate_grains()` 得到的共享表面 | 从第 465 行开始读取表面，由报告脚本用 Matplotlib 画出来 |

左、中图的绘图函数在 `src/kanapy/core/plotting.py`。颜色、视角、画面大小等展示设置由报告脚本调整。

右图的三角形也是 Kanapy 生成的。图上黑线的位置沿着这些三角形的边；把它们显示成深色细线，是报告脚本的选择。脚本没有另外划分三角形。中、右图隐藏了一角，方便看内部。

`ms` 是脚本起的变量名。`geometry` 是 Kanapy 对象自带的属性，`"Surface"` 是其中保存共享表面的数据项。

```python
surface = ms.geometry["Surface"]
```

这行只是把已算好的表面数据拿出来。里面有每个顶点的位置、每个三角形由哪三个点组成，以及三角形两边的晶粒编号。这行本身不重新算网格，也不画图。

三角形的顶点怎样得到：

1. 脚本设了 `ms.generate_grains(resolution=6)`。Kanapy 把 24 微米的盒子每个方向分成 6 格，格点间隔为 4 微米。每个小方格再拆成 6 个四面体，也就是四个角的小立体。
2. Kanapy 在这些小立体里算晶粒分界，得到切分后的面和角点。这一步用的是格点处计算结果在小立体里的插值，是对连续 APD 的近似。
3. 面已有三个角就直接保留。超过三个角时，取各个角的坐标平均值，补一个中间点，再连接周围相邻的角，拼成三角形。

所以顶点既可能是原格点，也可能是新算出的交点，或者后补的中间点。格子的粗细由脚本设置，具体顶点由 Kanapy 自动计算。不是所有点都等距，也不是每个三角形都另加了中间点。

对应源码：`src/kanapy/core/apd_geometry.py` 的 `build_grain_geometry()`，`src/kanapy/core/apd_mesh.py` 的背景格子和切分计算，以及 `src/kanapy/core/apd_boundary.py` 的 `triangulate()`。

Kanapy 还用这些三角片计算表面积和晶粒体积。三角形表示表面在三维绘图和网格软件里很常见，但也有四边形等其他做法。

## 接着讨论时要记住

第 1 页中图和右图，是对同一个 APD 的两种表示。右图不是把中图的阶梯外皮磨光得到的；表面三角形也不是最终的有限元体网格。

目前没有教授明确要求使用三角网格或显示三角形边线的原文依据。不能把这些选择归到教授的要求上。

第 1 页是为报告生成的说明例子。后面的点接触例子来自仓库里的可控测试，局部修复例子是人为构造的候选。这些都不是教授实际出错的那个 RVE。自动搜索修复方案、连续分区修复和最终力学验证仍是拟开展的工作。

本次逐页讨论只核对了第 1 页的来源和相关代码，没有改图、改算法或重新生成数据。下一次从用户当前的问题继续，讲清楚后再往后翻页。报告原文和中文讲稿可以作为定位材料，但不能当成教授原始要求的证据。

## 笔记本上怎么继续

仓库是 `git@github.com:JXRepo/Kanapy.git`，分支是 `jun`。

第一次下载：

```bash
git clone --branch jun git@github.com:JXRepo/Kanapy.git kanapy
cd kanapy
```

如果笔记本已有这个仓库，在仓库根目录更新：

```bash
git switch jun
git pull --ff-only origin jun
```

先读这份交接说明，再直接打开 `ICAMSbeamer/talk.pdf`。现有 PDF 和图片都已保存，继续看报告不用重新跑 Kanapy。

修改 `talk.tex` 后重新生成 PDF，需要安装 LaTeX、`latexmk` 和 `make`。然后运行：

```bash
cd ICAMSbeamer
make
```

只有需要改图时，才准备 Python 环境。仓库支持 Python 3.10 到 3.13。从仓库根目录运行，已有 `knpy` 环境时可以跳过创建：

```bash
conda env create -f environment.yml
conda activate knpy
python -m pip install .
python ICAMSbeamer/generate_figures.py --only context
```

`--only context` 只处理第 1 页这组三图；缓存可用时直接重画，缓存不存在时会重新生成这个 RVE。还有 `controls`、`local`、`periodic`；不加 `--only` 会处理全部图组。

`ICAMSbeamer/README.md` 中 `/home/users/xuejungs/.../python` 是办公室这台机器的路径。笔记本上激活自己的环境后用 `python`，不要照抄那个绝对路径。脚本会从当前仓库的 `src` 导入 Kanapy。

脚本会优先读取 `.pkl` 缓存。若换环境后缓存读不出来，先把缓存移到目录外备份，再运行 `--only context` 重新生成。更换依赖版本可能改变重新计算的结果；保留现有 PNG 和 PDF，就能保留当前这版报告。

## Git 状态和推送安排

已从本机 `.git/info/exclude` 删除 `/ICAMSbeamer/`。这份本地忽略配置不随 Git 提交，但 `ICAMSbeamer` 的文件现在可以正常添加。

用户已明确授权：把整个 `ICAMSbeamer` 和根目录这份 `HANDOFF.md` 一起提交，并推送到自己的 `origin/jun`。推送目标是 `JXRepo/Kanapy`，不是 `ICAMS/Kanapy`。

开始交接提交前，本地 `jun` 和刚获取的 `origin/jun` 都在 `19713e56`，已跟踪文件没有未提交修改。最终提交和远端状态以 `git log -1`、`git status` 及远端为准。
