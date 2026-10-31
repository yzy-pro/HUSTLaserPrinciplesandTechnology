# 光学谐振腔自再现模式的仿真计算

作者：陈恪瑾

学号：U202413925

班级：光电信息科学与工程202409班

报告见 [report/report.pdf](report/report.pdf)，LaTeX 源文件见 [report/report.tex](report/report.tex)。程序完成对称共焦腔的一维、二维 Fox–Li 迭代与单镜倾斜扫描。

## 直接运行

`bin/` 中提供 **Linux x86_64** 可执行文件，已在 Ubuntu 24.04 / glibc 2.39 环境验证，无需另外安装 Python、uv 或 Rust。其他操作系统使用源码构建；本包不包含 Windows exe。

解压后保留目录结构，在项目根目录执行：

```bash
./bin/confocal-fox-li all
```

程序打印结果目录 `results/runs/<时间戳>/`。将下列 `RUN` 换成该目录名称：

```bash
./bin/laser-tools verify results/runs/RUN
./bin/laser-tools plot results/runs/RUN
./bin/laser-tools prepare results/runs/RUN
```

`verify` 完成独立数值对照；`plot` 生成图表；`prepare` 把插图和数值表更新到 `assert/` 和 `report/generated/`。前两条命令保留原始数据，不覆盖已有的检查或绘图清单。已有报告的数据目录见 `report/manifest.json`，无需再次计算即可查看报告和保存结果。

单独计算及帮助：

```bash
./bin/confocal-fox-li baseline --mode 1d
./bin/confocal-fox-li baseline --mode 2d
./bin/confocal-fox-li tilt
./bin/confocal-fox-li validate
./bin/confocal-fox-li --help
./bin/laser-tools --help
```

所有参数集中在 `config/parameter.yaml`，长度单位为 m，角度列表单位为 mrad。命令行只选择计算任务和维度。程序根据可执行文件位置查找参数，整个目录可以移动，不依赖原机器的源码路径。

## 源码构建

安装 Rust 2024 工具链、uv 和 Linux 常用编译工具后执行：

```bash
./code/build.sh
```

脚本使用 `cargo build --release --locked` 编译数值程序，使用 uv 锁定 Python 依赖，并通过 PyInstaller 生成包含 Python 解释器及依赖的单文件程序，输出到 `bin/`。

也可以通过 uv 直接运行 Python 工具：

```bash
uv sync --locked
uv run --locked laser-tools --help
uv run --locked laser-tools plot results/runs/RUN
```

源码检查：

```bash
cargo test --locked --manifest-path code/Cargo.toml
cargo clippy --locked --manifest-path code/Cargo.toml --all-targets -- -D warnings
```

## 编译报告

安装含 ctex 和 Fandol 字体的 XeLaTeX：

```bash
cd report
xelatex -interaction=nonstopmode -halt-on-error report.tex
xelatex -interaction=nonstopmode -halt-on-error report.tex
```

插图统一位于 `assert/`，表格和结果数值位于 `report/generated/`。

## 提交文件

| 路径 | 内容 |
|---|---|
| `report/` | PDF、LaTeX 源文件、数值表和数据来源清单 |
| `assert/` | 最终报告插图和来源清单 |
| `bin/` | 数值程序与 Python 工具的可执行文件 |
| `code/src/` | 数值程序源码 |
| `code/python/` | 绘图、数值对照与报告整理源码 |
| `config/parameter.yaml` | 唯一参数文件 |
| `results/` | 报告对应的场数据、迭代记录、扫描和数值检查 |
| `code/Cargo.lock`、`uv.lock` | 依赖版本锁定文件 |
| `submission_checks.json` | 本次构建、检查和可执行文件信息 |

大倾角达到迭代上限的结果在报告和数据中标为未收敛，这是已完成实验的结果记录。

重新生成提交压缩包：

```bash
uv run --locked python code/package_submission.py
```

压缩包只包含上表中的最终文件、运行说明和构建脚本，不包含 Git 历史、虚拟环境、编译缓存或 LaTeX 辅助文件。
