"""陈恪瑾 / U202413925 / 光电信息科学与工程202409班."""
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="仿真结果绘图、数值对照和报告整理")
    parser.add_argument("command", choices=("plot", "verify", "prepare"))
    parser.add_argument("run_directory", help="results/runs/ 下的运行目录")
    args = parser.parse_args()
    if args.command == "plot":
        from . import plot_results as command
    elif args.command == "verify":
        from . import verify_python as command
    else:
        from . import prepare_report as command
    sys.argv = [f"laser-tools {args.command}", args.run_directory]
    command.main()


if __name__ == "__main__":
    main()
