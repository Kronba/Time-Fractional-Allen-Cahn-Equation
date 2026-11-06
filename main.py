"""入口：求解全部 (T, alpha) 组合 → 保存 .mat → 绘制演化图。

参数都在 config.py；运行方式见 README.md（或直接 ./run.sh）。
"""

import config
from scipy.io import savemat

from adi_solver import adi_solver
from plotting import plot_grid


def run_case(T: float, alpha: float) -> None:
    """计算一个 (T, alpha) 组合并保存终值解（总是重算、覆盖旧文件）。"""
    U, X, Y = adi_solver(alpha, T, config.TAU, config.N, config.NODES,
                         config.PARAMS, config.SEED)
    savemat(str(config.data_file(T, alpha)),
            {"U_final": U[:, :, -1], "X": X, "Y": Y,
             "T": T, "alpha": alpha, "seed": config.SEED,
             "params": config.PARAMS}, do_compression=True)
    print(f"saved T = {T:g}, alpha = {alpha:g}")


def main() -> None:
    (config.OUT_ROOT / "data").mkdir(parents=True, exist_ok=True)
    for T in config.T_LIST:
        for alpha in config.ALPHA_LIST:
            run_case(T, alpha)

    fig_file = config.figure_file()
    plot_grid(config.T_LIST, config.ALPHA_LIST, fig_file)
    print(f"finished: 数据在 {config.OUT_ROOT / 'data'}，图在 {fig_file}")


if __name__ == "__main__":
    main()
