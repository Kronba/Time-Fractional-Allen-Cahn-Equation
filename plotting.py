"""演化图绘制：T_list（行）× alpha_list（列）的子图网格。

每个面板从对应的 .mat 数据文件读取终值解，用 pcolormesh 绘制
（等价于 MATLAB 的 surf(..., 'EdgeColor', 'none') + view(2)）。
"""

import matplotlib

matplotlib.use("Agg")  # 无显示环境也能画图
import matplotlib.pyplot as plt  # noqa: E402
from scipy.io import loadmat  # noqa: E402

import config  # noqa: E402


def plot_grid(T_list, alpha_list, out_file) -> None:
    """按 T（行，上→下）× alpha（列，左→右）绘制演化图。"""
    letters = "abcdefghijklmnopqrstuvwxyz"
    rows, cols = len(T_list), len(alpha_list)

    # 每个子图约 250×267 px（沿袭 MATLAB 1000×800 版式）
    fig, axes = plt.subplots(rows, cols,
                             figsize=(250 * cols / 150, 267 * rows / 150),
                             squeeze=False)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.06,
                        wspace=0.12, hspace=0.34)  # 近似 tiledlayout compact

    for i, T in enumerate(T_list):
        for j, alpha in enumerate(alpha_list):
            ax = axes[i][j]
            k = i * cols + j
            path = config.data_file(T, alpha)
            if not path.exists():
                ax.axis("off")
                ax.set_title(f"({letters[k]}) missing", fontsize=8)
                continue

            dat = loadmat(str(path))
            ax.pcolormesh(dat["X"], dat["Y"], dat["U_final"],
                          shading="nearest")
            ax.set_title(f"({letters[k]}) T = {T:g}, "
                         f"$\\alpha$ = {alpha:g}", fontsize=8)
            ax.set_aspect("equal")
            ax.tick_params(labelsize=6)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=150)
    plt.close(fig)
