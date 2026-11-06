"""数值算例的全部可调参数（论文 Example 2：时间分数阶 Allen-Cahn，f = 0）。

修改本文件的常量后直接重跑 main.py：所有 (T, alpha) 组合都会全部重新
计算，图像也会重新生成。
"""

import math
from pathlib import Path

# ---- 求解参数 -----------------------------------------------------------
ALPHA_LIST = [1.0, 0.7, 0.5, 0.3]   # 分数阶阶数，对应图的列（左→右）
T_LIST = [0.0, 2.0, 8.0, 16.0]      # 终止时间，对应图的行（上→下）
N = 256                             # 每个方向的区间数（Nx = Ny）
TAU = 0.005                         # 时间步长
NODES = [0, 2 * math.pi]            # 空间区域 [0, 2π]^2
SEED = 10086                        # 随机初值 u0 = 0.05*(2*rand-1) 的种子
PARAMS = {"example": 2, "mu": 0.1 ** 2, "eps": 0.1, "p": 1}
# 求解器只读 PARAMS["mu"]；example/eps/p 仅随 .mat 保存作记录
# （p 是原 MATLAB 算例 1 精确解 sin(p·(x+y)) 里的频率，本算例不使用）。

# ---- 输出目录（固定在本目录下，与运行目录无关） ---------------------------
OUT_ROOT = Path(__file__).resolve().parent / "case5_results"


def data_file(T: float, alpha: float) -> Path:
    """(T, alpha) 组合对应的 .mat 数据文件路径。"""
    return OUT_ROOT / "data" / f"T{T:g}_alpha{alpha:g}.mat"


def figure_file() -> Path:
    """演化图输出路径（文件名按 行×列 自动命名）。"""
    return OUT_ROOT / "figures" / f"case5_{len(T_LIST)}x{len(ALPHA_LIST)}.png"
