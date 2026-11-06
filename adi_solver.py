"""FBDF2-ADI 求解器：零 Dirichlet 边界的时间分数阶 Allen-Cahn 方程。

    D_t^α u − μ·Δu + (u − u³) = 0,   (x, y) ∈ [0, 2π]²,  t ∈ (0, T]
    u(x, y, 0) = 0.05·(2·rand − 1)  （固定种子的均匀随机场，边界置 0）
    u|∂Ω = 0

算法要点（完整推导见 README.md）：
- 时间：分数阶 BDF2（FBDF2），历史卷积用初值修正恒等式实现；
- 空间：五点中心差分，ADI 分裂 + 混合导数修正（交叉项用一阶外推）；
- 简化：无 kappa 稳定化项；不做 u1 的起步计算（取 u¹ = u⁰）。
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from fractional_weights import fractional_weights


def adi_solver(alpha: float, T: float, tau: float, N: int,
               nodes, params: dict, seed: int):
    """计算一个 (T, alpha) 组合的数值解，返回 (U, X, Y)。

    U 形状 (N+1, N+1, Nt+1)，Nt = round(T/tau)；边界层恒为 0。
    """
    h = (nodes[1] - nodes[0]) / N
    x = np.linspace(nodes[0], nodes[1], N + 1)
    X, Y = np.meshgrid(x, x, indexing="ij")
    Nt = round(T / tau)

    # 1. FBDF2 权重 {ω_k} 与初值修正系数 b_n = Σ_{k≤n} ω_k
    omega = fractional_weights(alpha, Nt)
    b = np.cumsum(omega)
    tau_alpha = tau ** alpha
    a = params["mu"] * tau_alpha / omega[0]   # 隐式系数（无 kappa 项）

    # 2. (I − a·Dxx) 的三对角矩阵，按 LAPACK 带状格式存储；x、y 方向共用
    c = -a / h ** 2
    ab = np.vstack([np.full(N - 1, c),          # 上对角线
                    np.full(N - 1, 1 - 2 * c),  # 主对角线
                    np.full(N - 1, c)])         # 下对角线

    # 3. 初值：随机场（边界置 0）；u¹ = u⁰，不做起步计算
    rng = np.random.default_rng(seed)
    u0 = 0.05 * (2 * rng.random((N + 1, N + 1)) - 1)
    u0[0, :] = u0[-1, :] = u0[:, 0] = u0[:, -1] = 0.0

    U = np.zeros((N + 1, N + 1, Nt + 1))
    U[:, :, 0] = u0
    if Nt >= 1:
        U[:, :, 1] = u0

    Uin = U[1:N, 1:N, :]   # 内部点视图：边界层始终为 0，无需同步

    # 历史卷积用"时间轴在前的连续副本"做 BLAS 矩阵乘（直接法）
    Uflat = np.empty((Nt + 1, N - 1, N - 1))
    Uflat[0] = Uin[:, :, 0]
    if Nt >= 1:
        Uflat[1] = Uin[:, :, 1]

    # 4. 时间推进；第一个需要计算的时间层是 n = 2
    for n in range(2, Nt + 1):
        # (a) 历史卷积 Σ_{k=0}^{n-1} ω_{n-k} u^k
        weights = omega[1:n + 1][::-1]            # ω_n, …, ω_1
        history = (weights @ Uflat[:n].reshape(n, -1)).reshape(N - 1, N - 1)

        # (b) 初值修正 b_n u⁰
        init_corr = b[n] * Uin[:, :, 0]

        # (c) 非线性项二阶外推 2g(u^{n-1}) − g(u^{n-2})，g(u) = u − u³
        G = (2 * nonlinear_term(Uin[:, :, n - 1])
             - nonlinear_term(Uin[:, :, n - 2]))

        # (d) 右端项 + ADI 混合导数修正（交叉项一阶外推 a²·Dxy u^{n-1}）
        RHS = (-history + init_corr + tau_alpha * G) / omega[0]
        RHS = RHS + a ** 2 * cross_second_difference(U[:, :, n - 1], h)

        # (e) 两步带状三对角求解：(I − a·Dxx)U* = RHS，
        #     (I − a·Dyy)u^n = U*（零 Dirichlet 下中间层边界修正为 0）
        Ustar = solve_banded((1, 1), ab, RHS)
        Uin[:, :, n] = solve_banded((1, 1), ab, Ustar.T).T
        Uflat[n] = Uin[:, :, n]

    return U, X, Y


# --------------------------------------------------------------------------
# 辅助函数
# --------------------------------------------------------------------------

def cross_second_difference(U, h):
    """混合四阶差分 Dxy = Dx(Dy U)，用于 ADI 修正项。"""
    Dyy = (U[:, 2:] - 2 * U[:, 1:-1] + U[:, :-2]) / h ** 2
    return (Dyy[2:, :] - 2 * Dyy[1:-1, :] + Dyy[:-2, :]) / h ** 2


def nonlinear_term(U):
    """非线性项 g(u) = u − u³。"""
    return U - U ** 3
