"""FBDF2（分数阶 BDF2）卷积权重。

对生成函数 δ(ζ) = 3/2 − 2ζ + ζ²/2（经典 BDF2），[δ(ζ)]^α 的幂级数
系数 {ω_k} 给出离散 Caputo 导数的卷积核：

    D_τ^α u^n = τ^(−α) · Σ_{k=0}^n ω_{n−k} (u^k − u⁰)

递推公式由恒等式 δ·W' = α·δ'·W（W(ζ) = [δ(ζ)]^α）逐项比较系数得到，
推导过程见 README.md。
"""

import numpy as np


def fractional_weights(alpha: float, Nt: int):
    """返回 numpy 数组 [ω_0, ω_1, …, ω_Nt]，ω_k 为 [δ(ζ)]^α 的系数。

    d0 = 3/2, d1 = −2, d2 = 1/2 是 δ(ζ) 的系数；递推式

        n·d0·ω_n = [(α+1−n)·d1]·ω_{n−1} + [(2(α+1)−n)·d2]·ω_{n−2}

    对 n ≥ 1 成立（n = 1 时第二项为 0）。ω_0 = (3/2)^α。
    """
    d0, d1, d2 = 3.0 / 2.0, -2.0, 1.0 / 2.0
    omega = np.zeros(Nt + 1)
    omega[0] = d0 ** alpha

    for n in range(1, Nt + 1):
        term = ((alpha + 1) - n) * d1 * omega[n - 1]
        if n >= 2:
            term += (2 * (alpha + 1) - n) * d2 * omega[n - 2]
        omega[n] = term / (n * d0)
    return omega
