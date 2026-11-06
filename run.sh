#!/usr/bin/env bash
# 一键配置环境并运行（./run.sh；./run.sh setup 只装环境不运行）。
# 虚拟环境复用顺序：本目录 .venv → 上级目录 .venv → 新建。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ -x .venv/bin/python ]; then
    VENV=.venv
elif [ -x ../.venv/bin/python ]; then
    VENV=../.venv
else
    VENV=.venv
    echo "[run.sh] 创建虚拟环境 $VENV ..."
    python3 -m venv "$VENV" 2>/dev/null || python3 -m venv --without-pip "$VENV"
    if [ ! -x "$VENV/bin/pip" ]; then
        # 系统没有 ensurepip 时用 get-pip.py 引导
        echo "[run.sh] 引导 pip ..."
        curl -sSL https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py \
            || wget -qO /tmp/get-pip.py https://bootstrap.pypa.io/get-pip.py
        "$VENV/bin/python" /tmp/get-pip.py
        rm -f /tmp/get-pip.py
    fi
fi

echo "[run.sh] 安装依赖 ..."
"$VENV/bin/python" -m pip install -r requirements.txt

if [ "${1:-}" = "setup" ]; then
    echo "[run.sh] 环境就绪（仅配置，未运行）。"
    exit 0
fi

echo "[run.sh] 运行 main.py ..."
"$VENV/bin/python" main.py
