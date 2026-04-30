#!/usr/bin/env bash
set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/TonyZhang2004/Character_Table_of_Symmetric_Groups.git}"
REPO_DIR="${REPO_DIR:-$HOME/Character_Table_of_Symmetric_Groups}"
BRANCH="${BRANCH:-main}"
K="${K:-13}"
N="${N:-}"
MAX_MEMO_ENTRIES="${MAX_MEMO_ENTRIES:-5000000}"
LOG_INTERVAL="${LOG_INTERVAL:-10000}"
CHECKPOINT_INTERVAL="${CHECKPOINT_INTERVAL:-0}"
MEMO_FILE="${MEMO_FILE:-}"
RUN_TESTS="${RUN_TESTS:-1}"
DETACH="${DETACH:-1}"
FRESH="${FRESH:-0}"

if ! command -v git >/dev/null 2>&1; then
  echo "git is required. Install it first, then rerun this script." >&2
  exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required. Install it first, then rerun this script." >&2
  exit 1
fi

if [ -z "$N" ]; then
  N=$((K * (K + 1) / 2))
fi

if [ "$N" -le 78 ]; then
  echo "This VM helper is intended for triangular n larger than 78; got N=$N." >&2
  echo "Set K=13 for N=91, K=14 for N=105, etc." >&2
  exit 1
fi

OUTPUT_DIR="${OUTPUT_DIR:-$REPO_DIR/runs/staircase/S$N}"

if [ -d "$REPO_DIR/.git" ]; then
  echo "Updating existing repo at $REPO_DIR"
  git -C "$REPO_DIR" fetch origin "$BRANCH"
  git -C "$REPO_DIR" switch "$BRANCH"
  git -C "$REPO_DIR" pull --ff-only origin "$BRANCH"
else
  echo "Cloning $REPO_URL into $REPO_DIR"
  git clone --branch "$BRANCH" "$REPO_URL" "$REPO_DIR"
fi

cd "$REPO_DIR"

python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install numpy

if [ "$RUN_TESTS" = "1" ]; then
  PYTHONPYCACHEPREFIX=/tmp/mn_pycache .venv/bin/python -m unittest -v test_murnaghan_nakayama.py
fi

mkdir -p "$OUTPUT_DIR"

args=(
  run_staircase_column.py
  --n "$N"
  --output-dir "$OUTPUT_DIR"
  --max-memo-entries "$MAX_MEMO_ENTRIES"
  --log-interval "$LOG_INTERVAL"
  --checkpoint-interval "$CHECKPOINT_INTERVAL"
)

if [ -n "$MEMO_FILE" ]; then
  args+=(--memo-file "$MEMO_FILE")
fi

if [ "$FRESH" = "1" ]; then
  args+=(--fresh)
fi

echo "Triangular parameter K: $K"
echo "Symmetric group size N: $N"
echo "Output directory: $OUTPUT_DIR"
echo "CSV file: $OUTPUT_DIR/S${N}_staircase.csv"
echo "Progress file: $OUTPUT_DIR/S${N}_staircase.csv.progress.json"
echo "Log file: $OUTPUT_DIR/S${N}_staircase.log"

if [ "$DETACH" = "1" ]; then
  nohup .venv/bin/python "${args[@]}" > "$OUTPUT_DIR/S${N}_staircase.nohup.out" 2>&1 &
  pid="$!"
  echo "$pid" > "$OUTPUT_DIR/S${N}_staircase.pid"
  echo "Started detached staircase job with PID $pid"
  echo "Monitor with: tail -f $OUTPUT_DIR/S${N}_staircase.log"
else
  .venv/bin/python "${args[@]}"
fi
