#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

cd "$PROJECT_DIR"

python3 - "$TMP_DIR" <<'PY'
from pathlib import Path
import subprocess
import sys

root = Path.cwd()
outdir = Path(sys.argv[1])
tracked = subprocess.check_output(
    ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"]
).split(b"\0")
shell = []
python = []
for raw_path in tracked:
    if not raw_path:
        continue
    relative = raw_path.decode("utf-8", errors="surrogateescape")
    path = root / relative
    if not path.is_file():
        continue
    first_line = path.open("rb").readline(256)
    if relative.endswith(".sh") or (
        first_line.startswith(b"#!") and b"sh" in first_line
    ):
        shell.append(relative)
    if relative.endswith(".py") or (
        first_line.startswith(b"#!") and b"python" in first_line
    ):
        python.append(relative)

for name, paths in (("shell", shell), ("python", python)):
    if not paths:
        raise SystemExit(f"No {name} scripts found.")
    (outdir / f"{name}.txt").write_text("\n".join(paths) + "\n")
PY

mapfile -t shell_scripts < "$TMP_DIR/shell.txt"
mapfile -t python_scripts < "$TMP_DIR/python.txt"

printf 'Checking shell syntax (%d files)...\n' "${#shell_scripts[@]}"
bash -n "${shell_scripts[@]}"

printf 'Checking Python syntax (%d files)...\n' "${#python_scripts[@]}"
PYTHONPYCACHEPREFIX="$TMP_DIR/pycache" python3 -m py_compile \
    "${python_scripts[@]}"

printf 'Lint passed for %d shell and %d Python files.\n' \
    "${#shell_scripts[@]}" "${#python_scripts[@]}"
