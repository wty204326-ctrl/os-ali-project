# _activate_rna.sh — 定位并激活 conda 环境 rna
# 由 00_install_deps.sh / run_all.sh source
activate_rna() {
  if [ -n "${CONDA_PREFIX:-}" ] && [ "$(basename "$CONDA_PREFIX")" = "rna" ]; then
    return 0
  fi
  local act
  for act in \
    "${CONDA_EXE:+${CONDA_EXE%/*}/activate}" \
    "$HOME/miniforge3/bin/activate" \
    "$HOME/mambaforge/bin/activate" \
    "$HOME/miniconda3/bin/activate" \
    "$HOME/anaconda3/bin/activate"; do
    [ -n "$act" ] && [ -f "$act" ] || continue
    # conda activate 脚本依赖当前 shell
    # shellcheck disable=SC1090
    source "$act" rna && return 0
  done
  return 1
}
