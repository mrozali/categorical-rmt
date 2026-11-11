#!/bin/sh
# verify.sh: one-command check of the reproducible results of this repository.
#
# usage: ./verify.sh [quick|full]
#
#   quick   (default) The self-test of the S3 code, then every script in code/figure-reproduction
#           that runs from the stored data: Figs. 4, 7 and 9, Table 1, Figs. 3, 5 and 6, and the
#           quick checks of Secs. 2.4.3 and 2.5.2. About half a minute in all.
#   full    The same, then regenerates both S3 ensembles from their seeds and compares them level
#           by level with the stored files, about 9 minutes each on the tested Mac.
#
# Environment: PYTHON names the interpreter, with the packages of environment/requirements.txt.
# The default is .venv/bin/python if it exists, otherwise python3.
#
# Each check is reported OK or FAIL. A check is OK when its script exits normally and prints
# RESULT: PASS (or, for the self-test, "self-test passed"). The full output of every check goes to
# code/figure-reproduction/generated/verify_logs/. The exit status is the number of failed checks.

set -u
cd "$(dirname "$0")" || exit 99

MODE=${1:-quick}
case "$MODE" in
  quick|full) ;;
  -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
  *) echo "usage: ./verify.sh [quick|full]" >&2; exit 2 ;;
esac

if [ -z "${PYTHON:-}" ]; then
  if [ -x .venv/bin/python ]; then PYTHON=.venv/bin/python; else PYTHON=python3; fi
fi

FR=code/figure-reproduction
LOG=$FR/generated/verify_logs
mkdir -p "$LOG"
fails=0

run() {
  name=$1
  shift
  start=$(date +%s)
  if "$@" > "$LOG/$name.txt" 2>&1 && grep -q -E "RESULT: PASS|self-test passed" "$LOG/$name.txt"; then
    echo "OK    $name ($(( $(date +%s) - start )) s)"
  else
    echo "FAIL  $name (output in $LOG/$name.txt)"
    tail -n 3 "$LOG/$name.txt" | sed 's/^/      /'
    fails=$((fails + 1))
  fi
}

echo "verify.sh $MODE, interpreter $PYTHON"
run s3_selftest "$PYTHON" -B code/s3/python/s3_chain.py
for s in fig07_fib_sectors_L26 check_fib_isospectrality check_fib_generators_smallL \
         fig09_u1_threshold_check s3_from_stored_ensembles fig04_s3_sff_goe_binned; do
  run "$s" "$PYTHON" "$FR/$s.py"
done
if [ "$MODE" = full ]; then
  run s3_regeneration_june "$PYTHON" "$FR/check_s3_regeneration.py" june
  run s3_regeneration_nov "$PYTHON" "$FR/check_s3_regeneration.py" nov
fi

echo "$fails failed"
exit "$fails"
