#!/usr/bin/env bash
# Open the rtsc-targeted dashboard. From there: s = start, d = drain, x = stop now, q = quit.
#
#   ./start.sh                     # run "default" ($RTT_HOME/runs/default)
#   ./start.sh --run-id trial2     # another run
#   ./start.sh --rounds 0          # the start key runs rounds until stopped
#
# The dashboard is only a window onto the campaign: the work runs in a detached controller
# process, so closing this terminal (or losing an SSH session) does not stop it. Reopen with
# ./start.sh at any time. Why and how: docs/RUNNER.md §2-3.
set -euo pipefail
RTT_HOME="${RTT_HOME:-$HOME/.rtt}"
if [ ! -f "$RTT_HOME/env.sh" ]; then
  echo "Not installed yet. Run ./install.sh first." >&2
  exit 1
fi
# shellcheck disable=SC1091
source "$RTT_HOME/env.sh"
if ! rtt doctor > "$RTT_HOME/doctor.txt" 2>&1; then
  cat "$RTT_HOME/doctor.txt"
  read -r -p "Some required checks failed. Open the dashboard anyway? [y/N] " ans
  case "$ans" in y|Y|yes) ;; *) exit 1 ;; esac
fi
exec rtt ui "$@"
