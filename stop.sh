#!/usr/bin/env bash
# Stop a running campaign without the dashboard.
#   ./stop.sh            # stop now: terminate every job (they are retried on the next start)
#   ./stop.sh --drain    # finish running jobs, start nothing new
#   ./stop.sh --run-id X
#
# This writes one word into the run's `control` file; the controller reads it within seconds.
# "stop" kills every job's whole process tree (including QE's mpirun, which runs in its own
# session) and marks their ledger runs failed, so they are retried on the next start; an
# interrupted calculation is never used as a result. docs/RUNNER.md §3.
set -euo pipefail
RTT_HOME="${RTT_HOME:-$HOME/.rtt}"
# shellcheck disable=SC1091
source "$RTT_HOME/env.sh"
action=stop
args=()
for a in "$@"; do
  case "$a" in --drain) action=drain ;; *) args+=("$a") ;; esac
done
exec rtt campaign "$action" "${args[@]}"
