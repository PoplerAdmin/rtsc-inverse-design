"""Portable campaign runner: hardware detection, the background controller and the TUI.

* :mod:`rtt.campaign.hw` reads cores / RAM / GPU / disk and writes a ``machine.toml`` sized to
  this machine (``rtt init``).
* :mod:`rtt.campaign.home` resolves the install home (``RTT_HOME``, default ``~/.rtt``) and the
  run directory layout.
* :mod:`rtt.campaign.controller` is the long-running process behind ``rtt campaign start``: it
  runs the front end (S1-S3) and then keeps the CPU busy with S4-S6 jobs, each a child process
  in its own process group so that a stop request can end all of them.
* :mod:`rtt.campaign.ui` is the terminal dashboard (``rtt ui``, ``start.sh``).
"""
