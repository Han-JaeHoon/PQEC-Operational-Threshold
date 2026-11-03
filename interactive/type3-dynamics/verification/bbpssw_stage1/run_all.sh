#!/bin/sh
set -eu
cd "$(dirname "$0")"
python verify.py
python -m pytest -q test_bbpssw.py
