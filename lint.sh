#!/bin/bash
set -o errexit

python -m ruff check mhctools tests scripts

echo 'Passes ruff check'
