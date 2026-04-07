#!/bin/bash
set -o errexit

python -m ruff check mhctools tests

echo 'Passes ruff check'
