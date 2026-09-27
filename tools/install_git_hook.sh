#!/usr/bin/env bash
# Installs tools/pre-commit as this repository's git pre-commit hook. Run once after cloning:
#   bash tools/install_git_hook.sh
set -e
cd "$(git rev-parse --show-toplevel)"
cp tools/pre-commit .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
echo "Pre-commit hook installed: commits containing data files will be refused."
