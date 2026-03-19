#!/usr/bin/env bash
set -euo pipefail

# Install pre-commit hooks for AnnoQC development
# This script sets up git hooks to run formatting, linting, and tests before commits

echo "Installing pre-commit hooks for AnnoQC..."

# Check if pre-commit is installed
if ! command -v pre-commit >/dev/null 2>&1; then
    echo "pre-commit is not installed."
    echo ""
    echo "To install pre-commit, choose one of the following options:"
    echo ""
    echo "  Via pip:"
    echo "    pip install pre-commit"
    echo ""
    echo "  Via pixi (if using pixi environment):"
    echo "    pixi install pre-commit"
    echo ""
    exit 1
fi

# Install the hooks
pre-commit install

echo ""
echo "✓ Pre-commit hooks installed successfully!"
echo ""
echo "The following hooks will now run before each commit:"
echo "  - cargo fmt (code formatting)"
echo "  - cargo clippy (linting with -D warnings)"
echo "  - cargo test (run all tests)"
echo "  - trailing whitespace, end-of-file-fixer (general file checks)"
echo "  - check-yaml, check-toml (config file validation)"
echo ""
echo "To run all hooks manually on all files:"
echo "  pre-commit run --all-files"
echo ""
