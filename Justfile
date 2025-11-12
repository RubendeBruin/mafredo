# Justfile — project tasks (cookiecutter-aware)
#
# Default shell on Unix; PowerShell on Windows (rm is aliased to Remove-Item)
set shell := ["bash", "-cu"]
set windows-shell := ["powershell.exe", "-NoLogo", "-Command"]

# Run `just` with no arguments to show menu
default: help

# ---------------------------------------
# Setup
# ---------------------------------------

# Install the virtual environment and install the pre-commit hooks
install:
    @echo "🚀 Creating virtual environment using uv"
    @uv sync
    @uv run pre-commit install

install-dev:
    @echo "🚀 Creating virtual development environment using uv"
    @uv sync --dev
    @uv run pre-commit install

# ---------------------------------------
# Quality
# ---------------------------------------

# Run code quality tools
check:
    @echo "🚀 Checking lock file consistency with 'pyproject.toml'"
    @uv lock --locked
    @echo "🚀 Linting code: Running pre-commit"
    @uv run pre-commit run -a
    @echo "🚀 Static type checking: Running mypy"
    @uv run mypy
    @echo "🚀 Checking for obsolete dependencies: Running deptry"
    @uv run deptry .

# ---------------------------------------
# Tests
# ---------------------------------------

# Test the code with pytest
test:
    @echo "🚀 Testing code: Running pytest"
    @uv run python -m pytest

# Test the code with pytest and generate a coverage report
test-cov:
    @echo "🚀 Testing code with coverage: Running pytest"
    @uv run python -m pytest --cov --cov-config=pyproject.toml --cov-report=term-missing

# ---------------------------------------
# Build & Clean
# ---------------------------------------

# Clean build and examples
clean: clean-build

# Build wheel file
build: clean-build
    @echo "🚀 Creating wheel file"
    @uvx --from build pyproject-build --installer uv

# Clean build artifacts
clean-build:
    @echo "🧹 Removing build artifacts"
    @uv run python -c "import shutil; shutil.rmtree('dist', ignore_errors=True)"

# ---------------------------------------
# Publish
# ---------------------------------------

# Publish a release to PyPI
publish:
    @echo "🚀 Publishing to PyPI."
    @uvx twine upload --repository pypi dist/*

# Publish a release to TestPyPI
publish-test:
    @echo "🚀 Publishing to TestPyPI."
    @uvx twine upload --repository testpypi dist/*

# Build and publish to PyPI
build-and-publish: build publish

# Build and publish to TestPyPI
build-and-publish-test: build publish-test


# ---------------------------------------
# Docs
# ---------------------------------------

# Test if documentation can be built without warnings or errors
docs-test:
    @uv run mkdocs build -s

# Build and serve the documentation
docs:
    @uv run mkdocs serve


# ---------------------------------------
# Help / menu
# ---------------------------------------

# List all tasks and their descriptions
help:
    @just --list
