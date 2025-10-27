#!/bin/bash
# Release script for sboxanalyzer
# Usage: ./scripts/release.sh <version>
# Example: ./scripts/release.sh 1.0.2

set -e  # Exit on error

if [ -z "$1" ]; then
    echo "Usage: ./scripts/release.sh <version>"
    echo "Example: ./scripts/release.sh 1.0.2"
    exit 1
fi

VERSION=$1
TAG="v${VERSION}"

echo "Preparing release ${TAG}"

# Check if working directory is clean
if [ -n "$(git status --porcelain)" ]; then
    echo "ERROR: Working directory is not clean. Please commit or stash changes first."
    exit 1
fi

# Update version in pyproject.toml
echo "Updating version in pyproject.toml..."
sed -i.bak "s/^version = \".*\"/version = \"${VERSION}\"/" pyproject.toml
rm pyproject.toml.bak

# Update version in __init__.py
echo "Updating version in __init__.py..."
sed -i.bak "s/^__version__ = \".*\"/__version__ = \"${VERSION}\"/" sboxanalyzer/__init__.py
rm sboxanalyzer/__init__.py.bak

# Commit version changes
echo "Committing version changes..."
git add pyproject.toml sboxanalyzer/__init__.py
git commit -m "Bump version to ${VERSION}"

# Create and push tag
echo "Creating tag ${TAG}..."
git tag -a "${TAG}" -m "Release ${TAG}"

echo "SUCCESS: Version ${VERSION} prepared!"
echo ""
echo "Next steps:"
echo "1. Review the changes: git show"
echo "2. Push to GitHub: git push && git push --tags"
echo "3. GitHub Actions will automatically publish to PyPI"
echo ""
echo "Or to cancel: git reset --hard HEAD~1 && git tag -d ${TAG}"
