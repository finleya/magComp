#!/usr/bin/env bash
# Publish the reviewed site after running Rscript dev/build_articles.R.
set -euo pipefail
cd "$(dirname "$0")/.."

# Require a completed website build before replacing the published files.
test -f docs/index.html
test -f docs/.nojekyll

site_dir=$(mktemp -d)
trap 'rm -rf "$site_dir"' EXIT

git clone --quiet --single-branch --branch gh-pages \
  "$(git remote get-url origin)" "$site_dir"
rsync -a --delete --exclude=.git docs/ "$site_dir/"
git -C "$site_dir" add --all

if git -C "$site_dir" diff --cached --quiet; then
  echo "The published files already match docs/."
  exit 0
fi

git -C "$site_dir" \
  -c user.name="$(git config user.name)" \
  -c user.email="$(git config user.email)" \
  commit -m "Update website for $(git rev-parse --short HEAD)"
git -C "$site_dir" push origin gh-pages
