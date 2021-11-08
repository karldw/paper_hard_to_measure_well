#!/bin/bash

set -euf -o pipefail

# Get the directory the setup_hooks.sh script is in:
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEST=$DIR/../../.git/hooks

# See directions in http://mirrors.ctan.org/macros/latex/contrib/gitinfo2/gitinfo2.pdf

cp -pi "$DIR/post-checkout" "$DEST/"
cp -pi "$DIR/post-commit"   "$DEST/"
cp -pi "$DIR/post-merge"    "$DEST/"
chmod +x "$DEST/post-checkout"
chmod +x "$DEST/post-commit"
chmod +x "$DEST/post-merge"

