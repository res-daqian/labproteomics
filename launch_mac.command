#!/bin/bash
# Double-click to start the Lab Proteomics Workbench (macOS). First time: right-click > Open.
cd "$(dirname "$0")"
if ! command -v Rscript >/dev/null 2>&1; then
  export PATH="/usr/local/bin:/opt/homebrew/bin:/Library/Frameworks/R.framework/Resources/bin:$PATH"
fi
Rscript scripts/launch_app.R
