#!/bin/bash
# Double-click to start the workbench with the demo button pointing at the thesis data in golden_data/
# (Daqian's own Long-term retina runs; local only, never shipped). First time: right-click > Open.
cd "$(dirname "$0")"
if ! command -v Rscript >/dev/null 2>&1; then
  export PATH="/usr/local/bin:/opt/homebrew/bin:/Library/Frameworks/R.framework/Resources/bin:$PATH"
fi
Rscript scripts/launch_app.R golden_data
