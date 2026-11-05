# Start the workbench from this folder (used by launch_mac.command and launch_windows.bat).
# Optional argument: the folder whose files the demo button fills in (default demo_data/).
args <- commandArgs()
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script)) normalizePath(file.path(dirname(script), "..")) else getwd()
demo <- commandArgs(trailingOnly = TRUE)
demo <- if (length(demo) >= 1 && nzchar(demo[1])) demo[1] else "demo_data"
if (!grepl("^(/|[A-Za-z]:)", demo)) demo <- file.path(root, demo)
if (!dir.exists(demo)) stop("Demo folder not found: ", demo)
if (!requireNamespace("labproteomics", quietly = TRUE)) {
  message("labproteomics is not installed; loading it from source with pkgload.")
  pkgload::load_all(root, quiet = TRUE)
}
message("Demo folder: ", demo)
labproteomics::run_app(demo_dir = demo, runs_dir = file.path(root, "runs"), launch.browser = TRUE)
