# Start the workbench from this folder (used by launch_mac.command and launch_windows.bat).
args <- commandArgs()
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script)) normalizePath(file.path(dirname(script), "..")) else getwd()
if (!requireNamespace("labproteomics", quietly = TRUE)) {
  message("labproteomics is not installed; loading it from source with pkgload.")
  pkgload::load_all(root, quiet = TRUE)
  run_app <- labproteomics::run_app
}
labproteomics::run_app(demo_dir = file.path(root, "demo_data"), runs_dir = file.path(root, "runs"),
                       launch.browser = TRUE)
