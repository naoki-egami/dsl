.onAttach <- function(libname, pkgname) {
  pd <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pd$Version, " successfully loaded. See ?", pkgname, " for help. Note this is an early alpha release and backwards compatability may not be maintained.");
}
