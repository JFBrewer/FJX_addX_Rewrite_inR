# FAST-JX addX rewrite, version 5
# Reference implementation: Updating_FJX/UCI_fastJ_addX_73cx/addX_ACET.f
#
# PURPOSE
# -------
# Convert a high-resolution, wavelength-dependent photolysis table (FTBL) into
# the 18 FAST-JX photolysis bins used by GEOS-Chem / FAST-JX.  The calculation
# deliberately follows the order of the Fortran program:
#
#   1. Read the 77 Pratmo wavelength bins, 18-bin map, solar flux, and SRB.
#   2. Evaluate XNEW = cross section x quantum yield at every 0.05-nm point.
#   3. Flux-weight those values into 77 Pratmo bins.
#   4. Flux-weight the 77 Pratmo bins into 18 FAST-JX bins.
#   5. Optionally write the result in GEOS-Chem or FAST-JX fixed-width form.
#
# All required input files are passed explicitly: this script makes no hidden
# assumptions about the active species, except that FTBL has the standard
# FAST-JX three-case table layout.

source("FJX_Subroutines.R")

addx_v5 <- function(FTBL, FJX_Dir, output = NULL) {
  # `output` is intentionally NULL by default: callers can use the returned
  # tibble without creating a file.  When it is a path, GC_Format chooses the
  # legacy GEOS-Chem layout (TRUE) or FAST-JX-v73 a/b/c layout (FALSE).
  output_path <- output
  # Permit either a full table path or the table filename within FJX_Dir.
  if (!file.exists(FTBL)) {
    candidate <- file.path(FJX_Dir, FTBL)
    if (!file.exists(candidate)) stop("Cannot find FTBL: ", FTBL)
    FTBL <- candidate
  }
  # The first section of wavel-bins.dat defines the 77 parent (Pratmo) bins.
  # Do not read the historical 78th catch-all row: the Fortran code sets NB=77.
  PratmoBins <- readr::read_table(
    file.path(FJX_Dir, "wavel-bins.dat"), skip = 1, n_max = 77,
    col_names = c("PratmoBin", "Lambda_Start", "Lambda_End"),
    show_col_types = FALSE
  )
  # The last section maps Pratmo bins 16:77 directly to FAST-JX bins 5:18.
  JXBin_Num <- readr::read_table(
    file.path(FJX_Dir, "wavel-bins.dat"), skip = 90,
    col_names = c("PratmoBin", "FJX_Bin"), show_col_types = FALSE
  )
  # Flux is already integrated over each 0.05-nm solar microbin, exactly as
  # assumed by the original FBIN = FBIN + F(J) accumulation.
  SUSIM_Spectra <- readr::read_table(
    file.path(FJX_Dir, "solar-p05nm-UCI.dat"), skip = 2,
    col_names = c("SUSIM_Wavelength", "Flux"), show_col_types = FALSE
  )
  # SRB is the Schumann-Runge opacity-distribution matrix.  Its CSV rows are
  # FAST-JX bins and its columns are Pratmo bins 1:15.
  SRB <- readr::read_csv(
    file.path(FJX_Dir, "Bin_Assignment.csv"), skip = 1,
    show_col_types = FALSE
  )[, -1]

  # Keep these stages separate so they can be inspected independently when a
  # new cross-section table needs validation against a Fortran diagnostic run.
  xnew <- X_MeVK(SUSIM_Spectra, FTBL)
  pratmo <- Do_PrAtmo_Binning(xnew, PratmoBins)
  output <- Do_FJX_Binning(pratmo, JXBin_Num, SRB)
  output <- dplyr::left_join(output, read_fjx_pt(FTBL), by = "Altitude (km)")
  if (!is.null(output_path)) {
    write_geoschem_spec(output, FTBL, output_path)
  }
  output
}
