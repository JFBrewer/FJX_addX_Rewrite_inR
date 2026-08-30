# Direct, Fortran-faithful helpers for FAST-JX addX calculations.
#
# DATA FLOW
# ---------
# FTBL -> XNEW at each solar microbin -> 77 Pratmo means -> 18 FAST-JX means
#
# The names below intentionally separate the two kinds of bins:
# * "Pratmo" means one of the 77 wavelength intervals in wavel-bins.dat.
# * "FJX" means one of the final 18 FAST-JX photolysis bins.
# * "XNEW" is cross section x quantum yield in cm2, evaluated at a solar
#   wavelength and one of the three pressure/temperature cases.
#
# These functions avoid fuzzy joins.  A fuzzy join materializes many possible
# wavelength/bin pairs; the Fortran logic instead assigns each microbin once.

read_fjx_pt <- function(FTBL) {
  # Lines 4:6 in a standard FTBL define the three requested atmospheric cases.
  # [M] is stored in units of 1e19 molecules cm-3 and is converted here to the
  # physical units used by the Fortran program after its initialization call.
  raw <- suppressWarnings(readr::read_table(
    FTBL, skip = 3, n_max = 3,
    col_names = c("Variable", "13", "5", "0"), show_col_types = FALSE
  ))
  values <- as.matrix(raw[, 2:4])
  storage.mode(values) <- "double"

  tibble::tibble(
    `Altitude (km)` = as.numeric(colnames(values)),
    `p/hPa` = values[1, ],
    `[M]` = values[2, ] * 1e19,
    `T/K` = values[3, ]
  )
}

read_fjx_xs_table <- function(FTBL) {
  # Lines 13 onward contain: three interpolation temperatures, an NW/format
  # record, and then NW wavelength / cross-section / three-QY records.
  temperatures <- as.numeric(suppressWarnings(readr::read_table(
    FTBL, skip = 12, n_max = 1, col_names = c("13", "5", "0"),
    show_col_types = FALSE
  )))
  header <- suppressWarnings(readr::read_table(
    FTBL, skip = 13, n_max = 1, col_names = c("NW", "Format"),
    show_col_types = FALSE
  ))
  table <- suppressWarnings(readr::read_table(
    FTBL, skip = 14, n_max = header$NW,
    col_names = c("Wavelength", "Cross_Section", "QY_1", "QY_2", "QY_3"),
    show_col_types = FALSE
  ))

  list(temperatures = temperatures, table = table)
}

# Equivalent to the non-initialization branch of X_MeVK in addX_ACET.f.
# Unlike the prior R code, this returns XNEW for *all* solar grid points.
X_MeVK <- function(SUSIM_Spectra, FTBL) {
  # This function has a historical name because it was first validated against
  # MeVK.  It is generic: FTBL determines the species and quantum yields.
  pt <- read_fjx_pt(FTBL)
  xs <- read_fjx_xs_table(FTBL)
  reference <- xs$table
  qt3 <- xs$temperatures
  nw <- nrow(reference)
  if (nw < 2) stop("The cross-section table must contain at least two wavelengths.")

  # Fortran chooses W(I) < WW <= W(I+1), with IW limited to 1:(NW-1).
  # left.open=TRUE reproduces the strict lower / inclusive upper boundary.
  # Clamping IW and FW is important: unlike the earlier R rewrite, Fortran
  # evaluates every solar wavelength and holds the closest endpoint outside
  # the tabulated wavelength range.
  iw <- findInterval(SUSIM_Spectra$SUSIM_Wavelength,
                     reference$Wavelength, left.open = TRUE)
  iw <- pmin(nw - 1L, pmax(1L, iw))
  lower_w <- reference$Wavelength[iw]
  upper_w <- reference$Wavelength[iw + 1L]
  fw <- (SUSIM_Spectra$SUSIM_Wavelength - lower_w) / (upper_w - lower_w)
  fw <- pmin(1, pmax(0, fw))

  # Interpolate cross section once because it does not vary among the three
  # atmospheric cases.  Interpolate all three QY columns before choosing the
  # temperature pair needed for each case.
  xs_interp <- reference$Cross_Section[iw] + fw *
    (reference$Cross_Section[iw + 1L] - reference$Cross_Section[iw])
  qy_interp <- vapply(1:3, function(k) {
    qy <- reference[[paste0("QY_", k)]]
    qy[iw] + fw * (qy[iw + 1L] - qy[iw])
  }, numeric(nrow(SUSIM_Spectra)))

  purrr::map_dfr(seq_len(nrow(pt)), function(k) {
    # This is the Fortran XXT / IT / QFACT temperature interpolation.  The
    # target temperature is clamped to the table range, then QY is linearly
    # interpolated only within the appropriate adjacent temperature pair.
    xxt <- min(qt3[3], max(qt3[1], pt$`T/K`[k]))
    it <- if (xxt > qt3[2]) 2L else 1L
    qfact <- (xxt - qt3[it]) / (qt3[it + 1L] - qt3[it])
    qy <- qy_interp[, it] + qfact * (qy_interp[, it + 1L] - qy_interp[, it])

    tibble::tibble(
      SUSIM_Wavelength = SUSIM_Spectra$SUSIM_Wavelength,
      Flux = SUSIM_Spectra$Flux,
      `Altitude (km)` = pt$`Altitude (km)`[k],
      `p/hPa` = pt$`p/hPa`[k],
      `[M]` = pt$`[M]`[k],
      `T/K` = pt$`T/K`[k],
      XNEW = xs_interp * qy * 1e-20
    )
  })
}

assign_pratmo_bins <- function(wavelength, PratmoBins) {
  # Assign every solar wavelength to at most one Pratmo bin.  This matches
  # the nested Fortran searches: W(J) > WBIN(I) .and. W(J) <= WBIN(I+1).
  # A zero denotes a wavelength outside the 77 FAST-JX parent bins.
  bin <- findInterval(wavelength, PratmoBins$Lambda_Start, left.open = TRUE)
  valid <- bin >= 1L & bin <= nrow(PratmoBins)
  valid[valid] <- wavelength[valid] <= PratmoBins$Lambda_End[bin[valid]]
  bin[!valid] <- 0L
  bin
}

Do_PrAtmo_Binning <- function(Integrated_XSQY, PratmoBins) {
  # For every atmospheric case calculate the two Fortran work arrays:
  #   FBIN(I) = sum(F(J))
  #   ABIN(I) = sum(F(J) * XNEW(J)) / FBIN(I)
  # Here they are returned as PratmoFlux and PratmoJ, respectively.
  if (nrow(PratmoBins) != 77L) stop("FAST-JX requires exactly 77 Pratmo bins.")
  bin <- assign_pratmo_bins(Integrated_XSQY$SUSIM_Wavelength, PratmoBins)

  purrr::map_dfr(sort(unique(Integrated_XSQY$`Altitude (km)`)), function(altitude) {
    d <- Integrated_XSQY[Integrated_XSQY$`Altitude (km)` == altitude & bin > 0L, ]
    d_bin <- bin[Integrated_XSQY$`Altitude (km)` == altitude & bin > 0L]
    flux <- numeric(77)
    numerator <- numeric(77)
    # rowsum is the vectorized equivalent of the Fortran microbin loop.  The
    # explicit 77-element vectors retain zero-flux bins and their bin numbers.
    flux_values <- rowsum(d$Flux, d_bin, reorder = FALSE)
    numerator_values <- rowsum(d$Flux * d$XNEW, d_bin, reorder = FALSE)
    index <- as.integer(rownames(flux_values))
    flux[index] <- flux_values[, 1]
    numerator[index] <- numerator_values[, 1]

    tibble::tibble(
      `Altitude (km)` = altitude,
      PratmoBin = seq_len(77),
      PratmoFlux = flux,
      PratmoJ = ifelse(flux > 0, numerator / flux, 0)
    )
  })
}

Do_FJX_Binning <- function(Pratmo_Summary, JXBin_Num, SRB_raw) {
  # This is the second aggregation stage.  Bins 16:77 are mapped one-to-one
  # by IJX; the UV Schumann-Runge bins 1:15 are split among FAST-JX bins using
  # the SRB opacity-distribution weights.
  if (nrow(SRB_raw) > 18L || ncol(SRB_raw) != 15L) {
    stop("SRB must be stored as Fast-JX bins (rows) by Pratmo bins (15 columns).")
  }
  # Input CSV is rows=J and columns=I; Fortran accesses SRB(I,J).  Transpose
  # it once into the Fortran orientation and replace blank CSV cells by zero.
  SRB <- matrix(0, nrow = 15, ncol = 18)
  SRB[, seq_len(nrow(SRB_raw))] <- t(as.matrix(SRB_raw))
  SRB[is.na(SRB)] <- 0
  j_for_i <- integer(77)
  j_for_i[JXBin_Num$PratmoBin] <- JXBin_Num$FJX_Bin

  purrr::map_dfr(sort(unique(Pratmo_Summary$`Altitude (km)`)), function(altitude) {
    d <- Pratmo_Summary[Pratmo_Summary$`Altitude (km)` == altitude, ]
    d <- d[order(d$PratmoBin), ]
    flux <- d$PratmoFlux
    abin <- d$PratmoJ
    ffbin <- numeric(18)
    aabin <- numeric(18)

    # Direct mapping for the non-Schumann-Runge wavelength bins (I = 16:77).
    for (i in 16:77) {
      j <- j_for_i[i]
      ffbin[j] <- ffbin[j] + flux[i]
      aabin[j] <- aabin[j] + flux[i] * abin[i]
    }
    # Opacity-distribution split for the short-wavelength bins (I = 1:15).
    for (i in 1:15) {
      for (j in 1:18) {
        ffbin[j] <- ffbin[j] + flux[i] * SRB[i, j]
        aabin[j] <- aabin[j] + flux[i] * abin[i] * SRB[i, j]
      }
    }

    tibble::tibble(
      `Altitude (km)` = altitude,
      FJX_Bin = seq_len(18),
      FJX_Flux = ffbin,
      FJX_J = ifelse(ffbin > 0, aabin / ffbin, 0)
    )
  })
}

read_fjx_titles <- function(FTBL) {
  # The three 90-character titles immediately follow the T/K line.  FAST-JX
  # writes one title before each three-line v73 output block.
  lines <- readLines(FTBL, warn = FALSE)
  if (length(lines) < 9L) stop("FTBL does not contain the three FAST-JX title lines.")
  substr(lines[7:9], 1L, 90L)
}

# Write the final 18-bin results in the fixed-width layout produced by
# addX_ACET.f for FJX_spec.dat: title, then a/b/c lines holding six bins each.

# Write the legacy six-values-per-line layout used by GEOS-Chem FJX_spec.dat
# and shown in addX_MeVK.out before the title/a/b/c FAST-JX-v73 records.
write_geoschem_spec <- function(FJX_Output, FTBL, output_file) {
  # Legacy GEOS-Chem format: three 70-character numeric records per case.
  # This is the first output block in addX_MeVK.out (before v73 a/b/c records).
  required <- c("Altitude (km)", "FJX_Bin", "FJX_J", "p/hPa")
  if (!all(required %in% names(FJX_Output))) {
    stop("FJX_Output must contain: ", paste(required, collapse = ", "))
  }
  if (!all(table(FJX_Output$`Altitude (km)`) == 18L)) {
    stop("FJX_Output must contain exactly 18 Fast-JX bins for each case.")
  }

  pt <- read_fjx_pt(FTBL)
  titles <- read_fjx_titles(FTBL)
  species <- sprintf("%-6s", substr(titles[1], 1L, 6L))
  pressure_flag <- substr(readLines(FTBL, warn = FALSE)[3], 1L, 1L)
  if (pressure_flag != "p") pressure_flag <- " "
  format_values <- function(x) paste0(sprintf("%10.3E", x), collapse = "")

  lines <- character(0)
  for (k in seq_len(nrow(pt))) {
    altitude <- pt$`Altitude (km)`[k]
    d <- FJX_Output[FJX_Output$`Altitude (km)` == altitude, ]
    d <- d[order(d$FJX_Bin), ]
    if (!identical(as.integer(d$FJX_Bin), 1:18)) {
      stop("FJX bins must be numbered 1 through 18 for every case.")
    }

    pressure <- as.integer(round(pt$`p/hPa`[k]))
    first_prefix <- sprintf("%s%1s%3d", species, pressure_flag, pressure)
    lines <- c(
      lines,
      paste0(first_prefix, format_values(d$FJX_J[1:6])),
      paste0(strrep(" ", 10L), format_values(d$FJX_J[7:12])),
      paste0(strrep(" ", 10L), format_values(d$FJX_J[13:18])),
      ""
    )
  }
  writeLines(lines, output_file, useBytes = TRUE)
  invisible(normalizePath(output_file, winslash = "/", mustWork = FALSE))
}
