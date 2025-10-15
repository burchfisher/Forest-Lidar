# install.packages(c("lidR","terra"))  # if needed
suppressPackageStartupMessages({
  library(lidR)
  library(terra)
})

# ------------------------------------------------------------------------
# 🏔️  LIDAR DTM BUILDER (PER TILE) — unit-aware, Z→meters, EPSG-tagged outputs
# Created by Burch Fisher 10/15/25
# ------------------------------------------------------------------------
# What this script does
# 1) Reads a folder of LAS/LAZ tiles into a LAScatalog (no '*.laz' glob needed).
# 2) Filters to **class 2 (ground)** points only.
# 3) Processes **per source file** (chunk_size = 0) with a **buffer** for smooth
#    interpolation across tile edges.
# 4) Converts user-supplied spatial parameters (`res_m`, `buffer_m`) from meters
#    into the dataset’s **native XY units** (meters or US-survey-feet).
# 5) **Automatically skips tiles/chunks with no ground points** (e.g., all water)
#    so the run never crashes.
# 6) Builds per-tile DTMs using **TIN interpolation**, automatically retrying with
#    **IDW** if TIN fails for a given tile/chunk.
# 7) If the LAS CRS is in **US-survey-feet**, converts all Z (elevation) values to
#    **meters** before interpolation, so elevations are metric while XY remain native.
# 8) Optionally **reprojects** each finished DTM to a fully metric CRS you specify
#    in CONFIG (so both XY and Z end up in meters).
# 9) **Always appends an EPSG tag** to the output filenames:
#       "<tile>_DTM_3m_epsg2248.tif"  (native)
#       "<tile>_DTM_3m_epsg26917.tif" (reprojected)
#
# Units summary
# • Inputs (`res_m`, `buffer_m`) are given in meters and converted to native XY.
# • Elevations (Z) are meters if source CRS is US-ft (by default), otherwise native (see "force_z_m" varible).
# • If `output_crs` is set to a metric EPSG, both XY and Z will be metric in outputs.
# ------------------------------------------------------------------------

# =========================
# CONFIG (you set these in METERS)
# =========================
laz_dir    <- "/path/to/input/LAZ/"
out_dir    <- "/path/to/output/DTM"
res_m      <- 3           # pixel size (meters)
buffer_m   <- 50          # tile buffer (meters); 50–100 m recommended in steep terrain

# Reproject outputs? Set to NULL to keep native CRS, or to an EPSG string (e.g., "EPSG:26917")
output_crs <- NULL

# Force Z to meters when source CRS is US-ft?
force_z_m  <- TRUE

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# Helpers: unit detection & conversion
# =========================
USFT_TO_M <- 1200/3937         # 0.3048006096 m
M_TO_USFT <- 1 / USFT_TO_M

.get_proj_text <- function(x) {
  p <- tryCatch(projection(x), error = function(e) NULL)
  if (is.null(p)) "" else as.character(p)
}

detect_us_survey_foot <- function(x) {
  w <- .get_proj_text(x)
  if (!nzchar(w)) return(NA)
  if (grepl("US.?survey.?foot|Foot_US|US[-_ ]?ft|US[-_ ]?Survey[-_ ]?Foot", w, ignore.case = TRUE)) return(TRUE)
  if (grepl("metre|meter|\\bunit\\s*=\\s*metre", w, ignore.case = TRUE)) return(FALSE)
  NA
}

meters_to_native <- function(x_meters, las_or_ctg) {
  det <- detect_us_survey_foot(las_or_ctg)
  is_usft <- if (is.na(det)) FALSE else det
  if (is_usft) x_meters * M_TO_USFT else x_meters
}

# =========================
# Catalog setup (per-tile, buffered)
# =========================
ctg <- readLAScatalog(laz_dir)

opt_filter(ctg)           <- "-keep_class 2"   # ground only
opt_chunk_size(ctg)       <- 0                 # process per SOURCE FILE
opt_chunk_buffer(ctg)     <- meters_to_native(buffer_m, ctg)
opt_laz_compression(ctg)  <- TRUE
opt_progress(ctg)         <- TRUE

# Convert resolution (meters → native XY)
res_native <- meters_to_native(res_m, ctg)

# Warn if CRS missing
proj_txt <- .get_proj_text(ctg)
if (!nzchar(proj_txt)) {
  message("⚠️  Catalog has no CRS. Assuming meters for XY units. ",
          "If data are actually in US survey feet, set a CRS first.")
}

# ----------------------------
# Output template (always include an EPSG tag)
# ----------------------------
get_epsg_tag <- function(x) {
  tag <- tryCatch({
    code <- suppressWarnings(lidR::epsg(x))
    if (is.na(code)) "epsgNA" else paste0("epsg", code)
  }, error = function(e) "epsgNA")
  tag
}

src_epsg_tag <- get_epsg_tag(ctg)

# If reprojecting, prefer target EPSG in the filename tag
target_epsg_tag <- if (!is.null(output_crs) && grepl("^EPSG:", toupper(output_crs))) {
  tolower(gsub("EPSG:", "epsg", toupper(output_crs)))
} else {
  src_epsg_tag
}

opt_output_files(ctg) <- file.path(out_dir, sprintf("{*}_DTM_%dm_%s", res_m, target_epsg_tag))

# =========================
# Per-chunk DTM builder
# =========================
build_dtm_chunk <- function(las, res_native, force_z_m = TRUE, output_crs = NULL) {
  
  # Skip chunks with no ground points (catalog is already filtered to class 2)
  if (npoints(las) == 0) {
    e <- tryCatch(lidR::lasextent(las), error = function(...) NULL)
    if (!is.null(e)) {
      message(sprintf("⛵ No ground points in chunk (%.1f,%.1f)-(%.1f,%.1f). Skipping.",
                      e@xmin, e@ymin, e@xmax, e@ymax))
    } else {
      message("⛵ No ground points in chunk. Skipping.")
    }
    return(NULL)
  }
  
  # Convert Z to meters if source CRS uses US-survey-feet
  is_usft <- detect_us_survey_foot(las)
  if (isTRUE(is_usft) && isTRUE(force_z_m)) {
    message("↕️  Converting Z from US-survey-feet to meters for this chunk.")
    las@data$Z <- las@data$Z * USFT_TO_M
  }
  
  # Build DTM (TIN → fallback to IDW)
  r <- tryCatch(
    rasterize_terrain(las, res = res_native, algorithm = tin(), na.rm = TRUE),
    error = function(e) {
      message("TIN failed: ", conditionMessage(e), " → retrying with IDW")
      rasterize_terrain(las, res = res_native,
                        algorithm = knnidw(k = 10, p = 2), na.rm = TRUE)
    }
  )
  
  if (!inherits(r, "SpatRaster")) {
    message("⚠️  Could not build DTM even with IDW. Skipping.")
    return(NULL)
  }
  
  # Optional reprojection to metric (or any) target CRS
  if (!is.null(output_crs)) {
    r <- terra::project(r, output_crs, method = "bilinear")
  }
  
  r
}

# =========================
# Run (writes one DTM per source tile)
# =========================
invisible(
  catalog_map(
    ctg,
    build_dtm_chunk,
    res_native = res_native,
    force_z_m  = force_z_m,   # TRUE => Z in meters when source CRS is US-ft
    output_crs = output_crs   # NULL => keep native CRS; else reproject before writing
  )
)

# =========================
# Summary message
# =========================
if (is.null(output_crs)) {
  message(sprintf(
    "✅ Done. Per-tile DTMs written to: %s\n   Cell size: %.2f (native XY units)\n   Elevations (Z): meters if source CRS was US-ft; otherwise native.\n   Filenames include EPSG tag: %s",
    normalizePath(out_dir), res_native, src_epsg_tag
  ))
} else {
  message(sprintf(
    "✅ Done. Per-tile DTMs written to: %s\n   Reprojected to: %s (%s)\n   Cell size: %.2f (target CRS units)\n   Elevations (Z): meters.\n   Filenames include EPSG tag: %s",
    normalizePath(out_dir), output_crs, target_epsg_tag, res_m, target_epsg_tag
  ))
}
