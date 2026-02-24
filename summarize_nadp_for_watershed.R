library(macrosheds)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(httr)

sites <- c('ALBION', 'GREEN4', 'NAVAJO', 'w8', 'w9', 'W-9')

sheds <- ms_load_spatial_product(
  macrosheds_root = '~/ssd2/macrosheds_stuff/ms_test/',
  spatial_product = 'ws_boundary',
  site_codes = sites
)

# ---- configuration ----

nadp_dir <- 'nadp'
dir.create(nadp_dir, showWarnings = FALSE)

# NADP NTN annual raster base URL
# Variables: pH (lab pH), plus major ion concentrations (mg/L)
# Conductivity is not provided directly, but can be estimated from ions.
# Available from NADP's gridded data: https://nadp.slh.wisc.edu/maps-data/ntn-gradient-maps/

# Ion grids are precipitation-weighted mean concentrations in mg/L
nadp_vars <- c('pH', 'Ca', 'Mg', 'K', 'Na', 'Cl', 'NO3', 'NH4', 'SO4')
years <- 1985:2024  # adjust range as needed

# Equivalent weights (mg/meq) for converting mg/L to ueq/L
equiv_weights <- c(Ca = 20.04, Mg = 12.15, K = 39.10, Na = 22.99,
                   Cl = 35.45, NO3 = 62.00, NH4 = 18.04, SO4 = 48.03)

# Limiting equivalent conductances at 25C (uS/cm per ueq/L)
# Sources: CRC Handbook / Standard Methods
equiv_conductances <- c(Ca = 59.5, Mg = 53.1, K = 73.5, Na = 50.1,
                        Cl = 76.3, NO3 = 71.4, NH4 = 73.5, SO4 = 80.0)
# Note: H+ contributes too; will estimate from pH
H_conductance <- 349.8

# Base URL pattern for NADP NTN concentration grids (annual)
# URL pattern: https://nadp.slh.wisc.edu/filelib/maps/NTN/grids/{year}/{var}_conc_{year}.zip
# Zips contain raster data (ArcGrid or GeoTIFF)
base_url <- 'https://nadp.slh.wisc.edu/filelib/maps/NTN/grids'

# NADP uses these filename prefixes for concentration grids
nadp_var_prefixes <- c(pH  = 'pH',
                       Ca  = 'Ca',
                       Mg  = 'Mg',
                       K   = 'K',
                       Na  = 'Na',
                       Cl  = 'Cl',
                       NO3 = 'NO3',
                       NH4 = 'NH4',
                       SO4 = 'SO4')

# ---- prepare watershed geometries ----

# Reproject watersheds to WGS84 for consistency with NADP rasters
sheds_wgs84 <- st_transform(sheds, 4326)

# Get a bounding box with some buffer for cropping
sheds_bbox <- st_bbox(sheds_wgs84)

# ---- download NADP rasters ----

download_nadp_raster <- function(variable, year, dest_dir) {

    # Check if we already have an extracted raster for this variable/year
    extracted_dir <- file.path(dest_dir, paste0(variable, '_conc_', year))
    if(dir.exists(extracted_dir)) {
        rast_file <- find_raster_in_dir(extracted_dir)
        if(! is.null(rast_file)) {
            message('Already extracted: ', variable, ' ', year)
            return(rast_file)
        }
    }

    # Also check for a .tif directly
    tif_file <- file.path(dest_dir, paste0(variable, '_conc_', year, '.tif'))
    if(file.exists(tif_file)) {
        test <- try(rast(tif_file), silent = TRUE)
        if(! inherits(test, 'try-error')) {
            message('Already have: ', variable, ' ', year)
            return(tif_file)
        }
    }

    prefix <- nadp_var_prefixes[variable]
    zip_name <- paste0(prefix, '_conc_', year, '.zip')
    zip_file <- file.path(dest_dir, zip_name)

    # Try URL patterns (case variations in prefix)
    prefixes_to_try <- unique(c(prefix, tolower(prefix), toupper(prefix)))
    urls_to_try <- unlist(lapply(prefixes_to_try, function(p) {
        zn <- paste0(p, '_conc_', year, '.zip')
        c(
            paste0(base_url, '/', year, '/', zn),
            paste0(base_url, '/', year, '/', tolower(zn))
        )
    }))

    max_retries <- 4

    for(url in urls_to_try) {
        for(attempt in seq_len(max_retries)) {
            resp <- try(GET(url, write_disk(zip_file, overwrite = TRUE),
                            timeout(120)), silent = TRUE)

            if(! inherits(resp, 'try-error') && status_code(resp) == 200) {
                # Try to unzip
                exdir <- file.path(dest_dir, paste0(variable, '_conc_', year))
                unzipped <- try(unzip(zip_file, exdir = exdir), silent = TRUE)

                if(! inherits(unzipped, 'try-error') && length(unzipped) > 0) {
                    rast_file <- find_raster_in_dir(exdir)
                    if(! is.null(rast_file)) {
                        message('Downloaded and extracted: ', variable, ' ', year)
                        file.remove(zip_file)
                        Sys.sleep(2)  # polite delay between successful downloads
                        return(rast_file)
                    }
                }

                # Clean up failed extraction
                if(dir.exists(exdir)) unlink(exdir, recursive = TRUE)
                break  # got 200 but extraction failed; try next URL
            }

            if(file.exists(zip_file)) file.remove(zip_file)

            # If connection error, retry with backoff
            if(inherits(resp, 'try-error') ||
               status_code(resp) %in% c(429, 500, 502, 503, 504)) {
                wait <- 2^attempt + runif(1, 0, 1)
                message('  Retry ', attempt, '/', max_retries,
                        ' for ', variable, ' ', year,
                        ' (waiting ', round(wait, 1), 's)')
                Sys.sleep(wait)
            } else {
                break  # non-retryable HTTP error (e.g. 404); try next URL
            }
        }

        if(file.exists(zip_file)) file.remove(zip_file)
    }

    message('Could not download: ', variable, ' ', year)
    return(NA_character_)
}

find_raster_in_dir <- function(dir_path) {
    # Look for GeoTIFF first
    tifs <- list.files(dir_path, pattern = '\\.tif$', full.names = TRUE,
                       recursive = TRUE)
    if(length(tifs) > 0) return(tifs[1])

    # Look for ArcGrid (hdr.adf)
    adfs <- list.files(dir_path, pattern = 'hdr\\.adf$', full.names = TRUE,
                       recursive = TRUE, ignore.case = TRUE)
    if(length(adfs) > 0) {
        # Return the directory containing hdr.adf (that's the ArcGrid "file")
        return(dirname(adfs[1]))
    }

    # Look for .asc (Arc ASCII grid)
    ascs <- list.files(dir_path, pattern = '\\.asc$', full.names = TRUE,
                       recursive = TRUE)
    if(length(ascs) > 0) return(ascs[1])

    # Look for .img
    imgs <- list.files(dir_path, pattern = '\\.img$', full.names = TRUE,
                       recursive = TRUE)
    if(length(imgs) > 0) return(imgs[1])

    # Try to load anything terra can read
    all_files <- list.files(dir_path, full.names = TRUE, recursive = TRUE)
    for(f in all_files) {
        test <- try(rast(f), silent = TRUE)
        if(! inherits(test, 'try-error')) return(f)
    }

    return(NULL)
}

message('Downloading NADP NTN annual rasters...')
message('If automatic download fails, manually download rasters from:')
message('  https://nadp.slh.wisc.edu/filelib/maps/NTN/grids/{year}/')
message('Place zip files or extracted folders in the nadp/ directory')

download_log <- expand.grid(variable = nadp_vars, year = years,
                            stringsAsFactors = FALSE) %>%
    as_tibble()

download_log$filepath <- mapply(
    download_nadp_raster,
    variable = download_log$variable,
    year = download_log$year,
    dest_dir = nadp_dir
)

download_log <- download_log %>%
    filter(! is.na(filepath))

if(nrow(download_log) == 0) {
    stop('No NADP rasters were downloaded. Please download manually from:\n',
         '  https://nadp.slh.wisc.edu/filelib/maps/NTN/grids/{year}/\n',
         'and place them in nadp/ (as zips or extracted folders)\n',
         'Then re-run this script.')
}

message(nrow(download_log), ' rasters available for processing.')

# ---- also check for any manually placed rasters ----

# Check for .tif files and extracted directories
manual_tifs <- list.files(nadp_dir, pattern = '\\.tif$', full.names = TRUE)
manual_dirs <- list.dirs(nadp_dir, recursive = FALSE, full.names = TRUE)

manual_entries <- list()

if(length(manual_tifs) > 0) {
    manual_entries[['tifs']] <- tibble(filepath = manual_tifs) %>%
        mutate(
            basename = tools::file_path_sans_ext(basename(filepath)),
            variable = sub('_(conc_)?[0-9]{4}$', '', basename),
            year = as.integer(sub('.*_', '', basename))
        ) %>%
        filter(variable %in% nadp_vars, ! is.na(year)) %>%
        select(variable, year, filepath)
}

if(length(manual_dirs) > 0) {
    for(d in manual_dirs) {
        dname <- basename(d)
        # Parse variable_conc_year pattern
        yr <- as.integer(sub('.*_', '', dname))
        vr <- sub('_conc_[0-9]{4}$', '', dname)
        if(! is.na(yr) && vr %in% nadp_vars) {
            rast_file <- find_raster_in_dir(d)
            if(! is.null(rast_file)) {
                manual_entries[[d]] <- tibble(variable = vr, year = yr,
                                             filepath = rast_file)
            }
        }
    }
}

if(length(manual_entries) > 0) {
    manual_info <- bind_rows(manual_entries)
    download_log <- bind_rows(download_log, manual_info) %>%
        distinct(variable, year, .keep_all = TRUE)
}

# ---- extract raster values for each watershed ----

# Use exactextractr if available, otherwise terra::extract
use_exactextractr <- requireNamespace('exactextractr', quietly = TRUE)

if(use_exactextractr) {
    library(exactextractr)
    message('Using exactextractr for weighted extraction.')
} else {
    message('Using terra::extract. Install exactextractr for better accuracy.')
}

results <- list()

for(i in seq_len(nrow(download_log))) {

    row <- download_log[i, ]
    message('Processing: ', row$variable, ' ', row$year)

    r <- try(rast(row$filepath), silent = TRUE)
    if(inherits(r, 'try-error')) {
        message('  Skipping (invalid raster)')
        next
    }

    # Reproject watersheds to match raster CRS
    sheds_reproj <- st_transform(sheds_wgs84, crs(r))

    # Crop to reduce memory
    sheds_ext <- ext(vect(sheds_reproj))
    sheds_ext_buf <- sheds_ext + 0.5
    r_crop <- try(crop(r, sheds_ext_buf), silent = TRUE)
    if(inherits(r_crop, 'try-error')) r_crop <- r

    if(use_exactextractr) {
        vals <- exact_extract(r_crop, sheds_reproj, fun = 'mean')
    } else {
        ex <- terra::extract(r_crop, vect(sheds_reproj), fun = mean,
                             na.rm = TRUE, weights = TRUE)
        vals <- ex[, 2]
    }

    res <- tibble(
        site_code = sheds_wgs84$site_code,
        variable = row$variable,
        year = row$year,
        value = vals
    )

    results[[length(results) + 1]] <- res
}

# ---- combine and format results ----

nadp_summary <- bind_rows(results)

# Add descriptive labels
nadp_summary <- nadp_summary %>%
    mutate(
        variable_description = case_when(
            variable == 'pH'   ~ 'Precipitation-weighted mean pH',
            variable == 'Ca'   ~ 'Precipitation-weighted mean Ca (mg/L)',
            variable == 'Mg'   ~ 'Precipitation-weighted mean Mg (mg/L)',
            variable == 'K'    ~ 'Precipitation-weighted mean K (mg/L)',
            variable == 'Na'   ~ 'Precipitation-weighted mean Na (mg/L)',
            variable == 'Cl'   ~ 'Precipitation-weighted mean Cl (mg/L)',
            variable == 'NO3'  ~ 'Precipitation-weighted mean NO3 (mg/L)',
            variable == 'NH4'  ~ 'Precipitation-weighted mean NH4 (mg/L)',
            variable == 'SO4'  ~ 'Precipitation-weighted mean SO4 (mg/L)',
            TRUE ~ variable
        )
    )

# Pivot to wide format (one column per variable)
nadp_wide <- nadp_summary %>%
    select(site_code, year, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    arrange(site_code, year)

# ---- estimate conductivity from ions and pH ----

ion_cols <- intersect(names(equiv_weights), names(nadp_wide))

if(length(ion_cols) > 0 && 'pH' %in% names(nadp_wide)) {

    nadp_wide <- nadp_wide %>%
        mutate(
            # Convert each ion from mg/L to ueq/L, then multiply by equiv conductance
            Cond_est = {
                cond <- rep(0, n())
                for(ion in ion_cols) {
                    conc_mgl <- .data[[ion]]
                    conc_ueql <- conc_mgl / equiv_weights[ion] * 1000
                    cond <- cond + conc_ueql * equiv_conductances[ion] / 1000
                }
                # Add H+ contribution from pH
                H_ueql <- 10^(-pH) * 1e6  # mol/L -> ueq/L
                cond <- cond + H_ueql * H_conductance / 1000
                round(cond, 2)
            }
        )

    message('Estimated conductivity (Cond_est, uS/cm) computed from ion concentrations and pH.')

    # Also add to long format
    cond_long <- nadp_wide %>%
        select(site_code, year, Cond_est) %>%
        mutate(variable = 'Cond_est',
               variable_description = 'Estimated conductivity from ions + pH (uS/cm)') %>%
        rename(value = Cond_est)

    nadp_summary <- bind_rows(nadp_summary, cond_long)

} else {
    message('Not enough ion variables to estimate conductivity.')
}

# ---- save results ----

write.csv(nadp_summary, 'nadp_watershed_summary_long.csv', row.names = FALSE)
write.csv(nadp_wide, 'nadp_watershed_summary_wide.csv', row.names = FALSE)

message('\n---- Summary ----')
message('Sites: ', paste(unique(nadp_summary$site_code), collapse = ', '))
message('Variables: ', paste(sort(unique(nadp_summary$variable)), collapse = ', '))
message('Years: ', min(nadp_summary$year), '-', max(nadp_summary$year))
message('Results saved to:')
message('  nadp_watershed_summary_long.csv')
message('  nadp_watershed_summary_wide.csv')

# Print a quick summary
nadp_summary %>%
    group_by(site_code, variable) %>%
    summarize(
        n_years = n(),
        mean_value = mean(value, na.rm = TRUE),
        min_value = min(value, na.rm = TRUE),
        max_value = max(value, na.rm = TRUE),
        .groups = 'drop'
    ) %>%
    print(n = Inf)

