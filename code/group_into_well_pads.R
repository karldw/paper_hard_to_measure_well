suppressMessages(
  here::i_am("code/group_into_well_pads.R", uuid="d22f7910-5d5a-4280-ab7a-2685a1725e10")
)
source(here::here("code/shared_functions.r"))

INPUT_CRS   <- sf::st_crs(CRS_LONGLAT_int) # https://epsg.io/4326
WORKING_CRS <- sf::st_crs(CRS_PROJECT_int) # https://epsg.io/6350



#' Union intersecting geometries in x
#
#' @param x An simple features object (`sf`, `sfg`, or `sfg`)
#' @return The union of intersecting geometries
#' Specifically, union all of the intersecting groups in x
#' https://gis.stackexchange.com/a/323067
#' @examples
#' library(sf)
#' sq = function(pt, sz = 1) {
#'  st_polygon(list(rbind(c(pt - sz), c(pt[1] + sz, pt[2] - sz), c(pt + sz),
#'    c(pt[1] - sz, pt[2] + sz), c(pt - sz))))
#' }
#' x = st_sf(box = 1:6, st_sfc(sq(c(4.2,4.2)), sq(c(0,0)), sq(c(1, -0.8)),
#'   sq(c(0.5, 1.7)), sq(c(3,3)), sq(c(-3, -3))))
#' st_union_intersection(x)
#' @seealso [sf::st_union()]
#' @export
st_union_intersection <- function(x) {
  if (!requireNamespace("igraph", quietly=TRUE)) {
    stop("st_union_intersection requires the igraph package")
  }
  if (!requireNamespace("sf", quietly=TRUE)) {
    stop("st_union_intersection requires the sf package")
  }
  if (length(x) == 0) {
    return(x)
  }
  # Helper function:
  union_by_index <- function(idx, geom) {
    sf::st_union(geom[idx])
  }
  # doesn't matter if x is already a geometry
  x <- sf::st_geometry(x)
  unioned_list <- sf::st_intersects(x) %>%
    igraph::graph_from_adj_list() %>%
    igraph::components() %>%
    igraph::groups() %>%
    lapply(union_by_index, geom=x)
  do.call(c, unioned_list)
}



find_well_pads <- function(df, well_pad_number_start, verbose=FALSE) {
  if (verbose) {
    message("  ", paste0(sort(unique(df$aapg_geologic_province)), collapse=", "))
  }
  df_radius <- df %>%
    # Drop different entity_id at the exact same location
    # (distinct.sf already keeps geometry column by default, but for some
    # reason gives an error)
    dplyr::distinct(aapg_geologic_province, geometry) %>%
    dplyr::mutate(
      buffer_radius_m = dplyr::case_when(
        aapg_geologic_province == "SAN JOAQUIN BASIN" ~ units::as_units(20, "m"),
        TRUE ~ units::as_units(50, "m")
      ))
  well_pad_geometry <- df_radius %>%
    sf::st_geometry() %>%
    sf::st_buffer(dist=df_radius$buffer_radius_m) %>%
    st_union_intersection()
  new_well_pad_ids <- seq.int(
    from = well_pad_number_start,
    length.out = length(well_pad_geometry)
  )
  well_pad_info <- well_pad_geometry %>%
    sf::st_centroid() %>%
    geometry_to_lonlat() %>%
    dplyr::transmute(
      well_pad_lon = .data$longitude,
      well_pad_lat = .data$latitude,
      well_pad_id  = new_well_pad_ids,
    )
  well_pad_df <- sf::st_sf(well_pad_info, geometry = well_pad_geometry)

  out <- dplyr::select(df, entity_id) %>%
    sf::st_join(well_pad_df, left=FALSE) %>%
    sf::st_drop_geometry() %>%
    dplyr::as_tibble()
  out
}


assign_well_pads <- function(header_df, verbose=TRUE) {
  stopifnot(!anyDuplicated(header_df$entity_id))
  # My definition of a well pad, based on Alvarez et al. (2018) SI:
  # - a group of wells that are close together
  # - if the DI record is a lease (irrelevant for CA, NM, and post-1999 CO),
  #   assume it is its own well pad; don't group it with other records
  # - if more than two wells are nearby, take the union. (See the function
  #   st_union_intersection.)
  # NB: It would be nice to use operator, but there are a bunch of wells with
  # different operators reported at the same location.
  # Steps:
  # 1. Make into an sf dataframe
  # 2. Transform to a useful projection
  # 3. Buffer all points 50m, except the San Joaquin basin, where the buffer is 20m.
  # 4. Union the buffers and assign IDs
  # 5. Take the intersection of the original points with the buffers
  # 6. Do the above separately for each state or basin, then re-aggregate.
  # Note that the well-pad definitions are sensitive to the projection and
  # assumed well spacing (eg. 50m) for a small number of wells.

  # For records that are marked as "LEASE" (as opposed to "COM", "DRIP POINT",
  # "SWD", "UNIT", or "WELL"), don't do the geographic work, just label them
  # their own singleton well pad.
  singletons <- header_df %>%
    dplyr::group_by(aapg_geologic_province) %>%
    dplyr::filter(dplyr::n() == 1 | entity_type == "LEASE") %>%
    dplyr::ungroup() %>%
    dplyr::select(entity_id, surface_longitude_wgs84, surface_latitude_wgs84) %>%
    dplyr::mutate(well_pad_id = dplyr::row_number()) %>%
    dplyr::rename(
      well_pad_lon = surface_longitude_wgs84,
      well_pad_lat = surface_latitude_wgs84,
    )

  # For records that aren't LEASE or singleton, look for nearby wells within
  # groups defined by basin and operator_company_name. Note that there are a
  # small number of overlaps, where a well falls within a well pad with
  # different basin or operator_company_name. (This is true even if you only use
  # basin.)
  non_singleton <- header_df %>%
    dplyr::anti_join(singletons, by="entity_id") %>%
    lonlat_to_projected(c("surface_longitude_wgs84", "surface_latitude_wgs84"))

  # well_pads is the geometry of the well pads (a POLYGON buffered around the
  # points of each group of wells)
  non_singleton_lst <- non_singleton %>%
    # Split by basin so we can run different regions in parallel
    # Well pads will not span different subgroups of the group_by, so be careful
    # before adding grouping variables.
    # NB: we are *not* grouping by operator_company_name -- there are some wells
    # that are reported at exactly the same location, but with different
    # operators. We group these into a single well pad.
    dplyr::group_by(aapg_geologic_province) %>%
    dplyr::group_split()
  # We then also need to make sure we're not assigning duplicate well_pad_id
  # values in the different parallel processes, so calculate disjoint sets of
  # possible_well_pad_ids
  n_groups <- length(non_singleton_lst)
  subgroup_row_counts <- purrr::map_int(non_singleton_lst, nrow)
  well_pad_id_bounds <- c(nrow(singletons), subgroup_row_counts) %>% cumsum()
  well_pad_id_starts <- dplyr::lag(well_pad_id_bounds)[-1] + 1
  # possible_well_pad_ids isn't used, just checked
  possible_well_pad_ids <- purrr::map2(
      well_pad_id_starts,
      well_pad_id_bounds[-1],
      ~seq.int(from=.x, to=.y, by=1)
    ) %>%
    unlist() %>%
    c(singletons$well_pad_id)

  stopifnot(
    all(subgroup_row_counts) > 0,
    anyDuplicated(possible_well_pad_ids) == 0,
    length(non_singleton_lst) == length(well_pad_id_starts)
  )

  # Actually do the well pad creation, and re-add the singletons
  entity_pad_crosswalk <- furrr::future_map2_dfr(
      non_singleton_lst, well_pad_id_starts, find_well_pads,
      verbose=verbose,
      # we're not doing any RNG, but suppress warnings
      .options=furrr::furrr_options(seed=TRUE)
    ) %>%
    dplyr::bind_rows(singletons)
  # Note: singletons got assigned well_pad_id values 1:nrow(singletons), but it
  # doesn't matter that they're bound at the end here.

  out <- dplyr::inner_join(header_df, entity_pad_crosswalk, by="entity_id")
  stopifnot(
    anyDuplicated(header_df$entity_id) == 0,
    anyDuplicated(entity_pad_crosswalk$entity_id) == 0,
    setequal(header_df$entity_id, out$entity_id)
  )
  if (verbose) {
    n_wells <- nrow(header_df)
    n_unmatched <- n_wells - nrow(out)
    well_pad_well_counts <- dplyr::count(out, well_pad_id)
    n_wells_at_multi_well_pads <- dplyr::filter(well_pad_well_counts, n > 1) %>%
      dplyr::pull("n") %>%
      sum_() %>%
      fill_na(0)
    n_wells_at_single_well_pads <- dplyr::filter(well_pad_well_counts, n == 1) %>%
      dplyr::pull("n") %>%
      sum_() %>%
      fill_na(0)
    glue_message(
      "Of {n_wells} total wells, {n_wells_at_multi_well_pads} were matched to ",
      "multi-well pads, {n_wells_at_single_well_pads} were matched to single-",
      "well pads, and {n_unmatched} were unmatched for data quality reasons."
    )
  }
  out
}


parse_year_range <- function(year_range) {
  # Use a regex header_df to pull out the years
  # look for things like "1990-2018"
  stopifnot(length(year_range) == 1)
  years <- stringr::str_match(year_range, "^(\\d{4})-(\\d{4})$")[1, 2:3] %>%
    as.integer() %>%
    unclass()
  stopifnot(length(years) == 2, !anyNA(years))
  years
}


create_well_pad_crosswalk <- function(header_dir, output_file, year_range) {
  year_range %<>% parse_year_range()
  stopifnot(dir.exists(header_dir), length(output_file) == 1)
  header_df <- arrow::open_dataset(header_dir) %>%
    dplyr::filter(
      .data$first_prod_year >= year_range[1],
      .data$first_prod_year <= year_range[2],
      !is.na(aapg_geologic_province),
      aapg_geologic_province != "(N/A)",
      !is.na(surface_latitude_wgs84),
      !is.na(surface_longitude_wgs84),
    ) %>%
    dplyr::select(
      entity_id, aapg_geologic_province,
      surface_latitude_wgs84, surface_longitude_wgs84,
      # NA operator_company_name is allowed (currently coded as "(N/A)")
      operator_company_name, entity_type
    ) %>%
    dplyr::collect() %>%
    dplyr::distinct()  # distinct for multiple wells at the same entity and location

  well_pad_crosswalk <- assign_well_pads(header_df, verbose=FALSE) %>%
    add_distance_to_nearest_well_pad() %>%
    add_count_wells_nearby()
  arrow::write_parquet(well_pad_crosswalk, output_file)
}


add_distance_to_nearest_well_pad <- function(well_df) {
  df <- well_df %>%
    dplyr::select(well_pad_id, well_pad_lon, well_pad_lat) %>%
    dplyr::distinct() %>%
    lonlat_to_projected(c("well_pad_lon", "well_pad_lat"))

  stopifnot(
    anyDuplicated(df$well_pad_id) == 0L,
    nrow(df) > 0
  )
  if (nrow(df) == 1) {
    dists <- NA_real_
  } else {
    suppressMessages(
      neighbor_info <- nngeo::st_nn(df, df, k=2L, returnDist=TRUE, progress=FALSE)
    )
    # Since we're matching df with df, the nearest neighbor is always the same
    # well pad, so we ask for two neighbors to get the next nearest well pad.
    # Therefore, we take the max of these two distances (the distance to self is
    # always zero)
    dists <- purrr::map_dbl(neighbor_info$dist, max)
    stopifnot(
      all(purrr::map_dbl(neighbor_info$dist, min) == 0),
      all(dists >= 10) # centroids should be s
    )
  }

  df_dists <- df %>%
    sf::st_drop_geometry() %>%
    dplyr::select(well_pad_id) %>%
    # More precisely, this is the distance in meters between one well pad
    # centroid to the nearest well pad's centroid.
    dplyr::mutate(distance_to_nearest_pad_m = !!dists)
  out <- dplyr::inner_join(well_df, df_dists, by="well_pad_id")
  out
}


add_count_wells_nearby <- function(well_df) {
  # Work from most populous basin first, just so the parallel processing doesn't
  # wait too long on single straglers.
  basins <- well_df %>%
    dplyr::count(aapg_geologic_province) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    purrr::chuck("aapg_geologic_province")

  counts_df <- furrr::future_map_dfr(
    basins, count_wells_nearby_one_basin,
    well_df=well_df,
    .options=furrr::furrr_options(seed=TRUE)
  )
  stopifnot(anyDuplicated(counts_df$well_pad_id) == 0)
  out <- dplyr::inner_join(well_df, counts_df, by="well_pad_id")
  out
}


count_wells_nearby_one_basin <- function(aapg_geologic_province, well_df) {
  # Maybe useful for speed improvements:
  # SAN JOAQUIN BASIN took 2.95 mins
  # PERMIAN BASIN took 25.64 secs
  # DENVER JULESBURG took 23.74 secs
  # SAN JUAN took 8.05 secs
  # PICEANCE took 2.56 secs
  # LOS ANGELES BASIN took 2.31 secs
  # RATON took 0.89 secs
  # SACRAMENTO BASIN took 0.57 secs
  # VENTURA BASIN took 0.43 secs
  # SANTA MARIA BASIN took 0.69 secs
  # SALINAS BASIN took 0.28 secs
  # LAS ANIMAS ARCH took 0.21 secs
  # GREEN RIVER took 0.17 secs
  # ANADARKO BASIN took 0.16 secs
  # PARADOX took 0.15 secs
  # BRAVO DOME took 0.17 secs
  # NORTH PARK took 0.14 secs
  # EEL RIVER BASIN took 0.12 secs
  # HALF MOON BASIN took 0.13 secs
  # IMPERIAL VALLEY BASIN took 0.12 secs
  # NORTHERN COAST PRVC took 0.12 secs
  # SONOMA BASIN took 0.12 secs
  # ESTANCIA BASIN took 0.12 secs
  # TUCUMCARI BASIN took 0.12 secs

  individual_well_geom <- well_df %>%
    dplyr::filter(.data$aapg_geologic_province == !!aapg_geologic_province) %>%
    dplyr::select(entity_id, surface_longitude_wgs84, surface_latitude_wgs84) %>%
    lonlat_to_projected(c("surface_longitude_wgs84", "surface_latitude_wgs84"))

  well_pad_df <- well_df %>%
    dplyr::filter(.data$aapg_geologic_province == !!aapg_geologic_province) %>%
    dplyr::select(well_pad_id, well_pad_lon, well_pad_lat) %>%
    dplyr::distinct()

  nearby_dist_km <- 10

  stopifnot(
    length(aapg_geologic_province) == 1L,
    anyDuplicated(well_df[["entity_id"]]) == 0L,
    anyDuplicated(well_df[["well_pad_id"]]) != 0L,
    nrow(individual_well_geom) > 0,
    nrow(well_pad_df) > 0,
    nearby_dist_km > 0
  )

  well_pad_buffered <- well_pad_df %>%
    lonlat_to_projected(c("well_pad_lon", "well_pad_lat")) %>%
    sf::st_buffer(dist = nearby_dist_km * 1000) # default to 10km
  intersect_count <- sf::st_intersects(well_pad_buffered, individual_well_geom) %>%
    purrr::map_int(length)

  # Note: this count includes wells in the well pad.
  out <- well_pad_df %>%
    dplyr::select(well_pad_id) %>%
    dplyr::mutate(tot_count_wells_within_10km = !!intersect_count)

  out
}


test_well_pad_creation <- function() {
  header_df <- tibble::tibble(
    entity_id = 1:4,
    aapg_geologic_province = "SAN JOAQUIN BASIN",
    entity_type = "WELL",
    operator_company_name="AAAAA",
    # First two points with identical location, third point about 9m away,
    # and fourth point 4.3km away
    surface_longitude_wgs84 = c(-119.7514, -119.7514, -119.7515, -119.8),
    surface_latitude_wgs84 = 35.48126,
  )

  expected_results <- header_df %>%
    dplyr::mutate(well_pad_id = c(1, 1, 1, 2))
  actual_results <- assign_well_pads(header_df, verbose=FALSE) %>%
    dplyr::select(-well_pad_lon, -well_pad_lat) # don't test that part
  compare_results <- all.equal(actual_results, expected_results)
  if (!isTRUE(compare_results)) {
    print(tibble::tibble(
      expected_id = expected_results$well_pad_id,
      actual_id = actual_results$well_pad_id,
    ))
    print(actual_results)
    stop(compare_results)
  }
}


if (!exists("snakemake")) {
  message("This script is meant to be run with snakemake. Using a placeholder.")
  snakemake <- SnakemakePlaceholder(
    input = list(
      headers = "data/generated/production/well_headers",
      script = "code/group_into_well_pads.R"
    ),
    output = list(
      "data/generated/production/well_pad_crosswalk_1970-2018.parquet"
    ),
    wildcards = list(year_range = "1970-2018"),
    threads = 4L,
    resources = list(mem_mb = 10000L)
  )
}

# Set up resources (mem limit doesn't work on MacOS)
memory_limit(snakemake@resources[['mem_mb']])
future::plan("multicore", workers = snakemake@threads)

test_well_pad_creation()

create_well_pad_crosswalk(
  header_dir = snakemake@input$headers,
  output_file = snakemake@output[[1]],
  year_range = snakemake@wildcards$year_range
)
