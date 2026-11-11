## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE)

## ----flowchart, echo=FALSE, out.width="100%", fig.alt="Flowchart of a hellometry() call, from the input data through data_checker(), make_measurement_table(), full_estimation_table() for sizes then models, and predict_biomass(), to the returned list of three tables."----
knitr::include_graphics("hellometry_flowchart.svg")

## -----------------------------------------------------------------------------
# Load library
library(hellometry)

# We need to define a suite of levels to be used for estimation, i.e. columns in the data, 
# from lowest to coarsest resolution
level_vec  <-
  c("species", "genus", "family", "order")

# Read in the reference measurements, the real body sizes and body masses used to build the models
measurements <-
  bromeliad_inverts_measurements() %>%
  ## Rename columns to those expected by hellometry()
  dplyr::rename(size_col = body_size_mm,
                biomass_col = body_mass_mg,
                biomass_type = mass_type) %>%
  ## size_col holds both numbers and categories, so it must be character;
  ## biomass_col must be numeric for the allometric models
  dplyr::mutate(size_col = as.character(size_col),
                biomass_col = as.numeric(biomass_col))

# Have a look at the reference measurements
dplyr::glimpse(measurements)

# Read in the Trinidadian communities, the taxa we want size and biomass estimates for
communities <-
  trini_communities() %>%
  ## Rename the abundance column and add the columns hellometry() needs
  dplyr::rename(abundance = n) %>%
  ## Here we do not have any measurement for these invertebrates, so size is "unknown"
  ## and biomass NA, to be estimated. We use dry biomass because it tends to be more precise
  dplyr::mutate(size_col = "unknown",
                biomass_col = NA,
                biomass_type = "dry")

# Have a look at the communities to estimate
dplyr::glimpse(communities)

# Now we can use the package to estimate sizes and biomasses!

# Get the estimates
my_invertebrates_estimated <-
  hellometry(dats = communities, ## The data to be estimated
              level_vec  = level_vec , ## The taxonomic levels to be used
              biomass_type = "dry", ## Type of biomass to be estimated
              measurements = measurements) ## The reference measurements

# The result is a list of three elements
# - `data` - your data with estimated body sizes and biomasses, and column with
# information on the level at which were performed the estimation
dplyr::glimpse(my_invertebrates_estimated$data)
# - `size_estimates` - a tibble with all unique size estimates that were joined to your data
dplyr::glimpse(my_invertebrates_estimated$size_estimates)
# - `model_estimates` - a tibble with all unique allometric models that were joined to your data
dplyr::glimpse(my_invertebrates_estimated$model_estimates)


# But now let's say you do not really want to do estimation on the data, 
# but just see the results of all possible estimations. Fear not, 
# it is easy!
## First combine communities and reference measurements
my_invertebrates <-
  combine_measurements(dats = communities,
                       measurements = measurements)
## Then compile the measurement table
measurement_table <-
  make_measurement_table(dats = my_invertebrates, ## Data to be used
                         level_vec  = level_vec) ## Taxonomic levels to be used
## Get all possible size estimates
size_estimation_table <-
  full_estimation_table(level_vec  = level_vec ,
                        measurement_table = measurement_table,
                        what = "size_col")
### Have a glimpse of it
dplyr::glimpse(size_estimation_table)

## Get all possible models
model_estimation_table <-
  full_estimation_table(level_vec  = level_vec ,
                        measurement_table = dry_wet(measurement_table,
                                                        biomass_type = "dry"), ## Filter for dry biomass
                        what = "biomass_col")
### Have a glimpse of it
dplyr::glimpse(model_estimation_table)


## -----------------------------------------------------------------------------

# Read in data with traits
trait_data <-
  trini_communities(traits = TRUE)

# The fuzzy traits are every column after the taxonomy, i.e. from the first trait
# column (AS1) up to the last column of the table
trait_columns <-
  trait_data %>%
  ## Get the block of trait columns ..
  dplyr::select(AS1:dplyr::last_col()) %>%
  ## Keep just their names
  names()

# Rebuild a measurement table, this time asking it to keep every taxonomic level
# we want to match on (make_measurement_table() only keeps the levels we give it)
measurement_table <-
  make_measurement_table(dats = my_invertebrates,
                         level_vec = c("species", "genus", "subfamily",
                                       "family", "order", "class"))

# Here we want to match traits from different morphospecies, look for matches 
# from genus up to class
match_levels <-
  c("genus", "subfamily", "family", "order", "class")

# Build a lookup holding one representative trait profile per taxon, at 
# every level combined into a single long table that remembers the level 
# each profile came from
trait_lookup <-
  match_levels %>%
  purrr::map_dfr(.,
                 ~ trait_data %>%
                     ## Drop the taxa that are unnamed at this level
                     dplyr::filter(!is.na(.data[[.x]]), .data[[.x]] != "") %>%
                     ## One row per taxon name, carrying its trait profile
                     dplyr::distinct(key = .data[[.x]],
                                     dplyr::across(dplyr::all_of(trait_columns))) %>%
                     ## Keep a single profile should a taxon appear more than once
                     dplyr::distinct(key, .keep_all = TRUE) %>%
                     ## Tag which taxonomic level this profile belongs to
                     dplyr::mutate(match_level = .x))

# Now hand each measured taxon the trait profile of its most similar relative, i.e.
# the profile found at the finest taxonomic level at which the taxon matches
trait_assignment <-
  measurement_table %>%
  ## Tag each measurement row so we can return its assigned traits later
  dplyr::mutate(.row = dplyr::row_number()) %>%
  ## Keep only that tag and the taxonomy columns we match on
  dplyr::select(.row, dplyr::all_of(match_levels)) %>%
  ## Stretch to one row per measured taxon per taxonomic level
  tidyr::pivot_longer(dplyr::all_of(match_levels),
                      names_to = "match_level", values_to = "key") %>%
  ## Drop the levels where the taxon has no name
  dplyr::filter(!is.na(key), key != "") %>%
  ## Attach every trait profile that matches, at whatever level
  dplyr::inner_join(trait_lookup, 
                    by = c("match_level", "key")) %>%
  ## Rank the levels from finest to coarsest, following match_levels
  dplyr::mutate(.priority = match(match_level, 
                                  match_levels)) %>%
  ## For each measurement row, keep only the finest matching level
  dplyr::group_by(.row) %>%
  dplyr::slice_min(.priority, 
                   n = 1, 
                   with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  ## Return just the row tag and its assigned trait profile
  dplyr::select(.row, dplyr::all_of(trait_columns))

# Finally, join the assigned traits back onto the measurements, and keep only taxa
# that both received traits and are identified to family (our id_col below)
dats <-
  measurement_table %>%
  ## Recreate the same row tag to join on
  dplyr::mutate(.row = dplyr::row_number()) %>%
  ## Bring in each taxon's assigned trait profile
  dplyr::left_join(trait_assignment, by = ".row") %>%
  ## Drop the helper tag
  dplyr::select(-.row) %>%
  ## Keep taxa that got traits and have a family to identify them by
  dplyr::filter(!is.na(.data[[trait_columns[1]]]),
                !is.na(family), family != "")

## Get long format data of which taxa have similar traits
matched_traits <-
  matcher_of_traits(measurement_table = dats,
                    trait_columns = trait_columns,
                    id_col = "family") ## here taxa are identified by family
### See the output
head(matched_traits)

## Get clean list for full_estimation_table, includes matcher_of_traits()
## Here you can see which taxa were grouped together
clean_trait_list <-
  make_trait_table(measurement_table = dats,
                   trait_columns = trait_columns,
                   id_col = "family")
### See the output
clean_trait_list[[1]]

# The last two functions are called in full_estimation_table if you add the trait arguments
## Size estimation table
size_estimation_table_traits <-
  full_estimation_table(level_vec  = c(),
                        measurement_table = dats,
                        traits = TRUE, ## switch for traits
                        trait_columns = trait_columns, ## vector of fuzzy trait columns
                        id_col = "family", ## where to group
                        what = "size_col")
### See the output
head(size_estimation_table_traits)

## Model table
model_estimation_table_traits <-
  full_estimation_table(level_vec  = c(),
                        measurement_table = dry_wet(dats,
                                                    biomass_type = "dry"), ## Filter for dry biomass
                        traits = TRUE, ## traits = TRUE
                        trait_columns = trait_columns, ## vector of fuzzy trait columns
                        id_col = "family", ## column identifying taxa
                        what = "biomass_col")
### See the output
head(model_estimation_table_traits)

