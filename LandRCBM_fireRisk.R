
defineModule(sim, list(
  name = "LandRCBM_fireRisk",
  description = "",
  keywords = "",
  authors = c(
    person("Dominique", "Caron",  email = "dominique.caron@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Susan",     "Murray", email = "murray.e.susan@gmail.com",          role = c("ctb"))
  ),
  childModules = character(0),
  version = list(LandRCBM_fireRisk = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "LandRCBM_fireRisk.Rmd"),
  reqdPkgs = list("PredictiveEcology/SpaDES.core@development (>= 3.1.2.9005)", "terra", "data.table"),
  loadOrder = list(after = c("fireSense_burnProbability", "scfm_burnProbability")),
  parameters = bindrows(
    defineParameter(".runInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first fire-risk estimate should occur."),
    defineParameter(".runInterval", "numeric", default = 1, NA, NA,
                    desc = paste("Years between fire-risk estimates. `NA` estimates once, at",
                                 "`.runInitialTime`. Each estimate uses that year's ignition/spread",
                                 "covariates and model predictions, held fixed across all `nIterations` draws.")),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("burnProbability", objectClass = "SpatRaster", desc = "Pixel burn probability for the current year"),
    expectsInput("cohortDT",        objectClass = "data.table", desc = "Table of cohort attributes and carbon pools"),
    expectsInput("standDT",         objectClass = "data.table", desc = "Table of pixel attributes"),
    expectsInput("cbm_defaults_db", objectClass = "character",  desc = "Path to CBM defaults SQLite database")
  ),
  outputObjects = bindrows(
    createsOutput("fireSusceptibility", objectClass = "SpatRaster",
                  desc = "Raster of carbon emitted to the atmosphere if burned (t/ha)"),
    createsOutput("fireRisk", objectClass = "SpatRaster",
                  desc = "Raster of fire risk to carbon")
  )
))

doEvent.LandRCBM_fireRisk = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {

      # eventPriority must run after scfm_burnProbability's "estimate" event (priority 7.51)
      # so burnProbability is available, AND after CBM_core's "annual_carbonDynamics" event
      # (priority 8.5): LandRCBM_split3pools's "annualIncrements" event (priority 7, which runs
      # earlier in the same timestep) strips spatial_unit_id, mean_annual_temperature, species,
      # area, and sw_hw from cbm_vars$state when it rebuilds cohort groups for the new year;
      # those columns are only restored once CBM_core runs during annual_carbonDynamics.
      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, "LandRCBM_fireRisk", "calculateFireRisk", eventPriority = 8.6)

    },
    calculateFireRisk = {

      # calculate fire susceptibility
      pixelEmissions <- calculateFireSusceptibility(sim$standDT, sim$cohortDT, sim$cbm_defaults_db)

      # rasterize fire susceptibility
      sim$fireSusceptibility <- terra::rast(sim$burnProbability, names = "fireSusceptibility")
      terra::set.values(sim$fireSusceptibility, cells = pixelEmissions$pixelIndex, values = pixelEmissions$emissions)

      # calculate fire risk
      sim$fireRisk <- sim$burnProbability * sim$fireSusceptibility

      # schedule future event(s)
      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, "LandRCBM_fireRisk", "calculateFireRisk", eventPriority = 8.6)

    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

# Calculate pixel emissions if a fire were to occur
calculateFireSusceptibility <- function(standDT, cohortDT, cbm_defaults_db){

  emissionsByPool <- fireEmissionProportions(cbm_defaults_db)

  pools <- intersect(names(emissionsByPool)[-1], gsub("^pools\\.", "", names(cohortDT)))

  cohortPools <- cohortDT[, .SD, .SDcols = c("pixelIndex", paste0("pools.", pools))]
  cohortPools[standDT, spatial_unit_id := spatial_unit_id, on = "pixelIndex"]

  cohortPools[, emissions := rowSums(
    cohortPools[, .SD, .SDcols = paste0("pools.", pools)] *
      emissionsByPool[match(cohortPools$spatial_unit_id, emissionsByPool$spatial_unit_id), .SD, .SDcols = pools]
  )]

  cohortPools[ , .(emissions = sum(emissions)), by = "pixelIndex"]
}

# Get a summary of pool emission proportions for each spatial unit in the event of a wildfire
fireEmissionProportions <- function(cbm_defaults_db){

  cbmDBcon <- RSQLite::dbConnect(RSQLite::dbDriver("SQLite"), cbm_defaults_db, flags = RSQLite::SQLITE_RO)
  on.exit(RSQLite::dbDisconnect(cbmDBcon))
  cbmDBtables <- c("disturbance_type_tr", "disturbance_matrix_association", "disturbance_matrix_value", "pool")
  cbmDBtables <- setNames(cbmDBtables, cbmDBtables) |>
    lapply(RSQLite::dbReadTable, con = cbmDBcon) |>
    lapply(data.table::as.data.table)

  distTypeID <- cbmDBtables$disturbance_type_tr[name == "Wildfire", disturbance_type_id]
  distValues <- merge(
    cbmDBtables$disturbance_matrix_association[disturbance_type_id == distTypeID],
    cbmDBtables$disturbance_matrix_value,
    by = "disturbance_matrix_id")

  fireEmissions <- distValues[
    sink_pool_id %in% cbmDBtables$pool[code %in% c("CO2", "CO", "CH4"), id],
    .(proportion = sum(proportion)),
    by = c("spatial_unit_id", "source_pool_id")]

  fireEmissions[cbmDBtables$pool[, .(source_pool_id = id, code)], source_pool := code, on = "source_pool_id"]
  fireEmissions <- data.table::dcast(fireEmissions, spatial_unit_id ~ source_pool, value.var = "proportion")
  fireEmissions[is.na(fireEmissions)] <- 0

  return(fireEmissions)

}

.inputObjects <- function(sim) {

  return(invisible(sim))
}

