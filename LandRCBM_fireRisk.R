
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
    expectsInput("burnProbability", "SpatRaster", sourceURL = NA,
                 desc = paste("Per-pixel burn probability for this year, from either",
                              "`fireSense_burnProbability` or `scfm_burnProbability`.")),
    expectsInput("cbm_vars", "list", desc = paste("List of 5 data tables defining active cohorts in the current year:",
                                                  "key, parameters, pools, flux, and state.",
                                                  "This is created initially during the spinup and updated each year.")),
    expectsInput("cTransfers", "data.table", desc = "Carbon transfer values table with associated disturbance names and IDs."),
    expectsInput("rasterToMatch", "SpatRaster", desc = "Raster template used to rasterize the fire susceptibility output.")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "fireSusceptibility", objectClass = "SpatRaster", desc = "Raster of carbon emitted to the atmosphere if burned."),
    createsOutput(objectName = "fireRisk", objectClass = "SpatRaster", desc = "Raster of fire risk to carbon.")
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

      terra::compareGeom(sim$rasterToMatch, sim$burnProbability, stopOnError = TRUE)

      sim$fireSusceptibility <- calculateFireSusceptibility(sim$cbm_vars, sim$cTransfers, sim$rasterToMatch)

      sim$fireRisk <- sim$burnProbability * sim$fireSusceptibility

      # schedule future event(s)
      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, "LandRCBM_fireRisk", "calculateFireRisk", eventPriority = 8.6)

    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

calculateFireSusceptibility <- function(cbm_vars, cTransfers, rasterToMatch){
  # keep only rows needed in cTransfers: spatial units and transfer to the atmosphere
  spuids <- unique(cbm_vars$state$spatial_unit_id)
  transfers <- cTransfers[cTransfers$disturbance_type_id == 1 & cTransfers$spatial_unit_id %in% spuids, ]
  transfers <- transfers[transfers$sink_pool %in% c("CO2", "CO", "CH4"),]

  # apply each transfer to the cohort groups
  cohortGroupConsequences <- rep(0L, nrow(cbm_vars$pools))

  for (i in seq_len(nrow(transfers))){
    t <- transfers[i,]

    cohortGroups <- (cbm_vars$state$sw_hw == (t$sw_hw == "hw")) & cbm_vars$state$spatial_unit_id == t$spatial_unit_id

    cohortGroupConsequences <- cohortGroupConsequences + (cbm_vars$pools[, get(t$source_pool)] * t$proportion * cohortGroups)
  }

  # convert the cohort group consequences to pixels
  pixelConsequences <- data.table(pixelIndex = cbm_vars$key$pixelIndex, consequence = cohortGroupConsequences[cbm_vars$key$row_idx])
  pixelConsequences <- pixelConsequences[ , .(consequence = sum(consequence)), by = "pixelIndex"]

  # rasterize
  fireSusceptibility <- rast(rasterToMatch, names = "fireSusceptibility", vals = NA)
  fireSusceptibility[pixelConsequences$pixelIndex] <- pixelConsequences$consequence

  return(fireSusceptibility)
}

.inputObjects <- function(sim) {

  return(invisible(sim))
}

