## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "LandRCBM_fireRisk",
  description = "",
  keywords = "",
  authors = structure(list(list(given = "Dominique", family = "Caron", role = c("aut", "cre"), email = "dominique.caron@nrcan-rncan.gc.ca")), class = "person"),
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
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
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
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = "fireSusceptibility", objectClass = "SpatRaster", desc = "Raster of carbon emitted to the atmosphere if burned."),
    createsOutput(objectName = "fireRisk", objectClass = "SpatRaster", desc = "Raster of fire risk to carbon.")
  )
))

doEvent.LandRCBM_fireRisk = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {

      # do stuff for this event
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
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

