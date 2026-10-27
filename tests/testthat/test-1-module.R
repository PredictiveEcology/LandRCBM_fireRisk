
test_that("Module: LandRCBM_fireRisk", {

  # Set up project
  projectPath <- ifelse(basename(getwd()) == "testthat", dirname(dirname(getwd())), getwd())

  simInitInput <- SpaDES.project::setupProject(

    modules = c(
      "LandRCBM_fireRisk"
    ),
    times = list(start = 2000, end = 2000),
    paths = list(
      projectPath = projectPath,
      modulePath  = dirname(projectPath),
      packagePath = Sys.getenv("R_LIBS_USER"),
      inputPath   = file.path(tempdir(), "test_LandRCBM_fireRisk", "inputs"),
      cachePath   = file.path(tempdir(), "test_LandRCBM_fireRisk", "cache"),
      outputPath  = file.path(tempdir(), "test_LandRCBM_fireRisk", "outputs")
    ),

    require = c("data.table", "terra"),

    # Carbon transfer proportions
    cbm_defaults_db = {
      dbURL <- "https://raw.githubusercontent.com/cat-cfs/libcbm_py/main/libcbm/resources/cbm_defaults_db/cbm_defaults_v1.2.9300.391.db"
      destfile <- file.path(paths$inputPath, basename(dbURL))
      if (!file.exists(destfile)) download.file(dbURL, destfile, mode = "wb", quiet = TRUE)
      destfile
    },

    # 2 pixels with a burn probability of 50%
    # 1 SW, 1 HW
    # Set all pool values to 1
    burnProbability = terra::rast(matrix(c(0.5, 0.5), ncol = 2, nrow = 1), crs = "local"),

    standDT = data.table::data.table(pixelIndex = 1:2, spatial_unit_id = 27),
    cohortDT = {
      pools <- list(
        sw    = c("SoftwoodMerch", "SoftwoodFoliage", "SoftwoodOther",
                  "SoftwoodCoarseRoots", "SoftwoodFineRoots", "SoftwoodStemSnag", "SoftwoodBranchSnag"),
        hw    = c("HardwoodMerch", "HardwoodFoliage", "HardwoodOther",
                  "HardwoodCoarseRoots", "HardwoodFineRoots", "HardwoodStemSnag", "HardwoodBranchSnag"),
        other = c("AboveGroundVeryFastSoil", "BelowGroundVeryFastSoil", "AboveGroundFastSoil", "BelowGroundFastSoil",
                  "MediumSoil", "AboveGroundSlowSoil", "BelowGroundSlowSoil", "CO2", "CH4", "CO", "NO2", "Products")
      ) |> lapply(function(x) paste0("pools.", x))
      rbind(
        data.table::data.table(pixelIndex = 1)[, c(pools$other) := 1][, c(pools$sw) := 1][, c(pools$hw) := 0],
        data.table::data.table(pixelIndex = 2)[, c(pools$other) := 1][, c(pools$sw) := 0][, c(pools$hw) := 1]
      )
    }
  )

  # Run simInit
  simTestInit <- SpaDES.core::simInit2(simInitInput)
  expect_s4_class(simTestInit, "simList")

  # Run spades
  simTest <- SpaDES.core::spades(simTestInit)
  expect_s4_class(simTest, "simList")

  # Check outputs
  pixelEmissions <- c(
    3.775928, # SW
    4.005165  # HW
  )

  expect_true(terra::compareGeom(simTest$fireSusceptibility, simTest$burnProbability))
  expect_equal(terra::values(simTest$fireSusceptibility, mat = FALSE), pixelEmissions,
               tolerance = 0.000001, scale = 1)

  expect_true(terra::compareGeom(simTest$fireRisk, simTest$burnProbability))
  expect_equal(terra::values(simTest$fireRisk, mat = FALSE), pixelEmissions * 0.5,
               tolerance = 0.000001, scale = 1)

})


