
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
    cTransfers = readRDS(file.path(projectPath, "tests/testthat/testdata", "cTransfers.rds")),

    # 2 pixels with a burn probability of 50%
    # 1 SW, 1 HW
    burnProbability = terra::rast(matrix(c(0.5, 0.5), ncol = 2, nrow = 1), crs = "local"),
    rasterToMatch = burnProbability,

    # Set all pool values to 1
    cbm_vars = {

      pools <- unique(c(cTransfers$source_pool, cTransfers$sink_pool))

      list(
        key   = data.table::data.table(pixelIndex = 1:2, row_idx = 1:2),
        state = data.table::data.table(row_idx = 1:2, spatial_unit_id = 27, sw_hw = c(0, 1)),
        pools = data.table::data.table(row_idx = 1:2)[, c(pools) := 1]
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


