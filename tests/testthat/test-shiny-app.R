# Integration smoke test: drive the real Shiny app headlessly and confirm it
# launches and exposes its tabs. Guarded so it skips cleanly where the browser
# tooling (shinytest2 + chromote + Chrome) is unavailable -- e.g. local dev
# machines without Chrome -- and runs in CI where it is installed.

test_that("the app launches and exposes its core tabs", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  skip_if_not(nzchar(Sys.getenv("CHROMOTE_CHROME")) ||
                !is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
              "No Chrome/Chromium found for headless testing")

  app <- shinytest2::AppDriver$new(
    shiny::shinyApp(app_ui(), app_server),
    name = "rnaflow-launch",
    width = 1400, height = 900,
    load_timeout = 60 * 1000
  )
  on.exit(app$stop(), add = TRUE)

  html <- app$get_html("body")

  # Brand + stylesheet wiring loaded
  expect_match(html, "rnaflow-brand", fixed = TRUE)

  # All 13 tabs are present in the navbar
  for (tab in c("Data", "Volcano", "Explore", "Heatmap", "PCA", "QC",
                "Compare", "Enrichment", "Network", "Activity", "AI",
                "Project", "Report")) {
    expect_match(html, tab, fixed = TRUE)
  }

  # Fresh launch shows the getting-started guidance, not an error
  expect_match(html, "Getting started", fixed = TRUE)
})
