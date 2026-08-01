files <- c(
  "app.R", "README.md", "LICENSE",
  list.files("R", recursive = TRUE, full.names = TRUE),
  list.files("data/simulated", recursive = TRUE, full.names = TRUE),
  list.files("docs", recursive = TRUE, full.names = TRUE)
)

rsconnect::deployApp(
  appDir = ".",
  appName = "HEKBlueR",
  account = "caffeinated-code",
  server = "shinyapps.io",
  appFiles = files,
  forceUpdate = TRUE
)

