required_packages <- c(
  "here",
  "readxl",
  "knitr",
  "kableExtra"
)

installed <- rownames(installed.packages())
missing <- setdiff(required_packages, installed)

if (length(missing) > 0) {
  install.packages(missing, dependencies = TRUE)
}

invisible(lapply(required_packages, library, character.only = TRUE))