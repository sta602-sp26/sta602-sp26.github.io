deps <- function(path = ".") {
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
  }
  unique(renv::dependencies(path)$Package)
}

needed <- deps(".")
installed <- rownames(installed.packages())
missing <- setdiff(needed, installed)

if (length(missing)) {
  install.packages(missing)
}
# 
# quarto::quarto_render(".")
