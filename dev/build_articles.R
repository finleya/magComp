# Build the local website, article pages, and downloads.
# Run from the package root: Rscript dev/build_articles.R

Sys.setenv(TZ = "UTC")

unlink("docs/.nojekyll")
pkgdown::clean_site()
source("dev/article-data/build_downloads.R", chdir = TRUE)

pkgdown::init_site()
pkgdown::build_home()
pkgdown::build_reference(examples = FALSE)
pkgdown::build_news()
pkgdown::build_articles_index()

# pkgdown moves standalone Quarto HTML previews while building the site.
# Preserve those previews; these are not cached model fits or predictions.
preview_dir <- tempfile("magcomp-article-previews-")
dir.create(preview_dir)
previews <- list.files("vignettes/articles", pattern = "\\.html$", full.names = TRUE)
file.copy(previews, preview_dir)

# Build each nested Quarto article using pkgdown's public article API.
articles <- pkgdown::as_pkgdown(".")$vignettes$name

tryCatch({
  # Keep Quarto from copying linked standalone pages over pkgdown pages.
  unlink(previews)
  for (article in articles) {
    pkgdown::build_article(article, quiet = FALSE)
  }
}, finally = {
  file.copy(list.files(preview_dir, full.names = TRUE),
            "vignettes/articles", overwrite = TRUE)
  unlink(preview_dir, recursive = TRUE)
})

pkgdown::build_search()

# GitHub Pages serves these generated files without running Jekyll.
file.create("docs/.nojekyll")
