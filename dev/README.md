# Development resources

This directory contains resources for maintaining the package and its
website. It is included in the source repository and excluded from R
package builds. The website article sources in `vignettes/articles/` are
also excluded from package builds.

## Package checks

From the package root, build and check the source package:

```sh
R CMD build .
R CMD check --no-manual magComp_0.0.6.tar.gz
```

Routine tests are in `tests/testthat/`. Additional numerical comparisons
and simulation checks are described in [validation/README.md](validation/README.md).

## Website articles and data

The [article-data directories](article-data/) contain inputs, preparation
scripts, and population generators. Each example has a README describing
its data and how to reproduce them. The shared
[plot-sampling dataset](article-data/plot-sampling/README.md) supplies the
forest geometry and canopy-height covariates for the series.

Install the package before rendering articles. Rendering uses Quarto and
the R packages magComp, tidyverse, sf, knitr, and rmarkdown. The complete
website also uses pkgdown and the quarto R package.

Render one article from the package root with:

```sh
quarto render vignettes/articles/magnitude-plot-sampling.qmd
```

To refresh the downloadable data and build the website, run:

```sh
Rscript dev/build_articles.R
```

The builder uses `_pkgdown.yml` for navigation and writes the site to
`docs/`. It preserves standalone HTML previews beside the article sources.
Articles execute their analyses when rendered; fitted models and predictions
are not cached. Geometry preparation and optional population generation are
separate steps described in each dataset README.

## Website publication

The `master` branch holds the package and article sources. The `gh-pages`
branch holds only the generated website. In the repository's GitHub Pages
settings, select **Deploy from a branch**, **gh-pages**, and **/(root)**.

After building and reviewing `docs/index.html`, publish the site from the
package root with:

```sh
bash dev/publish_site.sh
```

This requires Git, rsync, and push access to the repository. It updates the
existing `gh-pages` branch using a temporary checkout. Publishing does not
rerun the article analyses; rerun the website builder when their sources or
the package change.

## Technical documentation

The [technical document](modeling_and_software_details/magComp.tex) describes
the model, priors, sampler, numerical methods, prediction, and aggregation.
To rebuild its PDF, run from `dev/modeling_and_software_details/`:

```sh
latexmk -pdf magComp.tex
```

Rebuild the function reference manual from the package root with:

```sh
R CMD Rd2pdf --no-preview --force --output=dev/manual/magComp-manual.pdf .
```
