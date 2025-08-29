## Estimating scale-dependent covariate responses using two-dimensional diffusion derived from the SPDE method
In this paper we introduce a novel method for estimating spatially weighted covariates based on applying local diffusion to the covariate internally in the model. We conduct simulation experiments to verify the model's ability to retrieve diffused covariate effects, and then apply this model to two case studies; the breeding bird survey in the northwest US and scientific trawl data from the eastern Bering Sea.

### Reproducing Results

To reproduce our results you can either:

1. Fork the repository, clone it, open a new RStudio project with version control, and paste the repo url

2. Download a zip and work locally on your computer

We use [`renv`](https://rstudio.github.io/renv/articles/renv.html) to manage package versions. Once you've downloaded the project, run `renv::restore()` in your current working directory. This will install the package versions we used when this repository was archived. Note that packages are installed in a stand-alone project library for this paper, and will not affect your installed R packages anywhere else! `renv` does *not* help with different versions of R. We used R version 4.3.2, and ran the analysis on a 24 GB Apple M2 Sequoia 15.6.1 laptop.

If you are more interested in a simple example of the steps to fit a diffusion model than reproducing our results in the paper, we made a short [vignette](https://maxlindmark.github.io/covariate_diffusion/) demonstrating a typical workflow with a minimal example on capelin.

### Repository structure

`analysis`: code to prepare data, run case studies (note these scripts source the data preparation scripts), run simulation experiment and make figures. The ebs case study code and the simulation experiment takes a few hours to run, while the bb case study is faster.

`data`: data for case studies, including shapefiles, bathymetry rasters, and species name tables

`docs`: .qmd and .html of simple vignette

`functions`: code to produce base maps

`results`: output and figures from case studies and experiments

`timing`: code to test timing of fitting diffused vs standard model, and to calculate \(\mathbf{D}^{-1}\mathbf{x}\) using either sparse LU decomposition or constructing the dense matrix

`tmb`: .cpp files for all models (case studies, simulation)

