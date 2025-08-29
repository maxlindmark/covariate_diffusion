# Modified from: https://github.com/pbs-assess/sdmTMB/blob/main/R/caic.R
# Zheng, N., Cadigan, N., & Thorson, J. T. (2024).
# A note on numerical evaluation of conditional Akaike information for
# nonlinear mixed-effects models (arXiv:2411.14185). arXiv.
# \doi{10.48550/arXiv.2411.14185}

#cAIC.sdmTMB <- function(object, what = c("cAIC", "EDF"), ...) {
cAIC.TMB <- function(obj, tmb_data, parlist, random, what = c("cAIC", "EDF")) {

  what <- tolower(what)
  what <- match.arg(what, choices = c("caic", "edf"))

  # if ("edf" %in% names(object) && what == "edf") {
  #   return(object$edf)
  # }

  # tmb_data <- object$tmb_data

  ## Ensure profile = NULL
  # if (is.null(object$control$profile)) {
  #   obj <- object$tmb_obj
  # } else {
  #   obj <- TMB::MakeADFun(
  #     data = tmb_data,
  #     parameters = object$parlist,
  #     map = object$tmb_map,
  #     random = object$tmb_random,
  #     DLL = "sdmTMB",
  #     profile = NULL #<
  #   )
  # }

  ## Make obj_new
  tmb_data$weights_i[] <- 0
  obj_new <- TMB::MakeADFun(
    data = tmb_data,
    parameters = parlist, #parameters = object$parlist,
    map = obj$env$map, #map = object$tmb_map,
    random = random, #random = object$tmb_random,
    DLL = obj$env$DLL, #DLL = "sdmTMB",
    profile = NULL
  )

  par <- obj$env$parList()
  parDataMode <- obj$env$last.par.best # FIXME: obj$env$last.par ??
  indx <- obj$env$lrandom()
  q <- sum(indx)
  p <- length(obj$par)

  ## use '-' for Hess because model returns negative loglikelihood
  #if (is.null(object$tmb_random)) {
  if (is.null(random) || length(random) == 0) {
    cli_inform(c("This model has no random effects.", "cAIC and EDF only apply to models with random effects."))
    return(invisible(NULL))
  }

  Hess_new <- -Matrix::Matrix(obj_new$env$f(parDataMode, order = 1, type = "ADGrad"), sparse = TRUE)
  Hess_new <- Hess_new[indx, indx] ## marginal precision matrix of REs

  ## Joint hessian etc
  Hess <- -Matrix::Matrix(obj$env$f(parDataMode, order = 1, type = "ADGrad"), sparse = TRUE)
  Hess <- Hess[indx, indx]
  negEDF <- Matrix::diag(Matrix::solve(Hess, Hess_new, sparse = FALSE))

  if (what == "caic") {
    jnll <- obj$env$f(parDataMode)
    cnll <- jnll - obj_new$env$f(parDataMode)
    cAIC_out <- 2 * cnll + 2 * (p + q) - 2 * sum(negEDF)
    return(cAIC_out)
  } else if (what == "edf") {
    # ## Figure out group for each random-effect coefficient
    # group <- names(object$last.par.best[obj$env$random])
    #
    # convert_bsmooth2names <- function(object, model = 1) {
    #   sn <- row.names(print_smooth_effects(object, m = model, silent = TRUE)$smooth_sds)
    #   sn <- gsub("^sd", "", sn)
    #   dms <- object$smoothers$sm_dims
    #   unlist(lapply(seq_along(dms), \(i) rep(sn[i], dms[i])))
    #
    # }
    # s_groups <- convert_bsmooth2names(object)
    # # smoothers always shared in delta models
    # if (is_delta(object)) s_groups <- c(paste0("1LP-", s_groups), paste0("2LP-", s_groups))
    # group[group == "b_smooth"] <- s_groups
    # group <- factor(group)
    #
    # ## Calculate total EDF by group
    # EDF <- tapply(negEDF, INDEX = group, FUN = length) - tapply(negEDF, INDEX = group, FUN = sum)
    # return(EDF)
  } else {
    cli_abort("Option not implemented")
  }
}
