

phylo_lm_util <- function(formula, data, tree, model, method, n.mvmorph.cores, boot = 0, ...) {
  dots <- list(...)
  dots_glm <- dots[names(dots) %in% names(formals(phylolm::phyloglm))]
  dots_lm <- dots[names(dots) %in% names(formals(phylolm::phylolm))]

  if (length(union(names(dots_glm), names(dots_lm))) != length(dots)) {
    warning("Some arguments in ... are not recognized.", call. = FALSE)
  }
  # we capture the first argument in the formula, to check whether it is binary
  x_var <- data[[all.vars(formula)[1]]]
  #NOW CHECK IF THE x_var (in this case, the response variable) is multivariate
  if (is.factor(x_var)) {
    # phyloglm need binary variables as 0,1 but I use factors
    data[all.vars(formula)[1]] <- as.numeric(x_var) - 1
    fun <- phylolm::phyloglm
    args <- c(list(formula = formula, data = data, phy = tree, method = method, boot = boot),
              dots_glm)
  } else {
    if (dim(x_var)[2] == 1){
      fun <- phylolm::phylolm
      args <- c(list(formula = formula, data = data, phy = tree, model = model, boot = boot),
                dots_lm)
    }
    else {
      res.temp <- mvMORPH::mvgls(formula = formula, data = data, tree = tree, model = model, method=c("PL-LOOCV") ,REML=TRUE)
      fun <- mvMORPH::manova.gls
      #TODO: add a check that n.mvmorph.cores is a positive integer
      args <- c(list(object = res.temp, test="Pillai", nbcores=n.mvmorph.cores, type = "I", nperm=1000))

    }
  }
  res <- do.call(phylopath:::quiet_safely(fun), args)
  if(exists('res.temp')){
    res$result$optpar <- res.temp$opt$par[2]
  }
  # Remove the call, since quiet_safely messes it up and it's annoying in printing
  res$result$call <- NULL

  return(res)
}


get_p2 <- function(m) {
  p <- switch(
    class(m),
    phylolm = stats::coef(phylolm::summary.phylolm(m))[, 'p.value'] %>% dplyr::last(),
    phyloglm = stats::coef(phylolm::summary.phyloglm(m))[, 'p.value'] %>% dplyr::last(),
    manova.mvgls = m$pvalue %>% dplyr::last()
  )
  if (p < .Machine$double.eps) p <- .Machine$double.eps
  return(p)
}


hiDsummary <- function(object, ...) {
  phylopath <- object
  stopifnot(inherits(phylopath, 'phylopath'))
  n <- nrow(phylopath$data[[1]])
  k <- sapply(phylopath$d_sep, nrow)
  q <- sapply(phylopath$model_set, function(m) nrow(m) + sum(m))
  C <- sapply(phylopath$d_sep, function(x) phylopath:::C_stat(x$p))
  p <- phylopath:::C_p(C, k)
  IC <- phylopath:::CICc(C, q, n)

  if (n <= max(q) + 1) {
    IC[n <= q] <- NA
    warning('CICc was not calculated for causal models where the number of parameters is equal to',
            'or larger than the number of species.')
  }

  d <- data.frame(
    model = names(phylopath$model_set), k = k, q = q, C = C, p = p,
    CICc = IC, stringsAsFactors = FALSE
  )
  d <- d[order(d$CICc), ]
  d$delta_CICc <- d$CICc - d$CICc[1]
  d$l <- phylopath:::l(d$delta_CICc)
  d$w <- phylopath:::w(d$l)
  class(d) <- c('phylopath_summary', 'data.frame')
  return(d)
}

HiDPPA <- function(model_set, data, tree, model = 'lambda', method = 'logistic_MPLE', n.mvmorph.cores = 1,
                   order = NULL, parallel = NULL, na.rm = TRUE, ...) {

  ###TODO: add a check for these model types
  # tmp <- check_models_data_tree(model_set, data, tree, na.rm)
  #model_set <- tmp$model_set
  #data <- tmp$data
  #tree <- tmp$tree

  if (is.null(order)) {
    order <- phylopath:::find_consensus_order(model_set)
  }
  formulas <- lapply(model_set, phylopath:::find_formulas, order)
  formulas <- purrr::map(
    formulas,
    ~purrr::map(.x, ~{attr(., ".Environment") <- NULL; .})
  )
  f_list <- unique(unlist(formulas))
  if (!is.null(parallel)) {
    cl <- parallel::makeCluster(
      min(c(parallel::detectCores() - 1, length(f_list))), parallel
    )
    parallel::clusterExport(cl, list('phylo_lm_util'), environment())
    on.exit(parallel::stopCluster(cl))
  } else {
    cl <- NULL
  }
  dsep_models_runs <- pbapply::pblapply(
    f_list,
    function(x, data, tree, model, method, ...) {
      phylo_lm_util(x, data, tree, model, method, n.mvmorph.cores, ...)
    },
    data = data, tree = tree, model = model, method = method, cl = cl)
  # Produce appropriate error if needed
  errors <- purrr::map(dsep_models_runs, 'error')
  purrr::map2(
    errors, f_list,
    ~if(!is.null(.x))
      stop(paste(
        'Fitting the following model:\n   ',
        Reduce(paste, deparse(.y)),
        '\nproduced this error:\n   ',
        .x
      ), call. = FALSE)
  )
  # Collect warnings as well, but save those for later.
  warnings <- purrr::map(dsep_models_runs, 'warning')
  warnings <- purrr::map2(
    warnings, f_list,
    ~if(!is.null(.x))
      paste(
        'Fitting the following model:\n   ',
        Reduce(paste, deparse(.y)),
        '\nproduced this/these warning(s):\n   ',
        .x
      )
  )
  warnings <- warnings[!sapply(warnings, is.null)]
  if (length(warnings) > 1) {
    warning('Some models produced warnings. Use `show_warnings()` to view them.')
  }

  # Collect models.
  dsep_models <- purrr::map(dsep_models_runs, 'result')
  dsep_models <- purrr::map(formulas, ~dsep_models[match(.x, f_list)])

  d_sep <- purrr::map2(
    formulas,
    dsep_models,
    ~tibble::tibble(
      d_sep = as.character(.x),
      p = purrr::map_dbl(.y, get_p2),
      phylo_par = purrr::map_dbl(.y, phylopath:::get_phylo_param),
      model = .y
    )
  )


  out <- list(
    d_sep = d_sep, model_set = model_set, data = data, tree = tree, model = model, method = method,
    dots = list(...), warnings = warnings
  )
  class(out) <- 'phylopath'
  return(out)
}
