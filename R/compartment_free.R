# Compartment-free models, source-name collisions and the variance floor.
#
# The cases translator#74 takes on, from Bill Denney's standalone nlmixr2
# translator for model-based meta-analysis (FeRx-NLME/ferx-r#721): a model()
# block with no linCmt() and no d/dt() whose error line reads an ordinary
# variable, `pred ~ add(sd)`. Every other statement already reaches
# [individual_parameters] through .parse_model_exprs(); what was missing is the
# prediction itself, which ferx takes as `y = <expr>` in [structural_model]
# (its "Form C" output expression, the same line an ODE readout uses).

# -- compartment-free lowering ------------------------------------------------

# ferx log-packs the Cholesky diagonal of a variance inside a box whose lower
# rail is exp(-12) = 6.144e-6 on the variance scale. Measured against ferx-r
# df7ffa7 (maxiter = 0, so the OFV is a function of the model text): an
# ESTIMATED omega written as 1e-8 is silently represented as 6.144e-6 and the
# fit reports that value, while the same omega with FIX is held at exactly
# 1e-8 and gives a different OFV (7932.16 against 7930.78). So the floor is for
# estimated variances only -- flooring a fixed one would change the model it was
# meant to preserve.
.VARIANCE_GUARD <- exp(-12)
# Where an estimated variance below the guard is started instead: clear of the
# rail, so the optimizer does not begin on the box edge. Bill Denney's value.
.VARIANCE_FLOOR <- 1e-5

# Move every estimated diagonal omega and kappa below ferx's guard up to
# .VARIANCE_FLOOR, and say so. The usual source is an nlmixr2 fit whose
# variance collapsed towards zero, handed to the translator as a start.
#
# A block omega's diagonal is reported but not moved: raising one diagonal
# element of a covariance matrix can make it indefinite, which is a worse
# model than the one the source wrote.
.floor_variances <- function(omegas, kappas) {
  warn <- character()
  msg <- function(kind, nm, v) paste0(
    "WARN  | ", kind, " ", nm, " starts at ", format(v, digits = 6),
    ", below the smallest variance ferx can represent while estimating it (",
    format(.VARIANCE_GUARD, digits = 4), ", exp(-12)) -- started at ",
    format(.VARIANCE_FLOOR), " instead. Left alone, ferx would silently move it to ",
    "its lower rail. A FIX variance is held exactly and is not moved.")
  for (i in seq_along(omegas)) {
    o <- omegas[[i]]
    if (isTRUE(o$fixed)) next
    if (identical(o$type, "block")) {
      k <- length(o$names)
      dpos <- cumsum(seq_len(k))           # diagonal positions in the lower triangle
      low <- which(o$values[dpos] < .VARIANCE_GUARD)
      for (j in low) warn <- c(warn, paste0(
        "WARN  | block omega element ", o$names[j], " starts at ",
        format(o$values[dpos[j]], digits = 6), ", below ferx's variance guard (",
        format(.VARIANCE_GUARD, digits = 4), "). It is NOT moved: raising one ",
        "diagonal of a covariance block can make it indefinite. Start it higher ",
        "in the source."))
      next
    }
    if (is.finite(o$values) && o$values < .VARIANCE_GUARD) {
      warn <- c(warn, msg("omega", o$names, o$values))
      omegas[[i]]$values <- .VARIANCE_FLOOR
    }
  }
  for (i in seq_along(kappas)) {
    k <- kappas[[i]]
    if (isTRUE(k$fixed)) next
    if (is.finite(k$value) && k$value < .VARIANCE_GUARD) {
      warn <- c(warn, msg("kappa", k$name, k$value))
      kappas[[i]]$value <- .VARIANCE_FLOOR
    }
  }
  list(omegas = omegas, kappas = kappas, warnings = warn)
}

# -- source-name collisions ---------------------------------------------------

# Every symbol a model() statement list assigns, every symbol it reads, and
# every ODE state it declares (`d/dt(x) <-`), at any depth -- including inside
# `if` arms, which .covariate_names() and .unmapped_symbols() do not descend
# into. A state is read on right-hand sides like any variable, so without its
# own set it would look like a free symbol, i.e. a data column.
.walk_names <- function(lst) {
  assigned <- character()
  used     <- character()
  states   <- character()
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    head <- as.character(e[[1]])[1L]
    if (identical(head, "{")) { for (x in as.list(e)[-1]) walk(x); return(invisible()) }
    if (identical(head, "if")) {
      used <<- c(used, .collect_symbols(e[[2]]))
      for (x in as.list(e)[-(1:2)]) walk(x)
      return(invisible())
    }
    if (.is_assignment(e)) {
      lhs <- e[[2]]
      if (is.symbol(lhs)) assigned <<- c(assigned, as.character(lhs))
      else if (.is_ddt_lhs(lhs)) states <<- c(states, .ddt_state(lhs))
      used <<- c(used, .collect_symbols(e[[3]]))
      return(invisible())
    }
    if (.is_tilde(e)) used <<- c(used, .collect_symbols(e))
  }
  for (e in lst) walk(e)
  list(assigned = unique(assigned), used = unique(used), states = unique(states))
}

# The data columns a model() statement list reads, as the source spells them
# (ferx binds a covariate to its column by exact header). Once a [covariates]
# block exists ferx expects every covariate in it and warns (W_PARSE) about the
# rest, so a block emitted for one column has to list them all.
.model_covariates <- function(lst, name_map) {
  nm <- .walk_names(lst)
  covs <- setdiff(nm$used, c(nm$assigned, nm$states, names(name_map)))
  covs[!toupper(covs) %in% .RESERVED_ODE_NAMES]
}

# Distinct source names that ferx would read as ONE name.
#
# The two kinds of name are compared differently, because ferx treats them
# differently. A model variable is emitted through .norm(),
# upper-cased with dots folded to `_`, so nlmixr2's `eff` and `EFF`, or `a.b`
# and `a_b`, are two variables in the source and one in the emitted file: the
# second assignment silently becomes a reassignment of the first. A data column
# is emitted verbatim and ferx binds it by its exact, case-sensitive header
# (E_MISSING_COVARIATE "(case-sensitive)"), so `wt` and `WT` are two columns to
# ferx as they are to the source, and are NOT a collision. A model name meets a
# column when its emitted spelling equals the column's header: a variable `dose`
# beside a column `DOSE` emits `DOSE = 3` and then reads it where the source
# read the column -- the covariate is gone, with no diagnostic from either side.
#
# Thetas, random effects and ODE states are not checked here: each has an owner
# that renames it with a message when it collides with a variable or a column
# (.deshadow_theta_names() and .uniquify_random_names(), both passed the column
# names, and .sanitise_state_names()). A state is still kept out of the column
# set -- it is read on right-hand sides like any variable. A model variable cannot be renamed the same way: two variables
# the source keeps apart cannot be told apart in the emitted text to decide
# which one a reference meant, and a column cannot be renamed at all.
#
# Runs for NONMEM too, where it should never fire: NM-TRAN is case-insensitive
# and has no dots, and nonmem2rx canonicalises (variables lower-case, data
# columns as $INPUT spells them -- `(wt/70)` arrives as `(WT/70)`). Checked
# across the bundled corpus, which raises nothing; a test pins that, so a check
# grown too broad shows up as a false refusal there rather than in a user's model.
.source_name_collisions <- function(lst, name_map) {
  nm    <- .walk_names(lst)
  model <- nm$assigned
  cols  <- setdiff(nm$used, c(model, nm$states, names(name_map)))
  key   <- .norm(model)
  groups <- list()
  for (k in unique(key)) {
    who <- c(model[key == k], cols[cols == k])
    if (length(who) > 1L) groups[[k]] <- who
  }
  warn <- character(); unsp <- character()
  for (k in names(groups)) {
    who  <- groups[[k]]
    what <- ifelse(who %in% model, "model variable", "data column")
    warn <- c(warn, paste0(
      "ERROR | ", paste0("'", who, "' (", what, ")", collapse = " and "),
      " are distinct in the source but are the same name to ferx: model names ",
      "are emitted upper-cased with dots as `_`, and a data column is read by ",
      "its exact header. The emitted file merges them into '", k,
      "', so a reference to one reads the other. Rename one of them in the source."))
    unsp <- c(unsp, paste0("source names ferx cannot tell apart: ",
                           paste(who, collapse = ", ")))
  }
  list(warnings = warn, unsupported = unsp)
}

# -- the prediction of a compartment-free model --------------------------------

# Strip `* sqrt(W)` from a prediction written for hand-built inverse-variance
# weighting, `f * sqrt(W)` or `sqrt(W) * f`, with W a bare symbol. Returns NULL
# when the expression does not have that shape.
.split_sqrt_weight <- function(e) {
  is_sqrt_sym <- function(x) is.call(x) && identical(x[[1]], as.name("sqrt")) &&
    length(x) == 2L && is.symbol(x[[2]])
  while (is.call(e) && identical(e[[1]], as.name("("))) e <- e[[2]]
  if (!is.call(e) || !identical(e[[1]], as.name("*")) || length(e) != 3L) return(NULL)
  if (is_sqrt_sym(e[[3]])) return(list(f = e[[2]], w = as.character(e[[3]][[2]])))
  if (is_sqrt_sym(e[[2]])) return(list(f = e[[3]], w = as.character(e[[2]][[2]])))
  NULL
}

# Lower a compartment-free model: the variable the error line reads becomes
# `y = ...` in [structural_model].
#
# Inline `y = <rhs>` when the prediction is assigned once, at top level, and
# nothing else reads it -- the shape the issue asks for and the one a reader
# expects. Otherwise keep the assignment(s) in [individual_parameters] and emit
# `y = <NAME>`, which ferx accepts for any individual parameter. The second
# shape is also forced when the right-hand side reads a kappa directly: ferx
# refuses a kappa in a Form C output expression ("cannot reference the IOV
# parameter"), and [individual_parameters] is the one block that may read one.
#
# `weighting = "weight"` turns `pred <- f * sqrt(W)` into `y = f` plus
# `DV ~ additive(SD) weight = 1/sqrt(W)`, for data whose DV has been divided
# back by sqrt(W) (ferx_weighted_data()). Same likelihood up to the data-only
# constant -sum(log W): ferx scores the weighted form on the natural scale. The
# weight is written as the expression rather than a precomputed `SE` column so
# the model reads exactly one column, the one the source already reads, and a
# dataset cannot carry an SE that disagrees with its W.
#
# Anything the weight form cannot represent faithfully aborts instead of falling
# back to the as-written model: the caller asked for a different data
# convention, so a silent fallback would score DV / sqrt(W) against
# f * sqrt(W).
.lower_compartment_free <- function(expr_out, kappa_names, param_names,
                                    weighting = "as_written") {
  lhs <- expr_out$tilde_lhs
  if (length(lhs) != 1L) return(NULL)
  ips <- expr_out$indiv_params
  kind <- vapply(ips, function(p) .stmt_kind(p, "assign"), "")
  lhs_of <- vapply(ips, function(p) if (identical(p$kind, "if")) "" else toupper(p$lhs), "")
  top <- which(kind != "if" & lhs_of == toupper(lhs))
  in_if <- any(vapply(ips[kind == "if"], function(p)
    toupper(lhs) %in% toupper(.ip_names(list(p))), logical(1)))
  if (length(top) == 0L && !in_if) return(NULL)

  # Who else reads the prediction. Statements are already in source order.
  others <- ips[setdiff(seq_along(ips), top)]
  read_elsewhere <- toupper(lhs) %in% toupper(.ip_refs(others)) ||
    (length(top) == 1L && toupper(lhs) %in%
       toupper(.collect_symbols(str2lang(ips[[top]]$rhs))))
  # Inlining moves the prediction from its own position to the end of the
  # model, past every later statement. That is the same value only if no later
  # statement reassigns something the prediction reads: `eff <- temax*DOSE;
  # pred <- E0 - eff; eff <- 0` predicts E0 - temax*DOSE, and `y = E0 - EFF`
  # would read the 0. Such a prediction stays where it is, as `y = PRED`.
  reassigned_after <- length(top) == 1L && top < length(ips) &&
    any(toupper(.collect_symbols(str2lang(ips[[top]]$rhs))) %in%
          toupper(.ip_names(ips[(top + 1L):length(ips)])))
  single <- length(top) == 1L && !in_if && !read_elsewhere && !reassigned_after

  warn <- character()
  weight_col <- NA_character_
  if (identical(weighting, "weight")) {
    why <- NULL
    if (!single)
      why <- paste0("'", lhs, "' is assigned more than once, inside an if, or read ",
                    "by another statement, or a name it reads is reassigned after ",
                    "it, so `* sqrt(W)` cannot be removed from one place")
    em <- expr_out$error_model
    if (is.null(why) && !(length(em) == 1L && identical(em[[1]]$type, "additive")))
      why <- "the error model is not a single additive one (ferx's `weight =` scales the additive loading only)"
    split <- if (is.null(why)) .split_sqrt_weight(str2lang(ips[[top]]$rhs))
    if (is.null(why) && is.null(split))
      why <- paste0("the prediction '", lhs, " = ", ips[[top]]$rhs,
                    "' is not written as `f * sqrt(W)` with W a data column")
    # .ip_names(), not the top-level left-hand sides: a W assigned only inside
    # an `if` arm is a model variable too, and taking it for a data column
    # declares a column that does not exist.
    if (is.null(why) && toupper(split$w) %in% toupper(c(param_names, .ip_names(ips))))
      why <- paste0("'", split$w, "' in sqrt(", split$w, ") is a model parameter or ",
                    "variable, not a data column; ferx's `weight =` may read data ",
                    "columns, TIME and thetas only")
    # A solver builtin is not a column either: `sqrt(TIME)` would declare TIME
    # in [covariates] and have ferx_weighted_data() divide DV by a column named
    # TIME that the model reads as the clock.
    if (is.null(why) && toupper(split$w) %in% .RESERVED_ODE_NAMES)
      why <- paste0("'", split$w, "' in sqrt(", split$w, ") is a ferx builtin (",
                    paste(.RESERVED_ODE_NAMES, collapse = ", "), "), not a data column")
    if (!is.null(why))
      cli::cli_abort(c(
        "{.code weighting = \"weight\"} cannot be applied: {why}.",
        i = "Translate with {.code weighting = \"as_written\"} (the default) instead."),
        call = NULL)
    weight_col <- split$w
    ips[[top]]$rhs <- paste(deparse(split$f, width.cutoff = 500L), collapse = " ")
    em[[1]]$weight <- paste0("1/sqrt(", weight_col, ")")
    expr_out$error_model <- em
    warn <- c(warn, paste0(
      "INFO  | weighting = \"weight\": the prediction is emitted without its ",
      "`* sqrt(", weight_col, ")` and the error model carries `weight = 1/sqrt(",
      weight_col, ")`. Fit it to data whose DV is divided by sqrt(", weight_col,
      ") -- ferx_weighted_data() does that -- and expect ferx's OFV to equal the ",
      "source's minus sum(log(", weight_col, "))."))
  }

  reads_kappa <- length(top) == 1L &&
    any(toupper(.collect_symbols(str2lang(ips[[top]]$rhs))) %in% toupper(kappa_names))
  if (single && !reads_kappa) {
    y <- ips[[top]]$rhs
    ips <- ips[-top]
  } else {
    y <- lhs
  }
  list(indiv_params = ips,
       structural   = list(type = "algebraic", y = y),
       error_model  = expr_out$error_model,
       weight_column = weight_col,
       warnings     = warn)
}
