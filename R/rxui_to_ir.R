#' Convert a rxode2 UI object to a ferx intermediate representation
#'
#' Accepts the rxUI S3 object returned by `rxode2::rxode2()`, `nonmem2rx::nonmem2rx()`,
#' or `monolix2rx::monolix2rx()` and converts it to a [new_ferx_ir()] ready for
#' [emit_ferx()].
#'
#' @param ui A rxUI S3 object (environment with `$iniDf` and `$lstExpr`).
#' @param source_format One of `"nonmem"`, `"nlmixr2"`, `"monolix"`, or `NA`.
#' @param source_file Path to the source file, or `NA`.
#' @param scaling_hint Named list mapping compartment number (as a character
#'   string) to the NONMEM `Sn=` scaling variable name, as returned by
#'   `.extract_nm_scaling()`. Used to emit a `[scaling]` block for ODE models
#'   where the observed compartment is scaled (e.g. `S2=V`). The default empty
#'   list disables scaling translation.
#'
#' @return A `ferx_ir` object.
#'
#' @seealso [new_ferx_ir()], [emit_ferx()], [to_ferx()]
#'
#' @importFrom stats setNames
#' @importFrom utils tail
#'
#' @examples
#' \dontrun{
#' ui <- rxode2::rxode2(function() {
#'   ini({ tvcl <- 0.134; eta.cl ~ 0.07; prop.err <- 0.01 })
#'   model({ cl <- tvcl * exp(eta.cl); linCmt() ~ prop(prop.err) })
#' })
#' ir <- rxui_to_ir(ui, source_format = "nlmixr2")
#' cat(emit_ferx(ir))
#' }
#' @export
rxui_to_ir <- function(ui, source_format = NA_character_, source_file = NA_character_,
                       scaling_hint = list()) {
  ini  <- ui$iniDf
  lst  <- ui$lstExpr
  warn <- character()
  unsp <- character()

  theta_out <- .extract_thetas(ini)
  warn      <- c(warn, theta_out$warnings)

  omega_out <- .extract_omegas(ini)
  # The flattening (ETA-coded IOV read as IIV) is a nonmem2rx behaviour, so the
  # warning -- which names nonmem2rx -- is only emitted for NONMEM sources.
  if (identical(source_format, "nonmem"))
    warn <- c(warn, .iov_flattening_warnings(omega_out$omegas))
  kappa_out <- .extract_kappas(ini)
  warn      <- c(warn, kappa_out$warnings)

  sigma_out <- .extract_sigmas(ini, tryCatch(ui$sigma, error = function(e) NULL))

  name_map  <- .norm_map_from_ini(ini)
  sigma_names_norm <- toupper(vapply(sigma_out$sigmas, function(s) s$name, ""))

  # Bind each sigma's SOURCE spelling to its emitted name. For a NONMEM source
  # sigma is not an iniDf row -- it lives in ui$sigma -- so .norm_map_from_ini()
  # never saw it and the eps symbol reached pass 3 exactly as written, matching
  # the declaration only while both were plain uppercase. See .extract_sigmas().
  for (i in seq_along(sigma_out$raw_names))
    name_map[[sigma_out$raw_names[i]]] <- sigma_out$sigmas[[i]]$name

  # A theta name and the name every reference normalises to must agree, or the
  # shadowing check below compares the wrong pair. They diverge whenever the
  # $THETA label differs from the iniDf key -- `; CL/F` gives name `t.CL` but
  # label `CL/F` -- so pin the map to the emitted name for every theta, not
  # only for the ones that end up renamed.
  theta_orig <- vapply(theta_out$thetas, function(t) t$name, "")
  for (i in seq_along(theta_orig)) name_map[[theta_out$raw_names[i]]] <- theta_orig[i]

  # Sanitise ODE state names to the ferx identifier grammar. States are not in
  # `iniDf`, so they never went through .norm() and were emitted exactly as
  # nonmem2rx spelled them -- one `c.RTOT` makes the entire file unparseable.
  # This runs BEFORE the de-shadow loop for two reasons: that loop reserves state
  # names so no theta is renamed onto one, and it must reserve the FINAL name;
  # and the odes it inspects are parsed from `name_map`, so the rename has to be
  # in the map by then or the two disagree about what the states are called.
  # Deliberately NOT the theta names, and the reason is narrower than it looks.
  # `.deshadow_theta_names()` does NOT resolve a theta/state clash -- it renames
  # only a theta that shadows an individual parameter or duplicates another
  # theta, and its `reserved` argument merely stops a replacement landing on a
  # state. Nothing renames a theta that simply shares a state's name.
  #
  # That is tolerable because ferx accepts the collision: a model with
  # `theta CENT` beside `states=[CENT]` validates clean, exact match or
  # case-differing (measured against ferx 0.2.0). Thetas are out of scope in
  # `[odes]`, so the state is what a d/dt right-hand side resolves to. Reserving
  # thetas here would rename states for no gain.
  #
  # It stops being tolerable at phase 4. `init()` and `[scaling]` DO have thetas
  # in scope, so once initial conditions are emitted a state sharing a theta's
  # name is silently shadowed there -- the defect-14 failure mode on a new pair.
  # Whoever adds init() owns that check.
  #
  # Shared with `reserved_base` below: the state sanitiser and the theta
  # de-shadower must agree about which random-effect names are spoken for, and
  # two hand-maintained copies of this list would silently drift apart the first
  # time a fourth kind of random effect is added.
  random_names <- c(unlist(lapply(omega_out$omegas, function(o) o$names)),
                    vapply(kappa_out$kappas, function(k) k$name, ""),
                    vapply(sigma_out$sigmas, function(s) s$name, ""))
  # Everything a state must not land on. Reserving too little here is not a
  # cosmetic problem: when a sanitised state collides with an assignment target,
  # that assignment is absorbed into `aux_vars`, dropped from
  # [individual_parameters], and its references resolve to the state instead --
  # `A_B = K * exp(ETA1)` beside `d/dt(A.B) = -A_B * A.B` emitted
  # `d/dt(A_B) = -A_B * A_B`, the amount squared, with the rate constant and its
  # IIV gone. ferx validates that clean (only W_UNUSED_PARAM), so nothing
  # downstream catches it.
  #
  # Assignment targets are taken raw from `lst` rather than from a parse, so they
  # are available before the model block is parsed. It over-approximates -- some
  # of those names are intermediates that get inlined away -- which only ever
  # costs a suffix on a state that would have collided with one.
  assigned_lhs <- character()
  for (e in lst) {
    if (!.is_assignment(e)) next
    if (is.symbol(e[[2]])) assigned_lhs <- c(assigned_lhs, .norm(as.character(e[[2]])))
  }
  state_out <- .sanitise_state_names(
    lst,
    taken = c(random_names,
              .covariate_names(lst, name_map),
              assigned_lhs,
              # ferx-core reserves these inside [odes] (RESERVED_ODE_NAMES); a
              # state named TIME would be shadowed by the integrator's own.
              .RESERVED_ODE_NAMES))

  # EVERY raw state name gets an entry, renamed or not. An identity entry is not
  # a no-op: `name_map` already holds every iniDf key, and .parse_model_exprs()
  # keeps extending its own copy as it walks, so a state whose raw name is also a
  # parameter key or an assignment target had its REFERENCES rewritten to that
  # other name while the declaration kept the state's own -- `d/dt(CENT) = -K *
  # VC` when a theta keyed CENT was labelled VC, and an ODE whose right-hand side
  # changed depending on whether `central <- 0` stood above or below it. Pinning
  # every state on every lookup map makes the decision taken here the one that
  # survives the walk.
  state_raw  <- .state_raw_names(lst)
  # `%in% names()`, not `is.null(map[[k]])`: `[[` on a list with a key it does
  # not hold is an error, not NULL.
  state_decl <- vapply(state_raw, function(k)
    if (k %in% names(state_out$map)) state_out$map[[k]] else k, "")

  # A state whose RAW name is also a model-parameter name is genuinely ambiguous
  # in the source: both deparse to the same symbol, and only scope separates the
  # two readings. Resolve it the way ferx does -- thetas and etas are out of
  # scope inside [odes], so the state wins there and the parameter wins
  # everywhere else -- and say so, because no .ferx file can record what the
  # source meant and the two readings are different models.
  #
  # The match is case-insensitive, which is how ferx compares them: measured
  # against ferx 0.3.0, `d/dt(central) = -central * (CENTRAL/V)` beside
  # `states=[central]` and `theta CENTRAL` validates clean and reports
  # `W_UNUSED_PARAM: theta 'CENTRAL' is ... not referenced in any model
  # expression` -- the engine read the ODE's CENTRAL as the state and the theta
  # went dead, turning the term into the amount squared over V. A case-sensitive
  # test missed that pair entirely, and the collision only becomes visible once
  # something drags the reference into [odes]: `kk <- CENTRAL/v` is an honest read
  # of the theta in [individual_parameters] scope, but its RHS "references a
  # state" by the same case-folded comparison, so the inliner moves the text into
  # the one block where that spelling means something else.
  ambiguous <- state_raw[toupper(state_raw) %in% toupper(names(name_map))]
  for (k in ambiguous) {
    warn <- c(warn, paste0(
      "ERROR | '", k, "' names both an ODE state and a model parameter in the ",
      "source. References were resolved by scope -- the state inside [odes], the ",
      "parameter everywhere else, which is how ferx reads them -- but that may ",
      "not be what the model meant. Rename one of the two in the source."))
    unsp <- c(unsp, paste0("state/parameter source-name collision: ", k))
  }

  # Two overlays, and the split IS the scope rule above. `pins` is applied to
  # every expression the parser normalises; `ode_pins` only to ODE right-hand
  # sides, because an ambiguous name means the parameter outside [odes] and the
  # state inside it. Neither is folded into `name_map`: that map is rebound as
  # the walk proceeds, and a decision taken here must not be.
  state_pins     <- state_decl[setdiff(state_raw, ambiguous)]
  state_ode_pins <- state_decl[ambiguous]
  state_arg      <- list(decl = state_decl, pins = state_pins,
                         ode_pins = state_ode_pins)
  renamed   <- if (length(state_out$map) > 0L) names(state_out$map) else character()
  state_map <- state_decl[renamed]
  warn <- c(warn, unname(state_out$warnings[renamed]))

  # De-shadow against the individual-parameter names the parser ACTUALLY
  # produces, by parsing first and parsing again with the corrected map. The
  # alternative -- predicting the parser's output -- has to re-implement its
  # filters by hand and drifts the moment they change; it over-predicted on
  # every bundled model, renaming thetas that shadowed nothing. .parse_model_exprs()
  # is pure (name_map is passed by value), so the throwaway pass is safe, and it
  # costs ~1% of a translation against the nonmem2rx parse that precedes it.
  # Iterate to a fixpoint. One pass is NOT enough: pass 3 drops a theta alias by
  # comparing names (`V <- V` is filtered as a self-assignment), so renaming a
  # theta can turn a filtered alias into an individual parameter that did not
  # exist in the previous parse -- and the surviving theta then shadows it. With
  # duplicate $THETA labels only one of the pair is renamed, which is exactly
  # how that happens. Re-parse and re-check until nothing new appears.
  reserved_base <- c(random_names, .covariate_names(lst, name_map))
  # The thetas referenced from [odes], recomputed every round: indices for the
  # carrier loop, source names for everything that reports on them.
  ode_theta_i <- integer()
  ode_theta   <- character()
  expr_out   <- .parse_model_exprs(lst, name_map, sigma_names_norm, state_arg)
  rename_why <- vector("list", length(theta_orig))
  for (round in 1:5) {
    # The linCmt passthrough invents an individual parameter for a fixed-effect
    # PK theta with no assignment of its own, so those names must be
    # de-shadowed too -- otherwise it emits the self-shadowing `V = V` this
    # function exists to prevent.
    cur_lhs <- vapply(expr_out$indiv_params, function(p) p$lhs, "")
    if (identical(expr_out$structural$type, "lincmt"))
      cur_lhs <- c(cur_lhs,
                   theta_orig[toupper(theta_orig) %in% .PK_CANDIDATES &
                              !toupper(theta_orig) %in% toupper(cur_lhs)])

    # A theta referenced from [odes], where thetas are out of scope, needs an
    # individual parameter to carry the value in. The test is whether the theta's
    # CURRENT EMITTED name appears in the emitted ODE text -- not its source name,
    # which is the individual parameter's name in every de-shadowed model:
    # `d/dt(ABS) = -KA * ABS` beside `theta TVKA` and `KA = TVKA * exp(ETA_KA)`
    # reads the parameter, and matching on the source name `KA` invented a second
    # carrier for a reference that was already correct.
    #
    # Recomputed each round rather than accumulated, for the same reason: a theta
    # matches on its source name only until de-shadowing renames it, and a set
    # carried over from that round keeps a match that is no longer true.
    #
    # Listing the name in `cur_lhs` is what triggers the rename. It is not an
    # individual parameter yet -- the next round re-parses and finds that it is,
    # which is why this must live inside the loop.
    #
    # A theta whose name is also a state name is excluded: in [odes] the state
    # wins, so the symbol is not a theta reference. Without the exclusion the ODE
    # still comes out right -- the rename moves the theta off the state's name and
    # the next round finds no match -- but the theta has been renamed for a
    # reference that was never to it, and reported as shadowing something it does
    # not shadow.
    #
    # The comparison is uppercased, and against BOTH spellings of every state.
    # `theta_orig` is `.norm()`ed and so always uppercase, while `state_decl` is
    # keyed by the RAW d/dt target -- lowercase in every nlmixr2 source -- so a
    # case-sensitive test against the keys alone never fired there: theta CENTRAL
    # beside state `central` was renamed to TVCENTRAL and reported as shadowing an
    # individual parameter that does not exist. The values matter as well as the
    # keys, because a sanitised state (`c.RTOT` -> `c_RTOT`) reaches [odes] under
    # its emitted name and it is that name a theta would collide with.
    state_syms <- toupper(c(names(state_decl), state_decl))
    ode_syms  <- .emitted_ode_symbols(expr_out$odes)
    theta_now <- vapply(theta_out$thetas, function(t) t$name, "")
    # An ODE symbol that spells a theta's SOURCE name is a second, separate case,
    # and the two are told apart by whether anything else declares that symbol.
    # nonmem2rx does not bind a theta referenced by its $THETA label -- `FLUX =
    # KTP*A(1)` for `(0,0.2) ; KTP` leaves a free `KTP` and records the theta in a
    # `rxmissingvars` placeholder -- so the symbol reaches [odes] declared nowhere
    # and has to be bound to a carrier. Where the same spelling IS declared it is
    # not a theta reference: `d/dt(ABS) = -KA * ABS` beside `KA = TVKA*exp(ETA_KA)`
    # reads the individual parameter, and treating it as a theta invented a second
    # carrier for a reference that was already correct.
    declared  <- toupper(c(cur_lhs, names(state_decl),
                           vapply(expr_out$odes, function(o) o$state, ""),
                           reserved_base))
    free_syms <- setdiff(ode_syms, declared)
    # Kept as INDICES, and converted to names only where names are what is
    # wanted. The predicate is per-theta, but two thetas can share a source name
    # -- that is exactly what a duplicate $THETA label is -- and collapsing the
    # result to `theta_orig[...]` threw away which of them matched. The carrier
    # loop then tested `theta_orig[i] %in% ode_theta` and let both through, so one
    # ODE reference produced two carriers: the one the reference resolves to, and
    # a dead `KTP = TVKTP` that ferx reports as `computed but never used`.
    ode_theta_i <- which((toupper(theta_now) %in% ode_syms |
                          toupper(theta_orig) %in% free_syms) &
                         !toupper(theta_orig) %in% state_syms)
    ode_theta   <- theta_orig[ode_theta_i]
    cur_lhs     <- c(cur_lhs,
                     ode_theta[!toupper(ode_theta) %in% toupper(cur_lhs)])
    desh <- .deshadow_theta_names(
      theta_names = vapply(theta_out$thetas, function(t) t$name, ""),
      indiv_names = cur_lhs,
      # States and covariates matter too: ferx resolves theta before both, so a
      # rename landing on either would reintroduce the shadowing on a new pair.
      # (The current theta names need not be listed -- .deshadow_theta_names()
      # already folds them into `taken`.)
      reserved    = c(reserved_base,
                      vapply(expr_out$odes, function(o) o$state, ""))
    )
    if (!any(!is.na(desh$map))) break
    for (i in seq_along(theta_out$thetas)) {
      if (is.na(desh$map[i])) next
      theta_out$thetas[[i]]$name <- desh$map[i]
      name_map[[theta_out$raw_names[i]]] <- desh$map[i]
      # Reasons accumulate; the message is written once, after the loop. A theta
      # can be renamed in more than one round (renaming one can reveal another
      # individual parameter), and reporting each hop would name intermediate
      # values that appear nowhere in the output -- or worse, name something
      # that ends up an individual parameter.
      rename_why[[i]] <- unique(c(rename_why[[i]], desh$reasons[[i]]))
    }
    expr_out <- .parse_model_exprs(lst, name_map, sigma_names_norm, state_arg)
  }

  # An [odes] reference is reported by the message the carrier emits, not as a
  # shadow. The rename is real, but at the point it happens nothing named KTP is
  # an individual parameter -- the trigger only PREDICTS one, and the prediction
  # comes true because the carrier is appended below. Calling that a shadow tells
  # the user their model had a collision it never had.
  parse_lhs <- toupper(vapply(expr_out$indiv_params, function(p) p$lhs, ""))
  ode_only  <- ode_theta[!toupper(ode_theta) %in% parse_lhs]

  for (i in seq_along(rename_why)) {
    if (length(rename_why[[i]]) == 0L) next
    final <- theta_out$thetas[[i]]$name
    if (any(rename_why[[i]] == "duplicate")) warn <- c(warn, paste0(
      "WARN  | two thetas are named '", theta_orig[i], "' (duplicate $THETA ",
      "label) -- renamed the later one to '", final, "'. ferx would have ",
      "resolved every reference to the first and silently ignored the second."))
    if (any(rename_why[[i]] == "shadow") && !theta_orig[i] %in% ode_only)
      warn <- c(warn, paste0(
        "INFO  | theta '", theta_orig[i], "' shares a name with an individual ",
        "parameter -- renamed to '", final, "' (in ferx a theta silently shadows ",
        "an identically named individual parameter)"))
    if (any(rename_why[[i]] == "builtin")) warn <- c(warn, paste0(
      "WARN  | theta '", theta_orig[i], "' collides with a ferx solver builtin (",
      paste(.RESERVED_ODE_NAMES, collapse = ", "), ") -- renamed to '", final,
      "'. ferx resolves the bare name to the builtin, so the theta would have ",
      "been declared and estimated but never read."))
  }

  # Define the passthroughs the [odes] scope now references. A theta cannot be
  # read from a d/dt right-hand side, so `KTP * CENT` is emitted as a reference to
  # an individual parameter named KTP and this is what defines it. Where the
  # source already supplies the assignment the entry is skipped: `KTP = THETA(3)`
  # arrives as the alias `KTP <- KTP`, which pass 3 drops as a self-assignment but
  # which de-shadowing turns into the surviving `KTP <- TVKTP` -- the same
  # passthrough, spelled by the model rather than by us. It is only the direct
  # reference (`DADT(1) = -THETA(3)*A(1)`, no $PK line) that has nothing to
  # convert, and that shape emitted a bare theta into [odes] until this ran.
  #
  # This runs BEFORE [scaling], and the order is load-bearing rather than
  # incidental. [scaling] resolves `S2 = VC` by looking the name up among the
  # theta names and then among the individual parameters, and a carrier moves it
  # from the first list to the second: the theta is renamed to TVVC, so `VC` no
  # longer matches there, and the parameter that answers to `VC` is the carrier.
  # Resolved first, the lookup found neither, `matched` came back NULL, and
  # [scaling] was dropped with no diagnostic -- the emitted model then predicts
  # amounts against concentration data and the engine validates it clean, which
  # is the S2=V failure CLAUDE.md warns about with the loud half removed.
  if (length(ode_theta_i) > 0) {
    theta_names <- vapply(theta_out$thetas, function(t) t$name, "")
    carrier     <- character()   # emitted theta name -> parameter carrying it
    for (i in seq_along(theta_names)) {
      nm <- theta_orig[i]
      if (!i %in% ode_theta_i) next
      lhs <- vapply(expr_out$indiv_params, function(p) p$lhs, "")
      rhs <- vapply(expr_out$indiv_params, function(p) p$rhs, "")
      # Reuse an existing parameter only when it is a PURE ALIAS of this theta.
      # Matching on the name alone is not enough and gets the arithmetic wrong:
      # `frac <- central/cl` in a model that later writes `cl <- cl*exp(eta.cl)`
      # reads the theta, so pointing that reference at the individual parameter
      # CL silently swaps in the IIV-applied value. Both forms parse.
      hit <- which(toupper(lhs) == toupper(nm) & rhs == theta_names[i])
      if (length(hit) > 0) {
        carrier[theta_names[i]] <- lhs[hit[1L]]
        next
      }
      # Otherwise define a new one, dodging everything already emitted -- including
      # the individual parameter that made the reuse invalid.
      taken <- c(lhs, theta_names, reserved_base,
                 vapply(expr_out$odes, function(o) o$state, ""))
      pt    <- .carrier_name(nm, theta_names[i], taken)
      expr_out$indiv_params <- c(expr_out$indiv_params,
                                 list(list(lhs = pt, rhs = theta_names[i])))
      carrier[theta_names[i]] <- pt
      warn <- c(warn, paste0(
        "INFO  | theta '", nm, "' is referenced from an ODE, where ferx cannot ",
        "resolve a theta -- emitting the individual parameter '", pt, " = ",
        theta_names[i], "' to carry its value"))
      # A numbered name means both preferred spellings were taken. Say so: the
      # number is positional, so it moves when another carrier appears before this
      # one, and anything indexing the emitted parameters by name breaks silently.
      if (grepl("_[0-9]+$", pt) && !grepl("_[0-9]+$", nm))
        warn <- c(warn, paste0(
          "WARN  | carrier for theta '", theta_names[i], "' had to be numbered ('",
          pt, "') because both '", nm, "' and '", theta_names[i], "_ODE' are ",
          "already in use. That name is positional and will change if another ",
          "carrier is added before it -- rename the colliding model variable."))
    }
    expr_out$odes <- .scope_odes_to_params(expr_out$odes, carrier)
  }

  structural <- expr_out$structural
  if (identical(structural$type, "ode")) {
    state_names <- vapply(expr_out$odes, function(o) o$state, "")
    obs_cmt     <- tryCatch(ui$central, error = function(e) NULL)
    # length must be checked, not just type: is.character() is TRUE for
    # character(0) and for a length-2 vector, and both then reached an `if`
    # that aborts ("argument is of length zero" / "the condition has length > 1")
    # on a model that previously translated.
    if (is.null(obs_cmt) || !is.character(obs_cmt) || length(obs_cmt) != 1L) {
      # state_names already carries the sanitised names, so the guess needs no
      # translation -- but ui$central is raw and does.
      obs_cmt <- tail(state_names, 1)
      warn <- c(warn, paste0("WARN  | obs_cmt could not be inferred -- guessed '",
                             obs_cmt, "', verify in [structural_model]"))
    } else if (obs_cmt %in% names(state_decl)) {
      # Through the state declaration map for the same reason the d/dt target is:
      # name_map would resolve a state name that happens to be a parameter key to
      # the parameter. `state_decl` covers every state, not only the renamed
      # ones, so an identity entry answers here rather than falling through.
      obs_cmt <- state_decl[[obs_cmt]]
    }
    structural$states  <- state_names
    structural$obs_cmt <- obs_cmt
  }

  scaling <- list()
  if (identical(structural$type, "ode") && length(scaling_hint) > 0L) {
    state_names_uc <- toupper(vapply(expr_out$odes, function(o) o$state, ""))
    # which() may match more than one state if two share an uppercased name;
    # take the first so the list [[ ]] index below is always scalar.
    obs_idx <- which(state_names_uc == toupper(structural$obs_cmt))[1L]
    if (!is.na(obs_idx)) {
      svar <- scaling_hint[[as.character(obs_idx)]]
      if (!is.null(svar)) {
        # .norm(), not a hand-inlined copy of it. This line used to spell out the
        # old dot-only rule; once .norm() gained the rest of the grammar the two
        # disagreed for any other illegal character, and a mismatch here does not
        # error -- `matched` comes back NULL and [scaling] is silently dropped,
        # which is precisely the S2=V failure CLAUDE.md warns about.
        norm_svar      <- .norm(svar)
        theta_names_uc <- toupper(vapply(theta_out$thetas, function(t) t$name, ""))
        indiv_lhs_uc   <- toupper(vapply(expr_out$indiv_params, function(p) p$lhs, ""))
        matched <- if (norm_svar %in% theta_names_uc) {
          vapply(theta_out$thetas, function(t) t$name, "")[match(norm_svar, theta_names_uc)]
        } else if (norm_svar %in% indiv_lhs_uc) {
          vapply(expr_out$indiv_params, function(p) p$lhs, "")[match(norm_svar, indiv_lhs_uc)]
        } else NULL
        if (!is.null(matched)) {
          scaling <- list(obs_scale = matched)
          warn    <- c(warn, paste0("INFO  | S", obs_idx, " = ", svar,
                                    " detected -- emitting [scaling] obs_scale = ", matched))
        }
      }
    }
  }

  lincmt_found <- identical(structural$type, "lincmt")
  if (lincmt_found) {
    # Fixed-effect PK params (theta with no ETA) are absent from indiv_params,
    # so the pk macro arg lookup misses them. Add passthrough entries so that
    # e.g. `V = THETA(3)` (no ETA) still produces `v=V` in the pk macro call.
    # The LHS is the PK name as written in the source (theta_orig) and the RHS
    # the possibly de-shadowed theta, so a passthrough reads `V = TVV` and never
    # the self-shadowing `V = V`.
    existing_lhs  <- toupper(vapply(expr_out$indiv_params, function(p) p$lhs, ""))
    theta_names   <- vapply(theta_out$thetas, function(t) t$name, "")
    for (i in seq_along(theta_names)) {
      pk_name <- theta_orig[i]
      if (toupper(pk_name) %in% .PK_CANDIDATES && !toupper(pk_name) %in% existing_lhs) {
        expr_out$indiv_params <- c(expr_out$indiv_params,
                                   list(list(lhs = pk_name, rhs = theta_names[i])))
        # Track the added name so a second theta normalising to the same PK
        # candidate does not append a duplicate passthrough entry.
        existing_lhs <- c(existing_lhs, toupper(pk_name))
      }
    }
    pk_out <- .infer_pk_macro(expr_out$indiv_params)
    warn   <- c(warn, pk_out$warnings)
    unsp   <- c(unsp, pk_out$unsupported)
    if (is.na(pk_out$pk_call)) {
      structural <- list()
    } else {
      structural <- list(type    = "pk_macro",
                         pk_call = pk_out$pk_call,
                         pk_args = pk_out$pk_args)
    }
  }

  # The invariant this whole dance exists to establish. ferx reports nothing
  # when it is violated, so assert it rather than trust the loop -- and do it
  # AFTER the linCmt passthrough, which appends to indiv_params and can add
  # the very `V = V` self-shadow the passthrough comment says it prevents.
  final_clash <- intersect(
    toupper(vapply(theta_out$thetas, function(t) t$name, "")),
    toupper(vapply(expr_out$indiv_params, function(p) p$lhs, "")))
  if (length(final_clash) > 0) {
    warn <- c(warn, paste0(
      "ERROR | could not give theta(s) ", paste(final_clash, collapse = ", "),
      " a name distinct from an individual parameter. In ferx a theta shadows ",
      "an identically named individual parameter silently, so this model would ",
      "fit with those parameters' individual definitions ignored."))
    unsp <- c(unsp, paste0("theta/individual-parameter name collision: ",
                           paste(final_clash, collapse = ", ")))
  }

  # A state and an individual parameter sharing a name is an E_PARSE in ferx,
  # but the translator never gets that far: the assignment is absorbed into
  # `aux_vars` (its RHS now "references a state"), dropped from
  # [individual_parameters], and inlined into itself until the depth cap. The
  # result is an ODE referencing a name nothing declares, with no diagnostic --
  # a loud dot-parse error turned into a silently deleted parameter. The state
  # sanitiser cannot prevent it (individual-parameter names are not known until
  # the model block is parsed), so assert it here, where they are.
  state_clash <- intersect(
    toupper(vapply(expr_out$odes, function(o) o$state, "")),
    toupper(vapply(expr_out$indiv_params, function(p) p$lhs, "")))
  if (length(state_clash) > 0) {
    warn <- c(warn, paste0(
      "ERROR | ", paste(state_clash, collapse = ", "), " names both an ODE state ",
      "and an individual parameter. ferx requires them to be distinct, and the ",
      "individual parameter is dropped from the output rather than emitted."))
    unsp <- c(unsp, paste0("state/individual-parameter name collision: ",
                           paste(state_clash, collapse = ", ")))
  }

  # Every name the emitted [odes] block references must be declared. This is the
  # one block where the check is unambiguous, which is why it is scoped to [odes]
  # and not applied everywhere: thetas and etas are NOT in scope here, so a bare
  # symbol can only be a state, an individual parameter, an ODE intermediate, a
  # covariate or a reserved name -- and all of those are known at this point.
  #
  # ferx does report this itself (`[odes]: RHS references undefined name(s)`), with
  # or without a dataset, but only where the engine runs: the fast PR job has no
  # ferx, and the phase-2 legality check next to this one tests the GRAMMAR, not
  # whether anything declares the name -- `KTP` is a perfectly legal identifier.
  # So the class of defect that leaks an undeclared name went unseen by every
  # engine-free tier. Two separate bugs reached the corpus that way (issue #6
  # defects 2 and 4).
  #
  # It is deliberately NOT extended to [individual_parameters]: thetas ARE in scope
  # there, so the set of legitimate names is much larger and a leftover carries far
  # less information. [odes] is where the constraint is tightest, and the set is
  # closed: ferx accepts declared states, individual parameters, ODE-block
  # intermediates and the reserved time variables, and NOTHING else. Not thetas,
  # not etas, not sigmas -- and not covariates either, which is the part that makes
  # this checkable at all. Measured against ferx 0.3.0, a covariate in an ODE RHS is
  # rejected with the remedy quoted below; the engine's advice is to pre-compute the
  # covariate-dependent term in [individual_parameters] and reference that.
  #
  # (0.2.0, the CI pin, could not be re-measured -- the local install was replaced
  # by 0.3.0 -- but no bundled model references a covariate from [odes], so nothing
  # that translated before is affected either way.)
  #
  # Because the set is closed, the covariate list is NOT consulted, and that is
  # deliberate rather than an omission: `.covariate_names()` defines a covariate as
  # a symbol nothing else binds, and rxode2's `ui$allCovs` does much the same, so
  # both classify a name the translator FAILED to bind as a legitimate covariate.
  # Measured: both call the unbound `CF` in qss_tmdd.mod and the unbound `KTP` in
  # ode_theta_ref.ctl covariates. An earlier version of this check consulted them,
  # passed the entire suite, and could not fire on either defect it was written for.
  if (length(expr_out$odes) > 0) {
    ode_declared <- toupper(c(
      vapply(expr_out$odes, function(o) o$state, ""),
      vapply(expr_out$indiv_params, function(p) p$lhs, ""),
      names(expr_out$odes_intermediates),   # NULL until phase 5 emits them
      # No function whitelist belongs here. `.collect_symbols()` recurses over
      # `expr[-1]`, so a call head is never collected -- `exp(-K*CENT)` yields
      # K and CENT and nothing else -- and a list of function names could only
      # ever whitelist ordinary identifiers that happen to share a spelling.
      # It did: with EXP/LOG/ABS/MIN/MAX/SIGN/... declared, an undeclared
      # covariate named MAX passed this check while WT was reported, and ferx
      # rejected the file we had just called clean. It blinded the scope_leak
      # arm equally, since that is computed from `ode_free`.
      .RESERVED_ODE_NAMES))
    ode_free <- setdiff(.emitted_ode_symbols(expr_out$odes), ode_declared)

    # Split the report by cause: the remedy differs. A theta, eta or sigma needs a
    # carrier; anything else is a name that resolves to nothing, which for a
    # covariate means pre-computing the term one block earlier.
    scope_leak <- intersect(ode_free, toupper(c(
      vapply(theta_out$thetas, function(t) t$name, ""), random_names)))
    unknown    <- setdiff(ode_free, scope_leak)

    if (length(scope_leak) > 0) {
      warn <- c(warn, paste0(
        "ERROR | [odes] references ", paste(scope_leak, collapse = ", "),
        ", which name a theta, eta or sigma. None of the three is in scope in a ",
        "ferx ODE block -- the value has to be carried in by an individual ",
        "parameter. ferx rejects this as E_PARSE."))
      unsp <- c(unsp, paste0("theta/eta/sigma referenced from [odes] without a ",
                             "carrier: ", paste(scope_leak, collapse = ", ")))
    }
    if (length(unknown) > 0) {
      warn <- c(warn, paste0(
        "ERROR | [odes] references ", paste(unknown, collapse = ", "),
        ", which nothing in the emitted model declares. A ferx ODE right-hand side ",
        "may only name declared states, individual parameters, ODE-block ",
        "intermediates and TIME/TAFD/TAD/MACHEPS. If one of these is a covariate, ",
        "the covariate-dependent term has to be pre-computed in ",
        "[individual_parameters] and referenced from here by that name."))
      unsp <- c(unsp, paste0("undeclared name referenced from [odes]: ",
                             paste(unknown, collapse = ", ")))
    }
  }

  # ferx-core checks RESERVED_ODE_NAMES against states, individual parameters
  # AND ODE intermediates (model_parser.rs, eq_ignore_ascii_case). The state
  # sanitiser covers the first and .deshadow_theta_names() renames a theta onto
  # a free name, but an individual parameter cannot be renamed the same way --
  # it carries the source's own name and [scaling] and [error_model] reference
  # it by that name -- so it is reported. ferx rejects the file, which is loud
  # but says nothing about which source variable is at fault; `TAD` and `T` are
  # ordinary $PK and $DES variable names in NONMEM.
  builtin_params <- intersect(
    toupper(vapply(expr_out$indiv_params, function(p) p$lhs, "")),
    .RESERVED_ODE_NAMES)
  if (length(builtin_params) > 0) {
    warn <- c(warn, paste0(
      "ERROR | ", paste(builtin_params, collapse = ", "), " names both an ",
      "individual parameter and a ferx solver builtin (",
      paste(.RESERVED_ODE_NAMES, collapse = ", "), "), which ferx-core rejects. ",
      "Rename the variable in the source model."))
    unsp <- c(unsp, paste0("individual parameter collides with a ferx builtin: ",
                           paste(builtin_params, collapse = ", ")))
  }

  # A covariate is the one name we must NOT sanitise: ferx resolves it against a
  # data column, case-sensitively, so any rewrite -- of case or of an illegal
  # character -- turns a working reference into E_MISSING_COVARIATE at fit time.
  # An illegal covariate name is therefore untranslatable rather than fixable,
  # and has to be said out loud. Computed after the de-shadow loop because a
  # rename can move a name into name_map and out of the covariate set.
  bad_covs <- Filter(Negate(.is_ferx_ident), .unmapped_symbols(lst, name_map))
  if (length(bad_covs) > 0) {
    warn <- c(warn, paste0(
      "ERROR | covariate reference(s) ", paste(bad_covs, collapse = ", "),
      " are not legal ferx identifiers. They cannot be renamed -- ferx matches ",
      "covariates to data columns by exact name -- so the data column has to be ",
      "renamed in both the dataset and the model."))
    unsp <- c(unsp, paste0("covariate name is not a legal ferx identifier: ",
                           paste(bad_covs, collapse = ", ")))
  }

  warn      <- c(warn, expr_out$warnings)
  unsp      <- c(unsp, expr_out$unsupported)

  if (!lincmt_found && length(structural) == 0 && length(expr_out$odes) == 0) {
    warn <- c(warn, "ERROR | No structural model detected -- [structural_model] section omitted")
    unsp <- c(unsp, "structural model (no linCmt() or d/dt() found in model block)")
  }

  fit_opts <- list(method = "focei", maxiter = 500L, covariance = TRUE)
  if (length(kappa_out$kappas) > 0)
    fit_opts$iov_column <- kappa_out$iov_column

  # CLAUDE.md: every translation warning is emitted at translation time so the
  # user sees it immediately, not only when they inspect result$warnings.
  .emit_warnings(warn)

  new_ferx_ir(
    source_format = source_format,
    source_file   = source_file,
    thetas        = theta_out$thetas,
    omegas        = omega_out$omegas,
    kappas        = kappa_out$kappas,
    sigmas        = sigma_out$sigmas,
    indiv_params  = expr_out$indiv_params,
    structural    = structural,
    odes          = expr_out$odes,
    error_model   = expr_out$error_model,
    scaling       = scaling,
    fit_options   = fit_opts,
    warnings      = warn,
    unsupported   = unsp,
    state_renames = state_map
  )
}

# -- name normalisation -------------------------------------------------------

# ferx identifiers are [A-Za-z_][A-Za-z0-9_]* (ferx-core model_parser.rs, the
# `identifier` rule). Anything else is an E_PARSE at every reference site --
# though NOT necessarily where it is declared: a dotted name is accepted inside
# `ode(states=[...])` and rejected the moment it is used, so "the states line
# parsed" is not evidence the name is legal.
.is_ferx_ident <- function(nm) grepl("^[A-Za-z_][A-Za-z0-9_]*$", nm)

# Map any name onto that grammar. Case is preserved -- callers that want the
# uppercase model convention go through .norm(); covariate references must NOT,
# because ferx matches data columns case-sensitively (see .covariate_names).
.ferx_ident <- function(nm) {
  x <- gsub("[^A-Za-z0-9_]", "_", as.character(nm))
  # NA has to go first: nzchar(NA) is TRUE and gsub/sub leave NA alone, so an NA
  # would slip past both guards below and out of a function documented to return
  # a legal identifier for any input.
  x[is.na(x)]   <- "X"
  # A leading digit cannot be fixed by substitution, and an empty name has
  # nothing to substitute. Both need a prefix rather than a replacement.
  x[!nzchar(x)] <- "X"
  sub("^([0-9])", "X_\\1", x)
}

# Strip the nonmem2rx t. prefix (theta) and e. prefix (effect eta) before normalising.
#
# .norm() is .ferx_ident() plus the uppercase model-name convention. It used to
# substitute only the dot, which covered every illegal character the bundled
# corpus happens to contain and none of the others.
.norm <- function(nm) toupper(.ferx_ident(nm))

.strip_prefix <- function(nm) sub("^[te][.]", "", nm)

.norm_map_from_ini <- function(ini) {
  nms <- unique(ini$name[!is.na(ini$name)])
  # Map the raw iniDf name (e.g. "t.TVCL") to the normalised ferx name ("TVCL").
  setNames(vapply(nms, function(nm) .norm(.strip_prefix(nm)), ""), nms)
}

# Recursively substitute known parameter names in an expression.
# Does NOT touch the function-name position of a call (call[[1]]).
.normalise_expr <- function(expr, map) {
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    return(as.name(if (nm %in% names(map)) map[[nm]] else nm))
  }
  if (!is.call(expr)) return(expr)
  as.call(c(list(expr[[1]]),
            lapply(as.list(expr[-1]), .normalise_expr, map = map)))
}

# Overwrite `map` with `pins`, unconditionally. State renames are folded in this
# way rather than written into `name_map` because that map is rebound while the
# statements are walked, and a name decided before the walk must not change
# during it.
.pin_names <- function(map, pins) {
  for (k in names(pins)) map[[k]] <- pins[[k]]
  map
}

# nonmem2rx bookkeeping placeholders (`rxmissingvars1 <- t.KTP`). Two passes have
# to agree that these are not model variables -- the pre-seed and the walk -- so
# the test lives in one place.
.is_rxmissingvars <- function(nm) grepl("^rxmissingvars[0-9]+$", tolower(nm))

# Route stored warnings to the console by their severity prefix. The text is
# passed as a cli value, never as a format string -- a model name containing
# braces would otherwise be evaluated as R code.
.emit_warnings <- function(warn) {
  for (w in warn) {
    if (startsWith(w, "ERROR") || startsWith(w, "WARN")) cli::cli_warn("{w}")
    else                                                 cli::cli_inform("{w}")
  }
  invisible(NULL)
}

# -- theta / individual-parameter de-shadowing --------------------------------

# ferx resolves a bare identifier as theta first, then eta, then individual
# parameter (ferx-core parse_atom). A theta whose name matches an individual
# parameter therefore SHADOWS it everywhere thetas are in scope --
# [individual_parameters], [scaling], [initial_conditions] -- so
#   CL = CL * exp(ETA_CL)
#   K20 = CL / V          # reads the THETA CL, not the line above
# leaves the first line dead and the model fits without that IIV. ferx emits no
# diagnostic for it. [odes] is unaffected (thetas are not in scope there) and so
# are `pk ...(cl=CL)` arguments (resolved by direct name lookup), which is why
# the same file can carry correct ODEs and a silently wrong derived parameter.
#
# nonmem2rx names thetas after the $PK parameters whenever $THETA carries
# labels, so the collision is the common case, not an edge case. Rename the
# theta, not the individual parameter: `theta TVCL` + `CL = TVCL * ...` is the
# idiomatic ferx form and leaves the model-facing names untouched.
#
# The rename is applied to `name_map` BEFORE the model expressions are parsed,
# never to the emitted strings afterwards. A textual pass cannot tell the two
# apart: in `cl <- t.CL * exp(eta2)` the RHS means the theta, but in the
# following `k20 <- cl/v` it means the individual parameter, and both deparse
# to the same token once normalised.

# PK parameter names the linCmt pk-macro arg lookup recognises. Shared with the
# passthrough in rxui_to_ir() so both agree on which fixed-effect thetas become
# individual parameters.
.PK_CANDIDATES <- c("CL", "V", "V1", "V2", "V3", "Q", "Q2", "Q3", "KA")

# Names the solver injects into every ferx model, and the only names an [odes]
# right-hand side may use without declaring them. ferx-core's RESERVED_ODE_NAMES
# (src/parser/model_parser.rs) is compared case-insensitively against states,
# individual parameters and ODE intermediates alike, so every producer of an
# emitted name has to agree on the list -- one copy, not one per call site.
#
# The two halves fail differently. A state or individual parameter that collides
# is an outright E_PARSE; a THETA that collides is not rejected at all, because
# ferx resolves the bare name to the builtin -- `KA = TIME * exp(ETA1)` makes KA
# read the integrator clock and leaves the estimated theta unreferenced, and the
# only diagnostic is a W_UNUSED_PARAM that says nothing about why.
.RESERVED_ODE_NAMES <- c("TIME", "T", "TAFD", "TAD", "MACHEPS")

# Names referenced by the model but never assigned and never declared in iniDf.
# ferx reads every such name as a covariate (a data column), and resolves theta
# before covariate -- so a theta renamed onto one shadows it exactly as it would
# an individual parameter.
.covariate_names <- function(lst, name_map) {
  assigned <- character()
  used     <- character()
  for (expr in lst) {
    if (!.is_assignment(expr)) {
      if (.is_tilde(expr)) used <- c(used, .collect_symbols(expr))
      next
    }
    lhs <- expr[[2]]
    if (is.symbol(lhs)) assigned <- c(assigned, .norm(as.character(lhs)))
    else if (.is_ddt_lhs(lhs)) assigned <- c(assigned, .norm(.ddt_state(lhs)))
    used <- c(used, .collect_symbols(expr[[3]]))
  }
  # name_map keys are the RAW iniDf names (`t.CL`, `e.ETA1`) and .norm() does not
  # strip those prefixes, so a reference to `t.CL` normalises to `T_CL` and would
  # be reported as a covariate unless the keys are treated as known too.
  known <- unique(c(assigned,
                    toupper(unlist(name_map, use.names = FALSE)),
                    .norm(names(name_map))))
  setdiff(unique(.norm(used)), known)
}

# Symbols that reach the output EXACTLY as the source spelled them.
#
# `.normalise_expr()` rewrites a symbol only if it is a key of `name_map`;
# everything else is emitted verbatim. That set is mostly covariates -- which is
# deliberate, since ferx matches data columns case-sensitively and rewriting one
# turns a working reference into E_MISSING_COVARIATE -- but it also contains any
# name the model simply never defines.
#
# Framing the legality check on THIS set rather than on the covariate set
# matters: `.covariate_names()` classifies by normalised name, so a raw `c.RTOT`
# alongside an eta named `c_RTOT` normalises to a known `C_RTOT` and is filtered
# out as "not a covariate" -- while `.normalise_expr()`, which matches raw keys,
# still emits `c.RTOT` verbatim into an unparseable file.
.unmapped_symbols <- function(lst, name_map) {
  used <- character()
  for (expr in lst) {
    if (.is_tilde(expr))          used <- c(used, .collect_symbols(expr))
    else if (.is_assignment(expr)) used <- c(used, .collect_symbols(expr[[3]]))
  }
  # Assignment targets are emitted through .norm(), so they are legal by
  # construction and are not free symbols. That is true regardless of where the
  # reference stands relative to the assignment only because
  # .parse_model_exprs() seeds every target before it walks; when the alias was
  # installed mid-walk, a reference ABOVE its assignment was emitted raw and
  # this subtraction hid it.
  assigned <- character()
  for (expr in lst) {
    if (!.is_assignment(expr)) next
    lhs <- expr[[2]]
    if (is.symbol(lhs))               assigned <- c(assigned, as.character(lhs))
    else if (.is_ddt_lhs(lhs))        assigned <- c(assigned, .ddt_state(lhs))
  }
  setdiff(unique(used), c(names(name_map), assigned))
}

# The raw `d/dt()` targets, in source order. One owner: `rxui_to_ir()` needs the
# same list the sanitiser works from -- it pins every state, not only the renamed
# ones -- and a second hand-rolled walk would be free to disagree about what
# counts as a state.
.state_raw_names <- function(lst) {
  raw <- character()
  for (expr in lst) {
    if (!.is_assignment(expr)) next
    if (.is_ddt_lhs(expr[[2]])) raw <- c(raw, .ddt_state(expr[[2]]))
  }
  unique(raw)
}

# Give every ODE state a name that is a legal ferx identifier and collides with
# nothing, and report the mapping so every REFERENCE can be rewritten with it.
#
# Renaming the declaration alone is not a fix: `nonmem2rx` prefixes `c.` onto a
# compartment whose name collides with a variable ($MODEL COMP=(RTOT) beside an
# $ERROR RTOT gives the state `c.RTOT`), and that name then appears in
# `ode(states=[...])`, in `obs_cmt=`, on the `d/dt` left-hand side, AND inlined
# into other ODE right-hand sides wherever the source wrote `A(3)`. The caller
# turns the returned map into the pins `.parse_model_exprs()` applies, so the
# declaration and every reference are rewritten from one decision.
#
# Only names that MUST change do change. `.ferx_ident()` preserves case rather
# than going through `.norm()`, because `depot`/`central` are already legal and
# renaming them to DEPOT/CENTRAL would rewrite names the user reads and indexes
# by, to no benefit. Collision detection is nonetheless case-INSENSITIVE: ferx
# rejects a state that matches an individual parameter or an ODE intermediate
# case-insensitively.
.sanitise_state_names <- function(lst, taken = character()) {
  raw  <- .state_raw_names(lst)
  map  <- list()
  warn <- character()
  # A state may legally keep its own name, so it must not be treated as taken by
  # itself -- seed `used` with everything EXCEPT the states, and add each state
  # as it is resolved.
  used <- taken

  # Two passes, already-legal names first. A state whose source name needs no
  # change has first claim on it; otherwise an illegal name that sanitises onto
  # the same spelling can take it first and displace a name that never needed to
  # move -- `A.B` becoming `A_B` renamed an existing `A_B` to `A_B_1`, churning
  # the one name in the pair that was fine.
  legal <- vapply(raw, .is_ferx_ident, logical(1))

  for (r in raw[c(which(legal), which(!legal))]) {
    cand <- .free_name(.ferx_ident(r), used)
    used <- c(used, cand)
    if (identical(cand, r)) next

    map[[r]] <- cand
    # Keyed by the raw state name so a caller that declines to apply one rename
    # can drop its message with it, rather than telling the user about a rename
    # that did not happen.
    warn[[r]] <- if (!.is_ferx_ident(r))
      paste0("INFO  | state '", r, "' is not a legal ferx identifier -- renamed",
             " to '", cand, "' (ferx names are letters, digits and underscore,",
             " not starting with a digit)")
    else
      paste0("INFO  | state '", r, "' renamed to '", cand,
             "' -- the source name collides with another emitted name")
  }

  # Only the renames. Callers that need the COMPLETE set of states -- and they
  # do, because inside [odes] a symbol naming a state must resolve to that state
  # even when it kept its own name -- call `.state_raw_names()`, which is the one
  # owner of "what counts as a state" and is what this function reads too.
  list(map = map, warnings = warn)
}

# Pick a name that is free in `taken`, comparing case-insensitively because
# that is how ferx compares them. `prefer` is tried first, in order; failing
# that, `base` is suffixed _1, _2, ... until something is free.
#
# `taken` is folded here rather than by the caller. The two users of this used
# to carry separate copies of the suffix search with OPPOSITE conventions --
# .free_theta_name() required an already-uppercased `taken`, .sanitise_state_names()
# uppercased its own -- which is exactly the sort of difference that survives a
# refactor and silently stops matching.
.free_name <- function(base, taken, prefer = character()) {
  taken <- toupper(taken)
  for (cand in c(prefer, base))
    if (!toupper(cand) %in% taken) return(cand)
  i <- 1L
  repeat {
    cand <- paste0(base, "_", i)
    if (!toupper(cand) %in% taken) return(cand)
    i <- i + 1L
  }
}

# Pick a free replacement for a shadowed theta.
.free_theta_name <- function(old, taken) {
  .free_name(old, taken, prefer = c(paste0("TV", old), paste0("THETA_", old)))
}

# Name the individual parameter that carries a theta's value into [odes].
#
# First choice is the name the source used (`KTP`), which is what ferx's own
# bundled examples do (`KM = TVKM`) and what keeps the emitted [odes] diffable
# against the source $DES.
#
# Second choice suffixes the theta's EMITTED name: `TVCL_ODE = TVCL`. Not the
# source name -- `CL_ODE` reads as "the CL used in the ODE", which is the
# individual value, and the entire defect class this carrier exists for is
# confusing the theta value with the IIV-applied one. The name has to make that
# distinction loud. It also needs no counter, because theta emitted names are
# already unique among thetas.
#
# `_1` is deliberately NOT the second choice even though `.free_name()` would
# supply it. That suffix already means "state disambiguated" (`central_1`) and is
# `.free_theta_name()`'s last resort, so a third meaning on the same spelling
# leaves a reader unable to tell a renamed compartment from a renamed theta from a
# carrier. `CL_1` and `V_1` are also plausible model variables in their own right.
# `base` is the suffixed form, not the source name, so that even the numbered last
# resort stays recognisable as a carrier: `TVCL_ODE_1`, never `CL_1`.
#
# The solver builtins are folded in HERE rather than left to the caller, for the
# reason `.free_name()` folds `taken`: the constraint has to travel with the
# function or the next call site forgets it. Without them the first choice is the
# theta's SOURCE name, which is precisely the name `.deshadow_theta_names()` has
# just renamed the theta away from when that name is a builtin -- a theta labelled
# TIME became `theta TVTIME` and then the carrier `TIME = TVTIME`, putting the
# collision back one block below where it was removed. ferx resolves the bare TIME
# to the integrator's clock, so the ODE term silently becomes time-dependent;
# `builtin_params` does report it, but as a source defect the user cannot act on,
# since no variable of that name exists in their model.
.carrier_name <- function(source_name, theta_name, taken) {
  .free_name(paste0(theta_name, "_ODE"), c(taken, .RESERVED_ODE_NAMES),
             prefer = source_name)
}

# Give every theta a name that is unique among thetas and does not shadow a
# predicted individual-parameter name. Both failures are silent in ferx: a
# shadowed individual parameter is written and never read, and a duplicate theta
# name resolves every reference to the first while the second sits dead.
#
# The result is keyed by theta INDEX, not by name, precisely because names can
# arrive duplicated -- a name-keyed map would collapse two distinct thetas onto
# one replacement and leave the duplication in place. `rxui_to_ir()` applies it
# to both the emitted name and `name_map`; renaming in one without the other
# leaves references pointing at a name nothing declares.
.deshadow_theta_names <- function(theta_names, indiv_names,
                                  reserved = character()) {
  n   <- length(theta_names)
  map <- rep(NA_character_, n)
  if (n == 0L) return(list(map = map, reasons = list(), warnings = character()))

  shadowed <- toupper(theta_names) %in% toupper(indiv_names)
  duped    <- duplicated(toupper(theta_names))
  # The third silent failure, and the quietest of the three: ferx resolves a bare
  # `TIME`/`T`/`TAD`/`TAFD`/`MACHEPS` to the value the solver injects, so a theta
  # of that name is not rejected -- it is simply never read, and every expression
  # that names it gets the clock instead. Renaming is safe here for the same
  # reason it is for a shadowing theta: this function owns theta naming, and
  # `theta TVTIME` + `KA = TVTIME * ...` is the same model under another label.
  builtin  <- toupper(theta_names) %in% .RESERVED_ODE_NAMES
  todo     <- which(shadowed | duped | builtin)
  if (length(todo) == 0L)
    return(list(map = map, reasons = vector("list", n), warnings = character()))

  taken   <- toupper(c(theta_names, indiv_names, reserved, .RESERVED_ODE_NAMES))
  reasons <- vector("list", n)
  for (i in todo) {
    new    <- .free_theta_name(theta_names[i], taken)
    taken  <- c(taken, toupper(new))
    map[i] <- new
    # Both conditions are recorded when both hold: they are different defects
    # with different consequences, and the shadowing one is the whole point of
    # this function. The caller turns these into prose once the final name is
    # known, since a theta may be renamed again in a later round.
    if (duped[i])    reasons[[i]] <- c(reasons[[i]], "duplicate")
    if (shadowed[i]) reasons[[i]] <- c(reasons[[i]], "shadow")
    if (builtin[i])  reasons[[i]] <- c(reasons[[i]], "builtin")
  }
  list(map = map, reasons = reasons, warnings = character())
}

# -- iniDf extractors ---------------------------------------------------------

.extract_thetas <- function(ini) {
  rows <- ini[!is.na(ini$ntheta) & is.na(ini$err), , drop = FALSE]
  thetas <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, ]
    # Use label only if it is a single valid identifier (no whitespace).
    # The label must be a legal ferx identifier, not merely whitespace-free:
    # `; CL/F` would otherwise be emitted as `theta CL/F(...)`, which the engine
    # cannot parse, and would diverge from the name every reference resolves to.
    lbl_raw <- if ("label" %in% names(row) && !is.na(row$label) &&
                   grepl("^[A-Za-z][A-Za-z0-9_.]*$", as.character(row$label)))
      as.character(row$label)
    else
      .strip_prefix(row$name)
    nm <- .norm(lbl_raw)
    list(name = nm, init = row$est, lower = row$lower, upper = row$upper,
         fixed = isTRUE(row$fix))
  })

  # raw_names keeps the iniDf key for each theta so a later rename can be
  # written back into the normalisation map, not just into the emitted name.
  list(thetas = thetas, raw_names = as.character(rows$name),
       warnings = character())
}

.extract_omegas <- function(ini) {
  iiv <- ini[!is.na(ini$neta1) & ini$condition == "id", , drop = FALSE]
  if (nrow(iiv) == 0) return(list(omegas = list()))

  off   <- iiv[iiv$neta1 != iiv$neta2, , drop = FALSE]
  diag  <- iiv[iiv$neta1 == iiv$neta2, , drop = FALSE]

  if (nrow(off) == 0) {
    omegas <- lapply(seq_len(nrow(diag)), function(i) {
      list(type = "diagonal", names = .norm(.strip_prefix(diag$name[i])), values = diag$est[i])
    })
    return(list(omegas = omegas))
  }

  blocks        <- .detect_blocks(off)
  block_eta_set <- unlist(blocks)
  omegas        <- list()

  for (bg in blocks) {
    lt    <- iiv[iiv$neta1 %in% bg & iiv$neta2 %in% bg, , drop = FALSE]
    lt    <- lt[order(lt$neta1, lt$neta2), ]
    nms   <- vapply(bg, function(e) {
      row <- lt[lt$neta1 == e & lt$neta2 == e, , drop = FALSE]
      .norm(.strip_prefix(row$name[1]))
    }, "")
    omegas <- c(omegas, list(list(type = "block", names = nms, values = lt$est)))
  }

  for (i in seq_len(nrow(diag))) {
    if (!diag$neta1[i] %in% block_eta_set)
      omegas <- c(omegas, list(
        list(type = "diagonal", names = .norm(.strip_prefix(diag$name[i])), values = diag$est[i])
      ))
  }

  list(omegas = omegas)
}

# Union-find over off-diagonal eta pairs; returns list of sorted integer vectors.
.detect_blocks <- function(off) {
  eta_set <- sort(unique(c(off$neta1, off$neta2)))
  parent  <- setNames(as.list(eta_set), as.character(eta_set))

  .find <- function(x) {
    while (!identical(parent[[as.character(x)]], x))
      x <- parent[[as.character(x)]]
    x
  }
  .union <- function(a, b) {
    ra <- .find(a); rb <- .find(b)
    if (!identical(ra, rb)) parent[[as.character(ra)]] <<- rb
  }

  for (i in seq_len(nrow(off))) .union(off$neta1[i], off$neta2[i])

  roots  <- vapply(eta_set, .find, 1.0)
  groups <- split(eta_set, roots)
  lapply(groups, sort)
}

.extract_kappas <- function(ini) {
  iov  <- ini[!is.na(ini$neta1) & ini$condition != "id", , drop = FALSE]
  warn <- character()
  if (nrow(iov) == 0) return(list(kappas = list(), iov_column = NULL, warnings = warn))

  diag <- iov[iov$neta1 == iov$neta2, , drop = FALSE]
  off  <- iov[iov$neta1 != iov$neta2, , drop = FALSE]

  if (nrow(off) > 0)
    warn <- c(warn, "WARN  | IOV block omega detected -- only diagonal kappas emitted")

  iov_col <- unique(iov$condition)
  if (length(iov_col) > 1)
    warn <- c(warn, paste0("WARN  | Multiple IOV condition columns: ",
                           paste(iov_col, collapse = ", "), " -- using first"))
  iov_col <- iov_col[1]

  kappas <- lapply(seq_len(nrow(diag)), function(i)
    list(name = .norm(.strip_prefix(diag$name[i])), value = diag$est[i])
  )
  list(kappas = kappas, iov_column = iov_col, warnings = warn)
}

# Detect etas that look like inter-occasion variability (IOV) but landed in the
# IIV [omega] block. NONMEM IOV coded as an extra ETA (`KAPPA = ETA(n)`, no
# `$OMEGA ... SAME`) is read as ordinary IIV by nonmem2rx, so the occasion
# structure is silently lost. Flag the conventional KAPPA*/IOV* names so the
# user can restore it as a `kappa` declaration. Returns one WARN string per
# matching eta (character(0) when none match). `omegas` is the list returned by
# .extract_omegas(); each element's `names` is one eta (diagonal) or several
# (block).
.iov_flattening_warnings <- function(omegas) {
  nms <- unlist(lapply(omegas, function(o) o$names), use.names = FALSE)
  iov <- unique(nms[grepl("^(KAPPA|IOV)", nms, ignore.case = TRUE)])
  vapply(iov, function(nm) paste0(
    "WARN  | ETA '", nm, "' looks like inter-occasion variability but was ",
    "emitted as IIV (nonmem2rx reads ETA-coded IOV as IIV). If this is IOV, ",
    "declare it as 'kappa ", nm, " ~ ...' and set 'iov_column' in [fit_options]."
  ), character(1), USE.NAMES = FALSE)
}

# `raw_names` is returned alongside, positionally matched to `sigmas`, because
# the emitted name and the name the MODEL TEXT spells are two different strings
# and both are needed. The caller binds one to the other in `name_map` so
# .normalise_expr() rewrites the eps reference; without that the declaration and
# the reference agree only by accident of both being plain uppercase, and the
# moment the name needs sanitising (`rownames(ui$sigma) = "eps.1"`) they diverge
# and pass 3 finds no sigma in the error assignment -- emitting `sigma EPS_1`
# and NO [error_model] block at all, a model with no residual error and no
# diagnostic. For the nlmixr2 branch the binding already exists (sigma is an
# iniDf row, so .norm_map_from_ini() covered it); seeding it again is a no-op.
.extract_sigmas <- function(ini, ui_sigma = NULL) {
  # nlmixr2 / rxode2 native: sigma rows appear in iniDf with err != NA.
  rows <- ini[!is.na(ini$err), , drop = FALSE]
  if (nrow(rows) > 0) {
    sigmas <- lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, ]
      list(name = .norm(.strip_prefix(row$name)), value = row$est, scale = "sd")
    })
    return(list(sigmas = sigmas, raw_names = as.character(rows$name)))
  }
  # nonmem2rx: sigma lives in the ui$sigma matrix (variance scale); convert to SD.
  if (!is.null(ui_sigma) && is.matrix(ui_sigma) && nrow(ui_sigma) > 0) {
    nms    <- rownames(ui_sigma)
    sigmas <- lapply(seq_along(nms), function(i)
      # .norm(), matching the nlmixr2 branch above. This was the last emitted-name
      # channel that only uppercased: the name comes from rownames(ui$sigma),
      # which is upstream-controlled, and it flows into both [parameters] and
      # [error_model].
      list(name  = .norm(nms[i]),
           value = sqrt(ui_sigma[i, i]),
           scale = "sd")
    )
    return(list(sigmas = sigmas, raw_names = nms))
  }
  list(sigmas = list(), raw_names = character())
}

# -- expression classifiers ---------------------------------------------------

.is_tilde <- function(expr) {
  is.call(expr) && identical(expr[[1]], quote(`~`))
}

.is_lincmt_tilde <- function(expr) {
  .is_tilde(expr) && is.call(expr[[2]]) &&
    identical(expr[[2]][[1]], as.name("linCmt"))
}

.is_assignment <- function(expr) {
  is.call(expr) && as.character(expr[[1]]) %in% c("<-", "=", "->", "<<-")
}

# Detect  d/dt(STATE) <- rhs   (assignment whose LHS is d/dt(...))
# nonmem2rx uses <- ; rxode2 native uses =
# Parsed by R as: `<-`(d/dt(STATE), rhs) where d/dt(STATE) = `/`(d, dt(STATE))
.is_ddt_lhs <- function(lhs) {
  is.call(lhs) &&
    identical(lhs[[1]], as.name("/")) &&
    length(lhs) >= 3 &&
    identical(lhs[[2]], as.name("d")) &&
    is.call(lhs[[3]]) &&
    identical(lhs[[3]][[1]], as.name("dt"))
}

# [[1L]] because as.character() on a call returns one element per part --
# `d/dt(f(x))` would otherwise yield c("f", "x"), which turns the state lookup's
# `if` into "the condition has length > 1" and invents a second state.
.ddt_state <- function(lhs) as.character(lhs[[3]][[2]])[[1L]]

# Point every theta reference in the emitted [odes] block at the individual
# parameter that carries its value.
#
# Applied to the FINAL right-hand sides rather than during the walk, because that
# is the only place every path converges. The d/dt line is normalised against a
# map that could carry the substitution, but an ODE-block intermediate that
# touches a state is inlined instead of emitted, and its text was normalised for
# a different context -- `ki <- KTP*CENT` reached [odes] as `TVKTP * CENT`, a
# bare theta, with the passthrough parameter sitting unreferenced beside it.
# Rewriting the emitted string closes both paths at once.
#
# `carrier` maps emitted theta name -> the individual parameter carrying it.
# Substitution is on the parsed expression, not the text: a regex for `TVKTP`
# also matches `TVKTP2`. It is re-parsed from `rhs` rather than taken from
# `rhs_expr`, which the inlining pass does not carry forward -- reading it gave
# every ODE the right-hand side `NULL`.
.scope_odes_to_params <- function(odes, carrier) {
  carrier <- carrier[names(carrier) != unname(carrier)]
  if (length(carrier) == 0L) return(odes)
  lapply(odes, function(o) {
    e <- tryCatch(str2lang(o$rhs), error = function(e) NULL)
    if (is.null(e)) return(o)
    e     <- .normalise_expr(e, as.list(carrier))
    # Only `rhs` is written. The inlining pass rebuilds each ode as
    # `list(state, rhs)` and its one consumer runs before this, so `rhs_expr` is
    # already gone by the time we get here -- setting it would put the field back
    # on carrier models and leave it NULL on every other one, and a phase-5
    # reader would get a correct expression from some models and nothing from
    # the rest.
    o$rhs <- paste(deparse(e, width.cutoff = 500L), collapse = " ")
    o
  })
}

# Every name the emitted [odes] block will reference, uppercased.
#
# Read off the FINAL rhs strings rather than the parsed expressions, because that
# text is exactly what ferx sees: intermediates have already been inlined, so a
# name that was inlined away is correctly absent and a name that only appears
# after inlining is correctly present.
#
# This is the set of names for which a theta is NOT in scope. Measured against
# ferx 0.2.0: a theta can be referenced from [individual_parameters], from
# `[scaling] y` and from `obs_scale`, but NOT from a d/dt right-hand side, an
# ODE-block intermediate, an `init()` expression, or a pk macro argument. The
# plan's phase-3 text lists [scaling] among the out-of-scope blocks; it is not.
#
# init() expressions belong in this set too and will join it when phase 4 starts
# emitting them -- they are parsed from the same odes list.
.emitted_ode_symbols <- function(odes) {
  if (length(odes) == 0L) return(character())
  syms <- unlist(lapply(odes, function(o) {
    e <- tryCatch(str2lang(o$rhs), error = function(e) NULL)
    if (is.null(e)) character() else .collect_symbols(e)
  }))
  unique(toupper(syms))
}

# Collect all symbol names (leaves) from an expression tree.
.collect_symbols <- function(expr) {
  if (is.symbol(expr)) return(as.character(expr))
  if (!is.call(expr))  return(character())
  unlist(lapply(as.list(expr[-1]), .collect_symbols))
}

# Recursively substitute aux-var symbols in an expression with their definitions.
# aux_map: named list mapping uppercase symbol name -> defining R expression.
# Used to inline $DES-internal intermediates (e.g. C2, EFF) into ODE RHS strings
# so that the emitted [odes] block references only thetas, etas, and states.
.inline_aux_vars <- function(expr, aux_map, depth = 0L) {
  if (depth > 30L) return(expr)
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    if (nm %in% names(aux_map))
      return(.inline_aux_vars(aux_map[[nm]], aux_map, depth + 1L))
    return(expr)
  }
  if (!is.call(expr)) return(expr)
  as.call(c(list(expr[[1]]),
            lapply(as.list(expr[-1]), .inline_aux_vars,
                   aux_map = aux_map, depth = depth + 1L)))
}

# -- expression parser --------------------------------------------------------

#' @param states List of three named character vectors describing the ODE
#'   states, all decided by `rxui_to_ir()` BEFORE this function runs and none of
#'   them folded into `name_map`. `decl` is raw name -> emitted name for every
#'   state, and answers the `d/dt` target. `pins` and `ode_pins` are the same
#'   entries split by scope and re-applied on top of `name_map` at every lookup.
#'
#'   The split matters because `name_map` is rebound as the walk proceeds --
#'   every ordinary assignment installs an alias -- and it already holds every
#'   iniDf key. Resolving a state through it therefore produced a name nobody
#'   chose: `d/dt(CENT) = -K * VC` when a theta keyed CENT was labelled VC, and
#'   an ODE whose right-hand side changed depending on whether `central <- 0`
#'   stood above or below it. `pins` covers the states whose raw name is
#'   unambiguous and is applied everywhere; `ode_pins` covers a state whose raw
#'   name is also a parameter name and is applied only to ODE right-hand sides,
#'   which is exactly ferx's own scoping -- thetas and etas are out of scope
#'   inside `[odes]`, so the state wins there and the parameter wins elsewhere.
#'
#'   Theta references are NOT resolved here -- see `.scope_odes_to_params()`,
#'   which rewrites them on the emitted text, because an ODE-block intermediate
#'   is inlined and never passes through this function's `d/dt` branch.
#' @noRd
.parse_model_exprs <- function(lst, name_map, sigma_names = character(),
                               states = list()) {
  state_decl     <- if (is.null(states$decl))     character() else states$decl
  state_pins     <- if (is.null(states$pins))     character() else states$pins
  state_ode_pins <- if (is.null(states$ode_pins)) character() else states$ode_pins

  # Pass 1: collect assignments; handle d/dt, linCmt, tilde directly.
  all_assigns  <- list()   # list(lhs_norm, rhs_norm, rhs_expr)
  odes         <- list()
  error_model  <- list()
  structural   <- list()
  warnings     <- character()
  unsupported  <- character()

  # Seed every ordinary assignment target before the walk, so a reference that
  # appears BEFORE its defining assignment resolves to the same emitted name as
  # one after it. The alias used to be installed only once the walk had passed
  # the assignment, so `d/dt(central) = -cl*central*f.rac` written above
  # `f.rac <- 0.5` emitted the illegal `f.rac` verbatim beside `F_RAC = 0.5` --
  # two spellings of one variable, an unparseable file, and nothing reported it
  # (.unmapped_symbols() subtracts assignment targets regardless of position).
  # Existing keys are never overwritten: where the name is also an iniDf key the
  # map already holds the binding that reference means, and the walk rebinds it
  # afterwards, which is what de-shadowing `cl <- cl * exp(eta.cl)` depends on.
  # nonmem2rx's rxmissingvars placeholders are skipped by the walk below, so
  # seeding them here would put a name in the map that nothing ever emits.
  for (expr in lst) {
    if (!.is_assignment(expr)) next
    lhs <- expr[[2]]
    if (!is.symbol(lhs)) next
    raw <- as.character(lhs)
    if (.is_rxmissingvars(raw)) next
    if (!raw %in% names(name_map)) name_map[raw] <- .norm(raw)
  }
  name_map <- .pin_names(name_map, state_pins)

  # Variables known to hold structural-model outputs (linCmt, ODE states).
  # Propagated forward; used in pass 2 to classify auxiliaries.
  aux_vars <- toupper(sigma_names)  # eps1, eps2, ...

  for (expr in lst) {
    # cmt() declarations from nonmem2rx -- skip silently
    if (is.call(expr) && identical(as.character(expr[[1]]), "cmt")) next

    # nonmem2rx bookkeeping. When a NONMEM model names a theta by its $THETA
    # label (`FLUX = KTP*A(1)` for `(0,0.2) ; KTP`), nonmem2rx does not bind the
    # symbol -- it emits `rxmissingvars1 <- t.KTP` to record that the theta was
    # referenced and left the reference itself dangling. The placeholder carries
    # no information the theta does not, and emitting it produced the meaningless
    # individual parameter `RXMISSINGVARS1 = TVKTP`. The dangling reference is
    # bound separately, by the [odes] carrier.
    if (.is_assignment(expr) && is.symbol(expr[[2]]) &&
        .is_rxmissingvars(as.character(expr[[2]]))) next

    if (.is_lincmt_tilde(expr)) {
      err_out     <- .parse_error_rhs(expr[[3]], name_map)
      error_model <- c(error_model,
                       list(list(dv = "DV", type = err_out$type,
                                 params = err_out$params)))
      warnings    <- c(warnings, err_out$warnings)
      if (!identical(structural$type, "ode"))
        structural <- list(type = "lincmt")
      next
    }

    if (.is_tilde(expr)) {
      err_out     <- .parse_error_rhs(expr[[3]], name_map)
      error_model <- c(error_model,
                       list(list(dv = "DV", type = err_out$type,
                                 params = err_out$params)))
      warnings    <- c(warnings, err_out$warnings)
      next
    }

    if (.is_assignment(expr)) {
      lhs_expr <- expr[[2]]

      # d/dt(STATE) <- rhs  or  d/dt(STATE) = rhs
      if (.is_ddt_lhs(lhs_expr)) {
        # The declaration has to be renamed with the references, so the d/dt
        # target is read through the rename map by hand -- .normalise_expr()
        # only rewrites symbols inside an expression.
        #
        # Through `state_decl`, NOT `name_map`: see the `states` parameter above.
        state_raw     <- .ddt_state(lhs_expr)
        state         <- if (state_raw %in% names(state_decl)) state_decl[[state_raw]]
                         else                                  state_raw
        # `ode_pins` on top of the map, and only here. Inside [odes] a bare name
        # that is both a state and a parameter is the state -- thetas and etas
        # are out of scope in that block -- while the same name in $PK means the
        # parameter, so the two scopes must not share one map. (`pins`, the
        # unambiguous states, is already folded in at entry.)
        rhs_expr_norm <- .normalise_expr(expr[[3]],
                                         .pin_names(name_map, state_ode_pins))
        rhs           <- paste(deparse(rhs_expr_norm, width.cutoff = 500L), collapse = " ")
        odes  <- c(odes, list(list(state = state, rhs = rhs, rhs_expr = rhs_expr_norm)))
        aux_vars <- c(aux_vars, toupper(state))  # ODE state vars are auxiliary
        if (!identical(structural$type, "ode"))
          structural <- list(type = "ode")
        next
      }

      # Skip non-symbol LHS (e.g. EFFECT(0) <- ..., f(ABS) <- ...)
      if (!is.symbol(lhs_expr)) next

      lhs_raw  <- as.character(lhs_expr)
      lhs_norm <- .norm(lhs_raw)
      rhs_expr <- expr[[3]]

      # rxlincmt1 <- linCmt()  -- nonmem2rx assignment form
      if (is.call(rhs_expr) && identical(rhs_expr[[1]], as.name("linCmt"))) {
        aux_vars <- c(aux_vars, lhs_norm)
        if (!identical(structural$type, "ode"))
          structural <- list(type = "lincmt")
        next
      }

      # Normalise the RHS BEFORE installing the alias: in `x <- f(x)` the RHS
      # refers to the previous binding of x, not the one being created. rxode2
      # sources can name a theta and the variable it defines identically
      # (`cl <- cl * exp(eta.cl)`), and installing the alias first rewrote that
      # theta reference into a self-reference the engine rejects. nonmem2rx
      # prefixes thetas (`t.CL`), which is why NONMEM never hit it.
      # A name with no mapping is left untouched by .normalise_expr(), so a
      # self-reference to a plain local (`k <- k * 2`, k not in iniDf) would
      # emit a bare lower-case `k` that is declared nowhere. Seed the alias for
      # that case only -- when lhs_raw IS a mapped name, the map already holds
      # the previous binding and must win.
      rhs_map <- name_map
      if (!lhs_raw %in% names(rhs_map)) rhs_map[lhs_raw] <- lhs_norm
      rhs_expr_norm <- .normalise_expr(rhs_expr, rhs_map)
      rhs_norm <- paste(deparse(rhs_expr_norm, width.cutoff = 500L), collapse = " ")
      # Update name_map so subsequent expressions see the alias -- unless the
      # name is a pinned state, whose emitted name was decided before the walk
      # and must survive it. `central <- 0` standing above `d/dt(central)` used
      # to rebind the state's key here, so the ODE referenced the individual
      # parameter instead of the compartment and the same model with that line
      # moved below it emitted a different equation.
      if (!lhs_raw %in% names(state_pins)) name_map[lhs_raw] <- lhs_norm

      # rhs_expr_norm is kept because name_map is time-varying: de-shadowing
      # rebinds a name mid-parse, so re-normalising the raw expression later
      # (pass 2b) would resolve it against the FINAL map and silently inline the
      # wrong variable into the ODEs.
      all_assigns <- c(all_assigns,
                       list(list(lhs = lhs_norm, rhs = rhs_norm,
                                 rhs_expr = rhs_expr,
                                 rhs_expr_norm = rhs_expr_norm)))
    }
  }

  # Pass 2: propagate aux_vars to fixpoint.
  # Any variable whose RHS contains an aux_var is itself auxiliary.
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (a in all_assigns) {
      if (a$lhs %in% aux_vars) next
      # The normalised form, not the raw one: a dotted local (`c.2`) uppercases
      # to `C.2` while aux_vars holds `C_2`, so the raw comparison misses it.
      syms <- toupper(.collect_symbols(a$rhs_expr_norm))
      if (any(syms %in% aux_vars)) {
        aux_vars <- c(aux_vars, a$lhs)
        changed  <- TRUE
      }
    }
  }

  # Pass 2b: inline aux-var definitions into ODE RHS strings.
  # $DES-internal intermediates (e.g. C2, EFF) are excluded from
  # [individual_parameters] but referenced in d/dt() expressions. Without
  # inlining, they appear as undefined names that ferx-core rejects at parse time.
  if (length(odes) > 0) {
    state_upper <- toupper(vapply(odes, function(o) o$state, ""))
    sigma_upper <- toupper(sigma_names)
    aux_map     <- list()
    for (a in all_assigns) {
      if (a$lhs %in% aux_vars &&
          !a$lhs %in% state_upper &&
          !a$lhs %in% sigma_upper)
        aux_map[[a$lhs]] <- a$rhs_expr_norm
    }
    if (length(aux_map) > 0) {
      odes <- lapply(odes, function(o) {
        inlined <- .inline_aux_vars(o$rhs_expr, aux_map)
        list(state = o$state,
             rhs   = paste(deparse(inlined, width.cutoff = 500L), collapse = " "))
      })
    } else {
      odes <- lapply(odes, function(o) list(state = o$state, rhs = o$rhs))
    }
  }

  # Pass 2c: collect RXM_* alias map for inline substitution.
  # nonmem2rx emits RXM_X = Y lines as internal IOV/eta copies. Collect them
  # all before Pass 3 (ordering in ui$lstExpr is not guaranteed) so that any
  # downstream indiv_param rhs that references RXM_X gets KAPPA_Y directly.
  rxm_map <- character()
  for (a in all_assigns) {
    if (grepl("^RXM_", a$lhs))
      rxm_map[[a$lhs]] <- a$rhs  # a$rhs is already normalised (e.g. "KAPPA_CL")
  }

  # Pass 3: classify each assignment into indiv_params or error_model.
  indiv_params <- list()
  for (a in all_assigns) {
    # Self-assignments arise from theta-alias resolution (tvcl <- t.TVCL -> TVCL <- TVCL).
    if (a$lhs == a$rhs) next

    # SCALE* vars are NONMEM-specific scaling intermediates.
    # RXINI* / RXF_* / RXM_* are nonmem2rx internal temporaries and IOV aliases.
    if (grepl("^SCALE\\d*$|^RXINI|^RXF_|^RXM_", a$lhs)) next

    if (a$lhs %in% aux_vars) {
      # Check if this is the error model assignment (RHS contains sigma vars).
      syms <- toupper(.collect_symbols(a$rhs_expr_norm))
      eps  <- intersect(syms, sigma_names)
      if (length(eps) > 0 && length(error_model) == 0) {
        # The normalised form, same as the detection above. Handing the raw
        # expression here made detection succeed and classification find no
        # sigma, emitting `DV ~ proportional()`.
        err  <- .classify_error_assignment(a$rhs_expr_norm, sigma_names)
        error_model <- c(error_model,
                         list(list(dv = "DV", type = err$type, params = err$params)))
      }
      next
    }

    # Inline RXM_* aliases so output references the real variable (e.g. KAPPA_CL).
    rhs_final <- a$rhs
    for (nm in names(rxm_map))
      rhs_final <- gsub(paste0("\\b", nm, "\\b"), rxm_map[[nm]], rhs_final, perl = TRUE)

    indiv_params <- c(indiv_params, list(list(lhs = a$lhs, rhs = rhs_final)))
  }

  list(
    indiv_params = indiv_params,
    odes         = odes,
    error_model  = error_model,
    structural   = structural,
    warnings     = warnings,
    unsupported  = unsupported
  )
}

# Classify an error expression (assignment RHS) into proportional / additive / combined.
# sigma_names: character vector of normalised sigma variable names (e.g. "EPS1").
.classify_error_assignment <- function(rhs_expr, sigma_names) {
  syms <- toupper(.collect_symbols(rhs_expr))
  eps  <- intersect(syms, sigma_names)

  if (length(eps) == 0)
    return(list(type = "proportional", params = character()))

  if (length(eps) >= 2)
    return(list(type = "combined", params = eps))

  # Single epsilon: multiplicative = proportional, additive = additive.
  fn <- tryCatch(as.character(rhs_expr[[1]]), error = function(e) "")
  type <- if (fn == "+") "additive" else "proportional"
  list(type = type, params = eps)
}

.parse_error_rhs <- function(rhs, name_map) {
  warn <- character()
  if (!is.call(rhs))
    return(list(type = "proportional", params = character(), warnings = warn))

  fn <- as.character(rhs[[1]])

  if (fn == "prop") {
    params <- .norm(.strip_prefix(as.character(rhs[[2]])))
    return(list(type = "proportional", params = params, warnings = warn))
  }
  if (fn == "add") {
    params <- .norm(.strip_prefix(as.character(rhs[[2]])))
    return(list(type = "additive", params = params, warnings = warn))
  }
  if (fn == "+") {
    lhs_fn <- tryCatch(as.character(rhs[[2]][[1]]), error = function(e) "")
    rhs_fn <- tryCatch(as.character(rhs[[3]][[1]]), error = function(e) "")
    if ((lhs_fn == "add" && rhs_fn == "prop") ||
        (lhs_fn == "prop" && rhs_fn == "add")) {
      add_node  <- if (lhs_fn == "add")  rhs[[2]] else rhs[[3]]
      prop_node <- if (lhs_fn == "prop") rhs[[2]] else rhs[[3]]
      params    <- c(.norm(.strip_prefix(as.character(prop_node[[2]]))),
                     .norm(.strip_prefix(as.character(add_node[[2]]))))
      return(list(type = "combined", params = params, warnings = warn))
    }
  }

  warn <- c(warn, "WARN  | complex $ERROR -- classified as proportional, verify")
  params <- tryCatch(.norm(.strip_prefix(as.character(rhs[[2]]))),
                     error = function(e) character())
  list(type = "proportional", params = params, warnings = warn)
}

# -- linCmt -> pk macro -------------------------------------------------------

.infer_pk_macro <- function(indiv_params) {
  lhs_lc   <- tolower(vapply(indiv_params, function(p) p$lhs, ""))
  lhs_uc   <- vapply(indiv_params, function(p) p$lhs, "")
  warn     <- character()
  unsp     <- character()

  # Detect model complexity by presence of q2/q3 (3-cpt) or q (2-cpt).
  has_ka  <- "ka"  %in% lhs_lc
  has_q   <- "q"   %in% lhs_lc
  has_q2  <- "q2"  %in% lhs_lc || "q3" %in% lhs_lc

  pk_call <- if (has_ka && has_q2) {
    "three_cpt_oral"
  } else if (has_ka && has_q) {
    "two_cpt_oral"
  } else if (has_ka) {
    "one_cpt_oral"
  } else if (has_q2) {
    "three_cpt_iv"
  } else if (has_q) {
    "two_cpt_iv"
  } else {
    "one_cpt_iv"
  }

  # For each ferx argument key, define an ordered list of candidate lhs_lc names.
  # nonmem2rx ADVAN4/TRANS4 (2cpt oral) names volumes v2/v3 instead of v1/v2.
  # nonmem2rx ADVAN3 (2cpt IV) names volumes v1/v2 directly.
  arg_aliases <- switch(pk_call,
    one_cpt_oral = list(cl = "cl", v = c("v", "v1", "v2"), ka = "ka"),
    one_cpt_iv   = list(cl = "cl", v = c("v", "v1", "v2")),
    # two_cpt_oral: try v1 first, then plain v (nlmixr2 alias), then v2 (NONMEM ADVAN4 TRANS4)
    # same pattern for peripheral: v2 -> v3 (NONMEM ADVAN4 TRANS4)
    two_cpt_oral = list(cl = "cl", v1 = c("v1", "v", "v2"), q = "q",
                        v2 = c("v2", "v3"), ka = "ka"),
    two_cpt_iv   = list(cl = "cl", v1 = c("v1", "v"), q = "q",
                        v2 = c("v2", "v3")),
    # 3-cpt: ferx slot Q (first inter-compartmental clearance) accepts arg name
    # `q` or `q2`; slot Q3 (second clearance) accepts ONLY `q3` (see ferx-core
    # PkParams::name_to_index). NONMEM ADVAN11 names them Q2 (first) / Q3
    # (second), so emit `q2=Q2` and `q3=Q3` -- emitting both clearances, not
    # dropping Q3.
    three_cpt_oral = list(cl = "cl", v1 = c("v1", "v"), q2 = c("q2", "q"),
                          v2 = "v2", q3 = "q3", v3 = c("v3", "v4"), ka = "ka"),
    three_cpt_iv   = list(cl = "cl", v1 = c("v1", "v"), q2 = c("q2", "q"),
                          v2 = "v2", q3 = "q3", v3 = c("v3", "v4")),
    list()
  )

  # A pk_call with no argument mapping (switch fell through to the default
  # empty list) is an unrecognised structure. Degrade gracefully: emit an
  # ERROR and omit the structural model rather than emitting an argument-less
  # macro call that ferx-core would reject at parse time. This keeps the
  # downstream is.na(pk_call) guard in rxui_to_ir() a live safety net.
  if (length(arg_aliases) == 0) {
    unsp <- c(unsp, paste0(pk_call, " (unsupported structure -- no argument mapping)"))
    warn <- c(warn, paste0("ERROR | ", pk_call,
                           " has no argument mapping -- structural model omitted"))
    return(list(pk_call = NA_character_, pk_args = list(),
                warnings = warn, unsupported = unsp))
  }

  # Optional args
  if ("f"    %in% lhs_lc) arg_aliases[["f"]]    <- "f"
  if ("alag" %in% lhs_lc || "lagtime" %in% lhs_lc || "tlag" %in% lhs_lc)
    arg_aliases[["alag"]] <- c("alag", "lagtime", "tlag")

  # Greedy matching -- each lhs_lc index used at most once.
  used_idx <- integer()
  pk_args  <- list()
  for (key in names(arg_aliases)) {
    candidates <- arg_aliases[[key]]
    for (cand in candidates) {
      idxs <- setdiff(which(lhs_lc == cand), used_idx)
      if (length(idxs) > 0) {
        used_idx       <- c(used_idx, idxs[1])
        pk_args[[key]] <- lhs_uc[idxs[1]]
        break
      }
    }
  }

  list(pk_call = pk_call, pk_args = pk_args, warnings = warn, unsupported = unsp)
}
