test_that("new_ferx_ir constructs with defaults", {
  ir <- new_ferx_ir()
  expect_s3_class(ir, "ferx_ir")
  expect_true(is.na(ir$source_format))
  expect_true(is.na(ir$source_file))
  expect_equal(ir$thetas, list())
  expect_equal(ir$omegas, list())
  expect_equal(ir$warnings, character())
})

test_that("new_ferx_ir stores supplied values", {
  ir <- new_ferx_ir(
    source_format = "nonmem",
    source_file   = "run001.ctl",
    thetas        = list(list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10)),
    omegas        = list(list(type = "diagonal", names = "ETA_CL", values = 0.07)),
    warnings      = c("WARN | something")
  )
  expect_equal(ir$source_format, "nonmem")
  expect_equal(ir$thetas[[1]]$name, "TVCL")
  expect_equal(ir$omegas[[1]]$type, "diagonal")
  expect_length(ir$warnings, 1L)
})

test_that("validate_ferx_ir accepts empty IR", {
  ir <- new_ferx_ir()
  expect_invisible(validate_ferx_ir(ir))
})

test_that("validate_ferx_ir accepts IR with valid structural type", {
  ir <- new_ferx_ir(structural = list(type = "pk_macro", pk_call = "one_cpt_oral", pk_args = list()))
  expect_invisible(validate_ferx_ir(ir))
})

test_that("validate_ferx_ir rejects structural with missing type", {
  ir <- new_ferx_ir(structural = list(pk_call = "one_cpt_oral"))
  expect_error(validate_ferx_ir(ir), "structural\\$type is missing")
})

test_that("validate_ferx_ir rejects unknown structural type", {
  ir <- new_ferx_ir(structural = list(type = "compartmental"))
  expect_error(validate_ferx_ir(ir), "pk_macro")
})

test_that("validate_ferx_ir rejects odes without structural type ode", {
  ir <- new_ferx_ir(
    structural = list(type = "pk_macro", pk_call = "one_cpt_oral", pk_args = list()),
    odes       = list(list(state = "depot", rhs = "-KA * depot"))
  )
  expect_error(validate_ferx_ir(ir), "structural\\$type is not")
})

test_that("validate_ferx_ir rejects non-ferx_ir input", {
  expect_error(validate_ferx_ir(list()), "ferx_ir")
})

test_that("validate_ferx_ir accepts ode structural with states and obs_cmt", {
  ir <- new_ferx_ir(
    structural = list(type = "ode", states = c("depot", "central"), obs_cmt = "central"),
    odes       = list(list(state = "depot", rhs = "-KA * depot"),
                      list(state = "central", rhs = "KA * depot - CL/V * central"))
  )
  expect_invisible(validate_ferx_ir(ir))
})

test_that("validate_ferx_ir rejects ode structural missing states", {
  ir <- new_ferx_ir(
    structural = list(type = "ode", obs_cmt = "central"),
    odes       = list(list(state = "central", rhs = "-CL/V * central"))
  )
  expect_error(validate_ferx_ir(ir), "structural\\$states")
})

test_that("validate_ferx_ir rejects ode structural missing obs_cmt", {
  ir <- new_ferx_ir(
    structural = list(type = "ode", states = c("depot", "central")),
    odes       = list(list(state = "depot", rhs = "-KA * depot"),
                      list(state = "central", rhs = "KA * depot"))
  )
  expect_error(validate_ferx_ir(ir), "structural\\$obs_cmt")
})

test_that("print.ferx_ir runs without error", {
  ir <- new_ferx_ir(
    source_format = "nonmem",
    source_file   = "run001.ctl",
    thetas        = list(list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10)),
    omegas        = list(list(type = "diagonal", names = "ETA_CL", values = 0.07)),
    structural    = list(type = "pk_macro", pk_call = "one_cpt_oral", pk_args = list()),
    warnings      = c("WARN | check theta bounds"),
    unsupported   = c("MIXTURE model")
  )
  expect_no_error(print(ir))
})

test_that("print.ferx_ir handles empty IR without error", {
  expect_no_error(print(new_ferx_ir()))
})

# -- emitted-name legality ----------------------------------------------------
#
# The structural half of the identifier guarantee. Legality used to be enforced
# only where each name was minted plus one corpus test over whichever models are
# bundled, so a name reaching the file through an unchecked channel was emitted
# verbatim into an unparseable file. These assert the guard fires per channel,
# so adding a channel without wiring it into .ir_declared_names() shows up here.

.pk_ir <- function(...) new_ferx_ir(
  structural = list(type = "pk_macro", pk_call = "one_cpt_oral",
                    pk_args = list(cl = "CL", v = "V")), ...)

test_that("a legal IR still validates and emits", {
  expect_silent(validate_ferx_ir(.pk_ir()))
  expect_type(emit_ferx(.pk_ir()), "character")
})

test_that("an illegal declared name is rejected in every channel", {
  channels <- list(
    theta       = .pk_ir(thetas = list(list(name = "1BAD", init = 1,
                                            lower = 0, upper = 10))),
    omega       = .pk_ir(omegas = list(list(type = "diagonal",
                                            names = "c.RTOT", values = 0.1))),
    kappa       = .pk_ir(kappas = list(list(name = "KAPPA-1", value = 0.1))),
    sigma       = .pk_ir(sigmas = list(list(name = "EPS.1", value = 0.1,
                                            scale = "sd"))),
    indiv_param = .pk_ir(indiv_params = list(list(lhs = "A B", rhs = "1"))),
    error_dv    = .pk_ir(error_model = list(list(dv = "D V", type = "proportional",
                                                 params = "PROP"))),
    error_param = .pk_ir(error_model = list(list(dv = "DV", type = "proportional",
                                                 params = "PROP.ERR"))),
    diffusion   = .pk_ir(diffusion = list(list(state = "9X", value = 0.1))),
    pk_arg      = new_ferx_ir(structural = list(type = "pk_macro",
                                                pk_call = "one_cpt_oral",
                                                pk_args = list(cl = "C.L"))),
    state       = new_ferx_ir(structural = list(type = "ode", states = "9CENT",
                                                obs_cmt = "9CENT"),
                              odes = list(list(state = "9CENT", rhs = "-K"))))
  for (ch in names(channels))
    expect_error(validate_ferx_ir(channels[[ch]]),
                 "not .*legal ferx identifier", info = ch)
})

test_that("the rejection names the channel, not just the name", {
  expect_error(validate_ferx_ir(.pk_ir(omegas = list(list(type = "diagonal",
                                                          names = "c.RTOT",
                                                          values = 0.1)))),
               "omega")
})

test_that("obs_cmt must be one of the declared states", {
  # Emitted, this produced only INFO-level warnings and an empty $unsupported --
  # a translation that looked clean and a file the engine refused with
  # "Observable compartment 'CENT' not in states".
  ir <- new_ferx_ir(structural = list(type = "ode", states = "c_CENT",
                                      obs_cmt = "CENT"),
                    odes = list(list(state = "c_CENT", rhs = "-K")))
  expect_error(validate_ferx_ir(ir), "obs_cmt.*not one of")
})

test_that("state_renames must be a named character vector", {
  # Unnamed, emit_ferx() rendered `# renamed: state  -> c_RTOT` -- provenance
  # with the source half blank, which is the half a reader holding only the
  # .ferx file cannot reconstruct.
  expect_error(validate_ferx_ir(.pk_ir(state_renames = c("c_RTOT", "A_B"))),
               "must be NAMED")
  expect_error(validate_ferx_ir(.pk_ir(state_renames = c(RTOT = "c_RTOT",
                                                         "A_B"))),
               "must be NAMED")
  expect_silent(validate_ferx_ir(.pk_ir(state_renames = c(RTOT = "c_RTOT"))))
})

test_that("the legality check reads declarations, not expression text", {
  # Covariates live in RHS expressions and must keep the data column's exact
  # spelling -- ferx matches them case-sensitively, so sanitising one turns a
  # working reference into E_MISSING_COVARIATE at fit time. They are reported as
  # untranslatable by the translator, not rewritten, so this guard must not fire
  # on them or it would abort the very models that need the diagnostic.
  ir <- .pk_ir(indiv_params = list(list(lhs = "CL", rhs = "TVCL * WT.KG")))
  expect_silent(validate_ferx_ir(ir))
})

test_that("an uncommented error_suggestion line is rejected", {
  # The suggestion is a reading of an expression we could NOT translate. One
  # uncommented line and the engine parses the guess as the translation --
  # the failure the field exists to prevent.
  ok <- new_ferx_ir(error_suggestion = c("# could not translate",
                                         "# [error_model]",
                                         "#   DV ~ additive(EPS1)"))
  expect_silent(validate_ferx_ir(ok))
  # Leading whitespace before the # is still a comment.
  expect_silent(validate_ferx_ir(new_ferx_ir(error_suggestion = "  # indented")))

  bad <- new_ferx_ir(error_suggestion = c("# could not translate",
                                          "[error_model]",
                                          "  DV ~ additive(EPS1)"))
  expect_error(validate_ferx_ir(bad), "only comment lines")
})

test_that("a readout and obs_scale cannot both be carried", {
  # Measured on ferx 0.3.0: the file validates clean and every prediction moves
  # by up to 11.37 units, because ferx applies obs_scale on TOP of the readout
  # expression rather than instead of it. No diagnostic from the engine, which
  # is why the guard belongs here. This is defect 15.
  ode <- function(scaling) new_ferx_ir(
    structural = list(type = "ode", states = "CENT", obs_cmt = "CENT"),
    odes = list(list(kind = "ddt", state = "CENT", rhs = "-CENT")),
    scaling = scaling)
  expect_error(validate_ferx_ir(ode(list(y = "CENT", obs_scale = "V"))), "double-scal")
  expect_error(validate_ferx_ir(ode(list(per_cmt = list(list(cmt = 1L, expr = "CENT")),
                                        obs_scale = "V"))), "double-scal")
  expect_error(validate_ferx_ir(
    ode(list(y = "CENT", per_cmt = list(list(cmt = 1L, expr = "CENT"))))),
    "both a plain `y` and a per-CMT")
  expect_silent(validate_ferx_ir(ode(list(y = "CENT"))))
  expect_silent(validate_ferx_ir(ode(list(obs_scale = "V"))))
})

test_that("obs_cmt is required for an ODE model only when no readout replaces it", {
  no_obs <- function(scaling) new_ferx_ir(
    structural = list(type = "ode", states = "CENT"),
    odes = list(list(kind = "ddt", state = "CENT", rhs = "-CENT")),
    scaling = scaling)
  expect_error(validate_ferx_ir(no_obs(list())), "obs_cmt must be a character scalar")
  expect_silent(validate_ferx_ir(no_obs(list(y = "CENT"))))
  expect_silent(validate_ferx_ir(no_obs(list(per_cmt = list(list(cmt = 1L, expr = "CENT"))))))
})

test_that("an error_model must be one form, and a selected chain must end in an else", {
  em <- function(...) validate_ferx_ir(new_ferx_ir(error_model = list(...)))
  expect_error(em(list(dv = "DV", type = "additive", params = "E1", cmt = 1L),
                  list(dv = "DV", type = "additive", params = "E2", cond = "F == 2")),
               "mixes per-CMT and covariate-selected")
  expect_error(em(list(dv = "DV", type = "additive", params = "E1", cmt = 1L),
                  list(dv = "DV", type = "additive", params = "E2")),
               "mixes per-CMT entries with unkeyed")
  # ferx requires the terminating bare `else` "so every observation maps to an
  # error model". The emitter derives it from a NULL cond, so a chain without one
  # would emit `else if (...) { ... }` and stop -- rejected by the engine.
  expect_error(em(list(dv = "DV", type = "additive", params = "E1", cond = "F == 1"),
                  list(dv = "DV", type = "additive", params = "E2", cond = "F == 2")),
               "must end with an unconditional entry")
  expect_error(em(list(dv = "DV", type = "additive", params = "E1"),
                  list(dv = "DV", type = "additive", params = "E2", cond = "F == 2")),
               "only the LAST covariate-selected")
  expect_silent(validate_ferx_ir(
    em(list(dv = "DV", type = "additive", params = "E1", cond = "F == 1"),
       list(dv = "DV", type = "additive", params = "E2"))))
})
