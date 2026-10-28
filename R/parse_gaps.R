# Control streams nonmem2rx cannot parse, classified by what the user wrote.
#
# nonmem2rx raises for a handful of ordinary NONMEM constructs, in rxode2's
# vocabulary and without saying the limit is upstream. Left alone that reaches the
# user as "Failed to parse ... as nonmem" with nothing to act on, and is not a
# gap: no $unsupported entry, so a model with one cannot appear in the sweep's
# `model -> gap` table. A message matching nothing here is NOT absorbed -- it keeps
# raising, so a new upstream break stays visible.

.nm_parse_gaps <- list(
  list(id      = "dae",
       pattern = "Differential Algebra Equations are not supported",
       feature = "ADVAN9 / ADVAN15 (differential-algebraic equations, $AES)",
       ferx    = paste0("ferx has no algebraic or implicit-equation facility ",
                        "either, so this needs a ferx format feature.")),
  list(id      = "mm",
       pattern = "Mi(chelis|chaelis) Ment[eo]n model translation not supported",
       feature = "ADVAN10 (Michaelis-Menten elimination)",
       ferx    = paste0("ferx can express it as an [odes] model (ferx-r bundles ",
                        "mm_oral.ferx); the limit is in nonmem2rx and rxode2, ",
                        "whose solved-system form has no Michaelis-Menten variant.")),
  list(id      = "solved_cmt",
       pattern = "can only request depot and central compartments for solved systems",
       feature = "a solved ADVAN model that reads a compartment amount other than depot or central",
       ferx    = paste0("ferx can express it as an [odes] model, where every state ",
                        "is addressable; the limit is in nonmem2rx and rxode2, whose ",
                        "solved-system form exposes only depot and central."))
)

# NULL when `msg` matches no known gap, otherwise list(id, feature, ferx).
.classify_nm_parse_error <- function(msg) {
  msg <- paste(msg, collapse = " ")
  for (g in .nm_parse_gaps)
    if (grepl(g$pattern, msg, ignore.case = TRUE)) return(g)
  NULL
}

# The ADVAN the control stream names, for the message. NA when unreadable.
.nm_advan_label <- function(ctl_file) {
  if (!is.character(ctl_file) || length(ctl_file) != 1L || !file.exists(ctl_file))
    return(NA_character_)
  body <- paste(.nm_block_lines(ctl_file, "SUB"), collapse = " ")
  m <- regmatches(body, regexpr("ADVAN[0-9]+", body, ignore.case = TRUE))
  if (length(m) == 0L) NA_character_ else toupper(m)
}

# The result handed back for a recognised gap: no model text (a comment-only
# file would read as a model), an ERROR warning and an $unsupported entry.
.nm_gap_result <- function(gap, msg, src_file, advan = NA_character_) {
  feature <- gap$feature
  if (!is.na(advan) && !grepl(advan, feature, fixed = TRUE))
    feature <- paste0(feature, " [", advan, "]")
  warn <- paste0(
    "ERROR | ", feature, " cannot be translated: nonmem2rx, the parser this ",
    "package reads NONMEM through, rejects it. ", gap$ferx,
    " nonmem2rx said: ", gsub("\\s+", " ", paste(msg, collapse = " ")))
  ir <- new_ferx_ir(source_format = "nonmem", source_file = src_file,
                    warnings = warn,
                    unsupported = paste0("unparseable by nonmem2rx: ", feature))
  new_ferx_translate_result(
    paste0("# No .ferx was produced: ", feature, " could not be parsed.\n"), ir)
}
