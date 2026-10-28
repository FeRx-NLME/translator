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

# -- identifiers rxode2 rejects ------------------------------------------------

# NONMEM is case-insensitive, and nonmem2rx lower-cases variables on the way in,
# so `DF = CL/V` reaches rxode2 as `df <- cl/v`, which its parser rejects (it is
# a distribution function there). The name is immaterial to the model, so it is
# renamed in the stream nonmem2rx reads rather than reported as a gap.
#
# The rejected set is asked of rxode2 itself rather than listed here: each name a
# record assigns is parsed as `NAME <- k*2` in a throwaway model, which takes
# milliseconds because it is a parse and not a compile.

.nm_code_record <- function(rec)
  grepl("^(PK|PRE[A-Z]*|ERR[A-Z]*|DES[A-Z]*|AES[A-Z]*)$", rec, ignore.case = TRUE)

# Per line: is it inside a record that holds code? Header lines count, since
# `$PK CL = 1` is legal.
.nm_code_lines <- function(lines) {
  in_code <- FALSE
  vapply(lines, function(l) {
    m <- regmatches(l, regexec("^\\s*\\$([A-Za-z]+)", l))[[1]]
    if (length(m) == 2L) in_code <<- .nm_code_record(m[2])
    in_code
  }, NA, USE.NAMES = FALSE)
}

# Names assigned in code records (upper-cased, unique), comments stripped.
.nm_assigned_names <- function(lines) {
  code <- sub(";.*$", "", lines[.nm_code_lines(lines)])
  m <- regmatches(code, gregexpr(
    "(?<![[:alnum:]_.$])[A-Za-z_][A-Za-z0-9_]*(?=\\s*=(?!=))", code, perl = TRUE))
  unique(toupper(unlist(m)))
}

# Does rxode2 accept `name` as an ordinary variable? Asked of the lower-cased
# name, because that is what nonmem2rx hands it (`DF` passes, `df` does not).
# NA when rxode2 cannot be asked, which callers must treat as "no rename".
.rx_accepts_name <- function(name) {
  if (!requireNamespace("rxode2", quietly = TRUE)) return(NA)
  name <- tolower(name)
  txt <- paste0("d/dt(rxprobe_s) <- -rxprobe_k*rxprobe_s\n",
                name, " <- rxprobe_k*2\n",
                "d/dt(rxprobe_t) <- ", name, "*rxprobe_s")
  tryCatch({
    utils::capture.output(suppressMessages(suppressWarnings(
      rxode2::rxModelVars(txt))))
    TRUE
  }, error = function(e) FALSE)
}

# Named character vector OLD -> NEW (upper-case) for every assigned name rxode2
# rejects. `T` and `TIME` are NONMEM's own names for time and are handled by
# nonmem2rx, so they are never renamed.
.nm_reserved_renames <- function(ctl_file) {
  if (!is.character(ctl_file) || length(ctl_file) != 1L || !file.exists(ctl_file))
    return(character())
  lines <- tryCatch(readLines(ctl_file, warn = FALSE), error = function(e) character())
  nms   <- setdiff(.nm_assigned_names(lines), c("T", "TIME"))
  if (length(nms) == 0L) return(character())
  bad <- nms[vapply(nms, function(n) identical(.rx_accepts_name(n), FALSE), NA)]
  if (length(bad) == 0L) return(character())
  taken <- unique(toupper(unlist(regmatches(
    lines, gregexpr("[A-Za-z_][A-Za-z0-9_]*", sub(";.*$", "", lines))))))
  new <- character()
  for (b in bad) {
    cand <- paste0(b, "_VAR")
    while (cand %in% c(taken, new)) cand <- paste0(cand, "X")
    new <- c(new, cand)
  }
  stats::setNames(new, bad)
}

# Rename whole-word identifiers in the code records only, outside comments;
# every other byte of the stream is left as it was.
.nm_rename_identifiers <- function(lines, renames) {
  if (length(renames) == 0L) return(lines)
  code <- .nm_code_lines(lines)
  for (i in which(code)) {
    parts <- strsplit(lines[i], ";", fixed = TRUE)[[1]]
    if (length(parts) == 0L) next
    body <- parts[1]
    for (old in names(renames))
      body <- gsub(paste0("(?<![[:alnum:]_.$])", old, "(?![[:alnum:]_])"),
                   renames[[old]], body, perl = TRUE, ignore.case = TRUE)
    rest <- if (grepl(";", lines[i], fixed = TRUE))
      substring(lines[i], nchar(parts[1]) + 1L) else ""
    lines[i] <- paste0(body, rest)
  }
  lines
}

# A temporary copy of `ctl_file` with `renames` applied, under the same file name.
.nm_renamed_copy <- function(ctl_file, renames) {
  d <- tempfile("ferxtr-rename-")
  dir.create(d)
  out <- file.path(d, basename(ctl_file))
  writeLines(.nm_rename_identifiers(readLines(ctl_file, warn = FALSE), renames), out)
  out
}

# A temporary copy of `ctl_file` with `renames` applied, under the same file name.
.nm_renamed_copy <- function(ctl_file, renames) {
  d <- tempfile("ferxtr-rename-")
  dir.create(d)
  out <- file.path(d, basename(ctl_file))
  writeLines(.nm_rename_identifiers(readLines(ctl_file, warn = FALSE), renames), out)
  out
}
