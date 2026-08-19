# Shared string and expression helpers

# Parse the NONMEM $DATA record for the dataset path, resolved relative to the
# control stream's own directory (NONMEM's convention). Returns NA_character_
# when there is no $DATA record or the file it names does not exist.
#
# Used to validate emitted .ferx WITH data. That distinction matters: model-only
# validation silently accepts any undefined name that could be a covariate --
# `y = ... UNDEFINED_NAME` reports ok = TRUE -- and only turns into
# E_MISSING_COVARIATE once a dataset is supplied.
.extract_nm_data_path <- function(ctl_file) {
  if (!is.character(ctl_file) || length(ctl_file) != 1L || !file.exists(ctl_file))
    return(NA_character_)
  lines <- tryCatch(readLines(ctl_file, warn = FALSE), error = function(e) character())
  idx   <- grep("^\\s*\\$DATA\\b", lines, ignore.case = TRUE)[1L]
  if (is.na(idx)) return(NA_character_)

  rest <- sub("^\\s*\\$DATA\\s*", "", lines[idx], ignore.case = TRUE)
  rest <- sub(";.*$", "", rest)          # strip a trailing NONMEM comment
  # A quoted path may contain spaces; an unquoted one ends at the first space,
  # which is where options such as IGNORE=@ begin.
  path <- if (grepl('^\\s*"', rest)) {
    sub('^\\s*"([^"]*)".*$', "\\1", rest)
  } else if (grepl("^\\s*'", rest)) {
    sub("^\\s*'([^']*)'.*$", "\\1", rest)
  } else {
    sub("^\\s*(\\S+).*$", "\\1", rest)
  }
  if (!nzchar(path)) return(NA_character_)

  resolved <- if (.is_abs_path(path)) path else file.path(dirname(ctl_file), path)
  # file.exists() is TRUE for directories, and `$DATA dat` beside a `dat/`
  # directory is an ordinary NONMEM layout -- pkpd_ir.mod declares exactly that.
  if (file.exists(resolved) && !dir.exists(resolved))
    normalizePath(resolved, mustWork = FALSE)
  else NA_character_
}

.is_abs_path <- function(p) grepl("^(/|~|[A-Za-z]:[\\\\/])", p)

# Parse NONMEM $PK block for Sn = varname compartment scaling assignments.
# Returns a named list mapping compartment number (as character) to variable name.
.extract_nm_scaling <- function(ctl_file) {
  lines <- tryCatch(readLines(ctl_file, warn = FALSE), error = function(e) character())
  if (length(lines) == 0L) return(list())
  # NONMEM block records (`$PK`, `$DES`, ...) may carry leading whitespace in
  # some tool outputs; tolerate it so an indented `$PK` is not silently missed.
  block_starts <- grep("^\\s*\\$", lines)
  pk_idx       <- grep("^\\s*\\$PK\\b", lines, ignore.case = TRUE)[1L]
  if (is.na(pk_idx)) return(list())
  next_after   <- block_starts[block_starts > pk_idx][1L]
  pk_lines     <- if (is.na(next_after)) lines[pk_idx:length(lines)]
                  else                   lines[pk_idx:(next_after - 1L)]
  result <- list()
  for (line in pk_lines) {
    m <- regmatches(line,
                    regexpr("^\\s*S(\\d+)\\s*=\\s*(\\w+)", line, perl = TRUE))
    if (length(m) > 0L && nchar(m) > 0L) {
      cmt_n        <- sub("^\\s*S(\\d+)\\s*=.*",      "\\1", m, perl = TRUE)
      var          <- sub("^\\s*S\\d+\\s*=\\s*(\\w+).*", "\\1", m, perl = TRUE)
      result[[cmt_n]] <- var
    }
  }
  result
}

# The compartment $MODEL marks DEFOBS, as a 1-based index into the COMP list.
#
# Read from the raw control stream for the same reason .extract_nm_scaling() is:
# nonmem2rx does not surface it. `ui$central` is NULL for both nonmem2rx and
# rxode2, so obs_cmt fell through to tail(state_names, 1) -- a positional guess
# that is right only when the observed compartment happens to be declared last.
# It is not a cosmetic guess: the same index selects the $PK scaling variable, so
# a model with `COMP=(CENT, DEFDOSE, DEFOBS)` before `COMP=(PERIPH)` and
# `S1 = V` got no [scaling] block at all -- the S2=V silent-divergence class
# CLAUDE.md warns about -- or, with `S2 = V2`, scaled the observation by the
# peripheral volume and validated clean.
#
# Returns NULL when $MODEL is absent or names no DEFOBS. Deliberately does NOT
# fall back to NONMEM's own default (compartment 1 when DEFOBS is omitted): the
# caller already has a documented guess with a warning attached, and replacing a
# loud guess with a quiet assumption about NONMEM semantics would trade a
# visible risk for an invisible one.
.extract_nm_defobs <- function(ctl_file) {
  lines <- tryCatch(readLines(ctl_file, warn = FALSE), error = function(e) character())
  if (length(lines) == 0L) return(NULL)
  block_starts <- grep("^\\s*\\$", lines)
  model_idx    <- grep("^\\s*\\$MODEL\\b", lines, ignore.case = TRUE)[1L]
  if (is.na(model_idx)) return(NULL)
  next_after   <- block_starts[block_starts > model_idx][1L]
  body <- if (is.na(next_after)) lines[model_idx:length(lines)]
          else                   lines[model_idx:(next_after - 1L)]
  # Strip comments before parsing: `COMP=(CENT) ; DEFOBS is below` must not
  # register a compartment attribute that is only mentioned in prose.
  body <- sub(";.*$", "", body)
  txt  <- paste(body, collapse = " ")
  # Each COMP declaration, in source order. Both `COMP=(...)` and `COMP (...)`
  # are legal NONMEM.
  m <- gregexpr("COMP\\s*=?\\s*\\(([^)]*)\\)", txt, ignore.case = TRUE, perl = TRUE)
  decls <- regmatches(txt, m)[[1]]
  if (length(decls) == 0L) return(NULL)
  for (i in seq_along(decls)) {
    inner <- sub("^COMP\\s*=?\\s*\\((.*)\\)$", "\\1", decls[i],
                 ignore.case = TRUE, perl = TRUE)
    toks  <- trimws(strsplit(inner, "[,[:space:]]+")[[1]])
    toks  <- toks[nzchar(toks)]
    if (length(toks) < 2L) next
    # NONMEM permits abbreviating an attribute to any unambiguous prefix, and
    # DEFOBSERVATION is the full spelling, so match a prefix rather than the
    # exact token. DEFDOSE shares "DEF", hence the four-character anchor.
    if (any(grepl("^DEFOBS", toks[-1], ignore.case = TRUE)))
      return(list(index = i, name = toks[1], n_comp = length(decls)))
  }
  NULL
}
