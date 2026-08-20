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
  pk_lines <- .nm_block_lines(ctl_file, "PK")
  if (length(pk_lines) == 0L) return(list())
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

# The body lines of one NONMEM block record, comments stripped.
#
# One owner, because two readers of the same control stream must agree about
# where a record ends. `.extract_nm_scaling()` decides which `S<n>` variable
# scales the observation and `.extract_nm_defobs()` decides which compartment
# number that `n` refers to, so a disagreement about block boundaries makes them
# number compartments differently -- and the result is an observation scaled by
# some other compartment's volume, which ferx validates clean.
#
# Record names may be abbreviated to any unambiguous prefix in NM-TRAN, hence
# the prefix match rather than the full word.
.nm_block_lines <- function(ctl_file, record) {
  lines <- tryCatch(readLines(ctl_file, warn = FALSE), error = function(e) character())
  if (length(lines) == 0L) return(character())
  # Block records may carry leading whitespace in some tool outputs; tolerate it
  # so an indented `$PK` is not silently missed.
  block_starts <- grep("^\\s*\\$", lines)
  idx <- grep(paste0("^\\s*\\$", record), lines, ignore.case = TRUE)[1L]
  if (is.na(idx)) return(character())
  nxt <- block_starts[block_starts > idx][1L]
  body <- if (is.na(nxt)) lines[idx:length(lines)] else lines[idx:(nxt - 1L)]
  # Strip comments before parsing: `COMP=(CENT) ; DEFOBS is below` must not
  # register an attribute that is only mentioned in prose.
  sub(";.*$", "", body)
}

# The compartment $MODEL marks DEFOBS, as a 1-based index into the COMP list.
#
# Read from the raw control stream for the same reason .extract_nm_scaling() is:
# nonmem2rx does not surface it. `ui$central` is NULL for both nonmem2rx and
# rxode2, so obs_cmt falls back to a positional guess without it. That index is
# also the NONMEM compartment number the `S<n>` scaling lookup uses, so getting
# it wrong loses the scaling silently.
#
# Returns NULL when $MODEL is absent or declares no compartments. When $MODEL
# exists but names no DEFOBS, `index` and `name` come back NA and only `comps`
# is populated -- callers already require a non-NA index, and the COMP list is
# worth returning on its own: it is the only way to check a NONMEM compartment
# NUMBER against the d/dt order it is about to index into. Deliberately does NOT
# fall back to NONMEM's own default (compartment 1 when DEFOBS is omitted): the
# caller has better evidence than that -- the compartment the DV expression
# names -- and a quiet assumption about NONMEM semantics would outrank it.
.extract_nm_defobs <- function(ctl_file) {
  # `$MOD` is a legal abbreviation of `$MODEL`.
  body <- .nm_block_lines(ctl_file, "MOD")
  if (length(body) == 0L) return(NULL)
  txt  <- paste(body, collapse = " ")
  # Every COMP declaration, in source order, in all the forms NM-TRAN accepts:
  # `COMP=(...)`, `COMP (...)`, `COMPARTMENT=(...)`, and the bare `COMP=NAME` /
  # `COMP NAME` forms, which carry no attributes but DO occupy a compartment
  # number -- miscounting them shifts every later ordinal.
  #
  # `\bCOMP` and not `COMP`: `NCOMP=4` is a different keyword on the same
  # record and must not be counted (there is no word boundary inside NCOMP, so
  # the anchor excludes it).
  re <- "\\bCOMP(ARTMENT)?\\s*=?\\s*(\\(([^)]*)\\)|[A-Za-z_][A-Za-z0-9_]*)"
  m  <- gregexpr(re, txt, ignore.case = TRUE, perl = TRUE)
  decls <- regmatches(txt, m)[[1]]
  if (length(decls) == 0L) return(NULL)
  comps  <- character(length(decls))
  defobs <- NA_integer_
  for (i in seq_along(decls)) {
    inner <- sub("^\\bCOMP(ARTMENT)?\\s*=?\\s*", "", decls[i],
                 ignore.case = TRUE, perl = TRUE)
    inner <- sub("^\\((.*)\\)$", "\\1", inner)
    toks  <- trimws(strsplit(inner, "[,[:space:]]+")[[1]])
    toks  <- toks[nzchar(toks)]
    if (length(toks) == 0L) next
    # Recorded for every compartment, attributes or not: a bare `COMP=NAME`
    # occupies an ordinal and naming it is the whole point of the list.
    comps[i] <- toks[1]
    if (length(toks) < 2L) next
    # NM-TRAN permits abbreviating an attribute to any unambiguous prefix, and
    # the full spelling is DEFOBSERVATION. The only other DEF* attribute is
    # DEFDOSE, which diverges at the fourth character, so DEFO is the shortest
    # unambiguous prefix -- `^DEFOBS` rejected the legal `DEFO` and `DEFOB` and
    # silently reverted the caller to its guess.
    if (is.na(defobs) && any(grepl("^DEFO", toks[-1], ignore.case = TRUE)))
      defobs <- i
  }
  list(index  = defobs,
       name   = if (is.na(defobs)) NA_character_ else comps[defobs],
       n_comp = length(decls),
       comps  = comps)
}
