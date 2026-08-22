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

# The compartment $MODEL marks DEFOBS, as a 1-based index into the COMP list --
# and, from the same parse, the compartment it marks DEFDOSE. The name says
# defobs and the return value says more than that: `comps`, `n_comp` and
# `defdose` come back too, because all four are read off the one COMP list and a
# second parser for the same record is how two readers start disagreeing about
# where it ends. See .nm_block_lines() below on that.
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
#
# `defdose` is NA on the same terms and for a different reason. NONMEM's default
# when no compartment is marked DEFDOSE is compartment 1, which is exactly what
# ferx does with a dose row it cannot resolve, so the two agree and there is
# nothing for a caller to report. Filling the NA in with 1 here would turn that
# agreement into a claim the control stream never made.
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
  comps   <- character(length(decls))
  defobs  <- NA_integer_
  defdose <- NA_integer_
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
    # DEFDOSE diverges from DEFOBSERVATION at the fourth character, so `DEFD` is
    # the shortest unambiguous prefix and the two anchors cannot cross-match.
    if (is.na(defdose) && any(grepl("^DEFD", toks[-1], ignore.case = TRUE)))
      defdose <- i
  }
  list(index   = defobs,
       name    = if (is.na(defobs)) NA_character_ else comps[defobs],
       defdose = defdose,
       n_comp  = length(decls),
       comps   = comps)
}

# Does $INPUT declare a CMT data item that NONMEM will actually read?
#
# TRUE / FALSE / NA, where NA means the question could not be answered (no
# $INPUT record, or one with no data items) and callers should stay silent
# rather than guess.
#
# This decides which rule assigns a dose row to a compartment. With a CMT item
# the data does, and NONMEM and ferx agree. Without one they diverge: NM-TRAN
# uses $MODEL's DEFDOSE compartment, while ferx-core reads
# `cmt_col.and_then(...).unwrap_or(1)` and uses compartment 1.
#
# `$INP` and not `$INPUT`: NM-TRAN accepts any unambiguous prefix, and `$IN`
# alone is ambiguous with `$INFN` while `$INP` is not.
#
# A data item may be written `CMT`, as a synonym pair `CMT=X` or `X=CMT`, or
# dropped with `CMT=DROP` / `CMT=SKIP`. A dropped item is NOT read by NONMEM, so
# it answers FALSE here -- for the dose-compartment question that is the right
# answer, because NM-TRAN falls back to DEFDOSE exactly as if the column were
# absent.
#
# What this canNOT answer is what the DATASET holds, and callers must not word a
# diagnostic as though it could. $INPUT names columns by POSITION; ferx binds
# them by CSV HEADER NAME and never reads $INPUT at all. So the two answers come
# apart in both directions:
#
#   - FALSE via `CMT=DROP` means NONMEM ignores a column that is physically
#     there, and ferx will happily read it. "ferx doses compartment 1" is false
#     for that spelling, and "add a CMT column" sends the user to add one that
#     already exists.
#   - TRUE means only that $INPUT named an item `CMT`. If the CSV header spells
#     that column something else, ferx still falls to compartment 1 and this
#     answers TRUE, so the divergence goes unreported. Tracked as #36.
#
# Settling either needs the dataset. `.extract_nm_data_path()` could supply it
# and a one-line header read would close the second case; until then the
# caller's remedy says what the dataset must CONTAIN rather than claiming to
# have looked.
.nm_input_has_cmt <- function(ctl_file) {
  body <- .nm_block_lines(ctl_file, "INP")
  if (length(body) == 0L) return(NA)
  # Items may continue onto following lines, which is why this reads the whole
  # record body rather than the `$INPUT` line.
  txt  <- sub("^\\s*[$]INP[A-Za-z]*", "", paste(body, collapse = " "),
              ignore.case = TRUE)
  # Whitespace around a synonym pair's `=` is collapsed FIRST. Splitting on
  # whitespace before parsing the `=` handed `CMT = DROP` back as three tokens,
  # of which the bare `CMT` answered TRUE -- so the DROP rule below applied to
  # `CMT=DROP` and not to `CMT = DROP`, and one spelling of one item silently
  # changed the answer.
  txt  <- gsub("[[:space:]]*=[[:space:]]*", "=", txt)
  toks <- trimws(strsplit(txt, "[,[:space:]]+")[[1]])
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0L) return(NA)
  any(vapply(toks, .nm_item_is_cmt, TRUE))
}

.nm_item_is_cmt <- function(tok) {
  parts <- toupper(trimws(strsplit(tok, "=")[[1]]))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) return(FALSE)
  if (any(parts %in% c("DROP", "SKIP"))) return(FALSE)
  any(parts == "CMT")
}
