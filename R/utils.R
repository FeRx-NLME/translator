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
