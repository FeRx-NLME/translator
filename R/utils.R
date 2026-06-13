# Shared string and expression helpers

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
