insert_line_breaks <- function(text, line_length) {
  words <- strsplit(text, " ")[[1]]
  lines <- character(0)
  current_line <- ""
  for (word in words) {
    if (nchar(current_line) + nchar(word) + 1 <= line_length) {
      if (current_line == "") {
        current_line <- word
      } else {
        current_line <- paste(current_line, word)
      }
    } else {
      lines <- c(lines, current_line)
      current_line <- word
    }
  }
  if (current_line != "") {
    lines <- c(lines, current_line)
  }
  return(paste(lines, collapse = "\n"))
}
# 
# text <- "Il faudra peut être raccourcir par la suite."
# formatted_text <- insert_line_breaks(text, 30)
# print(formatted_text)
