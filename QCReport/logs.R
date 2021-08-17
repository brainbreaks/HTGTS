log_input = function(input) {
  vals = shiny::reactiveValuesToList(input)
  vals_filter = !grepl("shinyActionButtonValue", sapply(vals, function(z) {paste(class(z), collapse="|")}))
  vals = vals[vals_filter]
  vals_pretty = sapply(vals, function(z) {
    if(is.data.frame(z) & "name" %in% colnames(z)) z = paste(z$name, collapse=", ")
    if(is.null(z)) z = "NULL"
    as.character(z)
  })

  log(names(vals_pretty), "=", vals_pretty, collapse="\n")
}