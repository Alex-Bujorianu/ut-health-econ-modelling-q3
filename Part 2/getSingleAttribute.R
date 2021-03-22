getSingleAttribute <- function(attribute, output, all=F) sapply(unique(output[,"name"]), function(entity) {
  list(
    if(all) {
      output[output[,"name"]==entity & output[,"key"]==attribute, "value"]
    } else {
      tail(output[output[,"name"]==entity & output[,"key"]==attribute, "value"], n=1)
    }
  )
})