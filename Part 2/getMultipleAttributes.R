getMultipleAttributes <- function(attributes, output) as.data.frame(sapply(attributes, function(attribute) as.numeric(getSingleAttribute(attribute, output))));
