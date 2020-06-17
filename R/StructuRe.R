#' Parse a STRUCTURE output file to a data frame or matrix.
#' 
#' @param file path to the STRUCTURE file.
#' @param matrix boolean; if true, returns only the matrix of inferred clusters.
#' @keywords STRUCTURE, parse
#' @export

parseStructure <- function(file, matrix = FALSE){
  # initial parse
  input_1 <- suppressWarnings(readLines(file(file)))
  input <- input_1[(grep("Inferred ancestry of individuals",input_1)+1):length(input_1)]
  input <- input[2:(grep("Estimated Allele Frequencies",input)-3)]
  # get pops
  k <- as.numeric(gsub(pattern = ".*([0-9]+) .*", replacement = "\\1", input_1[grep("populations assumed", input_1)]))
  # names
  nms <- gsub("  ", "", sapply(strsplit(x = input, split = "   ", fixed = TRUE), "[[", 2))
  # % miss
  miss <- gsub("\\(([0-9]+).*", "\\1", sapply(strsplit(x = input, split = "   ", fixed = TRUE), "[[", 3))
  # the group
  ass <- sapply(strsplit(x = input, split = "   ", fixed = TRUE), "[[", 4)
  group <- gsub("^ ([^ ]+).*", "\\1", x = ass)
  
  # get assignments
  foo<-strsplit(sapply(strsplit(input,":"),function(x) x[2])," ")
  # for quick plotting through barplot()
  if(matrix){
    mat <- sapply(foo,function(x) x[3:length(foo[[2]])])
    class(mat) <- c("PStructure", "matrix")
    return(mat)
  }
  
  foo.mat<-t(sapply(foo,function(x) x[3:length(foo[[2]])]))
  
  dat <- as.data.frame(cbind(nms, miss, group, foo.mat))
  
  # if k > 2
  n_IC <- dim(dat)[2] - 3
  col_IC <- paste("IC", 1:n_IC, sep = "_")
  
  colnames(dat) <- c("ID", "Missing", "Group", col_IC)
  class(dat) <- c("PStructure", "data.frame")
  attr(dat, "k") <- k
  
  return(dat)
}

#' Plot parsed STRUCTURE output.
#'
#' @param x R object that inherits PStructure (parsed STRUCTURE).
#' @param pop_names optional; supply a vector of population names.
#' @param col optional; group colours. Default is greyscale.
#' @param yAdjust optional; adjust the y position of x axis labels.
#' @keywords STRUCTURE, parse, plot
#' @export

plot.PStructure <- function(x, pop_names = NULL, col = NULL, yAdjust = -0.03){
  
  if(!is.null(pop_names) & length(pop_names) != attributes(x)$k) stop("Population names does not equal k.")
  if(!is.null(col) & length(col) != attributes(x)$k) stop("Colour vector must be same length as number of populations.")
  
  if(is.matrix(x)){
    barplot(x)
  } else {
    pmat <- t(x[grepl("IC", colnames(x))])
    colnames(pmat) <- paste(x$Group, x$ID, sep = "-")
    o_pmat <- pmat[,order(colnames(pmat))]
    # get groups
    groups <- table(sapply(strsplit(colnames(o_pmat), "-"), "[[", 1))
    prop_groups <- groups/sum(groups)
    # add plot
    par(xpd=NA)
    if(is.null(col)) barplot(o_pmat, border = 1, space = 0, axes = FALSE, xaxt='n')
    if(!is.null(col)) barplot(o_pmat, border = 1, space = 0, axes = FALSE, xaxt='n', col=col)
    # vector to add text
    lens <- c(0, cumsum(prop_groups)*l)
    lis <- list()
    for(i in 1:length(lens)){
      lis[[i]] <- c(lens[i], lens[i+1])
    }
    text_add <- sapply(lis, mean)
    # loop over groups to delimit 
    i <- 1
    for(text in text_add){
      if(!is.null(pop_names)){
        text(x = text, y = yAdjust, label = pop_names[i], font=3)
      } else if(is.null(pop_names)){
        text(x = text, y = yAdjust, label = paste("Pop", i, sep = "-"), font=3) 
      }
      i <- i + 1
    }
    
  }
  
}
