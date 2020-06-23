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
  
  if(identical(ass, group)){
    warning("No Pop data in initial input")
    grp <- match("Group", colnames(dat))
    dat <- dat[,-grp]
  }
  
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
#' @param pattern vector of length 2 describing the pattern and replacement for gsub() on the individual names to define groups.
#' the default regex takes the first letter from the vector of ID's.
#' @param order string of either byname (order by ID name), byq (order by q value of the first inferred cluster) or none, which orders the bars of the plot.
#' @keywords STRUCTURE, parse, plot
#' @export

plot.PStructure <- function(x, 
                            pop_names = NULL, 
                            col = NULL, 
                            yadjust = -0.03, 
                            pattern = c("(^.{1}).*", "\\1"),
                            order = "byname"){
  # fixme
  if(!is.null(pop_names) & !length(pop_names) %in% c(0, attributes(x)$k)) stop("Population names does not equal k.")
  if(!is.null(col) & !length(col) %in% c(0, attributes(x)$k)) stop("Colour vector must be same length as number of populations.")
  if(!is.character(order) & !order %in% c("byname", "byq", "none")) stop("Order must be one of byname, byq, or none.")
  
  if(is.matrix(x)){
    barplot(x)
  } else {
    pmat <- t(x[grepl("IC", colnames(x))])
    # fixme
    colnames(pmat) <- if(!is.null(x$Group)){
      paste(x$Group, x$ID, sep = "-")
    } else {
        x$ID
    }
    # order
    if(order == "byname"){
      o_pmat <- pmat[,order(colnames(pmat))]
    } else if(order == "byq"){
      o_pmat <- pmat[,order(as.numeric(pmat[1,]), colnames(pmat))]
    } else {
      o_pmat <- pmat
    }
   
    # get groups. FIXME
    if(!is.null(x$Group)){
    groups <- table(sapply(strsplit(colnames(o_pmat), "-"), "[[", 2))
    } else {
      ## some regex
      groups <- table(gsub(pattern = pattern[1], replacement = pattern[2], colnames(o_pmat)))
    }
    prop_groups <- groups/sum(groups)
    # add plot
    par(xpd=NA)
    if(is.null(col)) barplot(o_pmat, border = 1, space = 0, axes = FALSE, xaxt='n')
    if(!is.null(col)) barplot(o_pmat, border = 1, space = 0, axes = FALSE, xaxt='n', col=col)
    # vector to add text
    l <- dim(o_pmat)[2]
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
        text(x = text, y=yadjust, label = pop_names[i], font=1)
      } else if(is.null(pop_names)){
        text(x = text, y=yadjust, label = paste("Pop", i, sep = "-"), font=1) 
      }
      i <- i + 1
    }
    
  }
  
}

#' Print parsed STRUCTURE output.
#'
#' @param x R object that inherits PStructure (parsed STRUCTURE).
#' @keywords STRUCTURE, print
#' @export

print.PStructure <- function(x){
  to_print <- paste("Object of class PStructure with ",
        attributes(x)$k,
        " inferred clusters and ",
        dim(x)[1],
        " individuals.", sep = "")
  cat(to_print, "\n\n")
  no_rows_left <- if(dim(x)[1] > 5){
    dim(x)[1] - 5
  } else {
    NULL
  }
  print.data.frame(x[1:5,])
  to_print2 <- paste("  ...", no_rows_left, " rows truncated", sep = "")
  cat(to_print2)
}

