#' Make correlation vector plot
#' 
#' A function to visualise the differences between different found biclusters. 
#' Output is a matrix of plots. Each correlation vector is plotted against each 
#' other across the entire  measured gene set in the lower diagonal plots, and
#' a chosen gene set (e.g. mitochondrial) in the upper diagonal plots. The diagnal
#' plots themselves show the density plots of the entire measured and chosen
#' gene set. There are addition options to set the transparancy of the data
#' points and names of the correlation vectors.
#' 
#' @param cv.df A dataframe containing the correlation vectors of one or more patterns.
#' @param geneset.loc A gene set of interest (e.g. mitochondrial) to be plotted separately from rest of genes.
#' @param geneset.name Name of geneset (e.g. mitochondrial genes)
#' @param alpha1 Transparency level of non-gene set genes
#' @param alpha2 Transparency level of gene set genes
#' @param cnames Character vector containing names for the correlation vector
#' @return A plot of the correlation vectors
#' @example example_code/example_cvplot.R
#' @export


CVPlot <- function(cv.df,geneset.loc, geneset.name, alpha1 = 0.005,
                   alpha2 = 0.1, cnames=NULL){
    if(is.character(geneset.name)==FALSE){
        stop("geneset.name must be a single character value")
    }
  
    cv.num <- dim(cv.df)[2]
    if(length(cnames) != cv.num){
        colnames(cv.df)[seq_len(cv.num)] <- paste(rep("CV",cv.num),
                                              seq_len(cv.num),sep = "")
        }else{
        colnames(cv.df) <- cnames
        }
  
    mito_status <- rep(paste("Non", geneset.name),dim(cv.df)[1])
    mito_status[geneset.loc] <- geneset.name
    cv.df$Status <- mito_status
  
    status.loc <- (colnames(cv.df) == "Status")
  
    if(dim(cv.df)[2] == 2){
        a1 <- colnames(cv.df)[1]
        p <- ggplot(cv.df,aes(x = get(a1), col=Status)) + xlab(a1)
        return(p + geom_density() )
    }
  
    custom_cv_plot <- ggpairs(cv.df[,!status.loc],upper = "blank",
                              lower = "blank",
                              title = "",axisLabels ="show", legends=TRUE)
  
    col_cv <- hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0,
                      direction = 1)(2)
  
    H_plot_fun <- function(a,b){
        a1 <- colnames(cv.df)[a]
        b1 <- colnames(cv.df)[b]
        p <- ggplot(cv.df[-geneset.loc,],aes(x = get(a1), y = get(b1))) +
          xlab(a1) + ylab(b1)
        return(p + geom_point(size = 2,alpha = alpha1,col=col_cv[2]))
        }
  
    H_combn <- combn(cv.num,2)
  
    H_plot_fun2 <- function(i){
        p1 <- H_plot_fun(H_combn[1,i],H_combn[2,i])
        return(p1)
    }
    p_H_plots <- lapply(seq_len(dim(H_combn)[2]),H_plot_fun2)
  
    for(i in seq_len(dim(H_combn)[2])){
        custom_cv_plot <- putPlot(custom_cv_plot,p_H_plots[[i]],
                                  H_combn[2,i],H_combn[1,i])}
  
    H_m_plot_fun <- function(a,b){
        a1 <- colnames(cv.df)[a]
        b1 <- colnames(cv.df)[b]
        p <- ggplot(cv.df[geneset.loc,], aes(x = get(a1), y = get(b1))) +
          xlab(a1) + ylab(b1)
        return(p + geom_point(size = 2,alpha = alpha2,col=col_cv[1]))
    }
  
    H_m_plot_fun2 <- function(i){
        p1 <- H_m_plot_fun(H_combn[1,i],H_combn[2,i])
        return(p1)
    }
  
    p_mH_plots <- lapply(seq_len(dim(H_combn)[2]),H_m_plot_fun2)
  
    for(i in seq_len(dim(H_combn)[2])){
        custom_cv_plot <- putPlot(custom_cv_plot,p_mH_plots[[i]],
                                  H_combn[1,i],H_combn[2,i])}
  
  
    H_d_plot_fun <- function(a,l1=FALSE){
        a1 <- colnames(cv.df)[a]
        p <- ggplot(cv.df, aes(x = get(a1), col = Status)) +
          xlab(a1)
        if(l1 == FALSE){
            return(p + geom_density() + theme(legend.position="none"))}
        else{
            return(p + geom_density())}
    }
  
  
    for(i in seq_len(cv.num-1)){
        custom_cv_plot <- putPlot(custom_cv_plot,H_d_plot_fun(i),i,i)}
  
    custom_cv_plot <- putPlot(custom_cv_plot,H_d_plot_fun(cv.num,l1=FALSE),
                              cv.num,cv.num)
  
    return(custom_cv_plot)
}