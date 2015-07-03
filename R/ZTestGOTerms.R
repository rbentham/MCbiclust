ZTestGOTerms <- function(genes, gene.values, z.score = FALSE){
require(BSDA)
go.term.genes.loc <- which(genes %in% unique(unlist(GO_term_genes)))

genes.mean <- mean(gene.values[go.term.genes.loc])
genes.sd <-  sd(gene.values[go.term.genes.loc])

go.terms.num <- length(GO_term_genes)
p.values.go <- seq(length = go.terms.num)
if(z.score == TRUE) z.go <- seq(length = go.terms.num)

for(i in seq(length = go.terms.num)){
  genes.loc <- which(genes %in% GO_term_genes[[i]])
 
  if(length(genes.loc) > 2){
    z.test.GO <- z.test(genes[genes.loc], y=NULL, mu = genes.mean,sigma.x = genes.sd)
   p.values.go[i] <- z.test.GO$p.value
   
    if(z.score == TRUE) z.go[i] <- z.test.GO$statistic
   
  }else{
    p.values.go[i] <- 1
    if(z.score == TRUE) z.go[i] <- NA
    if(length(genes_loc)>0){
      mean.go[i] <- mean(gene.values[genes_loc])}
    else{
      mean.go[i] <- NA
    }
  }
}

  if(z.score == TRUE){
   return(z.go)
  } else {
    return(p.adjust(go.pvalues))
  }
}