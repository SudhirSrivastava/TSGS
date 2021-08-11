\name{featureSelect}
\alias{featureSelect}
\title{
Trait specific gene selection using support vector machine and genetic algorithm
}
\description{
This function gives the optimal set of informative genes based on RNA-Seq count data
}
\usage{
featureSelect(X, y, p = 20, n.iter = 5, alpha = 0.05, p.adj.method = "bonferroni")
}
\arguments{
  \item{X}{X is a G x N data frame of gene expression values (raw count data) where rows represent genes and columns represent samples. Each cell entry represents the read counts of of a gene in a sample (row names of X as gene names or gene ids)}
  \item{y}{y is a N x 1 numeric vector with entries 0 or 1 representing sample labels, where, 0/1 represents the sample label of samples for two conditions, e.g., 0 for Control and 1 for Case}
  \item{p}{Population size, by default 20}
  \item{n.iter}{The number of iterations, by default 5}
  \item{alpha}{The level of significance, by default 0.05}
  \item{p.adj.method}{Method of adjusting p-values, by default "bonferroni". The other methods available are "BH", "holm", "hochberg", "hommel", "BY".}
}
\value{
  \item{InformativeGenes}{List of informative genes selected}
  \item{LogCPM}{Log cpm data of informative genes}
  \item{DEA_Result}{Differential Expression Analysis Result of informative genes}
}
\author{
c(person("Md. Samir", "Farooqi", email = "ms.Farooqi@icar.gov.in", role = "aut"),
    person("K.K.", "Chaturvedi", email = "kk.Chaturvedi@icar.gov.in", role = "aut"),
    person("D.C.", "Mishra", email = "Dwijesh.Mishra@icar.gov.in", role = "aut"),
    person("Sudhir", "Srivastava", email = "Sudhir.Srivastava@icar.gov.in", role = c("cre","aut")))
}

\examples{
  filename <- system.file("extdata", "exampleData.csv", package = "TSGS")
  cdata <- read.csv(filename, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  X <- as.data.frame(cdata[-1,])
  y <- as.numeric(cdata[1,])
  set.seed(100)
  result <- featureSelect(X, y, 20, 2, 0.05, "bonferroni")
  gene_list <- result$InformativeGenes
  logcpm_data <- result$LogCPM
  dea_result <- result$DEA_Result
}
