# TSGS: Trait specific gene selection using support vector machine and genetic algorithm
Obtaining relevant set of trait specific genes from gene expression data is important for clinical diagnosis of disease and discovery of disease mechanisms in plants and animals. This process involves identification of relevant genes and removal of redundant genes as much as possible from a whole gene set. This package returns the trait specific gene set from the high dimensional RNA-seq count data by applying combination of two conventional machine learning algorithms, support vector machine (SVM) and genetic algorithm (GA). GA is used to control and optimize the subset of genes sent to the SVM for classification and evaluation. Genetic algorithm uses repeated learning steps and cross validation over number of possible solution and selects the best. The algorithm selects the set of genes based on a fitness function that is obtained via support vector machines. Using SVM as the classifier performance and the genetic algorithm for feature selection, a set of trait specific gene set is obtained.

# Authors:
    c(person("Md. Samir", "Farooqi", email = "ms.Farooqi@icar.gov.in", role = "aut"),
    person("K.K.", "Chaturvedi", email = "kk.Chaturvedi@icar.gov.in", role = "aut"),
    person("D.C.", "Mishra", email = "Dwijesh.Mishra@icar.gov.in", role = "aut"),
    person("Sudhir", "Srivastava", email = "Sudhir.Srivastava@icar.gov.in", role = c("cre","aut")))
