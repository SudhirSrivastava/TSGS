featureSelect <- function(X, y, p = 20, n.iter = 5, alpha = 0.05, p.adj.method = "bonferroni"){

  countData1 <- X
  geneNames1 <- rownames(X)
  Labels <- as.numeric(y)
  group <- unique(y)
  z <- factor(Labels, levels = group, labels = group)

  n1 <- length(which(z == group[1]))
  n2 <- length(which(z == group[2]))
  n <- n1+n2

  ## Filtering to remove low count reads
  LS <- colSums(countData1)
  LS.CPM <- LS/10^6
  t <- round(10/min(LS.CPM), 1) # Threshold

  y <- DGEList(counts = countData1, genes = geneNames1)
  keep <- rowSums(cpm(y) > t) >=2
  y <- y[keep, , keep.lib.sizes=FALSE]
  countData <- y$counts
  geneNames <- y$genes
  nGenes <- nrow(countData)

  y <- calcNormFactors(y)
  design <- model.matrix(~z)
  y <- estimateDisp(y, design = design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)

  res2 <- topTags(qlf, n=nGenes)
  res.tab <- res2$table
  ind1 <- which(res.tab$PValue < alpha)
  adj.pval <- p.adjust(res.tab$PValue, method = p.adj.method)
  res.tab$`Adjusted PValue` <- adj.pval
  ind <- which(adj.pval < alpha)
  geneNames.sel <- res.tab$genes[ind]

  res.f <- res.tab[ind,]
  log.counts <- cpm(y$counts, log = TRUE)
  countData.f <- log.counts[ind,] # Final log cpm data for feature selection

  data <- t(countData.f)
  s <- data.frame(data)

  t <- z
  eval_funct <- function(indices){

    evl_df <- cbind(s[,indices==1],t)
    evl_trng <- createDataPartition(evl_df$t, p=0.60,list = FALSE)
    evl_test <- evl_df[-evl_trng,]
    evl_train <- evl_df[evl_trng,]

    evl_svm <- train(t~.,data=evl_train,method="svmRadial",preProc=c("zv"),
                     trControl=trainControl(method = "cv",number = 5),savePredictions = "all")

    evl_cls <- predict(evl_svm,newdata=evl_test)
    evl_tbl <- confusionMatrix(evl_cls,evl_test$t)
    tbl <- evl_tbl$table
    tp <- tbl[1,1]
    tn <- tbl[2,2]
    t <- tbl[1,1]+tbl[1,2]+tbl[2,1]+tbl[2,2]
    result <- -(tn+t-tp)/((tn*t)-(tp*tn))
    return(result)
  }
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
    filter = obj$evaluations == minEval;
    bestObjectCount = sum(rep(1, obj$popSize)[filter]);
    if (bestObjectCount > 1) {
      bestSolution = obj$population[filter,][1,];
    } else {
      bestSolution = obj$population[filter,];
    }
    outputBest = paste(obj$iter, " #selected=", sum(bestSolution),
                       " Best (Error=", minEval, "): ", sep="");
    for (var in 1:length(bestSolution)) {
      outputBest = paste(outputBest,
                         bestSolution[var], " ",
                         sep="");
    }
    outputBest = paste(outputBest, "\n", sep="");

    cat(outputBest);
  }

  woppa <- rbga.bin(size=ncol(s),popSize=p,iters=n.iter, mutationChance=0.30, zeroToOneRatio=20,
                    evalFunc=eval_funct, verbose=TRUE, monitorFunc=monitor)

  bestSolution <- woppa$population[which.min(woppa$evaluations),]
  result <- cbind(data[,bestSolution==1],z)
  feature.selected <- res.f[bestSolution ==1,1]
  logcpm.feature.selected <- t(result)
  ind.m <- fmatch(as.character(feature.selected), as.character(res.f[,1]))
  result.pval <- res.f[ind.m, ]

  list(`InformativeGenes` = feature.selected,
       `LogCPM` = logcpm.feature.selected,
       `DEA_Result` = result.pval)
}
