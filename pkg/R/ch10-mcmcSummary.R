# mcmcSummary

mcmcSummary <- function(mat,rows=4,plot=TRUE)
  {
    d=dim(mat)
    p=d[2]
    message(paste("N =",d[1],"iterations"))
    print(summary(mat))
    message("Standard deviations:")
    print(apply(mat,2,sd))
    op=par(mfrow=c(rows,3))
    names=colnames(mat)
    for (i in 1:p) {
      plot(ts(mat[,i]),main=names[i],ylab="Value",xlab="Iteration")
      acf(mat[,i],lag.max=100,main=names[i])
      hist(mat[,i],30,main=names[i],xlab="Value",freq=FALSE)
    }
    par(op)
  }


# eof
