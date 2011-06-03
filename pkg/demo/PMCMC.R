# PMCMC.R

require(smfsb)

noiseSD=10

# first simulate some data
truth=simTs(c(x1=50,x2=100),0,20,2,stepLVc)
message("True states")
print(truth)
data2=truth+rnorm(prod(dim(truth)),0,noiseSD)
data2=as.timedData(data2)
data1=as.matrix(data2[,1])
colnames(data1)="x1"
message("Simulated data")
print(data2)

# now define the data likelihood functions
data1Lik <- function(x,t,y,log=TRUE,...)
{
	with(as.list(x),{
		return(dnorm(y,x1,noiseSD,log))
	})
}

data2Lik <- function(x,t,y,log=TRUE,...)
{
	ll=sum(dnorm(y,x,noiseSD,log=TRUE))
	if (log)
		return(ll)
	else
		return(exp(ll))
}

# now define a sampler for the prior on the initial state
simx0 <- function(N,t0,...)
{
	mat=cbind(rpois(N,50),rpois(N,100))
	colnames(mat)=c("x1","x2")
	mat
}

# create marginal log-likelihood functions, based on a particle filter
mLLik1=pfMLLik(100,simx0,0,stepLVc,data1Lik,data1)
mLLik2=pfMLLik(100,simx0,0,stepLVc,data2Lik,data2)

# Now create an MCMC algorithm...
iters=100
tune=0.01
thin=10
th=c(th1 = 1, th2 = 0.005, th3 = 0.6)

p=length(th)
ll=-1-99
thmat=matrix(0,nrow=iters,ncol=p)
colnames(thmat)=names(th)
for (i in 1:iters) {
	message(paste(i,""),appendLF=FALSE)
	for (j in 1:thin) {
		thprop=th*exp(rnorm(p,0,tune))
		llprop=mLLik1(thprop)
		if (log(runif(1)) < llprop - ll) {
			th=thprop
			ll=llprop
			}
		}
	thmat[i,]=th
	}
message("Done!")

# plot the results
mcmcsummaries <- function(mat,plot=TRUE)
  {
    d=dim(mat)
    p=d[2]
    message(paste("N =",d[1],"iterations"))
    print(summary(mat))
    message("Standard deviations:")
    print(apply(mat,2,sd))
    op=par(mfrow=c(p,3))
    names=colnames(mat)
    for (i in 1:p) {
      plot(ts(mat[,i]),main=names[i])
      acf(mat[,i],lag.max=100,main=names[i])
      hist(mat[,i],30,main=names[i])
    }
    par(op)
  }

mcmcsummaries(thmat)

# eof

