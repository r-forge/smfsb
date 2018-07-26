## PMCMC.R

## make sure the package is loaded
require(smfsb)

## load the reference data
data(LVdata)

## assume known measurement SD of 10
noiseSD=10

## now define the data likelihood functions
data1Lik <- function(x,t,y,log=TRUE,...)
{
	with(as.list(x),{
		return(dnorm(y,x1,noiseSD,log))
	})
}

data2Lik <- function(x,t,y,log=TRUE,...)
{
	ll=sum(dnorm(y,x,noiseSD,log=TRUE))
	if (log) ll else exp(ll)
}

data3Lik <- function(x,t,y,log=TRUE,...)
{
	ll=sum(dnorm(y,x,th["sd"],log=TRUE))
	if (log) ll else exp(ll)
}

## now define a sampler for the prior on the initial state
simx0 <- function(N,t0,...)
{
	mat=cbind(rpois(N,50),rpois(N,100))
	colnames(mat)=c("x1","x2")
	mat
}

LVdata=as.timedData(LVnoise10)
LVpreyData=as.timedData(LVpreyNoise10)
colnames(LVpreyData)=c("x1")

## create marginal log-likelihood functions, based on a particle filter
mLLik1=pfMLLik(100,simx0,0,stepLVc,data1Lik,LVpreyData)
mLLik2=pfMLLik(100,simx0,0,stepLVc,data2Lik,LVdata)
mLLik3=pfMLLik(100,simx0,0,stepLVc,data3Lik,LVdata)

## Now create an MCMC algorithm...
#th=c(th1 = 1, th2 = 0.005, th3 = 0.6)
th=c(th1 = 1, th2 = 0.005, th3 = 0.6, sd=10)
p = length(th)
rprop = function(th, tune=0.01) { th*exp(rnorm(p,0,tune)) }
thmat = metropolisHastings(th,mLLik3,rprop,iters=1000,thin=10)

## Dump MCMC output matrix to disk
save(thmat,file="LV.RData")

## Compute and plot some basic summaries
mcmcSummary(thmat,truth=th)

# eof

