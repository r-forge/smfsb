# demo for sampling a finite Markov chain

# load the package
library(smfsb)

# set up a Markov transition structure
P=matrix(c(0.9,0.1,0.2,0.8),ncol=2,byrow=TRUE)
pi0=c(0.5,0.5)

# generate a sample
samplepath=rfmc(200,P,pi0)

# analyse the sample
plot(samplepath)
summary(samplepath)
table(samplepath)
table(samplepath)/length(samplepath)

# now compute the exact stationary distribution...
e=eigen(t(P))$vectors[,1]
e/sum(e)


# eof

