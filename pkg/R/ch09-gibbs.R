# function

normgibbs<-function (N, n, a, b, cc, d, xbar, ssquared) 
{
        mat = matrix(ncol = 2, nrow = N)
        mu = cc
        tau = a/b
        mat[1, ] = c(mu, tau)
        for (i in 2:N) {
                muprec = n * tau + d
                mumean = (d * cc + n * tau * xbar)/muprec
                mu = rnorm(1, mumean, sqrt(1/muprec))
                taub = b + 0.5 * ((n - 1) * ssquared + n * (xbar - 
                        mu)^2)
                tau = rgamma(1, a + n/2, taub)
                mat[i, ] = c(mu, tau)
        }
        mat
}


# example

postmat=normgibbs(N=11000,n=15,a=3,b=11,cc=10,d=1/100,xbar=25,ssquared=20)
postmat=postmat[1001:11000,]
op=par(mfrow=c(3,3))
plot(postmat,col=1:10000)
plot(postmat,type="l")
plot.new()
plot(ts(postmat[,1]))
plot(ts(postmat[,2]))
plot(ts(sqrt(1/postmat[,2])))
hist(postmat[,1],40)
hist(postmat[,2],40)
hist(sqrt(1/postmat[,2]),40)
par(op)



# eof

