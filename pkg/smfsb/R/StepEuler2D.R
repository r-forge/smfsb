require(abind)

StepEuler2D <- function(N, d, dt=0.01) {
    S = t(N$Post-N$Pre)
    v = ncol(S)
    u = nrow(S)
    left <- function(a) {
        m = dim(a)[2]
        abind(a[,2:m,],a[,1,],along=2)
    }
    right <- function(a) {
        m = dim(a)[2]
        abind(a[,m,],a[,1:(m-1),],along=2)
    }
    down <- function(a) {
        n = dim(a)[3]
        abind(a[,,2:n],a[,,1],along=3)
    }
    up <- function(a) {
        n = dim(a)[3]
        abind(a[,,n],a[,,1:(n-1)],along=3)
    }
    laplacian <- function(a) left(a) + right(a) + up(a) + down(a) - 4*a
    rectify <- function(a) {
        a[a<0] = 0 # absorb at 0
        a
    }
    diffuse <- function(a) {
        m = dim(a)[2]
        n = dim(a)[3]
        a = a + d*laplacian(a)*dt
        a = rectify(a)
        a
    }
    return(function(x0, t0, deltat, ...) {
        x = x0
        t = t0
        m = dim(x0)[2]
        n = dim(x0)[3]
        termt = t0 + deltat
        repeat {
            x = diffuse(x)
            hr = apply(x, c(2,3), function(x){N$h(x, t, ...)})
            ## TODO: replace for loops with an apply command?
            for (i in 1:m) {
                for (j in 1:n) {
                    x[,i,j] = x[,i,j] + S %*% (hr[,i,j]*dt)
                }
            }
            x = rectify(x)
            t = t + dt
            if (t > termt)
                return(x)
        }
    })
}

