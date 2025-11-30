StepEuler1D <- function(N, d, dt=0.01) {
    S = t(N$Post-N$Pre)
    v = ncol(S)
    u = nrow(S)
    forward <- function(m) cbind(m[,2:ncol(m)],m[,1])
    back <- function(m) {
        n = ncol(m)
        cbind(m[,n],m[,1:(n-1)])
    }
    laplacian <- function(m) forward(m) + back(m) - 2*m
    rectify <- function(m) {
        m[m<0] = 0 # absorb at 0
        m
    }
    diffuse <- function(m) {
        n = ncol(m)
        m = m + d*laplacian(m)*dt
        m = rectify(m)
        m
    }
    return(function(x0, t0, deltat, ...) {
        x = x0
        t = t0
        n = ncol(x0)
        termt = t0 + deltat
        repeat {
            x = diffuse(x)
            hr = apply(x,2,function(x){N$h(x, t, ...)})
            x = x + S %*% (hr*dt)
            x = rectify(x)
            t = t + dt
            if (t > termt)
                return(x)
        }
    })
}


