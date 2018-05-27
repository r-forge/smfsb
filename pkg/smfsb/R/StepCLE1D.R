StepCLE1D <- function(N, d, dt=0.01) {
    S = t(N$Post-N$Pre)
    v = ncol(S)
    u = nrow(S)
    sdt = sqrt(dt)
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
        noise = matrix(rnorm(n*u,0,sdt),nrow=u)
        m = m + d*laplacian(m)*dt + sqrt(d)*(
            sqrt(m+forward(m))*noise -
            sqrt(m+back(m))*back(noise))
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
            dwt = matrix(rnorm(n*v,0,sdt),nrow=v)
            x = x + S %*% (hr*dt + sqrt(hr)*dwt)
            x = rectify(x)
            t = t + dt
            if (t > termt)
                return(x)
        }
    })
}


