extract_boot=function (w, smoothing = 1, mr.func = "wmr")
{
	require('doParallel')
    mr.func = match.fun(mr.func)
    mr.obs = mr.func(w, smoothing = smoothing)
    with(w, {
        nloc = length(x)
        nvars = dim(z)[3]
        nscales = length(y)
        exports = c("reps", "wmr", "lagMat", "regularize", "mr.func", 
            "Gauss", "mr.obs", "nscales", "smoothing")
        flibs = c("mvcwt")
        mr.obs$z.boot = foreach(i = 1:nscales, .combine = c, 
            .export = exports, .packages = flibs) %dopar% {
         rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, 
                    nloc)))
                  zp = z[, i, , drop = FALSE] * complex(argument = rphase)
                  mr.boot=as.vector(mr.func(list(x = x, y = y[i], z = zp), 
                    smoothing = smoothing)$z)
	res = foreach(j = 1:nloc, .combine = c) %dopar% {mr.boot[j]}
           	res
        	}
        	dim(mr.obs$z.boot) = c(length(x), length(y), 1)
	return(mr.obs)
    })
}
