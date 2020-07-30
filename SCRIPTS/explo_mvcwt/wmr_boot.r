wmr_boot=function (w, smoothing = 1, reps = 1000, mr.func = "wmr",output="quantile")  #Output in z.boot can be quantile, i.e. Pr(X<=obs) ; 'function', i.e. the ecdf from the surrogates
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
            mr.boot = foreach(j = 1:reps, .combine = cbind, .inorder = FALSE) %dopar% 
                {
                  rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, 
                    nloc)))
                  zp = z[, i, , drop = FALSE] * complex(argument = rphase)
                  as.vector(mr.func(list(x = x, y = y[i], z = zp), 
                    smoothing = smoothing)$z)
		}
		res = foreach(j = 1:nloc, .combine = c) %dopar% 
			{
				if(output=="quantile"){
                			ecdf(mr.boot[j, ])(mr.obs$z[j, i, ])
				}else if(output=="function"){
	                		ecdf(mr.boot[j, ])
				}
			}
           	res
        	}
	if(output=="quantile"){
        	dim(mr.obs$z.boot) = c(length(x), length(y), 1)
        }
	return(mr.obs)
    })
}
