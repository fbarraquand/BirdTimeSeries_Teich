mvcwt_test=function (x, y, scale.exp = 0.5, nscales = get.nscales(x), min.scale = get.min.scale(x), 
    max.scale = get.max.scale(x), scales = log2Bins(min.scale, 
        max.scale, nscales), loc = regularize(x), wave.fun = "Morlet") 
{
    s = 1
    wave.fun = match.fun(wave.fun)
    x = as.vector(unlist(x))
    lmat = lagMat(x, loc)
	print(lmat)
    y = matrix(unlist(y), nrow = length(x))
    w = foreach(s = scales) %dopar% {
        Conj(wave.fun(lmat/s)/s^scale.exp) %*% y
    }
    w = array(unlist(w), dim = c(length(loc), ncol(y), length(scales)))
    w = aperm(w, perm = c(1, 3, 2))
    structure(list(x = loc, y = scales, z = w), class = "mvcwt")
}

