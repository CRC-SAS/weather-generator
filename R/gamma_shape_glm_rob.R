gamma.shape.glm.rob <- function(object, it.lim = 10,
    eps.max = .Machine$double.eps^0.25,
    verbose = FALSE, ...)
{
    if(is.null(object$y)) object <- update(object, y = TRUE)
    y <- object$y
    A <- object$prior.weights
    if(is.null(A)) A <- rep(1, length(y))
    u <- object$fitted.values
    alpha <- 1/object$dispersion
    if(verbose) {
        message(gettextf("Initial estimate: %s", format(alpha)), domain = NA)
        utils::flush.console()
    }
    fixed <-  -y/u - log(u) + log(A) + 1 + log(y + (y == 0))
    eps <- 1
    itr <- 0
    while(abs(eps) > eps.max && (itr <- itr + 1) <= it.lim) {
        sc <- sum(A * (fixed + log(alpha) - digamma(A * alpha)))
        inf <- sum(A * (A * trigamma(A * alpha) - 1/alpha))
        alpha <- alpha + (eps <- sc/inf)
        if(verbose) {
            message(gettextf("Iter. %d Alpha: %s", itr, format(alpha)),
                domain = NA)
            utils::flush.console()
        }
    }
    if(itr > it.lim) warning("iteration limit reached")
    res <- list(alpha = alpha, SE = sqrt(1/inf))
    class(res) <- "gamma.shape"
    res
}
