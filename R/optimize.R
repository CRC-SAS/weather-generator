## least squares function for kriging parameters
LS <- function(variogram, dist_matrix, params) {
    M <- params[2] * exp((-dist_matrix)/params[3])
    diag(M) <- params[1] + params[2]
    return(sum(variogram - M)^2)
}

partially_apply_LS <- function(variogram, dist_matrix, base_p = c()) {
    new_ls <- function(p) {
        p <- c(base_p, p)
        LS(variogram, dist_matrix, p)
    }
    return(new_ls)
}

# Ver:
# http://www.win-vector.com/blog/2014/05/trimming-the-fat-from-glm-models-in-r/
# Es una fnción que borra del resultado de mgcv::gam, todo lo que no necesita predict!!
stripGlmLR = function(cm) {
    cm$y = c()
    cm$model = c()

    cm$residuals = c()
    cm$fitted.values = c()
    cm$effects = c()
    cm$qr$qr = c()
    cm$linear.predictors = c()
    cm$weights = c()
    cm$prior.weights = c()
    cm$data = c()


    cm$family$variance = c()
    cm$family$dev.resids = c()
    cm$family$aic = c()
    cm$family$validmu = c()
    cm$family$simulate = c()
    attr(cm$terms,".Environment") = c()
    attr(cm$formula,".Environment") = c()

    cm
}
