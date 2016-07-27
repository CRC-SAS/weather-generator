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
