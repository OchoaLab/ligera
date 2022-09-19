# wrapper around cPCG::cgsolve that does full matrices instead of just vectors, iteratively
# `transpose = TRUE` calculates t( cgsolve_mat( kinship, t(Y) ) ) but smartly
cgsolve_mat <- function(
                        kinship,
                        Y,
                        transpose = FALSE,
                        # cgsolve options
                        tol = 1e-15, # default 1e-6
                        maxIter = 1e6 # default 1e3
                        ) {
    # initialize matrix to fill, same dimensions as input (even for `transpose = TRUE`)
    Z <- matrix( NA, nrow = nrow(Y), ncol = ncol(Y) )
    if ( transpose ) {
        # fill it one row-turned-column at the time
        for ( k in 1 : nrow(Y) ) {
            if ( all( Y[ k, ] == 0 ) ) {
                # handle an annoying edge cases that is expected to occur in genetic data if there are fixed loci (residuals are zero there!)
                # sadly, in these cases cgsolve returns NaN instead of the obvious zero solution
                Z[ k, ] <- 0
            } else {
                Z[ k, ] <- drop( cPCG::cgsolve( kinship, Y[ k, ], tol = tol, maxIter = maxIter ) )
            }
        }
    } else {
        # fill it one column at the time
        for ( k in 1 : ncol(Y) ) {
            if ( all( Y[ , k ] == 0 ) ) {
                # handle an annoying edge cases that is expected to occur in genetic data if there are fixed loci (residuals are zero there!)
                # sadly, in these cases cgsolve returns NaN instead of the obvious zero solution
                Z[ , k ] <- 0
            } else {
                Z[ , k ] <- drop( cPCG::cgsolve( kinship, Y[ , k ], tol = tol, maxIter = maxIter ) )
            }
        }
    }
    return( Z )
}
