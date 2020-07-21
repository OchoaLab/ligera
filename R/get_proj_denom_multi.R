# @param trait_index The index for the trait in matrix Y (its column number), only needed if `trait_only == TRUE`
get_proj_denom_multi <- function(Z, Y, trait_only = TRUE, trait_index = 1) {
    # this is generalized calculation for trait, intercept, and possibly more covariates

    # make sure dimensions match
    if (
        nrow( Z ) != nrow( Y ) &&
        ncol( Z ) != ncol( Y )
    )
        stop('`Z` and `Y` dimensions must match!  Got dim(Z) = ', toString(dim(Z)), ', dim(Y) = ', toString(dim(Y)) )
    
    # this is the generalized "denom"
    # direct version, hoping this doesn't happen after fixed locus removal
    cov_mat <- solve( crossprod( Y, Z ) )
    ## # "try" matrix inversion and try to recover if there is colinearity
    ## cov_mat <- try( solve( crossprod( Y, Z ) ), silent = TRUE )
    ## while ( "try-error" %in% class( cov_mat ) ) {
    ##     ## # makes sense to throw away columns, one at the time
    ##     ## print( Y )
    ##     ## print( Z )
    ##     # keep all rows
    ##     # find index to remove
    ##     k <- ncol( Y )
    ##     if ( k == 1 ) {
    ##         # NOTE: if we were down to a single column, that's always invertible (unless the dot product is zero), very unlikely (let's actually die if that happens)
    ##         stop('Removed all columns trying to make `crossprod(Y, Z)` invertible!')
    ##     } else {
    ##         # remove columns now
    ##         Y <- Y[ , -k ]
    ##         Z <- Z[ , -k ]
    ##         # try again!
    ##         cov_mat <- try( solve( crossprod( Y, Z ) ) )
    ##     }
    ## }

    # this is the generalized projection matrix
    proj <- Z %*% cov_mat

    if (trait_only) {
        # here return values for the trait only, replicating old behavior (of `get_proj_denom`)

        # extract column that corresponds to trait (as determined by user)
        proj_trait <- proj[ , trait_index ]
        # and denominator is the single scalar entry of this matrix (there are also covariances), inverted!
        var_trait <- cov_mat[ trait_index, trait_index ]
        
        # done, return!
        return(
            list(
                proj = proj_trait,
                var = var_trait
            )
        )

    } else {
        # done, return full matrices (coefficients and variance calculations for all fixed covariates, including trait and intercept)
        return(
            list(
                proj = proj,
                cov_mat = cov_mat
            )
        )
    }
}
