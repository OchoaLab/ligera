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
    cov_mat <- solve( crossprod( Y, Z ) )

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
