# this version is used for tests only, not used in production
get_proj_denom_basic <- function(Z, trait) {
    # this is direct calculation for trait (including implicit intercept, coefficient not calculated; no other covariates)
    PhiInvy <- Z[ , 1 ]
    PhiInv1 <- Z[ , 2 ]
    
    # precompute quantities shared across loci
    PhiInv11 <- sum( PhiInv1 )
    PhiInv1y <- sum( PhiInvy )
    PhiInvyy <- drop( trait %*% PhiInvy )
    # a denominator that recurs
    denom <- PhiInvyy * PhiInv11 - PhiInv1y^2

    # the projection vector
    proj <- ( PhiInv11 * PhiInvy - PhiInv1y * PhiInv1 ) / denom

    # rephrase as variance, for direct comparison to a more complicated approach for covariates
    var_trait <- PhiInv11 / denom
    
    # done, return!
    return(
        list(
            proj = proj,
            var = var_trait
        )
    )
}
