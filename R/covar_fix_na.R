# this function fixes NA values, so the output covariates matrix has no NAs at all
# NAs are replaced by the mean value of each column
covar_fix_na <- function( covar ) {

    # if there are no NAs, return immediately
    if ( !anyNA( covar ) )
        return( covar )
    
    # navigate columns
    k_covars <- ncol( covar )
    for ( k in 1 : k_covars ) {
        # extract column vector
        covar_k <- covar[ , k ]
        # test if there are any NAs to bother with in this column
        if ( anyNA( covar_k ) ) {
            # perform replacement
            covar_k[ is.na( covar_k ) ] <- mean( covar_k, na.rm = TRUE )
            # store edited column vector back
            covar[ , k ] <- covar_k
        }
    }
    # done return modified matrix!
    return( covar )
}
