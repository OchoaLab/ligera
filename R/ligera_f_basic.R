# the original version without BEDMatrix or missingness support
# (for tests; do not export!)
# unlike main function, this one requires explicit kinship inverse to be provided, doesn't use regular kinship though
ligera_f_basic <- function(X, trait, kinship_inv, covar = NULL) {
    # need these dimensions
    m_loci <- nrow(X)
    n_ind <- ncol(X)
    
    # gather matrix of trait, intercept, and optional covariates
    Y <- cbind( trait, 1 )
    # add covariates, if present
    if ( !is.null( covar ) ) {
        # handle NAs now, so final Y has no missingness whatsoever
        covar <- covar_fix_na( covar )
        Y <- cbind( Y, covar )
    }
    # use kinship inverse
    Z <- kinship_inv %*% Y
    
    # proceed getting residuals for alt model first
    # calculate projection matrix
    # Y,Z,H have dimensions n x k
    H <- Z %*% solve( crossprod( Y, Z ) )
    # the coefficients are simply the genotypes projected!
    # (projection for trait coefficient only, alt model only)
    beta <- drop( X %*% H[ , 1 ] )
    # to get SSRs, calculate residuals first
    # (X %*% H) are coefficients, that times Y are fitted X values
    R <- X - tcrossprod( X %*% H, Y )
    # this order of multiplication is most scalable since k << n:
    # (m * n) * (n * k) * (k * n)
    # (X %*% H1) %*% t( Y1 )
    # m*n*k + m*k*n = 2*m*n*k
    # X %*% (H1 %*% t( Y1 ))
    # n^2*k + m*n^2 = n^2*(k+m)
    # then sums of residuals weighted by inverse kinship
    ssr1 <- rowSums( ( R %*% kinship_inv ) * R )
    # degrees of freedom of alt model
    df <- n_ind - ncol( Y )
    
    # repeat for null model now
    # just remove first column of these two
    Y <- Y[ , -1, drop = FALSE ]
    Z <- Z[ , -1, drop = FALSE ]
    H <- Z %*% solve( crossprod( Y, Z ) )
    R <- X - tcrossprod( X %*% H, Y )
    ssr0 <- rowSums( ( R %*% kinship_inv ) * R )
    
    # calculate F statistic now that all the parts are in place!
    f_stat <- (ssr0 - ssr1) / ssr1 * df
    
    # Get p-values for F test!!! (should be exact)
    pval <- stats::pf( f_stat, 1, df, lower.tail = FALSE )
    
    # done, return quantities of interest (nice table!)
    return(
        tibble::tibble(
            pval = pval,
            beta = beta,
            f_stat = f_stat,
            df = df
        )
    )
}
