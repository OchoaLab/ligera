# Full BOLT-like trick implemented!
# though only BEDMatrix makes sense in practice, for internal tests it makes more sense to admit regular kinship matrices as inputs too.
# computes inbreeding vector too, which is not used internally here but LIGERA needs it later (it's most convenient to compute it here instead of having to scan the genome again).
#
# MISSINGNESS: `Y` should be subsetted already (in LIGERA), but `X` can't be, so we subset that here using `indexes_ind`.
#
# my basic pure R implementation of the conjugate gradient method
# for educational purposes, as a stepping stone toward the BOLT-like trick
# here I use the kinship/popgen language so it makes sense to me
#
# finds a `Z` that satisfies: kinship %*% Z = Y
# same, but faster than: solve( kinship, Y )
# but here kinship is implicitly given by X
conj_grad_scan <- function(
                           X,
                           Y,
                           mean_kinship,
                           indexes_ind = NULL,
                           tol = 1e-15,
                           n_data_cut = 1e6, # block size limit for BEDMatrix
                           loci_on_cols = FALSE,
                           verbose = FALSE
                           ) {
    # validations
    if ( missing( X ) )
        stop( '`X` genotype matrix is required!' )
    if ( missing( Y ) )
        stop( '`Y` covariates matrix is required!' )
    if ( missing( mean_kinship ) )
        stop( '`mean_kinship` is required!' )

    # validate X and Y
    if ( !is.matrix(X) ) # TRUE for BEDMatrix
        stop( '`X` must be a matrix!' )
    if ( !is.matrix(Y) )
        stop( '`Y` must be a matrix!' )
    # force loci on columns here
    if ( 'BEDMatrix' %in% class( X ) )
        loci_on_cols <- TRUE

    # validate mean_kinship, must be numeric scalar
    if ( !is.numeric( mean_kinship ) )
        stop( '`mean_kinship` must be numeric!' )
    if ( length( mean_kinship ) != 1 )
        stop( '`mean_kinship` must be scalar!  Got length "', length( mean_kinship ), '"' )

    # get dimensions from Y
    n_ind <- nrow( Y )
    k_covars <- ncol( Y )
    # NOTE: X gets validated (repeatedly) in popkin_prod

    # starting point for solution is all zeroes
    # same dim as Y
    Z <- matrix(
        0,
        nrow = n_ind,
        ncol = k_covars
    )
    
    # Z_final stores converged columns only (unconverged are zeroes until they converge)
    Z_final <- Z
    
    # other internal matrices, which get updated as we go (columns drop as we converge)
    # residuals, initial values, updating as we progress
    # R <- Y - kinship %*% Z
    R <- Y
    # conjugate vectors (initial values)
    P <- R
    # residual norms
    Rn <- colSums( R^2 )
    # different covariate columns may converge at different times, let's keep track of that
    not_converged <- rep.int( TRUE, k_covars )

    # things that get computed by popkin_prod first time
    b <- NA # tells popkin_prod that we need it!
    inbr <- NULL

    if ( verbose )
        iter <- 0
    
    # start loop
    while ( any( not_converged ) ) {
        # P and R matrices are always non-converged subsets!

        # for info
        if ( verbose ) {
            iter <- iter + 1
            message( 'iter: ', iter )
            message( 'Rn: ', toString( Rn ) )
        }

        # NOTE: this is the slowest part!
        obj <- popkin_prod(
            X = X,
            P = P,
            mean_kinship = mean_kinship,
            b = b,
            indexes_ind = indexes_ind,
            n_data_cut = n_data_cut, # block size limit for BEDMatrix
            loci_on_cols = loci_on_cols
        )
        # always want KP
        KP <- obj$KP
        if ( is.na( b ) ) {
            # here we replace b and inbr with the value form the function
            b <- obj$b
            inbr <- obj$inbr
        }
        
        # now continue to update variables in the CG iterations, vectorized if possible
        
        # another vector of the same length
        alpha <- Rn / colSums(P * KP)
        # multiply matrices by alpha along rows (instead of columns, which is R default)
        alpha_mat <- matrix(
            alpha,
            nrow = n_ind,
            ncol = length( alpha ),
            byrow = TRUE
        )
        Z <- Z + P * alpha_mat
        R <- R - KP * alpha_mat
        # new residuals vector
        Rn1 <- colSums( R^2 )
        # take action if something has converged!
        new_converged <- Rn1 < tol
        if ( any( new_converged ) ) {
            # determine newly converged columns in terms of Z_final columns
            new_converged_in_final <- which( not_converged )[ new_converged ]
            # transfer converged columns from Z to Z_final
            Z_final[ , new_converged_in_final ] <- Z[ , new_converged ]
            # subset Z,P,R,Rn so unconverged columns are left only
            # matrices shouldn't drop to vectors (or matrix products die)
            Z <- Z[ , !new_converged, drop = FALSE ]
            P <- P[ , !new_converged, drop = FALSE ]
            R <- R[ , !new_converged, drop = FALSE ]
            Rn <- Rn[ !new_converged ]
            Rn1 <- Rn1[ !new_converged ]
            # update not_converged indicators
            not_converged[ new_converged_in_final ] <- FALSE
            # save a little bit of time in the last iteration by returning after this happens
            if ( !any( not_converged ) )
                break
        }
        # if there are still unconverged things, keep updating things
        # multiply `Rn1 / Rn` too across rows instead of columns
        P <- R + P * matrix(
                         Rn1 / Rn,
                         nrow = n_ind,
                         ncol = length( Rn1 ),
                         byrow = TRUE
                     )
        Rn <- Rn1
    }
    
    # after everything has converged, return the matrix of interest!
    return(
        list(
            Z = Z_final,
            inbr = inbr
        )
    )
}
