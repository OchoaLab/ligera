# Full BOLT-like trick implemented!
# though only BEDMatrix makes sense in practice, for internal tests it makes more sense to admit regular kinship matrices as inputs too.
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
#
# this version doesn't calculate b, inbr internally, requires external b instead
conj_grad_scan_bed <- function(
                               X,
                               Y,
                               b,
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
    if ( missing( b ) )
        stop( '`b` scalar is required!' )

    # validate X, Y, b
    if ( !is.matrix(X) ) # TRUE for BEDMatrix
        stop( '`X` must be a matrix!' )
    if ( !is.matrix(Y) )
        stop( '`Y` must be a matrix!' )
    if ( length( b ) != 1 )
        stop( '`b` must be scalar!  Got length ', length( b ) )
    if ( !is.numeric(b) )
        stop( '`b` must be numeric!' )

    # force loci on columns here
    if ( 'BEDMatrix' %in% class( X ) )
        loci_on_cols <- TRUE

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
        KP <- popkin_prod_bed(
            X = X,
            P = P,
            b = b,
            indexes_ind = indexes_ind,
            n_data_cut = n_data_cut, # block size limit for BEDMatrix
            loci_on_cols = loci_on_cols
        )
        
        # now continue to update variables in the CG iterations, vectorized if possible
        
        # another vector of the same length
        alpha <- Rn / colSums(P * KP)
        # sweep makes alpha multiply every row of P, KP (normal product is by columns)
        Z[ , not_converged ] <- Z[ , not_converged ] + sweep( P, 2, alpha, '*')
        R <- R - sweep( KP, 2, alpha, '*')
        # new residuals vector
        Rn1 <- colSums( R^2 )
        # take action if something has converged!
        new_converged <- Rn1 < tol
        if ( any( new_converged ) ) {
            # write to Z if needed
            # subset P,R,Rn so unconverged columns are left only
            still_not_converged <- Rn1 >= tol # columns of not_converged subset
            # if matrices drop to vectors, sweep complains (just below) :(
            P <- P[ , !new_converged, drop = FALSE ]
            R <- R[ , !new_converged, drop = FALSE ]
            Rn <- Rn[ !new_converged ]
            Rn1 <- Rn1[ !new_converged ]
            # update not_converged indicators
            not_converged[ which(not_converged)[ new_converged ] ] <- FALSE
            # save a little bit of time in the last iteration by returning after this happens
            if ( !any( not_converged ) )
                break
        }
        # if there are still unconverged things, keep updating things
        # have to "sweep" the `beta = Rn1 / Rn` too, to go across rows instead of columns
        P <- R + sweep( P, 2, Rn1 / Rn, '*')
        Rn <- Rn1
    }
    
    # after everything has converged, return the matrix of interest!
    return( Z )
}
