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
conj_grad_scan_bed_wcpp <- function(
                                    file,
                                    m_loci,
                                    n_ind, # number of individuals in BED file (may be less than kept individuals in Y or indexes_ind)
                                    Y,
                                    b,
                                    indexes_ind = NULL,
                                    tol = 1e-15
                                    ) {
    # validations
    if ( missing( file ) )
        stop( 'Genotype BED file `file` is required!' )
    if ( missing( Y ) )
        stop( '`Y` covariates matrix is required!' )
    if ( missing( b ) )
        stop( '`b` scalar is required!' )

    # validate Y, b
    if ( !is.matrix(Y) )
        stop( '`Y` must be a matrix!' )
    if ( length( b ) != 1 )
        stop( '`b` must be scalar!  Got length ', length( b ) )
    if ( !is.numeric(b) )
        stop( '`b` must be numeric!' )

    # get dimensions from Y
    n_ind_kept <- nrow( Y )
    k_covars <- ncol( Y )
    # NOTE: the BED file gets validated (repeatedly) in popkin_prod_bed_cpp
    # indexes_ind should match n_ind_kept...
    if ( !is.null( indexes_ind ) ) {
        if ( length( indexes_ind ) != n_ind )
            stop( 'Number of individuals disagrees between BED (', n_ind, ') and length of indexes_ind (', length( indexes_ind ), ')!' )
        if ( sum( indexes_ind ) != n_ind_kept )
            stop( 'Number of individuals kept disagrees between Y (', n_ind_kept, ') and indexes_ind (', sum( indexes_ind ), ')!' )
    }
    

    # starting point for solution is all zeroes
    # same dim as Y
    # Z will hold unconverged subset of data
    Z <- matrix(
        0,
        nrow = n_ind_kept,
        ncol = k_covars
    )
    # Z_final stores converged columns only (unconverged are zeroes until they converge)
    Z_final <- Z
    # other internal matrices, which get updated as we go (columns drop as we converge, like Z but not Z_final)
    # residuals, initial values, updating as we progress
    # R <- Y - kinship %*% Z
    R <- Y
    # conjugate vectors (initial values)
    P <- R
    # residual norms
    Rn <- colSums( R^2 )
    # different covariate columns may converge at different times, let's keep track of that
    # this corresponds to Z_final colums (not Z columns because they drop out as we go)
    not_converged <- rep.int( TRUE, k_covars )
    
    # start loop
    while ( any( not_converged ) ) {
        # Z, P and R matrices are always non-converged subsets!

        # NOTE: this is the slowest part!
        KP <- popkin_prod_bed_cpp(
            file,
            m_loci,
            n_ind, # number of individuals in BED file, not kept
            P,
            b,
            indexes_ind
        )
        
        # now continue to update variables in the CG iterations, vectorized if possible
        
        # another vector of the same length
        alpha <- Rn / colSums(P * KP)
        # sweep makes alpha multiply every row of P, KP (normal product is by columns)
        Z <- Z + sweep( P, 2, alpha, '*')
        R <- R - sweep( KP, 2, alpha, '*')
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
            # matrices shouldn't drop to vectors (sweep dies)
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
        # have to "sweep" the `beta = Rn1 / Rn` too, to go across rows instead of columns
        P <- R + sweep( P, 2, Rn1 / Rn, '*')
        Rn <- Rn1
    }
    
    # after everything has converged, return the matrix of interest!
    return( Z_final )
}
