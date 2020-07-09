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
                           loci_on_cols = FALSE
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

    # validate dimensions
    # these are pre-filters
    if ( loci_on_cols ) {
        n_ind <- nrow( X )
        m_loci <- ncol( X )
    } else {
        n_ind <- ncol( X )
        m_loci <- nrow( X )
    }
    # subset and validate
    if ( !is.null( indexes_ind ) ) {
        if ( is.logical( indexes_ind ) ) {
            if ( length( indexes_ind ) != n_ind )
                stop( 'Non-NULL `indexes_ind` length (', length( indexes_ind ), ') does not match number of individuals in genotype matrix (', n_ind, ')!' )
            # now replace n_ind with actual individuals we will work with
            n_ind <- sum( indexes_ind )
        } else {
            stop( 'Non-logical `indexes_ind` currently not supported!' )
        }
    }
    if ( nrow( Y ) != n_ind )
        stop( 'Number of rows of `Y` (', nrow( Y ), ') must match the number of individuals in the (possibly filtered) genotype matrix (', n_ind, ')!' )
    k_covars <- ncol( Y )

    # figure out how many loci to read at the time
    m_chunk <- floor( n_data_cut / n_ind )
    if ( m_chunk < 1 )
        m_chunk <- 1

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

    # compute b (scalar), the minimum B = crossprod(X-1)/m
    # since we're not forming B, we estimate b indirectly, as
    # b = ( 1 - 4 * mean( maf * (1 - maf) ) - mean_kinship ) / ( 1 - mean_kinship )
    # have to calculate it in first pass of genotypes scan, as a running sum
    b <- 0
    b_done <- FALSE
    # also need to estimate the inbreeding coefficient vector!
    inbr <- vector( 'numeric', n_ind )
    inbr_m <- vector( 'numeric', n_ind ) # normalization (number of non-NA loci per individual)
    # (inbr is done when b_done is TRUE as well)
    
    # start loop
    while ( any( not_converged ) ) {
        # P and R matrices are always non-converged subsets!
        
        # the corresponding multiplied version, approximately proportional to: kinship %*% P
        # will be a running average, so initialize to zeroes (same dim as P)
        KP <- matrix(
            0,
            nrow = n_ind,
            ncol = sum( not_converged )
        )
        
        # read the genotype matrix again from the start
        i_start <- 1
        while ( i_start < m_loci ) {
            i_end <- i_start + m_chunk - 1
            if ( i_end > m_loci )
                i_end <- m_loci
            # get submatrix
            X_chunk <-
                if ( !is.null( indexes_ind ) ) {
                    if ( loci_on_cols ) {
                        X[ indexes_ind, i_start : i_end, drop = FALSE ]
                    } else {
                        # transpose so below no other loci_on_cols dependence occurs
                        # (favors BEDMatrix orientation for speed when things get huge, so slower with true matrices in our ordinary layout, but meh)
                        t( X[ i_start : i_end, indexes_ind, drop = FALSE ] )
                    }
                } else {
                    if ( loci_on_cols ) {
                        X[ , i_start : i_end, drop = FALSE ]
                    } else {
                        # transpose so below no other loci_on_cols dependence occurs
                        # (favors BEDMatrix orientation for speed when things get huge, so slower with true matrices in our ordinary layout, but meh)
                        t( X[ i_start : i_end, , drop = FALSE ] )
                    }
                }
            
            # get allele frequencies and ultimately a key product, if it wasn't already calculated
            if ( !b_done ) {
                # get twice the allele frequencies
                # (vector of length m_chunk, approximately)
                x_bar <- colMeans( X_chunk, na.rm = TRUE )
                # compute term of interest and add to running sum
                # NOTE: since x_bar = 2 * maf, then this is as desired
                # x_bar * (2 - x_bar) = 4 * maf * ( 1 - maf )
                b <- b + sum( x_bar * (2 - x_bar) )
            }
            
            # manipulate for kinship-like estimation
            # "center" with -1
            X_chunk <- X_chunk - 1

            if ( !b_done ) {
                # do inbreeding coefficients now
                # here NAs can be handled correctly in reasonable memory
                inbr_m <- inbr_m + rowSums( !is.na( X_chunk ) )
                # now the actual inbreeding estimates
                # NOTE: here we're using the X-1 formulation, so this is really like B, not A
                inbr <- inbr + rowSums( X_chunk^2, na.rm = TRUE )
            }
            
            # TODO: how to we handle NAs properly here???
            # - if NAs are random (non-informative), their missingness affects scale randomly only, and absolute scale doesn't matter, so maybe meh?
            # now set NAs to zero (after centering step)
            X_chunk[ is.na( X_chunk ) ] <- 0

            # now compute desired matrix products, add to running sum
            # inner product has matrices of dimensions (X_chunk is transposed here)
            #   ( m_chunk * n_ind ) * ( n_ind * k_covar ) = ( m_chunk * k_covar )
            # outer product has matrices of dimensions
            #   ( n_ind * m_chunk ) * ( m_chunk * k_covar ) = ( n_ind * k_covar )
            # both products are O( n m k ), so they beat direct kinship estimation which has O( m n^2 )
            KP <- KP + ( X_chunk %*% crossprod( X_chunk, P ) )
            
            # update locus index for next iteration
            i_start <- i_end + 1
        }

        # now we complete the calculation of this key factor, and prevent estimating it again in untold further iterations
        if ( !b_done ) {
            # previous b was just sum(maf*(1-maf)), turn into average and incorporate into full equation
            b <- ( 1 - b / m_loci - mean_kinship ) / ( 1 - mean_kinship )

            # same with the inbreeding estimates, which need the previous b
            # note there's a conversion from self-kinship to inbreeding (2x-1)
            ## inbr <- 2 * ( inbr / inbr_m - b ) / ( 1 - b ) - 1
            inbr <- ( 2 * inbr / inbr_m - 1 - b ) / ( 1 - b )
            
            # prevent unneeded recalculations
            b_done <- TRUE
        }
        
        # complete final step in KP calculation
        # normalize
        KP <- KP / m_loci
        # and subtract b * colSums( P ), using sweep so it's along the rows
        KP <- sweep( KP, 2, b * colSums( P ), '-')
        # I don't renormalize by (1-b) as it doesn't matter, but it saves time to skip it (on huge matrices)
        # renormalize for troubleshooting only
        KP <- KP / ( 1 - b )
        
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
    return(
        list(
            Z = Z,
            inbr = inbr
        )
    )
}
