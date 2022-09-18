# calculates KP, and optionally b and inbr.
# KP uses b, so it must always be calculated (it is either calculated the first time or provided from an earlier iteration).
# inbr is not used except externally, and only an older version uses it, so we can skip that here.
# 
# @return If b = NA, a list with: KP, b, and optionally inbr.  Otherwise a list with KP only.
popkin_prod <- function (
                         X,
                         P,
                         mean_kinship = NA,
                         b = NA,
                         indexes_ind = NULL,
                         n_data_cut = 1e6, # block size limit for BEDMatrix
                         loci_on_cols = FALSE,
                         want_inbr = TRUE # so old code works
                         ) {
    # validations
    if ( missing( X ) )
        stop( '`X` genotype matrix is required!' )
    if ( missing( P ) )
        stop( '`P` matrix is required!' )

    # validate X and P
    if ( !is.matrix(X) ) # TRUE for BEDMatrix
        stop( '`X` must be a matrix!' )
    if ( !is.matrix(P) )
        stop( '`P` must be a matrix!' )
    # force loci on columns here
    if ( 'BEDMatrix' %in% class( X ) )
        loci_on_cols <- TRUE

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
    if ( nrow( P ) != n_ind )
        stop( 'Number of rows of `P` (', nrow( P ), ') must match the number of individuals in the (possibly filtered) genotype matrix (', n_ind, ')!' )
    k_covars <- ncol( P )
    
    # the corresponding multiplied version, approximately proportional to: kinship %*% P
    # will be a running average, so initialize to zeroes (same dim as P)
    KP <- matrix(
        0,
        nrow = n_ind,
        ncol = k_covars
    )

    # if b wasn't provided, we take that to mean that it has to be computed
    want_b <- is.na( b )
    
    if ( want_b ) {
        # need mean_kinship exclusively for b
        # validate mean_kinship, must be numeric scalar
        if ( length( mean_kinship ) != 1 )
            stop( '`mean_kinship` (required when `b` is NA) must be scalar!  Got length "', length( mean_kinship ), '"' )
        if ( is.na( mean_kinship ) )
            stop( '`mean_kinship` is required when `b` is NA!' )
        if ( !is.numeric( mean_kinship ) )
            stop( '`mean_kinship` (required when `b` is NA) must be numeric!' )
        
        # compute b (scalar), the minimum B = crossprod(X-1)/m
        # since we're not forming B, we estimate b indirectly, as
        # b = ( 1 - 4 * mean( maf * (1 - maf) ) - mean_kinship ) / ( 1 - mean_kinship )
        # have to calculate it in first pass of genotypes scan, as a running sum
        b <- 0
        if ( want_inbr ) {
            # also need to estimate the inbreeding coefficient vector!
            inbr <- vector( 'numeric', n_ind )
            inbr_m <- vector( 'numeric', n_ind ) # normalization (number of non-NA loci per individual)
        }
    }
    
    # figure out how many loci to read at the time
    m_chunk <- floor( n_data_cut / n_ind )
    if ( m_chunk < 1 )
        m_chunk <- 1
    
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
        if ( want_b ) {
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

        if ( want_b && want_inbr ) {
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
    if ( want_b ) {
        # previous b was just sum(maf*(1-maf)), turn into average and incorporate into full equation
        b <- ( 1 - b / m_loci - mean_kinship ) / ( 1 - mean_kinship )

        if ( want_inbr ) {
            # same with the inbreeding estimates, which need the previous b
            # note there's a conversion from self-kinship to inbreeding (2x-1)
            ## inbr <- 2 * ( inbr / inbr_m - b ) / ( 1 - b ) - 1
            inbr <- ( 2 * inbr / inbr_m - 1 - b ) / ( 1 - b )
        }
    }
    
    # complete final step in KP calculation
    # normalize
    KP <- KP / m_loci
    # and subtract b * colSums( P ) along the rows
    KP <- KP - matrix(
                   b * colSums( P ),
                   nrow = n_ind,
                   ncol = k_covars,
                   byrow = TRUE
               )
    # I thought renormalizing by (1-b) shouldn't matter, and it saves time to skip it (on huge matrices)
    # renormalize for unit tests though, actually it does matter elsewhere, surprisingly
    KP <- KP / ( 1 - b )

    # return everything we computed
    # KP is always calculated
    obj <- list( KP = KP )
    # these are optional, add only if we made them at all
    if ( want_b ) {
        obj$b <- b
        if ( want_inbr )
            obj$inbr <- inbr
    }
    return( obj )
}
