# this version just assumes that b and inbr are computed outside, never computes them in here

popkin_prod_bed <- function (
                             X,
                             P,
                             b,
                             indexes_ind = NULL,
                             n_data_cut = 1e6, # block size limit for BEDMatrix
                             loci_on_cols = FALSE
                             ) {
    # validations
    if ( missing( X ) )
        stop( '`X` genotype matrix is required!' )
    if ( missing( P ) )
        stop( '`P` matrix is required!' )
    if ( missing( b ) )
        stop( '`b` scalar is required!' )

    # validate X, P, b
    if ( !is.matrix(X) ) # TRUE for BEDMatrix
        stop( '`X` must be a matrix!' )
    if ( !is.matrix(P) )
        stop( '`P` must be a matrix!' )
    if ( length( b ) != 1 )
        stop( '`b` must be scalar!  Got length ', length( b ) )
    if ( !is.numeric(b) )
        stop( '`b` must be numeric!' )
    
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
        
        # manipulate for kinship-like estimation
        # "center" with -1
        X_chunk <- X_chunk - 1

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

    # complete final step in KP calculation
    # normalize
    KP <- KP / m_loci
    # replace sweep with matrix method
    KP <- KP - matrix((b * colSums( P )), dim(P)[1], length(b * colSums( P )), byrow = TRUE)
    # I thought renormalizing by (1-b) shouldn't matter, and it saves time to skip it (on huge matrices)
    # renormalize for unit tests though, actually it does matter elsewhere, surprisingly
    KP <- KP / ( 1 - b )

    # return KP only!
    return( KP ) 
}
