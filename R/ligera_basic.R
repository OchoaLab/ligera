# the original version without BEDMatrix or missingness support
# (for tests; do not export!)
# unlike main function, this one requires explicit kinship inverse to be provided (for extra sanity checks)
ligera_basic <- function(X, trait, kinship, kinship_inv, loci_on_cols = FALSE) {
    # some internal constants (preserving old tests)
    hetz <- TRUE
    hetz_indiv_inbr <- TRUE
    
    # override this for BEDMatrix
    if (class(X) == 'BEDMatrix') {
        loci_on_cols <- TRUE
    } else if (!is.matrix(X))
        stop('X has unsupported class: ', class(X))
    
    # need these dimensions
    if (loci_on_cols) {
        m_loci <- ncol(X)
        n_ind <- nrow(X)
    } else {
        m_loci <- nrow(X)
        n_ind <- ncol(X)
    }

    ##############################
    ### EFFECT SIZE ESTIMATION ###
    ##############################
    
    # precompute quantities shared across loci
    PhiInv11 <- sum( kinship_inv )
    PhiInv1y <- sum( drop( kinship_inv %*% trait ) )
    PhiInvyy <- drop( trait %*% kinship_inv %*% trait )
    # a denominator that recurs
    denom <- PhiInvyy * PhiInv11 - PhiInv1y^2

    # the projection vector
    proj <- ( PhiInv11 * drop(trait %*% kinship_inv) - PhiInv1y * rowSums(kinship_inv) ) / denom
    # the coefficients are simply the genotypes projected!
    beta <- drop( X %*% proj )
    
    ###########################
    ### VARIANCE ESTIMATION ###
    ###########################

    if (hetz) {
        # instead of estimating p_anc_hat first, here we estimate p*q per individual (that's what counting heterozygotes does), then average
        if ( hetz_indiv_inbr ) {
            # correct for inbreeding bias on a per-individual basis!
            # formulate as using weights (but these don't sum to one)
            weights_inbr <- 1 / (1 - popkin::inbr(kinship)) / 2 / n_ind

            # the desired "average"
            p_q <- drop( ( X == 1 ) %*% weights_inbr )
            
        } else {
            # bias in this case is given by FST (we correct), not mean_kinship (as for p_anc_hat version below):
            p_q <- rowMeans( X == 1 ) / 2 / ( 1 - popkin::fst(kinship) )
        }

        # in all hetz cases we need to regularize, as zero estimates are possible (likely even on the genome-wide level) and ruin inference completely
        # this is the crudest "laplace prior"-like version that prevents zeroes and is more conservative at rarer loci (which is good)
        # un-average the previous p_q by multiplying it by the sample size n_ind, then apply the correction and renormalize again.
        p_q <- ( 1 + n_ind * p_q ) / ( 2 + n_ind )
        
    } else {
        
        # to correct for a variance bias
        # will assume uniform weights!
        mean_kinship <- popkin::mean_kinship(kinship)
        
        # compute all p_anc_hat, vectorizing
        p_anc_hat <- rowMeans( X, na.rm = TRUE ) / 2
        
        # construct final estimate of p*q
        # includes (1 - mean_kinship) bias correction for this particular estimation approach
        p_q <- p_anc_hat * ( 1 - p_anc_hat ) / (1 - mean_kinship)
    }

    # construct final variance estimate of beta
    beta_std_dev <- sqrt( 4 * p_q * PhiInv11 / denom )
    
    ####################
    ### T-STATISTICS ###
    ####################
    
    # the test t-statistics
    t_stat <- beta / beta_std_dev
    
    # replace infinities with NAs
    t_stat[ is.infinite(t_stat) ] <- NA
    
    ################
    ### P-VALUES ###
    ################
    
    # Get naive p-values assuming a null t-distribution with the obvious degrees of freedom
    # so far I know this is wrong (the true tails are longer), may need to adjust the degrees of freedom (not yet known how exactly)
    # NOTE: two-sided test!
    pval <- 2 * stats::pt(-abs(t_stat), n_ind-1)

    # done, return quantities of interest (nice table!)
    return(
        tibble::tibble(
            pval = pval,
            beta = beta,
            beta_std_dev = beta_std_dev,
            p_q = p_q,
            t_stat = t_stat
        )
    )
}
