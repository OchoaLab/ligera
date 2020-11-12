#' @useDynLib ligera
#' @importFrom Rcpp sourceCpp
NULL

#' LIGERA2 BED: LIght GEnetic Robust Association main function
#'
#' This function performs the genetic association tests on every locus of a genotype matrix against a quantitative trait, implicitly computing the kinship matrix in a way that scales better than an explicit kinship estimate.
#' The function returns a tibble containing association statistics and several intermediates.
#' This optimized version requires the genotypes to be in a file in BED format.
#' 
#' Suppose there are `n` individuals and `m` loci.
#'
#' @param file The path to the BED file containing the genotypes, potentially excluding the BED extension.
#' @param m_loci The number of loci in the BED file.
#' @param n_ind The number of individuals in the BED file.
#' @param trait The length-`n` trait vector, which may be real valued and contain missing values.
#' @param mean_kinship An estimate of the mean kinship produced externally, to ensure internal estimates of kinship and inbreeding are unbiased.
#' @param covar An optional `n`-by-`K` matrix of `K` covariates, aligned with the individuals.
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param tol Tolerance value passed to conjugate gradient method solver.
#'
#' @return A tibble containing the following association statistics
#' 
#' - `pval`: The p-value of the association test
#' - `beta`: The estimated effect size coefficient for the trait vector at this locus
#' - `beta_std_dev`: The estimated coefficient variance of this locus (varies due to dependence on minor allele frequency)
#' - `p_q`: The allele variance estimate (estimate of `p*(1-p)`).  The number of heterozygotes, weighted by inbreeding coefficient, and with pseudocounts included, is used in this estimate (in other words, it does not equal MAF * ( 1 - MAF ), where MAF is the marginal allele frequency.
#' - `t_stat`: The test statistic, equal to `beta / beta_std_dev`.
#'
#' @examples
#' # MISSING SAMPLE BED FILE
#'
#' @seealso
#' The `popkin` package.
#' 
#' @export
ligera2_bed_matrix <- function(
                        file,
                        m_loci,
                        n_ind,
                        trait,
                        mean_kinship,
                        covar = NULL,
                        mem_factor = 0.7,
                        mem_lim = 1,
                        tol = 1e-15
                        ) {
    # - supports missingness in trait
    # TODO
    # - support true missingness in genotypes (inversion of matrix subsets, etc; right now only approximate)

    # informative errors when things are missing
    if ( missing( file ) )
        stop( 'Genotype BED file `file` is required!' )
    if ( missing( m_loci ) )
        stop( 'Number of loci `m_loci` is required!' )
    if ( missing( n_ind ) )
        stop( 'Number of individuals `n_ind` is required!' )
    if ( missing( trait ) )
        stop( '`trait` is required!' )
    if ( missing( mean_kinship ) )
        stop( '`mean_kinship` is required!' )

    # validate mean_kinship, must be numeric scalar
    if ( !is.numeric( mean_kinship ) )
        stop( '`mean_kinship` must be numeric!' )
    if ( length( mean_kinship ) != 1 )
        stop( '`mean_kinship` must be scalar!  Got length "', length( mean_kinship ), '"' )

    # load genotypes using BEDMatrix, for some parts that still use it
    X <- BEDMatrix::BEDMatrix( file, n = n_ind, p = m_loci )

    # NOTE: trait and X may have missing values
    
    # check dimensions of other items
    if ( length( trait ) != n_ind )
        stop('Number of individuals in `trait` (', length( trait ), ') does not match genotype matrix (', n_ind , ')')
    if ( !is.null( covar ) ) {
        if ( nrow( covar ) != n_ind )
            stop('Number of individuals in `covar` (', nrow( covar ), ') does not match genotype matrix (', n_ind , ')')
    }
    
    # handle missing values in trait
    # the good thing is that this is shared across loci, so comparably it's not so expensive
    # this NULL means there are no filters to apply
    indexes_ind <- NULL
    # differentiate between total individuals (in BED file) and indivdiuals kept after trait-NA filters (need both values)
    n_ind_kept <- n_ind
    if ( anyNA( trait ) ) {
        # indexes to keep (need to subset genotypes at load time)
        indexes_ind <- !is.na( trait )
        # subset trait
        trait <- trait[ indexes_ind ]
        # subset covariates, if present
        if ( !is.null( covar ) )
            covar <- covar[ indexes_ind, ]
        # reduce number of individuals, used in some calculations
        n_ind_kept <- length( trait )
        # NOTE: only genotypes are left to filter with indexes_ind
    }

    # for my cruder code, add BED extension explicitly if missing
    # based on genio:::add_ext
    if ( ! grepl( '\\.bed$', file) ) {
        # add extension if it wasn't there already
        file <- paste0(file, '.bed')
    }
    
    # calculate b and inbr with this fast code
    # (this scans the whole BED file once)
    obj <- get_b_inbr_bed_cpp(
        file,
        m_loci,
        n_ind,
        mean_kinship,
        indexes_ind
    )
    b <- obj$b
    inbr <- obj$inbr
    
    # only two things have to be solved, all vectors
    Y <- cbind( trait, 1 )
    # add covariates, if present
    if ( !is.null( covar ) ) {
        # handle NAs now, so final Y has no missingness whatsoever
        covar <- covar_fix_na( covar )
        Y <- cbind( Y, covar )
    }
    Z <- conj_grad_scan_bed_wcpp_matrix(
        file = file,
        m_loci = m_loci,
        n_ind = n_ind,
        Y = Y,
        b = b,
        indexes_ind = indexes_ind,
        tol = tol
    )

    # new way to abstract the rest of these
    obj <- get_proj_denom_multi( Z, Y ) # defaults hold: trait_only = TRUE, trait_index = 1
    proj <- obj$proj
    beta_var_fac <- obj$var
    
    ##############################
    ### EFFECT SIZE ESTIMATION ###
    ##############################

    # correct for inbreeding bias on a per-individual basis!
    # formulate as using weights (but these don't sum to one)
    weights_inbr <- 1 / ( 1 - inbr ) / 2
    
    # initialize output vectors
    # Do before get_mem_lim_m so free memory is accounted for properly
    beta <- vector('numeric', m_loci)
    p_q <- vector('numeric', m_loci)
    
    # this overcounts since there are logical branches, not all overlap, but meh seriously
    # as usual, this should be conservative
    # recall that ints count as 0.5, doubles as 1
    #
    # vec_m, int
    # # indexes_loci_chunk
    # vec_m, double
    # # drop( Xi %*% proj )
    # # p_q_i
    # # p_anc_hat_i
    #
    # vec_n # int
    # # n_ind_no_NA
    #
    # mat_m_n # all ints
    # # Xi
    # # M # unnamed
    # # ( Xi == 1 )
    # add one more double copy of Xi, this happens in matrix operations, are only temporary unnamed matrices
    
    # estimating total memory usage in bytes
    data <- popkin:::solve_m_mem_lim(
                         n = n_ind,
                         m = m_loci,
                         mat_m_n = 3,
                         vec_m = 3.5,
                         vec_n = 0.5,
                         mem = mem_lim,
                         mem_factor = mem_factor
                     )
    m_chunk <- data$m_chunk

    # navigate chunks
    i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while (TRUE) { # start an infinite loop, break inside as needed
        # this means all SNPs have been covered!
        if (i_chunk > m_loci)
            break
        
        # range of SNPs to extract in this chunk
        indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m_loci)

        if ( !is.null( indexes_ind ) ) {
            # individuals get filtered here (indexes_ind; required when there's missingness in trait)
            Xi <- t( X[ indexes_ind, indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
        } else {
            # keep all individuals (is this faster in that case?)
            Xi <- t( X[ , indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
        }
        
        # to have good averages, we need the number of non-NA individuals per row
        n_ind_no_NA <- rowSums( !is.na(Xi) )
        # now we can turn all NAs to zeroes (as ints, lower mem)
        Xi[ is.na(Xi) ] <- 0L
        
        # the coefficients are simply the genotypes projected!
        # adjust for the NAs? (not sure if this is reasonable or not yet)
        beta[ indexes_loci_chunk ] <- drop( Xi %*% proj ) * n_ind_kept / n_ind_no_NA
        
        ###########################
        ### VARIANCE ESTIMATION ###
        ###########################

        # instead of estimating p_anc_hat first, here we estimate p*q per individual (that's what counting heterozygotes does), then average
        # the desired "average"
        p_q_i <- drop( ( Xi == 1 ) %*% weights_inbr ) / n_ind_no_NA
        
        # in all hetz cases we need to regularize, as zero estimates are possible (likely even on the genome-wide level) and ruin inference completely
        # this is the crudest "laplace prior"-like version that prevents zeroes and is more conservative at rarer loci (which is good)
        # un-average the previous p_q by multiplying it by the sample size n_ind, then apply the correction and renormalize again.
        p_q_i <- ( 1 + n_ind_no_NA * p_q_i ) / ( 2 + n_ind_no_NA )
        
        # transfer this vector to main one
        p_q[ indexes_loci_chunk ] <- p_q_i
        
        # update starting point for next chunk! (overshoots at the end, that's ok)
        i_chunk <- i_chunk + m_chunk
    }

    # construct final variance estimate of beta
    beta_std_dev <- sqrt( 4 * p_q * beta_var_fac )
    
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
    pval <- 2 * stats::pt( -abs( t_stat ), n_ind_kept - 1 )

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
