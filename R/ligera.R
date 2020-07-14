#' LIGERA: LIght GEnetic Robust Association main function
#'
#' This function performs the genetic association tests on every locus of a genotype matrix against a quantitative trait, given a precomputed kinship matrix.
#' The function returns a tibble containing association statistics and several intermediates.
#' 
#' Suppose there are `n` individuals and `m` loci.
#'
#' @param X The `m`-by-`n` genotype matrix, containing dosage values in (0, 1, 2, NA) for the reference allele at each locus.
#' @param trait The length-`n` trait vector, which may be real valued and contain missing values.
#' @param kinship The `n`-by-`n` kinship matrix, estimated by other methods (i.e. the `popkin` package).
#' @param kinship_inv The optional matrix inverse of the kinship matrix.  Setting this parameter is not recommended, as internally a conjugate gradient method (`\link[cPCG]{cgsolve}`) is used to implicitly invert this matrix, which is much faster.  However, for very large numbers of traits without missingness and the same kinship matrix, inverting once might be faster.
#' @param covar An optional `n`-by-`K` matrix of `K` covariates, aligned with the individuals.
#' @param loci_on_cols If `TRUE`, `X` has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns.
#' If `X` is a BEDMatrix object, `loci_on_cols = TRUE` is set automatically.
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param tol Tolerance value passed to `\link[cPCG]{cgsolve}`.
#' @param maxIter Maximum number of iterations passed to `\link[cPCG]{cgsolve}`.
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
#' # Construct toy data
#' # genotype matrix
#' X <- matrix(
#'     c(0, 1, 2,
#'       1, 0, 1,
#'       1, 0, 2),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' trait <- 1 : 3
#' kinship <- diag( 3 ) / 2 # unstructured case
#'
#' tib <- ligera( X, trait, kinship )
#' tib
#'
#' @seealso
#' The `popkin` and `cPCG` packages.
#' 
#' @export
ligera <- function(
                   X,
                   trait,
                   kinship,
                   kinship_inv = NULL,
                   covar = NULL,
                   loci_on_cols = FALSE,
                   mem_factor = 0.7,
                   mem_lim = NA,
                   # cgsolve options
                   tol = 1e-15, # default 1e-6
                   maxIter = 1e6 # default 1e3
                   ) {
    # - supports missingness in trait (exact kinship matrix inverse in those cases)
    # TODO
    # - support true missingness in genotypes (inversion of matrix subsets, etc; right now only approximate)

    # some internal constants (preserving old tests)
    hetz <- TRUE
    hetz_indiv_inbr <- TRUE
    
    # informative errors when things are missing
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( trait ) )
        stop( '`trait` is required!' )
    if ( missing( kinship ) )
        stop( '`kinship` is required!' )
    # function from popkin validates further (includes square matrix test)
    popkin::validate_kinship( kinship )

    # and even further, as unexpected NAs are a pain
    if ( anyNA( kinship ) )
        stop( '`kinship` must not have any missing values!' )
    if ( !is.null( kinship_inv ) && anyNA( kinship_inv ) )
        stop( '`kinship_inv` must not have any missing values!' )
    # NOTE: trait and X may have missing values
    
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

    # check dimensions of other items
    if ( length( trait ) != n_ind )
        stop('Number of individuals in `trait` (', length( trait ), ') does not match genotype matrix (', n_ind , ')')
    if ( nrow( kinship ) != n_ind )
        stop('Number of individuals in `kinship` (', nrow( kinship ), ') does not match genotype matrix (', n_ind , ')')
    if ( !is.null( kinship_inv ) && nrow( kinship_inv ) != n_ind )
        stop('Number of individuals in `kinship_inv` (', nrow( kinship_inv ), ') does not match genotype matrix (', n_ind , ')')
    if ( !is.null( covar ) ) {
        if ( nrow( covar ) != n_ind )
            stop('Number of individuals in `covar` (', nrow( covar ), ') does not match genotype matrix (', n_ind , ')')
    }
    
    # update kinship, etc, if the trait has missing values
    # the good thing is that this is shared across loci, so comparably it's not so expensive
    # this NULL means there are no filters to apply
    indexes_ind <- NULL
    if ( anyNA( trait ) ) {
        # indexes to keep (need to subset genotypes at load time)
        indexes_ind <- !is.na( trait )
        # subset trait
        trait <- trait[ indexes_ind ]
        # subset kinship matrix
        kinship <- kinship[ indexes_ind, indexes_ind ]
        # force recomputing inverse of kinship matrix (see further below)
        kinship_inv <- NULL
        # subset covariates, if present
        if ( !is.null( covar ) )
            covar <- covar[ indexes_ind, ]
        # reduce number of individuals, used in some calculations
        n_ind <- length( trait )
        # NOTE: only genotypes are left to filter with indexes_ind
    }

    # gather matrix of trait, intercept, and optional covariates
    Y <- cbind( trait, 1 )
    # add covariates, if present
    if ( !is.null( covar ) )
        Y <- cbind( Y, covar )
    # compute inverse if needed
    if ( is.null( kinship_inv ) ) {
        # initialize matrix to fill
        Z <- matrix( NA, nrow = nrow(Y), ncol = ncol(Y) )
        for ( k in 1:ncol(Y) ) {
            Z[ , k ] <- drop( cPCG::cgsolve( kinship, Y[ , k ], tol = tol, maxIter = maxIter ) )
        }
    } else {
        # use kinship inverse if given
        Z <- kinship_inv %*% Y
    }
    # new way to abstract the rest of these
    obj <- get_proj_denom_multi( Z, Y )
    proj <- obj$proj
    beta_var_fac <- obj$var
    
    
    ##############################
    ### EFFECT SIZE ESTIMATION ###
    ##############################

    if ( hetz ) {
        if ( hetz_indiv_inbr ) {
            # correct for inbreeding bias on a per-individual basis!
            # formulate as using weights (but these don't sum to one)
            weights_inbr <- 1 / ( 1 - popkin::inbr(kinship) ) / 2
        } else {
            # the correction term is a scalar (same for all indvidiuals)
            weights_inbr <- 1 / ( 1 - popkin::fst(kinship) ) / 2
        }
    } else {
        # to correct for a variance bias
        # will assume uniform weights!
        mean_kinship <- popkin::mean_kinship(kinship)
    }
    
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
            if (loci_on_cols) {
                Xi <- t( X[ indexes_ind, indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
            } else {
                Xi <- X[ indexes_loci_chunk, indexes_ind, drop = FALSE ]
            }
        } else {
            # keep all individuals (is this faster in that case?)
            if (loci_on_cols) {
                Xi <- t( X[ , indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
            } else {
                Xi <- X[ indexes_loci_chunk, , drop = FALSE ]
            }
        }
        
        # to have good averages, we need the number of non-NA individuals per row
        n_ind_no_NA <- rowSums( !is.na(Xi) )
        # now we can turn all NAs to zeroes (as ints, lower mem)
        Xi[ is.na(Xi) ] <- 0L
        
        # the coefficients are simply the genotypes projected!
        # adjust for the NAs? (not sure if this is reasonable or not yet)
        beta[ indexes_loci_chunk ] <- drop( Xi %*% proj ) * n_ind / n_ind_no_NA
        
        ###########################
        ### VARIANCE ESTIMATION ###
        ###########################

        if (hetz) {
            # instead of estimating p_anc_hat first, here we estimate p*q per individual (that's what counting heterozygotes does), then average
            if ( hetz_indiv_inbr ) {
                # the desired "average"
                p_q_i <- drop( ( Xi == 1 ) %*% weights_inbr ) / n_ind_no_NA
            } else {
                # bias in this case is given by FST (we correct), not mean_kinship (as for p_anc_hat version below):
                p_q_i <- rowSums( Xi == 1 ) / n_ind_no_NA * weights_inbr 
            }

            # in all hetz cases we need to regularize, as zero estimates are possible (likely even on the genome-wide level) and ruin inference completely
            # this is the crudest "laplace prior"-like version that prevents zeroes and is more conservative at rarer loci (which is good)
            # un-average the previous p_q by multiplying it by the sample size n_ind, then apply the correction and renormalize again.
            p_q_i <- ( 1 + n_ind_no_NA * p_q_i ) / ( 2 + n_ind_no_NA )
            
        } else {
            # compute all p_anc_hat, vectorizing
            p_anc_hat_i <- rowSums( Xi ) / n_ind_no_NA / 2
            
            # construct final estimate of p*q
            # includes (1 - mean_kinship) bias correction for this particular estimation approach
            p_q_i <- p_anc_hat_i * ( 1 - p_anc_hat_i ) / (1 - mean_kinship)
        }
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
