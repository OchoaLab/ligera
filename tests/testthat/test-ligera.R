library(tibble)
library(popkin) # inbr

# simulate X (to share across tests)
# create a simple matrix with random valid data
n <- 100
m <- 1000
# create ancestral allele frequencies
p_anc <- runif( m )
# create genotypes
X <- rbinom( n * m, 2, p_anc )
# turn into matrix
X <- matrix( X, nrow = m, ncol = n )

# X with missing values
miss <- 0.1
X_miss <- X # copy first
X_miss[ sample( n * m, n * m * miss ) ] <- NA

# true kinship matrix in this unstructured case is I / 2
kinship <- diag( n ) / 2
mean_kinship <- mean( kinship )

# simulate phenotype (to share across tests)
trait <- rnorm( n )
# matrix version, for code that takes covariates
Y <- cbind( trait, 1 )
# get rid of names for comparison
dimnames( Y ) <- NULL

# trait with missingness
trait_miss <- trait # copy first
# repeat until we have at least one NA, otherwise this test doesn't serve its purpose
while ( !anyNA( trait_miss ) ) 
    trait_miss[ sample( n, n * miss ) ] <- NA
indexes_ind <- !is.na( trait_miss ) # to test for ind removals
Y_ind_rm <- Y[ indexes_ind , ]

# for ligera2 (full BOLT-like trick)
# construct exact kinship matrix estimate we'll test in trick version
# NOTE: this all requires that X have no NAs
x_bar <- rowMeans( X )
b <- (1 - mean( x_bar * ( 2 - x_bar ) ) - mean_kinship ) / ( 1 - mean_kinship )
kinship_est <- ( crossprod( X - 1 ) / m - b ) / ( 1 - b ) # here we do normalize properly for a plot
inbr_est <- inbr( kinship_est ) # for a comparison
# a basic validation
expect_true( !anyNA( inbr_est ) )

# version of b and inbr_est with ind removals (but no X missingness, first)
x_bar <- rowMeans( X[ , indexes_ind ] ) # overwrite x_bar, not reused elsewhere
b_ind_rm <- mean( x_bar * ( 2 - x_bar ) ) # compute in two parts for troubleshooting
b_ind_rm <- (1 - b_ind_rm - mean_kinship ) / ( 1 - mean_kinship )
kinship_est_ind_rm <- ( crossprod( X[ , indexes_ind ] - 1 ) / m - b_ind_rm ) / ( 1 - b_ind_rm ) # here we do normalize properly for a plot
inbr_est_ind_rm <- inbr( kinship_est_ind_rm ) # for a comparison

# true product for tests
KY <- kinship_est %*% Y    
# compare to vanilla `solve` (using actual inversion, which is least scalable solution)
Z <- solve( kinship_est, Y )
# a basic validation
expect_true( !anyNA( Z ) )

# versions with ind removals
# true product for tests
KY_ind_rm <- kinship_est_ind_rm %*% Y_ind_rm
# compare to vanilla `solve` (using actual inversion, which is least scalable solution)
Z_ind_rm <- solve( kinship_est_ind_rm, Y_ind_rm )

test_that("trait missingness is as desired", {
    # I had simulated this wrong, it was terrible, it deserves its own tests
    n_ind_kept <- sum( !is.na( trait_miss ) )
    # we should have repeated drawing above, until there was at least one NA, so this is a strict inequality
    expect_true( n_ind_kept < n )
    # check that the other items above agree
    expect_equal( n_ind_kept, sum( indexes_ind ) )
    expect_equal( n_ind_kept, nrow( Y_ind_rm ) )
    expect_equal( n_ind_kept, ncol( kinship_est_ind_rm ) )
    expect_equal( n_ind_kept, nrow( kinship_est_ind_rm ) )
    expect_equal( n_ind_kept, length( inbr_est_ind_rm ) )
    expect_equal( n_ind_kept, nrow( KY_ind_rm ) )
    expect_equal( n_ind_kept, nrow( Z_ind_rm ) )
    # other basic validations
    expect_true( !anyNA( b_ind_rm ) )
    expect_true( !anyNA( kinship_est_ind_rm ) )
    expect_true( !anyNA( inbr_est_ind_rm ) )
    expect_true( !anyNA( KY_ind_rm ) )
    expect_true( !anyNA( Z_ind_rm ) )
})

# missingness versions
x_bar_miss <- rowMeans( X_miss, na.rm = TRUE )
b_miss <- (1 - mean( x_bar_miss * ( 2 - x_bar_miss ), na.rm = TRUE ) - mean_kinship ) / ( 1 - mean_kinship )
expect_true( !is.na( b_miss ) )
# this one takes more steps
# missing values are just treated as zeroes internally, if missingness is random and uninformative then we're fine
Xc_miss <- X_miss - 1
Xc_miss[ is.na( Xc_miss ) ] <- 0
kinship_est_miss <- ( crossprod( Xc_miss ) / m - b_miss ) / ( 1 - b_miss )
# inbreding doesn't come from kinship under missigness, do things much more carefully actually!
# inbr_est_miss <- inbr( kinship_est_miss ) # WRONG in this implementation
inbr_est_miss <- colMeans( ( X_miss - 1 )^2, na.rm = TRUE ) # the first pass calculation
inbr_est_miss <- ( 2 * inbr_est_miss - 1 - b_miss ) / ( 1 - b_miss ) # this completes normalization
# a basic validation
expect_true( !anyNA( inbr_est_miss ) )
# true product for tests
KY_miss <- kinship_est_miss %*% Y    
# compare to vanilla `solve` (using actual inversion, which is least scalable solution)
Z_miss <- solve( kinship_est_miss, Y )
# a basic validation
expect_true( !anyNA( Z_miss ) )


test_that("ligera stops when needed", {
    # test that there are errors when crucial data is missing
    # there are 3 mandatory arguments
    # everything missing
    expect_error( ligera() )
    # singletons (2 things missing)
    expect_error( ligera( X ) )
    expect_error( ligera( trait = trait ) )
    expect_error( ligera( kinship = kinship ) )
    # pairs (1 thing missing)
    expect_error( ligera( X, trait ) )
    expect_error( ligera( X, kinship = kinship ) )
    expect_error( ligera( trait = trait, kinship = kinship ) )

    # other validations, when all three are present but some are wrong (wrong dimensions, lack of symmetry, etc)
    # here genotype has wrong number of individuals, compared to other two
    expect_error( ligera( X[ , -1 ], trait, kinship ) )
    # here trait is wrong length
    expect_error( ligera( X, trait[ -1 ], kinship ) )
    # here kinship is wrong dim
    expect_error( ligera( X, trait, kinship[ -1, -1 ] ) )
    # here kinship is not symmetric
    expect_error( ligera( X, trait, matrix(1:n^2, nrow = n) ) )
})

test_that("ligera runs on random data without missingness, matches basic version", {
    # NOTE: use <<- to remember variable globally (outside this scope)
    expect_silent( tib1 <<- ligera( X, trait, kinship ) )
    expect_true( is_tibble( tib1 ) )
    expect_equal( names( tib1 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib1 ), m )
    expect_true( !anyNA( tib1 ) )
    # range for some things
    expect_true( all( tib1$pval >= 0 ) )
    expect_true( all( tib1$pval <= 1 ) )
    expect_true( all( tib1$beta_std_dev > 0 ) )
    expect_true( all( tib1$p_q > 0 ) )
    #    expect_true( all( tib1$p_q < 1 ) ) # because of inbreeding weights, there is no upper limit technically

    # this is the basic version (requires non-missingness, so it can only be compared here)
    tib_basic <- ligera_basic( X, trait, kinship, solve(kinship) )
    expect_equal( tib1, tib_basic )
})

test_that("ligera runs on random data with missingness in X", {
    # NOTE: use <<- to remember variable globally (outside this scope)
    expect_silent( tib2 <<- ligera( X_miss, trait, kinship ) )
    expect_true( is_tibble( tib2 ) )
    expect_equal( names( tib2 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib2 ), m )
    expect_true( !anyNA( tib2 ) )
    # range for some things
    expect_true( all( tib2$pval >= 0 ) )
    expect_true( all( tib2$pval <= 1 ) )
    expect_true( all( tib2$beta_std_dev > 0 ) )
    expect_true( all( tib2$p_q > 0 ) )
})

test_that("ligera runs on random data with missingness in trait", {
    expect_silent( tib3 <- ligera( X, trait_miss, kinship ) )
    expect_true( is_tibble( tib3 ) )
    expect_equal( names( tib3 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib3 ), m )
    expect_true( !anyNA( tib3 ) )
    # range for some things
    expect_true( all( tib3$pval >= 0 ) )
    expect_true( all( tib3$pval <= 1 ) )
    expect_true( all( tib3$beta_std_dev > 0 ) )
    expect_true( all( tib3$p_q > 0 ) )
})

test_that("ligera runs on random data with missingness in X and trait", {
    expect_silent( tib4 <- ligera( X_miss, trait_miss, kinship ) )
    expect_true( is_tibble( tib4 ) )
    expect_equal( names( tib4 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib4 ), m )
    expect_true( !anyNA( tib4 ) )
    # range for some things
    expect_true( all( tib4$pval >= 0 ) )
    expect_true( all( tib4$pval <= 1 ) )
    expect_true( all( tib4$beta_std_dev > 0 ) )
    expect_true( all( tib4$p_q > 0 ) )
})

#################################
### LIGERA2 (Full BOLT trick) ###
#################################

test_that("popkin_prod stops when needed", {
    # NOTE either mean_kinship or b must be supplied (test all combos)
    # everything missing
    expect_error( popkin_prod() )
    # singletons (2 things missing)
    expect_error( popkin_prod( X = X ) )
    expect_error( popkin_prod( P = Y ) )
    expect_error( popkin_prod( mean_kinship = mean_kinship ) )
    expect_error( popkin_prod( b = b ) )
    # pairs (1 thing missing)
    expect_error( popkin_prod( X = X, P = Y ) )
    expect_error( popkin_prod( X = X, mean_kinship = mean_kinship ) )
    expect_error( popkin_prod( P = Y, mean_kinship = mean_kinship ) )
    expect_error( popkin_prod( X = X, b = b ) )
    expect_error( popkin_prod( P = Y, b = b ) )

    # other validations
    # X and Y must be matrices
    expect_error( popkin_prod( X = 1:10, P = Y, mean_kinship = mean_kinship ) )
    expect_error( popkin_prod( X = X, P = 1:10, mean_kinship = mean_kinship ) )
    # X and Y dimensions disagree
    expect_error( popkin_prod( X = X, P = Y[ , -1 ], mean_kinship = mean_kinship ) )
    # mean kinship is not scalar
    expect_error( popkin_prod( X = X, P = Y, mean_kinship = 1:10 ) )
})

test_that("popkin_prod matches direct product", {
    # for exact comparisons, use X without missingness

    # version with b present
    obj <- popkin_prod( X = X, P = Y, b = b ) # , mean_kinship = mean_kinship
    # compare to direct calculation
    expect_equal( obj$KP, KY )

    # version with missingness, b present
    obj <- popkin_prod( X = X_miss, P = Y, b = b_miss ) # , mean_kinship = mean_kinship
    # compare to direct calculation
    expect_equal( obj$KP, KY_miss )

    # version with b missing
    obj <- popkin_prod( X = X, P = Y, mean_kinship = mean_kinship )
    # compare to direct calculations
    expect_equal( obj$KP, KY )
    expect_equal( obj$b, b )
    expect_equal( obj$inbr, inbr_est )
})

test_that("conj_grad_scan stops when needed", {
    # everything missing
    expect_error( conj_grad_scan() )
    # singletons (2 things missing)
    expect_error( conj_grad_scan( X = X ) )
    expect_error( conj_grad_scan( Y = Y ) )
    expect_error( conj_grad_scan( mean_kinship = mean_kinship ) )
    # pairs (1 thing missing)
    expect_error( conj_grad_scan( X = X, Y = Y ) )
    expect_error( conj_grad_scan( X = X, mean_kinship = mean_kinship ) )
    expect_error( conj_grad_scan( Y = Y, mean_kinship = mean_kinship ) )

    # other validations
    # X and Y must be matrices
    expect_error( conj_grad_scan( X = 1:10, Y = Y, mean_kinship = mean_kinship ) )
    expect_error( conj_grad_scan( X = X, Y = 1:10, mean_kinship = mean_kinship ) )
    # X and Y dimensions disagree
    expect_error( conj_grad_scan( X = X, Y = Y[ , -1 ], mean_kinship = mean_kinship ) )
    # mean kinship is not scalar
    expect_error( conj_grad_scan( X = X, Y = Y, mean_kinship = 1:10 ) )
})

test_that("conj_grad_scan matches solve", {
    obj <- conj_grad_scan( X = X, Y = Y, mean_kinship = mean_kinship )
    # compare to precomputed Z from solve (considered ground truth)
    expect_equal( Z, obj$Z )
    # ditto direct computation
    expect_equal( inbr_est, obj$inbr )
    
    # here we just want to know that missingness in X doesn't result in missingness in outputs
    obj <- conj_grad_scan( X = X_miss, Y = Y, mean_kinship = mean_kinship )
    # compare to precomputed values
    expect_equal( Z_miss, obj$Z )
    expect_equal( inbr_est_miss, obj$inbr )

    # now test missingness in trait! (individuals removed)
    obj <- conj_grad_scan(
        X = X,
        Y = Y_ind_rm,
        mean_kinship = mean_kinship,
        indexes_ind = indexes_ind
    )
    # compare to precomputed values
    expect_equal( Z_ind_rm, obj$Z )
    expect_equal( inbr_est_ind_rm, obj$inbr ) # TODO
})

test_that("ligera2 stops when needed", {
    # test that there are errors when crucial data is missing
    # there are 3 mandatory arguments
    # everything missing
    expect_error( ligera2() )
    # singletons (2 things missing)
    expect_error( ligera2( X = X ) )
    expect_error( ligera2( trait = trait ) )
    expect_error( ligera2( mean_kinship = mean_kinship ) )
    # pairs (1 thing missing)
    expect_error( ligera2( X = X, trait = trait ) )
    expect_error( ligera2( X = X, mean_kinship = mean_kinship ) )
    expect_error( ligera2( trait = trait, mean_kinship = mean_kinship ) )

    # other validations, when all three are present but some are wrong (wrong dimensions, etc)
    # X must be matrix
    expect_error( ligera( X = 1:10, trait = trait, mean_kinship = mean_kinship ) )
    # here genotype has wrong number of individuals, compared to other two
    expect_error( ligera( X = X[ , -1 ], trait = trait, mean_kinship = mean_kinship ) )
    # here trait is wrong length
    expect_error( ligera( X = X, trait = trait[ -1 ], mean_kinship = mean_kinship ) )
    # mean_kinship is not scalar
    expect_error( ligera( X = X, trait = trait, mean_kinship = 1:10 ) )
})

test_that("ligera2 runs on random data without missingness, matches basic version", {
    # NOTE: use <<- to remember variable globally (outside this scope)
    expect_silent(
        tib5 <<- ligera2(
            X = X,
            trait = trait,
            mean_kinship = mean_kinship
        )
    )
    expect_true( is_tibble( tib5 ) )
    expect_equal( names( tib5 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib5 ), m )
    expect_true( !anyNA( tib5 ) )
    # range for some things
    expect_true( all( tib5$pval >= 0 ) )
    expect_true( all( tib5$pval <= 1 ) )
    expect_true( all( tib5$beta_std_dev > 0 ) )
    expect_true( all( tib5$p_q > 0 ) )
    #    expect_true( all( tib1$p_q < 1 ) ) # because of inbreeding weights, there is no upper limit technically

    # this is the basic version (requires non-missingness, so it can only be compared here)
    tib5_basic <- ligera_basic( X, trait, kinship_est, solve(kinship_est) )
    expect_equal( tib5, tib5_basic )
})

test_that("ligera2 runs on random data with missingness in X", {
    # NOTE: use <<- to remember variable globally (outside this scope)
    expect_silent( tib6 <<- ligera2( X = X_miss, trait = trait, mean_kinship = mean_kinship ) )
    expect_true( is_tibble( tib6 ) )
    expect_equal( names( tib6 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib6 ), m )
    expect_true( !anyNA( tib6 ) )
    # range for some things
    expect_true( all( tib6$pval >= 0 ) )
    expect_true( all( tib6$pval <= 1 ) )
    expect_true( all( tib6$beta_std_dev > 0 ) )
    expect_true( all( tib6$p_q > 0 ) )
})

test_that("ligera2 runs on random data with missingness in trait", {
    # NOTE: use <<- to remember variable globally (outside this scope)
    expect_silent( tib7 <<- ligera2( X = X, trait = trait_miss, mean_kinship = mean_kinship ) )
    expect_true( is_tibble( tib7 ) )
    expect_equal( names( tib7 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib7 ), m )
    expect_true( !anyNA( tib7 ) )
    # range for some things
    expect_true( all( tib7$pval >= 0 ) )
    expect_true( all( tib7$pval <= 1 ) )
    expect_true( all( tib7$beta_std_dev > 0 ) )
    expect_true( all( tib7$p_q > 0 ) )
})

test_that("ligera2 runs on random data with missingness in X and trait", {
    expect_silent( tib8 <- ligera2( X = X_miss, trait = trait_miss, mean_kinship = mean_kinship ) )
    expect_true( is_tibble( tib8 ) )
    expect_equal( names( tib8 ), c('pval', 'beta', 'beta_std_dev', 'p_q', 't_stat') )
    expect_equal( nrow( tib8 ), m )
    expect_true( !anyNA( tib8 ) )
    # range for some things
    expect_true( all( tib8$pval >= 0 ) )
    expect_true( all( tib8$pval <= 1 ) )
    expect_true( all( tib8$beta_std_dev > 0 ) )
    expect_true( all( tib8$p_q > 0 ) )
})

##############################################################
### LIGERA2_BED (Full BOLT trick, Rcpp optimized BED only) ###
##############################################################

test_that("popkin_prod_bed stops when needed", {
    # everything missing
    expect_error( popkin_prod_bed() )
    # singletons (2 things missing)
    expect_error( popkin_prod_bed( X = X ) )
    expect_error( popkin_prod_bed( P = Y ) )
    expect_error( popkin_prod_bed( b = b ) )
    # pairs (1 thing missing)
    expect_error( popkin_prod_bed( X = X, P = Y ) )
    expect_error( popkin_prod_bed( X = X, b = b ) )
    expect_error( popkin_prod_bed( P = Y, b = b ) )

    # other validations
    # X and Y must be matrices
    expect_error( popkin_prod_bed( X = 1:10, P = Y, b = b ) )
    expect_error( popkin_prod_bed( X = X, P = 1:10, b = b ) )
    # X and Y dimensions disagree
    expect_error( popkin_prod_bed( X = X, P = Y[ , -1 ], b = b ) )
    # b is not scalar
    expect_error( popkin_prod_bed( X = X, P = Y, b = 1:10 ) )
    # b is not numeric
    expect_error( popkin_prod_bed( X = X, P = Y, b = 'b' ) )
})

test_that("popkin_prod_bed matches direct product", {
    # for exact comparisons, use X without missingness

    # version with b present
    KP <- popkin_prod_bed( X = X, P = Y, b = b )
    # compare to direct calculation
    expect_equal( KP, KY )

    # version with missingness, b present
    KP <- popkin_prod_bed( X = X_miss, P = Y, b = b_miss )
    # compare to direct calculation
    expect_equal( KP, KY_miss )
})

test_that("conj_grad_scan_bed stops when needed", {
    # everything missing
    expect_error( conj_grad_scan_bed() )
    # singletons (2 things missing)
    expect_error( conj_grad_scan_bed( X = X ) )
    expect_error( conj_grad_scan_bed( Y = Y ) )
    expect_error( conj_grad_scan_bed( b = b ) )
    # pairs (1 thing missing)
    expect_error( conj_grad_scan_bed( X = X, Y = Y ) )
    expect_error( conj_grad_scan_bed( X = X, b = b ) )
    expect_error( conj_grad_scan_bed( Y = Y, b = b ) )

    # other validations
    # X and Y must be matrices
    expect_error( conj_grad_scan_bed( X = 1:10, Y = Y, b = b ) )
    expect_error( conj_grad_scan_bed( X = X, Y = 1:10, b = b ) )
    # X and Y dimensions disagree
    expect_error( conj_grad_scan_bed( X = X, Y = Y[ , -1 ], b = b ) )
    # b is not scalar
    expect_error( conj_grad_scan_bed( X = X, Y = Y, b = 1:10 ) )
    # b is not numeric
    expect_error( conj_grad_scan_bed( X = X, Y = Y, b = 'b' ) )
})

test_that("conj_grad_scan_bed matches solve", {
    # first without missingness
    Z_scan <- conj_grad_scan_bed( X = X, Y = Y, b = b )
    expect_equal( Z, Z_scan )
    
    # then with missingness
    Z_scan <- conj_grad_scan_bed( X = X_miss, Y = Y, b = b_miss )
    expect_equal( Z_miss, Z_scan )

    # now test missingness in trait! (individuals removed)
    Z_scan <- conj_grad_scan_bed(
        X = X,
        Y = Y_ind_rm,
        b = b_ind_rm,
        indexes_ind = indexes_ind
    )
    # compare to precomputed values
    expect_equal( Z_ind_rm, Z_scan )
})

# test BEDMatrix version, which also requires us to write the random data out (with genio)
if (
    suppressMessages( suppressWarnings( require(BEDMatrix) ) ) &&
    suppressMessages( suppressWarnings( require(genio) ) )
) {
    # write the random data to a temporary file
    # write the version with missingness separately
    name <- tempfile('delete-me-random-test') # output name without extensions!
    name_miss <- tempfile('delete-me-random-test-with-missingness') # output name without extensions!
    # write BED/BIM/FAM files
    write_plink( name, X, verbose = FALSE )
    write_plink( name_miss, X_miss, verbose = FALSE )
    
    # load back using BEDMatrix
    # load this way, passing dims, for speed and to avoid dimnames (complicates testing)
    X_BEDMatrix <- suppressMessages( suppressWarnings( BEDMatrix( name, n = n, p = m ) ) )
    X_miss_BEDMatrix <- suppressMessages( suppressWarnings( BEDMatrix( name_miss, n = n, p = m ) ) )
    
    ### ligera ###
    
    test_that("ligera runs correctly on BEDMatrix data, recovers R matrix outputs", {
        # this one compares to tib1
        expect_silent( tib1_BEDMatrix <- ligera( X_BEDMatrix, trait, kinship ) )
        expect_equal( tib1, tib1_BEDMatrix )
        
        # and this one compares to tib2
        expect_silent( tib2_BEDMatrix <- ligera( X_miss_BEDMatrix, trait, kinship ) )
        expect_equal( tib2, tib2_BEDMatrix )
    })

    ### ligera2 ###
    
    test_that("popkin_prod matches direct product", {
        # version with b present
        expect_silent(
            obj <- popkin_prod( X = X_BEDMatrix, P = Y, b = b )
        )
        expect_equal( KY, obj$KP )

        # version with missingness, b present
        expect_silent(
            obj <- popkin_prod( X = X_miss_BEDMatrix, P = Y, b = b_miss )
        )
        expect_equal( KY_miss, obj$KP )
    })
    
    test_that("conj_grad_scan runs correctly on BEDMatrix data, recovers R matrix outputs", {
        expect_silent(
            obj <- conj_grad_scan( X = X_BEDMatrix, Y = Y, mean_kinship = mean_kinship )
        )
        expect_equal( Z, obj$Z )
        expect_equal( inbr_est, obj$inbr )
        
        expect_silent(
            obj <- conj_grad_scan( X = X_miss_BEDMatrix, Y = Y, mean_kinship = mean_kinship )
        )
        expect_equal( Z_miss, obj$Z )
        expect_equal( inbr_est_miss, obj$inbr )

        # now test missingness in trait! (individuals removed)
        expect_silent(
            obj <- conj_grad_scan(
                X = X_BEDMatrix,
                Y = Y_ind_rm,
                mean_kinship = mean_kinship,
                indexes_ind = indexes_ind
            )
        )
        # compare to precomputed values
        expect_equal( Z_ind_rm, obj$Z )
        expect_equal( inbr_est_ind_rm, obj$inbr )
    })

    test_that("ligera2 runs correctly on BEDMatrix data, recovers R matrix outputs", {
        tib5_BEDMatrix <- ligera2( X = X_BEDMatrix, trait = trait, mean_kinship = mean_kinship )
        expect_equal( tib5, tib5_BEDMatrix )
        
        tib6_BEDMatrix <- ligera2( X = X_miss_BEDMatrix, trait = trait, mean_kinship = mean_kinship )
        expect_equal( tib6, tib6_BEDMatrix )
    })

    ### ligera2_bed ###

    file_bed <- paste0(name, '.bed')
    file_bed_miss <- paste0(name_miss, '.bed')
    
    test_that("get_b_inbr_bed works", {
        # version with no individuals removed
        # NULL = no individuals to remove
        expect_silent(
            obj <- get_b_inbr_bed_cpp( file_bed, m, n, mean_kinship, NULL )
        )
        expect_equal( class( obj ), 'list' )
        expect_equal( length( obj ), 2 )
        expect_equal( names( obj ), c('b', 'inbr') )
        expect_equal( obj$b, b )
        expect_equal( obj$inbr, inbr_est )

        # version with X missingness, no individuals removed
        # NULL = no individuals to remove
        expect_silent(
            obj <- get_b_inbr_bed_cpp( file_bed_miss, m, n, mean_kinship, NULL )
        )
        expect_equal( class( obj ), 'list' )
        expect_equal( length( obj ), 2 )
        expect_equal( names( obj ), c('b', 'inbr') )
        expect_equal( obj$b, b_miss )
        expect_equal( obj$inbr, inbr_est_miss )

        # version with individuals removed
        print( 'get_b_inbr_bed_cpp (ind_rm)...' )
        expect_silent(
            obj <- get_b_inbr_bed_cpp( file_bed, m, n, mean_kinship, indexes_ind )
        )
        expect_equal( class( obj ), 'list' )
        expect_equal( length( obj ), 2 )
        expect_equal( names( obj ), c('b', 'inbr') )
        expect_equal( obj$b, b_ind_rm ) # TODO: FAIL
        expect_equal( obj$inbr, inbr_est_ind_rm ) # TODO: FAIL
    })
    
    test_that("popkin_prod_bed matches direct product", {
        # version with b present
        expect_silent(
            KY_BM <- popkin_prod_bed( X = X_BEDMatrix, P = Y, b = b )
        )
        expect_equal( KY, KY_BM )

        # version with missingness, b present
        expect_silent(
            KY_miss_BM <- popkin_prod_bed( X = X_miss_BEDMatrix, P = Y, b = b_miss )
        )
        expect_equal( KY_miss, KY_miss_BM )

        # version with individuals removed
        expect_silent(
            KY_ind_rm_BM <- popkin_prod_bed(
                X = X_BEDMatrix,
                P = Y_ind_rm,
                b = b_ind_rm,
                indexes_ind = indexes_ind
            )
        )
        expect_equal( KY_ind_rm, KY_ind_rm_BM )
    })
    
    test_that("popkin_prod_bed_cpp matches direct product", {
        # version with b present
        # NULL = no individuals to remove
        expect_silent(
            KY_cpp <- popkin_prod_bed_cpp( file_bed, m, n, Y, b, NULL )
        )
        expect_equal( KY, KY_cpp )

        # version with missingness, b present
        # NULL = no individuals to remove
        expect_silent(
            KY_miss_cpp <- popkin_prod_bed_cpp( file_bed_miss, m, n, Y, b_miss, NULL )
        )
        expect_equal( KY_miss, KY_miss_cpp )

        # version with individuals removed
        print( 'popkin_prod_bed_cpp (ind_rm)...' )
        expect_silent(
            KY_ind_rm_cpp <- popkin_prod_bed_cpp( file_bed, m, n, Y_ind_rm, b_ind_rm, indexes_ind )
        )
        expect_equal( KY_ind_rm, KY_ind_rm_cpp )
    })
    
    test_that("conj_grad_scan_bed runs correctly on BEDMatrix", {
        # this one compares to Z
#        expect_silent(
            Z_BM <- conj_grad_scan_bed( X = X_BEDMatrix, Y = Y, b = b )
#        )
        expect_equal( Z, Z_BM )
        
        # and this one compares to Z_miss
#        expect_silent(
            Z_miss_BM <- conj_grad_scan_bed( X = X_miss_BEDMatrix, Y = Y, b = b_miss )
#        )
        expect_equal( Z_miss, Z_miss_BM )
        
        # now test missingness in trait! (individuals removed)
        print( 'conj_grad_scan_bed (ind_rm)...' )
        #        expect_silent(
        Z_ind_rm_BM <- conj_grad_scan_bed(
            X = X_BEDMatrix,
            Y = Y_ind_rm,
            b = b_ind_rm,
            indexes_ind = indexes_ind,
            verbose = TRUE
        )
        #        )
        # compare to precomputed values
        expect_equal( Z_ind_rm, Z_ind_rm_BM )
    })
    
    test_that("conj_grad_scan_bed_wcpp matches precomputed values", {
        # this one compares to Z
        expect_silent(
            Z_BM <- conj_grad_scan_bed_wcpp(
                file = file_bed,
                m_loci = m,
                n_ind = n,
                Y = Y,
                b = b
            )
        )
        expect_equal( Z, Z_BM )
        
        # and this one compares to Z_miss
        expect_silent(
            Z_miss_BM <- conj_grad_scan_bed_wcpp(
                file = file_bed_miss,
                m_loci = m,
                n_ind = n,
                Y = Y,
                b = b_miss
            )
        )
        expect_equal( Z_miss, Z_miss_BM )

        ## TODO: often not converging!
        ## # this one has individuals removed
        ## print( 'conj_grad_scan_bed_wcpp...' )
        ## #        expect_silent(
        ## Z_ind_rm_BM <- conj_grad_scan_bed_wcpp(
        ##     file = file_bed,
        ##     m_loci = m,
        ##     n_ind = n,
        ##     Y = Y_ind_rm,
        ##     b = b_ind_rm,
        ##     indexes_ind = indexes_ind,
        ##     verbose = TRUE
        ## )
        ## #        )
        ## expect_equal( Z_ind_rm, Z_ind_rm_BM )
    })

    test_that("ligera2_bed recovers R matrix outputs", {
        print('tib5_bed...')
        tib5_bed <- ligera2_bed(
            file = name,
            m_loci = m,
            n_ind = n,
            trait = trait,
            mean_kinship = mean_kinship #,
#            verbose = TRUE
        )
        expect_equal( tib5, tib5_bed )
        
        print('tib6_bed...')
        tib6_bed <- ligera2_bed(
            file = name_miss,
            m_loci = m,
            n_ind = n,
            trait = trait,
            mean_kinship = mean_kinship #,
#            verbose = TRUE
        )
        expect_equal( tib6, tib6_bed )

        # TODO: not converging :(
##         print('tib7_bed...')
##         tib7_bed <- ligera2_bed(
##             file = name,
##             m_loci = m,
##             n_ind = n,
##             trait = trait_miss,
##             mean_kinship = mean_kinship #,
## #            verbose = TRUE
##         )
##         expect_equal( tib7, tib7_bed )
    })

    # delete temporary files now
    delete_files_plink( name )
    delete_files_plink( name_miss )
}

