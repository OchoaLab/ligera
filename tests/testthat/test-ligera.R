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
trait_miss[ sample( n, n * miss ) ] <- NA

# for ligera2 (full BOLT-like trick)
# construct exact kinship matrix estimate we'll test in trick version
# NOTE: this all requires that X have no NAs
x_bar <- rowMeans( X )
b <- (1 - mean( x_bar * ( 2 - x_bar ) ) - mean_kinship ) / ( 1 - mean_kinship )
kinship_est <- ( crossprod( X - 1 ) / m - b ) / ( 1 - b ) # here we do normalize properly for a plot
inbr_est <- inbr( kinship_est ) # for a comparison 

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

test_that("conj_grad_scan matches solve, cPCG::cgsolve", {
    # for exact comparisons, use X without missingness
    
    # our method
    # NOTE: use <<- to remember variable globally (outside this scope)
    obj_scan1 <<- conj_grad_scan( X = X, Y = Y, mean_kinship = mean_kinship )
    Z <- obj_scan1$Z
    inbr_scan <- obj_scan1$inbr # estimate from here
    # compare this right away
    expect_equal( inbr_est, inbr_scan )
    # some basic validations
    expect_true( !anyNA( Z ) )
    expect_true( !anyNA( inbr_scan ) )

    # compare to vanilla `solve` (using actual inversion, which is least scalable solution)
    Zsolve <- solve( kinship_est, Y )
    # compare now!
    expect_equal( Z, Zsolve )
    
    if (
        suppressMessages( suppressWarnings( require(cPCG) ) )
    ) {
        # the other method, which uses existing kinship instead of using X (which scales more poorly)
        # need to tweak defaults so they agree to high accuracy
        tol <- 1e-16 # default 1e-6
        maxIter <- 1e6 # default 1e3
        # two things have to be solved separately (cgsolve doesn't do Y matrices)
        Zcg1 <- drop( cgsolve( kinship_est, Y[ , 1 ], tol = tol, maxIter = maxIter ) )
        Zcg2 <- drop( cgsolve( kinship_est, Y[ , 2 ], tol = tol, maxIter = maxIter ) )
        Zcg <- cbind( Zcg1, Zcg2 )
        # get rid of names for comparison
        dimnames( Zcg ) <- NULL

        # compare now!
        expect_equal( Z, Zcg )
    }
})

test_that("conj_grad_scan works with missingness in X", {
    # here we just want to know that missingness in X doesn't result in missingness in outputs
    
    # NOTE: use <<- to remember variable globally (outside this scope)
    obj_scan2 <<- conj_grad_scan( X = X_miss, Y = Y, mean_kinship = mean_kinship )
    # some basic validations
    expect_true( !anyNA( obj_scan2$Z ) )
    expect_true( !anyNA( obj_scan2$inbr ) )
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
    expect_silent( tib7 <- ligera2( X = X, trait = trait_miss, mean_kinship = mean_kinship ) )
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
    
    test_that("ligera runs correctly on BEDMatrix data, recovers R matrix outputs", {
        # this one compares to tib1
        expect_silent( tib1_BEDMatrix <- ligera( X_BEDMatrix, trait, kinship ) )
        expect_equal( tib1, tib1_BEDMatrix )
        
        # and this one compares to tib2
        expect_silent( tib2_BEDMatrix <- ligera( X_miss_BEDMatrix, trait, kinship ) )
        expect_equal( tib2, tib2_BEDMatrix )
    })
    
    test_that("conj_grad_scan runs correctly on BEDMatrix data, recovers R matrix outputs", {
        # this one compares to obj_scan1
        expect_silent(
            obj_scan1_BEDMatrix <- conj_grad_scan( X = X_BEDMatrix, Y = Y, mean_kinship = mean_kinship )
        )
        expect_equal( obj_scan1, obj_scan1_BEDMatrix )
        
        # and this one compares to obj_scan2
        expect_silent(
                obj_scan2_BEDMatrix <- conj_grad_scan( X = X_miss_BEDMatrix, Y = Y, mean_kinship = mean_kinship )
        )
        expect_equal( obj_scan2, obj_scan2_BEDMatrix )
    })

    test_that("ligera2 runs correctly on BEDMatrix data, recovers R matrix outputs", {
        tib5_BEDMatrix <- ligera2( X = X_BEDMatrix, trait = trait, mean_kinship = mean_kinship )
        expect_equal( tib5, tib5_BEDMatrix )
        
        tib6_BEDMatrix <- ligera2( X = X_miss_BEDMatrix, trait = trait, mean_kinship = mean_kinship )
        expect_equal( tib6, tib6_BEDMatrix )
    })

    # delete temporary files now
    delete_files_plink( name )
    delete_files_plink( name_miss )
}

