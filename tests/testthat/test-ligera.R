library(tibble)

# simulate X (to share across tests)
# create a simple matrix with random valid data
n <- 100
m <- 100
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

# simulate phenotype (to share across tests)
trait <- rnorm( n )

# trait with missingness
trait_miss <- trait # copy first
trait_miss[ sample( n, n * miss ) ] <- NA

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
    X_BEDMatrix <- suppressMessages( suppressWarnings( BEDMatrix( name ) ) )
    X_miss_BEDMatrix <- suppressMessages( suppressWarnings( BEDMatrix( name_miss ) ) )

    # this one compares to tib1
    expect_silent( tib1_BEDMatrix <- ligera( X_BEDMatrix, trait, kinship ) )
    expect_equal( tib1, tib1_BEDMatrix )

    # and this one compares to tib2
    expect_silent( tib2_BEDMatrix <- ligera( X_miss_BEDMatrix, trait, kinship ) )
    expect_equal( tib2, tib2_BEDMatrix )

    # delete temporary files now
    delete_files_plink( name )
    delete_files_plink( name_miss )
}
