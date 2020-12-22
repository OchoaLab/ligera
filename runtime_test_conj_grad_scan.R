# load packages
library(ligera)
library(BEDMatrix)
library(genio)
suppressMessages( library(dplyr) )
library(readr)
library(optparse)
library(Rcpp)
library(profvis)

# load necessary variables
name = "sample"
name_phen = "sample"
name_out = "sample_ligera2"
mean_kinship = 0.15

# remember to set working directory
# setwd("~/Documents/fall 2020/rotation1/ligera/scripts")
X <- BEDMatrix( name )
# read the full BIM and FAM tables too
bim <- read_bim( name )
fam <- read_fam( name )
# need dimensions for this version of the analysys
n_ind <- nrow( fam )
m_loci <- nrow( bim )

if ( !is.na( name_phen ) ) {
  # if the phenotype is in a separate PHEN file, load that
  phen <- read_phen( name_phen )
  # reorder phen to match fam, if needed
  phen <- phen[ match( fam$id, phen$id ), ]
  # sanity check
  # when phen has fewer individuals than fam, some are NAs after match above, so this takes care of those cases
  stopifnot( all(phen$id == fam$id, na.rm = TRUE) )
  # save trait now
  trait <- phen$pheno
} else {
  # if the trait of interest is in the fam file, extract it now
  trait <- fam$pheno
}

# tib <- ligera2_bed(
#   file = name,
#   m_loci = m_loci,
#   n_ind = n_ind,
#   trait = trait,
#   mean_kinship = mean_kinship)

####################################
Y <- cbind( trait, 1 )
tol = 1e-15
indexes_ind = NULL

# load cpp functions to the environment
sourceCpp("/Users/tiffanytu/Documents/fall 2020/rotation1/ligera/src/get_b_inbr_bed_cpp.cpp")
sourceCpp("/Users/tiffanytu/Documents/fall 2020/rotation1/ligera/src/popkin_prod_bed_cpp.cpp")

file = "sample.bed"

obj <- get_b_inbr_bed_cpp(
  file,
  m_loci,
  n_ind,
  mean_kinship,
  indexes_ind
)

b <- obj$b
inbr <- obj$inbr

# dissecting conj_grad_scan_bed_wcpp function

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
Z <- matrix(
  0,
  nrow = n_ind_kept,
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

profvis({
  # start loop
  while ( any( rep.int( TRUE, k_covars )) ) {
    # P and R matrices are always non-converged subsxets!
    
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
    #Z_sweep = sweep( P, 2, alpha, '*')
    #Z_trans = t( t(P) * alpha )
    Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
    Z[ , not_converged ] <- Z[ , not_converged ] + Z_matrix
    #R <- R - sweep( KP, 2, alpha, '*')
    #R <- R - t( t(KP) * alpha )
    R <- R - KP * matrix(alpha, dim(KP)[1], length(alpha), byrow = TRUE)
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
    #P <- R + sweep( P, 2, Rn1 / Rn, '*')
    #P <- R + t( t(P) * (Rn1 / Rn) )
    P <- R + P * matrix((Rn1/Rn), dim(P)[1], length((Rn1/Rn)), byrow = TRUE)
    Rn <- Rn1
  }
  
  # after everything has converged, return the matrix of interest!
  return( Z )
})

htmlwidgets::saveWidget(p, "profile.html")
profvis({conj_grad_scan_bed_wcpp_matrix()})

Z_sweep = Z
Z_matrix = Z
identical(Z_sweep, Z_matrix)

Z_rep = P * rep(alpha, each = nrow(P))
identical(Z_rep, Z_sweep)
Z_rep2 = P * rep(alpha, rep.int(nrow(P), length(alpha)))
identical(Z_rep, Z_rep2)
############ second method ##############
    # P and R matrices are always non-converged subsets!
    
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
    P_matrix <- P
    
    library(microbenchmark)
    microbenchmark(
    Z_sweep = sweep( P, 2, alpha, '*'),
    Z_trans = t( t(P) * alpha ),
    #Z_apply = t(apply(P,1 , '*', alpha)),
    Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
    #Z_rep = P * rep(alpha, each = nrow(P)),
    #Z_rep2 = P * rep(alpha, rep.int(nrow(P), length(alpha))),
    #Z_col = P * alpha[col(P)]
    )
    
    profvis({
        Z_sweep = sweep( P, 2, alpha, '*')
        # 2 - col wise/ 1 - row wise
        Z_trans = t( t(P) * alpha )
        Z_apply = t(apply(P,1 , '*', alpha))
        Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
    })
    
    # create matrix from vector before multiplying
   #P_matrix <- P
    Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
    identical(Z_sweep, Z_matrix)
    identical(P_matrix, P)
    Z_col = P * alpha[col(P)]
    identical(Z_sweep, Z_col)

    output = apply(P, 2, alpha)
    identical(output, Z_sweep)
#########################################
# memory
library(profmem)
profmem({
      Z_sweep = sweep( P, 2, alpha, '*')
    })

profmem({
  Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
}) 

profmem({
  Z_rep2 = P * rep(alpha, rep.int(nrow(P), length(alpha)))
}) 
#########################################
# increase dim of P
# trait initially at 100
trait_1 = rep(trait, 1000)
Y <- cbind( trait_1, 1, 1, 1, 1, 1, 1)

# get dimensions from Y
n_ind_kept <- nrow( Y )
k_covars <- ncol( Y )
Z <- matrix(
  0,
  nrow = n_ind_kept,
  ncol = k_covars
)
R <- Y
P <- R
Rn <- colSums( R^2 )
KP <- popkin_prod_bed_cpp(
  file,
  m_loci,
  n_ind, # number of individuals in BED file, not kept
  P,
  b,
  indexes_ind
)
alpha <- Rn / colSums(P * KP)
Z_sweep = sweep( P, 2, alpha, '*')

microbenchmark(
  Z_sweep = sweep( P, 2, alpha, '*'),
  Z_trans = t( t(P) * alpha ),
  #Z_apply = t(apply(P,1 , '*', alpha)),
  Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
  #Z_rep = P * rep(alpha, each = nrow(P)),
  #Z_rep2 = P * rep(alpha, rep.int(nrow(P), length(alpha))),
  #Z_col = P * alpha[col(P)]
)

profmem({
  Z_sweep = sweep( P, 2, alpha, '*')
})

profmem({
  Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
}) 

profmem({
  Z_rep2 = P * rep(alpha, rep.int(nrow(P), length(alpha)))
}) 

profmem({
  Z_trans = t( t(P) * alpha )
}) 
#########################################
# 12/15 testing
library(microbenchmark)
microbenchmark(
  Z_sweep = sweep( P, 2, alpha, '*'),
  Z_trans = t( t(P) * alpha ),
  Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE),
  times = 100L
  #setup <- {P <- runif(R); alpha <- Rn / colSums(P * KP)}
)

microbenchmark(
  Z_sweep = sweep( P, 2, alpha, '*'),
  Z_trans = t( t(P) * alpha ),
  Z_matrix = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE),
  times = 10000L,
  setup = {P <- matrix(runif(P), ncol = 2); alpha <- runif(Rn / colSums(P * KP))}
)
#########################################
# test conj_grad_scan functions
microbenchmark(
  conj_grad_scan_bed_wcpp_test(file, m_loci, n_ind, Y, b, indexes_ind = NULL, tol = 1e-15),
  conj_grad_scan_bed_wcpp_matrix_test(file, m_loci, n_ind, Y, b, indexes_ind = NULL, tol = 1e-15),
  #conj_grad_scan_bed_wcpp_trans(file, m_loci, n_ind, Y, b, indexes_ind = NULL, tol = 1e-15),
  times = 10L,
  setup = {P <- matrix(runif(P), ncol = 2); alpha <- runif(Rn / colSums(P * KP))}
)
Z1 <- Z
Z2 <- Z
microbenchmark(Z1[ , not_converged ] <- Z1[ , not_converged ] + sweep( P, 2, alpha, '*'), 
               Z2[ , not_converged ] <- Z2[ , not_converged ] + P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE),
               setup = {P <- matrix(runif(P), ncol = 2); alpha <- runif(Rn / colSums(P * KP))}
               )
microbenchmark(R <- R - sweep( KP, 2, alpha, '*'),
               R <- R - KP * matrix(alpha, dim(KP)[1], length(alpha), byrow = TRUE))
#########################################
sweep_method = function(P, alpha){
  Z = sweep( P, 2, alpha, '*')
  return(Z)
}
matrix_method = function(P, alpha){
  Z = P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
  return(Z)
}
test_sweep = sweep_method(P, alpha)
test_matrix = matrix_method(P, alpha)
identical(test_sweep, test_matrix) #TRUE
microbenchmark(sweep_method(P, alpha), matrix_method(P, alpha), times = 10L, 
               setup = {P <- matrix(runif(P), ncol = 2); 
               alpha <- runif(Rn / colSums(P * KP))})
microbenchmark(test_sweep = sweep_method(P, alpha), test_matrix = matrix_method(P, alpha), times = 10L, 
               setup = {P <- matrix(runif(P), ncol = 2); 
               alpha <- runif(Rn / colSums(P * KP))})
######################################
conj_grad_scan_bed_wcpp_test <- function(
  file,
  m_loci,
  n_ind, # number of individuals in BED file (may be less than kept individuals in Y or indexes_ind)
  Y,
  b,
  indexes_ind = NULL,
  tol = 1e-15
) {
  Z[ , not_converged ] <- Z[ , not_converged ] + sweep( P, 2, alpha, '*')
}

conj_grad_scan_bed_wcpp_matrix_test <- function(
  file,
  m_loci,
  n_ind, # number of individuals in BED file (may be less than kept individuals in Y or indexes_ind)
  Y,
  b,
  indexes_ind = NULL,
  tol = 1e-15
) {
  Z[ , not_converged ] <- Z[ , not_converged ] + P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
}

conj_grad_scan_bed_wcpp_test <- function(
  file,
  m_loci,
  n_ind, # number of individuals in BED file (may be less than kept individuals in Y or indexes_ind)
  Y,
  b,
  indexes_ind = NULL,
  tol = 1e-15
) {
  #validations
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
  Z <- matrix(
    0,
    nrow = n_ind_kept,
    ncol = k_covars
  )
  # # other internal matrices, which get updated as we go (columns drop as we converge)
  # # residuals, initial values, updating as we progress
  # # R <- Y - kinship %*% Z
  R <- Y
  # # conjugate vectors (initial values)
  P <- R
  # # residual norms
  Rn <- colSums( R^2 )
  # # different covariate columns may converge at different times, let's keep track of that
  # not_converged <- rep.int( TRUE, k_covars )

  # start loop
  # #while ( any( not_converged ) ) {
  # # P and R matrices are always non-converged subsets!
  # 
  # # NOTE: this is the slowest part!
  KP <- popkin_prod_bed_cpp(
    file,
    m_loci,
    n_ind, # number of individuals in BED file, not kept
    P,
    b,
    indexes_ind
  )
  # 
  # # now continue to update variables in the CG iterations, vectorized if possible
  # 
  # # another vector of the same length
   alpha <- Rn / colSums(P * KP)
  # sweep makes alpha multiply every row of P, KP (normal product is by columns)
  #sweep <- sweep( P, 2, alpha, '*')
  Z[ , not_converged ] <- Z[ , not_converged ] + sweep( P, 2, alpha, '*')
  R <- R - sweep( KP, 2, alpha, '*')
  # new residuals vector
   Rn1 <- colSums( R^2 )
  # # take action if something has converged!
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
  # # if there are still unconverged things, keep updating things
  # # have to "sweep" the `beta = Rn1 / Rn` too, to go across rows instead of columns
   P <- R + sweep( P, 2, Rn1 / Rn, '*')
   Rn <- Rn1
  # }
  
  # after everything has converged, return the matrix of interest!
  #message('sweep')
  return(Z)
}
################
conj_grad_scan_bed_wcpp_matrix_test <- function(
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
  Z <- matrix(
    0,
    nrow = n_ind_kept,
    ncol = k_covars
  )
  # # # other internal matrices, which get updated as we go (columns drop as we converge)
  # # # residuals, initial values, updating as we progress
  # # R <- Y - kinship %*% Z
  R <- Y
  # # conjugate vectors (initial values)
  P <- R
  # # residual norms
  Rn <- colSums( R^2 )
  # # different covariate columns may converge at different times, let's keep track of that
  # not_converged <- rep.int( TRUE, k_covars )

  # start loop
  # #while ( any( not_converged ) ) {
  # # P and R matrices are always non-converged subsets!
  # 
  # # NOTE: this is the slowest part!
  KP <- popkin_prod_bed_cpp(
    file,
    m_loci,
    n_ind, # number of individuals in BED file, not kept
    P,
    b,
    indexes_ind
  )
  # 
  # # now continue to update variables in the CG iterations, vectorized if possible
  # 
  # # another vector of the same length
   alpha <- Rn / colSums(P * KP)
  # sweep makes alpha multiply every row of P, KP (normal product is by columns)
  #alpha_matrix <- matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
  Z[ , not_converged ] <- Z[ , not_converged ] + P * matrix(alpha, dim(P)[1], length(alpha), byrow = TRUE)
  R <- R - KP * matrix(alpha, dim(KP)[1], length(alpha), byrow = TRUE)
  
  # new residuals vector
   Rn1 <- colSums( R^2 )
  # # take action if something has converged!
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
  # # if there are still unconverged things, keep updating things
  # # have to "sweep" the `beta = Rn1 / Rn` too, to go across rows instead of columns
   P <- R + P * matrix((Rn1/Rn), dim(P)[1], length((Rn1/Rn)), byrow = TRUE)
   Rn <- Rn1
  #}
  
  # after everything has converged, return the matrix of interest!
  #message('matrix')
  return(Z)
}

microbenchmark(
  conj_grad_scan_bed_wcpp_test(file, m_loci, n_ind, Y, b, indexes_ind = NULL, tol = 1e-15),
  conj_grad_scan_bed_wcpp_matrix_test(file, m_loci, n_ind, Y, b, indexes_ind = NULL, tol = 1e-15),
  #conj_grad_scan_bed_wcpp_trans(file, m_loci, n_ind, Y, b, indexes_ind = NULL, tol = 1e-15),
  times = 100L,
  setup = {P <- matrix(runif(P), ncol = 2); alpha <- runif(Rn / colSums(P * KP))}
)

#############################
