# this package, to perform the genetic association study
library(ligera)
# the package that estimates the population kinship matrix
library(popkin)
# a package that reads genotype matrices using low memory
library(BEDMatrix)
# a package that reads the BIM (variant annotations) and FAM (individual annotations) files, and reads and writes other genetics files
library(genio)
# for manipulating data.frames/tibbles
suppressMessages( library(dplyr) )
# to write the final report to a file
library(readr)
# for terminal options
library(optparse)

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--pheno", type = "character", default = NA, 
                help = "Base name for phenotype file (.PHEN).  If missing, phenotype in FAM table is used", metavar = "character"),
    make_option("--multi", action = "store_true", default = FALSE, 
                help = "Use multiscan version (slower but more powerful)"),
    make_option("--out", type = "character", default = NA, 
                help = "Base name for report output file (default same as --bfile)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
name_phen <- opt$pheno
name_out <- opt$out
multi <- opt$multi

# check for required values
if ( is.na( name ) )
    stop('Base name (--bfile) for input plink files is required!')
# if out is missing, set it to same as --file
if ( is.na( name_out ) )
    namme_out <- name
# if name_phen is missing, we read phenotype from FAM (see below)

message('BEDMatrix...')
# create a BEDMatrix object, which allows random access to the genotypes
X <- BEDMatrix( name )
# read the full BIM and FAM tables too
bim <- read_bim( name )
fam <- read_fam( name )

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

message('popkin...')
# to estimate the kinship matrix, it's best to group individuals into subpopulations (so the minimum value can be estimated through averaging).
# here we assume that these labels are in the first column of the FAM file
labs <- fam$fam
# estimate the kinship matrix
# (if there are no useful subpopulation labels, `labs` can be omitted)
kinship <- popkin( X, labs )
# NOTE: kinship estimates should always be inspected for correctness
# (visualize using `plot_popkin` from package `popkin`)
# to inspect later, save kinship matrix into an output file
write_grm( name_out, kinship, fam = fam )

message('ligera...')
# now we have all the parts we need, run LIGERA!
if (multi) {
    tib <- ligera_multi( X, trait, kinship )
} else {
    tib <- ligera( X, trait, kinship )
}
# can merge BIM into table to have a more complete report
tib <- bind_cols( bim, tib )
# write full table into an output file
file_out <- paste0( name_out, '.txt' )
write_tsv( tib, file_out )
