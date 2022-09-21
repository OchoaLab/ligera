# this package, to perform the genetic association study
library(ligera)
# a package that reads the BIM (variant annotations) and FAM (individual annotations) files, and reads and writes other genetics files
library(genio)
# for manipulating data.frames/tibbles
suppressMessages( library(dplyr) )
# to write the final report to a file
library(readr)
# for terminal options
library(optparse)

# hardcoded option for _f versions
V <- 1

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--pheno", type = "character", default = NA, 
                help = "Base name for phenotype file (.PHEN).  If missing, phenotype in FAM table is used", metavar = "character"),
    make_option("--mean_kinship", type = "double", default = NA, 
                help = "Mean kinship value", metavar = "double"),
    make_option("--multi", action = "store_true", default = FALSE, 
                help = "Use multiscan version (slower but more powerful)"),
    make_option("--fstat", action = "store_true", default = FALSE, 
                help = "Use f-statistic version (slower but more calibrated)"),
    make_option("--out", type = "character", default = NA, 
                help = "Base name for report output file (default same as --bfile)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
name_phen <- opt$pheno
name_out <- opt$out
mean_kinship <- opt$mean_kinship
multi <- opt$multi
fstat <- opt$fstat

# check for required values
if ( is.na( name ) )
    stop('Base name (--bfile) for input plink files is required!')
# if out is missing, set it to same as --file
if ( is.na( name_out ) )
    namme_out <- name
# if name_phen is missing, we read phenotype from FAM (see below)

# always need fam for phen matching
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

message('ligera...')
# now we have all the parts we need, run LIGERA!
if (multi) {
    if (fstat) {
        # this version reads bim/fam internally, outputs bim merged as below already
        tib <- ligera2_bed_f_multi(
            file = name,
            trait = trait,
            mean_kinship = mean_kinship,
            prune_ld = TRUE, # force this feature for real datasets
            V = V
        )
    } else {
        # this version reads bim/fam internally, outputs bim merged as below already
        tib <- ligera2_bed_multi(
            file = name,
            trait = trait,
            mean_kinship = mean_kinship,
            prune_ld = TRUE # force this feature for real datasets
        )
    }
} else {
    # read the full BIM table
    bim <- read_bim( name )
    # actual processing
    if (fstat) {
        tib <- ligera2_bed_f(
            file = name,
            m_loci = nrow( bim ),
            n_ind = nrow( fam ),
            trait = trait,
            mean_kinship = mean_kinship,
            V = V
        )
    } else {
        tib <- ligera2_bed(
            file = name,
            m_loci = nrow( bim ),
            n_ind = nrow( fam ),
            trait = trait,
            mean_kinship = mean_kinship
        )
    }
    # can merge BIM into table to have a more complete report
    tib <- bind_cols( bim, tib )
}
# write full table into an output file
file_out <- paste0( name_out, '.txt' )
write_tsv( tib, file_out )
