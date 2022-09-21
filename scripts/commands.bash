# sample commands, for running all variants on the sample data

# run ligera on sample data
time Rscript ligera.R --bfile sample --pheno sample --out sample_ligera
# 0m3.156s
time Rscript ligera.R --bfile sample --pheno sample --multi --out sample_ligera_multi
# 0m2.429s

# run ligera2 on the same sample data
# this version circumvents popkin, but we have to provide the mean kinship value instead
# here we use the true value for the simulation, but in reality we'd have to estimate it separately
# write a different output file
time Rscript ligera2.R --bfile sample --pheno sample --mean_kinship 0.15 --out sample_ligera2
# 0m5.059s 
time Rscript ligera2.R --bfile sample --pheno sample --mean_kinship 0.15 --multi --out sample_ligera2_multi
# 0m5.235s

# run ligera2_bed on the same sample data
# write yet a different output file
time Rscript ligera2_bed.R --bfile sample --pheno sample --mean_kinship 0.15 --out sample_ligera2_bed
# 0m3.645s
time Rscript ligera2_bed.R --bfile sample --pheno sample --mean_kinship 0.15 --multi --out sample_ligera2_bed_multi
# 0m3.887s

# only test with f-statistics
time Rscript ligera2_bed.R --bfile sample --pheno sample --mean_kinship 0.15 --fstat --out sample_ligera2_bed_f
# ran for too long, killed
time Rscript ligera2_bed.R --bfile sample --pheno sample --mean_kinship 0.15 --fstat --multi --out sample_ligera2_bed_f_multi
# ran for too long, killed

# the ligera and ligera2 outputs are slightly different because kinship estimates are slightly different, but estimates are highly correlated.
# the ligera2 and ligera2_bed outputs should be identical, but due to limited numerical accuracy there are small differences

# # this is R code that compares p-values
# library(readr)
# data1 <- read_tsv('sample_ligera.txt')
# data2 <- read_tsv('sample_ligera2.txt')
# data3 <- read_tsv('sample_ligera2_bed.txt')
# data1m <- read_tsv('sample_ligera_multi.txt')
# data2m <- read_tsv('sample_ligera2_multi.txt')
# data3m <- read_tsv('sample_ligera2_bed_multi.txt')
# # plot( data1$pval, data2$pval )
# # cor( data1$pval, data2$pval )
# # # [1] 0.9999781
# # cor( data2$pval, data3$pval )
# # # [1] 1
# mean( abs( data1$pval - data2$pval ) )
# # [1] 0.001345751
# mean( abs( data2$pval - data3$pval ) )
# # [1] 2.952728e-10

# # in this case all multi versions agreed with the non-multi versions, though with better sample sizes we would expect --multi to be more powerful (and these p-values would be different)
# mean( abs( data1$pval - data1m$pval ) )
# # [1] 0
# mean( abs( data2$pval - data2m$pval ) )
# # [1] 0
# mean( abs( data3$pval - data3m$pval ) )
# # [1] 0
