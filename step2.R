#!/usr/bin/Rscript --slave
argv <- commandArgs(TRUE)
data <- read.csv(argv[1], header = T, sep = ";")

print(data$Nb_single)

data$Nb_total_gene = data$Nb_single + data$Nb_dupli
data$Pct_single = data$Nb_single / data$Nb_total_gene
data$Pct_dupli = data$Nb_dupli / data$Nb_total_gene
data$Nb_total_gene_in_list = data$Nb_single_in_list + data$Nb_dupli_in_list
data$Pct_single_in_list = data$Nb_single_in_list / data$Nb_total_gene_in_list
data$Pct_dupli_in_list = data$Nb_dupli_in_list / data$Nb_total_gene_in_list

fct_chi2 <- function(species_list){
  for (i in species_list){
    
    ## test du chi2
    pct = c(data$Pct_single[data$Species == i ], data$Pct_dupli[data$Species == i])
    list = c(data$Nb_single_in_list[data$Species == i], data$Nb_dupli_in_list[data$Species == i])
    chi2_list = chisq.test(list, p = pct)
    
    data$chi2[data$Species == i] = chi2_list$statistic
    data$chi2_pv[data$Species == i] = chi2_list$p.value
    
    ## Ajustement de BH pour les valeurs de pvalue du test de chi2
    p = data$chi2_pv
    data$BH_chi2_pv = p.adjust(p, method = "BH", n = length(p))
    
    ## test hyperg??om??trique
    data$hypergeom_pv = phyper((data$Nb_dupli_in_list - 1), data$Nb_dupli,
                              data$Nb_single, data$Nb_total_gene_in_list, lower.tail = FALSE )
    
    ## Ajustement de BH pour les valeurs de pvalue du test hyperg??om??trique
    q = data$hypergeom_pv
    data$BH_hgeom_pv = p.adjust(q, method = "BH", n = length(q))
    
  }
  return (data)
}

species_list = data$Species
data = fct_chi2(species_list)


sortie <- write.table(data, argv[1], sep = ";", dec = ",", row.names = FALSE, col.names = TRUE)




