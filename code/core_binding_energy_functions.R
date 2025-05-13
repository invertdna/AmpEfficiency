
#Core and helper functions for estimating DNA proportions from amplicon sequencing data
#implementing binding-energy model to estimate amplification efficiency biases and correct for them
#RPK March/April 2025

#Trinucleotide Binding Energy Matrix (ΔΔG37 in kcal/mol), inferred from Hooyberghs et al 2009, SantaLucia and Hicks 2004
create_binding_energy_matrix <- function() {
  
  labels <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT")
  
  values <- matrix(c(0.0,2.4,2.0,2.2,2.0,5.3,4.9,5.1,2.4,4.9,4.5,4.7,2.2,5.1,4.7,4.9,2.4,5.3,4.9,5.1,5.3,8.2,7.8,8.0,4.9,7.8,7.4,7.6,5.1,8.0,7.6,7.8,2.0,4.9,4.5,4.7,4.9,7.8,7.4,7.6,4.5,7.4,7.0,7.2,4.7,7.6,7.2,7.4,2.2,5.1,4.7,4.9,5.1,8.0,7.6,7.8,4.7,7.6,7.2,7.4,4.9,7.8,7.4,7.6,
                     2.4,0.0,2.2,2.6,5.3,2.4,5.1,5.5,4.9,2.0,4.7,5.1,5.1,2.2,4.9,5.3,5.3,2.4,5.1,5.5,8.2,5.3,8.0,8.4,7.8,4.9,7.6,8.0,8.0,5.1,7.8,8.2,4.9,2.0,4.7,5.1,7.8,4.9,7.6,8.0,7.4,4.5,7.2,7.6,7.6,4.7,7.4,7.8,5.1,2.2,4.9,5.3,8.0,5.1,7.8,8.2,7.6,4.7,7.4,7.8,7.8,4.9,7.6,8.0,
                     2.2,2.2,0.0,2.2,5.1,5.1,2.4,5.1,4.7,4.7,2.0,4.7,4.9,4.9,2.2,4.9,5.1,5.1,2.4,5.1,8.0,8.0,5.3,8.0,7.6,7.6,4.9,7.6,7.8,7.8,5.1,7.8,4.7,4.7,2.0,4.7,7.6,7.6,4.9,7.6,7.2,7.2,4.5,7.2,7.4,7.4,4.7,7.4,4.9,4.9,2.2,4.9,7.8,7.8,5.1,7.8,7.4,7.4,4.7,7.4,7.6,7.6,4.9,7.6,
                     2.2,2.6,1.7,0.0,5.1,5.5,4.6,2.4,4.7,5.1,4.2,2.0,4.9,5.3,4.4,2.2,5.1,5.5,4.6,2.4,8.0,8.4,7.5,5.3,7.6,8.0,7.1,4.9,7.8,8.2,7.3,5.1,4.7,5.1,4.2,2.0,7.6,8.0,7.1,4.9,7.2,7.6,6.7,4.5,7.4,7.8,6.9,4.7,4.9,5.3,4.4,2.2,7.8,8.2,7.3,5.1,7.4,7.8,6.9,4.7,7.6,8.0,7.1,4.9,
                     2.0,5.3,4.9,5.1,0.0,2.4,2.0,2.2,2.2,5.1,4.7,4.9,2.6,5.5,5.1,5.3,5.3,8.2,7.8,8.0,2.4,5.3,4.9,5.1,5.1,8.0,7.6,7.8,5.5,8.4,8.0,8.2,4.9,7.8,7.4,7.6,2.0,4.9,4.5,4.7,4.7,7.6,7.2,7.4,5.1,8.0,7.6,7.8,5.1,8.0,7.6,7.8,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,
                     5.3,2.4,5.1,5.5,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,5.5,2.6,5.3,5.7,8.2,5.3,8.0,8.4,5.3,2.4,5.1,5.5,8.0,5.1,7.8,8.2,8.4,5.5,8.2,8.6,7.8,4.9,7.6,8.0,4.9,2.0,4.7,5.1,7.6,4.7,7.4,7.8,8.0,5.1,7.8,8.2,8.0,5.1,7.8,8.2,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,
                     5.1,5.1,2.4,5.1,2.2,2.2,0.0,2.2,4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,8.0,8.0,5.3,8.0,5.1,5.1,2.4,5.1,7.8,7.8,5.1,7.8,8.2,8.2,5.5,8.2,7.6,7.6,4.9,7.6,4.7,4.7,2.0,4.7,7.4,7.4,4.7,7.4,7.8,7.8,5.1,7.8,7.8,7.8,5.1,7.8,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,
                     5.1,5.5,4.6,2.4,2.2,2.6,1.7,0.0,4.9,5.3,4.4,2.2,5.3,5.7,4.8,2.6,8.0,8.4,7.5,5.3,5.1,5.5,4.6,2.4,7.8,8.2,7.3,5.1,8.2,8.6,7.7,5.5,7.6,8.0,7.1,4.9,4.7,5.1,4.2,2.0,7.4,7.8,6.9,4.7,7.8,8.2,7.3,5.1,7.8,8.2,7.3,5.1,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,
                     2.4,5.1,4.7,4.9,2.2,5.1,4.7,4.9,0.0,2.4,2.0,2.2,2.2,5.1,4.7,4.9,5.1,8.0,7.6,7.8,5.1,8.0,7.6,7.8,2.4,5.3,4.9,5.1,5.1,8.0,7.6,7.8,4.7,7.6,7.2,7.4,4.7,7.6,7.2,7.4,2.0,4.9,4.5,4.7,4.7,7.6,7.2,7.4,4.9,7.8,7.4,7.6,4.9,7.8,7.4,7.6,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,
                     5.1,2.2,4.9,5.3,5.1,2.2,4.9,5.3,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,8.0,5.1,7.8,8.2,8.0,5.1,7.8,8.2,5.3,2.4,5.1,5.5,8.0,5.1,7.8,8.2,7.6,4.7,7.4,7.8,7.6,4.7,7.4,7.8,4.9,2.0,4.7,5.1,7.6,4.7,7.4,7.8,7.8,4.9,7.6,8.0,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,
                     4.9,4.9,2.2,4.9,4.9,4.9,2.2,4.9,2.2,2.2,0.0,2.2,4.9,4.9,2.2,4.9,7.8,7.8,5.1,7.8,7.8,7.8,5.1,7.8,5.1,5.1,2.4,5.1,7.8,7.8,5.1,7.8,7.4,7.4,4.7,7.4,7.4,7.4,4.7,7.4,4.7,4.7,2.0,4.7,7.4,7.4,4.7,7.4,7.6,7.6,4.9,7.6,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,
                     4.9,5.3,4.4,2.2,4.9,5.3,4.4,2.2,2.2,2.6,1.7,0.0,4.9,5.3,4.4,2.2,7.8,8.2,7.3,5.1,7.8,8.2,7.3,5.1,5.1,5.5,4.6,2.4,7.8,8.2,7.3,5.1,7.4,7.8,6.9,4.7,7.4,7.8,6.9,4.7,4.7,5.1,4.2,2.0,7.4,7.8,6.9,4.7,7.6,8.0,7.1,4.9,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,
                     2.2,5.1,4.7,4.9,2.6,5.5,5.1,5.3,1.7,4.6,4.2,4.4,0.0,2.4,2.0,2.2,5.1,8.0,7.6,7.8,5.5,8.4,8.0,8.2,4.6,7.5,7.1,7.3,2.4,5.3,4.9,5.1,4.7,7.6,7.2,7.4,5.1,8.0,7.6,7.8,4.2,7.1,6.7,6.9,2.0,4.9,4.5,4.7,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,4.4,7.3,6.9,7.1,2.2,5.1,4.7,4.9,
                     5.1,2.2,4.9,5.3,5.5,2.6,5.3,5.7,4.6,1.7,4.4,4.8,2.4,0.0,2.2,2.6,8.0,5.1,7.8,8.2,8.4,5.5,8.2,8.6,7.5,4.6,7.3,7.7,5.3,2.4,5.1,5.5,7.6,4.7,7.4,7.8,8.0,5.1,7.8,8.2,7.1,4.2,6.9,7.3,4.9,2.0,4.7,5.1,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,7.3,4.4,7.1,7.5,5.1,2.2,4.9,5.3,
                     4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,4.4,4.4,1.7,4.4,2.2,2.2,0.0,2.2,7.8,7.8,5.1,7.8,8.2,8.2,5.5,8.2,7.3,7.3,4.6,7.3,5.1,5.1,2.4,5.1,7.4,7.4,4.7,7.4,7.8,7.8,5.1,7.8,6.9,6.9,4.2,6.9,4.7,4.7,2.0,4.7,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,7.1,7.1,4.4,7.1,4.9,4.9,2.2,4.9,
                     4.9,5.3,4.4,2.2,5.3,5.7,4.8,2.6,4.4,4.8,3.9,1.7,2.2,2.6,1.7,0.0,7.8,8.2,7.3,5.1,8.2,8.6,7.7,5.5,7.3,7.7,6.8,4.6,5.1,5.5,4.6,2.4,7.4,7.8,6.9,4.7,7.8,8.2,7.3,5.1,6.9,7.3,6.4,4.2,4.7,5.1,4.2,2.0,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,7.1,7.5,6.6,4.4,4.9,5.3,4.4,2.2,
                     2.4,5.3,4.9,5.1,5.3,8.2,7.8,8.0,4.9,7.8,7.4,7.6,5.1,8.0,7.6,7.8,0.0,2.4,2.0,2.2,2.4,5.3,4.9,5.1,2.0,4.9,4.5,4.7,2.2,5.1,4.7,4.9,2.2,5.1,4.7,4.9,5.1,8.0,7.6,7.8,4.7,7.6,7.2,7.4,4.9,7.8,7.4,7.6,2.6,5.5,5.1,5.3,5.5,8.4,8.0,8.2,5.1,8.0,7.6,7.8,5.3,8.2,7.8,8.0,
                     5.3,2.4,5.1,5.5,8.2,5.3,8.0,8.4,7.8,4.9,7.6,8.0,8.0,5.1,7.8,8.2,2.4,0.0,2.2,2.6,2.8,2.4,5.1,5.5,4.9,2.0,4.7,5.1,5.1,2.2,4.9,5.3,5.1,2.2,4.9,5.3,8.0,5.1,7.8,8.2,7.6,4.7,7.4,7.8,7.8,4.9,7.6,8.0,5.5,2.6,5.3,5.7,8.4,5.5,8.2,8.6,8.0,5.1,7.8,8.2,8.2,5.3,8.0,8.4,
                     5.1,5.1,2.4,5.1,8.0,8.0,5.3,8.0,7.6,7.6,4.9,7.6,7.8,7.8,5.1,7.8,2.2,2.2,0.0,2.2,5.1,5.1,2.4,5.1,4.7,4.7,2.0,4.7,4.9,4.9,2.2,4.9,4.9,4.9,2.2,4.9,7.8,7.8,5.1,7.8,7.4,7.4,4.7,7.4,7.6,7.6,4.9,7.6,5.3,5.3,2.6,5.3,8.2,8.2,5.5,8.2,7.8,7.8,5.1,7.8,8.0,8.0,5.3,8.0,
                     5.1,5.5,4.6,2.4,8.0,8.4,7.5,5.3,7.6,8.0,7.1,4.9,7.8,8.2,7.3,5.1,2.2,2.6,1.7,0.0,5.1,5.5,4.6,2.4,4.7,5.1,4.2,2.0,4.9,5.3,4.4,2.2,4.9,5.3,4.4,2.2,7.8,8.2,7.3,5.1,7.4,7.8,6.9,4.7,7.6,8.0,7.1,4.9,5.3,5.7,4.8,2.6,8.2,8.6,7.7,5.5,7.8,8.2,7.3,5.1,8.0,8.4,7.5,5.3,
                     5.3,8.2,7.8,8.0,2.4,5.3,4.9,5.1,5.1,8.0,7.6,7.8,5.5,8.4,8.0,8.2,2.4,2.8,4.9,5.1,0.0,2.4,2.0,2.2,2.2,3.0,4.7,4.9,2.6,3.0,5.1,5.3,5.1,8.0,7.6,7.8,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,5.5,8.4,8.0,8.2,2.6,5.5,5.1,5.3,5.3,8.2,7.8,8.0,5.7,8.6,8.2,8.4,
                     8.2,5.3,8.0,8.4,5.3,2.4,5.1,5.5,8.0,5.1,7.8,8.2,8.4,5.5,8.2,8.6,5.3,2.4,5.1,5.5,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,5.5,2.6,5.3,5.7,8.0,5.1,7.8,8.2,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,8.4,5.5,8.2,8.6,5.5,2.6,5.3,5.7,8.2,5.3,8.0,8.4,8.6,5.7,8.4,8.8,
                     8.0,8.0,5.3,8.0,5.1,5.1,2.4,5.1,7.8,7.8,5.1,7.8,8.2,8.2,5.5,8.2,5.1,5.1,2.4,5.1,2.2,2.2,0.0,2.2,4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,7.8,7.8,5.1,7.8,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,8.2,8.2,5.5,8.2,5.3,5.3,2.6,5.3,8.0,8.0,5.3,8.0,8.4,8.4,5.7,8.4,
                     8.0,8.4,7.5,5.3,5.1,5.5,4.6,2.4,7.8,8.2,7.3,5.1,8.2,8.6,7.7,5.5,5.1,5.5,4.6,2.4,2.2,2.6,1.7,0.0,4.9,5.3,4.4,2.2,5.3,5.7,4.8,2.6,7.8,8.2,7.3,5.1,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,8.2,8.6,7.7,5.5,5.3,5.7,4.8,2.6,8.0,8.4,7.5,5.3,8.4,8.8,7.9,5.7,
                     5.1,8.0,7.6,7.8,5.1,8.0,7.6,7.8,2.4,5.3,4.9,5.1,5.1,8.0,7.6,7.8,2.2,5.1,4.7,4.9,2.2,5.1,4.7,4.9,0.0,2.4,2.0,2.2,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,4.9,7.8,7.4,7.6,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,5.3,8.2,7.8,8.0,2.6,5.5,5.1,5.3,5.3,8.2,7.8,8.0,
                     8.0,5.1,7.8,8.2,8.0,5.1,7.8,8.2,5.3,2.4,5.1,5.5,8.0,5.1,7.8,8.2,5.1,2.2,4.9,5.3,3.0,2.2,4.9,5.3,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,8.2,5.3,8.0,8.4,5.5,2.6,5.3,5.7,8.2,5.3,8.0,8.4,
                     7.8,7.8,5.1,7.8,7.8,7.8,5.1,7.8,5.1,5.1,2.4,5.1,7.8,7.8,5.1,7.8,4.9,4.9,2.2,4.9,4.9,4.9,2.2,4.9,2.2,2.2,0.0,2.2,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,8.0,8.0,5.3,8.0,5.3,5.3,2.6,5.3,8.0,8.0,5.3,8.0,
                     7.8,8.2,7.3,5.1,7.8,8.2,7.3,5.1,5.1,5.5,4.6,2.4,7.8,8.2,7.3,5.1,4.9,5.3,4.4,2.2,4.9,5.3,4.4,2.2,2.2,2.6,1.7,0.0,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,8.0,8.4,7.5,5.3,5.3,5.7,4.8,2.6,8.0,8.4,7.5,5.3,
                     5.1,8.0,7.6,7.8,5.5,8.4,8.0,8.2,4.6,7.5,7.1,7.3,2.4,5.3,4.9,5.1,2.2,5.1,4.7,4.9,2.6,5.5,5.1,5.3,1.7,4.6,4.2,4.4,0.0,2.4,2.0,2.2,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,4.4,7.3,6.9,7.1,2.2,5.1,4.7,4.9,5.3,8.2,7.8,8.0,5.7,8.6,8.2,8.4,4.8,7.7,7.3,7.5,2.6,5.5,5.1,5.3,
                     8.0,5.1,7.8,8.2,8.4,5.5,8.2,8.6,7.5,4.6,7.3,7.7,5.3,2.4,5.1,5.5,5.1,2.2,4.9,5.3,3.0,2.6,5.3,5.7,4.6,1.7,4.4,4.8,2.4,0.0,2.2,2.6,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,7.3,4.4,7.1,7.5,5.1,2.2,4.9,5.3,8.2,5.3,8.0,8.4,8.6,5.7,8.4,8.8,7.7,4.8,7.5,7.9,5.5,2.6,5.3,5.7,
                     7.8,7.8,5.1,7.8,8.2,8.2,5.5,8.2,7.3,7.3,4.6,7.3,5.1,5.1,2.4,5.1,4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,4.4,4.4,1.7,4.4,2.2,2.2,0.0,2.2,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,7.1,7.1,4.4,7.1,4.9,4.9,2.2,4.9,8.0,8.0,5.3,8.0,8.4,8.4,5.7,8.4,7.5,7.5,4.8,7.5,5.3,5.3,2.6,5.3,
                     7.8,8.2,7.3,5.1,8.2,8.6,7.7,5.5,7.3,7.7,6.8,4.6,5.1,5.5,4.6,2.4,4.9,5.3,4.4,2.2,5.3,5.7,4.8,2.6,4.4,4.8,3.9,1.7,2.2,2.6,1.7,0.0,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,7.1,7.5,6.6,4.4,4.9,5.3,4.4,2.2,8.0,8.4,7.5,5.3,8.4,8.8,7.9,5.7,7.5,7.9,7.0,4.8,5.3,5.7,4.8,2.6,
                     2.2,5.1,4.7,4.9,5.1,8.0,7.6,7.8,4.7,7.6,7.2,7.4,4.9,7.8,7.4,7.6,2.2,5.1,4.7,4.9,5.1,8.0,7.6,7.8,4.7,7.6,7.2,7.4,4.9,7.8,7.4,7.6,0.0,2.4,2.0,2.2,2.4,5.3,4.9,5.1,2.0,4.9,4.5,4.7,2.2,5.1,4.7,4.9,2.2,5.1,4.7,4.9,5.1,8.0,7.6,7.8,4.7,7.6,7.2,7.4,4.9,7.8,7.4,7.6,
                     5.1,2.2,4.9,5.3,8.0,5.1,7.8,8.2,7.6,4.7,7.4,7.8,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,8.0,5.1,7.8,8.2,7.6,4.7,7.4,7.8,7.8,4.9,7.6,8.0,2.4,0.0,2.2,2.6,5.3,2.4,5.1,5.5,4.9,2.0,4.7,5.1,5.1,2.2,4.9,5.3,5.1,2.2,4.9,5.3,8.0,5.1,7.8,8.2,7.6,4.7,7.4,7.8,7.8,4.9,7.6,8.0,
                     4.9,4.9,2.2,4.9,7.8,7.8,5.1,7.8,7.4,7.4,4.7,7.4,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,7.8,7.8,5.1,7.8,7.4,7.4,4.7,7.4,7.6,7.6,4.9,7.6,2.2,2.2,0.0,2.2,5.1,5.1,2.4,5.1,2.5,4.7,2.0,4.7,4.9,4.9,2.2,4.9,4.9,4.9,2.2,4.9,7.8,7.8,5.1,7.8,7.4,7.4,4.7,7.4,7.6,7.6,4.9,7.6,
                     4.9,5.3,4.4,2.2,7.8,8.2,7.3,5.1,7.4,7.8,6.9,4.7,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,7.8,8.2,7.3,5.1,7.4,7.8,6.9,4.7,7.6,8.0,7.1,4.9,2.2,2.6,1.7,0.0,5.1,5.5,4.6,2.4,4.7,5.1,4.2,2.0,4.9,5.3,4.4,2.2,4.9,5.3,4.4,2.2,7.8,8.2,7.3,5.1,7.4,7.8,6.9,4.7,7.6,8.0,7.1,4.9,
                     5.1,8.0,7.6,7.8,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,5.1,8.0,7.6,7.8,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,2.4,5.3,4.9,5.1,0.0,2.4,2.0,2.2,2.2,5.1,4.7,4.9,2.6,5.5,5.1,5.3,5.1,8.0,7.6,7.8,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,
                     8.0,5.1,7.8,8.2,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,8.0,5.1,7.8,8.2,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,5.3,2.4,5.1,5.5,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,5.5,2.6,5.3,5.7,8.0,5.1,7.8,8.2,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,
                     7.8,7.8,5.1,7.8,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,7.8,7.8,5.1,7.8,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,5.1,5.1,2.4,5.1,2.2,2.2,0.0,2.2,4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,7.8,7.8,5.1,7.8,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,
                     7.8,8.2,7.3,5.1,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,7.8,8.2,7.3,5.1,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,5.1,5.5,4.6,2.4,2.2,2.6,1.7,0.0,4.9,5.3,4.4,2.2,5.3,5.7,4.8,2.6,7.8,8.2,7.3,5.1,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,
                     4.9,7.8,7.4,7.6,4.9,7.8,7.4,7.6,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,4.9,7.8,7.4,7.6,4.9,7.8,7.4,7.6,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,2.2,5.1,2.5,4.9,2.2,5.1,4.7,4.9,0.0,2.4,1.8,2.2,2.2,5.1,2.2,4.9,4.9,7.8,7.4,7.6,4.9,7.8,7.4,7.6,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,
                     7.8,4.9,7.6,8.0,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,7.8,4.9,7.6,8.0,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,5.1,2.2,4.9,5.3,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,
                     7.6,7.6,4.9,7.6,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,7.6,7.6,4.9,7.6,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,4.9,4.9,2.2,4.9,1.8,2.2,0.0,2.2,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,
                     7.6,8.0,7.1,4.9,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,7.6,8.0,7.1,4.9,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,4.9,5.3,4.4,2.2,2.2,2.6,1.7,0.0,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,
                     4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,4.4,7.3,6.9,7.1,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,4.4,7.3,6.9,7.1,2.2,5.1,4.7,4.9,2.2,5.1,4.7,4.9,2.6,5.5,5.1,5.3,1.7,4.6,4.2,4.4,0.0,2.4,2.0,2.2,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,4.4,7.3,6.9,7.1,2.2,5.1,4.7,4.9,
                     7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,7.3,4.4,7.1,7.5,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,7.3,4.4,7.1,7.5,5.1,2.2,4.9,5.3,5.1,2.2,4.9,5.3,5.5,2.6,5.3,5.7,4.6,1.7,4.4,4.8,2.4,0.0,2.2,2.6,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,7.3,4.4,7.1,7.5,5.1,2.2,4.9,5.3,
                     7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,7.1,7.1,4.4,7.1,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,7.1,7.1,4.4,7.1,4.9,4.9,2.2,4.9,4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,2.2,4.4,1.7,4.4,2.2,2.2,0.0,2.2,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,7.1,7.1,4.4,7.1,4.9,4.9,2.2,4.9,
                     7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,7.1,7.5,6.6,4.4,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,7.1,7.5,6.6,4.4,4.9,5.3,4.4,2.2,4.9,5.3,4.4,2.2,5.3,5.7,4.8,2.6,4.4,4.8,3.9,1.7,2.2,2.6,1.7,0.0,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,7.1,7.5,6.6,4.4,4.9,5.3,4.4,2.2,
                     2.2,5.1,4.7,4.9,5.1,8.0,7.6,7.8,4.7,7.6,7.2,7.4,4.9,7.8,7.4,7.6,2.6,5.5,5.1,5.3,5.5,8.4,8.0,8.2,5.1,8.0,7.6,7.8,5.3,8.2,7.8,8.0,1.7,4.6,4.2,4.4,4.6,7.5,7.1,7.3,4.2,7.1,6.7,6.9,4.4,7.3,6.9,7.1,0.0,2.4,2.0,2.2,2.4,5.3,4.9,5.1,2.0,4.9,4.5,4.7,2.2,5.1,4.7,4.9,
                     5.1,2.2,4.9,5.3,8.0,5.1,7.8,8.2,7.6,4.7,7.4,7.8,7.8,4.9,7.6,8.0,5.5,2.6,5.3,5.7,8.4,5.5,8.2,8.6,8.0,5.1,7.8,8.2,8.2,5.3,8.0,8.4,4.6,1.7,4.4,4.8,7.5,4.6,7.3,7.7,7.1,4.2,6.9,7.3,7.3,4.4,7.1,7.5,2.4,0.0,2.2,2.6,5.3,2.4,5.1,5.5,4.9,2.0,4.7,5.1,5.1,2.2,4.9,5.3,
                     4.9,4.9,2.2,4.9,7.8,7.8,5.1,7.8,7.4,7.4,4.7,7.4,7.6,7.6,4.9,7.6,5.3,5.3,2.6,5.3,8.2,8.2,5.5,8.2,7.8,7.8,5.1,7.8,8.0,8.0,5.3,8.0,4.4,4.4,1.7,4.4,7.3,7.3,4.6,7.3,6.9,6.9,4.2,6.9,7.1,7.1,4.4,7.1,2.2,2.2,0.0,2.2,5.1,5.1,2.4,5.1,4.7,4.7,2.0,4.7,4.9,4.9,2.2,4.9,
                     4.9,5.3,4.4,2.2,7.8,8.2,7.3,5.1,7.4,7.8,6.9,4.7,7.6,8.0,7.1,4.9,5.3,5.7,4.8,2.6,8.2,8.6,7.7,5.5,7.8,8.2,7.3,5.1,8.0,8.4,7.5,5.3,4.4,4.8,3.9,1.7,7.3,7.7,6.8,4.6,6.9,7.3,6.4,4.2,7.1,7.5,6.6,4.4,2.2,2.6,1.7,0.0,5.1,5.5,4.6,2.4,4.7,5.1,4.2,2.0,2.4,5.3,4.4,2.2,
                     5.1,8.0,7.6,7.8,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,5.5,8.4,8.0,8.2,2.6,5.5,5.1,5.3,5.3,8.2,7.8,8.0,5.7,8.6,8.2,8.4,4.6,7.5,7.1,7.3,1.7,4.6,4.2,4.4,4.4,7.3,6.9,7.1,4.8,7.7,7.3,7.5,2.4,5.3,4.9,5.1,0.0,2.4,2.0,2.2,2.2,5.1,4.7,4.9,2.6,5.5,5.1,5.3,
                     8.0,5.1,7.8,8.2,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,8.4,5.5,8.2,8.6,5.5,2.6,5.3,5.7,8.2,5.3,8.0,8.4,8.6,5.7,8.4,8.8,7.5,4.6,7.3,7.7,4.6,1.7,4.4,4.8,7.3,4.4,7.1,7.5,7.7,4.8,7.5,7.9,5.3,2.4,5.1,5.5,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,5.5,2.6,5.3,5.7,
                     7.8,7.8,5.1,7.8,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,8.2,8.2,5.5,8.2,5.3,5.3,2.6,5.3,8.0,8.0,5.3,8.0,8.4,8.4,5.7,8.4,7.3,7.3,4.6,7.3,4.4,4.4,1.7,4.4,7.1,7.1,4.4,7.1,7.5,7.5,4.8,7.5,5.1,5.1,2.4,5.1,2.2,2.2,0.0,2.2,4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,
                     7.8,8.2,7.3,5.1,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,8.2,8.6,7.7,5.5,5.3,5.7,4.8,2.6,8.0,8.4,7.5,5.3,8.4,8.8,7.9,5.7,7.3,7.7,6.8,4.6,4.4,4.8,3.9,1.7,7.1,7.5,6.6,4.4,7.5,7.9,7.0,4.8,5.1,5.5,4.6,2.4,2.2,2.6,1.7,0.0,4.9,5.3,4.4,2.2,2.5,5.7,4.8,2.6,
                     4.9,7.8,7.4,7.6,4.9,7.8,7.4,7.6,2.2,5.1,4.7,4.9,4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,5.3,8.2,7.8,8.0,2.6,5.5,5.1,5.3,5.3,8.2,7.8,8.0,4.4,7.3,6.9,7.1,4.4,7.3,6.9,7.1,1.7,4.6,4.2,4.4,4.4,7.3,6.9,7.1,2.2,5.1,4.7,4.9,2.2,5.1,4.7,4.9,0.0,2.4,2.0,2.2,2.2,5.1,4.7,4.9,
                     7.8,4.9,7.6,8.0,7.8,4.9,7.6,8.0,5.1,2.2,4.9,5.3,7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,8.2,5.3,8.0,8.4,5.5,2.6,5.3,5.7,8.2,5.3,8.0,8.4,7.3,4.4,7.1,7.5,7.3,4.4,7.1,7.5,4.6,1.7,4.4,4.8,7.3,4.4,7.1,7.5,5.1,2.2,4.9,5.3,5.1,2.2,4.9,5.3,2.4,0.0,2.2,2.6,5.1,2.2,4.9,5.3,
                     7.6,7.6,4.9,7.6,7.6,7.6,4.9,7.6,4.9,4.9,2.2,4.9,7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,8.0,8.0,5.3,8.0,5.3,5.3,2.6,5.3,8.0,8.0,5.3,8.0,7.1,7.1,4.4,7.1,7.1,7.1,4.4,7.1,4.4,4.4,1.7,4.4,7.1,7.1,4.4,7.1,4.9,4.9,2.2,4.9,4.9,4.9,2.2,4.9,2.2,2.2,0.0,2.2,4.9,4.9,2.2,4.9,
                     7.6,8.0,7.1,4.9,7.6,8.0,7.1,4.9,4.9,5.3,4.4,2.2,7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,8.0,8.4,7.5,5.3,5.3,5.7,4.8,2.6,8.0,8.4,7.5,5.3,7.1,7.5,6.6,4.4,7.1,7.5,6.6,4.4,4.4,4.8,3.9,1.7,7.1,7.5,6.6,4.4,4.9,5.3,4.4,2.2,4.9,5.3,4.4,2.2,2.2,2.6,1.7,0.0,2.2,5.3,4.4,2.2,
                     4.9,7.8,7.4,7.6,5.3,8.2,7.8,8.0,4.4,7.3,6.9,7.1,2.2,5.1,4.7,4.9,5.3,8.2,7.8,8.0,5.7,8.6,8.2,8.4,4.8,7.7,7.3,7.5,2.6,5.5,5.1,5.3,4.4,7.3,6.9,7.1,4.8,7.7,7.3,7.5,3.9,6.8,6.4,6.6,1.7,4.6,4.2,4.4,2.2,5.1,4.7,2.4,2.6,5.5,5.1,2.5,1.7,4.6,4.2,2.2,0.0,2.4,2.0,2.2,
                     7.8,4.9,7.6,8.0,8.2,5.3,8.0,8.4,7.3,4.4,7.1,7.5,5.1,2.2,4.9,5.3,8.2,5.3,8.0,8.4,8.6,5.7,8.4,8.8,7.7,4.8,7.5,7.9,5.5,2.6,5.3,5.7,7.3,4.4,7.1,7.5,7.7,4.8,7.5,7.9,6.8,3.9,6.6,7.0,4.6,1.7,4.4,4.8,5.1,2.2,4.9,5.3,5.5,2.6,5.3,5.7,4.6,1.7,4.4,4.8,2.4,0.0,2.2,2.6,
                     7.6,7.6,4.9,7.6,8.0,8.0,5.3,8.0,7.1,7.1,4.4,7.1,4.9,4.9,2.2,4.9,8.0,8.0,5.3,8.0,8.4,8.4,5.7,8.4,7.5,7.5,4.8,7.5,5.3,5.3,2.6,5.3,7.1,7.1,4.4,7.1,7.5,7.5,4.8,7.5,6.6,6.6,3.9,6.6,4.4,4.4,1.7,4.4,4.9,4.9,2.2,4.9,5.3,5.3,2.6,5.3,4.4,4.4,1.7,4.4,2.2,2.2,0.0,2.2,
                     7.6,8.0,7.1,4.9,8.0,8.4,7.5,5.3,7.1,7.5,6.6,4.4,4.9,5.3,4.4,2.2,8.0,8.4,7.5,5.3,8.4,8.8,7.9,5.7,7.5,7.9,7.0,4.8,5.3,5.7,4.8,2.6,7.1,7.5,6.6,4.4,7.5,7.9,7.0,4.8,6.6,7.0,6.1,3.9,4.4,4.8,3.9,1.7,4.9,5.3,4.4,2.2,5.3,5.7,4.8,2.6,4.4,4.8,3.9,1.7,2.2,2.6,1.7,0.0),
                   nrow = 64, ncol = 64, byrow = TRUE)
  
  # Assign row and column names
  rownames(values) <- labels
  colnames(values) <- labels
  
  # Return the matrix
  return(values)
}

# Reverse complement function using vectorized lookup
reverse_complement <- function(seq) {
  complements <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G",
                   "R" = "Y", "Y" = "R", "S" = "S", "W" = "W",
                   "K" = "M", "M" = "K", "B" = "V", "D" = "H",
                   "H" = "D", "V" = "B", "N" = "N")
  chars <- unlist(strsplit(seq, ""))
  rev_comp <- rev(sapply(chars, function(x) complements[x]))
  paste(rev_comp, collapse = "")
}

complement <- function(seq) {
  complements <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G",
                   "R" = "Y", "Y" = "R", "S" = "S", "W" = "W",
                   "K" = "M", "M" = "K", "B" = "V", "D" = "H",
                   "H" = "D", "V" = "B", "N" = "N")
  chars <- unlist(strsplit(seq, ""))
  comp <- sapply(chars, function(x) complements[x])
  paste(comp, collapse = "")
}

# IUPAC codes
iupac_codes <- list(
  "A" = c("A"), "C" = c("C"), "G" = c("G"), "T" = c("T"),
  "R" = c("A", "G"), "Y" = c("C", "T"), "S" = c("G", "C"),
  "W" = c("A", "T"), "K" = c("G", "T"), "M" = c("A", "C"),
  "B" = c("C", "G", "T"), "D" = c("A", "G", "T"),
  "H" = c("A", "C", "T"), "V" = c("A", "C", "G"), "N" = c("A", "C", "G", "T")
)

reverse <- function(seq){
  f <- rev(unlist(strsplit(seq, "")))
  paste(f, collapse = "")
}

expand_primer <- function(p1){
  
  iupac_ambig <- list(
    #"A" = c("A"), "C" = c("C"), "G" = c("G"), "T" = c("T"),
    "R" = c("A", "G"), "Y" = c("C", "T"), "S" = c("G", "C"),
    "W" = c("A", "T"), "K" = c("G", "T"), "M" = c("A", "C"),
    "B" = c("C", "G", "T"), "D" = c("A", "G", "T"),
    "H" = c("A", "C", "T"), "V" = c("A", "C", "G"), "N" = c("A", "C", "G", "T")
  )
  
  f <- strsplit(p1, split = "")[[1]]
  iupac_indices <- which(f %in% names(iupac_ambig))
  
  if(length(iupac_indices) == 0){
    return(p1)
  } else {
    iupac_combinations <- lapply(iupac_indices, function(i) {
      iupac_ambig[[f[i]]]
    })
    
    # Generate all combinations of the IUPAC codes
    all_combinations <- do.call(expand.grid, iupac_combinations)
    all_df <- t(replicate(nrow(all_combinations), f))
    
    for (i in 1:length(iupac_indices)) {
      all_df[, iupac_indices[i]] <- as.character(all_combinations[, i])
    }
    
    # Return the unique combinations
    return(apply(all_df, 1, function(x) paste(x, collapse = "")))
  }
}

substring_dataframe <- function(x, length_primer = 20, length_substring = 3) {
  # Input: vector x of character strings
  # Output: dataframe with each row corresponding to an element of x
  #         and columns containing substrings with varying start/end positions
  
  # Check if input is a character vector
  if (!is.character(x)) {
    stop("Input must be a character vector")
  }
  
  # Create start and end position vectors
  start_pos <- 1:(length_primer-length_substring+1)
  end_pos <- length_substring:length_primer
  
  # Create column names
  col_names <- paste0("substr_", start_pos, "_to_", end_pos)
  
  # Use mapply to efficiently generate all substrings
  result_list <- mapply(
    function(start, end) {
      substring(x, start, end)
    },
    start_pos, end_pos,
    SIMPLIFY = FALSE
  )
  
  # Convert list to data frame
  result_df <- do.call(cbind.data.frame, result_list)
  
  # Set column names
  names(result_df) <- col_names
  
  return(result_df)
}



weight_function <- function(x, primer_length, rev = FALSE, max_weight = 5, min_weight = 1){
  # x is the distance from the 3' end of primer
  # primer_length is the length of the primer
  # rev is a logical indicating if the primer being evaluated in reverse direction
  # returns a weight for the binding energy at position x
  a <- -(primer_length - 5)
  b <- 1
  
  if (rev){
    w <- primer_length - x
  } else {
    w <- x
  }
  # sigmoidal weighting function
  return(
    min_weight + ((max_weight - min_weight)*exp(a+b*w)/
           (1+exp(a+b*w)))
  ) 
}

get_triplet_binding_energy <- function(primer, template, binding_energy_mat = m, REV = FALSE){
  # NOTE: handles ambiguities in the primer, but not in the template
  # NOTE: primer and template both handled in the 5' - 3' orientation, such that the canonical match for primer AAA is template AAA in binding energy matrix
  # binding energy matrix `m` is a 64x64 matrix that has row names and column names corresponding to all possible triplets of nucleotides A,G,C,T
  
  p1 <- ifelse(REV == FALSE, primer, reverse_complement(primer))
  t1 <- template #ifelse(REV == FALSE, template, complement(template))
  
  p1 <- expand_primer(p1)
  pmat <- matrix(NA, nrow = length(p1), ncol = nchar(t1) - nchar(p1[1]) + 1) #matrix for temporarily storing values
  
  for (j in 1:length(p1)){ #to handle cases in which primer ambiguities mean a mixture of primer seqs
    
    length_p <- nchar(p1[j])
    length_t <- nchar(t1)
    p_triplets <- as.vector(mapply(substr, p1[j], 1:(length_p-2), 3:(length_p))) #get primer triplets
    t_potential_positions <- as.vector(mapply(substr, t1, 1:(length_t-length_p+1), length_p:length_t)) #index positions
    t_triplets <- substring_dataframe(t_potential_positions, length_primer = length_p, length_substring = 3) #get template triplets
    
    #get row/col indices of triplets, with which to efficiently find individual energies
    p_triplets_idx <- match(p_triplets,rownames(m))
    t_triplets_idx <- apply(t_triplets, MARGIN = 2, FUN = match, colnames(m))
    
    #matrix `t_triplets_idx` gives the indices of triplets in different potential primer-binding positions in rows, and the triplets implied by 
    #each starting position in columns. For a given starting position (row), the binding energy will be given by the rowSum
    #of binding energies considering the primer triplets and template triplets for that row
    
    #calculate binding energy, relative to perfect binding
    binding_raw_mat <- t(apply(t_triplets_idx, MARGIN = 1, FUN = function(x){diag(m[p_triplets_idx, x])}))
    
    #apply weighting function so as to penalize 3' mismatches to a greater extent
    binding_mat_weighted <- t(t(binding_raw_mat) * 
                                weight_function(1:(length_p-2), length_p-2, max_weight = 3, rev = REV)
    )
  
    pmat[j,] <- rowSums(binding_mat_weighted)
  }
  
  return(colMeans(pmat)) #return mean across primer variations
}
  
m <- create_binding_energy_matrix() #units kcal/mol at 37C

specificity_score <- function(g_vec){
  #for a vector of binding energies, sort and give (1 - ratio of probs of two lowest-energy scores), such that 
  #a specificity score of 0 means no difference between the two lowest scores, and 1 means unambiguous binding to lowest
  
  names(g_vec) <- 1:length(g_vec) #create index, storing info as names in vector
  g_vec <- sort(g_vec) #sort
  g_vec_mean <- mean(g_vec)
  g_vec_sd <- sd(g_vec)
  p1 <- pnorm(g_vec[1], g_vec_mean, g_vec_sd)
  p2 <- pnorm(g_vec[2], g_vec_mean, g_vec_sd)
  
  return(
    1 - ((p2 - p1)/p2)
  )
}

predict_binding <- function(forward_primer, reverse_primer, template){
  
  #check input
  Nidx <- which(!strsplit(template, split ="")[[1]] %in% c("A","G","C","T"))
  
  if (length(Nidx) > 0) {
    # cat("Non-standard base in template sequence (probably an N or similar); can't handle this right now")
    # return(NA)
    
    #alternatively, randomly insert nucleotides in those N positions
    template <- gsub("N", sample(c("A","G","C","T"), 1), template)
    foo <- strsplit(template, split = "")[[1]]
    foo[Nidx] <- sample(c("A","G","C","T"), length(Nidx), replace = TRUE)
    template <- paste(foo, collapse = "")
  } 
  
  
  
  #forward primer positions
  g_vec_fwd <- get_triplet_binding_energy(forward_primer, template, binding_energy_mat = m)
  fwd_specificity <- specificity_score(g_vec_fwd)
  
  #reverse primer positions
  g_vec_rev <- get_triplet_binding_energy(reverse_primer, template, binding_energy_mat = m, REV = TRUE)
  rev_specificity <- specificity_score(g_vec_rev)
  
  #output
  best_forward_binding_position <- which(g_vec_fwd == min(g_vec_fwd)) #keeps both in case of ties

  best_reverse_binding_position <- which(g_vec_rev == min(g_vec_rev))

  overall_binding_energy <- g_vec_fwd[best_forward_binding_position[1]] + g_vec_rev[best_reverse_binding_position[1]]
  amplicon_size_w_primers <- ((best_reverse_binding_position + nchar(reverse_primer)) - 
                                (best_forward_binding_position)) 
  amplicon_size_no_primers <- amplicon_size_w_primers - nchar(forward_primer) - nchar(reverse_primer)
  amplicon_w_primers <- substr(template, best_forward_binding_position[1], 
                               (best_reverse_binding_position[1] + nchar(reverse_primer) - 1)) #only keeps first, if there are multiple
  amplicon_no_primers <- substr(template, (best_forward_binding_position[1] + 1), 
                                (best_reverse_binding_position[1] - 1)) #only keeps first, if there are multiple
  
  return(list("best_forward_binding_position" = best_forward_binding_position,
              "forward_specificity" = unname(fwd_specificity),
              "reverse_specificity" = unname(rev_specificity),
              "warnings" = ifelse(unname(fwd_specificity) < 0.05 & unname(rev_specificity) < 0.05, #if second-best binding is more than 1/20th prob of best binding site, warn.
                                  "NONE", "Potential nonspecific binding"),
              "best_reverse_binding_position" = best_reverse_binding_position,
              "overall_binding_energy" = overall_binding_energy,
              "amplicon_size_w_primers" = ifelse(amplicon_size_w_primers > 0 & overall_binding_energy < 50, amplicon_size_w_primers, NA),
              "amplicon_size_no_primers" = ifelse(amplicon_size_no_primers > 0 & overall_binding_energy < 50, amplicon_size_no_primers, NA),
              "amplicon_w_primers" = ifelse(amplicon_size_w_primers > 0 & overall_binding_energy < 50, amplicon_w_primers, NA),
              "amplicon_no_primers" = ifelse(amplicon_size_no_primers > 0 & overall_binding_energy < 50, amplicon_no_primers, NA)))
}

#example
# predict_binding(
#   forward_primer = "GCCGGTAAAACTCGTGCCAGC",
#   reverse_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",
#   template = "GTTGGTAAATCTCGTGCCAGCCACCGCGGTTATACGAGAAACTCAAATTGAAAGTCTGCGGCGTAAAGGGTGATTATAGAAAATTTTACACATTAGAATTAAAAATTAATTAAGCTGTCGCACGCTTTTACTAATCAAAAACACAACAACGAAAGTGATTCTAAATAACTAATACACTTGAATTCACGACCCCTGAGACACAAACTGGGATTAGATACCCCACTATG"
# )

# predict_binding <- function(forward_primer, reverse_primer, template){
#   
#   #check input
#   Nidx <- which(!strsplit(template, split ="")[[1]] %in% c("A","G","C","T"))
#   
#   if (length(Nidx) > 0) {
#     # cat("Non-standard base in template sequence (probably an N or similar); can't handle this right now")
#     # return(NA)
# 
#     #alternatively, randomly insert nucleotides in those N positions
#      template <- gsub("N", sample(c("A","G","C","T"), 1), template)
#      foo <- strsplit(template, split = "")[[1]]
#      foo[Nidx] <- sample(c("A","G","C","T"), length(Nidx), replace = TRUE)
#      template <- paste(foo, collapse = "")
#       } 
#     
#   
#   
#   #iterate over forward primer positions
#   g_vec_fwd <- NA
#   for (pos in 1:(nchar(template)-nchar(forward_primer))) {
#     tmp <- substr(template, pos, (pos+nchar(forward_primer)-1))
#     g_vec_fwd[pos] <- get_binding_energy(forward_primer, tmp)
#   }
#   
#   #iterate over reverse primer positions, indexing from the leftward-most position of the binding position
#   #template in original 5' -> 3' orientation; rev primer in reverse complement
#   g_vec_rev <- NA
#   revcomp_reverse_primer <- reverse_complement(reverse_primer)
#   for (pos in 1:(nchar(template)-nchar(reverse_primer)+1)) {
#     tmp <- substr(template, pos, (pos+nchar(reverse_primer)-1))
#     g_vec_rev[pos] <- get_binding_energy(revcomp_reverse_primer, tmp, rev = TRUE)
#   }
#   
#   #initiation energy
#   deltaG_init = 3.4 #kcal/mol (or approximately +14.2 kJ/mol)
#   #SantaLucia gives initiation with G-C as 2.6, A-T as 1.7, and presumably these are per-primer. But then his fig 1 gives 1.03 as A-T initiation in example.
#   
#   
#   #output
#   best_forward_binding_position <- which(g_vec_fwd == min(g_vec_fwd)) #keeps both in case of ties
#   best_reverse_binding_position <- which(g_vec_rev == min(g_vec_rev))
#   overall_binding_energy <- g_vec_fwd[best_forward_binding_position[1]] + g_vec_rev[best_reverse_binding_position[1]] + deltaG_init
#   amplicon_size_w_primers <- ((best_reverse_binding_position + nchar(reverse_primer)) - 
#                                 (best_forward_binding_position)) 
#   amplicon_size_no_primers <- amplicon_size_w_primers - nchar(forward_primer) - nchar(reverse_primer)
#   amplicon_w_primers <- substr(template, best_forward_binding_position[1], 
#                                (best_reverse_binding_position[1] + nchar(reverse_primer) - 1)) #only keeps first, if there are multiple
#   amplicon_no_primers <- substr(template, (best_forward_binding_position[1] + 1), 
#                                 (best_reverse_binding_position[1] - 1)) #only keeps first, if there are multiple
#   
#   return(list("best_forward_binding_position" = best_forward_binding_position,
#               "best_reverse_binding_position" = best_reverse_binding_position,
#               "overall_binding_energy" = overall_binding_energy,
#               "amplicon_size_w_primers" = amplicon_size_w_primers,
#               "amplicon_size_no_primers" = amplicon_size_no_primers,
#               "amplicon_w_primers" = amplicon_w_primers,
#               "amplicon_no_primers" = amplicon_no_primers))
#   
#   }

organize_amplicons <- function(results_list){
  tmp <- lapply(X = 1:length(results_list),
         FUN = function(X){results_list[[X]]$amplicon_no_primers})
  return(unlist(tmp))
  }

find_identical <- function(results_list){
  #get amplicons from primer-binding result
  a <- organize_amplicons(results_list) 
  
  #is distinguishable? if 0, then not identical
  id_mat <- (outer(a,a,"==")) * 1
  
  #return index of identical sequences
  which(colSums(id_mat)>1)
}

amplicon_length_distribution <- function(results_list, minlength = 100, maxlength = 1000){
  f<-lapply(X = 1:length(b), FUN = function(X){b[[X]]$amplicon_size_no_primers})
  f <- unlist(f)
  f <- f[f>minlength & f<maxlength]
  return(f)
}
#amplicon_length_distribution(res, 100, 300) %>% hist()

#given alphas and a vector of reads, estimate true proportions
# estimate_proportions <- function(Nreads_vec, alphas, Npcr = 40){
#   ref_species <- which.max(alphas)
#   nu <- log(Nreads_vec/Nreads_vec[ref_species])
#   theta <- nu - (Npcr * alphas) 
#   p <- exp(theta)/sum(exp(theta))
#   return(p)
# }

get_alpha <- function(gamma, slope = -0.003, slope_se = 0.0006){
  #where gamma is the difference between observed deltaG and min deltaG for the primer/template pair
  return(list(
    est_alpha = gamma*slope,
    alpha_ci_lo = gamma*slope - 2*(sqrt(abs(gamma))*slope_se), #check SE calcs
    alpha_ci_hi = gamma*slope + 2*(sqrt(abs(gamma))*slope_se)
  ))
}
#get_alpha(gamma = -15)

estimate_proportions <- function(Nreads_vec, alphas, Npcr = 40, uncertainty = FALSE){
  require(dplyr)
  
  #select as reference species the one w the highest alpha and greatest number of reads
  ref_species <- data.frame(
    idx = 1:length(alphas),
    alphas, 
    Nreads_vec
  ) %>% 
    filter(Nreads_vec > 0) %>% 
    arrange(desc(alphas), desc(Nreads_vec)) %>% 
    pull(idx) %>% 
    nth(1)
  nu <- log(unlist(Nreads_vec)/unlist(Nreads_vec)[ref_species])
  p_df <- matrix(NA, nrow = 1000, ncol = length(alphas))
  if(!uncertainty){
    theta <- nu - (Npcr * alphas) 
    p_est <- exp(theta)/sum(exp(theta))
    return(list("p_est" = p_est, 
                "p_hi" = NA, 
                "p_lo" = NA))
    
  } else {
    for (i in 1:1000){
      tmp_alphas <- alphas + rnorm(length(alphas), 0, sd = 0.01) #adding estimated SD here
      tmp_alphas[ref_species] <- 0
      
      tmp_theta <- nu - (Npcr * tmp_alphas) 
      p_df[i,] <- exp(tmp_theta)/sum(exp(tmp_theta))
    }
    p_est <- colMeans(p_df)
    p_hi <- apply(p_df, MARGIN = 2, FUN = quantile, probs = 0.975)
    p_lo <- apply(p_df, MARGIN = 2, FUN = quantile, probs = 0.025)
    
    return(list("p_est" = p_est, 
                "p_hi" = p_hi, 
                "p_lo" = p_lo))
  }
}
# Nreads_vec <- c(5000, 300, 1000, 20, 700)
# estimate deltaG; find gammas
# gammas <-c(0, 15.526139, 24.810991, 41, 12)
# alphas <- get_alpha(gammas)
# estimate_proportions(Nreads_vec, alphas$est_alpha, Npcr = 40, uncertainty = FALSE)
# p <- estimate_proportions(Nreads_vec, alphas$est_alpha, Npcr = 40, uncertainty = TRUE)
# as.data.frame(p) %>% 
#   ggplot(aes(x = 1:5, y = p_est)) +
#   geom_point() +
#   geom_segment(aes(x = 1:5, xend = 1:5, y = p_lo, yend = p_hi))


# wrapper: given observations and info about species, predict correct proportions
correct_proportions <- function(Nreads_matrix, #m x n matrix of read counts, with species in rows and samples in columns
                                species_names, # vector of species names (in order of rows of Nreads_matrix)
                                species_templates, # vector of DNA template sequences against which to compare primers
                                forward_primer, # forward primer sequence, character vector
                                reverse_primer, # reverse primer sequence, character vector
                                Npcr = 40 # number of PCR cycles
){
  
  require(tidyverse)
  
  
  #get species binding energies for primer set
  
  binding_energies <- NA
  for (i in 1:length(species_names)){
    
    print(paste("Calculating binding energy for species", species_names[i]))
    
    binding_energies[i] <- predict_binding(forward_primer, reverse_primer, 
                                           template = as.character(species_templates[i]))$overall_binding_energy    
  }
  
  #calculate amp efficiency, given binding energies, relative to minimum binding energy observed
  
  gamma_vector <- binding_energies - min(binding_energies) #relative binding energy
  alpha_vector <- get_alpha(gamma_vector)$est_alpha #log-ratio amp efficiency; use means as point estimate
  
  # given alphas, estimate species' DNA proportions
  
  p_mat <- matrix(NA, nrow = nrow(Nreads_matrix), ncol = ncol(Nreads_matrix))
  p_mat_95ci_hi <- matrix(NA, nrow = nrow(Nreads_matrix), ncol = ncol(Nreads_matrix))
  p_mat_95ci_lo <- matrix(NA, nrow = nrow(Nreads_matrix), ncol = ncol(Nreads_matrix))
  
  
  for (i in 1:ncol(Nreads_matrix)){
    tmp <- estimate_proportions(Nreads_vec = as.vector(Nreads_matrix[,i]),
                                alphas = alpha_vector,
                                Npcr = Npcr,
                                uncertainty = T)
    
    p_mat[,i] <- tmp$p_est
    p_mat_95ci_hi[,i] <- tmp$p_hi
    p_mat_95ci_lo[,i] <- tmp$p_lo
  }
  
  tidy_estimates <- p_mat %>% 
    as.data.frame() %>%
    mutate(species = species_names) %>% 
    pivot_longer(cols = -species, names_to = "sample", values_to = "estimated_proportion") %>% 
    left_join(p_mat_95ci_lo %>% as.data.frame() %>% mutate(species = species_names) %>% 
                pivot_longer(cols = -species, names_to = "sample", values_to = "proportion_lo"), by = c("species", "sample")) %>% 
    left_join(p_mat_95ci_hi %>% as.data.frame() %>% mutate(species = species_names) %>% 
                pivot_longer(cols = -species, names_to = "sample", values_to = "proportion_hi"), by = c("species", "sample")) %>% 
    left_join(Nreads_matrix %>% 
                 as.data.frame() %>%
                 mutate(species = species_names) %>% 
                 pivot_longer(cols = -species, names_to = "sample", values_to = "observed_proportions") %>% 
                 group_by(sample) %>% 
                 mutate(observed_proportions = observed_proportions/sum(observed_proportions)), 
               by = c("species", "sample"))
  
  return(list(
    "estimated_proportions" = tidy_estimates,
    "binding_energies" = binding_energies,
    "est_amp_efficiencies" = alpha_vector
  ))
}




#TODO propagate uncertainty better from gamma through alpha to proportions 


#MFU primers
# forward_primer <- "GCCGGTAAAACTCGTGCCAGC"
# forward_primer_mm <- "GCCGGTAAAACTCGTAACAGC"
# reverse_primer <- "CATAGTGGGGTATCTAATCCCAGTTTG"
# #clupea
# template <- "GCCGGTAAAACTCGTGCCAGCCACCGCGGTTATACGAGAGACCCTAGTTGATATACTCGGCGTAAAGAGTGGTTATGGAAAACAAGCACTAAAGCCAAAGAGCCCTCAGGCCGTTATACGCACCCGGGGCCTCGAACCACTATCACGAAAGTAGCTTTACCCTCGCCCACCAGAACCCACGAGAGCTGGGACACAAACTGGGATTAGATACCCCACTATGCCCCGCCGTAAACTTAGATATATTAGTACAACAAATATCCGCCCGGGAACTACGAGCGCCAGCTTAAAACCCAAAGGACTTGGCGGTGCTTCAGACCCCCCTAGAGGAGCCTGTTCTAGAACCGATAACCCCCGTTCAACCTCACCACTCCTTGCCCCTCCCGCCTATATACCACCGTCGCCAGCTTACCCTGTGAAGGTACTACAGTAAGCAGAATGAGCATTCCTCAGAACGT"
# predict_binding(forward_primer_mm, reverse_primer, template)


# library(tidyverse)
library(here)
library(insect)


#example wrapper function
# taxbase <- read.csv(here("data/vert_mt_genomes/TaxonLookupTable.csv"))
# Nreads_matrix <- matrix(NA, nrow = 5, ncol = 10)
# for (i in 1:nrow(Nreads_matrix)){
#   Nreads_matrix[i,] <- rnbinom(ncol(Nreads_matrix), size = 1, mu = 1000)
# }
# species_names <- taxbase$species[1:5]
# species_templates <- readFASTA(here("data/vert_mt_genomes/vert_mt_genomes.fasta"), bin = F)[match(species_names, taxbase$species)]
# forward_primer <- "GCCGGTAAAACTCGTGCCAGC"
# reverse_primer <- "CATAGTGGGGTATCTAATCCCAGTTTG"
# 
# q <- correct_proportions(Nreads_matrix = Nreads_matrix,
#                          species_names = species_names,
#                          species_templates = species_templates,
#                          forward_primer = forward_primer,
#                          reverse_primer = reverse_primer,
#                          Npcr = 40
# )
# 
# #observed vs. estimated
# q$estimated_proportions %>% 
#   ggplot(aes(x = sample, y = estimated_proportion, color = species)) +
#   geom_point() +
#   geom_point(aes(x = sample, y = observed_proportions), color = "black") +
#   geom_segment(aes(x = sample, xend = sample, y = proportion_lo, yend = proportion_hi, color = species), alpha = 0.5) +
#   facet_wrap(~species)
# q$estimated_proportions %>% 
#   ggplot(aes(x = species, y = estimated_proportion, color = species)) +
#   geom_point() +
#   geom_point(aes(x = species, y = observed_proportions), color = "black") +
#   geom_segment(aes(x = species, xend = species, y = proportion_lo, yend = proportion_hi, color = species), alpha = 0.5) +
#   facet_wrap(~sample) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# 
