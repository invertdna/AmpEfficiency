
#Core and helper functions for estimating DNA proportions from amplicon sequencing data
#implementing binding-energy model to estimate amplification efficiency biases and correct for them
#RPK March/April 2025

#Dinucleotide Binding Energy Matrix (ΔG37 in kcal/mol)
create_binding_energy_matrix <- function() {
  dinucleotides <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
                     "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
  
  # Create empty 16x16 matrix with named rows and columns
  matrix_data <- matrix(0, nrow = 16, ncol = 16)
  rownames(matrix_data) <- dinucleotides
  colnames(matrix_data) <- dinucleotides
  
  # Fill in the matrix with the corrected ΔG37 values (kcal/mol)
  matrix_data["AA", "AA"] <- 0.91; matrix_data["AA", "AC"] <- 1.50; matrix_data["AA", "AG"] <- 1.20; matrix_data["AA", "AT"] <- 0.70
  matrix_data["AA", "CA"] <- 1.50; matrix_data["AA", "CC"] <- 1.50; matrix_data["AA", "CG"] <- 1.20; matrix_data["AA", "CT"] <- 0.70
  matrix_data["AA", "GA"] <- 1.20; matrix_data["AA", "GC"] <- 1.20; matrix_data["AA", "GG"] <- 1.20; matrix_data["AA", "GT"] <- 0.70
  matrix_data["AA", "TA"] <- 0.70; matrix_data["AA", "TC"] <- 0.70; matrix_data["AA", "TG"] <- 0.70; matrix_data["AA", "TT"] <- -1.00
  
  matrix_data["AC", "AA"] <- 1.50; matrix_data["AC", "AC"] <- 1.50; matrix_data["AC", "AG"] <- 0.70; matrix_data["AC", "AT"] <- 1.20
  matrix_data["AC", "CA"] <- 0.83; matrix_data["AC", "CC"] <- 1.50; matrix_data["AC", "CG"] <- 0.70; matrix_data["AC", "CT"] <- 1.20
  matrix_data["AC", "GA"] <- 1.20; matrix_data["AC", "GC"] <- 1.20; matrix_data["AC", "GG"] <- 0.70; matrix_data["AC", "GT"] <- 0.65
  matrix_data["AC", "TA"] <- 0.70; matrix_data["AC", "TC"] <- 0.70; matrix_data["AC", "TG"] <- -1.50; matrix_data["AC", "TT"] <- 0.70
  
  matrix_data["AG", "AA"] <- 1.20; matrix_data["AG", "AC"] <- 0.70; matrix_data["AG", "AG"] <- 1.50; matrix_data["AG", "AT"] <- 0.90
  matrix_data["AG", "CA"] <- 1.20; matrix_data["AG", "CC"] <- 0.70; matrix_data["AG", "CG"] <- 1.50; matrix_data["AG", "CT"] <- 0.90
  matrix_data["AG", "GA"] <- 0.08; matrix_data["AG", "GC"] <- 0.70; matrix_data["AG", "GG"] <- 1.20; matrix_data["AG", "GT"] <- 0.90
  matrix_data["AG", "TA"] <- 0.70; matrix_data["AG", "TC"] <- -1.50; matrix_data["AG", "TG"] <- 0.70; matrix_data["AG", "TT"] <- 0.40
  
  matrix_data["AT", "AA"] <- 0.70; matrix_data["AT", "AC"] <- 1.20; matrix_data["AT", "AG"] <- 0.90; matrix_data["AT", "AT"] <- 1.50
  matrix_data["AT", "CA"] <- 0.70; matrix_data["AT", "CC"] <- 1.20; matrix_data["AT", "CG"] <- 0.90; matrix_data["AT", "CT"] <- 1.50
  matrix_data["AT", "GA"] <- 0.70; matrix_data["AT", "GC"] <- 1.20; matrix_data["AT", "GG"] <- 0.90; matrix_data["AT", "GT"] <- 1.20
  matrix_data["AT", "TA"] <- -0.88; matrix_data["AT", "TC"] <- 0.70; matrix_data["AT", "TG"] <- 0.40; matrix_data["AT", "TT"] <- 0.70
  
  matrix_data["CA", "AA"] <- 1.50; matrix_data["CA", "AC"] <- 0.75; matrix_data["CA", "AG"] <- 1.20; matrix_data["CA", "AT"] <- 0.70
  matrix_data["CA", "CA"] <- 1.50; matrix_data["CA", "CC"] <- 1.50; matrix_data["CA", "CG"] <- 1.20; matrix_data["CA", "CT"] <- 0.70
  matrix_data["CA", "GA"] <- 0.70; matrix_data["CA", "GC"] <- 0.70; matrix_data["CA", "GG"] <- 0.70; matrix_data["CA", "GT"] <- -1.45
  matrix_data["CA", "TA"] <- 1.20; matrix_data["CA", "TC"] <- 1.20; matrix_data["CA", "TG"] <- 0.58; matrix_data["CA", "TT"] <- 0.70
  
  matrix_data["CC", "AA"] <- 1.50; matrix_data["CC", "AC"] <- 1.50; matrix_data["CC", "AG"] <- 0.70; matrix_data["CC", "AT"] <- 1.20
  matrix_data["CC", "CA"] <- 1.50; matrix_data["CC", "CC"] <- 1.05; matrix_data["CC", "CG"] <- 0.70; matrix_data["CC", "CT"] <- 1.20
  matrix_data["CC", "GA"] <- 0.70; matrix_data["CC", "GC"] <- 0.70; matrix_data["CC", "GG"] <- -1.50; matrix_data["CC", "GT"] <- 0.70
  matrix_data["CC", "TA"] <- 1.20; matrix_data["CC", "TC"] <- 1.20; matrix_data["CC", "TG"] <- 0.70; matrix_data["CC", "TT"] <- 1.20
  
  matrix_data["CG", "AA"] <- 1.20; matrix_data["CG", "AC"] <- 0.70; matrix_data["CG", "AG"] <- 1.50; matrix_data["CG", "AT"] <- 0.90
  matrix_data["CG", "CA"] <- 1.20; matrix_data["CG", "CC"] <- 0.70; matrix_data["CG", "CG"] <- 1.50; matrix_data["CG", "CT"] <- 0.90
  matrix_data["CG", "GA"] <- 0.70; matrix_data["CG", "GC"] <- -2.17; matrix_data["CG", "GG"] <- 0.70; matrix_data["CG", "GT"] <- 0.40
  matrix_data["CG", "TA"] <- 1.20; matrix_data["CG", "TC"] <- 0.70; matrix_data["CG", "TG"] <- 1.20; matrix_data["CG", "TT"] <- 0.90
  
  matrix_data["CT", "AA"] <- 0.70; matrix_data["CT", "AC"] <- 1.20; matrix_data["CT", "AG"] <- 0.90; matrix_data["CT", "AT"] <- 1.50
  matrix_data["CT", "CA"] <- 0.70; matrix_data["CT", "CC"] <- 1.20; matrix_data["CT", "CG"] <- 0.90; matrix_data["CT", "CT"] <- 1.50
  matrix_data["CT", "GA"] <- -1.28; matrix_data["CT", "GC"] <- 0.70; matrix_data["CT", "GG"] <- 0.40; matrix_data["CT", "GT"] <- 0.70
  matrix_data["CT", "TA"] <- 0.70; matrix_data["CT", "TC"] <- 0.72; matrix_data["CT", "TG"] <- 0.90; matrix_data["CT", "TT"] <- 1.20
  
  matrix_data["GA", "AA"] <- 1.20; matrix_data["GA", "AC"] <- 1.20; matrix_data["GA", "AG"] <- 0.25; matrix_data["GA", "AT"] <- 0.70
  matrix_data["GA", "CA"] <- 0.70; matrix_data["GA", "CC"] <- 0.70; matrix_data["GA", "CG"] <- 0.70; matrix_data["GA", "CT"] <- -1.30
  matrix_data["GA", "GA"] <- 1.50; matrix_data["GA", "GC"] <- 1.50; matrix_data["GA", "GG"] <- 1.20; matrix_data["GA", "GT"] <- 0.70
  matrix_data["GA", "TA"] <- 0.90; matrix_data["GA", "TC"] <- 0.90; matrix_data["GA", "TG"] <- 0.90; matrix_data["GA", "TT"] <- 0.40
  
  matrix_data["GC", "AA"] <- 1.20; matrix_data["GC", "AC"] <- 1.20; matrix_data["GC", "AG"] <- 0.70; matrix_data["GC", "AT"] <- 1.20
  matrix_data["GC", "CA"] <- 0.70; matrix_data["GC", "CC"] <- 0.70; matrix_data["GC", "CG"] <- -2.24; matrix_data["GC", "CT"] <- 0.70
  matrix_data["GC", "GA"] <- 1.50; matrix_data["GC", "GC"] <- 1.50; matrix_data["GC", "GG"] <- 0.70; matrix_data["GC", "GT"] <- 1.20
  matrix_data["GC", "TA"] <- 0.90; matrix_data["GC", "TC"] <- 0.90; matrix_data["GC", "TG"] <- 0.40; matrix_data["GC", "TT"] <- 0.90
  
  matrix_data["GG", "AA"] <- 1.20; matrix_data["GG", "AC"] <- 0.70; matrix_data["GG", "AG"] <- 1.20; matrix_data["GG", "AT"] <- 0.90
  matrix_data["GG", "CA"] <- 0.70; matrix_data["GG", "CC"] <- -1.84; matrix_data["GG", "CG"] <- 0.70; matrix_data["GG", "CT"] <- 0.40
  matrix_data["GG", "GA"] <- 1.20; matrix_data["GG", "GC"] <- 0.70; matrix_data["GG", "GG"] <- 0.59; matrix_data["GG", "GT"] <- 0.90
  matrix_data["GG", "TA"] <- 0.90; matrix_data["GG", "TC"] <- 0.40; matrix_data["GG", "TG"] <- 0.90; matrix_data["GG", "TT"] <- 1.07
  
  matrix_data["GT", "AA"] <- 0.70; matrix_data["GT", "AC"] <- 1.20; matrix_data["GT", "AG"] <- 0.90; matrix_data["GT", "AT"] <- 1.20
  matrix_data["GT", "CA"] <- -1.44; matrix_data["GT", "CC"] <- 0.70; matrix_data["GT", "CG"] <- 0.40; matrix_data["GT", "CT"] <- 0.70
  matrix_data["GT", "GA"] <- 0.70; matrix_data["GT", "GC"] <- 1.20; matrix_data["GT", "GG"] <- 0.90; matrix_data["GT", "GT"] <- 1.50
  matrix_data["GT", "TA"] <- 0.40; matrix_data["GT", "TC"] <- 0.90; matrix_data["GT", "TG"] <- 0.43; matrix_data["GT", "TT"] <- 0.90
  
  matrix_data["TA", "AA"] <- 0.70; matrix_data["TA", "AC"] <- 0.70; matrix_data["TA", "AG"] <- 0.70; matrix_data["TA", "AT"] <- -0.58
  matrix_data["TA", "CA"] <- 1.20; matrix_data["TA", "CC"] <- 1.20; matrix_data["TA", "CG"] <- 1.20; matrix_data["TA", "CT"] <- 0.70
  matrix_data["TA", "GA"] <- 0.90; matrix_data["TA", "GC"] <- 0.90; matrix_data["TA", "GG"] <- 0.90; matrix_data["TA", "GT"] <- 0.40
  matrix_data["TA", "TA"] <- 1.50; matrix_data["TA", "TC"] <- 1.50; matrix_data["TA", "TG"] <- 1.20; matrix_data["TA", "TT"] <- 0.70
  
  matrix_data["TC", "AA"] <- 0.70; matrix_data["TC", "AC"] <- 0.70; matrix_data["TC", "AG"] <- -1.50; matrix_data["TC", "AT"] <- 0.70
  matrix_data["TC", "CA"] <- 1.20; matrix_data["TC", "CC"] <- 1.20; matrix_data["TC", "CG"] <- 0.70; matrix_data["TC", "CT"] <- 0.88
  matrix_data["TC", "GA"] <- 0.90; matrix_data["TC", "GC"] <- 0.90; matrix_data["TC", "GG"] <- 0.40; matrix_data["TC", "GT"] <- 0.90
  matrix_data["TC", "TA"] <- 1.50; matrix_data["TC", "TC"] <- 1.50; matrix_data["TC", "TG"] <- 0.70; matrix_data["TC", "TT"] <- 1.20
  
  matrix_data["TG", "AA"] <- 0.70; matrix_data["TG", "AC"] <- -1.50; matrix_data["TG", "AG"] <- 0.70; matrix_data["TG", "AT"] <- 0.40
  matrix_data["TG", "CA"] <- 1.20; matrix_data["TG", "CC"] <- 0.70; matrix_data["TG", "CG"] <- 1.20; matrix_data["TG", "CT"] <- 0.90
  matrix_data["TG", "GA"] <- 0.90; matrix_data["TG", "GC"] <- 0.40; matrix_data["TG", "GG"] <- 0.90; matrix_data["TG", "GT"] <- 0.34
  matrix_data["TG", "TA"] <- 1.20; matrix_data["TG", "TC"] <- 0.70; matrix_data["TG", "TG"] <- 1.50; matrix_data["TG", "TT"] <- 0.90
  
  matrix_data["TT", "AA"] <- -1.50; matrix_data["TT", "AC"] <- 0.70; matrix_data["TT", "AG"] <- 0.40; matrix_data["TT", "AT"] <- 0.70
  matrix_data["TT", "CA"] <- 0.70; matrix_data["TT", "CC"] <- 1.20; matrix_data["TT", "CG"] <- 0.90; matrix_data["TT", "CT"] <- 1.20
  matrix_data["TT", "GA"] <- 0.40; matrix_data["TT", "GC"] <- 0.90; matrix_data["TT", "GG"] <- 0.98; matrix_data["TT", "GT"] <- 0.90
  matrix_data["TT", "TA"] <- 0.70; matrix_data["TT", "TC"] <- 1.20; matrix_data["TT", "TG"] <- 0.90; matrix_data["TT", "TT"] <- 0.86
  
  return(matrix_data)
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

weight_function <- function(x, primer_length, rev = FALSE, max_weight = 5){
  # x is the distance from the 3' end of primer
  # primer_length is the length of the primer
  # rev is a logical indicating if the primer being evaluated in reverse direction
  # returns a weight for the primer at position x
  a <- -(primer_length - 5)
  b <- 1
  
  if (rev){
    w <- primer_length - x
  } else {
    w <- x
  }
  # sigmoidal weighting function
  return(
    1 + ((max_weight-1)*exp(a+b*w)/
           (1+exp(a+b*w)))
  ) 
}

get_binding_energy <- function(p1, t1, rev = FALSE) {
  # NOTE: handles ambiguities in the primer, but not in the template
  
  comp_template <- complement(t1)
  p1 <- expand_primer(p1)
  pvec <- NA #vector for temporarily storing values
  for (j in 1:length(p1)){ #to handle cases in which primer ambiguities mean a mixture of primer seqs
    deltaG_raw <- 0
    deltaG_cumul <- 0
    deltaG_adjusted <- 0
    for (i in 1:(nchar(p1[j]) - 1)) {
      dinuc1 <- substr(p1[j], i, i + 1)
      dinuc2 <- substr(comp_template, i, i + 1)
      deltaG_raw <- m[dinuc1, dinuc2]
      deltaG_adjusted <- ifelse(deltaG_raw < 0, deltaG_raw, 
                                deltaG_raw * weight_function(i, nchar(p1[j]), rev))  #weighted to penalize mismatches toward the 3' end; assumes 5' - 3' orientation when rev = FALSE (default)
      deltaG_cumul <- sum(c(deltaG_cumul, deltaG_adjusted))
    }
    pvec[j] <- deltaG_cumul
  }
  return(mean(pvec)) #return mean across primer mixture, if applicable
  #note binding energy for initiation is added in the two-primer function, `predict_binding()` below; could put it in here.
  }
  
m <- create_binding_energy_matrix() #units kcal/mol at 37C

#TODO
# separate delta H and delta S components, so as to handle variable salt and MgCl2 conditions
# check against examples in SantaLucia 

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
    
  
  
  #iterate over forward primer positions
  g_vec_fwd <- NA
  for (pos in 1:(nchar(template)-nchar(forward_primer))) {
    tmp <- substr(template, pos, (pos+nchar(forward_primer)-1))
    g_vec_fwd[pos] <- get_binding_energy(forward_primer, tmp)
  }
  
  #iterate over reverse primer positions, indexing from the leftward-most position of the binding position
  #template in original 5' -> 3' orientation; rev primer in reverse complement
  g_vec_rev <- NA
  revcomp_reverse_primer <- reverse_complement(reverse_primer)
  for (pos in 1:(nchar(template)-nchar(reverse_primer)+1)) {
    tmp <- substr(template, pos, (pos+nchar(reverse_primer)-1))
    g_vec_rev[pos] <- get_binding_energy(revcomp_reverse_primer, tmp, rev = TRUE)
  }
  
  #initiation energy
  deltaG_init = 3.4 #kcal/mol (or approximately +14.2 kJ/mol)
  #SantaLucia gives initiation with G-C as 2.6, A-T as 1.7, and presumably these are per-primer. But then his fig 1 gives 1.03 as A-T initiation in example.
  
  
  #output
  best_forward_binding_position <- which(g_vec_fwd == min(g_vec_fwd)) #keeps both in case of ties
  best_reverse_binding_position <- which(g_vec_rev == min(g_vec_rev))
  overall_binding_energy <- g_vec_fwd[best_forward_binding_position[1]] + g_vec_rev[best_reverse_binding_position[1]] + deltaG_init
  amplicon_size_w_primers <- ((best_reverse_binding_position + nchar(reverse_primer)) - 
                                (best_forward_binding_position)) 
  amplicon_size_no_primers <- amplicon_size_w_primers - nchar(forward_primer) - nchar(reverse_primer)
  amplicon_w_primers <- substr(template, best_forward_binding_position[1], 
                               (best_reverse_binding_position[1] + nchar(reverse_primer) - 1)) #only keeps first, if there are multiple
  amplicon_no_primers <- substr(template, (best_forward_binding_position[1] + 1), 
                                (best_reverse_binding_position[1] - 1)) #only keeps first, if there are multiple
  
  return(list("best_forward_binding_position" = best_forward_binding_position,
              "best_reverse_binding_position" = best_reverse_binding_position,
              "overall_binding_energy" = overall_binding_energy,
              "amplicon_size_w_primers" = amplicon_size_w_primers,
              "amplicon_size_no_primers" = amplicon_size_no_primers,
              "amplicon_w_primers" = amplicon_w_primers,
              "amplicon_no_primers" = amplicon_no_primers))
  
  }

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
estimate_proportions <- function(Nreads_vec, alphas, Npcr = 40){
  ref_species <- which.max(alphas)
  nu <- log(Nreads_vec/Nreads_vec[ref_species])
  theta <- nu - (Npcr * alphas) 
  p <- exp(theta)/sum(exp(theta))
  return(p)
}

get_alpha <- function(gamma, slope = -0.0031, slope_se = 0.0006){
  #where gamma is the difference between observed deltaG and min deltaG for the primer/template pair
  return(list(
    est_alpha = gamma*slope,
    alpha_ci_lo = gamma*slope - 2*(sqrt(abs(gamma))*slope_se), #check SE calcs
    alpha_ci_hi = gamma*slope + 2*(sqrt(abs(gamma))*slope_se)
  ))
}
#get_alpha(gamma = -15)

estimate_proportions <- function(Nreads_vec, alphas, Npcr = 40, uncertainty = FALSE){
  ref_species <- which.max(alphas)
  nu <- log(Nreads_vec/Nreads_vec[ref_species])
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
    
    binding_energies[i] <- predict_binding(forward_primer, reverse_primer, template = as.character(species_templates[i]))$overall_binding_energy    
  }
  
  #calculate amp efficiency, given binding energies, relative to minimum binding energy observed
  
  gamma_vector <- binding_energies - min(binding_energies) #relative binding energy
  alpha_vector <- get_alpha(gamma_vector)$est_alpha #log-ratio amp efficiency; use means as point estimate
  
  # given alphas, estimate species' DNA proportions
  
  p_mat <- matrix(NA, nrow = nrow(Nreads_matrix), ncol = ncol(Nreads_matrix))
  p_mat_95ci_hi <- matrix(NA, nrow = nrow(Nreads_matrix), ncol = ncol(Nreads_matrix))
  p_mat_95ci_lo <- matrix(NA, nrow = nrow(Nreads_matrix), ncol = ncol(Nreads_matrix))
  
  
  for (i in 1:ncol(Nreads_matrix)){
    tmp <- estimate_proportions(Nreads_vec = Nreads_matrix[,i],
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
