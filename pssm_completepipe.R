#utilities
library(PSSMCOOL)
library(stringr)

setwd(r'(D:\Lab_work\Blast_work\output_pssm)')
pssm_directory <- r'(D:\Lab_work\Blast_work\output_pssm)'

pssm_files <- list.files(path = pssm_directory, pattern = "\\.pssm$",full.names = TRUE) 


file_names <- tools::file_path_sans_ext(basename(pssm_files))

combined <- data.frame()

col_names <- c()


Amino_vec <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")


#AAC header
aac <- paste0("AAC_PSSM_", Amino_vec)

# Create the AADP_AAC_PSSM part
aadp_aac <- paste0("AADP_AAC_PSSM_", Amino_vec)
# Create the AADP_DPC_PSSM part using expand.grid for pairwise combinations
aadp_dpc <- paste0("AADP_DPC_PSSM_", expand.grid(Amino_vec, Amino_vec)[, 1], "_", expand.grid(Amino_vec, Amino_vec)[, 2])
# AADP header
aadp <- c(aadp_aac, aadp_dpc)

# Create the AATP_TPC_AAC_PSSM part
aatp_tpc_aac <- paste0("AATP_TPC_PSSM_", Amino_vec)
# Create the AATP_TPC_DPC_PSSM part using expand.grid for pairwise combinations
aatp_tpc_dpc <- paste0("AATP_TPC_PSSM_", expand.grid(Amino_vec, Amino_vec)[, 1], "_", expand.grid(Amino_vec, Amino_vec)[, 2])
# AATP header
aatp <- c(aatp_tpc_aac, aatp_tpc_dpc)



# Generate combinations of blocks and amino acids
blocks <- 1:20
#AB header
ab <- paste0("Block", rep(blocks, each = length(Amino_vec)), "_", Amino_vec)

#DFPSSM header
dfpssm<- paste0("D_FPSSM_", Amino_vec)


#DP_PSSM header
dp_pssm <- c()

# Generate headers for positive and negative averages
for (aa in Amino_vec) {
  dp_pssm <- c(dp_pssm, paste0("DP_PSSM_PosAvg_", aa))
  dp_pssm <- c(dp_pssm, paste0("DP_PSSM_NegAvg_", aa))
}



# Generate headers for squared averages of differences
for (k in 1:5) {
  for (aa in Amino_vec) {
    dp_pssm <- c(dp_pssm, paste0("DP_PSSM_PosDiffSqAvg_k", k, "_", aa))
    dp_pssm <- c(dp_pssm, paste0("DP_PSSM_NegDiffSqAvg_k", k, "_", aa))
  }
}

#DPC_PSSM headers
dpc <- paste0("DPC_PSSM_", expand.grid(Amino_vec, Amino_vec)[, 1], "_", expand.grid(Amino_vec, Amino_vec)[, 2])

#edp
edp <- paste0("EDP_PSSM_", Amino_vec)

#eedp
eedp <- paste0("EEDP_PSSM_", expand.grid(Amino_vec, Amino_vec)[, 1], "_", expand.grid(Amino_vec, Amino_vec)[, 2])


#k_separated_bigrams_pssm
ksbi <- paste0("K_SEP_BIGRAM_PSSM_", expand.grid(Amino_vec, Amino_vec)[, 1], "_", expand.grid(Amino_vec, Amino_vec)[, 2])

#medp
medp_edp <- paste0("MEDP_PSSM_", Amino_vec)
# Create the MEDP_EEDP_PSSM part using expand.grid for pairwise combinations
medp_eedp <- paste0("MEDP_PSSM_", expand.grid(Amino_vec, Amino_vec)[, 1], "_", expand.grid(Amino_vec, Amino_vec)[, 2])
# MEDP header
medp <- c(medp_edp, medp_eedp)


fp <- paste0("PSE_PSSM_FPSSM_", Amino_vec)
# Generate headers for correlation factors
cor <- paste0("PSE_PSSM_CorrFact_1_", Amino_vec)
#PSE_PSSM headers
pse <- c(fp,cor)

#pssm_ac header
ac <- c()
for (g in 1:10) {  # g from 1 to 10
  for (aa in Amino_vec) {
    ac <- c(ac, paste0("PSSM_AC_g", g, "_", aa))
  }
}


#pssm_cc header
cc <- c()

# Set LG (the distance parameter)
lg <- 10  # You can change this value as needed

# Generate headers for PSSM_CC features
for (j1 in 1:20) {  # First amino acid column (1 to 20)
  for (j2 in 1:20) {  # Second amino acid column (1 to 20)
    if (j1 == j2) {
      next  # Skip the case where j1 == j2
    }
    for (g in 1:lg) {  # g from 1 to LG
      header <- paste0("PSSM_CC_g", g, "_", Amino_vec[j1], "_", Amino_vec[j2])
      cc <- c(cc, header)  # Append the header to the list
    }
  }
}

#pssm_composition
composition <- c()
for (g in 1:20) {  # g from 1 to 10
  for (aa in Amino_vec) {
    composition <- c(composition, paste0("PSSM_COMPOSITION_g", g, "_", aa))
  }
}


#rpm_pssm
rpm<-c()
for (p in 1:20) {  # g from 1 to 10
  for (aa in Amino_vec) {
    rpm <- c(rpm, paste0("RPM_Probe_",p , "_", aa))
  }
}

#rpssm

# Define amino acid groups and associated columns
grouped_columns <- list(
  list(c("F", "Y", "W"), c(5, 18, 19)),   # P1 = PF + PY + PW
  list(c("M", "L"), c(12, 11)),           # P2 = PM + PL
  list(c("I", "V"), c(9, 20)),            # P3 = PI + PV
  list(c("A", "T", "S"), c(1, 17, 16)),   # P4 = PA + PT + PS
  list(c("N", "H"), c(14, 8)),            # P5 = PN + PH
  list(c("Q", "E", "D"), c(15, 7, 4)),    # P6 = PQ + PE + PD
  list(c("R", "K"), c(2, 13)),            # P7 = PR + PK
  list(c("C"), c(6)),                     # P8 = PC
  list(c("G"), c(10)),                    # P9 = PG
  list(c("P"), c(3))                      # P10 = PP
)

# Step 1: Generate headers for feature vector of length 10 (Ds values)
headers_10 <- sapply(1:length(grouped_columns), function(s) {
  aas <- paste(grouped_columns[[s]][[1]], collapse = "+")
  paste0("RPSSM_", s, "_", aas)
})

# Step 2: Generate headers for feature vector of length 100 (Ds,t values)
headers_100 <- c()
for (s in 1:length(grouped_columns)) {
  aas_s <- paste(grouped_columns[[s]][[1]], collapse = "+")
  for (t in 1:length(grouped_columns)) {
    aas_t <- paste(grouped_columns[[t]][[1]], collapse = "+")
    headers_100 <- c(headers_100, paste0("RPSSM_", s, "_", aas_s, "_", t, "_", aas_t))
  }
}

# Combine headers
rpssm <- c(headers_10, headers_100)


#tpc_pssm
tpc <- paste0("TPC_PSSM_", expand.grid(Amino_vec, Amino_vec)[, 1], "_", expand.grid(Amino_vec, Amino_vec)[, 2]) 

#trigram
trigram <- c()
for (aa1 in Amino_vec) {
  for (aa2 in Amino_vec) {
    for (aa3 in Amino_vec) {
      trigram <- c(trigram, paste("PSSM_Triple", aa1, aa2, aa3, sep = "_"))
    }
  }
}



col_names <- c(aac,aadp,aatp,ab,dfpssm,dp_pssm,dpc,edp,eedp,ksbi,medp,pse,ac,cc,composition,rpm,rpssm,tpc,trigram)
length(col_names)



for (i in seq_len(length(pssm_files))){
  
  row_data <- c(aac_pssm(pssm_files[i]))
  row_data <- c(row_data,aadp_pssm(pssm_files[i]))
  row_data <- c(row_data,AATP_TPC(pssm_files[i])[[2]])
  row_data <- c(row_data,AB_PSSM(pssm_files[i]))
  row_data <- c(row_data,FPSSM(pssm_files[i])[[1]])
  row_data <- c(row_data,DP_PSSM(pssm_files[1]))
  row_data <- c(row_data,dpc_pssm(pssm_files[i]))
  row_data <- c(row_data,EDP_EEDP_MEDP(pssm_files[i])[[1]])
  row_data <- c(row_data,EDP_EEDP_MEDP(pssm_files[i])[[2]])
  row_data <- c(row_data,k_separated_bigrams_pssm(pssm_files[i]))
  row_data <- c(row_data,EDP_EEDP_MEDP(pssm_files[i])[[3]])
  row_data <- c(row_data,pse_pssm(pssm_files[i]))
  row_data <- c(row_data,pssm_ac(pssm_files[i]))
  row_data <- c(row_data,pssm_cc(pssm_files[i]))
  row_data <- c(row_data,pssm_composition(pssm_files[i]))
  row_data <- c(row_data,RPM_PSSM(pssm_files[i]))
  row_data <- c(row_data,rpssm(pssm_files[i]))
  row_data <- c(row_data,AATP_TPC(pssm_files[i])[[1]])
  row_data <- c(row_data,trigrame_pssm(pssm_files[i]))
  
  combined <- rbind(combined,as.data.frame(t(row_data)))
  
}

df <- data.frame(ID = file_names, combined)
colnames(df) <- c('ID',col_names)


write.csv(df, "outputR.csv", row.names = FALSE)


