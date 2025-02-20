




seq_files <- list.files(pattern = "seq_")

comb_seqs <- data.frame()
for (i in seq_files) {
	seqs <- read.csv(i)
	comb_seqs <- rbind(comb_seqs, seqs)
}

write.csv(comb_seqs, "sir_res.csv")
