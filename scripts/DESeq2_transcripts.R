log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# library("tximportData")
library("tximport")
library("readr")
library("DESeq2")


files <- snakemake@input[["counts"]]
names(files) <- snakemake@params[["names"]]
kallisto_counts <- tximport(files, type = "kallisto", txOut = TRUE)
# print(kallisto_counts$counts)

counts = kallisto_counts$counts

counts <- as.matrix(counts)
mode(counts) <- "integer"

# Loading samples information from file.
all_samples <- read.table(
    snakemake@input[["samples"]], header=TRUE,
    row.names="sample", check.names=FALSE
)

dir.create(snakemake@output[["results"]], showWarnings=FALSE)
# Looping through the samples
# conditions = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
# for (i in 1:9) {
# 	for (j in (i+1):10) {

conditions = c("A", "B") #JUST FOR TEST!!!!!!!!!!!!!!!!!!!
for (i in 1:1) {
	for (j in (i+1):2) {

		cnd1 = conditions[i]
		cnd2 = conditions[j]

		exp <- sprintf("%s-%s", cnd1, cnd2)
		base <- sprintf("%s", cnd1)
		other <- sprintf("%s", cnd2)

		# Slicing
		samples <- subset(all_samples, condition==base | condition==other)
		count <- counts[,c(row.names(samples))]


		# Calculating DESeq2
		dds <- DESeqDataSetFromMatrix(
	    		countData=count,
	    		colData=samples,
	    		design= ~condition
		)


		dds$condition <- relevel(dds$condition, ref=base)
		dds <- DESeq(dds)
		results = results(dds, contrast=c("condition", cnd1, cnd2))

	    	# Writing results to file
		fname <- paste(
			snakemake@output[["results"]],
			paste(exp, "csv", sep='.'), sep='/'
		)

		write.csv(
		        as.data.frame(results),
		        file=fname,
		        quote=FALSE
		)
	}
}
