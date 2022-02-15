library(data.table)

library(dplyr)

library(limma)
library(edgeR)



samples <- list.files("$PATH/chromosome_y_expression/requant/primaryassembly/quants_PAR_masked")
samples <- samples[grepl("fctx",samples)]

counts_PAR_mask <- as.data.frame(fread("$PATH/chromosome_y_expression/requant/primaryassembly/quants_PAR_masked/quants_PAR_masked_matrix.csv"))
rownames(counts_PAR_mask) <- counts_PAR_mask$'Geneid'
print(dim(counts_PAR_mask))
print(counts_PAR_mask[1:5,1:5])

counts_default_ref <- as.data.frame(fread("$PATH/chromosome_y_expression/requant/primaryassembly/quants_default_ref/quants_default_ref_matrix.csv"))
rownames(counts_default_ref) <- counts_default_ref$'Geneid'
print(dim(counts_default_ref))
print(counts_default_ref[1:5,1:5])

#get the sample names from the pheno file
covs <- fread("$PATH/chromosome_y_expression/FCTX_NABEC/sample_info_new_id.txt")
print(dim(covs))
print(head(covs))



print(length(covs$new_id))
print(length(samples))

dim((covs[which(covs$new_id %in% samples)]))

print(unique(covs$Gender))

print(dim(covs[covs$Gender=="male"]))
print(dim(covs[covs$Gender=="female"]))

male_covs <- covs[(which(covs$new_id %in% samples & covs$Gender == "male")),]
print(dim(male_covs))
print(head(male_covs))


head(male_covs)

#get male columns 
male_cols <- colnames(counts_PAR_mask)[(colnames(counts_PAR_mask) %in% male_covs$new_id)]
#select male columns from counts
counts_PAR_mask_male <- counts_PAR_mask[,male_cols]
print(dim(counts_PAR_mask_male))

#get male columns 
male_cols <- colnames(counts_default_ref)[(colnames(counts_default_ref) %in% male_covs$new_id)]
#select male columns from counts
counts_default_male <- counts_default_ref[,male_cols]
print(dim(counts_default_male))

print(head(counts_PAR_mask_male))
colnames(counts_PAR_mask_male) <- paste0(colnames(counts_PAR_mask_male),"_PAR_mask")
#print(head(counts_PAR_mask_male))

print(head(counts_default_male))
colnames(counts_default_male) <- paste0(colnames(counts_default_male),"_default")
#print(head(counts_default_male))



final_exp <- merge(counts_PAR_mask_male, counts_default_male, by=0, all= TRUE)
rownames(final_exp) <- final_exp$'Row.names'
print(dim(final_exp))
#print(head(final_exp))
final_exp <- final_exp[,2:length(colnames(final_exp))]




mask_covs <- male_covs
mask_covs$new_id <- paste0(mask_covs$new_id,"_PAR_mask")
mask_covs$masked <- "PAR_masked"
def_covs <- male_covs
def_covs$new_id <- paste0(def_covs$new_id,"_default")
def_covs$masked <- "default"

final_covs <- rbind(mask_covs, def_covs)
rownames(final_covs) <- final_covs$new_id
print(dim(final_covs))
print(head(final_covs))
print(tail(final_covs))

#remove rows with all zeros
print("how many after removing all zero rows")
print(dim(final_exp[rowSums(final_exp[])>0,]))
final_exp <- final_exp[rowSums(final_exp[])>0,]



#use dplyr select to reorder columns 
final_exp <- final_exp %>% dplyr::select(final_covs$new_id)
print(dim(final_exp))
print(dim(final_covs))

#check if colnames of cntTable match rownames of demogdiag metadata and if they are in the same order
print("check if data cols are same as meta data rows")
print(all(colnames(final_exp) == rownames(final_covs)))


##edgeR DE with samples unpaired
print("starting sample unpaired DE analysis")
dge <- DGEList(counts=final_exp, samples = final_covs, group = final_covs$masked)

design <- model.matrix(~0+group,data = dge$samples)

dge <- calcNormFactors(dge)

dge <- estimateDisp(dge, design, robust=TRUE)

fit <- glmQLFit(dge, design, robust=TRUE)
print(head(fit$coefficients))

contr.matrix <- makeContrasts(maskedvsdefault = groupPAR_masked-groupdefault,levels = colnames(design))
print(contr.matrix)

res <- glmQLFTest(fit, contrast=contr.matrix)

topTagGene <- topTags(res,n=Inf)
print(head(topTagGene$table))

write.csv(topTagGene$table,"$PATH/chromosome_y_expression/requant/primaryassembly/diff_exp_all_genes/diff_exp_edgeR_unpaired.csv",row.names = FALSE)

print("results for unpaired edgeR DE:")
is.de <- decideTestsDGE(res)
print(summary(is.de))


##edgeR DE with samples paired
print("starting sample paired DE analysis")
dge <- DGEList(counts=final_exp, samples = final_covs, group = final_covs$masked)

subject <- factor(dge$samples$SampleId)
design <- model.matrix(~subject+group,data = dge$samples)

dge <- calcNormFactors(dge)

dge <- estimateDisp(dge, design, robust=TRUE)

fit <- glmQLFit(dge, design, robust=TRUE)
print(head(fit$coefficients))

contr.matrix <- makeContrasts(maskedvsdefault = groupPAR_masked-groupdefault,levels = colnames(design))
print(contr.matrix)

res <- glmQLFTest(fit, contrast=contr.matrix)

topTagGene <- topTags(res,n=Inf)
print(head(topTagGene$table))

write.csv(topTagGene$table,"$PATH/chromosome_y_expression/requant/primaryassembly/diff_exp_all_genes/diff_exp_edgeR_paired.csv",row.names = FALSE)

print("results for paired edgeR DE:")
is.de <- decideTestsDGE(res)
print(summary(is.de))


