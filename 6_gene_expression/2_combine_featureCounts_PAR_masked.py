import pandas as pd
import os
counts = pd.DataFrame()
OUTDIR = "/PATH/primaryassembly"
WRKDIR = "/PATH/chrY"
NABECDIR = "/PATH/nabec"
#get the haplogroup data, which should have the sample ids and should all be male
haplos = pd.read_csv(f"{WRKDIR}/output_nabec/nabec_haplos.csv")
haplos['new_id'] = haplos['new_id']+'fctx'
print(haplos.head())


#get the samples that we have fastqs for
fastqs = os.listdir(f"{NABECDIR}/fastqs")
print(len(fastqs))
samples = list(set([s.replace("_R2.fastq.gz","").replace("_R1.fastq.gz","") for s in fastqs]))
print(len(samples))
print(samples[0:10])

#combine these to get the samples we want to run through salmon
samples_to_use = set(samples).intersection(set(haplos.new_id))
print(len(samples_to_use))

for s in samples_to_use:
    print(s)
    df = pd.read_table(f"{OUTDIR}/quants_PAR_masked/{s}/{s}.featureCounts.tsv",sep = "\s+",skiprows = 1)
    #print(df.shape)
    #print(df.head())
    df = df.iloc[:,[0,6]]
    df.columns = ["Geneid", s]
    if(len(counts.index)==0):
        counts = df
    else:
        counts = pd.merge(left = counts, right = df, left_on = "Geneid", right_on = "Geneid")
    
    #if(len(counts.columns) >4):
        #break;
counts.to_csv(f"{OUTDIR}/quants_PAR_masked/quants_PAR_masked_matrix.csv",index=None)