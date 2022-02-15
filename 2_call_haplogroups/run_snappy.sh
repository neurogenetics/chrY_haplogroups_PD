#!/bin/bash

module load python/2.7
module load plink

cd ../snappy
#python SNAPPY_v0.1.py --infile ../y_bfiles/chrY_hg19_final
#python SNAPPY_v0.1.py --infile ../y_male_only_bfiles/chrY_male_hemizygous_only_het_filter_hg19_final
#python SNAPPY_v0.1.py --infile ../y_ukbb/chrY_male_only
#python SNAPPY_v0.1.py --infile ../y_nabec_files/nabec_males_only_hg19_chrY
python SNAPPY_v0.1.py --infile ../y_neurox/neurox_chrY_male_only