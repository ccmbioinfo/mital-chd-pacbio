#!/bin/bash

source activate /hpf/largeprojects/ccmbio/mcouse/smital/envs

for f in `(find /hpf/largeprojects/smital/tcag2smital/pacbio/ -name '*hg38_sv.tsv*')`
do
        python3 annotate_filter_PacBio_pbsv.py -sv $f \
                -cmh /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread/annotate_SV/annotations/CMH.GRCh38.pbsv.bed \
                -chd /hpf/largeprojects/ccmbio/mcouse/smital/PacBio_CHD/CHD_tier_genes_with_class.tsv \
                -output /hpf/largeprojects/ccmbio/mcouse/smital/PacBio_CHD/annotated_output 
done
