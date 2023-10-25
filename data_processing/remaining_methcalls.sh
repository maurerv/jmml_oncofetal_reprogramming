#/bin/bash

BAMDIR="/omics/odcf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing"
OUTDIR="/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md"
BAMFILES=($(ls $BAMDIR/view-by-pid/*/*/paired/merged-alignment/*.bam | grep "FL\|FS\|JU"))

NPROC=15
mkdir -p ${OUTDIR}/Mbias
mkdir -p ${OUTDIR}//methylationCalls/
for file in "${BAMFILES[@]}"; do

    fname=$(basename $file)
    dirname=$(echo $fname | awk -F\. '{print $1}')
    bsub <<-EOT
        #!/bin/bash
        #BSUB -o $HOME/logs
        #BSUB $HOME/err
        #BSUB -R rusage[mem=20000]
        #BSUB -n $NPROC
        #BSUB -q long

        source $HOME/.bashrc
        conda activate jmml-thesis

        # Check for mbias
        # MethylDackel mbias -@ $NPROC \
        #     /omics/groups/OE0219/internal/genomes/Hsapiens/GRCh37/seq/hs37d5_PhiX_Lambda.fa \
        #     ${file} \
        #     ${OUTDIR}/Mbias/${fname}

        # Extract mbias adjusted methylation values
        mkdir -p ${OUTDIR}/methylationCalls/${dirname}/
        MethylDackel extract -@ 5 --OT 6,146,2,144 --OB 7,146,12,150  \
            /omics/groups/OE0219/internal/genomes/Hsapiens/GRCh37/seq/hs37d5_PhiX_Lambda.fa\
            ${file} \
            -o ${OUTDIR}/methylationCalls/${dirname}/${fname}

        # Replace chromosome names in bedGraph
        sed 's/^\([0-9].*\)/chr\1/g' ${OUTDIR}/methylationCalls/${dirname}/${fname}_CpG.bedGraph > \
            ${OUTDIR}/methylationCalls/${dirname}/${fname}.chr.bedGraph
EOT
done


#Check for mbias
# mkdir -p ${OUTDIR}/Mbias
# for file in "${BAMFILES[@]}"; do
#     fname=$(basename $file)
#     echo $fname
#     MethylDackel mbias -@ 10 \
#         /omics/groups/OE0219/internal/genomes/Hsapiens/GRCh37/seq/hs37d5_PhiX_Lambda.fa \
#         $file \
#         ${OUTDIR}/Mbias/${fname}
# done

#Extract methylation values with adjusting for mbias
mkdir -p ${OUTDIR}//methylationCalls/
for file in "${BAMFILES[@]}"; do
    fname=$(basename $file)
    dirname=$(echo $fname | awk -F\. '{print $1}')
    echo $dirname
    echo $fname
    mkdir -p ${OUTDIR}/methylationCalls/${dirname}/
    MethylDackel extract -@ 5 --OT 6,146,2,144 --OB 7,146,12,150  \
        /omics/groups/OE0219/internal/genomes/Hsapiens/GRCh37/seq/hs37d5_PhiX_Lambda.fa\
        $file \
        -o ${OUTDIR}/methylationCalls/${dirname}/${fname}
done

replace chromosome names
for file in ${OUTDIR}/methylationCalls/*/*.bedGraph; do
    awk '{print "chr" $0}' $file > $file.chr.bedgraph
    echo $file
done