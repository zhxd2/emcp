export PATH=/pub/software/eggnog-mapper-1.0.3/bin:$PATH
EMCP_HOME=/home/zhxd/software/emcp

# step1, run eggnog-mapper
#python2 /pub/software/eggnog-mapper-1.0.3/emapper.py -m diamond  -i GDDH13_1-1_prot.fasta -o my --cpu 20 

# step2, make orgDB
Rscript $EMCP_HOME/makeOrgPackageFromEmapper.R my.emapper.annotations

# step3, statictics
Rscript $EMCP_HOME/AnnoStat.R GDDH13_1-1_prot.fasta

# enrichment
Rscript $EMCP_HOME/enrich.R --de_result genes.counts.matrix.BLO_S2_vs_KID_S2.DESeq2.DE_results \
	--de_log2FoldChange 1 \
	--de_padj 0.05
