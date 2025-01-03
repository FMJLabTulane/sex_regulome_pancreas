##Analysis file created by Dr. Ruth Elgamal PhD
##Munge T2D sex-stratified GWAS meta-analysis
python /nfs/lab/ldsc/munge_sumstats.py \
--sumstats ~/LDSC/gwas/Mahajan.NatGenet2018b.T2D.{$SEX}.European.hg19.rsid.txt \
--a1 NEA \
--a2 EA \
--p Pvalue \
--snp rsid \
--N-col Neff \
--chunksize 500000 \
--frq EAF \
--maf-min 0.01 \
--signed-sumstats Beta,0 \
--merge-allele ~/LDSC/w_hm3.snplist \
--out ~/LDSC/gwas/Mahajan.NatGenet2018b.T2D.{$SEX}.European.hg19.rsid.ldsc

##Make functional annotations from cell type-specific snATACseq peaks from males and females separately
for annot in $(cat ~/LDSC/celltypes_w_sex.txt); do
mkdir ~/LDSC/celltype_sex_peaks/${annot}_hg19
for i in {1..22}; do
python ~/LDSC/make_annot.py \
--bed-file ~/Fahd/LDSC/celltype_sex_peaks/${annot}_peaks_hg19.bed \
--bimfile ~/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim \
--annot-file ~/LDSC/celltype_sex_peaks/${annot}_hg19/${annot}_hg19.${i}.annot.gz
done
done

##Run LD score regression on functional annotations
for annot in $(cat ~/LDSC/celltypes_w_sex.txt); do
for i in {1..22}; do
python ~/LDSC/ldsc.py \
--print-snps ~/LDSC/1000G_EUR_Phase3_baseline_snps/hm.${i}.snp \
--ld-wind-cm 1.0 --out ~/LDSC/celltype_sex_peaks/${annot}_hg19/${annot}_hg19.${i} \
--bfile ~/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} \
--thin-annot \
--annot ~/LDSC/celltype_sex_peaks/${annot}_hg19/${annot}_hg19.${i}.annot.gz \
--l2
done
done

##Run partitioned heritability
for sex in Mahajan.NatGenet2018b.T2D.FEMALE.European Mahajan.NatGenet2018b.T2D.MALE.European; do
for annot in $(cat ~/Fahd/LDSC/celltypes_w_sex.txt); do
python ~/LDSC/ldsc.py \
--h2 ~/LDSC/gwas/${sex}.hg19.rsid.ldsc.sumstats.gz \
--ref-ld-chr ~/LDSC/celltype_sex_peaks/${annot}_hg19/${annot}_hg19.,~/LDSC/1000G_EUR_Phase3_baseline/baseline. \
--out ~/LDSC/celltype_sex_results/${sex}_hg19.${annot} \
--overlap-annot  \
--frqfile-chr ~/LDSC/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ~/LDSC/weights_hm3_no_hla/weights. \
--print-coefficients
done
done
