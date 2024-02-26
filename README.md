# Add-seq: Chromatin accessibility profiling with small molecule intercalation and long-read sequencing
## Install Dependencies
```{bash}
conda install -f addseq.yml
```

## 0. Data Preprocessing
### 00. Basecalling and reads mapping
```{bash}
guppy_basecaller -i ~/data/fast5 -s ~/data/basecalled --bam_out -c ~/tools/guppy_basecaller_4.2.0/data/res_dna_r941_min_modbases-all-context_v001.cfg -x 'cuda:0' -r
```
```{bash}
minimap2 -x map-ont --d ref/sacCer3-ont.mmi ref/sacCer3.fa

minimap2 -t 32 -a --secondary=yes -x map-ont ref/sacCer3-ont.mmi basecalled/all_read.fastq | samtools view -S -b - | samtools sort - -o chrom.sorted.bam

samtools flagstat chrom.sorted.bam > chrom.sorted.flagstat

# Filter secondery and supplementary reads
samtools view -b -F SECONDARY,SUPPLEMENTARY chrom.sorted.bam | samtools sort -o chrom_pass.sorted.bam

samtools index chrom_pass.sorted.bam
```

### 01. Align signal events to reference
```{bash}
nanopolish index -d unbasecalled/ basecalled/all_read.fastq

nanopolish eventalign 
    --threads 16 --samples \
    --signal-index \
    --reads basecalled/all_read.fastq \
    --bam mapping/chrom_pass.sorted.bam \
    --genome ref/sacCer3.fa \
    --scale-events \
    --progress \
    --print-read-names > eventalign/chrom_pass.sorted.eventalign.txt
```
# Analyze Add-seq data with statistical analysis

All code for this is in scripts/addseqpaperfigscode.py

# Analyze Add-seq data with deep learning method

## Clone Add-seq neural network toolkit
```{bash}
git clone https://github.com/baigal628/NEMO.git
```

## 1. Train deep neural network model

```{bash}
python NEMO/nanopore_train_simple.py \
    --exp_id full \
    --device auto \
    --neg_data ctrl/sigalign/231128_unique.0_all_sig.tsv \
    --pos_data ctrl/sigalign/231128_unique.500_all_sig.tsv \
    --model_type resnet \
    --outpath train
```

## 2. Predict angelicin modification with deep neural network

## 20. Prepare signal alignment data
```{bash}
python3 NEMO/findNemo.py  \
    --mode init \
    --region all \
    --bam mapping/chrom_pass.sorted.bam \
    --genome ref/sacCer3.fa \
    --eventalign eventalign/chrom_pass.sorted.eventalign.txt \
    --outpath sigalign/ \
    --prefix chrom_pass.sorted

ll sigalign/
chrom_pass.sorted_all_readID.tsv
chrom_pass.sorted_all_sig.tsv
chrI_sig.tsv
chrII_sig.tsv
...
chrXVI_sig.tsv
```

## 21. Predict on specific loci
```{bash}
#Predict on CLN2 promoter region: chrXVI:66000-67550

python3 NEMO/findNemo.py  \
    --mode predict \
    --step 40 \
    --region chrXVI:66000-67550 \
    --bam mapping/chrom.sorted.bam \
    --genome ref/sacCer3.fa \
    --sigalign sigalign/chrXVI_sig.tsv \
    --outpath modPredict/ \
    --prefix resnet1d_step40_CLN2 \
    --readlist sigalign/chrom_pass.sorted_all_readID.tsv \
    --threads 32 \
    --load 50 \
    --weight train/best_models/resnet_0.001.pt

# Setting load=50 will save prediction file for every 50 reads
ll modPredict/
resnet1d_step40_CLN2_chrXVI:66000-67560_50_prediction.tsv
resnet1d_step40_CLN2_chrXVI:66000-67560_75_prediction.tsv

cat modPredict/resnet1d_step40_CLN2_chrXVI:66000-67560_*_prediction.tsv > modPredict/resnet1d_step40_CLN2_chrXVI:66000-67560_all_prediction.tsv

```

## 22. Predict on whole chromsome
```{bash}
python3 NEMO/findNemo.py  \
    --mode predict \
    --region chrV \
    --bam mapping/chrom_pass.sorted.bam \
    --genome ref/sacCer3.fa \
    --sigalign sigalign/chrV_sig.tsv \
    --outpath modPredict/ \
    --prefix resnet1d \
    --readlist sigalign/chrom_pass.sorted_all_readID.tsv \
    --threads 32 \
    --load 500 \
    --weight train/best_models/resnet_0.001.pt

ll modPredict/
resnet1d_chrV_500_prediction.tsv
resnet1d_chrX_1000_prediction.tsv
resnet1d_chrX_1500_prediction.tsv
resnet1d_chrX_2000_prediction.tsv
resnet1d_chrX_2100_prediction.tsv

cat modPredict/resnet1d_chrV_*_prediction > modPredict/resnet1d_chrV_all_prediction.tsv
```
