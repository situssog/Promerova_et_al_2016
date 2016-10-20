# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

module load bioinfo-tools
module load bowtie2
module load samtools

head -2 ../4513_197_1.fas > reference_4513_197_1.fasta
sed 's/^>.*/>reference_seq/g' reference_4513_197_1.fasta
samtools faidx reference_4513_197_1.fasta

bowtie2-build reference_4513_197_1.fasta reference_4513_197_1

sed -e 's/^--//g' ../4513_197_1_R1.fastq | sed -e '/^$/d' > 4513_197_1_R1.fastq
sed -e 's/^--//g' ../4513_197_1_R2.fastq | sed -e '/^$/d' > 4513_197_1_R2.fastq

bowtie2 -x reference_4513_197_1 -1 4513_197_1_R1.fastq -2 4513_197_1_R2.fastq -S 4513_197_1.sam

rm -rf 4513_197_1_R1.fastq 4513_197_1_R2.fastq


samtools view -b -S 4513_197_1.sam > 4513_197_1.bam
samtools sort 4513_197_1.bam 4513_197_1_sorted
samtools index 4513_197_1_sorted.bam 4513_197_1_sorted.bai

python /home/sergio/Downloads/Platypus_0.7.8/Platypus.py callVariants --bamFiles=4513_197_1_sorted.bam --refFile=reference_4513_197_1.fasta --output=4513_197_1_VariantCalls.vcf


