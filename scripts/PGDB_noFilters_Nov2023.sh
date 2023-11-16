#!/bin/bash

#SBATCH --job-name=${3}
#SBATCH --output=PGDB_${3}.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80Gb
#SBATCH --time=24:00:00
#SBATCH --partition=short

module load freebayes/1.3.4
module load gffcompare/0.11.5
module load stringtie/2.1.1
module load samtools/1.10


echo $1
echo $2
echo $3

   
cd /scratch/tsour.s/Sat_LSCC/RNAseq_data/$1

#sort gdc bam file by read name
if [[ ! -f $3.collate.bam ]]
then 
	srun samtools collate -o $3.collate.bam $2
	echo 'sorted bam'
fi 
 
#write paired end reads to 2 separate files
if [[ ! -f $3.paired_1.fq ]]
then 
	srun samtools fastq -1 $3.paired_1.fq -2 $3.paired_2.fq -n $3.collate.bam
	echo 'fasta conversion complete'
fi
    
# trim adapters, read quality filtering, make QC outputs [adapted from Spritz]
if [[ $3.paired.trim_1.fq ]]
then 
	srun /home/tsour.s/fastp -q 28 -i $3.paired_1.fq -I $3.paired_2.fq -o $3.paired.trim_1.fq -O $3.paired.trim_2.fq -h $3.fastp_QC_report.html -j $3.fastp_QC_report.json -w 16 -R $3.fastp_report --detect_adapter_for_pe
	echo "adapter trimming complete"
fi

#align using hisat2
if [[ ! -f $3.hisat2_summary.txt ]]
then 
	srun /home/tsour.s/hisat2-2.2.1/hisat2 -x /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic -1 $3.paired.trim_1.fq -2 $3.paired.trim_2.fq --dta -S $3.hisat2.sam --summary-file $3.hisat2_summary.txt -p 28
	echo 'hisat2 alignment complete'
fi
 
#convert to bam
if [[  -f $3.hisat2.sam ]]
then 
	srun samtools view -bS $3.hisat2.sam > $3.hisat2.bam
	srun samtools view -bS -h $3.hisat2.sam > $3.hisat2.header.bam #with header for customProDB
	echo 'bam conversion complete'
fi
   
#sort
if [[ ! -f $3.sorted.header.bam ]]
then 
	srun samtools sort $3.hisat2.bam -o $3.sorted.bam
	srun samtools sort $3.hisat2.header.bam -o $3.sorted.header.bam
	echo 'sorted bam'
fi
  
#create index
if [[ ! -f $3.sorted.header.bai ]]
then 
	srun samtools index -b $3.sorted.bam $3.sorted.bai
	srun samtools index -b $3.sorted.header.bam $3.sorted.header.bai
	echo 'index created'
fi

#assemble with string tie
if [[ ! -f $3.stringtie.gtf ]]
then 
	srun stringtie $3.sorted.bam -o $3.stringtie.gtf -G /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.100.gff3 -f 0.001 -c 1 -p 28
	echo 'stringtie assembly complete'
fi

#compare and annotate assembly with reference gff
if [[ ! -f $3.gffcmp ]]
then 
	srun gffcompare -r /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.100.gff3 -o $3.gffcmp $3.stringtie.gtf
	echo 'gffcompare complete'
fi
   
#filter for CDS (exons) only, and filter for all else and write these 2 transcript fastas to file
if [[ ! -f $3.gffcmp.annotated.transcript.gtf ]]
then 
	#srun /home/tsour.s/gffread/gffread $3.gffcmp.annotated.gtf --no-pseudo --force-exons -M -T -o $3.gffcmp.annotated.transcript.CDS.gtf
	#srun /home/tsour.s/gffread/gffread -w $3.CDS.DNA.fasta -g /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa $3.gffcmp.annotated.transcript.CDS.gtf
	srun /home/tsour.s/gffread/gffread $3.gffcmp.annotated.gtf --no-pseudo -M -T -o $3.gffcmp.annotated.transcript.gtf
	srun /home/tsour.s/gffread/gffread -w $3.all.DNA.fasta -g /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa $3.gffcmp.annotated.transcript.gtf
	echo 'genomic fastas written'
fi

#convert off compare output gtf to bed format, high confident and all but unknown/opposite strand/repeat
if [[ ! -f $3.gffcmp.annotated.all.bed  ]]
then 
	#srun python /home/tsour.s/scripts_templates/gffcompare_to_bed.py $3.gffcmp.annotated.gtf $3.gffcmp.annotated.HC.bed -C "=,c,k,m,n,j,e"
	srun python /home/tsour.s/scripts_templates/gffcompare_to_bed.py $3.gffcmp.annotated.gtf $3.gffcmp.annotated.all.bed
	echo 'bed conversion complete'
fi
   
#translate bed
if [[ ! -f $3.translation.fa ]]
then 
	#srun python /home/tsour.s/scripts_templates/translate_bed.py $3.gffcmp.annotated.HC.bed --twobit=/home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.2bit --cds --min_length 10 --reference GRCh38 --bed $3.CDS.translation.bed --fasta $3.CDS.translation.fa -v
	srun python /home/tsour.s/scripts_templates/translate_bed.py $3.gffcmp.annotated.all.bed --twobit=/home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.2bit --min_length 10 --reference GRCh38 --bed $3.translation.bed --fasta $3.translation.fa -v
	echo 'translation complete'
fi

#freeBayes variant calling
if [[ ! -f $3.vcf ]]
then 
	srun freebayes -b $3.sorted.bam -f /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa -C 1 -F 0.001  -v $3.vcf
	echo 'variant calling complete'
fi
  
#customProDB R script
if [[ ! -f $3.indel.fasta ]]
then 
	srun Rscript /home/tsour.s/scripts_templates/customProDB_ensembl_discovery.R $3.sorted.header.bam $3.vcf /home/tsour.s/genomes/customProDB_files/ensembl/ $3.
	echo 'customProDB complete'
   
fi

#merge fastas
if [[ ! -f $3.CDS.final.fasta ]]
then 
	srun python /home/tsour.s/scripts_templates/fasta_merge_files_and_filter_unique_sequences.py $3.All.final.fasta sequence '^>([^ ]+).*$' $3.indel.fasta $3.rpkm.fasta $3.SNV.fasta $3.translation.fa
	echo 'fasta merge complete'
fi


if [[ -f $3.All.final.fasta ]]
then
	mkdir files2transfer
	mv $3.vcf files2transfer/$3.vcf
	mv $3.translation.fa files2transfer/$3.translation.fa
	mv $3.translation.bed files2transfer/$3.translation.bed
	mv $3.stringtie.gtf files2transfer/$3.stringtie.gtf
	mv $3.rpkm.fasta files2transfer/$3.rpkm.fasta
	mv $3.indel.fasta files2transfer/$3.indel.fasta
	mv $3.hisat2_summary.txt files2transfer/$3.hisat2_summary.txt
	mv $3.gffcmp.tracking files2transfer/$3.gffcmp.tracking
	mv $3.gffcmp.loci files2transfer/$3.gffcmp.loci
	mv $3.gffcmp.annotated.transcript.gtf files2transfer/$3.gffcmp.annotated.transcript.gtf
	mv $3.gffcmp.annotated.gtf files2transfer/$3.gffcmp.annotated.gtf
	mv $3.gffcmp.annotated.all.bed files2transfer/$3.gffcmp.annotated.all.bed
	mv $3.gffcmp.$3.stringtie.gtf.tmap files2transfer/$3.gffcmp.$3.stringtie.gtf.tmap
	mv $3.gffcmp.$3.stringtie.gtf.refmap files2transfer/$3.gffcmp.$3.stringtie.gtf.refmap
	mv $3.gffcmp files2transfer/$3.gffcmp
	mv $3.fastp_QC_report.json files2transfer/$3.fastp_QC_report.json
	mv $3.fastp_QC_report.html files2transfer/$3.fastp_QC_report.html
	mv $3.all.DNA.fasta files2transfer/$3.all.DNA.fasta
	mv $3.SNV.fasta files2transfer/$3.SNV.fasta
	mv $3.SNV.tab files2transfer/$3.SNV.tab
	mv $3.All.final.fasta files2transfer/$3.All.final.fasta
	zip -r $3.files2transfer.zip files2transfer
fi





