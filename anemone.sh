#wget -r --no-parent ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP060/SRP060291
# the wget is not needed since all files can be downloaded directly from the fastq-dump command, and no FTP is encouraged by NCBI anymore. 



#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR181/SRR1819888/SRR1819888.sra
#fastq-dump --split-files SRR1952742
#fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1819888.sra

# concantenate them before or after the fastqc process ???? 

#fastqc SRR1952742_1.fastq --outdir=/mnt/dnarules/Anuj/
#fastqc SRR1952742_2.fastq --outdir=/mnt/dnarules/Anuj/

fastqc CGF1_R1_paired_trimmed.fastq

java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE CGF1_R1.fastq.gz CGF1_R2.fastq.gz CGF1_R1_paired_trimmed.fastq CGF1_R1_unpaired_trimmed.fastq CGF1_R2_paired_trimmed.fastq CGF1_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE CGF2_R1.fastq.gz CGF2_R2.fastq.gz CGF2_R1_paired_trimmed.fastq CGF2_R1_unpaired_trimmed.fastq CGF2_R2_paired_trimmed.fastq CGF2_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE CGF3_R1.fastq.gz CGF3_R2.fastq.gz CGF3_R1_paired_trimmed.fastq CGF3_R1_unpaired_trimmed.fastq CGF3_R2_paired_trimmed.fastq CGF3_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE CGT1_R1.fastq.gz CGT1_R2.fastq.gz CGT1_R1_paired_trimmed.fastq CGT1_R1_unpaired_trimmed.fastq CGT1_R2_paired_trimmed.fastq CGT1_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE CGT2_R1.fastq.gz CGT2_R2.fastq.gz CGT2_R1_paired_trimmed.fastq CGT2_R1_unpaired_trimmed.fastq CGT2_R2_paired_trimmed.fastq CGT2_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE CGT3_R1.fastq.gz CGT3_R2.fastq.gz CGT3_R1_paired_trimmed.fastq CGT3_R1_unpaired_trimmed.fastq CGT3_R2_paired_trimmed.fastq CGT3_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE EQF1_R1.fastq.gz EQF1_R2.fastq.gz EQF1_R1_paired_trimmed.fastq EQF1_R1_unpaired_trimmed.fastq EQF1_R2_paired_trimmed.fastq EQF1_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE EQF2_R1.fastq.gz EQF2_R2.fastq.gz EQF2_R1_paired_trimmed.fastq EQF2_R1_unpaired_trimmed.fastq EQF2_R2_paired_trimmed.fastq EQF2_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE EQF3_R1.fastq.gz EQF3_R2.fastq.gz EQF3_R1_paired_trimmed.fastq EQF3_R1_unpaired_trimmed.fastq EQF3_R2_paired_trimmed.fastq EQF3_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE EQT1_R1.fastq.gz EQT1_R2.fastq.gz EQT1_R1_paired_trimmed.fastq EQT1_R1_unpaired_trimmed.fastq EQT1_R2_paired_trimmed.fastq EQT1_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE EQT2_R1.fastq.gz EQT2_R2.fastq.gz EQT2_R1_paired_trimmed.fastq EQT2_R1_unpaired_trimmed.fastq EQT2_R2_paired_trimmed.fastq EQT2_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20
java -jar /home/dnarules/Downloads/Trimmomatic-0.33/trimmomatic-0.33.jar PE EQT3_R1.fastq.gz EQT3_R2.fastq.gz EQT3_R1_paired_trimmed.fastq EQT3_R1_unpaired_trimmed.fastq EQT3_R2_paired_trimmed.fastq EQT3_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/dnarules/Downloads/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20

fastqc *_paired_trimmed.fastq 

/home/dnarules/Trinity/Trinity --seqType fq --left CGF1_R1_paired_trimmed.fastq --right CGF1_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./CGF1/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left CGF2_R1_paired_trimmed.fastq --right CGF2_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./CGF2/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left CGF3_R1_paired_trimmed.fastq --right CGF3_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./CGF3/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left CGT1_R1_paired_trimmed.fastq --right CGT1_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./CGT1/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left CGT2_R1_paired_trimmed.fastq --right CGT2_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./CGT2/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left CGT3_R1_paired_trimmed.fastq --right CGT3_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./CGT3/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left EQF1_R1_paired_trimmed.fastq --right EQF1_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./EQF1/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left EQF2_R1_paired_trimmed.fastq --right EQF2_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./EQF2/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left EQF3_R1_paired_trimmed.fastq --right EQF3_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./EQF3/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left EQT1_R1_paired_trimmed.fastq --right EQT1_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./EQT1/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left EQT2_R1_paired_trimmed.fastq --right EQT2_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./EQT2/trinity"
/home/dnarules/Trinity/Trinity --seqType fq --left EQT3_R1_paired_trimmed.fastq --right EQT3_R2_paired_trimmed.fastq --CPU 18 --max_memory 50G --output "./EQT3/trinity"

/home/dnarules/Trinity/Trinity --seqType fq --max_memory 50G  \
         --left CGF1_R1_paired_trimmed.fastq,CGF2_R1_paired_trimmed.fastq,CGF3_R1_paired_trimmed.fastq,CGT1_R1_paired_trimmed.fastq,CGT2_R1_paired_trimmed.fastq,CGT3_R1_paired_trimmed.fastq \
         --right CGF1_R2_paired_trimmed.fastq,CGF2_R2_paired_trimmed.fastq,CGF3_R2_paired_trimmed.fastq,CGT1_R2_paired_trimmed.fastq,CGT2_R2_paired_trimmed.fastq,CGT3_R2_paired_trimmed.fastq \
         --CPU 32  
	 --output "./CG_trinity"


/home/dnarules/Trinity/Trinity --seqType fq --max_memory 50G  \
         --left EQF1_R1_paired_trimmed.fastq,EQF2_R1_paired_trimmed.fastq,EQF3_R1_paired_trimmed.fastq,EQT1_R1_paired_trimmed.fastq,EQT2_R1_paired_trimmed.fastq,EQT3_R1_paired_trimmed.fastq \
         --right EQF1_R2_paired_trimmed.fastq,EQF2_R2_paired_trimmed.fastq,EQF3_R2_paired_trimmed.fastq,EQT1_R2_paired_trimmed.fastq,EQT2_R2_paired_trimmed.fastq,EQT3_R2_paired_trimmed.fastq \
         --CPU 32  
	 --output "./EQ_trinity"

cd /home/dnarules/Anuj/ProjectAnemone/raw

/home/dnarules/Trinity/Trinity --seqType fq --max_memory 10G  \
         --left CGF1_R1_paired_trimmed.fastq,CGF2_R1_paired_trimmed.fastq,CGF3_R1_paired_trimmed.fastq \
         --right CGF1_R2_paired_trimmed.fastq,CGF2_R2_paired_trimmed.fastq,CGF3_R2_paired_trimmed.fastq \
         --CPU 8  \
	 --output "./CGF_trinity"


/home/dnarules/Trinity/Trinity --seqType fq --max_memory 10G  \
         --left CGT1_R1_paired_trimmed.fastq,CGT2_R1_paired_trimmed.fastq,CGT3_R1_paired_trimmed.fastq \
         --right CGT1_R2_paired_trimmed.fastq,CGT2_R2_paired_trimmed.fastq,CGT3_R2_paired_trimmed.fastq \
         --CPU 8  \
	 --output "./CGT_trinity"


/home/dnarules/Trinity/Trinity --seqType fq --max_memory 10G  \
         --left EQF1_R1_paired_trimmed.fastq,EQF2_R1_paired_trimmed.fastq,EQF3_R1_paired_trimmed.fastq \
         --right EQF1_R2_paired_trimmed.fastq,EQF2_R2_paired_trimmed.fastq,EQF3_R2_paired_trimmed.fastq \
         --CPU 8  \
	 --output "./EQF_trinity"

/home/dnarules/Trinity/Trinity --seqType fq --max_memory 10G  \
         --left EQT1_R1_paired_trimmed.fastq,EQT2_R1_paired_trimmed.fastq,EQT3_R1_paired_trimmed.fastq \
         --right EQT1_R2_paired_trimmed.fastq,EQT2_R2_paired_trimmed.fastq,EQT3_R2_paired_trimmed.fastq \
         --CPU 8  \
	 --output "./EQT_trinity"

/home/dnarules/Trinity/util/TrinityStats.pl /home/dnarules/Anuj/ProjectAnemone/trinity/CG_trinity/CG_Trinity.fasta
/home/dnarules/Trinity/util/TrinityStats.pl /home/dnarules/Anuj/ProjectAnemone/trinity/EQ_trinity/EQ_Trinity.fasta
/home/dnarules/Trinity/util/TrinityStats.pl /home/dnarules/Anuj/ProjectAnemone/trinity/CGT_assembled_RNA.fasta
/home/dnarules/Trinity/util/TrinityStats.pl /home/dnarules/Anuj/ProjectAnemone/trinity/CGF_assembled_RNA.fasta
/home/dnarules/Trinity/util/TrinityStats.pl /home/dnarules/Anuj/ProjectAnemone/trinity/EQT_assembled_RNA.fasta
/home/dnarules/Trinity/util/TrinityStats.pl /home/dnarules/Anuj/ProjectAnemone/trinity/EQF_assembled_RNA.fasta

#### need to do some ORF or CDS prediction here ... we do not need to do annotation. 
/home/dnarules/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /home/dnarules/Anuj/ProjectAnemone/trinity/CG_trinity/CG_Trinity.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.Predict -t /home/dnarules/Anuj/ProjectAnemone/trinity/CG_trinity/CG_Trinity.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /home/dnarules/Anuj/ProjectAnemone/trinity/EQ_trinity/EQ_Trinity.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.Predict -t /home/dnarules/Anuj/ProjectAnemone/trinity/EQ_trinity/EQ_Trinity.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /home/dnarules/Anuj/ProjectAnemone/trinity/CGT_assembled_RNA.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.Predict -t /home/dnarules/Anuj/ProjectAnemone/trinity/CGT_assembled_RNA.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /home/dnarules/Anuj/ProjectAnemone/trinity/CGF_assembled_RNA.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.Predict -t /home/dnarules/Anuj/ProjectAnemone/trinity/CGF_assembled_RNA.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /home/dnarules/Anuj/ProjectAnemone/trinity/EQT_assembled_RNA.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.Predict -t /home/dnarules/Anuj/ProjectAnemone/trinity/EQT_assembled_RNA.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /home/dnarules/Anuj/ProjectAnemone/trinity/EQF_assembled_RNA.fasta
/home/dnarules/TransDecoder-3.0.1/TransDecoder.Predict -t /home/dnarules/Anuj/ProjectAnemone/trinity/EQF_assembled_RNA.fasta

#stats from transcdecoder
grep ">" CGF_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "complete" CGF_assembled_RNA.fasta.transdecoder.pep | wc -l 
grep "internal" CGF_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "5prime" CGF_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "3prime" CGF_assembled_RNA.fasta.transdecoder.pep | wc -l 
grep ">" CGT_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "complete" CGT_assembled_RNA.fasta.transdecoder.pep | wc -l 
grep "internal" CGT_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "5prime" CGT_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "3prime" CGT_assembled_RNA.fasta.transdecoder.pep | wc -l 
grep ">" EQF_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "complete" EQF_assembled_RNA.fasta.transdecoder.pep | wc -l 
grep "internal" EQF_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "5prime" EQF_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "3prime" EQF_assembled_RNA.fasta.transdecoder.pep | wc -l 
grep ">" EQT_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "complete" EQT_assembled_RNA.fasta.transdecoder.pep | wc -l 
grep "internal" EQT_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "5prime" EQT_assembled_RNA.fasta.transdecoder.pep | wc -l
grep "3prime" EQT_assembled_RNA.fasta.transdecoder.pep | wc -l 

cd /home/Anuj/ProjectAnemone/trinotate/CG
/home/dnarules/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  CG_Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
blastx -query CG_Trinity.fasta -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
blastp -query CG_transdecoder.pep -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm CG_transdecoder.pep > pfam.log
signalp -f short -n signalp.out CG_transdecoder.pep
tmhmm --short < CG_transdecoder.pep > tmhmm.out
/home/dnarules/Trinity/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dnarules/Anuj/ProjectAnemone/trinity/CG_Trinity.fasta >  CG_Trinity.fasta.gene_trans_map
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite init --gene_trans_map CG_Trinity.fasta.gene_trans_map --transcript_fasta CG_Trinity.fasta --transdecoder_pep CG_transdecoder.pep
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_tmhmm tmhmm.out
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_signalp signalp.out
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > CG_trinotate_annotation_report.xls
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > CG_trinotate_annotation_report.tsv
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_hydra_blastx hydra.uniprot.blastx.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_nematostella_blastx nematostella.uniprot.blastx.outfmt6



cd /home/Anuj/ProjectAnemone/trinotate/CGF
/home/dnarules/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  CGF_Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
blastx -query CGF_Trinity.fasta -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
blastp -query CGF_transdecoder.pep -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm CGF_transdecoder.pep > pfam.log
signalp -f short -n signalp.out CGF_transdecoder.pep
tmhmm --short < CGF_transdecoder.pep > tmhmm.out
/home/dnarules/Trinity/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dnarules/Anuj/ProjectAnemone/trinity/CGF_Trinity.fasta >  CGF_Trinity.fasta.gene_trans_map
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite init --gene_trans_map CGF_Trinity.fasta.gene_trans_map --transcript_fasta CGF_Trinity.fasta --transdecoder_pep CGF_transdecoder.pep
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_tmhmm tmhmm.out
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_signalp signalp.out
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > CGF_trinotate_annotation_report.xls
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > CGF_trinotate_annotation_report.tsv
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6


cd /home/Anuj/ProjectAnemone/trinotate/CGT
/home/dnarules/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  CGT_Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
blastx -query CGT_Trinity.fasta -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
blastp -query CGT_transdecoder.pep -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm CGT_transdecoder.pep > pfam.log
signalp -f short -n signalp.out CGT_transdecoder.pep
tmhmm --short < CGT_transdecoder.pep > tmhmm.out
/home/dnarules/Trinity/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dnarules/Anuj/ProjectAnemone/trinity/CGT_Trinity.fasta >  CGT_Trinity.fasta.gene_trans_map
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite init --gene_trans_map CGT_Trinity.fasta.gene_trans_map --transcript_fasta CGT_Trinity.fasta --transdecoder_pep CGT_transdecoder.pep
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite LOAD_tmhmm tmhmm.out
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite LOAD_signalp signalp.out
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > CGT_trinotate_annotation_report.xls
/home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > CGT_trinotate_annotation_report.tsv

cd /home/Anuj/ProjectAnemone/trinotate/EQ
cd /home/dnarules/Anuj/ProjectAnemone/trinotate/EQ/
/home/dnarules/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  EQ_Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
blastx -query EQ_Trinity.fasta -db uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 > EQ_blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
blastp -query EQ_transdecoder.pep -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > EQ_blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm EQ_transdecoder.pep > pfam.log
/home/dnarules/signalp-4.1/signalp -f short -n signalp.out EQ_transdecoder.pep
/home/dnarules/tmhmm-2.0c/bin/tmhmm --short < EQ_transdecoder.pep > tmhmm.out
/home/dnarules/Trinity/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dnarules/Anuj/ProjectAnemone/trinity/EQ_Trinity.fasta >  EQ_Trinity.fasta.gene_trans_map
/home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite init --gene_trans_map EQ_Trinity.fasta.gene_trans_map --transcript_fasta EQ_Trinity.fasta --transdecoder_pep EQ_transdecoder.pep
/home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite LOAD_swissprot_blastp EQ_blastp.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite LOAD_swissprot_blastx EQ_blastx.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite LOAD_tmhmm tmhmm.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite LOAD_signalp signalp.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > EQ_trinotate_annotation_report.xls


cd /home/Anuj/ProjectAnemone/trinotate/EQF
/home/dnarules/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  EQF_Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
blastx -query EQF_Trinity.fasta -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastx.outfmt6
#blastx -query CG_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
#blastx -query CG_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
blastp -query EQF_transdecoder.pep -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm EQF_transdecoder.pep > pfam.log
signalp -f short -n signalp.out EQF_transdecoder.pep
tmhmm --short < EQF_transdecoder.pep > tmhmm.out
/home/dnarules/Trinity/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dnarules/Anuj/ProjectAnemone/trinity/EQF_Trinity.fasta >  EQF_Trinity.fasta.gene_trans_map
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite init --gene_trans_map EQF_Trinity.fasta.gene_trans_map --transcript_fasta EQF_Trinity.fasta --transdecoder_pep EQF_transdecoder.pep
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite LOAD_tmhmm tmhmm.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite LOAD_signalp signalp.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > EQF_trinotate_annotation_report.xls
/home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > EQF_trinotate_annotation_report.tsv

cd /home/Anuj/ProjectAnemone/trinotate/EQT
/home/dnarules/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  EQT_Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
blastx -query EQT_Trinity.fasta -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
blastp -query EQT_transdecoder.pep -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > swissprot.blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm CGT_transdecoder.pep > pfam.log
signalp -f short -n signalp.out CGT_transdecoder.pep
tmhmm --short < CGT_transdecoder.pep > tmhmm.out
/home/dnarules/Trinity/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/dnarules/Anuj/ProjectAnemone/trinity/EQT_Trinity.fasta >  EQT_Trinity.fasta.gene_trans_map
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite init --gene_trans_map EQT_Trinity.fasta.gene_trans_map --transcript_fasta EQT_Trinity.fasta --transdecoder_pep EQT_transdecoder.pep
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite LOAD_tmhmm tmhmm.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite LOAD_signalp signalp.out
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > EQT_trinotate_annotation_report.xls
/home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC > EQT_trinotate_annotation_report.tsv


cd /home/dnarules/Anuj/ProjectAnemone/trinotate/CG
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
blastx -query CG_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CG_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_hydra_blastx hydra.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate CG_Trinotate.sqlite LOAD_nematostella_blastx nematostella.uniprot.blastx.outfmt6


cd /home/dnarules/Anuj/ProjectAnemone/trinotate/CGF
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
blastx -query CGF_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CGF_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_hydra_blastx hydra.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate CGF_Trinotate.sqlite LOAD_nematostella_blastx nematostella.uniprot.blastx.outfmt6


cd /home/dnarules/Anuj/ProjectAnemone/trinotate/CGT
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
blastx -query CGT_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query CGT_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite LOAD_hydra_blastx hydra.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate CGT_Trinotate.sqlite LOAD_nematostella_blastx nematostella.uniprot.blastx.outfmt6


cd /home/dnarules/Anuj/ProjectAnemone/trinotate/EQ
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
blastx -query EQ_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query EQ_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite LOAD_hydra_blastx hydra.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate EQ_Trinotate.sqlite LOAD_nematostella_blastx nematostella.uniprot.blastx.outfmt6


cd /home/dnarules/Anuj/ProjectAnemone/trinotate/EQF
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
blastx -query EQF_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query EQF_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite LOAD_hydra_blastx hydra.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate EQF_Trinotate.sqlite LOAD_nematostella_blastx nematostella.uniprot.blastx.outfmt6


cd /home/dnarules/Anuj/ProjectAnemone/trinotate/EQT
makeblastdb -in uniprot_hydra.fasta -dbtype prot
makeblastdb -in uniprot_nematostella.fasta -dbtype prot
blastx -query EQT_Trinity.fasta -db uniprot_hydra.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > hydra.uniprot.blastx.outfmt6
blastx -query EQT_Trinity.fasta -db uniprot_nematostella.fasta -num_threads 24 -max_target_seqs 1 -outfmt 6 > nematostella.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite LOAD_hydra_blastx hydra.uniprot.blastx.outfmt6
home/dnarules/Trinotate-3.0.2/Trinotate EQT_Trinotate.sqlite LOAD_nematostella_blastx nematostella.uniprot.blastx.outfmt6


###################################CG##############################
#make sure the hmm and blast databases are there

#In this case use cut and not grep. 
#To extract the blastx report: 


To extract the blastp report:
cut -f 7 CG_trinotate_annotation_report.tsv > CG_swissprot_blastp.tsv

To extract the blastp report and To remove lines with . only
cut -f 7 CG_trinotate_annotation_report.tsv |  grep -v '^.$' > CG_swissprot_blastp.tsv

To count number of metazoans: 
grep "Metazoa" CG_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" CG_swissprot_blastp.tsv |  cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" CG_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of ion channels:
grep "ion" CG_swissprot_blastp.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" CG_swissprot_blastp.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l


To extract the blastx report:
cut -f 3 CG_trinotate_annotation_report.tsv > CG_swissprot_blastx.tsv

To extract the blastx report and To remove lines with . only
cut -f 3 CG_trinotate_annotation_report.tsv |  grep -v '^.$' > CG_swissprot_blastx.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' CG_swissprot_blastx.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" CG_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" CG_swissprot_blastx.tsv |  cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" CG_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of ion channels:
grep "ion" CG_swissprot_blastx.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" CG_swissprot_blastx.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l



########################################EQ #################
#In this case use cut and not grep. 

To extract the blastp report:
cut -f 7 EQ_trinotate_annotation_report.tsv > EQ_swissprot_blastp.tsv

To extract the blastp report and To remove lines with . only
cut -f 7 EQ_trinotate_annotation_report.tsv |  grep -v '^.$' > EQ_swissprot_blastp.tsv


To count number of metazoans: 
grep "Metazoa" EQ_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" EQ_swissprot_blastp.tsv |  cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" EQ_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of ion channels:
grep "ion" EQ_swissprot_blastp.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" EQ_swissprot_blastp.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l



To extract the blastx report: 
cut -f 3 EQ_trinotate_annotation_report.tsv > EQ_swissprot_blastx.tsv

To extract the blastx report and To remove lines with . only
cut -f 3 EQ_trinotate_annotation_report.tsv |  grep -v '^.$' > EQ_swissprot_blastx.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' EQ_swissprot_blastx.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" EQ_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" EQ_swissprot_blastx.tsv |  cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" EQ_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of ion channels:
grep "ion" EQ_swissprot_blastx.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" EQ_swissprot_blastx.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l



################FOR EQT##################33
#To count the number of 
#In this case use cut and not grep. 

To extract the blastp report:
cut -f 7 EQT_trinotate_annotation_report.tsv > EQT_swissprot_blastp.tsv

To extract the blastp report and To remove lines with . only
cut -f 7 EQT_trinotate_annotation_report.tsv |  grep -v '^.$' > EQT_swissprot_blastp.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' EQT_swissprot_blastp.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" EQT_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" EQT_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" EQT_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" EQT_swissprot_blastp.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" EQT_swissprot_blastp.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l



To extract the blastx report: 
cut -f 3 EQT_trinotate_annotation_report.tsv > EQT_swissprot_blastx.tsv

To extract the blastx report and To remove lines with . only
cut -f 3 EQT_trinotate_annotation_report.tsv |  grep -v '^.$' > EQT_swissprot_blastx.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' EQT_swissprot_blastx.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" EQT_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" EQT_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" EQT_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" EQT_swissprot_blastx.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" EQT_swissprot_blastx.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l


################FOR EQC##################33
#To count the number of 
#In this case use cut and not grep. 

To extract the blastp report:
cut -f 7 EQF_trinotate_annotation_report.tsv > EQF_swissprot_blastp.tsv

To extract the blastp report and To remove lines with . only
cut -f 7 EQF_trinotate_annotation_report.tsv |  grep -v '^.$' > EQF_swissprot_blastp.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' EQF_swissprot_blastp.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" EQF_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" EQF_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" EQF_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" EQF_swissprot_blastp.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" EQF_swissprot_blastp.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l



To extract the blastx report: 
cut -f 3 EQF_trinotate_annotation_report.tsv > EQF_swissprot_blastx.tsv

To extract the blastx report and To remove lines with . only
cut -f 3 EQF_trinotate_annotation_report.tsv |  grep -v '^.$' > EQF_swissprot_blastx.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' EQF_swissprot_blastx.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" EQF_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" EQF_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" EQF_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" EQF_swissprot_blastx.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" EQF_swissprot_blastx.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l

#####################################CGT##########################3
#In this case use cut and not grep. 
#To extract the blastx report: 


To extract the blastp report:
cut -f 7 CGT_trinotate_annotation_report.tsv > CGT_swissprot_blastp.tsv

To extract the blastp report and To remove lines with . only
cut -f 7 CGT_trinotate_annotation_report.tsv |  grep -v '^.$' > CGT_swissprot_blastp.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' CGT_swissprot_blastp.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" CGT_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" CGT_swissprot_blastp.tsv |  cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" CGT_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" CGT_swissprot_blastp.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" CGT_swissprot_blastp.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l


To extract the blastx report:
cut -f 3 CGT_trinotate_annotation_report.tsv > CGT_swissprot_blastx.tsv

To extract the blastx report and To remove lines with . only
cut -f 3 CGT_trinotate_annotation_report.tsv |  grep -v '^.$' > CGT_swissprot_blastx.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' CGT_swissprot_blastx.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" CGT_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" CGT_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" CGT_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" CGT_swissprot_blastx.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" CGT_swissprot_blastx.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l


#################################CGC#################333
#In this case use cut and not grep. 
#To extract the blastx report: 


To extract the blastp report:
cut -f 7 CGF_trinotate_annotation_report.tsv > CGF_swissprot_blastp.tsv

To extract the blastp report and To remove lines with . only
cut -f 7 CGF_trinotate_annotation_report.tsv |  grep -v '^.$' > CGF_swissprot_blastp.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' CGF_swissprot_blastp.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" CGF_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" CGF_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" CGF_swissprot_blastp.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" CGF_swissprot_blastp.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" CGF_swissprot_blastp.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l


To extract the blastx report:
cut -f 3 CGF_trinotate_annotation_report.tsv > CGF_swissprot_blastx.tsv

To extract the blastx report and To remove lines with . only
cut -f 3 CGF_trinotate_annotation_report.tsv |  grep -v '^.$' > CGF_swissprot_blastx.tsv

To count the number of unique All proteins: 
cut -f 1 -d '^' CGF_swissprot_blastx.tsv  | sort | uniq | wc -l 

To count number of metazoans: 
grep "Metazoa" CGF_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of cnidarians:
grep "Cnidaria" CGF_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "toxin" CGF_swissprot_blastx.tsv | cut -f 1 -d '^'  | sort | uniq | wc -l

To count the number of toxins:
grep "ion" CGF_swissprot_blastx.tsv | grep "channel" | cut -f 1 -d '^'  | sort | uniq | wc -l

Number of cnidarian toxins:
grep "Cnidaria" CGF_swissprot_blastx.tsv | grep "toxin" | cut -f 1 -d '^'  | sort | uniq | wc -l



#align the reads to the assembly and count the number of reads that are mapped etc. 


PATH=$PATH:/home/dnarules/express-1.5.1-linux_x86_64 ; export PATH
source ~/.bashrc

#build bowtie index
bowtie2-build EQ_Trinity.fasta EQ_Trinity.fasta
bowtie2-build CG_Trinity.fasta CG_Trinity.fasta
# get bowtie stats
bowtie2 -p 10 -q -x EQ_Trinity.fasta -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQF_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQF_R2_paired_trimmed.fastq 2>&1 1> /dev/null | tee align_stats.txt
bowtie2 -p 10 -q -x EQ_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQT_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQT_R2_paired_trimmed.fastq 2>&1 1> /dev/null | tee align_stats.txt
bowtie2 -p 10 -q -x CG_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/CGF_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/CGF_R2_paired_trimmed.fastq 2>&1 1> /dev/null | tee align_stats.txt
bowtie2 -p 10 -q -x CG_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/CGT_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/CGT_R2_paired_trimmed.fastq 2>&1 1> /dev/null | tee align_stats.txt
# run the bowtie mapping
bowtie2 -p 12 --time -x EQ_Trinity.fasta -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQF_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQF_R2_paired_trimmed.fastq -S /home/dnarules/Anuj/ProjectAnemone/bowtie/EQF/EQF.sam


bowtie2 -p 12 --time -x EQ_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQF1_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQF1_R2_paired_trimmed.fastq -S /home/dnarules/Anuj/ProjectAnemone/bowtie/EQF1.sam
bowtie2 -p 12 --time -x EQ_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQF2_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQF2_R2_paired_trimmed.fastq -S /home/dnarules/Anuj/ProjectAnemone/bowtie/EQF2.sam
bowtie2 -p 12 --time -x EQ_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQF3_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQF3_R2_paired_trimmed.fastq -S /home/dnarules/Anuj/ProjectAnemone/bowtie/EQF3.sam
bowtie2 -p 12 --time -x EQ_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQT1_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQT1_R2_paired_trimmed.fastq -S /home/dnarules/Anuj/ProjectAnemone/bowtie/EQT1.sam
bowtie2 -p 12 --time -x EQ_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQT2_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQT2_R2_paired_trimmed.fastq -S /home/dnarules/Anuj/ProjectAnemone/bowtie/EQT2.sam
bowtie2 -p 12 --time -x EQ_Trinity.fasta.bowtie2 -1 /home/dnarules/Anuj/ProjectAnemone/raw/EQT3_R1_paired_trimmed.fastq -2 /home/dnarules/Anuj/ProjectAnemone/raw/EQT3_R2_paired_trimmed.fastq -S /home/dnarules/Anuj/ProjectAnemone/bowtie/EQT3.sam


# change output from sam to bam
samtools view -b -S -o bowtie2.bam bowtie2.sam 
# do the read counting
express --rf-stranded -o /home/dnarules/Anuj/ProjectAnemone/express/EQF /home/dnarules/Anuj/ProjectAnemone/express/EQF/EQ_Trinity.fasta /home/dnarules/Anuj/ProjectAnemone/EQ/bowtie/EQF.bam

#extract the totcounts for each
cat ./EQT/EQT_express_outdir/EQT_results.xprs | cut -f 2,5 > EQT_transcripts_totCounts.tsv
cat ./EQF/EQF_express_outdir/EQF_results.xprs | cut -f 2,5 > EQF_transcripts_totCounts.tsv
cat ./CGT/CGT_express_outdir/CGT_results.xprs | cut -f 2,5 > CGT_transcripts_totCounts.tsv
cat ./CGF/CGF_express_outdir/CGF_results.xprs | cut -f 2,5 > CGF_transcripts_totCounts.tsv

#merge the totcounts into a single count matrix file in R perform the following
CGF_transcript_totcounts <- read.csv("/home/dnarules/Anuj/ProjectAnemone/countmatrix/CGF_transcripts_totCounts.tsv",sep="\t")
CGT_transcript_totcounts <- read.csv("/home/dnarules/Anuj/ProjectAnemone/countmatrix/CGT_transcripts_totCounts.tsv",sep="\t")
EQF_transcript_totcounts <- read.csv("/home/dnarules/Anuj/ProjectAnemone/countmatrix/EQF_transcripts_totCounts.tsv",sep="\t")
EQT_transcript_totcounts <- read.csv("/home/dnarules/Anuj/ProjectAnemone/countmatrix/EQT_transcripts_totCounts.tsv",sep="\t")
CG_countMatrix <- merge(CGF_transcript_totcounts,CGT_transcript_totcounts,by="target_id",suffixes=c("CGF","CGT"))
EQ_countMatrix <- merge(EQF_transcript_totcounts,EQT_transcript_totcounts,by="target_id",suffixes=c("EQF","EQT"))
countMatrix <- merge(CG_countMatrix,EQ_countMatrix,by="target_id",,suffixes=c("",""))	
write.csv(countMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/countmatrix/transcript_TotCount_Matrix.tsv")
write.csv(CG_countMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/countmatrix/CG_transcript_TotCount_Matrix.tsv")
write.csv(EQ_countMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/countmatrix/EQ_transcript_TotCount_Matrix.tsv")

#to get only the swissprotname and the transcript Id from the annotation matrix
cut -f 3 EQ_trinotate_annotation_report.tsv | cut -f 1 -d '^' > EQ_swissprot_blastx_OnlyProteinNames_annotation.tsv
cut -f 2 EQ_trinotate_annotation_report.tsv > EQ_swissprot_blastx_OnlyTranscriptNames_annotation.tsv
paste EQ_swissprot_blastx_OnlyTranscriptNames_annotation.tsv EQ_swissprot_blastx_OnlyProteinNames_annotation.tsv > EQ_swissprot_blastx_ProteinAndTranscript.tsv
cut -f 3 CG_trinotate_annotation_report.tsv | cut -f 1 -d '^' > CG_swissprot_blastx_OnlyProteinNames_annotation.tsv
cut -f 2 CG_trinotate_annotation_report.tsv > CG_swissprot_blastx_OnlyTranscriptNames_annotation.tsv
paste CG_swissprot_blastx_OnlyTranscriptNames_annotation.tsv CG_swissprot_blastx_OnlyProteinNames_annotation.tsv > CG_swissprot_blastx_ProteinAndTranscript.tsv

#to merge the annotation and the count matrix
CG_transcript_totcounts <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CG_transcript_TotCount_Matrix.tsv")
CG_transcript_annotation <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CG_swissprot_blastx_ProteinAndTranscript.tsv",sep="\t")
colnames(CG_transcript_totcounts)[2] <- "transcript_id"
CG_AnnotatedCountMatrix <- merge(CG_transcript_totcounts,CG_transcript_annotation,by="transcript_id")
CG_AnnotatedCountMatrix[,-2] -> CG_AnnotatedCountMatrix
data.frame(CG_AnnotatedCountMatrix[,1],CG_AnnotatedCountMatrix[,4],CG_AnnotatedCountMatrix[,2],CG_AnnotatedCountMatrix[,3]) -> CG_AnnotatedCountMatrix
colnames(CG_AnnotatedCountMatrix) <- c("transcript_id","sprot_top_blastx_hit","tot_countsCGF","tot_countsCGT")
# dont remove the dot's because it can give you the transcripts that are high expressed but not homologous to any other known protein
CG_AnnotatedCountMatrix <- subset(CG_AnnotatedCountMatrix, sprot_top_blastx_hit!=".")
write.table(CG_AnnotatedCountMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CG_AnnotatedCountMatrix.tsv")


EQ_transcript_totcounts <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQ_transcript_TotCount_Matrix.tsv")
EQ_transcript_annotation <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQ_swissprot_blastx_ProteinAndTranscript.tsv",sep="\t")
colnames(EQ_transcript_totcounts)[2] <- "transcript_id"
EQ_AnnotatedCountMatrix <- merge(EQ_transcript_totcounts,EQ_transcript_annotation,by="transcript_id")
EQ_AnnotatedCountMatrix[,-2] -> EQ_AnnotatedCountMatrix
data.frame(EQ_AnnotatedCountMatrix[,1],EQ_AnnotatedCountMatrix[,4],EQ_AnnotatedCountMatrix[,2],EQ_AnnotatedCountMatrix[,3]) -> EQ_AnnotatedCountMatrix
colnames(EQ_AnnotatedCountMatrix) <- c("transcript_id","sprot_top_blastx_hit","tot_countsEQF","tot_countsEQT")
# dont remove the dot's because it can give you the transcripts that are high expressed but not homologous to any other known protein
EQ_AnnotatedCountMatrix <- subset(EQ_AnnotatedCountMatrix, sprot_top_blastx_hit!=".")
write.table(EQ_AnnotatedCountMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQ_AnnotatedCountMatrix.tsv")


#########################################################
#repeat the whole thing for TPM and not CPM
#extract the TPM for each
cat ./EQT/EQT_express_outdir/EQT_results.xprs | cut -f 2,15 > EQT_transcripts_TPM.tsv
cat ./EQF/EQF_express_outdir/EQF_results.xprs | cut -f 2,15 > EQF_transcripts_TPM.tsv
cat ./CGT/CGT_express_outdir/CGT_results.xprs | cut -f 2,15 > CGT_transcripts_TPM.tsv
cat ./CGF/CGF_express_outdir/CGF_results.xprs | cut -f 2,15 > CGF_transcripts_TPM.tsv
#merge the totcounts into a single count matrix file in R perform the following
CGF_transcript_TPM <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CGF_transcripts_TPM.tsv",sep="\t")
CGT_transcript_TPM <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CGT_transcripts_TPM.tsv",sep="\t")
EQF_transcript_TPM <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQF_transcripts_TPM.tsv",sep="\t")
EQT_transcript_TPM <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQT_transcripts_TPM.tsv",sep="\t")
CG_countMatrix <- merge(CGF_transcript_TPM,CGT_transcript_TPM,by="target_id",suffixes=c("CGF","CGT"))
EQ_countMatrix <- merge(EQF_transcript_TPM,EQT_transcript_TPM,by="target_id",suffixes=c("EQF","EQT"))
countMatrix <- merge(CG_countMatrix,EQ_countMatrix,by="target_id",,suffixes=c("",""))	
write.csv(countMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/transcript_TPM_Matrix.tsv")
write.csv(CG_countMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CG_transcript_TPM_Matrix.tsv")
write.csv(EQ_countMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQ_transcript_TPM_Matrix.tsv")

#to get only the swissprotname and the transcript Id from the annottion matrix
cut -f 3 EQ_trinotate_annotation_report.tsv | cut -f 1 -d '^' > EQ_swissprot_blastx_OnlyProteinNames_annotation.tsv
cut -f 2 EQ_trinotate_annotation_report.tsv > EQ_swissprot_blastx_OnlyTranscriptNames_annotation.tsv
paste EQ_swissprot_blastx_OnlyTranscriptNames_annotation.tsv EQ_swissprot_blastx_OnlyProteinNames_annotation.tsv > EQ_swissprot_blastx_ProteinAndTranscript.tsv
cut -f 3 CG_trinotate_annotation_report.tsv | cut -f 1 -d '^' > CG_swissprot_blastx_OnlyProteinNames_annotation.tsv
cut -f 2 CG_trinotate_annotation_report.tsv > CG_swissprot_blastx_OnlyTranscriptNames_annotation.tsv
paste CG_swissprot_blastx_OnlyTranscriptNames_annotation.tsv CG_swissprot_blastx_OnlyProteinNames_annotation.tsv > CG_swissprot_blastx_ProteinAndTranscript.tsv

#to merge the annotation and the count matrix
CG_transcript_TPM <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CG_transcript_TPM_Matrix.tsv")
CG_transcript_annotation <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CG_swissprot_blastx_ProteinAndTranscript.tsv",sep="\t")
colnames(CG_transcript_TPM)[2] <- "transcript_id"
CG_AnnotatedCountMatrix <- merge(CG_transcript_TPM,CG_transcript_annotation,by="transcript_id")
CG_AnnotatedCountMatrix[,-2] -> CG_AnnotatedCountMatrix
data.frame(CG_AnnotatedCountMatrix[,1],CG_AnnotatedCountMatrix[,4],CG_AnnotatedCountMatrix[,2],CG_AnnotatedCountMatrix[,3]) -> CG_AnnotatedCountMatrix
colnames(CG_AnnotatedCountMatrix) <- c("transcript_id","sprot_top_blastx_hit","TPM_CGF","TPM_CGT")
# dont remove the dot's because it can give you the transcripts that are high expressed but not homologous to any other known protein
#CG_AnnotatedCountMatrix <- subset(CG_AnnotatedCountMatrix, sprot_top_blastx_hit!=".")
write.table(CG_AnnotatedCountMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/CG_AnnotatedCountMatrix_TPM.tsv")


EQ_transcript_TPM <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQ_transcript_TPM_Matrix.tsv")
EQ_transcript_annotation <- read.csv("/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQ_swissprot_blastx_ProteinAndTranscript.tsv",sep="\t")
colnames(EQ_transcript_TPM)[2] <- "transcript_id"
EQ_AnnotatedCountMatrix <- merge(EQ_transcript_TPM,EQ_transcript_annotation,by="transcript_id")
EQ_AnnotatedCountMatrix[,-2] -> EQ_AnnotatedCountMatrix
data.frame(EQ_AnnotatedCountMatrix[,1],EQ_AnnotatedCountMatrix[,4],EQ_AnnotatedCountMatrix[,2],EQ_AnnotatedCountMatrix[,3]) -> EQ_AnnotatedCountMatrix
colnames(EQ_AnnotatedCountMatrix) <- c("transcript_id","sprot_top_blastx_hit","TPM_EQF","TPM_EQT")
# dont remove the dot's because it can give you the transcripts that are high expressed but not homologous to any other known protein
#EQ_AnnotatedCountMatrix <- subset(EQ_AnnotatedCountMatrix, sprot_top_blastx_hit!=".")
write.table(EQ_AnnotatedCountMatrix,sep="\t",file="/home/dnarules/Anuj/ProjectAnemone/AnnotationAndCountMatrix/EQ_AnnotatedCountMatrix_TPM.tsv")


#download the toxins from the venomzone or uniprot and make a file of the IDS
cut -f 2 UniprotSodiumToxins.tab > UniprotSodiumToxinsProteinsOnly.tsv
# make a file of all the proteins you need for NaTx, KTx, CaTx, and other and ion channels
# search the ids in the annotatedcountmatrix and see which one are there.
grep -wFf UniprotSodiumToxinsProteinsOnly.tsv CG_AnnotatedCountMatrix_TPM.tsv
# note down the number of genes for that and 


#Differential expression USING TRINITY
########################################################
##########################################################3

# make the countmatrix for CGF CGT EQF and EQT
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left EQF_R1_paired_trimmed.fastq --right EQF_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left EQT_R1_paired_trimmed.fastq --right EQT_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left CGF_R1_paired_trimmed.fastq --right CGF_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left CGT_R1_paired_trimmed.fastq --right CGT_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir express_outdir
/home/dnarules/Trinity/util/abundance_estimates_to_matrix.pl --est_method express --out_prefix trans_counts --name_sample_by_basedir EQF/express_outdir/results.xprs EQT/express_outdir/results.xprs CGF/express_outdir/results.xprs CGT/express_outdir/results.xprs
/home/dnarules/Trinity/util/abundance_estimates_to_matrix.pl --est_method express --out_prefix trans_counts --name_sample_by_basedir EQF/express_outdir/results.xprs EQT/express_outdir/results.xprs CGF/express_outdir/results.xprs CGT/express_outdir/results.xprs

# make the countmatrix for CGF2 CGF2 CGF3 CGT1 CGT2 CGT3 EQF1 EQF2 EQF2 EQF3 and EQT1 EQT2 EQT3
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left ../raw/CGT1_R1_paired_trimmed.fastq --right ../raw/CGT1_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir CGT1express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left ../raw/CGT2_R1_paired_trimmed.fastq --right ../raw/CGT2_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir CGT2express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left ../raw/CGT3_R1_paired_trimmed.fastq --right ../raw/CGT3_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir CGT3express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left ../raw/CGF1_R1_paired_trimmed.fastq --right ../raw/CGF1_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir CGF1express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left ../raw/CGF2_R1_paired_trimmed.fastq --right ../raw/CGF2_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir CGF2express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts CG_Trinity.fasta --seqType fq --left ../raw/CGF3_R1_paired_trimmed.fastq --right ../raw/CGF3_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir CGF3express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left ../raw/EQT1_R1_paired_trimmed.fastq --right ../raw/EQT1_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir EQT1express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left ../raw/EQT2_R1_paired_trimmed.fastq --right ../raw/EQT2_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir EQT2express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left ../raw/EQT3_R1_paired_trimmed.fastq --right ../raw/EQT3_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir EQT3express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left ../raw/EQF1_R1_paired_trimmed.fastq --right ../raw/EQF1_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir EQF1express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left ../raw/EQF2_R1_paired_trimmed.fastq --right ../raw/EQF2_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir EQF2express_outdir
/home/dnarules/Trinity/util/align_and_estimate_abundance.pl --transcripts EQ_Trinity.fasta --seqType fq --left ../raw/EQF3_R1_paired_trimmed.fastq --right ../raw/EQF3_R2_paired_trimmed.fastq --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir EQF3express_outdir
/home/dnarules/Trinity/util/abundance_estimates_to_matrix.pl --est_method express --out_prefix trans_counts --name_sample_by_basedir EQF1express_outdir/results.xprs EQF2express_outdir/results.xprs EQF3express_outdir/results.xprs EQT1express_outdir/results.xprs EQT2express_outdir/results.xprs EQT3express_outdir/results.xprs
/home/dnarules/Trinity/util/abundance_estimates_to_matrix.pl --est_method express --out_prefix trans_counts --name_sample_by_basedir CGF1express_outdir/results.xprs CGF2express_outdir/results.xprs CGF3express_outdir/results.xprs CGT1express_outdir/results.xprs CGT2express_outdir/results.xprs CGT3express_outdir/results.xprs

# combine the count and the annotation matrix
/home/dnarules/Trinotate-3.0.2/util/Trinotate_get_feature_name_encoding_attributes.pl ./trinotate/CG/CG_trinotate_annotation_report.xls  > ./AnnotationAndCountMatrix2/CG_annot_feature_map.txt
/home/dnarules/Trinotate-3.0.2/util/Trinotate_get_feature_name_encoding_attributes.pl ./trinotate/EQ/EQ_trinotate_annotation_report.xls  > ./AnnotationAndCountMatrix2/EQ_annot_feature_map.txt
/home/dnarules/Trinity/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl ./AnnotationAndCountMatrix2/CG_trans_counts.counts.matrix ./AnnotationAndCountMatrix2/CG_annot_feature_map.txt > ./AnnotationAndCountMatrix2/CG_trans.counts.wAnnot.matrix
/home/dnarules/Trinity/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl ./AnnotationAndCountMatrix2/EQ_trans_counts.counts.matrix ./AnnotationAndCountMatrix2/EQ_annot_feature_map.txt > ./AnnotationAndCountMatrix2/EQ_trans.counts.wAnnot.matrix
/home/dnarules/Trinity/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl CG_trans_counts.TMM.EXPR.matrix CG_annot_feature_map.txt > GG_trans.TMM.EXPR.annotated.matrix

# ruun the diff expression part # use only the CG and EQ matrix
/home/dnarules/Trinity/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix ./edgeR_trinity/CG_trans.counts.wAnnot.matrix  --method edgeR  --samples_file ./edgeR_trinity/samplesCG.txt 
/home/dnarules/Trinity/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ./CG_trans.counts.wAnnot.matrix.CGF_vs_CGT.edgeR.count_matrix -P 1e-3 -C 2 --samples samplesCG.txt   --max_DE_genes_per_comparison 500  

# ruun the diff expression part # use only the CG and EQ matrix but only for the toxins
/home/dnarules/Trinity/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix ./CG_AnnotatedCountMatrix_totCounts.tsv  --method edgeR  --samples_file ./edgeR_trinity/samplesCG.txt 
/home/dnarules/Trinity/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ./CG_trans.counts.wAnnot.matrix.CGF_vs_CGT.edgeR.count_matrix -P 1e-3 -C 2 --samples samplesCG.txt   --max_DE_genes_per_comparison 500  


# run the gene set enrichment analysis part
/home/dnarules/Trinity/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ./CG_trans.TMM.EXPR.annotated.matrix -P 1e-3 -C 2 --samples samplesCG.txt  --max_genes_clust 100  -- also add in the gene set enrichment protocl 
/home/dnarules/Trinity/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R  diffExpr.P1e-3_C2.matrix.RData --Ptree 60


##############################################################
#############################################################

####################################
### EdgeR 
####################################
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("zebrafishRNASeq")
biocLite("edgeR")
biocLite("NBPSeq")
library(NBPSeq)
library(edgeR)

#####################################################################
### Analysis for CGFvsCGT ###
#####################################################################
setwd("~/Desktop/Anuj")
## reading the file from RUVg
CountMatrixCG <- read.csv("~/Desktop/Anuj/CG_AllCounts.csv", stringsAsFactors=FALSE)
View(CountMatrixCG)

targets <- read.csv("~/Desktop/Anuj/Design_matrixCG.csv")
Group <- factor(paste(targets$Treat,sep="."))
cbind(targets,Group=Group)
design <- model.matrix(~0+Group)
design
colnames(design) <- levels(Group)

#####################################################
### Making the count matrix with treatment and time
#####################################################
ModifiedCountMatrixCG <- DGEList(counts=CountMatrixCG[,2:7], group=Group,genes=CountMatrixCG[,1:2])
ModifiedCountMatrixCG
DispersedModifiedCountMatrixCG <- estimateDisp(ModifiedCountMatrixCG, design)
fitCG <- glmFit(DispersedModifiedCountMatrixCG, design)
my.contrasts <- makeContrasts(CGFvsCGT = (CGF)-(CGT),
                              levels=design)

####

lrt.CGFvsCGT <- glmLRT(fitCG, contrast=my.contrasts[,"CGFvsCGT"])
topTags(lrt.CGFvsCGT)
ArrangedlrtCGFvsCGT <- topTags(lrt.CGFvsCGT, n=500540)
GeneCG <- ArrangedlrtCGFvsCGT[[1]]$X
log2FoldChange <- ArrangedlrtCGFvsCGT[[1]]$logFC
pvalue <- ArrangedlrtCGFvsCGT[[1]]$PValue
padj <- ArrangedlrtCGFvsCGT[[1]]$FDR
resCG <- data.frame(GeneCG,log2FoldChange,pvalue,padj)
significantresCGFvsCGT<-subset(resCG,padj<0.05)
DiffExpGenesCGFvsCGT<-significantresCGFvsCGT$Gene

##############################
### Volcano plots 
##############################

# Make a basic volcano plot
with(resCG, plot(log2FoldChange, -log10(pvalue), col = "grey", pch=20, main="CGFvsCGT", xlim=c(-25,25),ylim=c(0,28)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resCG, padj<.05  & log2FoldChange>3), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resCG, padj<.05  & log2FoldChange<(-3)), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
################################################################################
### gives you the differntailly expressed genes in the whole datset 
#################################################################################

DiffExpGenesAllCG<-c(as.character(DiffExpGenesCGFvsCGT))
DiffExpGenesAllCGunique<-data.frame(unique(DiffExpGenesAllCG))
head(DiffExpGenesAllCGunique)
write.csv(DiffExpGenesAllCGunique,file="~/Desktop/Anuj/DiffExpGenesAll_CGFvsCGT.csv")

#####################################################################
### Analysis for EQFvsEQT ###
#####################################################################

setwd("~/Desktop/Anuj")
## reading the file from RUVg
CountMatrixEQ <- read.csv("~/Desktop/Anuj/EQ_AllCounts.csv", stringsAsFactors=FALSE)
View(CountMatrixEQ)

targets <- read.csv("~/Desktop/Anuj/Design_matrixEQ.csv")
Group <- factor(paste(targets$Treat,sep="."))
cbind(targets,Group=Group)
design <- model.matrix(~0+Group)
design
colnames(design) <- levels(Group)

#####################################################
### Making the count matrix with treatment and time
#####################################################
ModifiedCountMatrixEQ <- DGEList(counts=CountMatrixEQ[,2:7], group=Group,genes=CountMatrixEQ[,1:2])
ModifiedCountMatrixEQ
DispersedModifiedCountMatrixEQ <- estimateDisp(ModifiedCountMatrixEQ, design)
fitEQ <- glmFit(DispersedModifiedCountMatrixEQ, design)
my.contrasts <- makeContrasts(EQFvsEQT = (EQF)-(EQT),
                              levels=design)


lrt.EQFvsEQT <- glmLRT(fitEQ, contrast=my.contrasts[,"EQFvsEQT"])
topTags(lrt.EQFvsEQT)
ArrangedlrtEQFvsEQT <- topTags(lrt.EQFvsEQT, n=849496)
GeneEQ <- ArrangedlrtEQFvsEQT[[1]]$X
log2FoldChange <- ArrangedlrtEQFvsEQT[[1]]$logFC
pvalue <- ArrangedlrtEQFvsEQT[[1]]$PValue
padj <- ArrangedlrtEQFvsEQT[[1]]$FDR
resEQ <- data.frame(GeneEQ,log2FoldChange,pvalue,padj)
significantresEQFvsEQT<-subset(resEQ,padj<0.05)
DiffExpGenesEQFvsEQT<-significantresEQFvsEQT$Gene

##############################
### Volcano plots 
##############################

# Make a basic volcano plot
with(resEQ, plot(log2FoldChange, -log10(pvalue), col = "grey", pch=20, main="EQFvsEQT", xlim=c(-20,20),ylim=c(0,15)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resEQ, padj<.05  & log2FoldChange>3), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resEQ, padj<.05  & log2FoldChange<(-3)), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
################################################################################
### gives you the differntailly expressed genes in the whole datset 
#################################################################################

DiffExpGenesAllEQ<-c(as.character(DiffExpGenesEQFvsEQT))
DiffExpGenesAlluniqueEQ<-data.frame(unique(DiffExpGenesAllEQ))
head(DiffExpGenesAlluniqueEQ)
write.csv(DiffExpGenesAlluniqueEQ,file="~/Desktop/Anuj/DiffExpGenesAll_EQFvsEQT.csv")


##################################################################################################
###############################################################################################
########################################################################################
#######################33Top100GEnes in Heatmap #############################################

####################################
### EdgeR 
####################################
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("zebrafishRNASeq")
biocLite("edgeR")
biocLite("NBPSeq")
library(NBPSeq)
library(edgeR)

#####################################################################
### Analysis for CGFvsCGT ###
#####################################################################
setwd("~/Desktop/Anuj")
## reading the file from RUVg
CountMatrixCG <- read.csv("~/Desktop/Anuj/CG_AllCounts.csv", stringsAsFactors=FALSE)
View(CountMatrixCG)

targets <- read.csv("~/Desktop/Anuj/Design_matrixCG.csv")
Group <- factor(paste(targets$Treat,sep="."))
cbind(targets,Group=Group)
design <- model.matrix(~0+Group)
design
colnames(design) <- levels(Group)

#####################################################
### Making the count matrix with treatment and time
#####################################################
ModifiedCountMatrixCG <- DGEList(counts=CountMatrixCG[,2:7], group=Group,genes=CountMatrixCG[,1:2])
ModifiedCountMatrixCG
DispersedModifiedCountMatrixCG <- estimateDisp(ModifiedCountMatrixCG, design)
fitCG <- glmFit(DispersedModifiedCountMatrixCG, design)
my.contrasts <- makeContrasts(CGFvsCGT = (CGF)-(CGT),
                              levels=design)

####

lrt.CGFvsCGT <- glmLRT(fitCG, contrast=my.contrasts[,"CGFvsCGT"])
topTags(lrt.CGFvsCGT)
ArrangedlrtCGFvsCGT <- topTags(lrt.CGFvsCGT, n=100)
GeneCG <- ArrangedlrtCGFvsCGT[[1]]$X
log2FoldChange <- ArrangedlrtCGFvsCGT[[1]]$logFC
pvalue <- ArrangedlrtCGFvsCGT[[1]]$PValue
padj <- ArrangedlrtCGFvsCGT[[1]]$FDR
resCG <- data.frame(GeneCG,log2FoldChange,pvalue,padj)
significantresCGFvsCGT<-subset(resCG,padj<0.05)
DiffExpGenesCGFvsCGT<-significantresCGFvsCGT$Gene

################################################################################
### gives you the differntailly expressed genes in the whole datset 
#################################################################################

DiffExpGenesAllCG<-c(as.character(DiffExpGenesCGFvsCGT))
DiffExpGenesAllunique<-data.frame(unique(DiffExpGenesAllCG))
head(DiffExpGenesAllunique)
write.csv(DiffExpGenesAllunique,file="~/Desktop/Anuj/DiffExpGenesTop100_CGFvsCGT.csv")

##################### heat map for CGFvsCGT ##########################

library(gplots)
library(RColorBrewer) 
setwd("~/Desktop/Anuj/CG")

Top100DEGsCG <- read.csv("~/Desktop/Anuj/CG/AverageTop100DEGs_CG.csv", stringsAsFactors=FALSE)
Top100DEGsCGNames <- Top100DEGsCG[,1]
Top100DEGsCG = data.frame(Top100DEGsCG, row.names =1 )
head(Top100DEGsCG)
Top100DEGsCG_matrix = log(Top100DEGsCG[1:2]+1, 2)
head(Top100DEGsCG)
FinalMatrix4Heatmap_CG <- as.matrix(Top100DEGsCG_matrix)
head(FinalMatrix4Heatmap_CG)
#Top100DEGsCGGenes <-Top100DEGsCG[,2:3]

rownames(Top100DEGsCG) <- Top100DEGsCGNames
head(Top100DEGsCG)

library(heatmap.plus)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
cols <- palette(brewer.pal(5, "Dark2"))

tiff ("HeatmapTop100_CGT.tiff", height = 7, width = 7, units = 'in', res = 1000, compression = 'lzw')
heatmap.2(FinalMatrix4Heatmap_CG, trace="none", Colv  = F, dendrogram = c("row"), col=hmcol,sepcolor="black", cexCol = .9, cexRow = .6, key = TRUE, keysize = 0.8, key.xlab = "log2FoldChange", key.title = NA,  key.ylab  = NA, density.info=c("none") , scale="row", symkey=FALSE)

dev.off()

##################### heat map for CG Toxins ##########################

library(gplots)
library(RColorBrewer) 
setwd("~/Desktop/Anuj/CG/CG_Toxin")

ToxinsDEGsCG <- read.csv("~/Desktop/Anuj/CG/CG_Toxin/Average_DEGs4toxin_CG.csv", stringsAsFactors=FALSE)
View(ToxinsDEGsCG)
ToxinsDEGsCGNames <- ToxinsDEGsCG[,1]
ToxinsDEGsCG = data.frame(ToxinsDEGsCG, row.names =1 )
head(ToxinsDEGsCG)
ToxinsDEGsCG_matrix = log(ToxinsDEGsCG[1:2]+1, 2)
head(ToxinsDEGsCG_matrix)
FinalMatrix4ToxinsDEGsCG <- as.matrix(ToxinsDEGsCG_matrix)
head(FinalMatrix4ToxinsDEGsCG)

library(heatmap.plus)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
cols <- palette(brewer.pal(5, "Dark2"))

tiff ("HeatmapToxin_CG.tiff", height = 7, width = 7, units = 'in', res = 1000, compression = 'lzw')
heatmap.2(FinalMatrix4ToxinsDEGsCG, trace="none", Colv  = F, dendrogram = c("row"), col=hmcol,sepcolor="black", cexCol = .9, cexRow = .6, key = TRUE, keysize = 0.8, key.xlab = "log2FoldChange", key.title = NA,  key.ylab  = NA, density.info=c("none") , scale="row", symkey=FALSE)

dev.off()

##################### heat map for CG NeuroToxins ##########################

library(gplots)
library(RColorBrewer) 
setwd("~/Desktop/Anuj/CG/CG_neurotoxin")
library(heatmap.plus)

NeuroToxinsDEGsCG <- read.csv("~/Desktop/Anuj/CG/CG_neurotoxin/Average_DEGs4Neurotoxin_CG.csv", stringsAsFactors=FALSE)
View(NeuroToxinsDEGsCG)
NeuroToxinsDEGsCGNames <- NeuroToxinsDEGsCG[,1]
NeuroToxinsDEGsCG = data.frame(NeuroToxinsDEGsCG, row.names =1 )
head(NeuroToxinsDEGsCG)
NeuroToxinsDEGsCG_matrix = log(NeuroToxinsDEGsCG[1:2]+1, 2)
head(NeuroToxinsDEGsCG_matrix)
FinalMatrix4NeuroToxinsDEGsCG <- as.matrix(NeuroToxinsDEGsCG_matrix)
head(FinalMatrix4NeuroToxinsDEGsCG)


hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
cols <- palette(brewer.pal(5, "Dark2"))

tiff ("Heatmap4NeuroToxin_CG.tiff", height = 7, width = 7, units = 'in', res = 1000, compression = 'lzw')
heatmap.2(FinalMatrix4NeuroToxinsDEGsCG, trace="none", Colv  = F, dendrogram = c("row"), col=hmcol,sepcolor="black", cexCol = .9, cexRow = .6, key = TRUE, keysize = 0.8, key.xlab = "log2FoldChange", key.title = NA,  key.ylab  = NA, density.info=c("none") , scale="row", symkey=FALSE)

dev.off()
#####################################################################
### Analysis for EQFvsEQT ###
#####################################################################

setwd("~/Desktop/Anuj")
## reading the file from RUVg
CountMatrixEQ <- read.csv("~/Desktop/Anuj/EQ_AllCounts.csv", stringsAsFactors=FALSE)
View(CountMatrixEQ)

targets <- read.csv("~/Desktop/Anuj/Design_matrixEQ.csv")
Group <- factor(paste(targets$Treat,sep="."))
cbind(targets,Group=Group)
design <- model.matrix(~0+Group)
design
colnames(design) <- levels(Group)

#####################################################
### Making the count matrix with treatment and time
#####################################################
ModifiedCountMatrixEQ <- DGEList(counts=CountMatrixEQ[,2:7], group=Group,genes=CountMatrixEQ[,1:2])
ModifiedCountMatrixEQ
DispersedModifiedCountMatrixEQ <- estimateDisp(ModifiedCountMatrixEQ, design)
fitEQ <- glmFit(DispersedModifiedCountMatrixEQ, design)
my.contrasts <- makeContrasts(EQFvsEQT = (EQF)-(EQT),
                              levels=design)


lrt.EQFvsEQT <- glmLRT(fitEQ, contrast=my.contrasts[,"EQFvsEQT"])
topTags(lrt.EQFvsEQT)
ArrangedlrtEQFvsEQT <- topTags(lrt.EQFvsEQT, n=100)
GeneEQ <- ArrangedlrtEQFvsEQT[[1]]$X
log2FoldChange <- ArrangedlrtEQFvsEQT[[1]]$logFC
pvalue <- ArrangedlrtEQFvsEQT[[1]]$PValue
padj <- ArrangedlrtEQFvsEQT[[1]]$FDR
resEQ <- data.frame(GeneEQ,log2FoldChange,pvalue,padj)
significantresEQFvsEQT<-subset(resEQ,padj<0.05)
DiffExpGenesEQFvsEQT<-significantresEQFvsEQT$Gene

##############################################################################
### gives you the differntailly expressed genes in the whole datset 
#################################################################################

DiffExpGenesAllEQ<-c(as.character(DiffExpGenesEQFvsEQT))
DiffExpGenesAllunique<-data.frame(unique(DiffExpGenesAllEQ))
head(DiffExpGenesAllunique)
write.csv(DiffExpGenesAllunique,file="~/Desktop/Anuj/DiffExpGenesTop100_EQFvsEQT.csv")

##################### heat map for EQFvsEQT ##########################

library(gplots)
library(RColorBrewer) 
setwd("~/Desktop/Anuj/EQ")

Top100DEGsEQ <- read.csv("~/Desktop/Anuj/EQ/AverageTop100DEGs_EQ.csv", stringsAsFactors=FALSE)
View(Top100DEGsEQ)
Top100DEGsEQNames <- Top100DEGsEQ[,1]
Top100DEGsEQ = data.frame(Top100DEGsEQ, row.names =1 )
head(Top100DEGsEQ)
Top100DEGsEQ_matrix = log(Top100DEGsEQ[1:2]+1, 2)
head(Top100DEGsEQ)
FinalMatrix4Heatmap_EQ <- as.matrix(Top100DEGsEQ_matrix)
head(FinalMatrix4Heatmap_EQ)
#Top100DEGsEQGenes <-Top100DEGsEQ[,2:3]

rownames(Top100DEGsEQ) <- Top100DEGsEQNames
head(Top100DEGsEQ)
#heatmap(FinalMatrix4Heatmap_EQ, col = rev(heat.colors(256)),Colv = NA, cexRow = 0.6, keep.dendro = TRUE)
library(heatmap.plus)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
cols <- palette(brewer.pal(5, "Dark2"))

tiff ("HeatmapTop100_EQ.tiff", height = 7, width = 7, units = 'in', res = 1000, compression = 'lzw')
heatmap.2(FinalMatrix4Heatmap_EQ, trace="none", Colv  = F, dendrogram = c("row"), col=hmcol,sepcolor="black", cexCol = .9, cexRow = .6, key = TRUE, keysize = 0.8, key.xlab = "log2FoldChange", key.title = NA,  key.ylab  = NA, density.info=c("none") , scale="row", symkey=FALSE)

dev.off()

##################### heat map for EQ Toxins ##########################

library(gplots)
library(RColorBrewer) 
setwd("~/Desktop/Anuj/EQ/EQ_Toxin")

ToxinsDEGsEQ <- read.csv("~/Desktop/Anuj/EQ/EQ_Toxin/Average_DEGs4toxin_EQ.csv", stringsAsFactors=FALSE)
View(ToxinsDEGsEQ)
ToxinsDEGsEQNames <- ToxinsDEGsEQ[,1]
ToxinsDEGsEQ = data.frame(ToxinsDEGsEQ, row.names =1 )
head(ToxinsDEGsEQ)
ToxinsDEGsEQ_matrix = log(ToxinsDEGsEQ[1:2]+1, 2)
head(ToxinsDEGsEQ_matrix)
FinalMatrix4ToxinsDEGsEQ <- as.matrix(ToxinsDEGsEQ_matrix)
head(FinalMatrix4ToxinsDEGsEQ)

library(heatmap.plus)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
cols <- palette(brewer.pal(5, "Dark2"))

tiff ("HeatmapToxin_EQ.tiff", height = 7, width = 7, units = 'in', res = 1000, compression = 'lzw')
heatmap.2(FinalMatrix4ToxinsDEGsEQ, trace="none", Colv  = F, dendrogram = c("row"), col=hmcol,sepcolor="black", cexCol = .9, cexRow = .6, key = TRUE, keysize = 0.8, key.xlab = "log2FoldChange", key.title = NA,  key.ylab  = NA, density.info=c("none") , scale="row", symkey=FALSE)

dev.off()







