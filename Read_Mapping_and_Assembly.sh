

#Author: Oumie Kuyateh 


#Date written: 2023


#Aim: I aim to quantify viral infections and investigate host-virus interaction using the transcriptomes of adult Drosophila that were initially sequenced #for other purposes, but which were also coincidentally infected with natural Drosophila viruses. 
Because these constitute entirely natural infections, #this study could help reveal the true nature of host gene expression in response to native viruses. 


#Data: I obtained nine large RNA-seq datasets of D. melanogaster from the European Nucleotide Archive (European Nucleotide Archive, 2022) that each #comprised 130 to over 900 sequencing libraries. 
The Project codes PRJNA258012, PRJNA261333, PRJNA281652, PRJNA305983, PRJNA325107, PRJNA483441 #PRJNA518903, PRJNA527284, and PRJNA75285 were used. 


#STAR Version 2.7.3a
#DIAMOND Version 2.0.11
#Trinity Version 2.9.1

##This is a bash code that shows how the datasets were mapped to the D. melanogaster genome and denovo builds by Trinity carried out. This was followed by #a DIAMOND BLASTX search to find virus sequences. 

###Read each line containing the Aspera link to download the datasets. The Aspera links are stored in Aspera_links.txt 

while read line;

####do the follwing

do

####Donwload file

	ascp -QT -l 300m -P33001 -i /localdisk/home/s1790958/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:$line  ./ 

####Get the SRR number only from the Aspera link and call it stem

	stem=$(echo $line | sed 's/.*\///g' | sed 's/\.fastq.gz*//g')

###print the stem on the screen

	echo $stem

###Make a while loop to test if the fastq file is broken and continue removing and downloading it until it is no longer broken
###The while loop will keep running until it is no longer true, i.e. the file is not broken

	while !  gunzip  -v -t ${stem}.fastq.gz
	do
		rm ${stem}.fastq.gz
		ascp -QT -l 300m -P33001 -i /localdisk/home/s1790958/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:$line  ./ 
	done

###If file is not broken run STAR

###Map each read to the annotated Drosophila genome using "STAR" but don't sort the bam file that contains the mapped reads 
###so that we can sort it using "samtools" and create a stats file

	STAR --runThreadN 12 --genomeDir /mnt/drive2-3tb/shi_analysis_real/D.melanogaster/Drosophila_STAR_Reference --outReadsUnmapped Fastx --outSAMtype BAM 	  Unsorted --readFilesCommand zcat --readFilesIn ${stem}.fastq.gz --outFileNamePrefix ${stem}.Dmel.Yanzhu.

###Sort the bam file that contains reads that mapped/aligned to Dmel using "samtools sort" commands. 
###The "sort" command positions the alignemts based on the order they occur in the chromosome

	samtools sort ${stem}.Dmel.Yanzhu.Aligned.out.bam -@ 6 -m 6G -O bam -o ${stem}.Dmel.Yanzhu.Aligned.sortedByCoord.out.bam 

###Then index the sorted bam file which allows one or a program to go through files that overlap a genomic region quickly.
###You can still get stats from an unindexed file but it will be slower

	samtools index ${stem}.Dmel.Yanzhu.Aligned.sortedByCoord.out.bam  

###Now generate a stats file using the "samtools idxstats" command

	samtools idxstats ${stem}.Dmel.Yanzhu.Aligned.sortedByCoord.out.bam  > ${stem}.Dmel.Yanzhu.stats

###Now move the Dmel stats file to another directory called "Dmel_Mapped_Stats_Yanzhu"

	mv ${stem}.Dmel.Yanzhu.stats /mnt/drive2-3tb/Yanzhu_analysis/D.melanogaster_Yanzhu/Dmel_Mapped_Stats_Yanzhu

###Now generate a count matrix using the "featureCounts" command 

	featureCounts -T 12 -a /mnt/drive2-3tb/shi_analysis_real/D.melanogaster/DrosophilaReferenceFromFlybase/*gtf -o /mnt/drive2-3tb/Yanzhu_analysis/${stem}.featureCounts.Yanzhu.count ${stem}.Dmel.Yanzhu.Aligned.sortedByCoord.out.bam 

###Move the count file and summary file from feature coutns to their designated directories using "mv" command

	mv ${stem}.featureCounts.Yanzhu.count /mnt/drive2-3tb/Yanzhu_analysis/D.melanogaster_Yanzhu/Dmel_featureCounts_Yanzhu/Dmel_featureCounts_Yanzhu_Count_matrix

	mv ${stem}.featureCounts.Yanzhu.count.summary /mnt/drive2-3tb/Yanzhu_analysis/D.melanogaster_Yanzhu/Dmel_featureCounts_Yanzhu/Dmel_featureCounts_Yanzhu_Count_summary

###Rename the reads that did not map to Dmel from "*.Unmapped.out.mate1" to "*.Unmapped.fastq*

	mv ${stem}.Dmel.Yanzhu.Unmapped.out.mate1 ${stem}.Dmel.Yanzhu.Unmapped.fastq

###Now map the reads the Non_Dmel reads to the Bowtie database (bowtie index) file that contains the common non-virus microbiomw
###in Dmel, D. simulans and xxxxx 

	bowtie2  -p 18  --un ${stem}.not_holobiont.fastq -x /mnt/drive2-3tb/Yanzhu_analysis/Holobiont_Bowtie_Index_Database/HOLOBIONT -U ${stem}.Dmel.Yanzhu.Unmapped.fastq 2> ${stem}.holobiont.log | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > ${stem}.holobiont.bam

###Now index it and generate a stats file using samtools

	samtools index ${stem}.holobiont.bam
	samtools idxstats ${stem}.holobiont.bam > ${stem}.holobiont.stats

###Use trinity to build the reads. Trinity brings out a lot of files but the only one I'm interested in is the trinity.fas file

	Trinity --seqType fq --max_memory 10G --CPU 12 --single ${stem}.not_holobiont.fastq   --output /mnt/drive2-3tb/Yanzhu_analysis/Denovo_Build_Yanzhu/${stem}.Dmel.Yanzhu.Trinity

###Look into the file in the directory containing the trinity.fas file using the "cat" command and use "fasta_formatter" 
###to unwrap the sequences and replace the word "TRINITY" in  each sequence name in the fasta file to the "stem" name using the "sed" command.
###This will enable us to know where each transcript-like sequence is coming from and save this into a file called "stem.TrinityBuild.fas"

	cat /mnt/drive2-3tb/Yanzhu_analysis/Denovo_Build_Yanzhu/${stem}.Dmel.Yanzhu.Trinity/Trinity.fasta | fasta_formatter | sed  "s/TRINITY/${stem}/g" > ${stem}.TrinityBuild.fas

###Now delete the directory that was obtained from trinity using the "rm -r" command since we no longer need it

	rm -r /mnt/drive2-3tb/Yanzhu_analysis/Denovo_Build_Yanzhu/${stem}.Dmel.Yanzhu.Trinity

###Now only keep seqences in the "stem.TrinityBuild.fas" that are longer than 500 bases using the "grep" command and save this into 
###stem.Trinity500Plus.fas


	grep "[AGCTN]\{500,\}" -B 1 --no-group-separator ${stem}.TrinityBuild.fas > ${stem}.Trinity500Plus.fas

###Now remove near identical sequences using the "cd-hit-est" command and save it in a file called "stem.Trinity500Plus.filtered.fas"
###I have explained some of the options I used in the command below

	cd-hit-est -c 0.97 -i ${stem}.Trinity500Plus.fas  -p 1 -d 0 -b 3 -T 16 -o ${stem}.Trinity500Plus.filtered.fas

###Now do a diamond search of thr refseq database using the diamond command as explained below.
###Diamond is a sequence aligner and works just like NCBI-BLAST but its up to 500 to 20,000 X faster.It is suitable 
###for protein-protein search and DNA-protein search on long and short reads including contigs (e.g. from TRINITY) and assemblies
###diamond 


	diamond blastx --frameshift 15 --top 5 --evalue 1e-20 --threads 16 --db /localdisk/BLAST/NCBI_databases/refseq_protein_and_viruses_diamond.dmnd --outfmt 6 qseqid sseqid pident length mismatch evalue bitscore sscinames sskingdoms --query ${stem}.Trinity500Plus.filtered.fas >  ${stem}.blastx.tsv

###Now zip the fastq file that did not map to the holobiont and move it to  a directory call non-holobiont

	gzip ${stem}.not_holobiont.fastq 

	mv  ${stem}.not_holobiont.fastq.gz  /mnt/drive2-3tb/Yanzhu_analysis/Holobiont_Yanzhu/Non_holobiont

###Move the holobiont stats to a directory called Holobiont_stats_files

	mv  ${stem}.holobiont.stats     /mnt/drive2-3tb/Yanzhu_analysis/Holobiont_Yanzhu/Holobiont/Holobiont_stats_files

###Move the holobiont bam file to the Holobiont_bam_files directory

	mv ${stem}.holobiont.bam  /mnt/drive2-3tb/Yanzhu_analysis/Holobiont_Yanzhu/Holobiont/Holobiont_bam_files

###Move the tsv file from the diamond blastx to the Diamond_Blast_Yanzhu_output directory

	mv ${stem}.blastx.tsv  /mnt/drive2-3tb/Yanzhu_analysis/Diamond_Blast_Yanzhu_output

###Move the raw reads obtained from the ENA database to the Yanzhu_Raw_Fastq_reads directory

	mv ${stem}.fastq.gz      /mnt/drive2-3tb/Yanzhu_analysis/Yanzhu_Raw_Fastq_reads

###Move the fasta file obtained from the trinity builder to Trinity_fasta_files directory

	mv ${stem}.TrinityBuild.fas    /mnt/drive2-3tb/Yanzhu_analysis/Denovo_Build_Yanzhu/Trinity_Processes_Fasta_Yanzhu

#zip the file that did not map to Dmel genome when the STAR mapping was used and move it to Dmel_Yanzhu_Unmapped_Reads directory

	gzip ${stem}.Dmel.Yanzhu.Unmapped.fastq 

	mv ${stem}.Dmel.Yanzhu.Unmapped.fastq.gz  /mnt/drive2-3tb/Yanzhu_analysis/D.melanogaster_Yanzhu/Dmel_Yanzhu_Unmapped_Reads

###Move log.out files and the aligned BAM files obtained after mapping to Dmel to the Dmel_Yanzhu_Log.final_files 
### and Dmel_Yanzhu_Soorted_Bam_files directories respectively


	mv *.Log.final.out /mnt/drive2-3tb/Yanzhu_analysis/D.melanogaster_Yanzhu/Dmel_Yanzhu_Log.final_files

	mv *Aligned.sortedByCoord.out.bam /mnt/drive2-3tb/Yanzhu_analysis/D.melanogaster_Yanzhu/Dmel_Yanzhu_Soorted_Bam_files

###Move all the fasta files from the post trinity

	mv *.fas* /mnt/drive2-3tb/Yanzhu_analysis/Denovo_Build_Yanzhu/Trinity_Processed_Fasta_Yanzhu


###Remove the remaining files that I don't need

	rm *.out*   *.holobiont.log *bai*


###Close the while loop and tell it where to read the files

done < Aspera_links.txt




