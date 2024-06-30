
# Download modkit

rm modkit_v0.3.0_centos7_x86_64.tar.gz
wget https://github.com/nanoporetech/modkit/releases/download/v0.3.0/modkit_v0.3.0_centos7_x86_64.tar.gz
tar xvf modkit_v0.3.0_centos7_x86_64.tar.gz

# Run modkit

for f in data/*.bam; do dist/modkit pileup --preset traditional --with-header --ref human_g1k_v37.fasta $f modkit_pileup/$(basename $f).bed -t 36; done

# Make a regions.bed file from Table 2
# 4       111100841       113912749
# 4       113915113       113954280
# X       29020059        33038509
# X       29570408        33426916
# 14      21501001        21879481
# 14      21502976        21849620
# 16      3688503 3689781
# 16      3762460 3899656

# Add 10 Mb padding and merge regions
samtools faidx human_g1k_v37.fasta
bedtools slop -i regions.bed -b 10000000 -g human_g1k_v37.fasta.fai > regions.padded.bed
bedtools merge -i regions.padded.bed > regions.padded.merged.bed


echo "ANK2\nIL1RAPL1\nCHD8\nCREBBP" > names.txt
paste regions.padded.merged.bed names.txt > regions.padded.merged.named.bed


echo "ANK2\nANK2\nIL1RAPL1\nIL1RAPL1\nCHD8\nCHD8\nCREBBP\nCREBBP" > names2.txt
paste regions.bed names2.txt > regions.named.bed

# Extract only relevant regions 

for f in modkit_pileup/*bam.bed;do 
  bedtools intersect -a $f -b regions.padded.merged.bed > $f.regions.padded.merged.bed &
done

# Convert to DSS format

for f in modkit_pileup/*bam.bed.regions.padded.merged.bed;do
  awk -v OFS='\t' '{print $1,$2,$10,$12}' $f > $f.dds
done
