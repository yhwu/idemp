idemp
=====

Barcode demultiplex for Illumina I1, R1, R2 fastq.gz files. 

Only for typical Illumina runs, where the barcode sequence reads are saved in the I1_*.fastq.gz files, and the first fields of the sequence names are exactly same for I1, R1, and/or R2 fastq.gz files. This little program works as follows:

1. compare sequence names in the I1, R1, R2 read files;
2. read in barcode sequence from the I1 file;
3. read in barcode and sample id table;
4. calculate minimum edit distance among the designed barcodes, nEd;
5. find exact match for each barcode sequence read;
6. for reads not fully matched, calculate the minimum edit distance with the designed barcodes;
7. assign a barcode sequence if the edit distance is smaller than a cutoff;
8. by default, the cutoff=min(n, nEd), n(n=1) can be set by user.

#### Functions to be added
1. Trim R1 reads to get rid of barcode ends.


#### Compile and test

1. git clone https://github.com/yhwu/idemp
2. cd idemp
3. make
4. make test

#### Usage
```
Usage:
   idemp -b code -I1 I1 -R1 R1 -R2 R2 -m n -o folder

Options:
   code    barcode file, each line contains barcode\tid
   I1      barcode fastq file, text or gzipped
   R1      read1 fastq file, text or gzipped
   R2      read2 fastq file, text or gzipped, optional
   n       allowed base mismatches, optional, default=1
   folder  output folder, optional, default=.

Output:
   folder/R1.id.fastq.gz   #reads assigned to ids
   folder/R2.id.fastq.gz   #reads assigned to ids
   folder/I1.id            #read name to id
   folder/I1.id.stat       #barcode base error stat
```
