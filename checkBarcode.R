## This script checks whether the reads in each splitted file belong to the same barcode
## and whether the filename contains the corresponding sample id name.
## Note, this script is a memory hog. It is only intended to check the idemp code.
##source("http://bioconductor.org/biocLite.R")
##biocLite("ShortRead")
require(ShortRead)

#### Files original and decoded ####
barcodeFile="barcode_sample.txt"
I1File="~/projects/141104_M03249_0008_000000000-ABY6R/Data/Intensities/BaseCalls/Undetermined_S0_L001_I1_001.fastq.gz"
R1File="~/projects/141104_M03249_0008_000000000-ABY6R/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz" 
R2File="~/projects/141104_M03249_0008_000000000-ABY6R/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz"
decodedFile=list.files("./", pattern="*.fastq.gz$")

#### Read in barcode table, sequence names and barcodes ####
firstLine=readLines("barcode_sample.txt",1)
barcodeTable=read.table(barcodeFile, header=grepl("arcode", firstLine))
### Read in seqquence names and barcodes ###
barcodeReads=readFastq(I1File)
I1Names=unname(sapply(as.character(id(barcodeReads)), function(x) (strsplit(x," ")[[1]][1] ) ) )
I1Barcodes=as.character(sread(barcodeReads))

#### Loop through decoded files ####
for( i in 1:length(decodedFile) ) {
  cat(i, " of ", length(decodedFile), "\n")
  cat(decodedFile[i],"\n")
  decodedReads=readFastq(decodedFile[i])
  seqNames=unname(sapply(as.character(id(decodedReads)), function(x) (strsplit(x," ")[[1]][1] ) ) )
  seqBarcodes=as.character(sread(decodedReads))
  
  idx=match(seqNames, I1Names)
  #table(I1Barcodes[idx])
  #sort(table(I1Barcodes[idx]), decreasing=T)
  bctable=sort(table(I1Barcodes[idx]), decreasing=T)
  if ( length(bctable)> 1000 ) {
    message(length(bctable), " barcodes found out of ", length(seqNames), " reads.") 
  }
  if ( length( grep(names(bctable)[1], barcodeTable[,1]) )==0 ) {
    message("most frequent barcode ",names(bctable)[1], " not matched.") 
    next
  }
  if ( length(bctable)> 1000 ) next
  
  cat(bctable,"\n")
  codeMax=names(bctable)[1]
  cat(codeMax,"\n")
  if( sort(table(I1Barcodes[idx]), decreasing=T)[1] < length(idx)*0.9 ) {
    message("suspicious decoding\n")
    message(sort(table(I1Barcodes[idx]), decreasing=T),"\n","Total reads:",length(idx),"\n")
  }
  sampleid = as.character(barcodeTable[barcodeTable[,1]==codeMax,2][1])
  print( sampleid )
  print( grepl(sampleid, decodedFile[i]) )  
  if ( ! grepl(sampleid, decodedFile[i]) ) 
    message("SampleID not in file name\n", sampleid, "\n", decodedFile[i], "\n")
  ##if ( i>10 ) break; 
}

q("no")
