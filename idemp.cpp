/**** system headers ****/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
#include <zlib.h>
//#include <sys/types.h>
//#include <cstdlib>
using namespace std;

/**** user headers ****/
#include "functions.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void check_are_read_names_same(string I1File, string R1File, string R2File, 
				 string& readBarcode)
{
readBarcode="";
  
  vector<string> inpFile(0);
  inpFile.push_back(I1File);
  inpFile.push_back(R1File);
  if ( R2File != "" ) inpFile.push_back(R2File);
  
  vector<gzFile> inpFH(inpFile.size(), Z_NULL);
  vector<kseq_t *> seq(inpFile.size());
  for(size_t i=0; i<inpFile.size(); ++i) {
    inpFH[i] = gzopen(inpFile[i].c_str(), "r");
if(inpFH[i]==Z_NULL){ cerr << "Cannot open " << inpFile[i] << endl; exit(1); }
seq[i] = kseq_init( inpFH[i] );
  }
  
  bool isReadNameSame=true;
  bool readDone=false;
  vector<int> l(inpFile.size());
  vector<string> sname(inpFile.size());
  size_t icount=0;
  while( ! readDone ) {
    readDone=true;
    for(size_t i=0; i<inpFile.size(); ++i) {
      l[i] = kseq_read( seq[i] );
      sname[i] = string(seq[i]->name.s);
      if ( l[i]>=0 ) readDone=false;
    }
    if ( readDone ) break;  // done reading all reads 
    ++icount;
    readBarcode += string(seq[0]->seq.s) + "\t";
    for(size_t i=1; i<inpFile.size(); ++i) {
      if ( sname[i] == sname[i-1] ) continue; 
      isReadNameSame=false;
      cerr << "Read names are different at read " << icount << " "
	   << "between files\n" << inpFile[i-1] << " " << inpFile[i] << endl
	   << sname[i-1] << "\n" << sname[i] << endl;
      break;
    }
  }
  
  cerr << "Read names are same:\t" << isReadNameSame << endl;
  
  for(size_t i=0; i<inpFile.size(); ++i) {
    gzclose( inpFH[i] );
    kseq_destroy( seq[i] );
  }
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main_usage() 
{
  cerr << "Usage:\n"
       << "   idemp -b code -I1 I1 -R1 R1 -R2 R2 -m n -o folder\n"
       << "\nOptions:\n"
       << "   code    barcode file, each line contains barcode\\tid \n"
       << "   I1      barcode fastq file, text or gzipped\n"
       << "   R1      read1 fastq file, text or gzipped\n"
       << "   R2      read2 fastq file, text or gzipped, optional\n"
       << "   n       allowed base mismatches, optional, default=1\n"
       << "   folder  output folder, optional, default=.\n"
       << "\nOutput:\n"
       << "   folder/R1.id.fastq.gz   #reads assigned to ids\n"
       << "   folder/R2.id.fastq.gz   #reads assigned to ids\n"
       << "   folder/I1.id            #read name to id\n"
       << "   folder/I1.id.stat       #barcode base error stat\n"
       << endl;
  return 0;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main(int argc, char* argv[])
{
  if ( argc < 2 ) return main_usage();
  
  int res;
  string mycommand="";
  string barcodeFile="", I1File="", R1File="", R2File="", folder=".";
  int nMismatch=1; 
  long int BUFFERSIZE=(int)600E6; // buffer size to hold reads, two buffers used
  
  /*  read in parameters
   *  barcodeFile: each line contains barcode TAB sampleid, ^# lines ignored
   *  I1File: index reads file in fastq format, text or gzipped
   *  R1File: read1 fastq file
   *  R2File: read2 fastq file
   *  folder: output folder, default to current .
   */
  vector<string> ARGV(0);
  for(int i=0;i<argc;++i) ARGV.push_back(string(argv[i]));
  for(int i=1;i<(int)ARGV.size();++i) {
#define _next2 ARGV[i]=""; ARGV[i+1]=""; continue;
#define _next1 ARGV[i]=""; continue;
    if ( ARGV[i]=="-b" ) { barcodeFile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-I1" ) { I1File=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-R1" ) { R1File=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-R2" ) { R2File=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-m" ) { nMismatch=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-o" ) { folder=ARGV[i+1]; _next2; }
  }
  cerr << "barcode:\t" << barcodeFile << "\n"
       << "Index reads:\t" << I1File << "\n"
       << "Read1 reads:\t" << R1File << "\n"
       << "Read2 reads:\t" << R2File << "\n"
       << "Output folder:\t" << folder << "\n"
       << endl;
  
  string cmd="mkdir -p " + folder;
  res = system( cmd.c_str() );
  
  /*  1. read in barcode and sample id from barcodeFile 
   *     barcode: vector
   *     sampleid: vector
   *     barcode and sampleid are row matched
   *  2. check edit distance between barcodes
   *  3. check conflicks
   */
  vector<string> barcode(0), sampleid(0);
  ifstream FIN(barcodeFile.c_str());
  if(!FIN){ cerr << "Cannot open " << barcodeFile << endl; exit(1); }
  string tmps1, tmps2;
  while( FIN >> tmps1 >> tmps2 ) {
    if ( tmps1[1] == '#' ) continue;
    if ( tmps1.find("arcode") != string::npos ) continue;
    barcode.push_back(tmps1);
    sampleid.push_back(tmps2);
  }
  FIN.close();
  if ( barcode.empty() ) {
    cerr << "no barcode found in file " << barcodeFile << endl;
    return 1;
  }
  for(size_t i=0; i<barcode.size(); ++i) 
    cerr << barcode[i] << "\t" << sampleid[i] << endl;
  cerr << "barcodes:\t" << barcode.size() << endl << endl;
  /* check edit distance beyween barcodes */
  int minEditDistance=barcode[1].length();
  vector<string> miniEdbarcodes(0);
  cerr << "Pairwise barcode edit distance:" << endl;
  for(size_t i=0; i<barcode.size(); ++i) {
    for(size_t j=0; j<=i; ++j) {
      int ed=edit_distance(barcode[i], barcode[j]);
      if ( j==0 ) cerr << ed;
      else cerr << " " << ed ;
      if ( minEditDistance >= ed && i>j ) {
	if (  minEditDistance > ed ) miniEdbarcodes.clear();
	minEditDistance = ed;
	miniEdbarcodes.push_back( barcode[i] + ":" +  barcode[j] );
      }
      if ( ed==0 && i<j ) {
	if ( sampleid[i] != sampleid[j] ) {
	  cerr << "one barcode points to two sampleids:\n"
	       << barcode[i] << "\t" << sampleid[i] << "\t" << sampleid[j]
	       << endl;
	}
      }
      
    }
    cerr << endl;
  }
  cerr << "\nClosest barcodes, editDistance=" << minEditDistance << endl;
  for(size_t i=0; i<miniEdbarcodes.size(); ++i) cerr << miniEdbarcodes[i] << endl;
  cerr << endl;
  /* check if any sampleids are same */
  for(size_t i=0; i<sampleid.size(); ++i) {
    for(size_t j=i+1; j<sampleid.size(); ++j) {
      if ( sampleid[i] == sampleid[j] ) {
	if ( barcode[i] != barcode[j] ) {
	  cerr << "two barcodes point to the same sampleid:\n"
	       << barcode[i] << "\t" << barcode[j] << "\t" << sampleid[i]
	       << endl;
	}
      }
    }
  }
  
  /* check if the sequence names are same
   * read in sequence barcodes into a single string
   * index the barcodes delimiters
   */
  string readBarcode=""; readBarcode.reserve(BUFFERSIZE+1);
  vector<size_t> readBarcode_qpos(0); readBarcode_qpos.push_back(0);
  check_are_read_names_same(I1File, R1File, R2File, readBarcode);
  for(size_t i=0; i<readBarcode.size(); ++i) 
    if ( readBarcode[i]=='\t' ) readBarcode_qpos.push_back(i+1);
  
  /* this is the main decoding part
   * if complete match is found return
   * if not, if only one single base mutation return 
   * if not, check edit distance
   * Note: more than 80% barcodes should perfectly match
   * Return: 
   *     readBarcodeIdx match index 
   *     readBarcodeMis mismatches, edit distance
   */
  vector<int> readBarcodeIdx( readBarcode_qpos.size()-1, barcode.size() );
  vector<int> readBarcodeMis( readBarcode_qpos.size()-1, barcode[0].size() );
  for(size_t i=0; i<readBarcode_qpos.size()-1; ++i) {
    string query=readBarcode.substr(readBarcode_qpos[i], 
				    readBarcode_qpos[i+1]-readBarcode_qpos[i]-1 );
    // cerr << query << endl;
    // string seqname=readName.substr(readName_qpos[i], 
    //				   readName_qpos[i+1]-readName_qpos[i]-1 );
    
    if ( (i+1)%1000000==0 ) cerr << i+1 << endl;
    if ( query.size() != barcode[0].size() ) continue;
    
    // check exact match
    bool ismatched=false;
    for(size_t j=0; j<barcode.size(); ++j) {
      if ( query != barcode[j] ) continue;
      ismatched=true;
      readBarcodeIdx[i]=j;
      readBarcodeMis[i]=0;    
      break;
    }
    if ( ismatched ) continue;
    
    int minMismatch=query.size();
    // check single base mutation
    for(size_t j=0; j<barcode.size(); ++j) {
      int mismatch=0;
      for(size_t k=0; k<query.size(); ++k) mismatch+=( query[k]!= barcode[j][k] );
      if ( mismatch < minMismatch ) {
	minMismatch=mismatch;
	readBarcodeIdx[i]=j;
	readBarcodeMis[i] = minMismatch;    
	if ( minMismatch<=1 ) continue;
      }
    }
    //if ( minMismatch<=2 ) continue;
    
    // check deletion or insertation of single base
    for(size_t j=0; j<barcode.size(); ++j) {
      int mismatch=edit_distance(query, barcode[j]);
      if ( mismatch < minMismatch ) {
	minMismatch=mismatch;
	readBarcodeIdx[i]=j;
	readBarcodeMis[i]=minMismatch;    
      }
    }
    
  }
  cerr << "Done matching barcodes" << endl;
  
  /* write sequence barcodes assignment
   * in the format "#barcode\tbarcodeIdx\tdesigned_barcode\tmismatch\n"
   * the assignments are determined by
   *     if misMatch > minEditDistance between barcodes undecoded
   *     if misMatch > nMismatch  undecoded
   *     minEditDistance is calculated above
   *     nMismatch defaults to 1 allowing for base error; can be set by user
   */
  string tmps=I1File+".decode";
  size_t p=tmps.rfind('/');
  string fileName= p==string::npos ? tmps: tmps.substr(p+1);
  fileName=folder+"/"+fileName;
  ofstream FOUT(fileName.c_str());
  if ( !FOUT ) { cerr << "Can't open " << fileName << endl; exit(0); }
  FOUT << "#seq\tbarcodeIdx\tbarcode\tmismatch\n";
  for(size_t i=0; i<readBarcode_qpos.size()-1; ++i) {
    string query=readBarcode.substr(readBarcode_qpos[i], 
				    readBarcode_qpos[i+1]-readBarcode_qpos[i]-1 );
    //string seqname=readName.substr(readName_qpos[i], 
    //				   readName_qpos[i+1]-readName_qpos[i]-1 );
    
    /*
      cerr << i+1 << "\t|" << query << "|\t"
      //<< seqname << "|\t"
      << readBarcodeIdx[i] << "\t"
	 << barcode[ readBarcodeIdx[i] ] << "\t"
	 << readBarcodeMis[i] << endl;
    */
    //FOUT << seqname << "\t" 
    if ( readBarcodeMis[i] <= nMismatch){
    FOUT << query << "\t"
	 << readBarcodeIdx[i] << "\t"
	 << barcode[ readBarcodeIdx[i] ] << "\t"
	 << readBarcodeMis[i] << "\n";}
    
    if ( readBarcodeMis[i] > minEditDistance ) readBarcodeIdx[i]=barcode.size();
    if ( readBarcodeMis[i] > nMismatch ) readBarcodeIdx[i]=barcode.size();
  }
  FOUT.close();
  
  /* check which barcodes were used, create files to be written
   * check number of reads for each barcode/file
   * Note: the last index [barcode.size()] is for undecoded
   */
  vector<bool> barcodeWrite(barcode.size()+1, false);
  vector<size_t> readsPerBarcode(barcode.size()+1, 0);
  for(size_t i=0; i<readBarcodeIdx.size(); ++i) {
    barcodeWrite[ readBarcodeIdx[i] ] = true;
    readsPerBarcode[ readBarcodeIdx[i] ] +=1;
  }
  fileName += ".stat";
  FOUT.open(fileName.c_str());
  if ( !FOUT ) { cerr << "Can't open " << fileName << endl; exit(0); }
  FOUT << "#barcode\tsampleid\tfreq\n";
  for(size_t i=0; i<barcode.size(); ++i) 
    FOUT << barcode[i] << "\t" << sampleid[i] << "\t" << readsPerBarcode[i] << endl;
  FOUT << "undecoded\tunsigned\t" << readsPerBarcode[barcode.size()] << endl;
  FOUT << "All\tTotal\t" << readBarcodeIdx.size() << endl;
  FOUT.close();
  
  /* Loop through read1 and read2 files
   * write to decoded files
   */
  vector<string> inputReadFile(1, R1File);
  if ( R2File !="" ) inputReadFile.push_back(R2File) ;
  
  gzFile fp;
  kseq_t *seq;
  int l;
  for(int iFile=0; iFile<inputReadFile.size(); ++iFile) {
    
    /* set file names and file handles 
     * note: OS may limit number of open files, probably 512 or 1024.
     *       May need to rewrite this part for large number of barcodes. 
     */
    vector<string> outputNames(barcode.size()+1, "");
    vector<gzFile> outputfp(barcode.size()+1, Z_NULL);
    vector<size_t> readsWrittenPerBarcode(barcode.size()+1, 0);
    tmps=inputReadFile[iFile];
    p=tmps.rfind('/');
    fileName= p==string::npos ? tmps: tmps.substr(p+1);
    fileName=folder+"/"+fileName;
    for(size_t i=0; i<barcode.size(); ++i) {
      if ( ! barcodeWrite[i] )  continue;
      outputNames[i] = fileName + "_" + sampleid[i] + ".fastq.gz";
    }
    outputNames[barcode.size()]=fileName + "_unsigned.fastq.gz";
    for(size_t i=0; i<=barcode.size(); ++i) {
      if ( ! barcodeWrite[i] )  continue;
      cerr << i << "\t"
	   << ( i<barcode.size() ? barcode[i] : "undecoded" ) << "\t"
	   << outputNames[i] << endl;
      outputfp[i] = gzopen(outputNames[i].c_str(), "wb");
      if ( outputfp[i] == Z_NULL ) {
	cerr << "Can't open " << outputNames[i] << endl;
	exit(0);
      }
    }
    
    int l;
    fp = gzopen(inputReadFile[iFile].c_str(), "r");
    if(fp==Z_NULL) { 
      cerr << "Can't open " << inputReadFile[iFile] << endl; 
      exit(1); 
    }
    seq = kseq_init(fp);
    size_t icount=0;
    while ((l = kseq_read(seq)) >= 0) {
      //string seqname = string(seq->name.s);
      
      int ibc = readBarcodeIdx[icount];
      if ( ibc<0 ) ibc=barcode.size();
      
      tmps = "@"+string(seq->name.s);
      if (seq->comment.l) tmps += " " + string(seq->comment.s);
      tmps += "\n";
      tmps += string(seq->seq.s)+"\n+\n";
      if (seq->qual.l) tmps += string(seq->qual.s)+"\n";
      gzwrite( outputfp[ibc], &tmps[0], tmps.size() );
      ++readsWrittenPerBarcode[ ibc ];
      
      ++icount;
      if ( icount%1000000==0 ) cerr << icount << endl;
    }
    
    gzclose(fp); //close input
    for(size_t i=0; i<=barcode.size(); ++i) gzclose(outputfp[i]); //close output
    for( size_t i=0; i<barcode.size()+1; ++i) 
      if ( readsWrittenPerBarcode[i] != readsPerBarcode[i] ) {
	cerr << "Something wrong:\n"
	     << readsPerBarcode[i] << " reads should be written to " 
	     << outputNames[i] << "\n"
	     << readsWrittenPerBarcode[i] << " reads were written\n"
	     << endl;
      } // double check number of reads
    
  }
  
  kseq_destroy(seq);
  return 0;
} 
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

