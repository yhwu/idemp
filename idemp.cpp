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
#include <unordered_map>
//#include <sys/types.h>
#include <cstdlib>
//#include <stdlib>
using namespace std;

/**** user headers ****/
#include "functions.h"
#include "idemp_func.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main_usage() {
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
int main(int argc, char* argv[]) {
    if (argc < 2) return main_usage();

    int res;
    string mycommand = "";
    long int BUFFERSIZE = (int) 600E6; // buffer size to hold reads, two buffers used

    /*  read in parameters
     *  barcodeFile: each line contains barcode TAB sampleid, ^# lines ignored
     *  I1File: index reads file in fastq format, text or gzipped
     *  R1File: read1 fastq file
     *  R2File: read2 fastq file
     *  folder: output folder, default to current .
     */
    IdempArgument arg = parse_arguments(argc, argv);
    cerr << "barcode:\t" << arg.barcodeFile << "\n"
            << "Index reads:\t" << arg.I1File << "\n"
            << "Read1 reads:\t" << arg.R1File << "\n"
            << "Read2 reads:\t" << arg.R2File << "\n"
            << "Output folder:\t" << arg.outputFolder << "\n"
            << endl;

    string cmd = "mkdir -p " + arg.outputFolder;
    res = system(cmd.c_str());

    /*  1. read in barcode and sample id from barcodeFile 
     *     barcode: vector
     *     sampleid: vector
     *     barcode and sampleid are row matched
     *  2. check edit distance between barcodes
     *  3. check conflicks
     */
    vector<string> barcode(0), sampleid(0);
    read_barcode_sampleid(arg.barcodeFile, barcode, sampleid);

    for (size_t i = 0; i < barcode.size(); ++i) cerr << barcode[i] << "\t" << sampleid[i] << endl;
    cerr << "barcodes:\t" << barcode.size() << endl << endl;

    /* check edit distance beyween barcodes */
    int minEditDistance = min_edit_distance(barcode, sampleid);



    /* check if the sequence names are same
     * index the barcodes delimiters
     */
    string readBarcode = "";
    readBarcode.reserve(BUFFERSIZE + 1);
    check_are_read_names_same(arg.I1File, arg.R1File, arg.R2File, readBarcode);

    vector<size_t> readBarcode_qpos(0);
    readBarcode_qpos.push_back(0);
    for (size_t i = 0; i < readBarcode.size(); ++i)
        if (readBarcode[i] == '\t') readBarcode_qpos.push_back(i + 1);

    /* this is the main decoding part
     * if complete match is found return
     * if not, if only one single base mutation return 
     * if not, check edit distance
     * Note: more than 80% barcodes should perfectly match
     * Note: last readBarcode_qpos points to string end
     * Return: 
     *     readBarcodeIdx match index 
     *     readBarcodeMis mismatches, edit distance
     */
    vector<int> readBarcodeIdx(readBarcode_qpos.size() - 1, barcode.size());
    vector<int> readBarcodeMis(readBarcode_qpos.size() - 1, barcode[0].size());

    map_readbarcode(readBarcode, readBarcode_qpos, barcode, minEditDistance, readBarcodeIdx, readBarcodeMis);
    cerr << "Done matching barcodes" << endl;

    /* write sequence barcodes assignment
     * in the format "#barcode\tbarcodeIdx\tdesigned_barcode\tmismatch\n"
     * the assignments are determined by
     *     if misMatch > minEditDistance between barcodes undecoded
     *     if misMatch > nMismatch  undecoded
     *     minEditDistance is calculated above
     *     nMismatch defaults to 1 allowing for base error; can be set by user
     */
    string tmps = arg.I1File + ".decode";
    size_t p = tmps.rfind('/');
    string fileName = p == string::npos ? tmps : tmps.substr(p + 1);
    fileName = arg.outputFolder + "/" + fileName;
    ofstream FOUT(fileName.c_str());
    if (!FOUT) {
        cerr << "Can't open " << fileName << endl;
        exit(0);
    }
    FOUT << "#seq\tbarcodeIdx\tbarcode\tmismatch\n";
    for (size_t i = 0; i < readBarcode_qpos.size() - 1; ++i) {

        string query = readBarcode.substr(readBarcode_qpos[i],
                readBarcode_qpos[i + 1] - readBarcode_qpos[i] - 1);

        cerr << query << '\t' << readBarcodeMis[i] << endl;
        if (readBarcodeMis[i] <= arg.nMismatch) {
            FOUT << query << "\t"
                    << readBarcodeIdx[i] << "\t"
                    << barcode[ readBarcodeIdx[i] ] << "\t"
                    << readBarcodeMis[i] << "\n";
        }

        if (readBarcodeMis[i] >= minEditDistance) readBarcodeIdx[i] = barcode.size();
        if (readBarcodeMis[i] > arg.nMismatch) readBarcodeIdx[i] = barcode.size();
    }
    FOUT.close();

    /* check which barcodes were used, create files to be written
     * check number of reads for each barcode/file
     * Note: the last index [barcode.size()] is for undecoded
     */
    vector<bool> barcodeWrite(barcode.size() + 1, false);
    vector<size_t> readsPerBarcode(barcode.size() + 1, 0);
    for (size_t i = 0; i < readBarcodeIdx.size(); ++i) {
        barcodeWrite[ readBarcodeIdx[i] ] = true;
        readsPerBarcode[ readBarcodeIdx[i] ] += 1;
    }
    fileName += ".stat";
    FOUT.open(fileName.c_str());
    if (!FOUT) {
        cerr << "Can't open " << fileName << endl;
        exit(0);
    }
    FOUT << "#barcode\tsampleid\tfreq\n";
    for (size_t i = 0; i < barcode.size(); ++i)
        FOUT << barcode[i] << "\t" << sampleid[i] << "\t" << readsPerBarcode[i] << endl;
    FOUT << "undecoded\tunsigned\t" << readsPerBarcode[barcode.size()] << endl;
    FOUT << "All\tTotal\t" << readBarcodeIdx.size() << endl;
    FOUT.close();

    /* Loop through read1 and read2 files
     * write to decoded files
     */
    vector<string> inputReadFile(1, arg.R1File);
    if (arg.R2File != "") inputReadFile.push_back(arg.R2File);

    gzFile fp;
    kseq_t *seq;
    int l;
    for (int iFile = 0; iFile < inputReadFile.size(); ++iFile) {

        /* set file names and file handles 
         * note: OS may limit number of open files, probably 512 or 1024.
         *       May need to rewrite this part for large number of barcodes. 
         */
        vector<string> outputNames(barcode.size() + 1, "");
        vector<gzFile> outputfp(barcode.size() + 1, Z_NULL);
        vector<size_t> readsWrittenPerBarcode(barcode.size() + 1, 0);
        tmps = inputReadFile[iFile];
        p = tmps.rfind('/');
        fileName = p == string::npos ? tmps : tmps.substr(p + 1);
        fileName = arg.outputFolder + "/" + fileName;
        for (size_t i = 0; i < barcode.size(); ++i) {
            if (!barcodeWrite[i]) continue;
            outputNames[i] = fileName + "_" + sampleid[i] + ".fastq.gz";
        }
        outputNames[barcode.size()] = fileName + "_unsigned.fastq.gz";
        for (size_t i = 0; i <= barcode.size(); ++i) {
            if (!barcodeWrite[i]) continue;
            cerr << i << "\t"
                    << (i < barcode.size() ? barcode[i] : "undecoded") << "\t"
                    << outputNames[i] << endl;
            outputfp[i] = gzopen(outputNames[i].c_str(), "wb");
            if (outputfp[i] == Z_NULL) {
                cerr << "Can't open " << outputNames[i] << endl;
                exit(0);
            }
        }

        int l;
        fp = gzopen(inputReadFile[iFile].c_str(), "r");
        if (fp == Z_NULL) {
            cerr << "Can't open " << inputReadFile[iFile] << endl;
            exit(1);
        }
        seq = kseq_init(fp);
        size_t icount = 0;
        while ((l = kseq_read(seq)) >= 0) {
            //string seqname = string(seq->name.s);

            int ibc = readBarcodeIdx[icount];
            if (ibc < 0) ibc = barcode.size();

            tmps = "@" + string(seq->name.s);
            if (seq->comment.l) tmps += " " + string(seq->comment.s);
            tmps += "\n";
            tmps += string(seq->seq.s) + "\n+\n";
            if (seq->qual.l) tmps += string(seq->qual.s) + "\n";
            gzwrite(outputfp[ibc], &tmps[0], tmps.size());
            ++readsWrittenPerBarcode[ ibc ];

            ++icount;
            if (icount % 1000000 == 0) cerr << icount << endl;
        }

        gzclose(fp); //close input
        for (size_t i = 0; i <= barcode.size(); ++i) gzclose(outputfp[i]); //close output
        for (size_t i = 0; i < barcode.size() + 1; ++i)
            if (readsWrittenPerBarcode[i] != readsPerBarcode[i]) {
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

