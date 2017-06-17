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
#include <cstdlib>
#include <stdexcept>
using namespace std;

/**** user headers ****/
#include "idemp_func.h"
#include "functions.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

/* check that R1 and R2 have the same read names, read in I1 reads */
void check_are_read_names_same(string I1File, string R1File, string R2File,
        string& readBarcode) {
    readBarcode = "";

    vector<string> inpFile(0);
    inpFile.push_back(I1File);
    inpFile.push_back(R1File);
    if (R2File != "") inpFile.push_back(R2File);

    vector<gzFile> inpFH(inpFile.size(), Z_NULL);
    vector<kseq_t *> seq(inpFile.size());
    for (size_t i = 0; i < inpFile.size(); ++i) {
        inpFH[i] = gzopen(inpFile[i].c_str(), "r");
        if (inpFH[i] == Z_NULL) {
            cerr << "Cannot open " << inpFile[i] << endl;
            exit(1);
        }
        seq[i] = kseq_init(inpFH[i]);
    }

    bool isReadNameSame = true;
    bool readDone = false;
    vector<int> l(inpFile.size());
    vector<string> sname(inpFile.size());
    size_t icount = 0;
    while (!readDone) {
        readDone = true;
        for (size_t i = 0; i < inpFile.size(); ++i) {
            l[i] = kseq_read(seq[i]);
            sname[i] = string(seq[i]->name.s);
            if (l[i] >= 0) readDone = false;
        }
        if (readDone) break; // done reading all reads 
        ++icount;
        readBarcode += string(seq[0]->seq.s) + "\t";
        for (size_t i = 1; i < inpFile.size(); ++i) {
            if (sname[i] == sname[i - 1]) continue;
            isReadNameSame = false;
            cerr << "Read names are different at read " << icount << " "
                    << "between files\n" << inpFile[i - 1] << " " << inpFile[i] << endl
                    << sname[i - 1] << "\n" << sname[i] << endl;
            break;
        }
    }

    if ( ! isReadNameSame ) {
        cerr << "Read names are not same, exit" << endl;
        exit(1);
    }
    cerr << "Read names are same:\t" << isReadNameSame << endl;

    for (size_t i = 0; i < inpFile.size(); ++i) {
        gzclose(inpFH[i]);
        kseq_destroy(seq[i]);
    }
}


/* parse input */
IdempArgument parse_arguments(int argc, char* argv[]) {
    std::unordered_map<std::string, std::string> args;
    std::unordered_map<std::string, std::string>::iterator got;
    vector<string> ARGV(0);
    for (int i = 0; i < argc; ++i) ARGV.push_back(string(argv[i]));
    for (int i = 1; i < (int) ARGV.size(); ++i) {
        if (ARGV[i] == "-b") { args["barcodeFile"] = ARGV[i + 1]; continue;}
        if (ARGV[i] == "-I1") { args["I1File"] = ARGV[i + 1]; continue;}
        if (ARGV[i] == "-R1") { args["R1File"] = ARGV[i + 1]; continue;}
        if (ARGV[i] == "-R2") { args["R2File"] = ARGV[i + 1]; continue; }
        if (ARGV[i] == "-m") { args["nMismatch"] = ARGV[i + 1]; continue;}
        if (ARGV[i] == "-o") { args["outputFolder"] = ARGV[i + 1]; continue;}
    }
    
    IdempArgument arg;
    
    got = args.find("barcodeFile");
    std::string barcodeFile = got == args.end() ? "" : args["barcodeFile"];

    if (args.find("barcodeFile") == args.end()) throw std::invalid_argument("need barcodeFile, -b barcodeFile");
    if (args.find("I1File") == args.end()) throw std::invalid_argument( "need I1File, -I1 I1File" );
    if (args.find("R1File") == args.end()) throw std::invalid_argument( "need R1File, -R1 R1File" );
    if (args.find("R2File") == args.end()) throw std::invalid_argument( "need R2File, -R2 R1File" );

    if (args.find("nMismatch") != args.end()) arg.nMismatch = std::stoi(args["nMismatch"]);
    if (args.find("outputFolder") != args.end() ) arg.outputFolder = args["outputFolder"];
    arg.barcodeFile = args["barcodeFile"]; 
    arg.I1File = args["I1File"];
    arg.R1File = args["R1File"];
    arg.R2File = args["R2File"];
        
    return arg;
}
