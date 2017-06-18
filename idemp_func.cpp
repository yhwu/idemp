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
//using namespace std;

/**** user headers ****/
#include "idemp_func.h"
#include "functions.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


void map_readbarcode(string & readBarcode, vector<size_t> & readBarcode_qpos, vector<string> & barcode,
        int minEditDistance,
        vector<int> & readBarcodeIdx, vector<int> & readBarcodeMis) {

    std::unordered_map<std::string, int> barcodeDict;
    std::unordered_map<std::string, int>::iterator got;
    for (size_t j = 0; j < barcode.size(); ++j) {
        got = barcodeDict.find(barcode[j]);
        if (got != barcodeDict.end()) throw::invalid_argument("dupilcated barcode: " + barcode[j]);
        barcodeDict[barcode[j]] = j;
    }


    for (size_t i = 0; i < readBarcode_qpos.size() - 1; ++i) {

        string query = readBarcode.substr(readBarcode_qpos[i],
                readBarcode_qpos[i + 1] - readBarcode_qpos[i] - 1);
        readBarcodeIdx[i] = barcode.size();
        readBarcodeMis[i] = barcode[0].length();

        if ((i + 1) % 1000000 == 0) cerr << i + 1 << endl;
        if (query.size() != barcode[0].size()) {
            cerr << "Warning " << query << " does not have same length with barcodes" << endl;
        };


        // check exact match with barcode dictionary; most read barcodes should match 
        bool ismatched = false;
        if (barcodeDict.find(query) != barcodeDict.end()) {
            readBarcodeIdx[i] = barcodeDict[query];
            readBarcodeMis[i] = 0;
            ismatched = true;
        }
        if (ismatched) continue;

        int minMismatch = barcode[0].size();

        // check single base mutation
        for (size_t j = 0; j < barcode.size(); ++j) {
            int mismatch = 0;
            for (size_t k = 0; k < query.size(); ++k) {
                mismatch += (query[k] != barcode[j][k]);
                if (mismatch >= minEditDistance) {
                    mismatch = query.size();
                    break;
                }
            }
            if (mismatch < minMismatch) {
                minMismatch = mismatch;
                readBarcodeIdx[i] = j;
                readBarcodeMis[i] = minMismatch;
                if (minMismatch <= 1) break;
            }
        }
        if (minMismatch <= 1) continue; // next best

        // check edit distance for indel
        for (size_t j = 0; j < barcode.size(); ++j) {
            int mismatch = edit_distance(query, barcode[j]);
            if (mismatch < minMismatch) {
                minMismatch = mismatch;
                readBarcodeIdx[i] = j;
                readBarcodeMis[i] = minMismatch;
            }
        }
    }

}

int min_edit_distance(vector<string> & barcode, vector<string> & sampleid) {
    int minEditDistance = barcode[0].length();

    vector<string> miniEdbarcodes(0);
    cerr << "Pairwise barcode edit distance:" << endl;
    for (size_t i = 0; i < barcode.size(); ++i) {
        for (size_t j = 0; j <= i; ++j) {
            int ed = edit_distance(barcode[i], barcode[j]);
            if (j == 0) cerr << ed;
            else cerr << " " << ed;
            if (minEditDistance >= ed && i > j) {
                if (minEditDistance > ed) miniEdbarcodes.clear();
                minEditDistance = ed;
                miniEdbarcodes.push_back(barcode[i] + ":" + barcode[j]);
            }
            if ((ed == 0) &&
                    (j < i) &&
                    (sampleid[i] != sampleid[j])) {
                cerr << "one barcode points to two sampleids:\n"
                        << barcode[i] << "\t" << sampleid[i] << "\t" << sampleid[j]
                        << endl;
                throw std::invalid_argument("duplicated sample ids");
            }
        }
        cerr << endl;
    }

    cerr << "\nClosest barcodes, editDistance=" << minEditDistance << endl;
    for (size_t i = 0; i < miniEdbarcodes.size(); ++i) cerr << miniEdbarcodes[i] << endl;
    cerr << endl;

    /* check if any sampleids are same */
    for (size_t i = 0; i < sampleid.size(); ++i) {
        for (size_t j = i + 1; j < sampleid.size(); ++j) {
            if (sampleid[i] == sampleid[j]) {
                if (barcode[i] != barcode[j]) {
                    cerr << "two barcodes point to the same sampleid:\n"
                            << barcode[i] << "\t" << barcode[j] << "\t" << sampleid[i]
                            << endl;
                    throw std::invalid_argument("duplicated barcodes for one sampleid");
                }
            }
        }
    }

    return minEditDistance;
}

void read_barcode_sampleid(string barcodeFile,
        vector<string> & barcode, vector<string> & sampleid) {

    ifstream FIN(barcodeFile.c_str());
    if (!FIN) {
        cerr << "Cannot open " << barcodeFile << endl;
        exit(1);
    }
    string tmps1, tmps2;
    while (FIN >> tmps1 >> tmps2) {
        if (tmps1[0] == '#') continue;
        if (tmps1.find("arcode") != string::npos) continue;
        barcode.push_back(tmps1);
        sampleid.push_back(tmps2);
    }
    FIN.close();

    if (barcode.empty()) throw std::invalid_argument("no barcode found");

}

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

    if (!isReadNameSame) {
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
        if (ARGV[i] == "-b") {
            args["barcodeFile"] = ARGV[i + 1];
            continue;
        }
        if (ARGV[i] == "-I1") {
            args["I1File"] = ARGV[i + 1];
            continue;
        }
        if (ARGV[i] == "-R1") {
            args["R1File"] = ARGV[i + 1];
            continue;
        }
        if (ARGV[i] == "-R2") {
            args["R2File"] = ARGV[i + 1];
            continue;
        }
        if (ARGV[i] == "-m") {
            args["nMismatch"] = ARGV[i + 1];
            continue;
        }
        if (ARGV[i] == "-o") {
            args["outputFolder"] = ARGV[i + 1];
            continue;
        }
    }

    IdempArgument arg;

    got = args.find("barcodeFile");
    std::string barcodeFile = got == args.end() ? "" : args["barcodeFile"];

    if (args.find("barcodeFile") == args.end()) throw std::invalid_argument("need barcodeFile, -b barcodeFile");
    if (args.find("I1File") == args.end()) throw std::invalid_argument("need I1File, -I1 I1File");
    if (args.find("R1File") == args.end()) throw std::invalid_argument("need R1File, -R1 R1File");
    if (args.find("R2File") == args.end()) throw std::invalid_argument("need R2File, -R2 R1File");

    if (args.find("nMismatch") != args.end()) arg.nMismatch = std::stoi(args["nMismatch"]);
    if (args.find("outputFolder") != args.end()) arg.outputFolder = args["outputFolder"];
    arg.barcodeFile = args["barcodeFile"];
    arg.I1File = args["I1File"];
    arg.R1File = args["R1File"];
    arg.R2File = args["R2File"];

    return arg;
}
