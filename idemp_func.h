#ifndef IDEMP_H
#define IDEMP_H

#include <exception>
#include <sstream>
using namespace std;
#include <vector>
#include <stdio.h>
#include <inttypes.h> 

class IdempArgument {
  public:
    string I1File = "", 
            R1File = "", 
            R2File = "", 
            barcodeFile = "", 
            outputFolder = "decoded";
    int nMismatch=1;
};

IdempArgument parse_arguments(int argc, char* argv[]);

void read_barcode_sampleid(string barcodeFile, 
        vector<string> & barcode, vector<string> & sampleid);

void check_are_read_names_same(string I1File, string R1File, string R2File,
        string& readBarcode);

#endif /* IDEMP_H */

