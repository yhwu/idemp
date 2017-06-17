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

#endif /* IDEMP_H */

