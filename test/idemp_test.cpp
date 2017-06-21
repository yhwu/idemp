#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include "../idemp_func.h"
#include "../functions.h"
/*
 * Simple C++ Test Suite
 */

void test_parse_arguments() {
    int argc = 9;
    const char* myargv[] = {"a.out", "-b", "barcodeFile", "-I1", "I1File", "-R1", "R1File", "-R2", "R2File"};
    char ** argv = const_cast<char**>(myargv);
    IdempArgument arg = parse_arguments(argc, argv);
    assert ( (arg.I1File == "I1File") && (arg.R1File == "R1File")  && (arg.R2File == "R2File") && 
            (arg.barcodeFile == "barcodeFile") && (arg.nMismatch==1));
}


void test_edit_distance() {
    assert ( edit_distance("GACTTTCCCTCG", "GACTTTCCCTCG")==0 );
    assert ( edit_distance("GACTTTCCCTCG", "GACCTTCCCTCG")==1 );
    assert ( edit_distance("GACTTTCCCTCG", "GACATTTCCCTC")==2 );
    assert ( edit_distance("GACTTTCCCTCG", "GACTTCCCTCGG")==2 );
    assert ( edit_distance("NNNNNNNNNNNN", "GACTTCCCTCGG")==12 );
}


int main(int argc, char** argv) {
    std::cout << "test_edit_distance (idemp_test)" << std::endl;
    test_edit_distance();
    std::cout << "test_edit_distance (idemp_test) pass" << std::endl;

    std::cout << "test_parse_arguments (idemp_test)" << std::endl;
    test_parse_arguments();
    std::cout << "test_parse_arguments (idemp_test) pass" << std::endl;

    return (EXIT_SUCCESS);
}

