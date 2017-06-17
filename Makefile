CXXFLAGS = -std=gnu++11 -O3

DEBUG ?= 1
ifeq ($(DEBUG), 1)
    CXXFLAGS += -g -DDEBUG
    CFLAGS += -g -DDEBUG
else
    CFLAGS = -DNDEBUG
endif


idemp: idemp.o functions.o
	$(CXX) $(CXXFLAGS) $? -o $@ -lz
	
idemp.o: idemp.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

idemp_func.o: idemp_func.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

functions.o: functions.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

test/idemp_test.o: test/idemp_test.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

clean:
	rm -f *.o idemp
	rm -f test/*.o test/idemp_test

test:   idemp barcode_sample.txt I1.fastq R1.fastq R2.fastq
	./idemp -b barcode_sample.txt -I1 I1.fastq -R1 R1.fastq -R2 R2.fastq -o test

unittest:   test/idemp_test.o idemp_func.o functions.o
	$(CXX) $? -o ./test/idemp_test -lz $(CXXFLAGS)
	./test/idemp_test
	