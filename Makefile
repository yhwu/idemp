idemp: idemp.o functions.o
	g++ -o idemp idemp.o functions.o -lz

idemp.o: idemp.cpp
	g++ -c -O3 idemp.cpp
functions.o: functions.cpp
	g++ -c -O3 functions.cpp

clean:
	rm *.o idemp

test:   idemp barcode_sample.txt I1.fastq R1.fastq R2.fastq
	./idemp -b barcode_sample.txt -I1 I1.fastq -R1 R1.fastq -R2 R2.fastq
