test: H_Border.o h_Coord.o Inside.o Partition.o RWsgy.o DataTran.o testTDFWI.o testmain.o
	mpic++ -o test H_Border.o h_Coord.o Inside.o Partition.o  RWsgy.o DataTran.o testTDFWI.o testmain.o

H_Border.o: Partition.h H_Border.cpp
	mpic++ -c -std=c++11 H_Border.cpp

h_Coord.o: Partition.h h_Coord.cpp
	mpic++ -c -std=c++11 h_Coord.cpp

Inside.o: Partition.h Inside.cpp
	mpic++ -c -std=c++11 Inside.cpp

Partition.o: Partition.cpp Partition.h H_Border.cpp h_Coord.cpp Inside.cpp
	mpic++ -c -std=c++11 Partition.cpp

RWsgy.o: RWsgy.cpp RWsgy.h Partition.h
	mpic++ -c -std=c++11 RWsgy.cpp

DataTran.o: DataTran.cpp DataTran.h Partition.h 
	mpic++ -c -std=c++11 DataTran.cpp

testTDFWI.o: testTDFWI.cpp testTDFWI.h Partition.h
	mpic++ -c -std=c++11 testTDFWI.cpp

testmain.o: testmain.cpp 
	mpic++ -c -std=c++11 testmain.cpp
