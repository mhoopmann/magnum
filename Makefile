#Set these variables if needed
C = gcc
CC = g++
FLAGS = -O3 -static -std=c++11 -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC

#Paths to supporting software
MSTOOLKITPATH = ../MSToolkit
HARDKLORPATH = ../Hardklor
PEPXMLPATH = ../NeoPepXMLParser

#Do not touch these variables
LIBPATH = -L$(MSTOOLKITPATH) -L$(HARDKLORPATH) -L$(PEPXMLPATH)
LIBS = -lmstoolkitlite -lhardklor -lneopepxmlparser -lpthread
INCLUDE = -I$(MSTOOLKITPATH)/include -I$(HARDKLORPATH) -I$(PEPXMLPATH)


#Do not touch these variables
MAGNUM = MagnumManager.o MParams.o MAnalysis.o MData.o MDB.o MLog.o MPrecursor.o MSpectrum.o MIons.o MIonSet.o MTopPeps.o Threading.o CometDecoys.o


#Make statements
magnum : Magnum.cpp $(MAGNUM)
	$(CC) $(FLAGS) $(INCLUDE) $(MAGNUM) Magnum.cpp $(LIBPATH) $(LIBS) -o magnum

clean:
	rm *.o magnum


#Magnum objects
MParams.o : MParams.cpp
	$(CC) $(FLAGS) $(INCLUDE) MParams.cpp -c

MAnalysis.o : MAnalysis.cpp
	$(CC) $(FLAGS) $(INCLUDE) MAnalysis.cpp -c

MData.o : MData.cpp
	$(CC) $(FLAGS) $(INCLUDE) MData.cpp -c

MDB.o : MDB.cpp
	$(CC) $(FLAGS) $(INCLUDE) MDB.cpp -c

MLog.o : MLog.cpp
	$(CC) $(FLAGS) $(INCLUDE) MLog.cpp -c

MagnumManager.o : MagnumManager.cpp
	$(CC) $(FLAGS) $(INCLUDE) MagnumManager.cpp -c

MPrecursor.o : MPrecursor.cpp
	$(CC) $(FLAGS) $(INCLUDE) MPrecursor.cpp -c

MSpectrum.o : MSpectrum.cpp
	$(CC) $(FLAGS) $(INCLUDE) MSpectrum.cpp -c

MIons.o : MIons.cpp
	$(CC) $(FLAGS) $(INCLUDE) MIons.cpp -c

MIonSet.o : MIonSet.cpp
	$(CC) $(FLAGS) $(INCLUDE) MIonSet.cpp -c

MTopPeps.o : MTopPeps.cpp
	$(CC) $(FLAGS) $(INCLUDE) MTopPeps.cpp -c

Threading.o : Threading.cpp
	$(CC) $(FLAGS) $(INCLUDE) Threading.cpp -c

CometDecoys.o : CometDecoys.cpp
	$(CC) $(FLAGS) $(INCLUDE) CometDecoys.cpp -c
