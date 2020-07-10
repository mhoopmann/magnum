/*
Copyright 2018, Michael R. Hoopmann, Institute for Systems Biology

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef _MANALYSIS_H
#define _MANALYSIS_H

#include "MDB.h"
#include "MData.h"
#include "MIons.h"
#include "Threading.h"
#include "ThreadPool.h"
#include "CometDecoys.h"

//=============================
// Structures for threading
//=============================
struct mAnalysisStruct {
  bool*       bKIonsMem;    //Pointer to the memory manager array to mark memory is in use
  Mutex*      mutex;        //Pointer to a mutex for protecting memory
  mPeptide*   pep;
  int         pepIndex;
  mAnalysisStruct(Mutex* m, mPeptide* p, int i){
    mutex=m;
    pep=p;
    pepIndex=i;
  }
  ~mAnalysisStruct(){
    //Mark that memory is not being used, but do not delete it here.
    Threading::LockMutex(*mutex);
    if(bKIonsMem!=NULL) *bKIonsMem=false;
    bKIonsMem=NULL;
    Threading::UnlockMutex(*mutex);
    mutex=NULL;   //release mutex
    pep=NULL;
  }
};

class MAnalysis{
public:

  //Constructors & Destructors
  MAnalysis  (mParams& p, MDatabase* d, MData* dat);
  ~MAnalysis ();

  //Master Functions
  bool doPeptideAnalysis   ();
  bool doEValueAnalysis();
  bool doEValuePrecalc();

private:

  //Thread-start functions
  static void analyzePeptideProc(mAnalysisStruct* s); 
  static void analyzeEValueProc(MSpectrum* s);
  static void analyzeEValuePrecalcProc(MSpectrum* s);

  //Analysis functions
  static bool analyzePeptide(mPeptide* p, int pepIndex, int iIndex);

  //Private Functions
  bool         allocateMemory          (int threads);
  static bool  analyzeSinglets         (mPeptide& pep, int index, int iIndex);
  static void  checkXLMotif            (int motifA, char* motifB, vector<int>& v);
  void         deallocateMemory        (int threads);
  static void  scoreSingletSpectra     (int index, int sIndex, double mass, int len, int pep, char k, double minMass, double maxMass, int iIndex);
  static void  scoreSpectra            (vector<int>& index, int sIndex, int len, double modMass, int pep1, int pep2, int k1, int k2, int link, int iIndex);
  static float magnumScoring           (int specIndex, double modMass, int sIndex, int iIndex, int z=0);
  static void  setBinList              (kMatchSet* m, int iIndex, int charge, double preMass, mPepMod* mods, char modLen);

  //Data Members
  static bool*      bKIonsManager;
  static MDatabase* db;
  static MIons*     ions;
  static double     maxMass;
  static double     minMass;
  static mParams    params;
  static MData*     spec;
  static bool*      adductSites;
  static bool**     scanBuffer;

  static int        numIonSeries;

  static int* pepMassSize;
  static double**   pepMass;
  static bool** pepBin;
  static int* pepBinSize;
  static int skipCount;
  static int nonSkipCount;

  static MDecoys decoys;

  static bool scoreSingletSpectra2(int index, int sIndex, double mass, double xlMass, int counterMotif, int len, int pep, char k, double minMass, int iIndex, char linkSite, int linkIndex);

  static Mutex  mutexKIonsManager; 
  static Mutex* mutexSpecScore; //these signal PSM list reads/additions/deleteions
  static Mutex** mutexSingletScore; //these signal singlet list reads/additions/deletions

  //Utilities
  static int compareD           (const void *p1,const void *p2);
  static int comparePeptideBMass(const void *p1,const void *p2);
  
};

#endif
