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

#ifndef _MDATA_H
#define _MDATA_H

#include "MDB.h"
#include "MIons.h"
#include "MLog.h"
#include "MParams.h"
#include "MPrecursor.h"
#include "MSpectrum.h"
#include "MSReader.h"
#include "NeoPepXMLParser.h"
#include <deque>
#include <iostream>
#include "CometDecoys.h"
#include "Threading.h"
#include "ThreadPool.h"

//=============================
// Structures for threading
//=============================
typedef struct mMS2struct{
  MSToolkit::Spectrum* s;
  MSpectrum* pls;
  int state;
  bool thread;
  mMS2struct(MSToolkit::Spectrum* sp, mParams* params){
    s = sp;
    pls = new MSpectrum(*params);
    state = 0;
    thread = false;
  }
  ~mMS2struct(){
    delete s;
    pls = NULL;  //this needs to be deleted elsewhere
  }
} mMS2struct;

class MData {
public:

  MData();
  MData(mParams* p);
  ~MData();

  MSpectrum& operator[ ](const int& i);
  MSpectrum& at(const int& i);
  MSpectrum* getSpectrum(const int& i);

  bool      createDiag        (FILE*& f);
  NeoPepXMLParser* createPepXML(std::string& str, MDatabase& db);
  bool      createPercolator  (FILE*& f, FILE*& f2);
  bool      createTXT         (FILE*& f);
  void      diagSinglet       ();
  void      exportPepXML      (NeoPepXMLParser*& p, std::vector<mResults>& r);
  void      exportPercolator  (FILE*& f, std::vector<mResults>& r);
  void      exportTXT         (FILE*& f, std::vector<mResults>& r);
  bool      getBoundaries     (double mass1, double mass2, std::vector<int>& index, bool* buffer);
  bool      getBoundaries2    (double mass, double prec, std::vector<int>& index, bool* buffer);
  double    getMaxMass        ();
  double    getMinMass        ();
  void      outputDiagnostics (FILE* f, MSpectrum& s, MDatabase& db);
  bool      outputResults     (MDatabase& db);
  void      processPSM        (MSpectrum& s, mScoreCard3& sc, mResults& r);
  void      processSpectrumInfo (MSpectrum& s, mResults& r);
  bool      readSpectra       ();
  void      setAdductSites    (std::string s);
  void      setLog            (MLog* c);
  void      setParams         (MParams* p);
  void      setVersion        (const char* v);
  int       size              ();

  bool*     getAdductSites    ();

private:

  //Data Members
  bool* bScans;
  char               version[32];
  std::vector<MSpectrum*>  spec;
  std::vector<mMass>      massList;
  static mParams*           params;
  MParams*           parObj;
  MIons              aa;
  static MLog*              mlog;
  int                pepXMLindex;

  //Common memory to be shared by all threads during spectral processing
  static bool* memoryPool;
  static double** tempRawData;
  static double** tmpFastXcorrData;
  static float**  fastXcorrData;
  static Mutex    mutexMemoryPool;
  static mPreprocessStruct** preProcess;

  //Optimized file loading structures for spectral processing
  static std::deque<MSToolkit::Spectrum*> dMS1;
  static std::vector<MSToolkit::Spectrum*> vMS1Buffer;
  static Mutex mutexLockMS1;
  static CHardklor2** h;
  static CHardklorSetting hs;
  static Mutex* mutexHardklor;
  static CAveragine** averagine;
  static CMercury8** mercury;
  static CModelLibrary* models;
  static bool* bHardklor;
  static int maxPrecursorMass;

  bool adductSite[128];

  //static void xCorrProc(MSpectrum* s);

  //spectral processing functions
  static void averageScansCentroid(std::vector<MSToolkit::Spectrum*>& s, MSToolkit::Spectrum& avg, double min, double max);
  static int  findPeak(MSToolkit::Spectrum* s, double mass);
  static int  findPeak(MSToolkit::Spectrum* s, double mass, double prec);
  static void formatMS2(MSToolkit::Spectrum* s, MSpectrum* pls);
  void initHardklor();
  void memoryAllocate();
  void memoryFree();
  static void processMS2(mMS2struct* s);
  static int  processPrecursor(mMS2struct* s, int tIndex);
  void releaseHardklor();

  //Utilities
  static void        centroid(MSToolkit::Spectrum* s, MSpectrum* out, double resolution, int instrument = 0);
  static void        collapseSpectrum(MSpectrum& s);
  static int  compareInt        (const void *p1, const void *p2);
  static int  compareMassList   (const void *p1, const void *p2);
  static int compareScanBinRev2(const void *p1, const void *p2);
  static bool compareSpecPoint(const mSpecPoint& p1, const mSpecPoint& p2){ return p1.mass<p2.mass; }
  static int         getCharge(MSpectrum& s, int index, int next);
  static double      polynomialBestFit (std::vector<double>& x, std::vector<double>& y, std::vector<double>& coeff, int degree=2);
  bool        processPath       (const char* in_path, char* out_path);
  std::string      processPeptide(mPeptide& pep, std::vector<mPepMod>& mod, int site, double massA, MDatabase& db);

};


#endif

