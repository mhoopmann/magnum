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
#include <iostream>
#include "CometDecoys.h"
#include "Threading.h"
#include "ThreadPool.h"

using namespace MSToolkit;

class MData {
public:

  MData();
  MData(mParams* p);
  ~MData();

  MSpectrum& operator[ ](const int& i);
  MSpectrum& at(const int& i);
  MSpectrum* getSpectrum(const int& i);

  bool      createDiag        (FILE*& f);
  NeoPepXMLParser* createPepXML(string& str, MDatabase& db);
  bool      createPercolator  (FILE*& f);
  bool      createTXT         (FILE*& f);
  void      diagSinglet       ();
  void      exportPepXML      (NeoPepXMLParser*& p, vector<mResults>& r);
  void      exportPercolator  (FILE*& f, vector<mResults>& r);
  void      exportTXT         (FILE*& f, vector<mResults>& r);
  bool      getBoundaries     (double mass1, double mass2, vector<int>& index, bool* buffer);
  bool      getBoundaries2    (double mass, double prec, vector<int>& index, bool* buffer);
  double    getMaxMass        ();
  double    getMinMass        ();
  bool      mapPrecursors     ();
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
  void      xCorr             ();

  bool*     getAdductSites    ();

private:

  //Data Members
  bool* bScans;
  char               version[32];
  vector<MSpectrum>  spec;
  vector<mMass>      massList;
  mParams*           params;
  MParams*           parObj;
  MIons              aa;
  MLog*              mlog;
  int                pepXMLindex;

  bool adductSite[128];

  static void xCorrProc(MSpectrum* s);

  //Utilities
  void        centroid          (Spectrum& s, Spectrum& out, double resolution, int instrument=0);
  void        collapseSpectrum  (Spectrum& s);
  static int  compareInt        (const void *p1, const void *p2);
  static int  compareMassList   (const void *p1, const void *p2);
  int         getCharge         (Spectrum& s, int index, int next);
  double      polynomialBestFit (vector<double>& x, vector<double>& y, vector<double>& coeff, int degree=2);
  bool        processPath       (const char* in_path, char* out_path);
  string      processPeptide(mPeptide& pep, vector<mPepMod>& mod, int site, double massA, MDatabase& db);

};


#endif

