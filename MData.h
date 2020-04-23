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
#include "MParams.h"
#include "MPrecursor.h"
#include "MSpectrum.h"
#include "MSReader.h"
#include "pepXMLWriter.h"
#include <iostream>
#include "CometDecoys.h"
#include "Threading.h"
#include "ThreadPool.h"

using namespace MSToolkit;

/*
#ifdef _MSC_VER
#include <direct.h>
#define getcwd _getcwd
#define slashdir '\\'
#else
#define slashdir '/'
#endif
*/

class MData {
public:

  MData();
  MData(mParams* p);
  ~MData();

  MSpectrum& operator[ ](const int& i);
  MSpectrum& at(const int& i);
  MSpectrum* getSpectrum(const int& i);

  void      diagSinglet       ();
  bool      getBoundaries     (double mass1, double mass2, vector<int>& index, bool* buffer);
  bool      getBoundaries2    (double mass, double prec, vector<int>& index, bool* buffer);
  double    getMaxMass        ();
  double    getMinMass        ();
  bool      mapPrecursors     ();
  void      outputDiagnostics (FILE* f, MSpectrum& s, MDatabase& db);
  int       outputPepXML      (PXWSpectrumQuery& sq, MDatabase& db, kResults& r);
  void      outputPepXML2     (PXWSpectrumQuery& sq, int shIndex, kResults& r);
  bool      outputPercolator  (FILE* f, MDatabase& db, kResults& r, int count);
  bool      outputResults     (MDatabase& db, MParams& par);
  bool      readSpectra       ();
  void      setAdductSites    (char* s);
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
  MIons              aa;

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
  string      processPeptide(mPeptide& pep, vector<mPepMod>* mod, int site, double massA, MDatabase& db);

};


#endif

