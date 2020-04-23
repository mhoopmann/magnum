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

#ifndef _MPRECURSOR_H
#define _MPRECURSOR_H

#include "MStructs.h"
#include "MSpectrum.h"

#include "CHardklor2.h"
#include "CModelLibrary.h"
#include "CHardklorSetting.h"
#include "CHardklorVariant.h"
#include "MSObject.h"
#include "MSReader.h"
#include <cmath>
#include <deque>

#define GAUSSCONST 5.5451774444795623

typedef struct mRTProfile{
  int   scan;
  float rt;
} mRTProfile;

typedef struct mScanBin{
  int   index;
  float intensity;
} mScanBin;

class MPrecursor {
public:

  MPrecursor(mParams* p);
  ~MPrecursor();

  bool  estimatePrecursor (MSpectrum& s);
  int   getSpecRange      (MSpectrum& pls);
  bool  setFile           ();

  std::vector<int> preCharges;

private:

  //Functions
  void    averageScansCentroid(std::vector<MSToolkit::Spectrum*>& s, MSToolkit::Spectrum& avg, double min, double max);
  void    centroid(MSToolkit::Spectrum& s, MSToolkit::Spectrum& out, double resolution, int instrument = 0);
  int     findPeak(MSToolkit::Spectrum& s, double mass);
  int     findPeak(MSToolkit::Spectrum& s, double mass, double prec);
  double  polynomialBestFit     (std::vector<double>& x, std::vector<double>& y, std::vector<double>& coeff, int degree=2);

  //Data Members
  char      fileName[256];
  int       lastScan;
  MSToolkit::MSObject  mso;
  MSToolkit::MSReader  msr;
  MSToolkit::MSReader  msw;
  MSToolkit::Spectrum  avg;
  MSToolkit::Spectrum  cent;
  MSToolkit::Spectrum  spec;
  
  CAveragine*      averagine;
	CMercury8*       mercury;
  CModelLibrary*   models;
  CHardklor2*      h;
  CHardklorSetting hs;
  CHardklorVariant hv;

  std::deque<MSToolkit::Spectrum>*  centBuf;

  mParams*   params;

  //Utilities
  static int    comparePLPScanNum (const void *p1, const void *p2);
  static int    compareScanBinRev (const void *p1, const void *p2);

};

#endif
