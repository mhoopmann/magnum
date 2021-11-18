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

#ifndef _MSPECTRUM_H
#define _MSPECTRUM_H

#include <cmath>
#include <list>
#include <vector>
#include "MHistogram.h"
#include "MStructs.h"
#include "MTopPeps.h"
#include "CometDecoys.h"

#define HISTOSZ 152

typedef struct sHistoPep {
  int pepIndex;
  int topScore;
} sHistoPep;

class MSpectrum {

public:

  //Constructors & Destructors
  MSpectrum(const int& i, const double& bs, const double& os, const int& th);
  MSpectrum(mParams& p);
  MSpectrum(const MSpectrum& p);
  ~MSpectrum();

  //Operators
  MSpectrum&  operator=(const MSpectrum& p);
  mSpecPoint&  operator[](const int& i);

  //Data Members
  mSparseMatrix*  xCorrSparseArray;
  int             xCorrSparseArraySize;
  char**          kojakSparseArray;
  int             kojakBins;
  
  double lowScore;

  int cc;
  int sc;

  //Accessors
  double              getBinOffset          ();
  int                 getCharge             ();
  bool                getInstrumentPrecursor();
  double              getInvBinSize         ();
  float               getMaxIntensity       ();
  double              getMZ                 ();
  mPrecursor&         getPrecursor          (int i);
  mPrecursor*         getPrecursor2         (int i);
  float               getRTime              ();
  int                 getScanNumber         ();
  mScoreCard&         getScoreCard          (int i);
  int                 getSingletCount       ();
  mScoreCard&         getSingletScoreCard   (int i);
  MTopPeps*           getTopPeps            (int index);
  int                 size                  ();
  int                 sizePrecursor         ();
  
  mScoreCard*    singletFirst;   //pointer to start of linked list
  mScoreCard*    singletLast;    //pointer to end of linked list
  int            singletMax;

  int hpSize;
  sHistoPep* hp;
  int histogram[HISTOSZ];
  int histogramCount;
  int histoMaxIndex;

  //** temporary
  //int hX[60][HISTOSZ];
  //int hXCount[60];
  //void tHistogram(double score, int len);
  //void exportHisto();
  //**

  MHistogram** mHisto;
  MDecoys* decoys;
  double computeE(double score, int len);
  bool ionSeries[6];
  int maxHistogramCount;
  double minAdductMass;
  int maxPepLen;
  double bigMonoMass;
  int bigZ;

  //Modifiers
  void addPoint               (mSpecPoint& s);
  void addPrecursor           (mPrecursor& p, int sz);
  void clear                  ();
  void erasePrecursor         (int i);
  void setCharge              (int i);
  void setInstrumentPrecursor (bool b);
  void setMaxIntensity        (float f);
  void setMZ                  (double d);
  void setPrecursor           (double d, int i);
  void setRTime               (float f);
  void setScanNumber          (int i);

  //Functions
  bool calcEValue(mParams* params, MDecoys& decoys);
  bool checkReporterIon(double mz,mParams* params);
  void  checkScore(mScoreCard& s, int th);
  void  checkSingletScore (mScoreCard& s);
  bool generateXcorrDecoys(mParams* params, MDecoys& decoys);
  bool generateXcorrDecoys2(int maxPepLen);
  bool generateXcorrDecoys3(int minP, int maxP, int depth);
  void linearRegression(double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr);
  void linearRegression2(double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared);
  void linearRegression3(double& slope, double& intercept, double& rSquared);
  void linearRegression4(int* h, int sz, double& slope, double& intercept, double& rSquared);
  double makeXCorrB(int decoyIndex, double modMass, int maxZ, int len, int offset=0);
  double makeXCorrY(int decoyIndex, double modMass, int maxZ, int len, int offset=0);
  void  resetSingletList  ();
  void  shortResults(std::vector<mScoreCard2>& v);
  void  shortResults2(std::vector<mScoreCard3>& v);
  void  sortMZ            ();
  void  xCorrScore        ();

private:

  //Data members
  double                binOffset;
  double                binSize;
  int                   charge;
  bool                  instrumentPrecursor;
  double                invBinSize;
  float                 maxIntensity;
  double                mz;
  std::vector<mPrecursor>*   precursor;
  float                 rTime;
  int                   scanNumber;
  int                   singletCount;
  
  
  std::vector<MTopPeps>*     singlets;
  std::vector<mSpecPoint>*   spec;
  mScoreCard            topHit[20];
  int                   xCorrArraySize;
  

  //Functions
  void BinIons      (mPreprocessStruct *pPre);
  void MakeCorrData (double *pdTempRawData, mPreprocessStruct *pPre, double scale);

  void kojakXCorr   ();

  //Utilities
  static int compareIntensity (const void *p1,const void *p2);
  static int compareMZ        (const void *p1,const void *p2);



};

#endif
