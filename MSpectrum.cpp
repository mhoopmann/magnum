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

#include "MSpectrum.h"
#include <iostream>

using namespace std;

/*============================
  Constructors & Destructors
============================*/
//MSpectrum::MSpectrum(const int& i, const double& bs, const double& os, const int& th){
//  binOffset=os;
//  binSize=bs;
//  instrumentPrecursor=false;
//  invBinSize=1.0/binSize;
//  charge = 0; 
//  maxIntensity=0;
//  mz = 0;
//  precursor = new vector<mPrecursor>;
//  singlets = new vector<MTopPeps>;
//  spec = new vector<mSpecPoint>;
//  scanNumber = 0;
//  rTime = 0;
//  xCorrArraySize=0;
//  xCorrSparseArraySize=0;
//  xCorrSparseArray=NULL;
//  bigMonoMass=0;
//  bigZ=0;
//  
//  singletCount=0;
//  singletFirst=NULL;
//  singletLast=NULL;
//  singletMax=i;
//
//  kojakSparseArray=NULL;
//  kojakBins=0;
//
//  //singletList=NULL;
//  //singletBins=0;
//
//  lowScore=0;
//
//  hpSize=th;
//  hp=new sHistoPep[hpSize];
//  for(int j=0;j<hpSize;j++){
//    hp[j].pepIndex=-1;
//    hp[j].topScore=0;
//  }
//  for(int j=0;j<HISTOSZ;j++) histogram[j]=0;
//  histogramCount=0;
//  histoMaxIndex=0;
//
//  //** temporary
//  //for(int a=0;a<60;a++){
//  //  for (int j = 0; j<HISTOSZ; j++) hX[a][j] = 0;
//  //  hXCount[a] = 0;
//  //}
//  //**
//
//  cc=0;
//  sc=0;
//
//  peakCounts=0;
//}

MSpectrum::MSpectrum(mParams& p){
  //params->topCount, params->binSize, params->binOffset, params->threads
  binOffset = p.binOffset;
  binSize = p.binSize;
  instrumentPrecursor = false;
  invBinSize = 1.0 / binSize;
  charge = 0;
  maxIntensity = 0;
  mz = 0;
  precursor = new vector<mPrecursor>;
  singlets = new vector<MTopPeps>;
  spec = new vector<mSpecPoint>;
  scanNumber = 0;
  rTime = 0;
  xCorrArraySize = 0;
  xCorrSparseArraySize = 0;
  xCorrSparseArray = NULL;
  bigMonoMass=0;
  bigZ=0;

  singletCount = 0;
  singletFirst = NULL;
  singletLast = NULL;
  singletMax = p.topCount;

  kojakSparseArray = NULL;
  kojakBins = 0;

  lowScore = 0;

  maxHistogramCount=3000;
  minAdductMass=p.minAdductMass;
  maxPepLen=p.maxPepLen;
  mHisto=new MHistogram*[maxPepLen+1];
  for(int j=0;j<p.maxPepLen+1;j++)mHisto[j]=NULL;

  for(int j=0;j<6;j++) ionSeries[j]=p.ionSeries[j];

  hpSize = p.threads;
  hp = new sHistoPep[hpSize];
  for (int j = 0; j<hpSize; j++){
    hp[j].pepIndex = -1;
    hp[j].topScore = 0;
  }
  for (int j = 0; j<HISTOSZ; j++) histogram[j] = 0;
  histogramCount = 0;
  histoMaxIndex = 0;

  //** temporary
  //for (int a = 0; a<60; a++){
  //  for (int j = 0; j<HISTOSZ; j++) hX[a][j] = 0;
  //  hXCount[a] = 0;
  //}
  //**

  cc = 0;
  sc = 0;

  peakCounts=0;
  maxX = (int)((p.maxPepMass + p.maxAdductMass + 100.0) / p.binSize);
}

MSpectrum::MSpectrum(const MSpectrum& p){
  unsigned int i;
  int j;
  spec = new vector<mSpecPoint>(*p.spec);
  for(i=0;i<20;i++) topHit[i]=p.topHit[i];
  precursor = new vector<mPrecursor>(*p.precursor);
  singlets = new vector<MTopPeps>(*p.singlets);

  binOffset = p.binOffset;
  binSize = p.binSize;
  instrumentPrecursor = p.instrumentPrecursor;
  invBinSize = p.invBinSize;
  charge= p.charge;
  maxIntensity = p.maxIntensity;
  mz = p.mz;
  scanNumber = p.scanNumber;
  rTime = p.rTime;
  xCorrArraySize = p.xCorrArraySize;
  xCorrSparseArraySize = p.xCorrSparseArraySize;
  lowScore=p.lowScore;
  for(i = 0; i<HISTOSZ; i++) histogram[i] = p.histogram[i];
  histogramCount = p.histogramCount;
  histoMaxIndex = p.histoMaxIndex;
  bigMonoMass=p.bigMonoMass;
  bigZ=p.bigZ;

  decoys=p.decoys;
  maxHistogramCount = p.maxHistogramCount;
  minAdductMass = p.minAdductMass;
  maxPepLen = p.maxPepLen;
  mHisto = new MHistogram*[maxPepLen+1];
  for (int j = 0; j<p.maxPepLen+1; j++)mHisto[j] = NULL;
  for (int j = 0; j<6; j++) ionSeries[j] = p.ionSeries[j];

  cc=p.cc;
  sc=p.sc;
  peakCounts=p.peakCounts;

  hpSize=p.hpSize;
  hp=new sHistoPep[hpSize];
  for(j=0;j<hpSize;j++)hp[j]=p.hp[j];

  singletCount=p.singletCount;
  singletMax=p.singletMax;
  singletFirst=NULL;
  singletLast=NULL;
  mScoreCard* sc=NULL;
  mScoreCard* tmp=p.singletFirst;
  if(tmp!=NULL) {
    singletFirst=new mScoreCard(*tmp);
    sc=singletFirst;
    tmp=tmp->next;
    while(tmp!=NULL){
      sc->next=new mScoreCard(*tmp);
      sc->next->prev=sc;
      sc=sc->next;
      tmp=tmp->next;
    }
    singletLast=sc;
  }

  if(p.xCorrSparseArray==NULL){
    xCorrSparseArray=NULL;
  } else {
    xCorrSparseArray = (mSparseMatrix *)calloc((size_t)xCorrSparseArraySize, (size_t)sizeof(mSparseMatrix));
    for(j=0;j<xCorrSparseArraySize;j++) xCorrSparseArray[j]=p.xCorrSparseArray[j];
  }

  kojakBins=p.kojakBins;
  if(p.kojakSparseArray==NULL){
    kojakSparseArray=NULL;
  } else {
    for(j=0;j<kojakBins;j++){
      if(p.kojakSparseArray[j]==NULL){
        kojakSparseArray[j]=NULL;
      } else {
        kojakSparseArray[j] = new char[(int)invBinSize+1];
        for(i=0;i<(unsigned int)invBinSize+1;i++) kojakSparseArray[j][i]=p.kojakSparseArray[j][i];
      }
    }
  }

  //** temporary
  //for (int a = 0; a<60; a++){
  //  for (int j = 0; j<HISTOSZ; j++) hX[a][j] = p.hX[a][j];
  //  hXCount[a] = p.hXCount[a];
  //}
  //**
  maxX=p.maxX;
}
  
MSpectrum::~MSpectrum(){
  delete spec;
  delete precursor;
  delete singlets;
  if(xCorrSparseArray!=NULL) free(xCorrSparseArray);

  while(singletFirst!=NULL){
    mScoreCard* tmp=singletFirst;
    singletFirst=singletFirst->next;
    delete tmp;
  }
  singletLast=NULL;

  int j;
  if(kojakSparseArray!=NULL){
    for(j=0;j<kojakBins;j++){
      if(kojakSparseArray[j]!=NULL) delete [] kojakSparseArray[j];
    }
    delete [] kojakSparseArray;
  }

  delete [] hp;

  for (j = 0; j<maxPepLen+1; j++){
    if (mHisto[j] != NULL) delete mHisto[j];
  }
  delete[] mHisto;

  decoys=NULL;
}


/*============================
  Operators
============================*/
MSpectrum& MSpectrum::operator=(const MSpectrum& p){
  if(this!=&p){
    delete spec;
    spec = new vector<mSpecPoint>(*p.spec);
    for(int i=0;i<20;i++) topHit[i]=p.topHit[i];
    delete precursor;
    precursor = new vector<mPrecursor>(*p.precursor);
    delete singlets;
    singlets = new vector<MTopPeps>(*p.singlets);

    binOffset = p.binOffset;
    binSize = p.binSize;
    instrumentPrecursor = p.instrumentPrecursor;
    invBinSize = p.invBinSize;
    charge = p.charge;
    maxIntensity = p.maxIntensity;
    mz = p.mz;
    scanNumber = p.scanNumber;
    rTime = p.rTime;
    xCorrArraySize = p.xCorrArraySize;
    xCorrSparseArraySize = p.xCorrSparseArraySize;
    lowScore = p.lowScore;
    bigMonoMass = p.bigMonoMass;
    bigZ=p.bigZ;

    for (int i = 0; i<HISTOSZ; i++) histogram[i] = p.histogram[i];
    histogramCount = p.histogramCount;
    histoMaxIndex = p.histoMaxIndex;

    for(int i=0;i<maxPepLen+1;i++){
      if(mHisto[i]!=NULL) delete mHisto[i];
    }
    delete [] mHisto;
    decoys = p.decoys;
    maxHistogramCount = p.maxHistogramCount;
    minAdductMass = p.minAdductMass;
    maxPepLen = p.maxPepLen;
    mHisto = new MHistogram*[maxPepLen+1];
    for (int j = 0; j<p.maxPepLen+1; j++)mHisto[j] = NULL;
    for (int j = 0; j<6; j++) ionSeries[j] = p.ionSeries[j];

    cc = p.cc;
    sc = p.sc;
    peakCounts=p.peakCounts;
    maxX=p.maxX;

    delete [] hp;
    hpSize = p.hpSize;
    hp = new sHistoPep[hpSize];
    for (int j = 0; j<hpSize; j++)hp[j] = p.hp[j];

    singletCount=p.singletCount;
    singletMax=p.singletMax;
    singletFirst=NULL;
    singletLast=NULL;
    mScoreCard* sc=NULL;
    mScoreCard* tmp=p.singletFirst;
    if(tmp!=NULL) {
      singletFirst=new mScoreCard(*tmp);
      sc=singletFirst;
      tmp=tmp->next;
      while(tmp!=NULL){
        sc->next=new mScoreCard(*tmp);
        sc->next->prev=sc;
        sc=sc->next;
        tmp=tmp->next;
      }
      singletLast=sc;
    }

    if(xCorrSparseArray!=NULL) free(xCorrSparseArray);
    if(p.xCorrSparseArray==NULL){
      xCorrSparseArray=NULL;
    } else {
      xCorrSparseArray = (mSparseMatrix *)calloc((size_t)xCorrSparseArraySize, (size_t)sizeof(mSparseMatrix));
      for(int j=0;j<xCorrSparseArraySize;j++) xCorrSparseArray[j]=p.xCorrSparseArray[j];
    }
    
    if(kojakSparseArray!=NULL){
      for(int j=0;j<kojakBins;j++){
        if(kojakSparseArray[j]!=NULL) delete [] kojakSparseArray[j];
      }
      delete [] kojakSparseArray;
    }
    kojakBins=p.kojakBins;
    if(p.kojakSparseArray==NULL){
      kojakSparseArray=NULL;
    } else {
      for(int j=0;j<kojakBins;j++){
        if(p.kojakSparseArray[j]==NULL){
          kojakSparseArray[j]=NULL;
        } else {
          kojakSparseArray[j] = new char[(int)invBinSize+1];
          for(int i=0;i<(int)invBinSize+1;i++) kojakSparseArray[j][i]=p.kojakSparseArray[j][i];
        }
      }
    }
  }
  return *this;
}

mSpecPoint& MSpectrum::operator [](const int &i){
  return spec->at(i);
}


/*============================
  Accessors
============================*/
double MSpectrum::getBinOffset(){
  return binOffset;
}

int MSpectrum::getCharge(){
  return charge;
}

bool MSpectrum::getInstrumentPrecursor(){
  return instrumentPrecursor;
}

double MSpectrum::getInvBinSize(){
  return invBinSize;
}

float MSpectrum::getMaxIntensity(){
  return maxIntensity;
}

double MSpectrum::getMZ(){
  return mz;
}

mPrecursor& MSpectrum::getPrecursor(int i){
  return precursor->at(i);
}

mPrecursor* MSpectrum::getPrecursor2(int i){
  return &precursor->at(i);
}

float MSpectrum::getRTime(){
  return rTime;
}

int MSpectrum::getScanNumber(){
  return scanNumber;
}

mScoreCard& MSpectrum::getScoreCard(int i){
  return topHit[i];
}

int MSpectrum::getSingletCount(){
  return singletCount;
}

mScoreCard& MSpectrum::getSingletScoreCard(int i){
  if(i>=singletCount) return *singletLast;
  mScoreCard* sc=singletFirst;
  int j=0;
  while(j<i){
    if(sc->next==NULL) break;
    sc=sc->next;
    j++;
  }
  return *sc;
}

MTopPeps* MSpectrum::getTopPeps(int i){
  return &singlets->at(i);
}

int MSpectrum::size(){
  return (int)spec->size();
}

int MSpectrum::sizePrecursor(){
  return (int)precursor->size();
}


/*============================
  Modifiers
============================*/
void MSpectrum::addPoint(mSpecPoint& s){
  spec->push_back(s);
}

void MSpectrum::addPrecursor(mPrecursor& p, int sz){
  precursor->push_back(p);
  MTopPeps tp;
  tp.peptideMax=sz;
  //tp.resetSingletList(p.monoMass);
  singlets->push_back(tp);
  if(p.monoMass>bigMonoMass) {
    bigMonoMass=p.monoMass;
    bigZ=p.charge;
  }
}

void MSpectrum::clear(){
  spec->clear();
  precursor->clear();
  singlets->clear();
  bigMonoMass=0;
  bigZ=0;
}

void MSpectrum::clearPrecursors(){
  precursor->clear();
  singlets->clear();
}

void MSpectrum::erasePrecursor(int i){
  precursor->erase(precursor->begin()+i);
  singlets->erase(singlets->begin()+i);
}

void MSpectrum::setCharge(int i){
  charge=i;
}

void MSpectrum::setInstrumentPrecursor(bool b){
  instrumentPrecursor=b;
}

void MSpectrum::setMaxIntensity(float f){
  maxIntensity=f;
}

void MSpectrum::setMZ(double d){
  mz=d;
}

void MSpectrum::setNativeID(string s){
  nativeID = s;
}

void MSpectrum::setRTime(float f){
  rTime=f;
}

void MSpectrum::setScanNumber(int i){
  scanNumber=i;
}

/*============================
  Functions
============================*/
bool MSpectrum::checkReporterIon(double mz, mParams* params){
  int key = (int)mz;
  if (key >= kojakBins) return false;
  if (kojakSparseArray[key] == NULL) return false;
  int pos = (int)((mz - key)*invBinSize);
  if(kojakSparseArray[key][pos]>params->rIonThreshold) return true;
  return false;
}

void MSpectrum::checkScore(mScoreCard& s, int th){
  unsigned int i;
  unsigned int j;

  //edge case for "reversible" cross-links: check if already matches top hit identically
  //note that such duplications still occur below the top score, but shouldn't influence the final result to the user
  int k=0;
  while(k<20 && s.eVal==topHit[k].eVal){
    if(s.pep==topHit[k].pep && s.site==topHit[k].site){
      if(s.mods->size()==topHit[k].mods->size()){
        for(i=0;i<s.mods->size();i++){
          if(s.mods->at(i).mass!=topHit[k].mods->at(i).mass || s.mods->at(i).pos!=topHit[k].mods->at(i).pos) break;
        }
        if(i==s.mods->size()) return;
      }
    }
    k++;
  }

  for(i=0;i<20;i++){
    if(s.eVal < topHit[i].eVal) {
      for(j=19;j>i;j--) {
        topHit[j]=topHit[j-1];
      }
      topHit[i] = s;
      lowScore=topHit[19].eVal;
      return;
    }
  }
}

//void MSpectrum::checkSingletScore(mScoreCard& s){
//
//  mScoreCard* sc;
//  mScoreCard* cur;
//  
//  //If list is empty, add the score card
//  if(singletCount==0){
//    singletFirst=new mScoreCard(s);
//    singletLast=singletFirst;
//    singletCount++;
//    return;
//  }  
//
//  //check if we can just add to the end
//  if(s.simpleScore<singletLast->simpleScore){
//    //check if we need to store the singlet
//    if(singletCount==singletMax) return;
//
//    singletLast->next=new mScoreCard(s);
//    singletLast->next->prev=singletLast;
//    singletLast=singletLast->next;
//    singletCount++;
//    return;
//  }
//  
//  //check if it goes in the front
//  if(s.simpleScore>=singletFirst->simpleScore){
//
//    singletFirst->prev=new mScoreCard(s);
//    singletFirst->prev->next=singletFirst;
//    singletFirst=singletFirst->prev;
//
//    //add to singlet list
//    if(singletCount<singletMax) {
//      singletCount++;
//    } else {
//      cur=singletLast;
//      singletLast=singletLast->prev;
//      singletLast->next=NULL;
//      delete cur;
//    }
//    return;
//  }
//
//
//  //scan to find insertion point
//  cur = singletFirst->next;
//  int i=1;
//  while(s.simpleScore < cur->simpleScore){
//    i++;
//    cur=cur->next;
//  }
//
//  sc=new mScoreCard(s);
//  sc->prev=cur->prev;
//  sc->next=cur;
//  cur->prev->next=sc;
//  cur->prev=sc;
//  if(sc->prev==NULL) singletFirst=sc;
//
//  if(singletCount<singletMax) {
//    singletCount++;
//  } else {
//    cur=singletLast;
//    singletLast=singletLast->prev;
//    singletLast->next=NULL;
//    delete cur;
//  }
//
//}

double MSpectrum::computeE(double score, int len){
  //if we haven't computed a histogram yet, do it now
  if(mHisto[len]==NULL){
    //cout << "Compute Histo: " << scanNumber << "\t" << len << endl;
    //Clear histogram
    for (int j = 0; j<HISTOSZ; j++) histogram[j] = 0;
    histogramCount = 0;

    if (!generateXcorrDecoys2(len)) {
      //handle failure
      cout << "Fail genXCorr2: " << scanNumber << "\t" << len << endl;
      exit(1);
    }

    mHisto[len]=new MHistogram();
    linearRegression3(mHisto[len]->slope, mHisto[len]->intercept, mHisto[len]->rSq);
    //cout << "linReg3 success" << endl;
    mHisto[len]->slope*=10;
    
    //** temporary
    //cout << "\nDecoy Histogram: " << len << "\t" << histogramCount << endl;
    //for(int i=0;i<HISTOSZ;i++) cout << i << "\t" << histogram[i] << endl;
    //cout << mHisto[len]->slope << "\t" <<  mHisto[len]->intercept << "\t" << mHisto[len]->rSq << endl;
    //**
  }

  if (mHisto[len]->slope>=0){
    //handle bad slopes
    //cout << "Fail regression: " << scanNumber << "\t" << len << endl;
    //cout << histogramCount << endl;
    //for (int j = 0; j<HISTOSZ; j++) cout << j << "\t" << histogram[j] << endl;
    //cout << score << "\t" << len << "\tslope=" << mHisto[len]->slope << endl;
    return 1000;
    //exit(1);
  }

  //cout << "Compute E: " << scanNumber << "\t" << len << endl;
  return pow(10.0, mHisto[len]->slope * score + mHisto[len]->intercept);

}

void MSpectrum::resetSingletList(){
  /*
  size_t j;
  double max;
  if (singletList != NULL){
    for (j = 0; j<singletBins; j++){
      if (singletList[j] != NULL) delete singletList[j];
    }
    delete[] singletList;
  }
  max=precursor->at(0).monoMass;
  for(j=1;j<precursor->size();j++){
    if (precursor->at(j).monoMass>max) max = precursor->at(j).monoMass;
  }
  singletBins=(int)(max/10+1);
  singletList = new list<mSingletScoreCard*>*[singletBins];
  for (j = 0; j<singletBins; j++) singletList[j] = NULL;
  */
}

void MSpectrum::sortMZ(){
  qsort(&spec->at(0),spec->size(),sizeof(mSpecPoint),compareMZ);
}

//void MSpectrum::xCorrScore(){
//  kojakXCorr();
//}


/*============================
  Private Functions
============================*/
void MSpectrum::kojakXCorr(double* pdTempRawData, double* pdTmpFastXcorrData, float* pfFastXcorrData, mPreprocessStruct*& pPre){
  int i;
  int j;
  int iTmp;
  double dTmp;
  double dSum;

  pPre->iHighestIon = 0;
  pPre->dHighestIntensity = 0;
  BinIons(pPre);

  memset(pdTempRawData, 0, xCorrArraySize*sizeof(double));
  memset(pdTmpFastXcorrData, 0, xCorrArraySize*sizeof(double));
  memset(pfFastXcorrData, 0, xCorrArraySize*sizeof(float));
  kojakSparseArray = new char*[kojakBins];
  for (i = 0; i<kojakBins; i++) kojakSparseArray[i] = NULL;


  // Create data for correlation analysis.
  MakeCorrData(pdTempRawData, pPre, 50.0);

  // Make fast xcorr spectrum.
  double dm = 1.0 / 150;
  mSpecPoint *pdCorrelationData = pPre->pdCorrelationData;
  dSum = 0.0;
  for (i = 0; i<75; i++) dSum += pdCorrelationData[i].intensity;
  for (i = 75; i < xCorrArraySize + 75; i++) {
    if (i<xCorrArraySize && pdCorrelationData[i].intensity>0) dSum += pdCorrelationData[i].intensity;
    if (i >= 151 && pdCorrelationData[i - 151].intensity>0) dSum -= pdCorrelationData[i - 151].intensity;
    pdTmpFastXcorrData[i - 75] = (dSum - pdCorrelationData[i - 75].intensity)* dm;
  }

  xCorrSparseArraySize = 1;

  //double dTmp0 = pdCorrelationData[0].intensity - pdTmpFastXcorrData[0];
  //double dTmp1 = pdCorrelationData[1].intensity - pdTmpFastXcorrData[1];
  //double dTmp2 = pdCorrelationData[2].intensity - pdTmpFastXcorrData[2];
  //pfFastXcorrData[0] = (float)(dTmp0 + dTmp1*0.5);
  //pfFastXcorrData[1] = (float)(dTmp1 + (dTmp0 + dTmp2)*0.5);
  //for (i = 2; i<xCorrArraySize - 1; i++){
  //  dTmp0 = dTmp1;
  //  dTmp1 = dTmp2;
  //  dTmp2 = pdCorrelationData[i + 1].intensity - pdTmpFastXcorrData[i + 1];
  //  pfFastXcorrData[i] = (float)(dTmp1 + (dTmp0 + dTmp2)*0.5);
  //}
  //pfFastXcorrData[xCorrArraySize - 1] = (float)(dTmp2 + dTmp1*0.5);
  pfFastXcorrData[0] = 0;
  //double dTmp;
  for (i = 1;i < xCorrArraySize-1;i++) {
    dTmp = pdCorrelationData[i].intensity - pdTmpFastXcorrData[i];
    pfFastXcorrData[i] = (float)dTmp;

    //TODO: Parameterize this
    // Allow user to set flanking peaks
    if (true) {
      iTmp = i - 1;
      pfFastXcorrData[i] += (float)((pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp]) * 0.5);

      iTmp = i + 1;
      pfFastXcorrData[i] += (float)((pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp]) * 0.5);
    }
  }

  //MH: Fill sparse matrix
  for (i = 0; i<xCorrArraySize; i++){
    if (pfFastXcorrData[i]>0.5 || pfFastXcorrData[i]<-0.5){

      dTmp = binSize*i;
      iTmp = (int)dTmp;
      
      if (kojakSparseArray[iTmp] == NULL) {
        kojakSparseArray[iTmp] = new char[(int)invBinSize + 1];
        for (j = 0; j<(int)invBinSize + 1; j++) kojakSparseArray[iTmp][j] = 0;
      }
      j = (int)((dTmp - iTmp)*invBinSize+0.5);

      if (pfFastXcorrData[i]>127) kojakSparseArray[iTmp][j] = 127;
      else if (pfFastXcorrData[i]<-128) kojakSparseArray[iTmp][j] = -128;
      else if (pfFastXcorrData[i]>0) kojakSparseArray[iTmp][j] = (char)(pfFastXcorrData[i] + 0.5);
      else kojakSparseArray[iTmp][j] = (char)(pfFastXcorrData[i] - 0.5);
    }
  }

}

void MSpectrum::BinIons(mPreprocessStruct *pPre) {
  int i;
  unsigned int j;
  double dPrecursor;
  double dIon;
  double dIntensity;
  mSpecPoint *pdCorrelationData = pPre->pdCorrelationData;

  // Just need to pad iArraySize by 75.
  dPrecursor=0;
  for(j=0;j<precursor->size();j++){
    if(precursor->at(j).monoMass>dPrecursor) dPrecursor=precursor->at(j).monoMass;
  }
  xCorrArraySize = (int)((dPrecursor + 100.0) / binSize);
  if(xCorrArraySize>pPre->iMaxXCorrArraySize) xCorrArraySize=pPre->iMaxXCorrArraySize;
  kojakBins = (int)(spec->at(spec->size()-1).mass+100.0);

  memset(pdCorrelationData, 0, xCorrArraySize*sizeof(mSpecPoint));

  for(i=0;i<(int)spec->size();i++){
    dIon = spec->at(i).mass;
    dIntensity = spec->at(i).intensity;   

    if (dIntensity > 0.0) {
      if (dIon < (dPrecursor + 50.0)) {

        //#define BIN(dMass) (int)(dMass*invBinSize + binOffset)
        int iBinIon = (int)(dIon*invBinSize+binOffset);
        dIntensity = sqrt(dIntensity);
        if (iBinIon > pPre->iHighestIon) pPre->iHighestIon = iBinIon;

        if ((iBinIon < xCorrArraySize) && (dIntensity > pPre->pdCorrelationData[iBinIon].intensity)) {
          if (dIntensity > pPre->pdCorrelationData[iBinIon].intensity) {
            pPre->pdCorrelationData[iBinIon].intensity = (float)dIntensity;
            pPre->pdCorrelationData[iBinIon].mass = dIon;
          }
          if (pPre->pdCorrelationData[iBinIon].intensity > pPre->dHighestIntensity) pPre->dHighestIntensity = pPre->pdCorrelationData[iBinIon].intensity;    
        }
      }
    }
  }

  //Clear spectrum data that we no longer need
  spec->clear();

}

// pdTempRawData now holds raw data, pdCorrelationData is windowed data.
void MSpectrum::MakeCorrData(double *pdTempRawData, mPreprocessStruct *pPre, double scale){
  int  i;
  int  ii;
  int  iBin;
  int  iWindowSize;
  int  iNumWindows=10;
  double dMaxWindowInten;
  double dMaxOverallInten;
  double dTmp1;
  double dTmp2;

  dMaxOverallInten = 0.0;

  // Normalize maximum intensity to 100.
  dTmp1 = 1.0;
  if (pPre->dHighestIntensity > 0.000001) dTmp1 = 100.0 / pPre->dHighestIntensity;

  for (i=0; i < xCorrArraySize; i++) {
    pdTempRawData[i] = pPre->pdCorrelationData[i].intensity*dTmp1;
    pPre->pdCorrelationData[i].intensity=0.0;
    if (dMaxOverallInten < pdTempRawData[i]) dMaxOverallInten = pdTempRawData[i];
  }

  iWindowSize = (int) ceil( (double)(pPre->iHighestIon)/iNumWindows);

  for (i=0; i<iNumWindows; i++){
    dMaxWindowInten = 0.0;
    for (ii=0; ii<iWindowSize; ii++) {   // Find max inten. in window.
      iBin = i*iWindowSize+ii;
      if (iBin < xCorrArraySize) {
        if (pdTempRawData[iBin] > dMaxWindowInten)dMaxWindowInten = pdTempRawData[iBin];
      }
    }

    if (dMaxWindowInten > 0.0) {
      dTmp1 = scale / dMaxWindowInten;
      dTmp2 = 0.05 * dMaxOverallInten;

      for (ii=0; ii<iWindowSize; ii++){    // Normalize to max inten. in window.      
        iBin = i*iWindowSize+ii;
        if (iBin < xCorrArraySize){
          if (pdTempRawData[iBin] > dTmp2) pPre->pdCorrelationData[iBin].intensity = (float)(pdTempRawData[iBin]*dTmp1);
        }
      }
    }
  }

}

//from Comet
bool MSpectrum::calcEValue(mParams* params, MDecoys& decoys) {
  int i;
  int iLoopCount;
  int iMaxCorr;
  int iStartCorr;
  int iNextCorr;
  double dSlope;
  double dIntercept;
  double dRSq;
  if(topHit[0].simpleScore==0) return true;

  if (histogramCount < DECOY_SIZE) {
    if (!generateXcorrDecoys(params,decoys)) return false;
  }

  //linearRegression(dSlope, dIntercept, iMaxCorr, iStartCorr, iNextCorr);
  linearRegression2(dSlope, dIntercept, iMaxCorr, iStartCorr, iNextCorr,dRSq);
  histoMaxIndex = iMaxCorr;

  dSlope *= 10.0;

  iLoopCount=20; //score all e-values among top hits?
  for (i = 0; i<iLoopCount; i++) {
    if(topHit[i].simpleScore==0) break; //score all e-values among top hits?
    if (dSlope >= 0.0) {
      topHit[i].eVal=1e12;
    } else {
      topHit[i].eVal = pow(10.0, dSlope * topHit[i].simpleScore + dIntercept);
    }
  }

  return true;
}

//from Comet
// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
bool MSpectrum::generateXcorrDecoys(mParams* params, MDecoys& decoys) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  size_t r=0;
  int key;
  int pos;
  double dBion;
  double dYion;
  double dXcorr;
  double dFragmentIonMass = 0.0;

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = DECOY_SIZE - histogramCount;
  int seed = (scanNumber*histogramCount);
  if (seed<0) seed = -seed;
  seed = seed % DECOY_SIZE; //don't always start at the top, but not random either; remains reproducible across threads
  int decoyIndex;

  j = 0;
  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % DECOY_SIZE;

    //iterate over precursors
    r++;
    if (r >= precursor->size()) r = 0;
    maxZ = precursor->at(r).charge;
    if (maxZ>4) maxZ = 4;

    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      dBion = decoys.decoyIons[decoyIndex].pdIonsN[j];
      dYion = decoys.decoyIons[decoyIndex].pdIonsC[j];
      if(dBion>precursor->at(r).monoMass && dYion>precursor->at(r).monoMass) break; //stop when fragment ion masses exceed precursor mass

      for (n = 0; n<6; n++) {
        if(!params->ionSeries[n]) continue;
        switch (n) {
        case 0: dFragmentIonMass = dBion - 27.9949141; break;
        case 1: dFragmentIonMass = dBion; break;
        case 2: dFragmentIonMass = dBion + 17.026547; break;
        case 3: dFragmentIonMass = dYion + 25.9792649; break;
        case 4: dFragmentIonMass = dYion; break;
        case 5: dFragmentIonMass = dYion - 16.0187224; break;
        }
        if(dFragmentIonMass>precursor->at(r).monoMass) continue;

        for (z=1;z<maxZ;z++) {
          mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
          key = (int)mz;
          if(key>=kojakBins) break;
          if(kojakSparseArray[key]==NULL) continue;
          pos = (int)((mz - key)*invBinSize);
          dXcorr += kojakSparseArray[key][pos];
        }
      }
    }

    if(dXcorr<=0.0) dXcorr=0.0;
    k = (int)(dXcorr*0.05 + 0.5);  // 0.05=0.005*10; see MAnalysis::mangnumScoring
    if (k < 0) k = 0;
    else if (k >= HISTOSZ) k = HISTOSZ - 1;
    histogram[k] += 1;
    histogramCount++;
  }

  return true;
}

//from Comet
// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
bool MSpectrum::generateXcorrDecoys2(int maxPepLen) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  size_t r = 0;
  int key;
  int pos;
  double dBion[3]; //0=N-term, 1=C-term, 2=center
  double dYion[3];
  double dXcorr[3];
  double dFragmentIonMass = 0.0;
  double oMass;
  bool bAdduct;
  int ac;


  

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = DECOY_SIZE;
  int seed = scanNumber % DECOY_SIZE; //don't always start at the top, but not random either; remains reproducible across threads
  int decoyIndex;

  //cout << "generateXcorrDecoys2 " << iLoopMax << "\t" << precursor->size() << endl;

  j = 0;
  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr[0]=dXcorr[1]=dXcorr[2] = 0;
    decoyIndex = (seed + i) % DECOY_SIZE;

    //iterate over precursors
    r++;
    if (r >= precursor->size()) r = 0;
    maxZ = precursor->at(r).charge;
    if (maxZ>4) maxZ = 4;
    
    //if we need to add an adduct mass, then do so here
    bAdduct=true;
    ac=3;
    oMass=precursor->at(r).monoMass-110.5822*maxPepLen;
    //cout << "oMass: " << oMass << "\t" << minAdductMass << "\t" << precursor->at(r).monoMass << endl;
    if(oMass<minAdductMass) {
      oMass=0;
      bAdduct=false;
      ac=1;
    }

    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      if(j==maxPepLen) break;    //stop when fragment ions exceed peptide length
      dBion[0] = decoys->decoyIons[decoyIndex].pdIonsN[j] + oMass;
      dYion[0] = decoys->decoyIons[decoyIndex].pdIonsC[j];
      if(bAdduct){
        dBion[1] = decoys->decoyIons[decoyIndex].pdIonsN[j];
        dYion[1] = decoys->decoyIons[decoyIndex].pdIonsC[j] + oMass;
        if(j<maxPepLen/2) {
          dBion[2]=dBion[1];
          dYion[2]=dYion[0];
        } else {
          dBion[2]=dBion[0];
          dYion[2]=dYion[1];
        }
      }

      for(int x=0;x<ac;x++){
        for (n = 0; n<6; n++) {
          if (!ionSeries[n]) continue;
          switch (n) {
          case 0: dFragmentIonMass = dBion[x] - 27.9949141; break;
          case 1: dFragmentIonMass = dBion[x]; break;
          case 2: dFragmentIonMass = dBion[x] + 17.026547; break;
          case 3: dFragmentIonMass = dYion[x] + 25.9792649; break;
          case 4: dFragmentIonMass = dYion[x]; break;
          case 5: dFragmentIonMass = dYion[x] - 16.0187224; break;
          }
          if (dFragmentIonMass>precursor->at(r).monoMass) continue;

          for (z = 1; z<maxZ; z++) {
            mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
            mz = binSize * (int)(mz*invBinSize + binOffset);
            key = (int)mz;
            if (key >= kojakBins) break;
            if (key<0 || kojakSparseArray[key] == NULL) continue;
            pos = (int)((mz - key)*invBinSize);
            dXcorr[x] += kojakSparseArray[key][pos];
          }
        }
      }
    }

    for(int x=0;x<ac;x++){
      if (dXcorr[x] <= 0.0) dXcorr[x] = 0.0;
      k = (int)(dXcorr[x]*0.05 + 0.5);  // 0.05=0.005*10; see MAnalysis::mangnumScoring
      if (k < 0) k = 0;
      else if (k >= HISTOSZ) k = HISTOSZ - 1;
      histogram[k] += 1;
      histogramCount++;
      if(histogramCount==maxHistogramCount) break;
    }
    if (histogramCount == maxHistogramCount) break;
  }

  //cout << "success" << endl;
  return true;
}

//from Comet
// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
bool MSpectrum::generateXcorrDecoys3(int minP, int maxP, int depth) {
  //int i;
  //int n;
  //int j;
  //int k;
  //int maxZ;
  //int z;
  //size_t r = 0;
  //int key;
  //int pos;
  //double dBion[3]; //0=N-term, 1=C-term, 2=center
  //double dYion[3];
  //double dXcorr[3];
  //double dFragmentIonMass = 0.0;
  //double oMass;
  //bool bAdduct;
  //int ac;


  //iterate from 0 to maxlen
  //compute yn;
  //if n>minlen/2, function(n+1 to n*2,+offset)
  //if n>minlen, function(b+offset);
  //compute bn;
  //if n>minlex, function(y+offset);

  //function(n + 1 to n * 2, +offset)
  //function( b+offset, upto n);
  //compute yn+offset
  //compute bn;  <-this is redundant for all lengths

  int histoXCount[MAX_DECOY_PEP_LEN];
  int histoX[MAX_DECOY_PEP_LEN][HISTOSZ];
  for(int a=0;a<MAX_DECOY_PEP_LEN;a++){
    histoXCount[a]=0;
    for(int b=0;b<HISTOSZ;b++){
      histoX[a][b]=0;
    }
  }

  //int histoXCount[51];
  //int histoX[51][HISTOSZ];
  //for (int a = 0; a<51; a++){
  //  histoXCount[a] = 0;
  //  for (int b = 0; b<HISTOSZ; b++){
  //    histoX[a][b] = 0;
  //  }
  //}

  int decoyIndex=0;
  double monoMass=bigMonoMass;
  double ionB,ionY;
  double dFragmentIonMass;
  double xcorr,xcorr2B, xcorr2Y;
  double m;
  int maxZ=bigZ;
  if (maxZ>4) maxZ = 4;
  int key,pos;
  double xcorrY,xcorrB;
  for(int x=0;x<depth/3;x++){
    
    //cout << "MonoMass: " << monoMass << endl;
    decoyIndex=x;
    xcorrY=0;
    xcorrB=0;
    for(int a=0;a<maxP-1;a++){  

      //the unmodified ions
      ionB = decoys->decoyIons[decoyIndex].pdIonsN[a];
      //cout << "B" << a + 1 << "\t" << ionB << endl;
      ionY = decoys->decoyIons[decoyIndex].pdIonsC[a];
      //cout << "Y" << a+1 << "\t" << ionY << endl;
      for (int b = 3; b<6; b++) {
        if (!ionSeries[b]) continue;
        switch (b) {
        case 0: dFragmentIonMass = ionB - 27.9949141; break;
        case 1: dFragmentIonMass = ionB; break;
        case 2: dFragmentIonMass = ionB + 17.026547; break;
        case 3: dFragmentIonMass = ionY + 25.9792649; break;
        case 4: dFragmentIonMass = ionY; break;
        case 5: dFragmentIonMass = ionY - 16.0187224; break;
        }
        //if (dFragmentIonMass>monoMass) continue;

        for (int z = 1; z<maxZ; z++) {
          m = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          m = binSize * (int)(m*invBinSize + binOffset);
          key = (int)m;
          if (key >= kojakBins) break;
          if (key<0 || kojakSparseArray[key] == NULL) continue;
          pos = (int)((m - key)*invBinSize);
          //cout << (int)kojakSparseArray[key][pos] << endl;
          if (b<3) xcorrB += kojakSparseArray[key][pos];
          else xcorrY += kojakSparseArray[key][pos];
        }
      }
      if(a>(minP-3)) {
        //Modification on N-term
        double massOffsetB = monoMass - decoys->decoyIons[decoyIndex].pdIonsN[a+1]-17.00273965; //17.00328821;
        //cout << "Start B: " << massOffsetB << "\t" << xcorrY << endl;
        xcorr2Y=xcorrY;
        double massOffsetY = monoMass - decoys->decoyIons[decoyIndex].pdIonsC[a + 1] + 1.00782503;
        //cout << "Start Y: " << massOffsetY << "\t" << xcorrB << endl;
        xcorr2B = xcorrB;
        for (int b = 0; b<6; b++) {
          if (!ionSeries[b]) continue;
          switch (b) {
          case 0: dFragmentIonMass = massOffsetB - 27.9949141; break;
          case 1: dFragmentIonMass = massOffsetB; break;
          case 2: dFragmentIonMass = massOffsetB + 17.026547; break;
          case 3: dFragmentIonMass = massOffsetY - 27.9949141; break;
          case 4: dFragmentIonMass = massOffsetY; break;
          case 5: dFragmentIonMass = massOffsetY + 17.026547; break;
          }
          if(b<3) xcorr2Y += makeXCorrB(decoyIndex, dFragmentIonMass, maxZ, a);
          else xcorr2B += makeXCorrY(decoyIndex, dFragmentIonMass, maxZ, a);
        }
        //get key
        //cout << "TotY: " << xcorr2Y << "\t" << a+2 << endl;
        if (xcorr2Y <= 0.0) xcorr2Y = 0.0;
        int k = (int)(xcorr2Y * 0.05 + 0.5);  // 0.05=0.005*10; see MAnalysis::mangnumScoring
        if (k < 0) k = 0;
        else if (k >= HISTOSZ) k = HISTOSZ - 1;
        histoX[a+2][k]++;  //N-terminal open mod
        histoXCount[a+2]++;
        //cout << "TotB: " << xcorr2B << "\t" << a + 2 << endl;
        if (xcorr2B <= 0.0) xcorr2B = 0.0;
        k = (int)(xcorr2B * 0.05 + 0.5);  // 0.05=0.005*10; see MAnalysis::mangnumScoring
        if (k < 0) k = 0;
        else if (k >= HISTOSZ) k = HISTOSZ - 1;
        histoX[a + 2][k]++;  //N-terminal open mod
        histoXCount[a + 2]++;
      }

      //Middle mod - for all peptides
      if(a>=(minP-2)/2 && a<=(maxP-2)/2){ 
        xcorr=xcorrB+xcorrY;
        double massOffsetB = monoMass - decoys->decoyIons[decoyIndex].pdIonsN[a*2 + 1] - 17.00273965;
        double massOffsetY = monoMass - decoys->decoyIons[decoyIndex].pdIonsC[a*2 + 1] + 1.00782503;
        for (int b = 0; b<6; b++) {
          if (!ionSeries[b]) continue;
          switch (b) {
          case 0: dFragmentIonMass = massOffsetB - 27.9949141; break;
          case 1: dFragmentIonMass = massOffsetB; break;
          case 2: dFragmentIonMass = massOffsetB + 17.026547; break;
          case 3: dFragmentIonMass = massOffsetY - 27.9949141; break;
          case 4: dFragmentIonMass = massOffsetY; break;
          case 5: dFragmentIonMass = massOffsetY + 17.026547; break;
          }
          if (b<3) xcorr += makeXCorrB(decoyIndex, dFragmentIonMass, maxZ, a*2,a+1);
          else xcorr += makeXCorrY(decoyIndex, dFragmentIonMass, maxZ, a*2,a+1);
        }
        //cout << "TotM: " << xcorr << "\t" << a*2 + 2 << endl;
        if (xcorr <= 0.0) xcorr = 0.0;
        int k = (int)(xcorr * 0.05 + 0.5);  // 0.05=0.005*10; see MAnalysis::mangnumScoring
        if (k < 0) k = 0;
        else if (k >= HISTOSZ) k = HISTOSZ - 1;
        histoX[a*2 + 2][k]++;  //N-terminal open mod
        histoXCount[a*2 + 2]++;
        
        if((a*2+1)<maxP-1){
          xcorr = xcorrB + xcorrY;
          double massOffsetB = monoMass - decoys->decoyIons[decoyIndex].pdIonsN[a * 2 + 2] - 17.00273965;
          double massOffsetY = monoMass - decoys->decoyIons[decoyIndex].pdIonsC[a * 2 + 2] + 1.00782503;
          for (int b = 0; b<6; b++) {
            if (!ionSeries[b]) continue;
            switch (b) {
            case 0: dFragmentIonMass = massOffsetB - 27.9949141; break;
            case 1: dFragmentIonMass = massOffsetB; break;
            case 2: dFragmentIonMass = massOffsetB + 17.026547; break;
            case 3: dFragmentIonMass = massOffsetY - 27.9949141; break;
            case 4: dFragmentIonMass = massOffsetY; break;
            case 5: dFragmentIonMass = massOffsetY + 17.026547; break;
            }
            if (b<3) xcorr += makeXCorrB(decoyIndex, dFragmentIonMass, maxZ, a * 2+1, a+1);
            else xcorr += makeXCorrY(decoyIndex, dFragmentIonMass, maxZ, a * 2+1, a+1);
          }
          //cout << "TotM: " << xcorr << "\t" << a * 2 + 3 << endl;
          if (xcorr <= 0.0) xcorr = 0.0;
          int k = (int)(xcorr * 0.05 + 0.5);  // 0.05=0.005*10; see MAnalysis::mangnumScoring
          if (k < 0) k = 0;
          else if (k >= HISTOSZ) k = HISTOSZ - 1;
          histoX[a * 2 + 3][k]++;  //N-terminal open mod
          histoXCount[a * 2 + 3]++;
        }
      }

    }

  }
          
  for(int a=minP;a<maxP+1;a++){
    mHisto[a] = new MHistogram();
    linearRegression4(histoX[a],histoXCount[a],mHisto[a]->slope, mHisto[a]->intercept, mHisto[a]->rSq);
    mHisto[a]->slope *= 10;

    //** temporary
    //cout << "\nDecoy Histogram: " << a << "\t" << histoXCount[a] << endl;
    //for (int i = 0; i<HISTOSZ; i++) cout << i << "\t" << histoX[a][i] << endl;
    //cout << mHisto[a]->slope << "\t" << mHisto[a]->intercept << "\t" << mHisto[a]->rSq << endl;
    //**

  }

  return true;
}

//from Comet
//Drastically simplifying this process. Future steps would be to use a 
//decoy library index to rapidly calculate these values, should the memory be available.
bool MSpectrum::generateXcorrDecoys4(int minP, int maxP, int depth) {
  //cout << "generateXcorrDecoys4: " << scanNumber << endl;
  int histoXCount[MAX_DECOY_PEP_LEN];
  int histoX[MAX_DECOY_PEP_LEN][HISTOSZ];
  for (int a = 0;a < MAX_DECOY_PEP_LEN;a++) {
    histoXCount[a] = 0;
    for (int b = 0;b < HISTOSZ;b++) {
      histoX[a][b] = 0;
    }
  }

  int decoyIndex = 0;
  double monoMass = bigMonoMass;
  double ionB, ionY;
  double dFragmentIonMass;
  double xcorr=0;
  int k;
  double m;
  int maxZ = bigZ;
  if (maxZ > 4) maxZ = 4;
  int key, pos;

  for (int x = 0;x < depth;x++) {

    //cout << "MonoMass: " << monoMass << "\t" << "Depth: " << x << endl;
    decoyIndex = x;
    xcorr = 0;
    for (int a = 0;a < maxP;a++) {

      //the unmodified ions
      ionB = decoys->decoyIons[decoyIndex].pdIonsN[a];
      //cout << "B" << a + 1 << "\t" << ionB << endl;
      ionY = decoys->decoyIons[decoyIndex].pdIonsC[a];
      //cout << "Y" << a+1 << "\t" << ionY << endl;
      for (int b = 0; b < 6; b++) {
        if (!ionSeries[b]) continue;
        switch (b) {
        case 0: dFragmentIonMass = ionB - 27.9949141; break;
        case 1: dFragmentIonMass = ionB; break;
        case 2: dFragmentIonMass = ionB + 17.026547; break;
        case 3: dFragmentIonMass = ionY + 25.9792649; break;
        case 4: dFragmentIonMass = ionY; break;
        case 5: dFragmentIonMass = ionY - 16.0187224; break;
        }
        //if (dFragmentIonMass>monoMass) continue;

        for (int z = 1; z < maxZ; z++) {
          m = (dFragmentIonMass + (z - 1) * 1.007276466) / z;
          m = binSize * (int)(m * invBinSize + binOffset);
          key = (int)m;
          if (key >= kojakBins) break;
          if (key < 0 || kojakSparseArray[key] == NULL) continue;
          pos = (int)((m - key) * invBinSize+0.5);
          //cout << (int)kojakSparseArray[key][pos] << endl;
          xcorr += kojakSparseArray[key][pos];
        }
      }

      if (xcorr <= 0.0) k = 0;
      else k = (int)(xcorr * 0.05 + 0.5);  // 0.05=0.005*10; see MAnalysis::mangnumScoring
      if (k < 0) k = 0;
      else if (k >= HISTOSZ) k = HISTOSZ - 1;
      histoX[a+2][k]++;   //the first ion results from a peptide of len=2
      histoXCount[a+2]++;

    }

  }

  for (int a = minP;a < maxP + 1;a++) {
    mHisto[a] = new MHistogram();
    linearRegression4(histoX[a], histoXCount[a], mHisto[a]->slope, mHisto[a]->intercept, mHisto[a]->rSq);
    mHisto[a]->slope *= 10;

    //** temporary
    //cout << "\nDecoy Histogram: " << a << "\t" << histoXCount[a] << endl;
    //for (int i = 0; i<HISTOSZ; i++) cout << i << "\t" << histoX[a][i] << endl;
    //cout << mHisto[a]->slope << "\t" << mHisto[a]->intercept << "\t" << mHisto[a]->rSq << endl;
    //**

  }

  return true;
}

double MSpectrum::makeXCorrB(int decoyIndex, double modMass, int maxZ, int len, int offset){
 
  double xcorr=0;
  double m;
  int key;
  int pos;
  double b;
  for(int a=offset;a<len+1;a++){
    b = modMass + decoys->decoyIons[decoyIndex].pdIonsN[a];
    if(b<0) continue;
    //if (b>maxMass) break; //stop when fragment ion masses exceed precursor mass
    //cout << "xb" << a+1 << "\t" << b << endl;
    for (int z = 1; z<maxZ; z++) {
      m = (b + (z - 1)*1.007276466) / z;
      m = binSize * (int)(m*invBinSize + binOffset);
      key = (int)m;
      if (key >= kojakBins) break;
      if (kojakSparseArray[key] == NULL) continue;
      pos = (int)((m - key)*invBinSize);
      xcorr += kojakSparseArray[key][pos];
    }
  }
  //cout << "xB: " << xcorr << endl;
  return xcorr;
}

double MSpectrum::makeXCorrY(int decoyIndex, double modMass, int maxZ, int len, int offset){

  double xcorr = 0;
  double m;
  int key;
  int pos;
  double b;
  for (int a = offset; a<len + 1; a++){
    b = modMass + decoys->decoyIons[decoyIndex].pdIonsC[a];
    if(b<0) continue;
    //cout << "xy" << a+1 << "\t" << b << endl;
    //if (b>maxMass) break; //stop when fragment ion masses exceed precursor mass
    for (int z = 1; z<maxZ; z++) {
      m = (b + (z - 1)*1.007276466) / z;
      m = binSize * (int)(m*invBinSize + binOffset);
      key = (int)m;
      if (key >= kojakBins) break;
      if (kojakSparseArray[key] == NULL) continue;
      pos = (int)((m - key)*invBinSize);
      xcorr += kojakSparseArray[key][pos];
    }
  }
  //cout << "xY: " << xcorr << endl;
  return xcorr;
}

//from Comet
void MSpectrum::linearRegression(double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  // Find maximum correlation score index.
  for (i = HISTOSZ- 2; i >= 0; i--) {
    if (histogram[i] > 0)  break;
  }
  iMaxCorr = i;

  iNextCorr = 0;
  for (i = 0; i<iMaxCorr; i++)  {
    if (histogram[i] == 0)  {
      // register iNextCorr if there's a histo value of 0 consecutively
      if (histogram[i + 1] == 0 || i + 1 == iMaxCorr) {
        if (i>0) iNextCorr = i - 1;
        break;
      }
    }
  }

  if (i == iMaxCorr) {
    iNextCorr = iMaxCorr;
    if (iMaxCorr>12) iNextCorr = iMaxCorr - 2;
  }

  // Create cummulative distribution function from iNextCorr down, skipping the outliers.
  dCummulative[iNextCorr] = histogram[iNextCorr];
  for (i = iNextCorr - 1; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histogram[i];
    if (histogram[i + 1] == 0) dCummulative[i + 1] = 0.0;
  }

  // log10
  for (i = iNextCorr; i >= 0; i--)  {
    histogram[i] = (int)dCummulative[i];  // First store cummulative in histogram.
    dCummulative[i] = log10(dCummulative[i]);
  }

  iStartCorr = 0;
  if (iNextCorr >= 30) iStartCorr = (int)(iNextCorr - iNextCorr*0.25);
  else if (iNextCorr >= 15) iStartCorr = (int)(iNextCorr - iNextCorr*0.5);

  Mx = My = a = b = 0.0;
  while (iStartCorr >= 0) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (histogram[i] > 0) {
        SumY += (float)dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }

    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (dCummulative[i] > 0) {
        dx = i - Mx;
        dy = dCummulative[i] - My;
        Sx += dx*dx;
        Sxy += dx*dy;
      }
    }

    if (Sx > 0) b = Sxy / Sx;   // slope
    else b = 0;

    if (b < 0.0) break;
    else iStartCorr--;
  }

  a = My - b*Mx;  // y-intercept

  slope = b;
  intercept = a;
  iMaxXcorr = iMaxCorr;
  iStartXcorr = iStartCorr;
  iNextXcorr = iNextCorr;
}

void MSpectrum::linearRegression2(double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.
  double SST, SSR;
  double rsq;
  double bestRSQ;
  double bestSlope;
  double bestInt;

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int bestNC;
  int bestStart;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  // Find maximum correlation score index.
  for (i = HISTOSZ - 2; i >= 0; i--) {
    if (histogram[i] > 0)  break;
  }
  iMaxCorr = i;

  //bail now if there is no width to the distribution
  if (iMaxCorr<3) {
    slope = 0;
    intercept = 0;
    iMaxXcorr = 0;
    iStartXcorr = 0;
    iNextXcorr = 0;
    rSquared = 0;
    return;
  }

  //More aggressive version summing everything below the max
  dCummulative[iMaxCorr - 1] = histogram[iMaxCorr - 1];
  for (i = iMaxCorr - 2; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histogram[i];
  }

  //get middle-ish datapoint as seed. Using count/10.
  for (i = 0; i<iMaxCorr; i++){
    if (dCummulative[i] < histogramCount / 10) break;
  }
  if (i >= (iMaxCorr - 1)) iNextCorr = iMaxCorr - 2;
  else iNextCorr = i;

  // log10...and stomp all over the original...hard to troubleshoot later
  for (i = iMaxCorr - 1; i >= 0; i--)  {
    histogram[i] = (int)dCummulative[i];
    if(dCummulative[i]>0) dCummulative[i] = log10(dCummulative[i]);
  }

  iStartCorr = iNextCorr - 1;
  iNextCorr++;

  bool bRight = false; // which direction to add datapoint from
  bestRSQ = 0;
  bestNC = 0;
  bestSlope = 0;
  bestInt = 0;
  rsq = Mx = My = a = b = 0.0;
  while (true) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (histogram[i] > 0) {
        SumY += dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }
    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    SST = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dx = i - Mx;
      dy = dCummulative[i] - My;
      Sx += dx*dx;
      Sxy += dx*dy;
      SST += dy*dy;
    }
    b = Sxy / Sx;
    a = My - b*Mx;  // y-intercept

    //MH: compute R2
    SSR = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dy = dCummulative[i] - (b*i + a);
      SSR += (dy*dy);
    }
    rsq = 1 - SSR / SST;

    if (rsq>0.95 || rsq>bestRSQ){
      if (rsq>bestRSQ || iNextCorr - iStartCorr + 1<8){ //keep better RSQ only if more than 8 datapoints, otherwise keep every RSQ below 8 datapoints
        bestRSQ = rsq;
        bestNC = iNextCorr;
        bestSlope = b;
        bestInt = a;
        bestStart = iStartCorr;
      }
      if (bRight){
        if (iNextCorr<(iMaxCorr - 1) && dCummulative[iNextCorr]>0) iNextCorr++;
        else if (iStartCorr>0) iStartCorr--;
        else break;
      } else {
        if (iStartCorr>0) iStartCorr--;
        else if (iNextCorr<(iMaxCorr - 1) && dCummulative[iNextCorr]>0) iNextCorr++;
        else break;
      }
      bRight = !bRight;
    } else {
      break;
    }
  }

  slope = bestSlope;
  intercept = bestInt;
  iMaxXcorr = iMaxCorr;
  iStartXcorr = bestStart;
  iNextXcorr = bestNC;
  rSquared = bestRSQ;

}

void MSpectrum::linearRegression3(double& slope, double& intercept, double& rSquared) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.
  double SST, SSR;
  double rsq;
  double bestRSQ;
  double bestSlope;
  double bestInt;

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  //cout << "x " << histogramCount << endl;
  //for(i=0;i<HISTOSZ;i++) cout << "x " << i << "\t" << histogram[i] << endl;

  // Find maximum correlation score index.
  for (i = HISTOSZ - 2; i >= 0; i--) {
    if (histogram[i] > 0)  break;
  }
  iMaxCorr = i;

  //bail now if there is no width to the distribution
  if (iMaxCorr<2) {
    slope = 0;
    intercept = 0;
    rSquared = 0;
    return;
  }

  //More aggressive version summing everything
  dCummulative[iMaxCorr] = histogram[iMaxCorr];
  for (i = iMaxCorr - 1; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histogram[i];
  }

  //get middle-ish datapoint as seed. Using count/10.
  for (i = 0; i<iMaxCorr; i++){
    if (dCummulative[i] < histogramCount / 10) break;
  }
  if (i >= (iMaxCorr)) iNextCorr = iMaxCorr - 1;
  else iNextCorr = i;

  // log10...and stomp all over the original...hard to troubleshoot later
  for (i = iMaxCorr; i >= 0; i--)  {
    histogram[i] = (int)dCummulative[i];
    if (dCummulative[i]>0) dCummulative[i] = log10(dCummulative[i]);
  }

  iStartCorr = iNextCorr - 1;
  iNextCorr++;

  bool bRight = false; // which direction to add datapoint from
  bestRSQ = 0;
  bestSlope = 0;
  bestInt = 0;
  rsq = Mx = My = a = b = 0.0;
  while (true) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (histogram[i] > 0) {
        SumY += dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }
    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    SST = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dx = i - Mx;
      dy = dCummulative[i] - My;
      Sx += dx*dx;
      Sxy += dx*dy;
      SST += dy*dy;
    }
    b = Sxy / Sx;
    a = My - b*Mx;  // y-intercept

    //MH: compute R2
    SSR = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dy = dCummulative[i] - (b*i + a);
      SSR += (dy*dy);
    }
    rsq = 1 - SSR / SST;

    if (rsq>0.95 || rsq>bestRSQ){
      if (rsq>bestRSQ || iNextCorr - iStartCorr + 1<8){ //keep better RSQ only if more than 8 datapoints, otherwise keep every RSQ below 8 datapoints
        bestRSQ = rsq;
        bestSlope = b;
        bestInt = a;
      }
      if (bRight){
        if (iNextCorr<iMaxCorr && dCummulative[iNextCorr]>0) iNextCorr++;
        else if (iStartCorr>0) iStartCorr--;
        else break;
      } else {
        if (iStartCorr>0) iStartCorr--;
        else if (iNextCorr<iMaxCorr && dCummulative[iNextCorr]>0) iNextCorr++;
        else break;
      }
      bRight = !bRight;
    } else {
      break;
    }
  }

  slope = bestSlope;
  intercept = bestInt;
  rSquared = bestRSQ;

}

void MSpectrum::linearRegression4(int* h, int sz, double& slope, double& intercept, double& rSquared) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.
  double SST, SSR;
  double rsq;
  double bestRSQ;
  double bestSlope;
  double bestInt;

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  //cout << "x " << histogramCount << endl;
  //for(i=0;i<HISTOSZ;i++) cout << "x " << i << "\t" << histogram[i] << endl;

  // Find maximum correlation score index.
  for (i = HISTOSZ - 2; i >= 0; i--) {
    if (h[i] > 0)  break;
  }
  iMaxCorr = i;

  //bail now if there is no width to the distribution
  if (iMaxCorr<2) {
    slope = 0;
    intercept = 0;
    rSquared = 0;
    return;
  }

  //More aggressive version summing everything
  dCummulative[iMaxCorr] = h[iMaxCorr];
  for (i = iMaxCorr - 1; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + h[i];
  }

  //get middle-ish datapoint as seed. Using count/10.
  for (i = 0; i<iMaxCorr; i++){
    if (dCummulative[i] < sz / 10) break;
  }
  if (i >= (iMaxCorr)) iNextCorr = iMaxCorr - 1;
  else iNextCorr = i;

  // log10...and stomp all over the original...hard to troubleshoot later
  for (i = iMaxCorr; i >= 0; i--)  {
    h[i] = (int)dCummulative[i];
    if (dCummulative[i]>0) dCummulative[i] = log10(dCummulative[i]);
  }

  iStartCorr = iNextCorr - 1;
  iNextCorr++;

  bool bRight = false; // which direction to add datapoint from
  bestRSQ = 0;
  bestSlope = 0;
  bestInt = 0;
  rsq = Mx = My = a = b = 0.0;
  while (true) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (h[i] > 0) {
        SumY += dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }
    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    SST = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dx = i - Mx;
      dy = dCummulative[i] - My;
      Sx += dx*dx;
      Sxy += dx*dy;
      SST += dy*dy;
    }
    b = Sxy / Sx;
    a = My - b*Mx;  // y-intercept

    //MH: compute R2
    SSR = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dy = dCummulative[i] - (b*i + a);
      SSR += (dy*dy);
    }
    rsq = 1 - SSR / SST;

    if (rsq>0.95 || rsq>bestRSQ){
      if (rsq>bestRSQ || iNextCorr - iStartCorr + 1<8){ //keep better RSQ only if more than 8 datapoints, otherwise keep every RSQ below 8 datapoints
        bestRSQ = rsq;
        bestSlope = b;
        bestInt = a;
      }
      if (bRight){
        if (iNextCorr<iMaxCorr && dCummulative[iNextCorr]>0) iNextCorr++;
        else if (iStartCorr>0) iStartCorr--;
        else break;
      } else {
        if (iStartCorr>0) iStartCorr--;
        else if (iNextCorr<iMaxCorr && dCummulative[iNextCorr]>0) iNextCorr++;
        else break;
      }
      bRight = !bRight;
    } else {
      break;
    }
  }

  slope = bestSlope;
  intercept = bestInt;
  rSquared = bestRSQ;

}

bool MSpectrum::matchMods(mPepMod2& v1, vector<mPepMod>& v2){
  //make simplified list of mods
  vector<mSimpleMod> v;
  size_t b;
  for(size_t a=0; a<v1.mods.size();a++){
    for(b=0;b<v.size();b++){
      if(fabs(v[b].mass-v1.mods[a].mass)<0.001){
        v[b].count++;
        break;
      }
    }
    if(b==v.size()){
      mSimpleMod m;
      m.count=1;
      m.mass=v1.mods[a].mass;
      v.push_back(m);
    }
  }

  //check other mod set for matches
  for(size_t a=0;a<v2.size();a++){
    for(b=0;b<v.size();b++){
      if (fabs(v[b].mass - v2[a].mass)<0.001){
        v[b].count--;
        break;
      }
    }
    if (b == v.size()) return false;
  }

  //verify same number and type of mods
  for(b = 0; b<v.size(); b++){
    if(v[b].count!=0) return false;
  }
  return true;
}

void MSpectrum::shortResults(std::vector<mScoreCard2>& v){
  //cout << "shortResults: ";
  v.clear();
  int count=0;
  int scoreIndex=0;
  mScoreCard tmpSC = getScoreCard(scoreIndex);

  //iterate over all results
  while(tmpSC.simpleScore>0){
    count++;

    //check if this peptide seen in another configuration
    size_t a;
    for(a=0;a<v.size();a++){
      if(v[a].eVal==tmpSC.eVal && v[a].simpleScore==tmpSC.simpleScore && v[a].pep==tmpSC.pep && v[a].mods.size()==tmpSC.mods->size()){
        //size_t b;
        //for(b=0;b<v[a].mods.size();b++){
        //  if(v[a].mods[b].mass!=tmpSC.mods->at(b).mass) break;
        //  if (v[a].mods[b].pos != tmpSC.mods->at(b).pos) break;
        //  if (v[a].mods[b].term != tmpSC.mods->at(b).term) break;
        //}
        //if(b==v[a].mods.size()){ //identical peptides
        size_t b;
        for(b=0;b<v[a].sites.size();b++){
          if(tmpSC.site==v[a].sites[b]) break;
        }
        if(b==v[a].sites.size()) v[a].sites.push_back(tmpSC.site);
        /*} else {
          mScoreCard2 sc;
          sc.conFrag=tmpSC.conFrag;
          sc.eVal=tmpSC.eVal;
          sc.mass=tmpSC.mass;
          sc.massA=tmpSC.massA;
          sc.match=tmpSC.match;
          sc.pep=tmpSC.pep;
          sc.precursor=tmpSC.precursor;
          sc.simpleScore=tmpSC.simpleScore;
          sc.sites.push_back(tmpSC.site);
          for(size_t c=0; c<tmpSC.mods->size();c++) sc.mods.push_back(tmpSC.mods->at(c));
          v.push_back(sc);
        }*/
        break;
      }
    }

    //add unique results
    if (a==v.size()){
      mScoreCard2 sc;
      sc.conFrag = tmpSC.conFrag;
      sc.eVal = tmpSC.eVal;
      sc.mass = tmpSC.mass;
      sc.massA = tmpSC.massA;
      sc.match = tmpSC.match;
      sc.pep = tmpSC.pep;
      sc.precursor = tmpSC.precursor;
      sc.simpleScore = tmpSC.simpleScore;
      sc.sites.push_back(tmpSC.site);
      for (size_t c = 0; c<tmpSC.mods->size(); c++) sc.mods.push_back(tmpSC.mods->at(c));
      v.push_back(sc);
    }

    if(++scoreIndex<20) tmpSC = getScoreCard(scoreIndex);
    else break;
  }

  //cout << count << endl;

}

void MSpectrum::shortResults2(std::vector<mScoreCard3>& v){
  //cout << "shortResults: ";
  v.clear();
  int count = 0;
  int scoreIndex = 0;
  mScoreCard tmpSC = getScoreCard(scoreIndex);

  //iterate over all results
  while (tmpSC.simpleScore>0){
    count++;

    //check if this peptide seen in another configuration
    size_t a;
    for (a = 0; a<v.size(); a++){
      if (v[a].eVal == tmpSC.eVal && v[a].simpleScore == tmpSC.simpleScore && v[a].pep == tmpSC.pep && v[a].modCount == (int)tmpSC.mods->size()){
        if(matchMods(v[a].mSet[0],*tmpSC.mods)){ //mods must also be of the same type, not just same number
          //alternate version of existing peptide, so append the novel mod positions
          mPepMod2 pm;
          for (size_t c = 0; c<tmpSC.mods->size(); c++) {
            pm.mods.push_back(tmpSC.mods->at(c));
          }
          v[a].mSet.push_back(pm);
          v[a].aSites.push_back(tmpSC.site);
          break;
        }
      }
    }

    //add unique results
    if (a == v.size()){
      mScoreCard3 sc;
      sc.conFrag = tmpSC.conFrag;
      sc.eVal = tmpSC.eVal;
      sc.mass = tmpSC.mass;
      sc.massA = tmpSC.massA;
      sc.match = tmpSC.match;
      sc.pep = tmpSC.pep;
      sc.precursor = tmpSC.precursor;
      sc.simpleScore = tmpSC.simpleScore;
      sc.modCount=(int)tmpSC.mods->size();
      sc.aSites.push_back(tmpSC.site);
      mPepMod2 pm;
      for (size_t c = 0; c<tmpSC.mods->size(); c++) {
        pm.mods.push_back(tmpSC.mods->at(c));
      }
      sc.mSet.push_back(pm);
      v.push_back(sc);
    }

    if (++scoreIndex<20) tmpSC = getScoreCard(scoreIndex);
    else break;
  }

}

//void MSpectrum::makeXCorrB(int index, double massOffset, double maxMass, int maxZ, int maxIon){
//  double xcorr=0;
//  double ion;
//  double m;
//  int key,pos;
//
//  for (int a = 0; a<maxIon; a++) { // iterate through required # of ions
//    ion = decoys->decoyIons[index].pdIonsN[a]+massOffset;
//    if (ion>maxMass) break; //stop when fragment ion masses exceed precursor mass
//
//    for (int z = 1; z<maxZ; z++) {
//      m = (ion + (z - 1)*1.007276466) / z;
//      m = binSize * (int)(m*invBinSize + binOffset);
//      key = (int)m;
//      if (key >= kojakBins) break;
//      if (kojakSparseArray[key] == NULL) continue;
//      pos = (int)((m - key)*invBinSize);
//      xcorr += kojakSparseArray[key][pos];
//    }
//  }
//
//
//
//  //iterate from 0 to maxlen
//    //compute yn;
//    //if n>minlen/2, function(n+1 to n*2,+offset)
//    //if n>minlen, function(b+offset);
//    //compute bn;
//    //if n>minlex, function(y+offset);
//
//  //function(n + 1 to n * 2, +offset)
//    //function( b+offset, upto n);
//    //compute yn+offset
//    //compute bn;  <-this is redundant for all lengths
//}


/*============================
  Utilities
============================*/
int MSpectrum::compareIntensity(const void *p1, const void *p2){
  const mSpecPoint d1 = *(mSpecPoint *)p1;
  const mSpecPoint d2 = *(mSpecPoint *)p2;
  if(d1.intensity<d2.intensity) return -1;
  else if(d1.intensity>d2.intensity) return 1;
  else return 0;
}

int MSpectrum::compareMZ(const void *p1, const void *p2){
  const mSpecPoint d1 = *(mSpecPoint *)p1;
  const mSpecPoint d2 = *(mSpecPoint *)p2;
  if(d1.mass<d2.mass) return -1;
  else if(d1.mass>d2.mass) return 1;
  else return 0;
}


//** temporary
//void MSpectrum::tHistogram(double score, int len){
//  int k = (int)(score * 10 + 0.5);
//  if (k < 0) k = 0;
//  else if (k >= HISTOSZ) k = HISTOSZ - 1;
//  hX[len][k] += 1;
//  hXCount[len]++;
//}
//
//void MSpectrum::exportHisto(){
//
//  for(int k=10;k<30;k++){
//    for (int j = 0; j<HISTOSZ; j++) histogram[j] = hX[k][j];
//    histogramCount = hXCount[k];
//
//    MHistogram mHistoX;
//    linearRegression3(mHistoX.slope, mHistoX.intercept, mHistoX.rSq);
//    cout << "\nReal Histogram: " << k << "\t" << histogramCount << endl;
//    for (int i = 0; i<HISTOSZ; i++) cout << i << "\t" << histogram[i] << endl;
//    cout << mHistoX.slope << "\t" << mHistoX.intercept << "\t" << mHistoX.rSq << endl;;
//  }
//}
//**
