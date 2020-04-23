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
MSpectrum::MSpectrum(const int& i, const double& bs, const double& os){
  binOffset=os;
  binSize=bs;
  instrumentPrecursor=false;
  invBinSize=1.0/binSize;
  charge = 0; 
  maxIntensity=0;
  mz = 0;
  precursor = new vector<mPrecursor>;
  singlets = new vector<MTopPeps>;
  spec = new vector<mSpecPoint>;
  scanNumber = 0;
  rTime = 0;
  xCorrArraySize=0;
  xCorrSparseArraySize=0;
  xCorrSparseArray=NULL;
  
  singletCount=0;
  singletFirst=NULL;
  singletLast=NULL;
  singletMax=i;

  kojakSparseArray=NULL;
  kojakBins=0;

  //singletList=NULL;
  //singletBins=0;

  lowScore=0;

  for(int j=0;j<HISTOSZ;j++) histogram[j]=0;
  histogramCount=0;
  histoMaxIndex=0;

  cc=0;
  sc=0;
}

MSpectrum::MSpectrum(const MSpectrum& p){
  unsigned int i;
  int j;
  spec = new vector<mSpecPoint>;
  for(i=0;i<p.spec->size();i++) spec->push_back(p.spec->at(i));
  for(i=0;i<20;i++) topHit[i]=p.topHit[i];
  precursor = new vector<mPrecursor>;
  for(i=0;i<p.precursor->size();i++) precursor->push_back(p.precursor->at(i));
  singlets = new vector<MTopPeps>;
  for (i = 0; i<p.singlets->size(); i++) singlets->push_back(p.singlets->at(i));

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

  cc=p.cc;
  sc=p.sc;

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
}


/*============================
  Operators
============================*/
MSpectrum& MSpectrum::operator=(const MSpectrum& p){
  if(this!=&p){
    unsigned int i;
    int j;
    delete spec;
    spec = new vector<mSpecPoint>;
    for(i=0;i<p.spec->size();i++) spec->push_back(p.spec->at(i));
    for(i=0;i<20;i++) topHit[i]=p.topHit[i];
    delete precursor;
    precursor = new vector<mPrecursor>;
    for(i=0;i<p.precursor->size();i++) precursor->push_back(p.precursor->at(i));
    delete singlets;
    singlets = new vector<MTopPeps>;
    for (i = 0; i<p.singlets->size(); i++) singlets->push_back(p.singlets->at(i));

    binOffset = p.binOffset;
    binSize = p.binSize;
    instrumentPrecursor = p.instrumentPrecursor;
    invBinSize = p.invBinSize;
    charge = p.charge;
    maxIntensity = p.maxIntensity;
    mz = p.charge;
    scanNumber = p.scanNumber;
    rTime = p.rTime;
    xCorrArraySize = p.xCorrArraySize;
    xCorrSparseArraySize = p.xCorrSparseArraySize;
    lowScore = p.lowScore;

    for (i = 0; i<HISTOSZ; i++) histogram[i] = p.histogram[i];
    histogramCount = p.histogramCount;
    histoMaxIndex = p.histoMaxIndex;

    cc = p.cc;
    sc = p.sc;

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
      for(j=0;j<xCorrSparseArraySize;j++) xCorrSparseArray[j]=p.xCorrSparseArray[j];
    }
    
    if(kojakSparseArray!=NULL){
      for(j=0;j<kojakBins;j++){
        if(kojakSparseArray[j]!=NULL) delete [] kojakSparseArray[j];
      }
      delete [] kojakSparseArray;
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
}

void MSpectrum::clear(){
  spec->clear();
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

void MSpectrum::checkScore(mScoreCard& s){
  unsigned int i;
  unsigned int j;

  //char str[256];
  //sprintf(str,"%d.txt",scanNumber);
  //FILE* f=fopen(str,"at");
  //fprintf(f,"%.4f\n",s.simpleScore);
  //fclose(f);

  int index;
  index = (int)(s.simpleScore * 10.0 + 0.5);
  if (index >= HISTOSZ) index = HISTOSZ - 1;
  histogram[index]++;
  histogramCount++;

  //edge case for "reversible" cross-links: check if already matches top hit identically
  //note that such duplications still occur below the top score, but shouldn't influence the final result to the user
  int k=0;
  while(k<20 && s.simpleScore==topHit[k].simpleScore){
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
    if(s.simpleScore > topHit[i].simpleScore) {
      for(j=19;j>i;j--) {
        topHit[j]=topHit[j-1];
      }
      topHit[i] = s;
      lowScore=topHit[19].simpleScore;
      return;
    }
  }
}

void MSpectrum::checkSingletScore(mScoreCard& s){

  mScoreCard* sc;
  mScoreCard* cur;
  
  //If list is empty, add the score card
  if(singletCount==0){
    singletFirst=new mScoreCard(s);
    singletLast=singletFirst;
    singletCount++;
    return;
  }  

  //check if we can just add to the end
  if(s.simpleScore<singletLast->simpleScore){
    //check if we need to store the singlet
    if(singletCount==singletMax) return;

    singletLast->next=new mScoreCard(s);
    singletLast->next->prev=singletLast;
    singletLast=singletLast->next;
    singletCount++;
    return;
  }
  
  //check if it goes in the front
  if(s.simpleScore>=singletFirst->simpleScore){

    singletFirst->prev=new mScoreCard(s);
    singletFirst->prev->next=singletFirst;
    singletFirst=singletFirst->prev;

    //add to singlet list
    if(singletCount<singletMax) {
      singletCount++;
    } else {
      cur=singletLast;
      singletLast=singletLast->prev;
      singletLast->next=NULL;
      delete cur;
    }
    return;
  }


  //scan to find insertion point
  cur = singletFirst->next;
  int i=1;
  while(s.simpleScore < cur->simpleScore){
    i++;
    cur=cur->next;
  }

  sc=new mScoreCard(s);
  sc->prev=cur->prev;
  sc->next=cur;
  cur->prev->next=sc;
  cur->prev=sc;
  if(sc->prev==NULL) singletFirst=sc;

  if(singletCount<singletMax) {
    singletCount++;
  } else {
    cur=singletLast;
    singletLast=singletLast->prev;
    singletLast->next=NULL;
    delete cur;
  }

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

void MSpectrum::xCorrScore(){
  kojakXCorr();
}


/*============================
  Private Functions
============================*/
void MSpectrum::kojakXCorr(){
  int i;
  int j;
  int iTmp;
  double dTmp;
  double dSum;
  double *pdTempRawData;
  double *pdTmpFastXcorrData;
  float  *pfFastXcorrData;
  mPreprocessStruct pPre;

  pPre.iHighestIon = 0;
  pPre.dHighestIntensity = 0;

  BinIons(&pPre);
  //cout << scanNumber << ": " << kojakBins << "\t" << xCorrArraySize << "\t" << invBinSize << "\t" << (int)invBinSize+1 << endl;
  kojakSparseArray=new char*[kojakBins];
  for(i=0;i<kojakBins;i++) kojakSparseArray[i]=NULL;

  pdTempRawData = (double *)calloc((size_t)xCorrArraySize, (size_t)sizeof(double));
  if (pdTempRawData == NULL) {
    fprintf(stderr, " Error - calloc(pdTempRawData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  pdTmpFastXcorrData = (double *)calloc((size_t)xCorrArraySize, (size_t)sizeof(double));
  if (pdTmpFastXcorrData == NULL) {
    fprintf(stderr, " Error - calloc(pdTmpFastXcorrData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  pfFastXcorrData = (float *)calloc((size_t)xCorrArraySize, (size_t)sizeof(float));
  if (pfFastXcorrData == NULL) {
    fprintf(stderr, " Error - calloc(pfFastXcorrData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  // Create data for correlation analysis.
  MakeCorrData(pdTempRawData, &pPre, 50.0);

  // Make fast xcorr spectrum.
  dSum=0.0;
  for (i=0; i<75; i++) dSum += pPre.pdCorrelationData[i].intensity;
  for (i=75; i < xCorrArraySize +75; i++) {
    if (i<xCorrArraySize) dSum += pPre.pdCorrelationData[i].intensity;
    if (i>=151) dSum -= pPre.pdCorrelationData[i-151].intensity;
    pdTmpFastXcorrData[i-75] = (dSum - pPre.pdCorrelationData[i-75].intensity)* 0.0066666667;
  }

  xCorrSparseArraySize=1;
  for (i=0; i<xCorrArraySize; i++) {
    dTmp = pPre.pdCorrelationData[i].intensity - pdTmpFastXcorrData[i];
    pfFastXcorrData[i] = (float)dTmp;

    // Add flanking peaks if used
    iTmp = i-1;
    if (iTmp >= 0) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp])*0.5);

    iTmp = i+1;
    if (iTmp < xCorrArraySize) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp])*0.5);

  }
  free(pdTmpFastXcorrData);

  //MH: Fill sparse matrix
  for(i=0;i<xCorrArraySize;i++){
    if(pfFastXcorrData[i]>0.5 || pfFastXcorrData[i]<-0.5){

      //Fill in missing masses as a result of adding flanking peaks
      if(pPre.pdCorrelationData[i].mass==0){
        j=1;
        while(true){
          if( (i+j)<xCorrArraySize){
            if(pPre.pdCorrelationData[i+j].mass>0){
              pPre.pdCorrelationData[i].mass=pPre.pdCorrelationData[i+j].mass-j*binSize;
              break;
            }
          }
          if( (i-j)>-1){
            if(pPre.pdCorrelationData[i-j].mass>0){
              pPre.pdCorrelationData[i].mass=pPre.pdCorrelationData[i-j].mass+j*binSize;
              break;
            }
          }
          j++;
        }
      }

      //convert i to sparse array key
      //dTmp=pPre.pdCorrelationData[i].mass+binSize*binOffset;
      dTmp=binSize*i;
      iTmp=(int)dTmp;
      //cout << i << "\t" << pfFastXcorrData[i] << "\t" << dTmp << "\t" << iTmp << endl;
      if(kojakSparseArray[iTmp]==NULL) {
        kojakSparseArray[iTmp]=new char[(int)invBinSize+1];
        for(j=0;j<(int)invBinSize+1;j++) kojakSparseArray[iTmp][j]=0;
      }
      j=(int)((dTmp-iTmp)*invBinSize/*+0.5*/);
      //cout << (dTmp-iTmp) << "\t" << (dTmp-iTmp)*invBinSize/*+0.5*/ << endl;
      //cout << j << endl;
      //if( j>(int)invBinSize) {
      //  cout << "ERROR!" << endl;
      //  exit(0);
      //}
      if(pfFastXcorrData[i]>127) kojakSparseArray[iTmp][j]=127;
      else if(pfFastXcorrData[i]<-128) kojakSparseArray[iTmp][j]=-128;
      else if(pfFastXcorrData[i]>0) kojakSparseArray[iTmp][j]=(char)(pfFastXcorrData[i]+0.5);
      else kojakSparseArray[iTmp][j]=(char)(pfFastXcorrData[i]-0.5);
      //cout << i << "\t" << iTmp << "\t" << j << "\t" << (int)kojakSparseArray[iTmp][j] << endl;
    }
  }

  /*
  if(scanNumber==11368){
    for(i=0;i<kojakBins;i++){
      if(kojakSparseArray[i]==NULL) {
        cout << i << "\tNULL" << endl;
        continue;
      }
      for(j=0;j<(int)invBinSize+1;j++){
        cout << i << "\t" << j << "\t" << (int)kojakSparseArray[i][j] << endl;
      }
    }
  }
  */

  //exit(1);
  free(pPre.pdCorrelationData);
  free(pfFastXcorrData);
  free(pdTempRawData);

}

void MSpectrum::BinIons(mPreprocessStruct *pPre) {
  int i;
  unsigned int j;
  double dPrecursor;
  double dIon;
  double dIntensity;

  // Just need to pad iArraySize by 75.
  dPrecursor=0;
  for(j=0;j<precursor->size();j++){
    if(precursor->at(j).monoMass>dPrecursor) dPrecursor=precursor->at(j).monoMass;
  }
  xCorrArraySize = (int)((dPrecursor + 100.0) / binSize);
  kojakBins = (int)(spec->at(spec->size()-1).mass+100.0);

  pPre->pdCorrelationData = (mSpecPoint *)calloc(xCorrArraySize, (size_t)sizeof(mSpecPoint));
  if (pPre->pdCorrelationData == NULL) {
    fprintf(stderr, " Error - calloc(pdCorrelationData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  i = 0;
  while(true) {
    if (i >= (int)spec->size()) break;

    dIon = spec->at(i).mass;
    dIntensity = spec->at(i).intensity;   
    i++;

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

  if (histogramCount < DECOY_SIZE) {
    if (!generateXcorrDecoys(params,decoys)) return false;
  }

  linearRegression(dSlope, dIntercept, iMaxCorr, iStartCorr, iNextCorr);
  histoMaxIndex = iMaxCorr;

  dSlope *= 10.0;

  iLoopCount=20; //score all e-values among top hits?
  for (i = 0; i<iLoopCount; i++) {
    if(topHit[i].simpleScore==0) break; //score all e-values among top hits?
    if (dSlope >= 0.0) {
      topHit[i].eVal=999.0;
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
  int r;
  int key;
  int pos;
  double dBion;
  double dYion;
  double dXcorr;
  double dFragmentIonMass = 0.0;

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = DECOY_SIZE - histogramCount;
  int seed = rand()%3000;
  int decoyIndex;

  j = 0;
  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex=(seed+i)%3000;

    //grab precursor at random
    r=rand()%precursor->size();
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
  }

  return true;
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
