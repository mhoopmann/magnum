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

#ifndef _MSTRUCTS_H
#define _MSTRUCTS_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <iostream>
#include <exception>

//FASTA database structure
typedef struct mDB{
  std::string description; //FASTA description (header - after first space)
  std::string name;        //FASTA name (header - before first space)
  std::string sequence;    //FASTA sequence
} mDB;

typedef struct mFile{
  std::string input;
  std::string base;
  std::string ext;
} kFile;

//structure holds peptide mappings to database
typedef struct mPepMap{
  int index;  //protein index
  unsigned short start;  //first aa
  unsigned short stop;   //last aa
} mPepMap;

//Peptide reference to an entry in pldbDB
typedef struct mPeptide{
  bool cTerm;
  bool nTerm;
  char xlSites;
  double mass;            //monoisotopic, zero mass
  std::vector<mPepMap>* map;   //array of mappings where peptides appear in more than one place
  mPeptide(){
    cTerm=false;
    nTerm=false;
    xlSites=0;
    mass=0;
    map = new std::vector<mPepMap>;
  }
  mPeptide(const mPeptide& m){
    cTerm=m.cTerm;
    nTerm=m.nTerm;
    xlSites=m.xlSites;
    mass=m.mass;
    map = new std::vector<mPepMap>(*m.map);
  }
  ~mPeptide(){
    delete map;
  }
  mPeptide& operator=(const mPeptide& m){
    if(this!=&m){
      cTerm = m.cTerm;
      nTerm = m.nTerm;
      xlSites = m.xlSites;
      mass=m.mass;
      delete map;
      map = new std::vector<mPepMap>(*m.map);
    }
    return (*this);
  }
} mPeptide;

typedef struct kPeptideB{
  double  mass;
  int     index;
  bool    linkable;
} kPeptideB;

//For sorting peptide lists
typedef struct mPepSort{
  int index;        //peptide array index
  std::string sequence;  //peptide sequence
} mPepSort;

typedef struct mMass {
  bool    xl;
  int     index;
  double  mass;
} mMass;

typedef struct mSparseMatrix{
  int   bin;
  float fIntensity;
} mSparseMatrix;

typedef struct mEnzymeRules{
  bool cutC[128];
  bool cutN[128];
  bool exceptC[128];
  bool exceptN[128];
} mEnzymeRules;

typedef struct mParams {
  int     eValDepth;
  int     instrument;     //0=Orbi, 1=FTICR
  int     isotopeError;
  int     maxMods;
  int     maxPeaks;
  int     maxPepLen;
  int     minPeaks;
  int     minPepLen;
  int     miscleave;
  int     ms1Centroid;
  int     ms2Centroid;
  int     ms1Resolution;
  int     ms2Resolution;
  int     preferPrecursor;
  int     setA;
  int     setB;
  int     specProcess;
  int     threads;
  int     topCount;
  int     truncate;
  bool    exportPepXML;
  bool    exportPercolator;
  bool    ionSeries[6];
  bool    precursorRefinement;
  bool    xcorr;
  double  binOffset;
  double  binSize;
  double  maxPepMass;
  double  minPepMass;
  double  maxAdductMass;
  double  minAdductMass;
  double  percVersion;
  double  ppmPrecursor;
  double  rIonThreshold;
  std::string    adductSites;
  std::string     dbFile;
  std::string    decoy;
  std::string     enzyme;
  std::string     enzymeName;
  std::string     ext;
  std::string     inFile;  //true input file with full path
  std::string     msFile;   //input file parameter from confic
  std::string     outFile;  //true output file with full path
  std::string     resPath;
  std::string     dbPath;
  std::string     msBase;
  std::vector<mMass>    aaMass;
  std::vector<int>      diag;
  std::vector<mMass>    mods;
  std::vector<mMass>    fMods;
  std::vector<double>   rIons;
  mParams(){
    eValDepth=3000;
    instrument=1;
    isotopeError=1;
    maxMods=0;
    maxPeaks=0;
    maxPepLen=50;
    minPeaks=20;
    minPepLen=6;
    miscleave=2;
    ms1Centroid=0;
    ms2Centroid=0;
    ms1Resolution=60000;
    ms2Resolution=15000;
    preferPrecursor=1;
    setA=0;
    setB=0;
    specProcess=0;
    threads=1;
    topCount=250;
    truncate=0;
    exportPepXML=false;
    exportPercolator=false;
    ionSeries[0]=false; //a-ions
    ionSeries[1]=true;  //b-ions
    ionSeries[2]=false; //c-ions
    ionSeries[3]=false; //x-ions
    ionSeries[4]=true;  //y-ions
    ionSeries[5]=false; //z-ions
    precursorRefinement=true;
    xcorr=false;
    binSize=0.03;
    binOffset=0.0;
    maxPepMass=4000.0;
    minPepMass=500.0;
    maxAdductMass=500.0;
    minAdductMass=10.0;
    percVersion=2.04;
    ppmPrecursor=25.0;
    rIonThreshold=10;
    decoy="random";
    enzyme="[KR]|{P}";
    enzymeName="Trypsin";
  }
} mParams;

typedef struct mSpecPoint{
  double mass;
  float intensity;
} mSpecPoint;

typedef struct mPreprocessStruct { //adapted from Comet
   int iHighestIon;
   double dHighestIntensity;
   mSpecPoint *pdCorrelationData;
} mPreprocessStruct;

typedef struct mPepMod{
  bool term;
  char pos;
  double mass;
} mPepMod;

typedef struct mScoreCard{
  char    precursor;
  char    site;
  int     conFrag;
  int     match;
  int     pep;
  float   simpleScore;
  double  mass;
  double  massA;
  double  eVal;
  std::vector<mPepMod>* mods;
  mScoreCard*  next;
  mScoreCard*  prev;
  mScoreCard(){
    precursor=0;
    site=0;
    match=0;
    conFrag=0;
    pep=0;
    simpleScore=0;
    mass=0;
    massA=0;
    eVal=1000;
    mods = new std::vector<mPepMod>;
    next = NULL;
    prev = NULL;
  }
  mScoreCard(const mScoreCard& p){
    precursor=p.precursor;
    site=p.site;
    match=p.match;
    conFrag=p.conFrag;
    pep=p.pep;
    simpleScore=p.simpleScore;
    mass=p.mass;
    massA=p.massA;
    eVal=p.eVal;
    mods = new std::vector<mPepMod>(*p.mods);
    next = NULL;
    prev = NULL;
  }
  ~mScoreCard(){
    delete mods;
  }
  mScoreCard& operator=(const mScoreCard& p){
    if(this!=&p){
      precursor=p.precursor;
      site=p.site;
      match = p.match;
      conFrag = p.conFrag;
      pep=p.pep;
      simpleScore=p.simpleScore;
      mass=p.mass;
      massA=p.massA;
      eVal=p.eVal;
      delete mods;
      mods = new std::vector<mPepMod>(*p.mods);
      next = NULL;
      prev = NULL;
    }
    return *this;
  }
} mScoreCard;

typedef struct mScoreCard2{
  char    precursor;
  int     conFrag;
  int     match;
  int     pep;
  float   simpleScore;
  double  mass;
  double  massA;
  double  eVal;
  std::vector<mPepMod> mods;
  std::vector<char> sites;
} mScoreCard2;

typedef struct mPepMod2{
  std::vector<mPepMod> mods;
} mPepMod2;

typedef struct mScoreCard3{
  char    precursor;
  int     conFrag;
  int     match;
  int     pep;
  int     modCount;
  float   simpleScore;
  double  mass;
  double  massA;
  double  eVal;
  std::vector<mPepMod2> mSet;
  std::vector<char> aSites;
} mScoreCard3;



typedef struct mPrecursor{ 
  int     charge;
  double  corr;
  char    label;
  double  monoMass;
  mPrecursor(){
    charge=0;
    corr=0;
    label=0;
    monoMass=0;
  }
} mPrecursor;

typedef struct rMods{
  double mass;
  char pos;
  bool term;
  bool variable;
  bool adduct;
} rMods;

typedef struct rMods2{
  std::vector<rMods> mods;
} rMods2;

typedef struct mProtRes{
  std::string protein;
  char prevAA;
  char nextAA;
  int startPos;
} mProtRes;

typedef struct mResults{
  bool    decoy;
  int     charge;
  int     psmID;
  int     scanNumber;
  float   rTimeSec;
  double  eValue;
  double  monoMass;
  double  ppm;
  double  psmMass;
  double  scoreDelta;
  double  scoreMagnum;
  double  selectedMZ;
  std::vector<std::string>  modPeptide;
  std::string  mods;
  std::string  openMod;
  std::string  peptide;
  std::string  percPeptide;
  std::vector<rMods2> vMods;
  std::vector<mProtRes> proteins;
  std::vector<double> rIon;
  mResults(){
    decoy=false;
    charge=0;
    psmID=0;
    scanNumber=0;
    rTimeSec=0;
    eValue=0;
    monoMass=0;
    ppm=0;
    scoreDelta=0;
    scoreMagnum=0;
    selectedMZ=0;
  }
} mResults;

typedef struct kResults{
  bool    decoy;
  bool    linkable1;
  bool    linkable2;
  bool    cTerm1; //peptide contains protein c-terminus
  bool    nTerm1; //peptide contains protein n-terminus
  bool    cTerm2; //peptide contains protein c-terminus
  bool    nTerm2; //peptide contains protein n-terminus
  bool    n15Pep1;
  bool    n15Pep2;
  int     charge;
  int     link1;
  int     link2;
  int     pep1;
  int     pep2;
  int     rank;
  int     rankA;
  int     rankB;
  int     scanNumber;
  int     type;
  float   rTime;
  double  eVal;
  double  hk;
  double  massA;
  double  massB;
  double  obsMass;
  double  ppm;
  double  psmMass;
  double  score;
  double  scoreA;
  double  scoreB;
  double  scoreDelta;
  double  scorePepDif;
  double  xlMass;
  std::string  modPeptide1;
  std::string  modPeptide2;
  std::string  peptide1;
  std::string  peptide2;
  std::string  xlLabel;
  std::string  rIon;
  std::vector<mPepMod> mods1;
  std::vector<mPepMod> mods2;
} kResults;

typedef struct kCKey{
  int key;
  int pos;
} kCKey;

typedef struct kMatchSet{
  int sz;
  int charge;
  kCKey** a;
  kCKey** b;
  kCKey** c;
  kCKey** x;
  kCKey** y;
  kCKey** z;
  kMatchSet(){
    sz=0;
    charge=0;
    a=NULL;
    b=NULL;
    c=NULL;
    x=NULL;
    y=NULL;
    z=NULL;
  }
  ~kMatchSet(){
    deAllocate();
  }
  void allocate(int s, int ch){
    deAllocate();
    sz=s;
    charge=ch;
    //kCKey k;
    //k.key=0;
    //k.pos=0;
    a = new kCKey*[charge];
    b = new kCKey*[charge];
    c = new kCKey*[charge];
    x = new kCKey*[charge];
    y = new kCKey*[charge];
    z = new kCKey*[charge];
    for(int i=0;i<charge;i++){
      a[i] = new kCKey[sz];
      b[i] = new kCKey[sz];
      c[i] = new kCKey[sz];
      x[i] = new kCKey[sz];
      y[i] = new kCKey[sz];
      z[i] = new kCKey[sz];
      //for(int j=0;j<sz;j++){
      //  a[i][j]=k;
      //  b[i][j]=k;
      //  c[i][j]=k;
      //  x[i][j]=k;
      //  y[i][j]=k;
      //  z[i][j]=k;
      //}
    }
  }
  void deAllocate(){
    int i;
    if(a!=NULL){
      for(i=0;i<charge;i++) delete [] a[i];
      delete [] a;
    }
    if(b!=NULL){
      for(i=0;i<charge;i++) delete [] b[i];
      delete [] b;
    }
    if(c!=NULL){
      for(i=0;i<charge;i++) delete [] c[i];
      delete [] c;
    }
    if(x!=NULL){
      for(i=0;i<charge;i++) delete [] x[i];
      delete [] x;
    }
    if(y!=NULL){
      for(i=0;i<charge;i++) delete [] y[i];
      delete [] y;
    }
    if(z!=NULL){
      for(i=0;i<charge;i++) delete [] z[i];
      delete [] z;
    }
    a=NULL;
    b=NULL;
    c=NULL;
    x=NULL;
    y=NULL;
    z=NULL;
  }

} kMatchSet;

typedef struct kXLMotif {
  std::string motif;
  int xlIndex[10];
  int counterMotif[10];
} kXLMotif;

typedef struct kXLTarget{
  size_t linkerID;
  bool target[128];
} kXLTarget;

#endif
