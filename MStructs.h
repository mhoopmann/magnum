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

#define MAX_PRECURSOR 32

//FASTA database structure
typedef struct mDB{
  bool decoy;
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

typedef struct mAtom {
  int index = 0;
  int count = 0;
}mAtom;

typedef struct mParams {
  int     atomicProcessing = 0;
  int     eValDepth = 5000;
  int     instrument = 0;     //0=Orbi, 1=FTICR
  int     isotopeError = 3;
  int     maxMods = 2;
  int     maxPeaks = 0;
  int     maxPepLen = 50;
  int     minPeaks = 20;
  int     minPepLen = 6;
  int     miscleave = 2;
  int     ms1Centroid = 1;
  int     ms2Centroid = 1;
  int     ms1Resolution = 60000;
  int     ms2Resolution = 15000;
  int     preferPrecursor = 2;
  int     setA = 0;
  int     setB = 0;
  int     specProcess = 1;
  int     threads = 1;
  int     topCount = 5;
  int     truncate = 0;
  bool    buildDecoy = false;
  bool    exportPepXML = true;
  bool    exportPercolator = false;
  bool    ionSeries[6] = { false,true,false,false,true,false };
  bool    precursorRefinement = true;
  bool    splitPercolator = false;
  bool    xcorr = false;
  double  binOffset = 0.0;
  double  binSize = 0.03;
  double  maxPepMass = 4000.0;
  double  minPepMass = 500.0;
  double  maxAdductMass = 500.0;
  double  minAdductMass = 10.0;
  double  percVersion = 2.04;
  double  ppmPrecursor = 25.0;
  double  rIonThreshold = 10;
  std::string     adductSites;
  std::string     dbFile;
  std::string     decoy = "random";
  std::string     enzyme = "[KR]|{P}";
  std::string     enzymeName = "Trypsin";
  std::string     ext;
  std::string     inFile;  //true input file with full path
  std::string     inFileNoExt;  //for pepXML
  std::string     msFile;   //input file parameter from config
  std::string     outFile;  //true output file with full path
  std::string     resPath;
  std::string     dbPath;
  std::string     msBase;
  std::vector<mMass>    aaMass;
  std::vector<mAtom>    atomSig;
  std::vector<int>      diag;
  std::vector<mMass>    mods;
  std::vector<mMass>    fMods;
  std::vector<double>   rIons;
} mParams;

typedef struct mSpecPoint{
  double mass;
  float intensity;
} mSpecPoint;

typedef struct mPreprocessStruct { //adapted from Comet
  int iMaxXCorrArraySize;
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
  char    type;  //0=selected peak, 1=instrument predicted, 2=hardklor predicted
  char    offset; //offset amount (rounded)
  mPrecursor(){
    charge=0;
    corr=0;
    label=0;
    monoMass=0;
    type=0;
    offset=0;
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
  double  openModMass;
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
    openModMass=0;
    monoMass=0;
    ppm=0;
    scoreDelta=0;
    scoreMagnum=0;
    selectedMZ=0;
  }
} mResults;

typedef struct mSimpleMod{
  double mass;
  int count;
} mSimpleMod;


typedef struct sLink2{
  size_t nextNode=SIZE_MAX;
  size_t nextIndex = SIZE_MAX;
  size_t pepNum=0;
  size_t parentPepNum=0;
  double score[MAX_PRECURSOR]{0};    //make an array large enough to hold all precursors? Is 10 enough?
  double scoreNL[MAX_PRECURSOR]{0};  //make an array large enough to hold all precursors? Is 10 enough?
} sLink2;


typedef struct sNode2{
  bool visit=false;
  double mass=0;
  double lastMass=0;
  double score[4]{};       //always same? is it necessary to track charge states?
  double scoreAlt[MAX_PRECURSOR][4]{0};    //make an array large enough to hold all precursors? Is 10 enough?
  double scoreAltNL[4]{0};  //always same?
  //int match[4]{};
  //int matchAlt[10][4]{};
  //int matchAltNL[4]{};
  std::vector<sLink2> start;
  std::vector<sLink2> next;
  //sNode2(){
  //  visit = false;
  //  mass = 0;
  //}
  int id=0; //for diagnostics only
} sNode2;

typedef struct sPepModSet{
  std::vector<mPepMod> mods;
  mPepMod& operator[](const size_t& index) { return mods[index]; }
} sPepModSet;

typedef struct sPrecursor{
  int index;
  double monomass;
  int maxZ;
} sPrecursor;

typedef struct sScoreSet {
  double scores[MAX_PRECURSOR]{};  //never more than 10 precursors
  //int match[10]{};
  bool scored=false;
  //sScoreSet(){
  //  scored = false;
  //}
} sScoreSet;

typedef struct sDIndex{
  size_t a;
  size_t b;
} sDIndex;

typedef struct sIPep {
  double pepMass=0;
  std::vector<size_t> pepIndex;
} sIPep;

typedef struct sIPeak {
  double mass=0;
  std::vector<sIPep> index;
} sIPeak;

typedef struct sScoreSet2 {
  double scoreP[MAX_PRECURSOR]{};  //never more than 10 precursors
  double score = 0;
} sScoreSet2;

#endif
