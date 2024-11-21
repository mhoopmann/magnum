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

#include "MAnalysis.h"

using namespace std;

bool*       MAnalysis::bKIonsManager;
MDatabase*  MAnalysis::db;
MIons*      MAnalysis::ions;
double      MAnalysis::maxMass;
double      MAnalysis::minMass;
Mutex       MAnalysis::mutexKIonsManager;
Mutex*      MAnalysis::mutexSpecScore;
Mutex**     MAnalysis::mutexSingletScore;
mParams     MAnalysis::params;
MData*      MAnalysis::spec;
bool*       MAnalysis::adductSites;
bool**      MAnalysis::scanBuffer;

int         MAnalysis::numIonSeries;

int*        MAnalysis::pepMassSize;
double**    MAnalysis::pepMass;
bool**      MAnalysis::pepBin;
int*        MAnalysis::pepBinSize;

int MAnalysis::skipCount;
int MAnalysis::nonSkipCount;

MDecoys MAnalysis::decoys;

double MAnalysis::dummy[10]{};
int MAnalysis::dummyM[10]{};
int* MAnalysis::maxZ2;
size_t* MAnalysis::bufSize2;

//bool MAnalysis::bEcho;
//int MAnalysis::sCounter;
bool MAnalysis::echo;

/*============================
  Constructors & Destructors
============================*/
MAnalysis::MAnalysis(mParams& p, MDatabase* d, MData* dat){
  unsigned int i;
  int j,k;
  
  //bEcho=false;

  //Assign pointers and structures
  params=p;
  db=d;
  spec=dat;
  adductSites = spec->getAdductSites();

  //Do memory allocations and initialization
  bKIonsManager=NULL;
  ions=NULL;
  allocateMemory(params.threads);
  for(j=0;j<params.threads;j++){
    for(i=0;i<params.fMods.size();i++) ions[j].addFixedMod((char)params.fMods[i].index,params.fMods[i].mass);
    for(i=0;i<params.mods.size();i++) ions[j].addMod((char)params.mods[i].index,params.mods[i].xl,params.mods[i].mass);
    for(i=0;i<params.aaMass.size();i++) ions[j].setAAMass((char)params.aaMass[i].index, params.aaMass[i].mass);
    ions[j].setMaxModCount(params.maxMods);
  }

  //Initalize variables
  maxMass = spec->getMaxMass()+0.25;
  minMass = spec->getMinMass()-0.25;

  numIonSeries=0;
  for(i=0;i<6;i++){
    if(params.ionSeries[i]) numIonSeries++;
  }

  //Create mutexes
  Threading::CreateMutex(&mutexKIonsManager);
  mutexSingletScore = new Mutex*[spec->size()];
  mutexSpecScore = new Mutex[spec->size()];
  for(j=0;j<spec->size();j++){
    Threading::CreateMutex(&mutexSpecScore[j]);
    mutexSingletScore[j] = new Mutex[spec->at(j).sizePrecursor()];
    for(k=0;k<spec->at(j).sizePrecursor();k++){
      Threading::CreateMutex(&mutexSingletScore[j][k]);
    }
  }

  //xCorrCount=0;
  echo=false;
}

MAnalysis::~MAnalysis(){
  int i,j;

  //Destroy mutexes
  Threading::DestroyMutex(mutexKIonsManager);
  for(i=0;i<spec->size();i++){
    Threading::DestroyMutex(mutexSpecScore[i]);
    for(j=0;j<spec->at(i).sizePrecursor();j++){
      Threading::DestroyMutex(mutexSingletScore[i][j]);
    }
    delete [] mutexSingletScore[i];
  }
  delete [] mutexSingletScore;
  delete [] mutexSpecScore;

  //Deallocate memory and release pointers
  deallocateMemory(params.threads);
  db=NULL;
  spec=NULL;
  adductSites=NULL;
  
}

//============================
//  Public Functions
//============================
bool MAnalysis::doPeptideAnalysis(){
  size_t i;
  int iPercent;
  int iTmp;
  vector<mPeptide>* p;
  vector<int> index;
  vector<mPepMod> mods;

  //mScoreCard sc;

  ThreadPool<mAnalysisStruct*>* threadPool = new ThreadPool<mAnalysisStruct*>(analyzePeptideProc,params.threads,params.threads,1);

  //Set progress meter
  iPercent=0;
  printf("%2d%%",iPercent);
  fflush(stdout);

  //Set which list of peptides to search (with and without internal lysine)
  p=db->getPeptideList();

  //Iterate the peptide for the first pass
  for(i=0;i<p->size();i++){

    threadPool->WaitForQueuedParams();

    mAnalysisStruct* a = new mAnalysisStruct(&mutexKIonsManager,&p->at(i),(int)i);
    threadPool->Launch(a);

    //Update progress meter
    iTmp=(int)((double)i/p->size()*100);
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }
  }

  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  //clean up memory & release pointers
  delete threadPool;
  threadPool=NULL;
  p=NULL;

  return true;
}

bool MAnalysis::doEValuePrecalc(){
  int i;
  int iPercent;
  int iTmp;

  ThreadPool<MSpectrum*>* threadPool = new ThreadPool<MSpectrum*>(analyzeEValuePrecalcProc, params.threads, params.threads, 1);

  //Set progress meter
  iPercent = 0;
  printf("%2d%%", iPercent);
  fflush(stdout);

  //Iterate the peptide for the first pass
  for (i = 0; i<spec->size(); i++){

    threadPool->WaitForQueuedParams();

    MSpectrum* a = &spec->at(i);
    threadPool->Launch(a);

    //Update progress meter
    iTmp = (int)((double)i / spec->size() * 100);
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }
  }

  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  //clean up memory & release pointers
  delete threadPool;
  threadPool = NULL;

  return true;
}

//============================
//  Private Functions
//============================

//============================
//  Thread-Start Functions
//============================

//These functions fire off when a thread starts. They pass the variables to for
//each thread-specific analysis to the appropriate function.
void MAnalysis::analyzePeptideProc(mAnalysisStruct* s){
  int i;
  Threading::LockMutex(mutexKIonsManager);
  for(i=0;i<params.threads;i++){
    if(!bKIonsManager[i]){
      bKIonsManager[i]=true;
      break;
    }
  }
  Threading::UnlockMutex(mutexKIonsManager);
  if(i==params.threads){
    cout << "Error in KAnalysis::analyzePeptidesProc" << endl;
    exit(-1);
  }
  s->bKIonsMem = &bKIonsManager[i];
  analyzePeptide(s->pep,s->pepIndex,i);
  delete s;
  s=NULL;
}

void MAnalysis::analyzeEValuePrecalcProc(MSpectrum* s){
  s->generateXcorrDecoys3(params.minPepLen, db->getMaxPepLen(s->bigMonoMass),params.eValDepth);
  s = NULL;
}

//============================
//  Analysis Functions
//============================

//Analyzes all single peptides. Also analyzes cross-linked peptides when in full search mode, 
//or stage 1 of relaxed mode analysis
bool MAnalysis::analyzePeptide(mPeptide* p, int pepIndex, int iIndex){

  vector<int> index;
  vector<mPepMod> mods;

  char str[256];
  db->getPeptideSeq(p->map->at(0).index,p->map->at(0).start,p->map->at(0).stop,str);
  //if(strcmp(str,"LVRPEVDVMCTAFHDNEETFLK")==0) {
  if (strcmp(str, "NFPPSQDASGDLYTTSSQLTLPATQCLAGK") == 0) {
    cout << "\n" << str << "\t" << p->mass << "\t" << pepIndex << endl;
    echo=true;
  } else echo=true; //else return true;
  //if(strcmp(str,"EFNAETFTFHADICTLSEK")==0) cout <<"Peptide: " << pepIndex << "." << endl;
  //cout << str << "\t" << p->mass << "\t" << pepIndex << endl;

  //Set the peptide, calc the ions, and score it against the spectra
  int len = (p->map->at(0).stop - p->map->at(0).start) + 1;

  //skip peptide if its smallest possible mass with modifications is more than largest precursor.
  if(p->mass>spec->getMaxMass()+1) {
    if(echo) cout << "Too big. Next" << endl;
    echo=false;
    return true;
  }

  ions[iIndex].setPeptide(&db->at(p->map->at(0).index).sequence[p->map->at(0).start],p->map->at(0).stop-p->map->at(0).start+1,p->mass,p->nTerm,p->cTerm);
  ions[iIndex].buildModIons2(false); //having this here is bad if there are lots of mods and few spectra
  if (echo) {
    cout << "set and build" << endl;

    cout << "Forward: " << endl;
    for (size_t a = 0; a < ions[iIndex].vPeaks.size(); a++) {
      cout << a << "\t" << ions[iIndex].vPeaks[a].mass;
      for (size_t b = 0; b < ions[iIndex].vPeaks[a].index.size(); b++) {
        cout << "," << ions[iIndex].vPeaks[a].index[b].pepMass << "(";
        for (size_t c = 0; c < ions[iIndex].vPeaks[a].index[b].pepIndex.size(); c++) {
          cout << "." << ions[iIndex].vPeaks[a].index[b].pepIndex[c];
        }
        cout << ")";
      }
      cout << endl;
    }
    cout << "Reverse: " << endl;
    for (size_t a = 0; a < ions[iIndex].vPeaksRev.size(); a++) {
      cout << a << "\t" << ions[iIndex].vPeaksRev[a].mass;
      for (size_t b = 0; b < ions[iIndex].vPeaksRev[a].index.size(); b++) {
        cout << "," << ions[iIndex].vPeaksRev[a].index[b].pepMass << "(";
        for (size_t c = 0; c < ions[iIndex].vPeaksRev[a].index[b].pepIndex.size(); c++) {
          cout << "." << ions[iIndex].vPeaksRev[a].index[b].pepIndex[c];
        }
        cout << ")";
      }
      cout << endl;
    }
  }

  //Check peptide without open modifications
  double lastMass=0;
  if(echo) cout << "ions[iIndex].pepCount = " << ions[iIndex].pepCount << endl;
  for(size_t j=0;j<ions[iIndex].pepCount;j++){
    if(echo) cout << j << " lastmass = " << lastMass << endl;
    if(ions[iIndex].pepMass[j]<=lastMass) continue; //skip peptide variants already searched
    if(spec->getBoundaries2(ions[iIndex].pepMass[j],params.ppmPrecursor,index,scanBuffer[iIndex])){
      scoreSpectra2(index, ions[iIndex].pepMass[j], len, pepIndex, iIndex);
    }
    lastMass= ions[iIndex].pepMass[j];
  }
  
  if (p->xlSites == 0) {
    if(echo) cout << "No sites" << endl;
    echo=false;
    return true;
  }

  //Search for open modifications on peptide as well if it has sites where modification can bind
  analyzeSinglets(*p, pepIndex, iIndex);

  //cout << "Done: " << str << endl;

  echo=false;
  return true;
  
  
  //
  //ions[iIndex].buildIons();
  //ions[iIndex].modIonsRec2(0,-1,0,0,false);

  //for(j=0;j<ions[iIndex].size();j++){
  //  bt= spec->getBoundaries2(ions[iIndex][j].mass, params.ppmPrecursor, index, scanBuffer[iIndex]);
  //  if(bt) scoreSpectra(index,j,len,ions[iIndex][j].difMass,pepIndex,-1,-1,-1,-1,iIndex);
  // }

  //if(p->xlSites==0) return true;

  ////Search for open modifications on peptide as well if it has sites where modification can bind
  //analyzeSinglets(*p,pepIndex,iIndex);

  //return true;
}

bool MAnalysis::analyzeSinglets(mPeptide& pep, int index, int iIndex) {

  //get the peptide sequence
  string pepSeq;
  db->getPeptideSeq(pep, pepSeq);
  if (echo) cout << "analyzeSinglets(): " << pepSeq << endl;

  //if(index!=903) return false;

  //Build our peptide
  int len = (pep.map->at(0).stop - pep.map->at(0).start) + 1;
  ions[iIndex].setPeptide(&db->at(pep.map->at(0).index).sequence[pep.map->at(0).start], len, pep.mass, pep.nTerm, pep.cTerm);
  ions[iIndex].buildModIons2(); //It would be more efficient to do this after spec->getBoundaries below. Must use alternative way to compute min and max peptides.
  if (echo) {
    cout << "Forward: " << endl;
    for (size_t a = 0; a < ions[iIndex].vPeaks.size(); a++) {
      cout << a << "\t" << ions[iIndex].vPeaks[a].mass;
      for (size_t b = 0; b < ions[iIndex].vPeaks[a].index.size(); b++) {
        cout << "," << ions[iIndex].vPeaks[a].index[b].pepMass << "(";
        for (size_t c = 0; c < ions[iIndex].vPeaks[a].index[b].pepIndex.size(); c++) {
          cout << "." << ions[iIndex].vPeaks[a].index[b].pepIndex[c];
        }
        cout << ")";
      }
      cout << endl;
    }
    cout << "Reverse: " << endl;
    for (size_t a = 0; a < ions[iIndex].vPeaksRev.size(); a++) {
      cout << a << "\t" << ions[iIndex].vPeaksRev[a].mass;
      for (size_t b = 0; b < ions[iIndex].vPeaksRev[a].index.size(); b++) {
        cout << "," << ions[iIndex].vPeaksRev[a].index[b].pepMass << "(";
        for (size_t c = 0; c < ions[iIndex].vPeaksRev[a].index[b].pepIndex.size(); c++) {
          cout << "." << ions[iIndex].vPeaksRev[a].index[b].pepIndex[c];
        }
        cout << ")";
      }
      cout << endl;
    }
  }
                                
  //cout << "BUILD DONE: " << ions[iIndex].pepCount << endl;

  /*if (index == 103) {
    cout << "Bions:" << endl;
    vector<sNode2>* peaks = ions[iIndex].peaks;
    for (size_t a = 0; a < peaks->size(); a++) {
      cout << a << " " << peaks->at(a).id << "\t" << peaks->at(a).mass << "\t" << peaks->at(a).next.size() << endl;
      for (size_t b = 0; b < peaks->at(a).next.size(); b++) {
        cout << "  N\t" << b << "\t" << peaks->at(a).next[b].nextIndex << "\t" << peaks->at(a).next[b].nextNode << "\t" << peaks->at(a).next[b].pepNum << endl;
      }
      for (size_t b = 0; b < peaks->at(a).start.size(); b++) {
        cout << "  S\t" << b << "\t" << peaks->at(a).start[b].nextIndex << "\t" << peaks->at(a).start[b].nextNode << "\t" << peaks->at(a).start[b].pepNum << endl;
      }
    }

    cout << "Yions:"<<endl;
    vector<sNode2>* peaksRev=ions[iIndex].peaksRev;
    for (size_t a = 0; a < peaksRev->size(); a++) {
      cout << a << " " << peaksRev->at(a).id << "\t" << peaksRev->at(a).mass << "\t" << peaksRev->at(a).next.size() << endl;
      for (size_t b = 0; b < peaksRev->at(a).next.size(); b++) {
        cout << "  N\t" << b << "\tNindex: " << peaksRev->at(a).next[b].nextIndex << "\tNnode: " << peaksRev->at(a).next[b].nextNode << "\tPep: " << peaksRev->at(a).next[b].pepNum << endl;
      }
      for (size_t b = 0; b < peaksRev->at(a).start.size(); b++) {
        cout << "  S\t" << b << "\tSindex: " << peaksRev->at(a).start[b].nextIndex << "\tSnode: " << peaksRev->at(a).start[b].nextNode << "\tPep: " << peaksRev->at(a).start[b].pepNum << "\tParentPep: " << peaksRev->at(a).start[b].parentPepNum << endl;
      }
    }
  }*/

  /*cout << endl;
  for(size_t a=0;a<ions[iIndex].peaks->size();a++){
    cout << ions[iIndex].peaks->at(a).mass << "\t" << ions[iIndex].peaks->at(a).next.size() << "\t" << ions[iIndex].peaks->at(a).start.size() << endl;
  }*/

  //get all spectra that might contain this peptide and adduct
  //Set Mass boundaries
  double minMass = ions[iIndex].pepMassMin + params.minAdductMass;
  double maxMass = ions[iIndex].pepMassMax + params.maxAdductMass;
  minMass-=(minMass/1000000*params.ppmPrecursor); //is this necessary because all adduct results end in 0ppm error?
  maxMass+=(maxMass/1000000*params.ppmPrecursor);
 
  vector<int> scanIndex;
  //cout << "Boundaries" << endl;
  if (!spec->getBoundaries(minMass, maxMass, scanIndex, scanBuffer[iIndex])) return true;


  //cout << "onward" << endl;
  for (size_t j = 0; j<scanIndex.size(); j++){
    scoreSingletSpectra2(scanIndex[j], pep.mass, len, index, minMass, maxMass, iIndex);
  }

  //cout << "Done Singlets" << endl;

  return true;
}

/*============================
  Private Functions
============================*/
bool MAnalysis::allocateMemory(int threads){
  bKIonsManager = new bool[threads];
  ions = new MIons[threads];
  scanBuffer = new bool*[threads];
  for(int i=0;i<threads;i++) {
    bKIonsManager[i]=false;
    scanBuffer[i] = new bool[spec->size()];
    for(int j=0;j<128;j++){
      ions[i].site[j]=adductSites[j];
    }
  }
  maxZ2=new int[threads];
  bufSize2=new size_t[threads];
  return true;
}

void MAnalysis::deallocateMemory(int threads){
  delete [] bKIonsManager;
  delete [] ions;
  for (int i = 0; i < threads; i++){
    delete[] scanBuffer[i];
  }
  delete[] scanBuffer;
  delete[] maxZ2;
  delete[] bufSize2;
}

//This function is way out of date. Particularly the mutexes and how to deal with multiple precursors.
//void MAnalysis::scoreSingletSpectra(int index, int sIndex, double mass, int len, int pep, char k, double minMass, double maxMass, int iIndex, bool bSiteless){
//  //cout << "scoreSingletSpectra()" << endl;
//  mScoreCard sc;
//  MIonSet* iset;
//  mPepMod mod;
//  double score=0;
//  int i,j;
//  int precI;
//  int match;
//  int conFrag;
//
//  MSpectrum* s=spec->getSpectrum(index);
//  mPrecursor* p=NULL;
//  MTopPeps* tp;
//  int sz=s->sizePrecursor();
//  double topScore=0;
//  int topMatch=0;
//  int topConFrag=0;
//
//  int code;
//  for(i=0;i<sz;i++){
//    p=s->getPrecursor2(i);
//    if (p->corr<-4) code = 2;
//    else if (p->corr<0)code = 3;
//    else if (p->corr == 0)code = 2;
//    else code = 1;
//    if(code==1) break;
//
//  }
//  if(i==sz) {
//    i=0;
//    p = s->getPrecursor2(i);
//  }
//
//  if((p->monoMass-mass)<params.minAdductMass) score=0;  //this could be narrowed down to user-defined precursor tolerance.
//  else if((p->monoMass-mass)>params.maxAdductMass) score=0;
//  else if(bSiteless) score = magnumScoring(index, 0, sIndex, iIndex, match, conFrag, p->charge); //score peptide without open mod (i.e. scores peptide without localization)
//
//  if(score>0){
//    topScore = score;
//    topMatch = match;
//    topConFrag = conFrag;
//    precI = i;
//    sc.simpleScore = (float)score;
//    sc.pep = pep;
//    sc.mass = mass;
//    sc.massA = p->monoMass - mass;
//    sc.precursor = i;
//    sc.site = -99;
//    sc.mods->clear();
//    iset = ions[iIndex].at(sIndex);
//    if (iset->difMass != 0){
//      for (j = 0; j<ions[iIndex].getIonCount(); j++) {
//        if (iset->mods[j] != 0){
//          if (j == 0 && iset->modNTerm) mod.term = true;
//          else if (j == ions[iIndex].getIonCount() - 1 && iset->modCTerm) mod.term = true;
//          else mod.term = false;
//          mod.pos = (char)j;
//          mod.mass = iset->mods[j];
//          sc.mods->push_back(mod);
//        }
//      }
//    }
//  }
//
//  for(i=0;i<sz;i++){
//    p=s->getPrecursor2(i);
//    //cout << i << " of " << sz << "\t" << p->monoMass << "\t" << minMass << "\t" << maxMass << "\t" << mass << endl;
//    if(p->monoMass<minMass) continue;
//    if(p->monoMass>maxMass) continue;
//    if ((p->monoMass - mass)>params.maxAdductMass) continue;
//    if ((p->monoMass - mass)<params.minAdductMass) continue;
//    //cout << "Before magnumScoring" << endl;
//    score=magnumScoring(index,p->monoMass-mass,sIndex,iIndex,match,conFrag,p->charge);  //open mod with localization
//    //cout << score << endl;
//    if(score==0) continue;
//    else if(score>topScore) { //replace the previous peptide scores, if this version of the peptide scores better.
//      topScore=score;
//      topMatch=match;
//      topConFrag=conFrag;
//      precI=i;
//      sc.simpleScore = (float)score;
//      sc.pep = pep;
//      sc.mass = mass;
//      sc.massA = p->monoMass - mass;
//      sc.precursor = i;
//      sc.site = k;
//      sc.mods->clear();
//      iset = ions[iIndex].at(sIndex);
//      if (iset->difMass != 0){
//        for (j = 0; j<ions[iIndex].getIonCount(); j++) {
//          if (iset->mods[j] != 0){
//            if (j == 0 && iset->modNTerm) mod.term = true;
//            else if (j == ions[iIndex].getIonCount() - 1 && iset->modCTerm) mod.term = true;
//            else mod.term = false;
//            mod.pos = (char)j;
//            mod.mass = iset->mods[j];
//            sc.mods->push_back(mod);
//          }
//        }
//      }
//    }
//  }
//
//  if(topScore>0){
//    //cout << "Topper " << topScore << endl;
//    double ev = 1000;
//    Threading::LockMutex(mutexSpecScore[index]);
//    ev = s->computeE(topScore, len);
//
//    //** temporary
//    //s->tHistogram(topScore, len);
//    //**
//
//    //cout << "DoneE " << ev << endl;
//    Threading::UnlockMutex(mutexSpecScore[index]);
//    sc.eVal=ev;
//    sc.match=topMatch;
//    sc.conFrag=topConFrag;
//
//    tp = s->getTopPeps(precI);
//    //cout << "GotTopPeps" << endl;
//    Threading::LockMutex(mutexSingletScore[index][precI]);
//    tp->checkPeptideScore(sc);
//    //cout << "CheckPepScore" << endl;
//    Threading::UnlockMutex(mutexSingletScore[index][precI]);
//  
//    Threading::LockMutex(mutexSpecScore[index]);
//    s->checkScore(sc,iIndex);
//    //cout << "CheckScore" << endl;
//    Threading::UnlockMutex(mutexSpecScore[index]);
//  } 
//
//  //** temporary
//  //else {
//  //  Threading::LockMutex(mutexSpecScore[index]);
//  //  s->tHistogram(0, len);
//  //  Threading::UnlockMutex(mutexSpecScore[index]);
//  //}
//  //**
//
//}

void MAnalysis::scoreSingletSpectra2(int index, double mass, int len, int pep, double minMass, double maxMass, int iIndex){
  mScoreCard sc;
  double score = 0;
  int precI;
  int match=0;
  int conFrag=0;

  MSpectrum* s = spec->getSpectrum(index);
  mPrecursor* p;
  MTopPeps* tp;
  int sz = s->sizePrecursor();
  double topScore = 0;
  int topMatch = 0;
  int topConFrag = 0;

  //score all peptides against all appropriate precursors
  vector<sPrecursor> pre;
  maxZ2[iIndex] = 1;
  //double maxPre=0;
  //double minPre=100000;
  if(echo)cout << "scoreSingletSpectra2 " << s->getScanNumber() << "\t" << sz << endl;
  //bEcho=true; //comment this out
  if(sz>MAX_PRECURSOR){
    cout << "WARNING: " << s->getScanNumber() << " has too many candidate precursor ions. Analysis limited to " << MAX_PRECURSOR << " precursors." << endl;
  }
  for(int i=0;i<sz;i++){
    p = s->getPrecursor2(i);
    if (p->monoMass<minMass) continue;
    if (p->monoMass>maxMass) continue;
    //if ((p->monoMass - mass)>params.maxAdductMass) continue; //Not sure here, peptides have multiple masses
    //if ((p->monoMass - mass)<params.minAdductMass) continue;
    sPrecursor pr;
    pr.index=i;
    pr.monomass=p->monoMass;
    pr.maxZ=p->charge-1;
    if(pr.maxZ<1) pr.maxZ=1;
    if(pr.maxZ>3) pr.maxZ=3;
    if (pr.maxZ>maxZ2[iIndex]) maxZ2[iIndex] = pr.maxZ;
    //if(pr.monomass>maxPre) maxPre=pr.monomass; //adjust later for ppm error?
    //if(pr.monomass<minPre) minPre=pr.monomass;
    pre.push_back(pr);
    if(pre.size()==MAX_PRECURSOR) break;
  }
  if(pre.size()==0) {
    cout << "WTF" << endl;
    exit(1);
  }

  bufSize2[iIndex] = sizeof(double)*pre.size();
  size_t pepCount=ions[iIndex].pepCount;
  //cout << "First pepCount: " << pepCount << "\t" << sizeof(double) * pre.size() << "\t" << pre.size() << endl;
  //for(size_t a=0;a<pre.size();a++){
  //  cout << pre[a].monomass << "\t" << pre[a].maxZ << "\t" << pre[a].index << endl;
  //}
  size_t preCount=pre.size();
  sScoreSet2* pScores = new sScoreSet2[pepCount];
  if (echo) cout << "Before Score9-f: " << ions[iIndex].vPeaks.size() << endl;
  score9(s, ions[iIndex].vPeaks, pScores, &pre, iIndex);

  sScoreSet2* pScores3 = new sScoreSet2[pepCount];
  if(echo) cout << "Before Score9: " << ions[iIndex].vPeaksRev.size() << endl;
  score9(s,ions[iIndex].vPeaksRev,pScores3, &pre, iIndex);

  //keep only the best score(s).
  vector<sDIndex> vTop;
  topScore=0;
  size_t minMods=100;
  for (size_t a = 0; a<pepCount; a++){
    //cout << "Peptide: " << a << "\t" << ions[iIndex].pepMass[a] << "\t" << ions[iIndex].pepLinks[a] << "\t" << ions[iIndex].pepMods[a].mods.size() << "\t" << pScores[a].scores << endl;
    double topPreScore=0;
    size_t topPreIndex=0;
    for (size_t b = 0; b<preCount; b++){

      double massA = pre[b].monomass - ions[iIndex].pepMass[a];
      if(massA<params.minAdductMass || massA>params.maxAdductMass) continue; //skip adducts outside our bounds

      double trueScore=pScores[a].score+pScores[a].scoreP[b]+pScores3[a].score+pScores3[a].scoreP[b];
      if(trueScore>topPreScore){
      //if (pScores[a].scores[b]+pScores2[a].scores[b] > topPreScore){
        topPreScore=trueScore;
        //topPreScore=pScores[a].scores[b] + pScores2[a].scores[b];
        topPreIndex=b;
      }
    }

    if (topPreScore > topScore){
      topScore = topPreScore;
      //topScore = pScores[a].scores[topPreIndex]+pScores2[a].scores[topPreIndex];
      vTop.clear();
      sDIndex di;
      di.a=a;
      di.b = topPreIndex;
      vTop.push_back(di);
      minMods = ions[iIndex].pepMods[a].mods.size();
    } else if (topPreScore == topScore){
      sDIndex di;
      di.a = a;
      di.b = topPreIndex;
      vTop.push_back(di);
      if (ions[iIndex].pepMods[a].mods.size()<minMods) minMods = ions[iIndex].pepMods[a].mods.size();
    }
   
  }

  //TODO: Combine all ambiguous localizations here
  if(topScore>0){
    for(size_t a=0;a<vTop.size();a++){
      double score= (pScores[vTop[a].a].score+pScores[vTop[a].a].scoreP[vTop[a].b] + pScores3[vTop[a].a].score+ pScores3[vTop[a].a].scoreP[vTop[a].b]) * 0.005;
      //double score = (pScores[vTop[a].a].scores[vTop[a].b]+ pScores2[vTop[a].a].scores[vTop[a].b]) *0.005;
      if (ions[iIndex].pepMods[vTop[a].a].mods.size()>minMods) continue; //skip modified peptides that are explained with fewer modifications

      topScore = score;
      //topMatch = pScores[vTop[a].a].match[vTop[a].b];
      topConFrag = conFrag;
      precI = pre[vTop[a].b].index;
      sc.simpleScore = (float)score;
      sc.pep = pep;
      sc.mass = ions[iIndex].pepMass[vTop[a].a];
      sc.massA = pre[vTop[a].b].monomass - ions[iIndex].pepMass[vTop[a].a];
      sc.precursor = pre[vTop[a].b].index;
      sc.site = ions[iIndex].pepLinks[vTop[a].a];
      sc.mods->clear();
      for (size_t c = 0; c<ions[iIndex].pepMods[vTop[a].a].mods.size(); c++){
        sc.mods->push_back(ions[iIndex].pepMods[vTop[a].a][c]);
      }

      double ev = 1000;
      Threading::LockMutex(mutexSpecScore[index]);
      ev = s->computeE(topScore, len);
      Threading::UnlockMutex(mutexSpecScore[index]);
      sc.eVal = ev;
      sc.match = 0;//topMatch;
      sc.conFrag = topConFrag;

      tp = s->getTopPeps(precI);
      Threading::LockMutex(mutexSingletScore[index][precI]);
      tp->checkPeptideScore(sc);
      Threading::UnlockMutex(mutexSingletScore[index][precI]);

      Threading::LockMutex(mutexSpecScore[index]);
      s->checkScore(sc, iIndex);
      Threading::UnlockMutex(mutexSpecScore[index]);
    }
  }

  delete [] pScores;
  //delete [] pScores2;
  delete [] pScores3;
  //peaks=NULL;

}

//could speed up this function significantly (>5%) by replacing the calls to magnumScoring2 with the code inside magnumScoring2
//void MAnalysis::score6(MSpectrum* s, sNode2* node, sLink2* link, double* score, double* scoreNL, int depth, sScoreSet* v, vector<sPrecursor>* pre, int iIndex/*, double minMass, double maxMass*/) {
//  cout << "Peptide: " << link->pepNum << "\tNode: " << node->id << "\tmass: " << node->mass << "\tsite: " << ions[iIndex].pepLinks[link->pepNum];
//
//  if (node->mass > 0) {
//
//    if (!node->visit) {
//      for (int b = 1; b <= maxZ2[iIndex]; b++){
//        double mz = (node->mass + 1.007276466*b)/b;
//        node->score[b] = node->score[b - 1] + magnumScoring2(s, mz);//score forward
//        //cout << "ScoreA: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->score[b] << endl;
//        mz = (ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466*b)/ b;
//        node->scoreAltNL[b] = node->scoreAltNL[b - 1] + magnumScoring2(s, mz);//score reverse without link
//        //cout << "ScoreA-ALtNL: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->scoreAltNL[b] << endl;
//      }
//      cout << "\tScore: " << node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//      cout << "\tScoreAltNL: " << ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 << "," << node->scoreAltNL[1] << "," << node->scoreAltNL[2];
//    } else {
//      cout << "\tCarryover-Score: " << node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//      cout << "\tCarryover-ScoreAltNL: " << ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 << "," << node->scoreAltNL[1] << "," << node->scoreAltNL[2];
//    }
//
//    for (size_t a = 0; a<pre->size(); a++){
//      link->score[a] = score[a] + node->score[pre->at(a).maxZ];
//      link->scoreNL[a] = scoreNL[a] + node->score[pre->at(a).maxZ] + node->scoreAltNL[pre->at(a).maxZ];
//      cout << "\tDepth: " << depth << "," << ions[iIndex].maxLink;
//      if (depth < ions[iIndex].maxLink) {  
//        if (!node->visit) {
//          for (int b = 1; b <= maxZ2[iIndex]; b++) {
//            double mz = (pre->at(a).monomass - node->mass + 1.007276466*b) / b;
//            node->scoreAlt[a][b] = node->scoreAlt[a][b - 1] + magnumScoring2(s, mz); //score reverse after precursor subtraction
//            //cout << "Pre: " << a << "\tScoreA-Alt: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->scoreAlt[a][b] << endl;
//          }
//          cout << "\tScoreAlt: " << pre->at(a).monomass - node->mass + 1.007276466 << "," << node->scoreAlt[a][1] << "," << node->scoreAlt[a][2];
//        } else {
//          cout << "\tCarryover-Alt: " << pre->at(a).monomass - node->mass + 1.007276466 << "," << node->scoreAlt[a][1] << "," << node->scoreAlt[a][2];
//        }
//        link->score[a] += node->scoreAlt[a][pre->at(a).maxZ];
//        //cout << "Pre: " << a << " scoreA final: " << link->score[a] << endl;
//      }
//    }
//    cout << "\tLinkScore: " << link->score[0] << "," << link->score[1];
//    cout << "\tLinkScoreNL: " << link->scoreNL[0] << "," << link->scoreNL[1];
//
//  } else {
//
//    if (!node->visit) {
//      for (int b = 1; b <= maxZ2[iIndex]; b++) {
//        double mz = (ions[iIndex].pepMass[link->pepNum] + node->mass + 1.007276466*b) / b;
//        node->score[b] = node->score[b - 1] + magnumScoring2(s, mz);//score forward
//        //cout << "ScoreB: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->score[b] << endl;
//      }
//      cout << "\tScore: " << ions[iIndex].pepMass[link->pepNum] + node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//    } else {
//      cout << "\tCarryover-Score: " << ions[iIndex].pepMass[link->pepNum] + node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//    }
//    for (size_t a = 0; a < pre->size(); a++) {
//      if (!node->visit) {
//        for (int b = 1; b <= maxZ2[iIndex]; b++) {
//          double mz = (pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466*b) / b;
//          node->scoreAlt[a][b] = node->scoreAlt[a][b - 1] + magnumScoring2(s, mz); //score after precursor subtraction
//          //cout << "ScoreB-Alt: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->scoreAlt[a][b] << "\t" << ions[iIndex].pepMass[link->pepNum] << "\t" << link->pepNum << "\t" << pre->at(a).monomass  << "\t" << node->mass << endl;
//        }
//        cout << "\tScoreAlt: " << pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 << "," << node->scoreAlt[a][1] << "," << node->scoreAlt[a][2];
//      } else {
//        cout << "\tCarryover-Alt: " << pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 << "," << node->scoreAlt[a][1] << "," << node->scoreAlt[a][2];
//      }
//      link->score[a] = score[a] + node->scoreAlt[a][pre->at(a).maxZ] + node->score[pre->at(a).maxZ];
//      link->scoreNL[a] = scoreNL[a] + node->scoreAlt[a][pre->at(a).maxZ];
//    }
//    cout << "\tLinkScore: " << link->score[0] << "," << link->score[1];
//    cout << "\tLinkScoreNL: " << link->scoreNL[0] << "," << link->scoreNL[1];
//
//  }
//
//  node->visit = true;
//  cout << endl;
//  if (link->nextNode == SIZE_MAX) {
//    //cout << "In memswap" << endl;
//    if (ions[iIndex].pepLinks[link->pepNum] < 0) {
//      //cout << ions[iIndex].pepLinks[link->pepNum] << " no localization" << endl;
//      memcpy(v[link->pepNum].scores, link->scoreNL, bufSize2[iIndex]);
//      cout << link->pepNum << " scores now " << v[link->pepNum].scores[0] << "," << v[link->pepNum].scores[1] << endl;
//      //cout << "Copy NL scores for: " << link->pepNum << "\t" << bufSize2[iIndex] << endl;
//    } else {
//      memcpy(v[link->pepNum].scores, link->score, bufSize2[iIndex]);
//      cout << link->pepNum << " scores now " << v[link->pepNum].scores[0] << "," << v[link->pepNum].scores[1] << endl;
//    }
//  } else score6(s, &ions[iIndex].peaks->at(link->nextNode), &ions[iIndex].peaks->at(link->nextNode).next[link->nextIndex], link->score, link->scoreNL, depth + 1, v, pre, iIndex/*, minMass, maxMass*/);
//
//  for (size_t a = 0; a < node->start.size(); a++) {
//    //cout << "Iterate from " << link->nextNode << "\t" << a << "\t" << node->start[a].pepNum << "\t" << link->pepNum << "\tScored: " << v[node->start[a].pepNum].scored << endl;
//    if (v[node->start[a].pepNum].scored) continue; //maybe create a structure and flag instead?
//    if (node->start[a].pepNum < link->pepNum) continue;
//    //if (ions[iIndex].pepMass[node->start[a].pepNum] <minMass || ions[iIndex].pepMass[node->start[a].pepNum]>maxMass) continue; //interesting!! skip peptides of too high mass - doesn't work for low masses (other than fewer iterations)
//    v[node->start[a].pepNum].scored = true;
//    //cout << "Next node: " << node->start[a].nextNode << endl;
//    if (node->start[a].nextNode == SIZE_MAX) {
//      cout << "Obscure score position: " << ions[iIndex].pepLinks[node->start[a].pepNum] << endl;
//      if (ions[iIndex].pepLinks[node->start[a].pepNum] < 0)  memcpy(v[node->start[a].pepNum].scores, link->scoreNL, bufSize2[iIndex]);
//      else memcpy(v[node->start[a].pepNum].scores, link->score, bufSize2[iIndex]);
//    } else score6(s, &ions[iIndex].peaks->at(node->start[a].nextNode), &ions[iIndex].peaks->at(node->start[a].nextNode).next[node->start[a].nextIndex], link->score, link->scoreNL, depth + 1, v, pre, iIndex/*, minMass, maxMass*/);
//  }
//  cout << "Done " << node->id << endl;
//
//}

//could speed up this function significantly (>5%) by replacing the calls to magnumScoring2 with the code inside magnumScoring2
//void MAnalysis::score7(MSpectrum* s, vector<sNode2>* peakSet, sNode2* node, sLink2* link, double* score, double* scoreNL, int depth, sScoreSet* v, vector<sPrecursor>* pre, int iIndex/*, double minMass, double maxMass*/) {
//  //sCounter++;
//  //cout << "score7 : " << sCounter << endl;
//  //if(bEcho) cout << "Peptide: " << link->pepNum << "\tNode: " << node->id << "\tmass: " << node->mass << "\tsite: " << ions[iIndex].pepLinks[link->pepNum] << "\t";
//
//  if (node->mass > 0) {
//
//    if (!node->visit) {
//      //cout << link->pepNum << "," << node->mass;
//      for (int b = 1; b <= maxZ2[iIndex]; b++) {
//        double mz = (node->mass + 1.007276466 * b) / b;
//        node->score[b] = node->score[b - 1] + magnumScoring2(s, mz);//score forward
//        //cout << "\t" << mz << "," << (int)magnumScoring2(s, mz);
//        //cout << "ScoreA: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->score[b] << endl;
//      }
//      //cout << endl;
//      //if (bEcho) cout << "\tScore: " << node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//    } else {
//      //if (bEcho) cout << "\tCarryover-Score: " << node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//    }
//
//    for (size_t a = 0; a < pre->size(); a++) {
//      link->score[a] = score[a] + node->score[pre->at(a).maxZ];
//      link->scoreNL[a] = scoreNL[a] + node->score[pre->at(a).maxZ];// + node->scoreAltNL[pre->at(a).maxZ];
//      //cout << "\tDepth: " << depth << "," << ions[iIndex].maxLink;
//    }
//    //cout << "\tLinkScore: " << link->score[0] << "," << link->score[1];
//    //cout << "\tLinkScoreNL: " << link->scoreNL[0] << "," << link->scoreNL[1];
//
//  } else {
//
//    //if (!node->visit) {
//    //  for (int b = 1; b <= maxZ2[iIndex]; b++) {
//    //    double mz = (ions[iIndex].pepMass[link->pepNum] + node->mass + 1.007276466 * b) / b;
//    //    node->score[b] = node->score[b - 1] + magnumScoring2(s, mz);//score forward
//    //    cout << "ScoreB: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->score[b] << endl;
//    //  }
//    //  cout << "\tScore: " << ions[iIndex].pepMass[link->pepNum] + node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//    //} else {
//    //  cout << "\tCarryover-Score: " << ions[iIndex].pepMass[link->pepNum] + node->mass + 1.007276466 << "," << node->score[1] << "," << node->score[2];
//    //}
//    for (size_t a = 0; a < pre->size(); a++) {
//      if (!node->visit) {
//        //cout << link->pepNum << "," << node->mass;
//        for (int b = 1; b <= maxZ2[iIndex]; b++) {
//          double mz = (pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 * b) / b;
//          node->scoreAlt[a][b] = node->scoreAlt[a][b - 1] + magnumScoring2(s, mz); //score after precursor subtraction
//          //cout << "\t" << mz << "," << (int)magnumScoring2(s, mz);
//          //cout << "ScoreB-Alt: " << mz << " " << b << "\t" << (int)magnumScoring2(s, mz) << "\t" << node->scoreAlt[a][b] << "\t" << ions[iIndex].pepMass[link->pepNum] << "\t" << link->pepNum << "\t" << pre->at(a).monomass  << "\t" << node->mass << endl;
//        }
//        //cout << endl;
//        //if (bEcho) cout << "\tScoreAlt: " << pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 << "," << node->scoreAlt[a][1] << "," << node->scoreAlt[a][2];
//      } else {
//        //if (bEcho) cout << "\tCarryover-Alt: " << pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 << "," << node->scoreAlt[a][1] << "," << node->scoreAlt[a][2];
//      }
//      link->score[a] = score[a] + node->scoreAlt[a][pre->at(a).maxZ] + node->score[pre->at(a).maxZ];
//      link->scoreNL[a] = scoreNL[a] + node->scoreAlt[a][pre->at(a).maxZ];
//    }
//
//    //cout << "\tLinkScore: " << link->score[0] << "," << link->score[1];
//    //cout << "\tLinkScoreNL: " << link->scoreNL[0] << "," << link->scoreNL[1];
//
//  }
//
//  node->visit = true;
//
//  //if (bEcho) cout << endl;
//  if (link->nextNode == SIZE_MAX) {
//    //cout << "In memswap" << endl;
//    if (ions[iIndex].pepLinks[link->pepNum] < 0) {
//      //cout << ions[iIndex].pepLinks[link->pepNum] << " no localization" << endl;
//      memcpy(v[link->pepNum].scores, link->scoreNL, bufSize2[iIndex]);
//      //cout << link->pepNum << " scores now " << v[link->pepNum].scores[0] << "," << v[link->pepNum].scores[1] << endl;
//      //cout << "Copy NL scores for: " << link->pepNum << "\t" << bufSize2[iIndex] << endl;
//    } else {
//      memcpy(v[link->pepNum].scores, link->score, bufSize2[iIndex]);
//      //cout << link->pepNum << " scores now " << v[link->pepNum].scores[0] << "," << v[link->pepNum].scores[1] << endl;
//    }
//  } else {
//    //cout << "Continuing: " << link->nextNode << "\tindex: " << link->nextIndex << endl;
//    score7(s, peakSet, &peakSet->at(link->nextNode), &peakSet->at(link->nextNode).next[link->nextIndex], link->score, link->scoreNL, depth + 1, v, pre, iIndex/*, minMass, maxMass*/);
//    //cout << "Gniunitnoc: " << link->nextNode << "\tindex: " << link->nextIndex << endl;
//  }
//
//  for (size_t a = 0; a < node->start.size(); a++) { //Change this to only do starts at first time to node?
//    //cout << "Iterate from " << node->id << "\t" << a << "\tNewPep: " << node->start[a].pepNum << "\tOldPep: " << link->pepNum << "\tScored: " << v[node->start[a].pepNum].scored << endl;
//    if (v[node->start[a].pepNum].scored) continue; //maybe create a structure and flag instead?
//    //if (node->start[a].pepNum < link->pepNum) continue;
//    //if (ions[iIndex].pepMass[node->start[a].pepNum] <minMass || ions[iIndex].pepMass[node->start[a].pepNum]>maxMass) continue; //interesting!! skip peptides of too high mass - doesn't work for low masses (other than fewer iterations)
//    if(node->start[a].parentPepNum!=link->pepNum) continue; //skip start if it comes from a different peptide...must still iterate all starts (slow)
//    v[node->start[a].pepNum].scored = true;
//    if (node->start[a].nextNode == SIZE_MAX) {
//      cout << "Obscure score position: " << ions[iIndex].pepLinks[node->start[a].pepNum] << endl;
//      if (ions[iIndex].pepLinks[node->start[a].pepNum] < 0)  memcpy(v[node->start[a].pepNum].scores, link->scoreNL, bufSize2[iIndex]);
//      else memcpy(v[node->start[a].pepNum].scores, link->score, bufSize2[iIndex]);
//    } else {
//      //if (bEcho) cout << "Starting: " << node->start[a].pepNum << "\tnode: " << node->start[a].nextNode << "\tindex: " << node->start[a].nextIndex << endl;
//      score7(s, peakSet, &peakSet->at(node->start[a].nextNode), &peakSet->at(node->start[a].nextNode).next[node->start[a].nextIndex], link->score, link->scoreNL, depth + 1, v, pre, iIndex/*, minMass, maxMass*/);
//      //if(bEcho) cout << "Gnitrats: " << node->start[a].pepNum << "\tnode: " << node->start[a].nextNode << "\tindex: " << node->start[a].nextIndex << endl;
//    }
//  }
//
//  //if (bEcho) cout << "Done " << node->id << endl;
//  //sCounter--;
//  //cout << "end Score7: " << sCounter << endl;
//
//}

//could speed up this function significantly (>5%) by replacing the calls to magnumScoring2 with the code inside magnumScoring2
//void MAnalysis::score8(MSpectrum* s, vector<sNode2>* peakSet, sNode2* node, sLink2* link, sScoreSet* v, vector<sPrecursor>* pre, int iIndex) {
//  //cout << "score8: " << node->id << " pepNum:" << link->pepNum;
//  //First visit to node means it needs to be scored
//  if(!node->visit){
//    if (node->mass > 0) {
//      for (int b = 1; b <= maxZ2[iIndex]; b++) {
//        double mz = (node->mass + 1.007276466 * b) / b;
//        node->score[b] = node->score[b - 1] + magnumScoring2(s, mz);//score forward
//        //cout << "\tM:"<<mz<<","<<node->score[b];
//      }
//    } else {
//      node->lastMass= pre->at(0).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass;
//      for (size_t a = 0; a < pre->size(); a++) {
//        for (int b = 1; b <= maxZ2[iIndex]; b++) {
//          double mz = (pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 * b) / b;
//          node->scoreAlt[a][b] = node->scoreAlt[a][b - 1] + magnumScoring2(s, mz); //score after precursor subtraction
//          //cout << "\tA:" << mz <<","<< node->scoreAlt[a][b];
//        }
//        //cout << " Pre:" << pre->at(a).monomass;
//        //cout << " Ion:" << ions[iIndex].pepMass[link->pepNum];
//      }
//    }
//    node->visit = true;
//
//  } else { 
//
//    //but repeat visits may need to be rescored if the peptide mass has changed (due to differential modifications)...
//    if (node->mass <0) {
//      double mass= pre->at(0).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass - node->lastMass;
//      if(mass>0.01 || mass<0.01){
//        for (size_t a = 0; a < pre->size(); a++) {
//          for (int b = 1; b <= maxZ2[iIndex]; b++) {
//            double mz = (pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 * b) / b;
//            node->scoreAlt[a][b] = node->scoreAlt[a][b - 1] + magnumScoring2(s, mz); //score after precursor subtraction  
//            //cout << "\tR:" << mz << "," << node->scoreAlt[a][b];
//          }
//        }
//      }
//    }
//  }
//  //cout << endl;
//
//  //Add the node score to the running tally for the peptide
//  if(node->mass>0){
//    for(size_t a=0;a<pre->size();a++) v[link->pepNum].scores[a]+=node->score[maxZ2[iIndex]];
//  } else {
//    for (size_t a = 0; a < pre->size(); a++) v[link->pepNum].scores[a] += node->scoreAlt[a][maxZ2[iIndex]];
//  }
//  //cout << "  pepnum:" << link->pepNum << "\t" << v[link->pepNum].scores[0] << endl;
//  
//
//
//  //Branch to next peptide, if necessary
//  for (size_t a = 0; a < node->start.size(); a++) {
//    //if(echo) cout << "Node: " << node->id << " Check: " << node->start[a].pepNum << endl;
//    if (v[node->start[a].pepNum].scored) continue; //maybe create a structure and flag instead?
//    if (node->start[a].parentPepNum != link->pepNum) continue; //skip start if it comes from a different peptide...must still iterate all starts (slow)
//    v[node->start[a].pepNum].scored = true;
//    //if (echo) cout << "Node: " << node->id << " Check: " << node->start[a].pepNum << " scored." << endl;
//
//    //copy over scores to the new peptide
//    memcpy(v[peakSet->at(node->start[a].nextNode).next[node->start[a].nextIndex].pepNum].scores, v[link->pepNum].scores, bufSize2[iIndex]);
//    //for(size_t x=0;x<bufSize2[iIndex];x++){
//    //  cout << "sc: " << x << "\t" << v[peakSet->at(node->start[a].nextNode).next[node->start[a].nextIndex].pepNum].scores[x] << endl;
//    //}
//    score8(s, peakSet, &peakSet->at(node->start[a].nextNode), &peakSet->at(node->start[a].nextNode).next[node->start[a].nextIndex], v, pre, iIndex);
//  }
//
//    //Continue to next amino acid in this peptide
//  if (link->nextNode != SIZE_MAX) {
//    score8(s, peakSet, &peakSet->at(link->nextNode), &peakSet->at(link->nextNode).next[link->nextIndex], v, pre, iIndex);
//  }
//
//}
//
////backup until checked in
////void MAnalysis::score9(MSpectrum* s, vector<sIPeak>& peakSet, sScoreSet2* ss, vector<sPrecursor>* pre, int iIndex) {
////  for (size_t a = 0; a < peakSet.size(); a++) {
////    if (peakSet[a].mass > 0) {
////      double score = 0;
////      if (echo) cout << a << " score9 A: " << peakSet[a].mass << endl;
////      for (int b = 1; b <= maxZ2[iIndex]; b++) {
////        double mz = (peakSet[a].mass + 1.007276466 * b) / b;
////        score += magnumScoring2(s, mz);
////      }
////      //if(score!=0){
////      for (size_t b = 0; b < peakSet[a].pepIndex.size(); b++) {
////        if (echo) cout << "Add " << score << " to " << peakSet[a].pepIndex[b] << endl;
////        ss[peakSet[a].pepIndex[b]].score += score;
////      }
////      //}
////    } else {
////      //if two peptides have the same precursor mass, then scoring is repeated.
////      //TODO: see if we can avoid the repeat...
////      for (size_t b = 0; b < peakSet[a].pepIndex.size(); b++) {
////        for (size_t c = 0; c < pre->size(); c++) {
////          double score = 0;
////          if (echo) cout << a << " score9 b: " << peakSet[a].mass << " pep:" << peakSet[a].pepIndex[b] << " pre: " << c << " becomes " << pre->at(c).monomass - ions[iIndex].pepMass[peakSet[a].pepIndex[b]] - peakSet[a].mass << endl;
////          for (int d = 1; d <= maxZ2[iIndex]; d++) {
////            double mz = (pre->at(c).monomass - ions[iIndex].pepMass[peakSet[a].pepIndex[b]] - peakSet[a].mass + 1.007276466 * d) / d;
////            if (mz < 0) continue;
////            score += magnumScoring2(s, mz);
////          }
////          //if(echo) cout << "Add " << score << " to " << peakSet[a].pepIndex[b] << endl;
////          ss[peakSet[a].pepIndex[b]].scoreP[c] += score;
////        }
////      }
////    }
////  }
////}

void MAnalysis::score9(MSpectrum* s, vector<sIPeak>& peakSet, sScoreSet2* ss, vector<sPrecursor>* pre, int iIndex){
  for(size_t a=0;a<peakSet.size();a++){
    if (peakSet[a].mass > 0) {
      double score=0;
      //if(echo) cout << a << " score9 A: " << peakSet[a].mass << endl;
      for (int b = 1; b <= maxZ2[iIndex]; b++) {
        double mz = (peakSet[a].mass + 1.007276466 * b) / b;
        score+=magnumScoring2(s, mz);
      }
      //if(score!=0){
        for (size_t b = 0; b < peakSet[a].index[0].pepIndex.size(); b++) {
          //if(echo) cout << "Add " << score << " to " << peakSet[a].index[0].pepIndex[b] << endl;
          ss[peakSet[a].index[0].pepIndex[b]].score += score;
        }
      //}
    } else {
      //if two peptides have the same precursor mass, then scoring is repeated.
      //TODO: see if we can avoid the repeat...
      for(size_t b=0;b< peakSet[a].index.size();b++){
        for (size_t c = 0; c < pre->size(); c++) {
          double score = 0;
          //if(echo) cout << a << " score9 b: " << peakSet[a].mass << " pep:"<< peakSet[a].index[b].pepMass << " pre: " << c << " becomes " << pre->at(c).monomass- peakSet[a].index[b].pepMass - peakSet[a].mass << endl;
          for (int d = 1; d <= maxZ2[iIndex]; d++) {
            double mz = (pre->at(c).monomass - peakSet[a].index[b].pepMass - peakSet[a].mass + 1.007276466 *d) /d;
            if(mz<0) continue;
            score += magnumScoring2(s, mz);
          }
          //if(score!=0){
          for (size_t d = 0; d < peakSet[a].index[b].pepIndex.size(); d++) {
            //if (echo) cout << "Add " << score << " to " << peakSet[a].index[b].pepIndex[d] << endl;
            ss[peakSet[a].index[b].pepIndex[d]].scoreP[c] += score;
          }
          //}
        }
      }
    }
  }
}

//TODO: Fix inefficiencies. Right now all precursor variations are scored, including those that are the wrong mass.
void MAnalysis::score9solo(MSpectrum* s, vector<sIPeak>& peakSet, sScoreSet2* ss, int maxZ, int iIndex) {
  for (size_t a = 0;a < peakSet.size();a++) {
    double score = 0;
    //if (echo) cout << a << " score9solo A: " << peakSet[a].mass << endl;
    for (int b = 1; b <= maxZ; b++) {
      double mz = (peakSet[a].mass + 1.007276466 * b) / b;
      //cout << " charge " << b << "\t" << mz << endl;
      score += magnumScoring2(s, mz);
    }
    for (size_t b = 0; b < peakSet[a].index[0].pepIndex.size(); b++) {
      //if (echo) cout << "Add " << score << " to " << peakSet[a].index[0].pepIndex[b] << endl;
      ss[peakSet[a].index[0].pepIndex[b]].score += score;
    }
  }
}

//could speed up this function significantly (>5%) by replacing the calls to magnumScoring2 with the code inside magnumScoring2
//void MAnalysis::score6solo(MSpectrum* s, sNode2* node, sLink2* link, double* score, double* scoreNL, sScoreSet* v, sPrecursor& pre, int iIndex, double maxMass) {
//
//  if (!node->visit) {
//    for (int b = 1; b <= pre.maxZ; b++) {
//      double mz = (node->mass + 1.007276466 * b) / b;
//      node->score[b] = node->score[b - 1] + magnumScoring2(s, mz);//score forward
//      mz = (ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 * b) / b;
//      node->scoreAltNL[b] = node->scoreAltNL[b - 1] + magnumScoring2(s, mz);//score reverse without link
//    }
//  }
//
//  link->score[0] = score[0] + node->score[pre.maxZ];
//  link->scoreNL[0] = scoreNL[0] + node->score[pre.maxZ] + node->scoreAltNL[pre.maxZ];
//  node->visit = true;
//
//  if (link->nextNode == SIZE_MAX) v[link->pepNum].scores[0]=link->scoreNL[0];
//  else score6solo(s, &ions[iIndex].peaks->at(link->nextNode), &ions[iIndex].peaks->at(link->nextNode).next[link->nextIndex], link->score, link->scoreNL, v, pre, iIndex, maxMass);
//
//  for (size_t a = 0; a < node->start.size(); a++) {
//    if (v[node->start[a].pepNum].scored) continue; //maybe create a structure and flag instead?
//    if (node->start[a].pepNum < link->pepNum) continue;
//    //if (ions[iIndex].pepMass[node->start[a].pepNum]>maxMass) continue; //interesting!! skip peptides of too high mass - doesn't work for low masses (other than fewer iterations)
//    v[node->start[a].pepNum].scored = true;
//    if (node->start[a].nextNode == SIZE_MAX)  v[link->pepNum].scores[0] = link->scoreNL[0];
//    else score6solo (s, &ions[iIndex].peaks->at(node->start[a].nextNode), &ions[iIndex].peaks->at(node->start[a].nextNode).next[node->start[a].nextIndex], link->score, link->scoreNL, v, pre, iIndex, maxMass);
//  }
//
//}

//could speed up this function significantly (>5%) by replacing the calls to magnumScoring2 with the code inside magnumScoring2
//void MAnalysis::score7solo(MSpectrum* s, vector<sNode2>* peakSet, sNode2* node, sLink2* link, double* score, double* scoreNL, sScoreSet* v, sPrecursor& pre, int iIndex, double maxMass) {
//
//  if (!node->visit) {
//    for (int b = 1; b <= pre.maxZ; b++) {
//      double mz = (node->mass + 1.007276466 * b) / b;
//      node->score[b] = node->score[b - 1] + magnumScoring2(s, mz);//score forward
//      //mz = (ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466 * b) / b;
//      //node->scoreAltNL[b] = node->scoreAltNL[b - 1] + magnumScoring2(s, mz);//score reverse without link
//    }
//  }
//
//  link->score[0] = score[0] + node->score[pre.maxZ];
//  link->scoreNL[0] = scoreNL[0] + node->score[pre.maxZ];// + node->scoreAltNL[pre.maxZ];
//  node->visit = true;
//
//  if (link->nextNode == SIZE_MAX) v[link->pepNum].scores[0] = link->scoreNL[0];
//  else score7solo(s, peakSet, &peakSet->at(link->nextNode), &peakSet->at(link->nextNode).next[link->nextIndex], link->score, link->scoreNL, v, pre, iIndex, maxMass);
//
//  for (size_t a = 0; a < node->start.size(); a++) {
//    if (v[node->start[a].pepNum].scored) continue; //maybe create a structure and flag instead?
//    //if (node->start[a].pepNum < link->pepNum) continue;
//    //if (ions[iIndex].pepMass[node->start[a].pepNum]>maxMass) continue; //interesting!! skip peptides of too high mass - doesn't work for low masses (other than fewer iterations)
//    if (node->start[a].parentPepNum != link->pepNum) continue;
//    v[node->start[a].pepNum].scored = true;
//    if (node->start[a].nextNode == SIZE_MAX)  v[link->pepNum].scores[0] = link->scoreNL[0];
//    else score7solo(s, peakSet, &peakSet->at(node->start[a].nextNode), &peakSet->at(node->start[a].nextNode).next[node->start[a].nextIndex], link->score, link->scoreNL, v, pre, iIndex, maxMass);
//  }
//
//}

//could speed up this function significantly (>5%) by replacing the calls to magnumScoring2 with the code inside magnumScoring2
//this version tries to count matched peaks, but it is currently buggy...
//void MAnalysis::score7(MSpectrum* s, sNode2* node, sLink2* link, double* score, double* scoreNL, int* match, int* matchNL, int depth, sScoreSet* v, vector<sPrecursor>* pre, int iIndex, int maxZ, size_t bufSize, size_t bufSizeM/*, double minMass, double maxMass*/) {
//
//  if (node->mass > 0) {
//
//    if (!node->visit) {
//      for (int b = 1; b <= maxZ; b++){
//        double mz = (node->mass + 1.007276466*b) / b;
//        char c = magnumScoring2(s, mz);//score forward
//        node->score[b] = node->score[b - 1] + c;
//        if (c>5) node->match[b] = node->match[b - 1] + 1;
//        mz = (ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466*b) / b;
//        c = magnumScoring2(s, mz);//score reverse without link
//        node->scoreAltNL[b] = node->scoreAltNL[b - 1] + c;
//        if (c>5) node->matchAltNL[b] = node->matchAltNL[b - 1] + 1;
//      }
//    }
//
//    for (size_t a = 0; a<pre->size(); a++){
//      link->score[a] = score[a] + node->score[pre->at(a).maxZ];
//      link->scoreNL[a] = scoreNL[a] + node->score[pre->at(a).maxZ] + node->scoreAltNL[pre->at(a).maxZ];
//      link->match[a] = match[a] + node->match[pre->at(a).maxZ];
//      link->matchNL[a] = matchNL[a] + node->match[pre->at(a).maxZ] + node->matchAltNL[pre->at(a).maxZ];
//
//      if (depth < ions[iIndex].maxLink) {
//        if (!node->visit) {
//          for (int b = 1; b <= maxZ; b++) {
//            double mz = (pre->at(a).monomass - node->mass + 1.007276466*b) / b;
//            char c = magnumScoring2(s, mz); //score reverse after precursor subtraction
//            node->scoreAlt[a][b] = node->scoreAlt[a][b - 1] + c;
//            if (c>5) node->matchAlt[a][b] = node->matchAlt[a][b - 1] + 1;
//          }
//        }
//        link->score[a] += node->scoreAlt[a][pre->at(a).maxZ];
//        link->match[a] += node->matchAlt[a][pre->at(a).maxZ];
//      }
//    }
//
//  } else {
//
//    if (!node->visit) {
//      for (int b = 1; b <= maxZ; b++) {
//        double mz = (ions[iIndex].pepMass[link->pepNum] + node->mass + 1.007276466*b) / b;
//        char c = magnumScoring2(s, mz);//score forward
//        node->score[b] = node->score[b - 1] + c;
//        if (c>5) node->match[b] = node->match[b - 1] + 1;
//      }
//    }
//    for (size_t a = 0; a < pre->size(); a++) {
//      if (!node->visit) {
//        for (int b = 1; b <= maxZ; b++) {
//          double mz = (pre->at(a).monomass - ions[iIndex].pepMass[link->pepNum] - node->mass + 1.007276466*b) / b;
//          char c = magnumScoring2(s, mz); //score after precursor subtraction  
//          node->scoreAlt[a][b] = node->scoreAlt[a][b - 1] + c;
//          if (c>5) node->matchAlt[a][b] = node->matchAlt[a][b - 1] + 1;
//        }
//      }
//      link->score[a] = score[a] + node->scoreAlt[a][pre->at(a).maxZ] + node->score[pre->at(a).maxZ];
//      link->scoreNL[a] = scoreNL[a] + node->scoreAlt[a][pre->at(a).maxZ];
//      link->match[a] = match[a] + node->matchAlt[a][pre->at(a).maxZ] + node->match[pre->at(a).maxZ];
//      link->matchNL[a] = matchNL[a] + node->matchAlt[a][pre->at(a).maxZ];
//    }
//
//  }
//
//  node->visit = true;
//  if (link->nextNode == SIZE_MAX) {
//    if (ions[iIndex].pepLinks[link->pepNum] < 0) {
//      memcpy(v[link->pepNum].scores, link->scoreNL, bufSize);
//      memcpy(v[link->pepNum].match, link->matchNL, bufSizeM);
//    } else {
//      memcpy(v[link->pepNum].scores, link->score, bufSize);
//      memcpy(v[link->pepNum].match, link->match, bufSizeM);
//    }
//  } else score7(s, &ions[iIndex].peaks->at(link->nextNode), &ions[iIndex].peaks->at(link->nextNode).next[link->nextIndex], link->score, link->scoreNL, link->match, link->matchNL, depth + 1, v, pre, iIndex, maxZ, bufSize, bufSizeM/*, minMass, maxMass*/);
//
//  for (size_t a = 0; a < node->start.size(); a++) {
//    if (v[node->start[a].pepNum].scored) continue; //maybe create a structure and flag instead?
//    if (node->start[a].pepNum < link->pepNum) continue;
//    //if (ions[iIndex].pepMass[node->start[a].pepNum] <minMass || ions[iIndex].pepMass[node->start[a].pepNum]>maxMass) continue; //interesting!! skip peptides of too high mass - doesn't work for low masses (other than fewer iterations)
//    v[node->start[a].pepNum].scored = true;
//    if (node->start[a].nextNode == SIZE_MAX) {
//      if (ions[iIndex].pepLinks[node->start[a].pepNum] < 0) {
//        memcpy(v[node->start[a].pepNum].scores, link->scoreNL, bufSize);
//        memcpy(v[node->start[a].pepNum].match, link->matchNL, bufSizeM);
//      } else {
//        memcpy(v[node->start[a].pepNum].scores, link->score, bufSize);
//        memcpy(v[node->start[a].pepNum].match, link->match, bufSizeM);
//      }
//    } else score7(s, &ions[iIndex].peaks->at(node->start[a].nextNode), &ions[iIndex].peaks->at(node->start[a].nextNode).next[node->start[a].nextIndex], link->score, link->scoreNL, link->match, link->matchNL, depth + 1, v, pre, iIndex, maxZ, bufSize, bufSizeM/*, minMass, maxMass*/);
//  }
//
//}

//void MAnalysis::scoreSpectra(vector<int>& index, int sIndex, int len, double modMass, int pep1, int pep2, int k1, int k2, int link, int iIndex){
//  unsigned int a;
//  int i,z,ps;
//  mScoreCard sc;
//  mPepMod mod;
//  double mass = ions[iIndex][sIndex].mass;
//  mPrecursor* p=NULL;
//  MTopPeps* tp=NULL;
//
//  //score spectra
//  for(a=0;a<index.size();a++){
//    
//    //find the specific precursor mass in this spectrum to identify the charge state
//    z=0;
//    for (ps = 0; ps<spec->at(index[a]).sizePrecursor(); ps++){
//      p = spec->at(index[a]).getPrecursor2(ps);
//      double ppm = (p->monoMass - mass) / mass*1e6;
//      if (ppm<params.ppmPrecursor && ppm>-params.ppmPrecursor){
//        z = p->charge;
//        tp = spec->at(index[a]).getTopPeps(ps);
//        break;
//      }
//    }
//    
//    sc.simpleScore=magnumScoring(index[a],modMass,sIndex,iIndex,sc.match,sc.conFrag,z);
//    if(sc.simpleScore==0) {
//
//      //Threading::LockMutex(mutexSpecScore[index[a]]);
//      //if (spec->at(index[a]).hp[iIndex].pepIndex != pep1){ //first score
//      //  spec->at(index[a]).hp[iIndex].pepIndex = pep1;
//      //  spec->at(index[a]).hp[iIndex].topScore = 0;
//      //  spec->at(index[a]).histogram[0]++;
//      //  spec->at(index[a]).histogramCount++;
//      //} else {
//      //  //do nothing, can't be a better score than what is already there.
//      //}
//      //Threading::UnlockMutex(mutexSpecScore[index[a]]);
//
//      continue;
//    }
//
//    double ev = 1000;
//    Threading::LockMutex(mutexSpecScore[index[a]]);
//    ev = spec->at(index[a]).computeE(sc.simpleScore, len);
//    Threading::UnlockMutex(mutexSpecScore[index[a]]);
//
//    sc.eVal=ev;
//    sc.mods->clear();
//    sc.site=k1;
//    sc.mass=mass;
//    sc.massA=0;
//    sc.pep=pep1;
//    sc.precursor=(char)ps;
//    if(ions[iIndex][sIndex].difMass!=0){
//      for(i=0;i<ions[iIndex].getPeptideLen();i++) {
//        if(ions[iIndex][sIndex].mods[i]!=0){
//          if (i == 0 && ions[iIndex][sIndex].modNTerm) mod.term = true;
//          else if (i == ions[iIndex].getIonCount() - 1 && ions[iIndex][sIndex].modCTerm) mod.term = true;
//          else mod.term = false;
//          mod.pos=(char)i;
//          mod.mass=ions[iIndex][sIndex].mods[i];
//          sc.mods->push_back(mod);
//        }
//      }
//    }
//
//    Threading::LockMutex(mutexSingletScore[index[a]][ps]);
//    tp->checkPeptideScore(sc);
//    Threading::UnlockMutex(mutexSingletScore[index[a]][ps]);
//
//    Threading::LockMutex(mutexSpecScore[index[a]]);
//    spec->at(index[a]).checkScore(sc,iIndex);
//    Threading::UnlockMutex(mutexSpecScore[index[a]]);
//
//  }
//
//  p=NULL;
//  tp=NULL;
//
//}

void MAnalysis::scoreSpectra2(vector<int>& index, double mass, int len, int pep1, int iIndex) {
  unsigned int a;
  int ps;
  mScoreCard sc;
  mPrecursor* p = NULL;
  MTopPeps* tp = NULL;
  MSpectrum* s=NULL;

  if(echo) cout << "scoreSpectra2 " << mass << " " << pep1 << endl;

  //score spectra
  for (a = 0; a < index.size(); a++) {
    s=spec->getSpectrum(index[a]);
    if(echo) cout << s->getScanNumber() << endl;

    //find the specific precursor mass in this spectrum to identify the charge state
    sPrecursor pre;
    pre.index=-1;
    for (int i = 0; i < s->sizePrecursor(); i++) {
      p = s->getPrecursor2(i);
      double ppm = (p->monoMass - mass) / mass * 1e6;
      if (ppm<params.ppmPrecursor && ppm>-params.ppmPrecursor) {
        pre.index=i;
        pre.monomass=p->monoMass;
        pre.maxZ = p->charge - 1;
        if (pre.maxZ < 1) pre.maxZ = 1;
        if (pre.maxZ > 3) pre.maxZ = 3;
        tp = s->getTopPeps(i);
        ps=i;
        break;
      }
    }

    //bufSize2[iIndex] = sizeof(double);
    size_t pepCount = ions[iIndex].pepCount;
    if(echo) cout << "pepCount = " << pepCount << endl;

    //TODO: FIX TO ONLY SCORE AROUND ONCE
    sScoreSet2* pScores = new sScoreSet2[pepCount];
    sScoreSet2* pScores2 = new sScoreSet2[pepCount];
    score9solo(s, ions[iIndex].vPeaks, pScores, pre.maxZ, iIndex);
    score9solo(s, ions[iIndex].vPeaksRev, pScores2, pre.maxZ, iIndex);

    vector<sDIndex> vTop;
    size_t minMods = 100;
    sc.simpleScore = 0;
    for (size_t b = 0; b < pepCount; b++) {
      if(echo) cout << b << " of " << pepCount << "\t" << mass << "\t" << ions[iIndex].pepMass[b] << endl;
      double sumScore=pScores[b].score+pScores2[b].score;
      if(ions[iIndex].pepMass[b]!=mass) continue;
      if (sumScore<=0) continue;
      if (sumScore > sc.simpleScore) {
        sc.simpleScore = (float)sumScore;
        vTop.clear();
        sDIndex di;
        di.a = b;
        di.b = pre.index;
        vTop.push_back(di);
        minMods = ions[iIndex].pepMods[b].mods.size();
      } else if (sumScore == sc.simpleScore) {
        sDIndex di;
        di.a = b;
        di.b = pre.index;
        vTop.push_back(di);
        if (ions[iIndex].pepMods[b].mods.size() < minMods) minMods = ions[iIndex].pepMods[b].mods.size();
      }
    }
    //cout << "simpleScore: " << sc.simpleScore << endl;
    if (sc.simpleScore == 0) {
      //peaks=NULL;
      delete [] pScores;
      delete[] pScores2;
      continue;
    }

    sc.simpleScore *= 0.005;
    double ev = 1000;
    Threading::LockMutex(mutexSpecScore[index[a]]);
    ev = spec->at(index[a]).computeE(sc.simpleScore, len);
    Threading::UnlockMutex(mutexSpecScore[index[a]]);
    
    for (size_t b = 0; b < vTop.size(); b++) {
      if (ions[iIndex].pepMods[vTop[b].a].mods.size() > minMods) continue; //skip modified peptides that are explained with fewer modifications
      sc.eVal = ev;
      sc.site = -1;
      sc.mass = mass;
      sc.massA = 0;
      sc.pep = pep1;
      sc.precursor = (char)pre.index;
      sc.mods->clear();
      for (size_t c = 0; c < ions[iIndex].pepMods[vTop[b].a].mods.size(); c++) {
        sc.mods->push_back(ions[iIndex].pepMods[vTop[b].a][c]);
      }

      Threading::LockMutex(mutexSingletScore[index[a]][ps]);
      tp->checkPeptideScore(sc);
      Threading::UnlockMutex(mutexSingletScore[index[a]][ps]);

      Threading::LockMutex(mutexSpecScore[index[a]]);
      spec->at(index[a]).checkScore(sc, iIndex);
      Threading::UnlockMutex(mutexSpecScore[index[a]]);
    }
    
    delete[] pScores;
    delete[] pScores2;
    //peaks=NULL;
  }

  p = NULL;
  tp = NULL;
  s=NULL;

}

//An alternative score uses the XCorr metric from the Comet algorithm
//This version allows for fast scoring when the cross-linked mass is added.
//float MAnalysis::magnumScoring(int specIndex, double modMass, int sIndex, int iIndex, int& match, int& conFrag, int z) { 
//
//  MSpectrum* s=spec->getSpectrum(specIndex);
//  MIonSet* ki=ions[iIndex].at(sIndex);
//
//  double dXcorr=0.0;
//  double invBinSize=s->getInvBinSize();
//  double binOffset=params.binOffset;
//  double dif;
//  double mz;
//
//  int ionCount=ions[iIndex].getIonCount();
//  int maxCharge=z;
//  if(maxCharge<1) maxCharge=s->getCharge();  
//
//  int i,j,k;
//  int key;
//  int pos;
//  int con;
//  match=0;
//  conFrag=0;
//
//  //Assign ion series
//  double***  ionSeries;
//  ionSeries=new double**[numIonSeries];
//  k=0;
//  if(params.ionSeries[0]) ionSeries[k++]=ki->aIons;
//  if(params.ionSeries[1]) ionSeries[k++]=ki->bIons;
//  if(params.ionSeries[2]) ionSeries[k++]=ki->cIons;
//  if(params.ionSeries[3]) ionSeries[k++]=ki->xIons;
//  if(params.ionSeries[4]) ionSeries[k++]=ki->yIons;
//  if(params.ionSeries[5]) ionSeries[k++]=ki->zIons;
//
//  //The number of fragment ion series to analyze is PrecursorCharge-1
//  //However, don't analyze past the 3+ series
//  if(maxCharge>4) maxCharge=4;
//
//  bool hardStop=false;
//
//  //Iterate all series
//  for(k=1;k<maxCharge;k++){
//
//    dif=modMass/k;
//
//    //Iterate through pfFastXcorrData
//    for(j=0;j<numIonSeries;j++){
//      con=0;
//
//      for(i=0;i<ionCount;i++){
//
//        //get key
//        if(ionSeries[j][k][i]<0) {
//          mz = params.binSize * (int)((dif-ionSeries[j][k][i])*invBinSize+binOffset);
//          if(mz<0) {
//            hardStop=true;
//            break;
//          }
//          key = (int)mz;
//          if (key >= s->kojakBins) {
//            if(con>conFrag) conFrag=con;
//            con=0;
//            break;
//          }
//          if (s->kojakSparseArray[key] == NULL) {
//            if (con>conFrag) conFrag = con;
//            con = 0;
//            continue;
//          }
//          pos = (int)((mz - key)*invBinSize);
//          dXcorr += s->kojakSparseArray[key][pos];
//          if (s->kojakSparseArray[key][pos]>5) {
//            match++;
//            con++;
//          } else {
//            if (con>conFrag) conFrag = con;
//            con = 0;
//          }
//
//        } else {
//          mz = params.binSize * (int)(ionSeries[j][k][i]*invBinSize+binOffset);
//          key = (int)mz;
//          if(key>=s->kojakBins) {
//            if (con>conFrag) conFrag = con;
//            con = 0;
//            break;
//          }
//          if(s->kojakSparseArray[key]==NULL) {
//            if (con>conFrag) conFrag = con;
//            con = 0;
//            continue;
//          }
//          pos = (int)((mz-key)*invBinSize);
//          dXcorr += s->kojakSparseArray[key][pos];
//          if (s->kojakSparseArray[key][pos]>5) {
//            match++;
//            con++;
//          } else {
//            if (con>conFrag) conFrag = con;
//            con = 0;
//          }
//        }
//      }
//      if(hardStop) break;
//      if(con>conFrag) conFrag=con;
//    }
//    if (hardStop) break;
//  }
//  
//
//  //Scale score appropriately
//  if(hardStop || dXcorr <= 0.0) dXcorr=0.0;
//  else dXcorr *= 0.005;
//
//  //Clean up memory
//  k = 0;
//  if (params.ionSeries[0]) ionSeries[k++] = NULL;
//  if (params.ionSeries[1]) ionSeries[k++] = NULL;
//  if (params.ionSeries[2]) ionSeries[k++] = NULL;
//  if (params.ionSeries[3]) ionSeries[k++] = NULL;
//  if (params.ionSeries[4]) ionSeries[k++] = NULL;
//  if (params.ionSeries[5]) ionSeries[k++] = NULL;
//
//  delete[] ionSeries;
//
//  return float(dXcorr);
//}

//An alternative score uses the XCorr metric from the Comet algorithm
//This version allows for fast scoring when the cross-linked mass is added.
char MAnalysis::magnumScoring2(MSpectrum* s, double mass) {
  //if(echo) cout << s->getScanNumber() << "\t" << mass << "\t" << s->kojakBins;
  double invBinSize = s->getInvBinSize();
  double mz = params.binSize * (int)(mass * invBinSize + params.binOffset);
  int key = (int)mz;
  //if(echo) cout << " key:" << key << " mz:"<<mz;
  if (key >= s->kojakBins) {
    //if(echo) cout << ": 0" << endl;
    return 0;
  }
  if (s->kojakSparseArray[key] == NULL) {
    //if(echo) cout << ": 0" << endl;
    return 0;
  }
  int pos = (int)((mz - key)*invBinSize);
  //if(echo) cout << ": " << (int)s->kojakSparseArray[key][pos] << endl;
  return s->kojakSparseArray[key][pos];
}


/*============================
  Utilities
============================*/
int MAnalysis::compareD(const void *p1, const void *p2){
  const double d1 = *(double *)p1;
  const double d2 = *(double *)p2;
  if(d1<d2) return -1;
  else if(d1>d2) return 1;
  else return 0;
}

