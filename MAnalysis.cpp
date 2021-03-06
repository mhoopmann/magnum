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

/*============================
  Constructors & Destructors
============================*/
MAnalysis::MAnalysis(mParams& p, MDatabase* d, MData* dat){
  unsigned int i;
  int j,k;
  
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

  int j;
  bool bt;
  vector<int> index;
  vector<mPepMod> mods;

  //char str[256];
  //db->getPeptideSeq(p->map->at(0).index,p->map->at(0).start,p->map->at(0).stop,str);
  //cout << str << "\t" << p->mass << endl;
  //Set the peptide, calc the ions, and score it against the spectra
  int len = (p->map->at(0).stop - p->map->at(0).start) + 1;
  ions[iIndex].setPeptide(&db->at(p->map->at(0).index).sequence[p->map->at(0).start],p->map->at(0).stop-p->map->at(0).start+1,p->mass,p->nTerm,p->cTerm);
  ions[iIndex].buildIons();
  ions[iIndex].modIonsRec2(0,-1,0,0,false);

  for(j=0;j<ions[iIndex].size();j++){
    bt=spec->getBoundaries2(ions[iIndex][j].mass,params.ppmPrecursor,index,scanBuffer[iIndex]);
    if(bt) scoreSpectra(index,j,len,ions[iIndex][j].difMass,pepIndex,-1,-1,-1,-1,iIndex);
   }

  if(p->xlSites==0) return true;

  //Crosslinked peptides must also search singlets with reciprocol mass on each lysine
  analyzeSinglets(*p,pepIndex,iIndex);

  return true;
}

bool MAnalysis::analyzeSinglets(mPeptide& pep, int index, int iIndex){

  int i;
  size_t j;
  int k;
  int len;
  double minMass;
  double maxMass;
  vector<int> scanIndex;
  string pepSeq;

  //get the peptide sequence
  db->getPeptideSeq(pep,pepSeq);

  //Set Mass boundaries
  minMass = pep.mass + params.minAdductMass;
  maxMass = pep.mass + params.maxAdductMass;
  minMass-=(minMass/1000000*params.ppmPrecursor);
  maxMass+=(maxMass/1000000*params.ppmPrecursor);

  //Find mod mass as difference between precursor and peptide
  
  len=(pep.map->at(0).stop-pep.map->at(0).start)+1;
  ions[iIndex].setPeptide(&db->at(pep.map->at(0).index).sequence[pep.map->at(0).start], len, pep.mass, pep.nTerm, pep.cTerm);
  
  //Iterate every adduct site
  for(k=0;k<len;k++){

    if (k == len - 1 && pep.cTerm){
      if (!adductSites['c'])  continue;
    } else if (k==len-1) {
      if(!adductSites['%']) continue;
    } else if (!adductSites[pepSeq[k]]) {
      if (k == 0 && pep.nTerm){
        if(!adductSites['n']) continue;
      } else if (k == 0) {
        if(!adductSites['$']) continue;
      } else {
        continue;
      }
    }
  
    //build fragment ions and score against all potential spectra
    ions[iIndex].reset();
    ions[iIndex].buildModIons(k);
    ions[iIndex].modIonsRec2(0,k,0,0,true);

    //iterate through all ion sets
    for(i=0;i<ions[iIndex].size();i++){

      //Iterate all spectra from (peptide mass + minimum mass) to (peptide mass + maximum mass)
      if (!spec->getBoundaries(minMass + ions[iIndex][i].difMass, maxMass + ions[iIndex][i].difMass, scanIndex, scanBuffer[iIndex])) continue;

      for(j=0;j<scanIndex.size();j++){
        //cout << "before score" << endl;
        scoreSingletSpectra(scanIndex[j], i, ions[iIndex][i].mass, len, index, (char)k, minMass + ions[iIndex][i].difMass, maxMass + ions[iIndex][i].difMass, iIndex);
        //cout << "after score" << endl;
      }
    }

  }

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
  }
  return true;
}

void MAnalysis::deallocateMemory(int threads){
  delete [] bKIonsManager;
  delete [] ions;
  for (int i = 0; i < threads; i++){
    delete[] scanBuffer[i];
  }
  delete[] scanBuffer;
}

//This function is way out of date. Particularly the mutexes and how to deal with multiple precursors.
void MAnalysis::scoreSingletSpectra(int index, int sIndex, double mass, int len, int pep, char k, double minMass, double maxMass, int iIndex){
  //cout << "scoreSingletSpectra()" << endl;
  mScoreCard sc;
  MIonSet* iset;
  mPepMod mod;
  double score=0;
  int i,j;
  int precI;
  int match;
  int conFrag;

  MSpectrum* s=spec->getSpectrum(index);
  mPrecursor* p;
  MTopPeps* tp;
  int sz=s->sizePrecursor();
  double topScore=0;
  int topMatch=0;
  int topConFrag=0;

  int code;
  for(i=0;i<sz;i++){
    p=s->getPrecursor2(i);
    if (p->corr<-4) code = 2;
    else if (p->corr<0)code = 3;
    else if (p->corr == 0)code = 2;
    else code = 1;
    if(code==1) break;

  }
  if(i==sz) {
    i=0;
    p = s->getPrecursor2(i);
  }
  if((p->monoMass-mass)<params.minAdductMass) score=0;  //this could be narrowed down to user-defined precursor tolerance.
  else if((p->monoMass-mass)>params.maxAdductMass) score=0;
  else score = magnumScoring(index, 0, sIndex, iIndex, match, conFrag, p->charge); //score peptide without open mod (i.e. scores peptide without localization)
  if(score>0){
    topScore = score;
    topMatch = match;
    topConFrag = conFrag;
    precI = i;
    sc.simpleScore = (float)score;
    sc.pep = pep;
    sc.mass = mass;
    sc.massA = p->monoMass - mass;
    sc.precursor = i;
    sc.site = -99;
    sc.mods->clear();
    iset = ions[iIndex].at(sIndex);
    if (iset->difMass != 0){
      for (j = 0; j<ions[iIndex].getIonCount(); j++) {
        if (iset->mods[j] != 0){
          if (j == 0 && iset->modNTerm) mod.term = true;
          else if (j == ions[iIndex].getIonCount() - 1 && iset->modCTerm) mod.term = true;
          else mod.term = false;
          mod.pos = (char)j;
          mod.mass = iset->mods[j];
          sc.mods->push_back(mod);
        }
      }
    }
  }

  for(i=0;i<sz;i++){
    p=s->getPrecursor2(i);
    //cout << i << " of " << sz << "\t" << p->monoMass << "\t" << minMass << "\t" << maxMass << "\t" << mass << endl;
    if(p->monoMass<minMass) continue;
    if(p->monoMass>maxMass) continue;
    if ((p->monoMass - mass)>params.maxAdductMass) continue;
    if ((p->monoMass - mass)<params.minAdductMass) continue;
    //cout << "Before magnumScoring" << endl;
    score=magnumScoring(index,p->monoMass-mass,sIndex,iIndex,match,conFrag,p->charge);  //open mod with localization
    //cout << score << endl;
    if(score==0) continue;
    else if(score>topScore) { //replace the previous peptide scores, if this version of the peptide scores better.
      topScore=score;
      topMatch=match;
      topConFrag=conFrag;
      precI=i;
      sc.simpleScore = (float)score;
      sc.pep = pep;
      sc.mass = mass;
      sc.massA = p->monoMass - mass;
      sc.precursor = i;
      sc.site = k;
      sc.mods->clear();
      iset = ions[iIndex].at(sIndex);
      if (iset->difMass != 0){
        for (j = 0; j<ions[iIndex].getIonCount(); j++) {
          if (iset->mods[j] != 0){
            if (j == 0 && iset->modNTerm) mod.term = true;
            else if (j == ions[iIndex].getIonCount() - 1 && iset->modCTerm) mod.term = true;
            else mod.term = false;
            mod.pos = (char)j;
            mod.mass = iset->mods[j];
            sc.mods->push_back(mod);
          }
        }
      }
    }
  }

  if(topScore>0){
    //cout << "Topper " << topScore << endl;
    double ev = 1000;
    Threading::LockMutex(mutexSpecScore[index]);
    ev = s->computeE(topScore, len);
    //cout << "DoneE " << ev << endl;
    Threading::UnlockMutex(mutexSpecScore[index]);
    sc.eVal=ev;
    sc.match=topMatch;
    sc.conFrag=topConFrag;

    tp = s->getTopPeps(precI);
    //cout << "GotTopPeps" << endl;
    Threading::LockMutex(mutexSingletScore[index][precI]);
    tp->checkPeptideScore(sc);
    //cout << "CheckPepScore" << endl;
    Threading::UnlockMutex(mutexSingletScore[index][precI]);
  
    Threading::LockMutex(mutexSpecScore[index]);
    s->checkScore(sc,iIndex);
    //cout << "CheckScore" << endl;
    Threading::UnlockMutex(mutexSpecScore[index]);
  }

}

void MAnalysis::scoreSpectra(vector<int>& index, int sIndex, int len, double modMass, int pep1, int pep2, int k1, int k2, int link, int iIndex){
  unsigned int a;
  int i,z,ps;
  mScoreCard sc;
  mPepMod mod;
  double mass = ions[iIndex][sIndex].mass;
  mPrecursor* p=NULL;
  MTopPeps* tp=NULL;

  //score spectra
  for(a=0;a<index.size();a++){
    
    //find the specific precursor mass in this spectrum to identify the charge state
    z=0;
    for (ps = 0; ps<spec->at(index[a]).sizePrecursor(); ps++){
      p = spec->at(index[a]).getPrecursor2(ps);
      double ppm = (p->monoMass - mass) / mass*1e6;
      if (ppm<params.ppmPrecursor && ppm>-params.ppmPrecursor){
        z = p->charge;
        tp = spec->at(index[a]).getTopPeps(ps);
        break;
      }
    }
    
    sc.simpleScore=magnumScoring(index[a],modMass,sIndex,iIndex,sc.match,sc.conFrag,z);
    if(sc.simpleScore==0) {

      //Threading::LockMutex(mutexSpecScore[index[a]]);
      //if (spec->at(index[a]).hp[iIndex].pepIndex != pep1){ //first score
      //  spec->at(index[a]).hp[iIndex].pepIndex = pep1;
      //  spec->at(index[a]).hp[iIndex].topScore = 0;
      //  spec->at(index[a]).histogram[0]++;
      //  spec->at(index[a]).histogramCount++;
      //} else {
      //  //do nothing, can't be a better score than what is already there.
      //}
      //Threading::UnlockMutex(mutexSpecScore[index[a]]);

      continue;
    }

    double ev = 1000;
    Threading::LockMutex(mutexSpecScore[index[a]]);
    ev = spec->at(index[a]).computeE(sc.simpleScore, len);
    Threading::UnlockMutex(mutexSpecScore[index[a]]);

    sc.eVal=ev;
    sc.mods->clear();
    sc.site=k1;
    sc.mass=mass;
    sc.massA=0;
    sc.pep=pep1;
    sc.precursor=(char)ps;
    if(ions[iIndex][sIndex].difMass!=0){
      for(i=0;i<ions[iIndex].getPeptideLen();i++) {
        if(ions[iIndex][sIndex].mods[i]!=0){
          if (i == 0 && ions[iIndex][sIndex].modNTerm) mod.term = true;
          else if (i == ions[iIndex].getIonCount() - 1 && ions[iIndex][sIndex].modCTerm) mod.term = true;
          else mod.term = false;
          mod.pos=(char)i;
          mod.mass=ions[iIndex][sIndex].mods[i];
          sc.mods->push_back(mod);
        }
      }
    }

    Threading::LockMutex(mutexSingletScore[index[a]][ps]);
    tp->checkPeptideScore(sc);
    Threading::UnlockMutex(mutexSingletScore[index[a]][ps]);

    Threading::LockMutex(mutexSpecScore[index[a]]);
    spec->at(index[a]).checkScore(sc,iIndex);
    Threading::UnlockMutex(mutexSpecScore[index[a]]);

  }

  p=NULL;
  tp=NULL;

}

//An alternative score uses the XCorr metric from the Comet algorithm
//This version allows for fast scoring when the cross-linked mass is added.
float MAnalysis::magnumScoring(int specIndex, double modMass, int sIndex, int iIndex, int& match, int& conFrag, int z) { 

  MSpectrum* s=spec->getSpectrum(specIndex);
  MIonSet* ki=ions[iIndex].at(sIndex);

  double dXcorr=0.0;
  double invBinSize=s->getInvBinSize();
  double binOffset=params.binOffset;
  double dif;
  double mz;

  int ionCount=ions[iIndex].getIonCount();
  int maxCharge=z;
  if(maxCharge<1) maxCharge=s->getCharge();  

  int i,j,k;
  int key;
  int pos;
  int con;
  match=0;
  conFrag=0;

  //Assign ion series
  double***  ionSeries;
  ionSeries=new double**[numIonSeries];
  k=0;
  if(params.ionSeries[0]) ionSeries[k++]=ki->aIons;
  if(params.ionSeries[1]) ionSeries[k++]=ki->bIons;
  if(params.ionSeries[2]) ionSeries[k++]=ki->cIons;
  if(params.ionSeries[3]) ionSeries[k++]=ki->xIons;
  if(params.ionSeries[4]) ionSeries[k++]=ki->yIons;
  if(params.ionSeries[5]) ionSeries[k++]=ki->zIons;

  //The number of fragment ion series to analyze is PrecursorCharge-1
  //However, don't analyze past the 3+ series
  if(maxCharge>4) maxCharge=4;

  bool hardStop=false;

  //Iterate all series
  for(k=1;k<maxCharge;k++){

    dif=modMass/k;

    //Iterate through pfFastXcorrData
    for(j=0;j<numIonSeries;j++){
      con=0;

      for(i=0;i<ionCount;i++){

        //get key
        if(ionSeries[j][k][i]<0) {
          mz = params.binSize * (int)((dif-ionSeries[j][k][i])*invBinSize+binOffset);
          if(mz<0) {
            hardStop=true;
            break;
          }
          key = (int)mz;
          if (key >= s->kojakBins) {
            if(con>conFrag) conFrag=con;
            con=0;
            break;
          }
          if (s->kojakSparseArray[key] == NULL) {
            if (con>conFrag) conFrag = con;
            con = 0;
            continue;
          }
          pos = (int)((mz - key)*invBinSize);
          dXcorr += s->kojakSparseArray[key][pos];
          if (s->kojakSparseArray[key][pos]>5) {
            match++;
            con++;
          } else {
            if (con>conFrag) conFrag = con;
            con = 0;
          }

        } else {
          mz = params.binSize * (int)(ionSeries[j][k][i]*invBinSize+binOffset);
          key = (int)mz;
          if(key>=s->kojakBins) {
            if (con>conFrag) conFrag = con;
            con = 0;
            break;
          }
          if(s->kojakSparseArray[key]==NULL) {
            if (con>conFrag) conFrag = con;
            con = 0;
            continue;
          }
          pos = (int)((mz-key)*invBinSize);
          dXcorr += s->kojakSparseArray[key][pos];
          if (s->kojakSparseArray[key][pos]>5) {
            match++;
            con++;
          } else {
            if (con>conFrag) conFrag = con;
            con = 0;
          }
        }
      }
      if(hardStop) break;
      if(con>conFrag) conFrag=con;
    }
    if (hardStop) break;
  }
  

  //Scale score appropriately
  if(hardStop || dXcorr <= 0.0) dXcorr=0.0;
  else dXcorr *= 0.005;

  //Clean up memory
  k = 0;
  if (params.ionSeries[0]) ionSeries[k++] = NULL;
  if (params.ionSeries[1]) ionSeries[k++] = NULL;
  if (params.ionSeries[2]) ionSeries[k++] = NULL;
  if (params.ionSeries[3]) ionSeries[k++] = NULL;
  if (params.ionSeries[4]) ionSeries[k++] = NULL;
  if (params.ionSeries[5]) ionSeries[k++] = NULL;

  delete[] ionSeries;

  return float(dXcorr);
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

