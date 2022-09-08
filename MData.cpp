/*
Copyright 2014, Michael R. Hoopmann, Institute for Systems Biology

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

#include "MData.h"

using namespace std;
using namespace MSToolkit;

Mutex MData::mutexMemoryPool;
bool* MData::memoryPool;
double** MData::tempRawData;
double** MData::tmpFastXcorrData;
float**  MData::fastXcorrData;
mPreprocessStruct** MData::preProcess;
mParams* MData::params;

MLog* MData::mlog;

deque<Spectrum*> MData::dMS1;
vector<Spectrum*> MData::vMS1Buffer;
Mutex MData::mutexLockMS1;
CHardklor2** MData::h;
CAveragine** MData::averagine;
CMercury8** MData::mercury;
CModelLibrary* MData::models;
Mutex* MData::mutexHardklor;
CHardklorSetting MData::hs;
bool* MData::bHardklor;

int MData::maxPrecursorMass;

/*============================
  Constructors
============================*/
MData::MData(){
  bScans=NULL;
  params=NULL;
  for(int i=0;i<128;i++) adductSite[i]=false;
  pepXMLindex=0;
}

MData::MData(mParams* p){
  bScans=NULL;
  params=p;
  size_t i;
  for(i=0;i<p->fMods.size();i++) aa.addFixedMod((char)p->fMods[i].index,p->fMods[i].mass);
  for (i = 0; i<128; i++) adductSite[i] = false;
  pepXMLindex = 0;
}

MData::~MData(){
  params=NULL;
  parObj=NULL;
  if(bScans!=NULL) delete[] bScans;
}


/*============================
  Operators
============================*/
MSpectrum& MData::operator [](const int& i){
  return *spec[i];
}


/*============================
  Functions
============================*/
bool MData::createDiag(FILE*& f){
  string fName=params->outFile;
  fName+=".diag.xml";
  f = fopen(fName.c_str(), "wt");
  if (f == NULL) {
    cout << "ERROR: Cannot open: " << fName << endl;
    return false;
  }
  return true;
}

NeoPepXMLParser* MData::createPepXML(string& str, MDatabase& db){
  str = params->outFile;
  str += ".pep.xml";
  FILE* f = fopen(str.c_str(), "wt");
  if (f == NULL) {
    cout << "ERROR: Cannot open: " << str << endl;
    return NULL;
  }
  fclose(f);
  
  NeoPepXMLParser* p = new NeoPepXMLParser();
  CnpxMSMSPipelineAnalysis pa;
  time_t timeNow;
  time(&timeNow);
  tm local_tm = *localtime(&timeNow);
  pa.date.date.year = local_tm.tm_year+1900;
  pa.date.date.month = local_tm.tm_mon + 1;
  pa.date.date.day = local_tm.tm_mday;
  pa.date.time.hour = local_tm.tm_hour;
  pa.date.time.minute = local_tm.tm_min;
  pa.date.time.second = local_tm.tm_sec;
  pa.summary_xml=str;

  CnpxMSMSRunSummary rs;
  rs.base_name=params->inFileNoExt; //params->inFile.substr(0,params->inFile.find_last_of('.'));
  rs.raw_data_type="raw";
  rs.raw_data=params->ext;

  CnpxSampleEnzyme se;
  se.name = params->enzymeName;
  CnpxSpecificity sp;
  mEnzymeRules enzyme = db.getEnzymeRules();
  for (int i = 65; i < 90; i++){
    if (enzyme.cutC[i]) sp.cut += (char)i;
    if (enzyme.exceptN[i]) sp.no_cut += (char)i;
  }
  if (sp.cut.size()>0) sp.sense = "C";
  else {
    sp.sense = "N";
    for (int i = 65; i < 90; i++){
      if (enzyme.cutN[i]) sp.cut += (char)i;
      if (enzyme.exceptC[i]) sp.no_cut += (char)i;
    }
  }
  se.specificity.push_back(sp);
  rs.sample_enzyme.push_back(se);
  
  CnpxSearchSummary ss;
  ss.base_name=params->outFile;
  ss.search_engine="Magnum";
  ss.search_engine_version=version;
  ss.precursor_mass_type="monoisotopic";
  ss.fragment_mass_type="monoisotopic";
  ss.search_id=1;
  for(size_t a=0;a<params->mods.size();a++){
    CnpxAminoAcidModification aam;
    aam.aminoacid=(char)params->mods[a].index;
    aam.massdiff=params->mods[a].mass;
    aam.mass = aam.massdiff + aa.getAAMass(params->mods[a].index);
    aam.variable="Y";
    ss.aminoacid_modification.push_back(aam);
  }
  for (size_t a = 0; a<params->fMods.size(); a++){
    CnpxAminoAcidModification aam;
    aam.aminoacid = (char)params->fMods[a].index;
    aam.massdiff = params->fMods[a].mass;
    aam.mass = aa.getAAMass(params->fMods[a].index);
    aam.variable = "N";
    ss.aminoacid_modification.push_back(aam);
  }

  for (size_t a = 0; a<parObj->xmlParams.size(); a++){
    CnpxParameter p;
    p.name=parObj->xmlParams[a].name;
    p.value=parObj->xmlParams[a].value;
    ss.parameter.push_back(p);
  }

  CnpxSearchDatabase sdb;
  sdb.local_path=params->dbPath;
  sdb.type="AA";
  
  ss.search_database.push_back(sdb);
  rs.search_summary.push_back(ss);
  pa.msms_run_summary.push_back(rs);
  p->msms_pipeline_analysis.push_back(pa);
  pepXMLindex=0;
  return p;
}

bool MData::createPercolator(FILE*& f, FILE*& f2){
  string sPercName = params->outFile;
  if(params->splitPercolator){
    string file1= sPercName+ ".standard.perc.txt";
    f = fopen(file1.c_str(), "wt");
    if (f == NULL) {
      cout << "ERROR: Cannot open: " << file1 << endl;
      return false;
    }
    fprintf(f, "SpecId\tLabel\tscannr\tMScore\tDScore");
    fprintf(f, "\tnLog10Eval\tCharge1\tCharge2\tCharge3");
    fprintf(f, "\tCharge4\tCharge5\tCharge6plus");
    fprintf(f, "\tMass\tPPM\tLength");
    fprintf(f, "\tPeptide\tProteins\n");

    string file2 = sPercName + ".open.perc.txt";
    f2 = fopen(file2.c_str(), "wt");
    if (f2 == NULL) {
      cout << "ERROR: Cannot open: " << file2 << endl;
      return false;
    }
    fprintf(f2, "SpecId\tLabel\tscannr\tMScore\tDScore");
    fprintf(f2, "\tnLog10Eval\tCharge1\tCharge2\tCharge3");
    fprintf(f2, "\tCharge4\tCharge5\tCharge6plus");
    fprintf(f2, "\tMass\tPPM\tLength");
    fprintf(f2, "\tPeptide\tProteins\n");


  } else {

    sPercName += ".perc.txt";
    f = fopen(sPercName.c_str(), "wt");
    if (f == NULL) {
      cout << "ERROR: Cannot open: " << sPercName << endl;
      return false;
    }
    fprintf(f, "SpecId\tLabel\tscannr\tMScore\tDScore");
    fprintf(f, "\tnLog10Eval\tCharge1\tCharge2\tCharge3");
    fprintf(f, "\tCharge4\tCharge5\tCharge6plus");
    fprintf(f, "\tMass\tPPM\tLength");
    fprintf(f, "\tPeptide\tProteins\n");

  }
  return true;
}

bool MData::createTXT(FILE*& f){
  string sTXTName = params->outFile;
  sTXTName += ".magnum.txt";
  f = fopen(sTXTName.c_str(), "wt");
  if (f == NULL) {
    cout << "ERROR: Cannot open: " << sTXTName << endl;
    return false;
  }
  fprintf(f,"ScanNumber\tPSMID\tRTsec\tSelectedMZ\tMonoisotopicMass");
  fprintf(f,"\tCharge\tPSMmass\tPPM\tEvalue");
  fprintf(f,"\tMscore\tDscore\tReporterIons");
  fprintf(f,"\tModifications\tOpenMod\tPeptide");
  fprintf(f,"\tPeptideA\tDecoy\tProteins\n");
  return true;
}

bool* MData::getAdductSites(){
  return &adductSite[0];
}

MSpectrum* MData::getSpectrum(const int& i){
  return spec[i];
}

void MData::initHardklor(){
  CHardklorVariant hv;
  hv.clear();
  hs.winSize = 10;
  hs.peptide = 4;
  hs.sn = 0;
  hs.depth = 3;
  hs.minCharge = 2;
  hs.maxCharge = 8;
  hs.algorithm = Version2;
  if (params->instrument == 1) hs.msType = FTICR;
  else hs.msType = OrbiTrap;
  hs.res400 = params->ms1Resolution;
  hs.corr = 0.875;
  hs.centroid = true;

  strcpy(hs.inFile, "PLTmp.ms1");

  hs.fileFormat = ms1;

  CHardklorVariant hkv;
  vector<CHardklorVariant> pepVariants;
  pepVariants.clear();
  pepVariants.push_back(hkv);

  averagine = new CAveragine*[params->threads]();
  mercury = new CMercury8*[params->threads]();
  h = new CHardklor2*[params->threads]();
  mutexHardklor = new Mutex[params->threads]();
  bHardklor = new bool[params->threads]();
  for (int a = 0; a<params->threads; a++){
    averagine[a] = new CAveragine(NULL, NULL);
    mercury[a] = new CMercury8(NULL);
    if (a == 0) models = new CModelLibrary(averagine[a], mercury[a]);
    h[a] = new CHardklor2(averagine[a], mercury[a], models);
    h[a]->Echo(false);
    h[a]->SetResultsToMemory(true);
    Threading::CreateMutex(&mutexHardklor[a]);
    bHardklor[a] = false;
  }
  models->eraseLibrary();
  models->buildLibrary(2, 8, pepVariants);

}

MSpectrum& MData::at(const int& i){
  return *spec[i];
}

//Places centroided peaks in bins, then finds average. Has potential for error when accuracy drifts into neighboring bins.
void MData::averageScansCentroid(vector<Spectrum*>& s, Spectrum& avg, double min, double max){
  unsigned int i;
  int j, k;
  double binWidth;
  double offset = -1.0;
  float* bin;
  float intensity;
  int binCount;
  int* pos;
  double  lowMZ = -1.0;
  mScanBin sb;
  vector<mScanBin> topList;

  avg.clear();

  //if vector is just one scan, then simply copy the peaks
  if (s.size() == 1){
    int index = findPeak(s[0], min);
    if (s[0]->at(index).mz<min) index++;
    while (index<s[0]->size() && s[0]->at(index).mz<max){
      avg.add(s[0]->at(index++));
    }
    return;
  }

  pos = new int[s.size()];

  //Set some really small bin width related to the mz region being summed.
  binWidth = 0.0001;

  binCount = (int)((max - min) / binWidth + 1);
  bin = new float[binCount];
  for (j = 0; j<binCount; j++) bin[j] = 0;

  //align all spectra to point closest to min and set offset
  for (i = 0; i<s.size(); i++)  {
    pos[i] = findPeak(s[i], min);
    if (s[i]->at(pos[i]).mz<min) pos[i]++;
    if (offset<0) offset = s[i]->at(pos[i]).mz;
    else if (s[i]->at(pos[i]).mz<offset) offset = s[i]->at(pos[i]).mz;
  }

  //Iterate all spectra and add peaks to bins
  for (i = 0; i<s.size(); i++) {
    while (pos[i]<s[i]->size() && s[i]->at(pos[i]).mz<max){
      j = (int)((s[i]->at(pos[i]).mz - offset) / binWidth);
      if (j<0){
        pos[i]++;
        continue;
      }
      if (j >= binCount) break;
      bin[j] += s[i]->at(pos[i]).intensity;
      pos[i]++;
    }
  }

  //Unsure of current efficiency. Finds bin of tallest peak, then combines with neighboring bins
  //to produce an average. Thus summing of neighboring bins allows flexibility when dealing with mass accuracy drift.
  //Using larger bins has the same effect, but perhaps fewer significant digits in the final result.
  for (j = 0; j<binCount; j++){
    if (bin[j]>0) {
      sb.index = j;
      sb.intensity = bin[j];
      topList.push_back(sb);
    }
  }
  if (topList.size()>0) {
    qsort(&topList[0], topList.size(), sizeof(mScanBin), compareScanBinRev2);
    for (i = 0; i<topList.size(); i++){
      if (bin[topList[i].index] == 0) continue;
      intensity = 0;
      j = topList[i].index - 50; //This will sum bins within 0.05 m/z... might be too wide...
      k = topList[i].index + 51;
      if (j<0) j = 0;
      if (k>binCount) k = binCount;
      for (j = j; j<k; j++){
        intensity += bin[j];
        bin[j] = 0;
      }
      intensity /= s.size();
      avg.add(offset + topList[i].index*binWidth, intensity);
    }
    avg.sortMZ();
  }

  //clean up memory
  delete[] bin;
  delete[] pos;

}

void MData::diagSinglet(){
  int oddCount = 0;
  double maxScore;
  int bigCount = 0;
  int twoCount = 0;
  double bigScore = 0;
  double bigMass;
  double unknownMass;
  for (size_t b = 0; b<spec.size(); b++){
    //if(spec[b].getScoreCard(0).simpleScore==0) continue;
    maxScore = 0;
    unknownMass = 0;
    for (int q = 0; q<spec[b]->sizePrecursor(); q++){
      if (spec[b]->getTopPeps(q)->peptideCount == 0) continue;
      if (spec[b]->getTopPeps(q)->peptideFirst->simpleScore>spec[b]->getScoreCard(0).simpleScore){
        if (spec[b]->getTopPeps(q)->peptideFirst->simpleScore>maxScore) {
          maxScore = spec[b]->getTopPeps(q)->peptideFirst->simpleScore;
          unknownMass = spec[b]->getPrecursor(q).monoMass - spec[b]->getTopPeps(q)->peptideFirst->mass;
        }
      }
    }
    if (maxScore>0) oddCount++;
    if (maxScore>3) bigCount++;
    if (maxScore>3 && spec[b]->getScoreCard(0).simpleScore>0 && maxScore>spec[b]->getScoreCard(0).simpleScore * 2) twoCount++;
    if (maxScore>bigScore) {
      bigScore = maxScore;
      bigMass = unknownMass;
    }
  }
  cout << " Diagnostics:" << endl;
  cout << "  Odd counts " << oddCount << " of " << spec.size() << endl;
  cout << "  Big counts " << bigCount << " of " << oddCount << endl;
  cout << "  Two fold " << twoCount << " of " << bigCount << endl;
  cout << "  Biggest Score: " << bigScore << "  Mass: " << bigMass << endl;
  cout << endl;
}

void MData::exportPercolator(FILE*& f, vector<mResults>& r){
  for (size_t a = 0; a<r.size(); a++){
    if (r[a].decoy) {
      fprintf(f, "D-%d-%.2f", r[a].scanNumber, r[a].rTimeSec / 60);
      if(a>0) fprintf(f,"-%d",(int)a+1);
      fprintf(f, "\t-1");
    } else {
      fprintf(f, "T-%d-%.2f", r[a].scanNumber, r[a].rTimeSec / 60);
      if (a>0) fprintf(f, "-%d", (int)a + 1);
      fprintf(f, "\t1");
    }
    fprintf(f, "\t%d", r[a].scanNumber);
    fprintf(f, "\t%.4lf", r[a].scoreMagnum);
    fprintf(f, "\t%.4lf", r[a].scoreDelta);
    fprintf(f, "\t%.6lf", -log10(r[a].eValue));
    int z=r[a].charge;
    if(z>6) z=6;
    for(int b=1;b<z;b++) fprintf(f,"\t0");
    fprintf(f,"\t1");
    for (int b = z + 1; b<7; b++) fprintf(f, "\t0");
    fprintf(f, "\t%.4lf", r[a].psmMass);
    fprintf(f, "\t%.4lf", r[a].ppm);
    fprintf(f, "\t%d", (int)r[a].peptide.size());
    fprintf(f, "\t-.%s.-", r[a].modPeptide[0].c_str());
    for(size_t b=0;b<r[a].proteins.size();b++){
      fprintf(f, "\t%s", r[a].proteins[b].protein.c_str());
    }
    fprintf(f, "\n");
  }
}

void MData::exportTXT(FILE*& f, vector<mResults>& r){
  for(size_t a=0;a<r.size();a++){
    fprintf(f,"%d",r[a].scanNumber);
    fprintf(f, "\t%d", r[a].psmID);
    fprintf(f, "\t%.2f", r[a].rTimeSec);
    fprintf(f,"\t%.6lf",r[a].selectedMZ);
    fprintf(f, "\t%.6lf", r[a].monoMass);
    fprintf(f,"\t%d",r[a].charge);
    fprintf(f, "\t%.6lf", r[a].psmMass);
    fprintf(f, "\t%.4lf", r[a].ppm);
    fprintf(f, "\t%.4e", r[a].eValue);
    fprintf(f, "\t%.2lf", r[a].scoreMagnum);
    fprintf(f, "\t%.2lf", r[a].scoreDelta);
    if(r[a].rIon.size()==0) fprintf(f,"\t");
    else {
      fprintf(f,"\t%.3lf",r[a].rIon[0]);
      for(size_t b=1;b<r[a].rIon.size();b++) fprintf(f,";%.3lf",r[a].rIon[b]);
    }
    fprintf(f, "\t%s", r[a].mods.c_str());
    fprintf(f, "\t%s", r[a].openMod.c_str());
    fprintf(f, "\t%s", r[a].peptide.c_str());
    fprintf(f, "\t%s", r[a].modPeptide[0].c_str());
    if(r[a].decoy) fprintf(f, "\t1");
    else fprintf(f,"\t0");
    fprintf(f, "\t%s", r[a].proteins[0].protein.c_str());
    for(size_t b=1;b<r[a].proteins.size();b++){
      fprintf(f, ";%s", r[a].proteins[b].protein.c_str());
    }
    fprintf(f,"\n");
  }
}

void MData::exportPepXML(NeoPepXMLParser*& p, vector<mResults>& r){
  if(r.size()==0) return;
  
  CnpxSpectrumQuery s;
  char str[256];
  sprintf(str,"%s.%d.%d.%d",params->msBase.c_str(),r[0].scanNumber,r[0].scanNumber,r[0].charge);
  s.spectrum=str;
  s.start_scan=r[0].scanNumber;
  s.end_scan=r[0].scanNumber;
  s.precursor_neutral_mass=r[0].monoMass;
  s.assumed_charge=r[0].charge;
  s.retention_time_sec=r[0].rTimeSec;
  s.index=++pepXMLindex;
  
  CnpxSearchResult sr;
  for (size_t a = 0; a<r.size(); a++){
    for(size_t x=0;x<r[a].vMods.size();x++){ //to export EVERY alternate arrangement of mods for this PSM
      CnpxSearchHit sh;
      sh.hit_rank=1;
      sh.peptide=r[a].peptide;
      sh.protein=r[a].proteins[0].protein;
      sh.peptide_prev_aa=r[a].proteins[0].prevAA;
      sh.peptide_next_aa=r[a].proteins[0].nextAA;
      sh.peptide_start_pos = r[a].proteins[0].startPos;
      for(size_t b=1;b<r[a].proteins.size();b++){
        CnpxAlternativeProtein ap;
        ap.protein=r[a].proteins[b].protein;
        ap.peptide_prev_aa = r[a].proteins[b].prevAA;
        ap.peptide_next_aa = r[a].proteins[b].nextAA;
        ap.peptide_start_pos = r[a].proteins[b].startPos;
        sh.alternative_protein.push_back(ap);
      }
      sh.num_tot_proteins=(int)r[a].proteins.size();
      sh.calc_neutral_pep_mass=r[a].psmMass;

      //Add scores
      sprintf(str,"%.4lf",r[a].scoreMagnum);
      sh.addSearchScore("Mscore",string(str));
      sprintf(str, "%.4lf", r[a].scoreDelta);
      sh.addSearchScore("Dscore", string(str));
      sprintf(str, "%.4lf", r[a].ppm);
      sh.addSearchScore("ppm_error", string(str));
      sprintf(str, "%.3e", r[a].eValue);
      sh.addSearchScore("Evalue", string(str));
      for(size_t b=0;b<r[a].rIon.size(); b++){
        sprintf(str, "%.4lf", r[a].rIon[b]);
        sh.addSearchScore("reporter_ion", string(str));
      }

      //check modifications
      CnpxModificationInfo mi;
      //check n-term
      for(size_t b=0;b<r[a].vMods[x].mods.size();b++){
        if (r[a].vMods[x].mods[b].term && r[a].vMods[x].mods[b].pos == 0 && !r[a].vMods[x].mods[b].adduct){
          sprintf(str, "n[%d]", (int)(r[a].vMods[x].mods[b].mass + 0.5));
          mi.modified_peptide+=str;
          mi.mod_nterm_mass = r[a].vMods[x].mods[b].mass + 1.00782503;
          break;
        }
      }
      //check all amino acids
      for(size_t b=0;b<r[a].peptide.size();b++){
        mi.modified_peptide+=r[a].peptide[b];
        //check fixed mods
        if(aa.getFixedModMass(r[a].peptide[b])!=0){
          char c = r[a].peptide[b];
          sprintf(str, "[%d]", (int)(aa.getAAMass(c) + 0.5));
          mi.modified_peptide += str;
          CnpxModAminoAcidMass mam;
          mam.position = (int)b + 1;
          mam.mass = aa.getAAMass(r[a].peptide[b]);
          mam.staticMass = aa.getFixedModMass(c);
          mam.source = "param";
          mi.mod_aminoacid_mass.push_back(mam);
        }
        //check variable mods
        for (size_t c = 0; c<r[a].vMods[x].mods.size(); c++){
          if (r[a].vMods[x].mods[c].pos == b && !r[a].vMods[x].mods[c].term){
            sprintf(str, "[%d]", (int)(r[a].vMods[x].mods[c].mass + 0.5));
            mi.modified_peptide += str;
            CnpxModAminoAcidMass mam;
            mam.position = (int)r[a].vMods[x].mods[c].pos + 1;
            mam.mass = aa.getAAMass(r[a].peptide[b]) + r[a].vMods[x].mods[c].mass;
            if (r[a].vMods[x].mods[c].variable) mam.variable = r[a].vMods[x].mods[c].mass;
            else mam.staticMass = r[a].vMods[x].mods[c].mass;
            if (r[a].vMods[x].mods[c].adduct) mam.source = "adduct";
            else mam.source="param";
            mi.mod_aminoacid_mass.push_back(mam);
            break;
          }
        }
      }
      //check c-term
      for (size_t b = 0; b<r[a].vMods[x].mods.size(); b++){
        if (r[a].vMods[x].mods[b].term && r[a].vMods[x].mods[b].pos > 0 && !r[a].vMods[x].mods[b].adduct){
          sprintf(str, "c[%d]", (int)(r[a].vMods[x].mods[b].mass + 0.5));
          mi.modified_peptide += str;
          mi.mod_cterm_mass = r[a].vMods[x].mods[b].mass + 17.00273963;
          break;
        }
      }
      if(mi.mod_aminoacid_mass.size()>0 || mi.mod_cterm_mass!=0 || mi.mod_nterm_mass!=0) sh.modification_info.push_back(mi);

      //check for open modifications (adducts)
      for (size_t b = 0; b<r[a].vMods[x].mods.size(); b++){
        if (r[a].vMods[x].mods[b].adduct && r[a].vMods[x].mods[b].term){
          sh.massdiff = r[a].vMods[x].mods[b].mass;
        }
      }

      sr.search_hit.push_back(sh);
    }
  }
  s.search_result.push_back(sr);

  p->msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.push_back(s);
}

//Returns closet mz value to desired point
int MData::findPeak(Spectrum* s, double mass){
  int sz = s->size();
  int lower = 0;
  int mid = sz / 2;
  int upper = sz;


  //binary search to closest mass
  while (s->at(mid).mz != mass){
    if (lower >= upper) break;
    if (mass<s->at(mid).mz){
      upper = mid - 1;
      mid = (lower + upper) / 2;
    } else {
      lower = mid + 1;
      mid = (lower + upper) / 2;
    }
    if (mid == sz) {
      mid--;
      break;
    }
  }

  //Check that mass is closest
  if (mid>0 && fabs(s->at(mid - 1).mz - mass)<fabs(s->at(mid).mz - mass)) return mid - 1;
  if (mid<s->size() - 1 && fabs(s->at(mid + 1).mz - mass)<fabs(s->at(mid).mz - mass)) return mid + 1;
  return mid;

}

//Returns point within precision or -1 if doesn't exist
int MData::findPeak(Spectrum* s, double mass, double prec){
  int sz = s->size();
  int lower = 0;
  int mid = sz / 2;
  int upper = sz;

  double minMass = mass - (mass / 1000000 * prec);
  double maxMass = mass + (mass / 1000000 * prec);

  //binary search to closest mass
  while (s->at(mid).mz<minMass || s->at(mid).mz>maxMass){
    if (lower >= upper) break;
    if (mass<s->at(mid).mz){
      upper = mid - 1;
      mid = (lower + upper) / 2;
    } else {
      lower = mid + 1;
      mid = (lower + upper) / 2;
    }
    if (mid == sz) {
      mid--;
      break;
    }
  }

  //Check that mass is correct
  if (s->at(mid).mz>minMass && s->at(mid).mz<maxMass) return mid;

  return -1;
}

void MData::formatMS2(MSToolkit::Spectrum* s, MSpectrum* pls){
  char nStr[256];
  string sStr;

  pls->setRTime(s->getRTime());
  pls->setScanNumber(s->getScanNumber());
  s->getNativeID(nStr, 256);
  sStr = nStr;
  pls->setNativeID(sStr);

  bool doCentroid = false;
  switch (s->getCentroidStatus()){
  case 0:
    if (params->ms2Centroid) {
      char tmpStr[256];
      sprintf(tmpStr, "Magnum parameter indicates MS/MS data are centroid, but spectrum %d labeled as profile.", s->getScanNumber());
      mlog->addError(string(tmpStr));
    } else doCentroid = true;
    break;
  case 1:
    if (!params->ms2Centroid) {
      mlog->addWarning(0, "Spectrum is labeled as centroid, but Kojak parameter indicates data are profile. Ignoring Kojak parameter.");
    }
    break;
  default:
    if (!params->ms2Centroid) doCentroid = true;
    break;
  }

  //If not centroided, do so now.
  int totalPeaks = 0;
  if (doCentroid){
    centroid(s, pls, params->ms2Resolution, params->instrument);
    totalPeaks += pls->size();
  } else {
    mSpecPoint sp;
    for (int i = 0; i<s->size(); i++){
      sp.mass = s->at(i).mz;
      sp.intensity = s->at(i).intensity;
      pls->addPoint(sp);
    }
  }

  ////remove precursor if requested
  //if (params->removePrecursor>0){
  //  double pMin = s->getMZ() - params->removePrecursor;
  //  double pMax = s->getMZ() + params->removePrecursor;
  //  for (int i = 0; i<pls->size(); i++){
  //    if ((*pls)[i].mass>pMin && (*pls)[i].mass<pMax) (*pls)[i].intensity = 0;
  //  }
  //}

  //Collapse the isotope peaks
  int collapsedPeaks = 0;
  if (params->specProcess == 1 && pls->size()>1) {
    collapseSpectrum(*pls);
    collapsedPeaks += pls->size();
  }

  //If user limits number of peaks to analyze, sort by intensity and take top N
  if (params->maxPeaks>0 && pls->size()>params->maxPeaks){
    if (pls->size()>1) pls->sortIntensityRev();
    vector<mSpecPoint> v;
    for (int i = 0; i<params->maxPeaks; i++) v.push_back((*pls)[i]);
    pls->clear();
    for (int i = 0; i<params->maxPeaks; i++) pls->addPoint(v[i]);
    pls->sortMZ();
    pls->setMaxIntensity(v[0].intensity);
  } else {
    float max = 0;
    for (int i = 0; i<pls->size(); i++){
      if ((*pls)[i].intensity>max) max = (*pls)[i].intensity;
    }
    pls->setMaxIntensity(max);
  }

  //Get any additional information user requested
  pls->setCharge(s->getCharge());
  pls->setMZ(s->getMZ());
  if (params->preferPrecursor>0){
    if (s->getMonoMZ()>0 && s->getCharge()>0){
      mPrecursor pre;
      pre.monoMass = s->getMonoMZ()*s->getCharge() - s->getCharge()*1.007276466;
      pre.charge = s->getCharge();
      pre.corr = 0;
      pls->addPrecursor(pre, params->topCount);
      for (int px = 1; px <= params->isotopeError; px++){
        if (px == 4) break;
        pre.monoMass -= 1.00335483;
        pre.corr -= 0.1;
        pls->addPrecursor(pre, params->topCount);
      }
      pls->setInstrumentPrecursor(true);
    }
  }

}

bool MData::getBoundaries(double mass1, double mass2, vector<int>& index, bool* buffer){
  int sz=(int)massList.size();

  if(mass1>massList[sz-1].mass) return false;

  int lower=0;
  int mid=sz/2;
  int upper=sz;
  int i;
  int low;
  int high;

  vector<int> v;

  //binary search to closest mass
  while(massList[mid].mass!=mass1){
		if(lower>=upper) break;
    if(mass1<massList[mid].mass){
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
		if(mid==sz) {
			mid--;
			break;
		}
	}

  //Adjust if needed
  if(massList[mid].mass<mass1) mid++;
  if(mid>0 && massList[mid-1].mass>mass1) mid--;
  if(massList[mid].mass>mass2) return false;
  if(mid==sz) return false;
  low=mid;

  //binary search to next mass
  lower=0;
  mid=sz/2;
  upper=sz;
  while(massList[mid].mass!=mass2){
		if(lower>=upper) break;
    if(mass2<massList[mid].mass){
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
		if(mid==sz) {
			mid--;
			break;
		}
  }

  //Adjust if needed
  if(massList[mid].mass>mass2) mid--;
  if(mid<sz-1 && massList[mid+1].mass<mass2) mid++;
  if(mid<0) return false;
  high=mid;

  sz=(int)spec.size();
  memset(buffer,false,sz);
  for (i = low; i <= high; i++) buffer[massList[i].index]=true;
  index.clear();
  for (i = 0; i < sz; i++) {
    if(buffer[i]) index.push_back(i);
  }
  return true;

}

//Get the list of spectrum array indexes to search based on desired mass
bool MData::getBoundaries2(double mass, double prec, vector<int>& index, bool* buffer){
  int sz=(int)massList.size();
  int lower=0;
  int mid=sz/2;
  int upper=sz;
	int i;

  vector<int> v;

  double minMass = mass - (mass/1000000*prec);
  double maxMass = mass + (mass/1000000*prec);

  //binary search to closest mass
  while(massList[mid].mass<minMass || massList[mid].mass>maxMass){
		if(lower>=upper) break;
    if(mass<massList[mid].mass){
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
		if(mid==sz) {
			mid--;
			break;
		}
	}

	//Check that mass is correct
	if(massList[mid].mass<minMass || massList[mid].mass>maxMass) return false;

  v.push_back(massList[mid].index);

	//check left 
  i=mid;
	while(i>0){
		i--;
    if(massList[i].mass<minMass) break;
    v.push_back(massList[i].index);
	}

	//check right
	i=mid;
	while(i<(sz-1)){
		i++;
    if(massList[i].mass>maxMass) break;
    v.push_back(massList[i].index);
	}

  sz = (int)spec.size();
  memset(buffer, false, sz);
  for (i = 0; i < (int)v.size(); i++) buffer[v[i]] = true;
  index.clear();
  for (i = 0; i < sz; i++) {
    if (buffer[i]) index.push_back(i);
  }
  return true;

  //Sort indexes and copy to final array, removing duplicates.
  //This may be a potentially slow step and should be profiled
  qsort(&v[0],v.size(),sizeof(int),compareInt);
  index.clear();
  index.push_back(v[0]);
  mid=(int)v.size();
  for(i=1;i<mid;i++){
    if(v[i]!=v[i-1]) index.push_back(v[i]);
  }

	return true;

}

double MData::getMaxMass(){
  if(massList.size()==0) return 0;
  else return massList[massList.size()-1].mass;
}

double MData::getMinMass(){
  if(massList.size()==0) return 0;
  else return massList[0].mass;
}

void MData::memoryAllocate(){
  //find largest possible array for a spectrum
  int threads = params->threads;
  double td = params->maxPepMass + params->maxAdductMass + 1;
  maxPrecursorMass = (int)td;
  int xCorrArraySize = (int)((params->maxPepMass + params->maxAdductMass + 100.0) / params->binSize);

  //Mark all arrays as available
  memoryPool = new bool[threads];
  for (int a = 0; a<threads; a++) memoryPool[a] = false;

  //Allocate arrays
  tempRawData = new double*[threads]();
  for (int a = 0; a<threads; a++) tempRawData[a] = new double[xCorrArraySize]();

  tmpFastXcorrData = new double*[threads]();
  for (int a = 0; a<threads; a++) tmpFastXcorrData[a] = new double[xCorrArraySize]();

  fastXcorrData = new float*[threads]();
  for (int a = 0; a<threads; a++) fastXcorrData[a] = new float[xCorrArraySize]();

  preProcess = new mPreprocessStruct*[threads]();
  for (int a = 0; a<threads; a++) {
    preProcess[a] = new mPreprocessStruct();
    preProcess[a]->pdCorrelationData = new mSpecPoint[xCorrArraySize]();
    preProcess[a]->iMaxXCorrArraySize=xCorrArraySize;
  }

  //Create mutex
  Threading::CreateMutex(&mutexMemoryPool);
}

void MData::memoryFree(){
  delete[] memoryPool;
  for (int a = 0; a<params->threads; a++){
    delete[] tempRawData[a];
    delete[] tmpFastXcorrData[a];
    delete[] fastXcorrData[a];
    delete[] preProcess[a]->pdCorrelationData;
    delete preProcess[a];
  }
  delete[] tempRawData;
  delete[] tmpFastXcorrData;
  delete[] fastXcorrData;
  delete[] preProcess;

  //Destroy mutexes
  Threading::DestroyMutex(mutexMemoryPool);
}

void MData::outputDiagnostics(FILE* f, MSpectrum& s, MDatabase& db){
  size_t i,x;
  int j,k;
  int code;
  char strs[256];
  char st[32];
  string pep1,tmp;
  mPeptide pep;
  mPrecursor* p;
  MTopPeps* tp;
  mScoreCard* sc;
  mScoreCard psm;
  
  fprintf(f, " <scan id=\"%d\">\n", s.getScanNumber());
  fprintf(f, "  <precursor_list size=\"%d\">\n", s.sizePrecursor());
  
  for (j = 0; j<s.sizePrecursor(); j++) {
    p=s.getPrecursor2(j);
    if(p->corr<-4) code=2;
    else if (p->corr<0)code=3;
    else if(p->corr==0)code=2;
    else code=1;
    fprintf(f, "   <precursor mass=\"%.4lf\" charge=\"%d\" type=\"%d\" hk_corr=\"%.4lf\">\n", p->monoMass, p->charge, code, p->corr);

    tp = s.getTopPeps(j);
    sc = tp->peptideFirst;
    k=1;
    while (sc != NULL){
      fprintf(f,"    <peptide rank=\"%d\" sequence=\"",k++);
      db.getPeptideSeq(db.getPeptideList()->at(sc->pep).map->at(0).index, db.getPeptideList()->at(sc->pep).map->at(0).start, db.getPeptideList()->at(sc->pep).map->at(0).stop, strs);
      for (i = 0; i<strlen(strs); i++){
        fprintf(f, "%c", strs[i]);
        for (x = 0; x<sc->mods->size(); x++){
          if (sc->mods->at(x).pos == char(i)) fprintf(f, "[%.2lf]", sc->mods->at(x).mass);
        }
        if (char(i) == sc->site) fprintf(f, "[x]");
      }
      fprintf(f, "\" adduct_site=\"%d\" evalue=\"%.3e\" score=\"%.4lf\" mass=\"%.4lf\" matches=\"%d\" con_fragments=\"%d\" additional_mass=\"%.4lf\"/>\n", (int)sc->site+1, sc->eVal, sc->simpleScore, sc->mass, sc->match, sc->conFrag, sc->massA);
      sc = sc->next;
    }
    fprintf(f,"   </precursor>\n");
  }
  fprintf(f,"  </precursor_list>\n");

  k=0;
  for(j=0;j<20;j++){
    if(s.getScoreCard(j).simpleScore>0) k++;
    else break;
  }
  fprintf(f,"  <results_list size=\"%d\">\n",k);
  for (j = 0; j<k; j++){
    fprintf(f,"   <result rank=\"%d\" ",j+1);
    psm = s.getScoreCard(j);
    pep = db.getPeptide(psm.pep);
    db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, strs);
    pep1.clear();
    if (pep.nTerm && aa.getFixedModMass('$') != 0) {
      sprintf(st, "[%.2lf]", aa.getFixedModMass('$'));
      pep1+=st;
    }
    for (i = 0; i<strlen(strs); i++){
      pep1+=strs[i];
      for (x = 0; x<psm.mods->size(); x++){
        if (psm.mods->at(x).pos == (char)i) {
          sprintf(st, "[%.2lf]", psm.mods->at(x).mass);
          pep1+=st;
        }
      }
      if ((int)i == psm.site) pep1+="[x]";
    }
    if (pep.cTerm && aa.getFixedModMass('%') != 0) {
      sprintf(st, "[%.2lf]", aa.getFixedModMass('%'));
      pep1+=st;
    }

    fprintf(f, "sequence=\"%s\" score=\"%.4lf\" evalue=\"%.3e\" mass=\"%.4lf\"",&pep1[0], psm.simpleScore, psm.eVal, psm.mass);
    fprintf(f, " additional_mass=\"%.4lf\" precursor=\"%d\"/>\n", psm.massA, (int)psm.precursor);
  }
  fprintf(f, "  </results_list>\n");
  /*
  fprintf(f,"  <histogram count=\"%d\">\n",s.histogramCount);
  for(j=0;j<HISTOSZ;j++){
    fprintf(f,"   <bin id=\"%d\" value=\"%d\"/>\n",j,s.histogram[j]);
  }
  fprintf(f,"  </histogram>\n");
  */
  fprintf(f, " </scan>\n");
 
}

bool MData::outputResults(MDatabase& db){
  FILE* fTXT=NULL;
  FILE* fPerc=NULL;
  FILE* fPerc2 = NULL;
  FILE* fDiag=NULL;
  string fXML;
  NeoPepXMLParser* pepxml;

  //Create output files
  if(!createTXT(fTXT)) return false;
  if(params->diag.size()>0) {
    if(!createDiag(fDiag)) return false;
  }
  if(params->exportPepXML){
    pepxml=createPepXML(fXML,db);
    if(pepxml==NULL) return false;
  }
  if(params->exportPercolator){
    if(!createPercolator(fPerc,fPerc2)) return false;
  }
  
  //iterate over all spectra to create exportable results
  //Output top score for each spectrum
  //Must iterate through all possible precursors for that spectrum
  for (size_t a = 0; a<spec.size(); a++) {

    //** temporary
    //spec[a].exportHisto();
    //**

    //Check if we need to output diagnostic information
    bool bDiag = false;
    if (params->diag.size()>0){
      if (params->diag[0] == -1) bDiag = true;
      else {
        for (size_t b = 0; b<params->diag.size(); b++){
          if (spec[a]->getScanNumber() == params->diag[b]){
            bDiag = true;
            break;
          }
        }
      }
    }
    if (bDiag) outputDiagnostics(fDiag, *spec[a], db);

    //Default score card
    vector<mResults> vRes;
    vector<mScoreCard3> shorts;
    spec[a]->shortResults2(shorts);
    //condenseResults(shorts);
    if (shorts.size() > 0) { //figure out how many results need reporting
      double topScore=shorts[0].eVal;
      size_t scoreIndex=0;
      while(shorts[0].simpleScore>0 && shorts[scoreIndex].eVal==topScore){
        mResults res;

        res.psmID=(int)scoreIndex+1;
        processSpectrumInfo(*spec[a],res);
        processPSM(*spec[a],shorts[scoreIndex],res);

        //delta score
        size_t n = scoreIndex + 1;
        while (n<shorts.size()){
          if (shorts[n].simpleScore == 0) break;
          if (shorts[n].eVal == topScore) {
            n++;
            continue;
          }
          if (shorts[n].precursor != shorts[scoreIndex].precursor) {
            n++;
            continue;
          }
          res.scoreDelta=res.scoreMagnum-shorts[n].simpleScore;
          break;
        }
        if(res.scoreDelta==0 && n>0) res.scoreDelta=res.scoreMagnum-shorts[n-1].simpleScore;

        //peptide sequence
        mPeptide pep = db.getPeptide(shorts[scoreIndex].pep);
        db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, res.peptide);
        for(size_t b=0;b<shorts[scoreIndex].mSet.size();b++){
          string spep=processPeptide(pep, shorts[scoreIndex].mSet[b].mods, (int)shorts[scoreIndex].aSites[b], shorts[scoreIndex].massA, db);
          res.modPeptide.push_back(spep);
        }

        //Determine if target or decoy
        bool bTarget = false;
        for (size_t b = 0; b<pep.map->size(); b++) if (db[pep.map->at(b).index].name.find(params->decoy) == string::npos) bTarget = true;
        if (bTarget) res.decoy = false;
        else res.decoy = true;

        //proteins
        for (size_t b = 0; b<pep.map->size(); b++){
          mProtRes pr;
          pr.protein = db[pep.map->at(b).index].name;
          if(pep.map->at(b).start==0) pr.prevAA='-';
          else pr.prevAA = db[pep.map->at(b).index].sequence[pep.map->at(b).start-1];
          if (pep.map->at(b).stop + 1 >= db[pep.map->at(b).index].sequence.size()) pr.nextAA = '-';
          else pr.nextAA = db[pep.map->at(b).index].sequence[pep.map->at(b).stop + 1];
          pr.startPos = pep.map->at(b).start+1;
          res.proteins.push_back(pr);
        }

        //Check for Reporter Ions
        for (size_t b = 0; b<params->rIons.size(); b++){
          if (spec[a]->checkReporterIon(params->rIons[b], params)){
            res.rIon.push_back(params->rIons[b]);
          }
        }

        vRes.push_back(res);
        scoreIndex++;
        if(scoreIndex==shorts.size()) break;
      }
    }
    //Uncomment this to report spectra with null results
    //if(vRes.size()==0){
    //  mResults res;
    //  vRes.push_back(res);
    //}

    exportTXT(fTXT, vRes);
    if(params->exportPepXML) exportPepXML(pepxml,vRes);
    if(params->exportPercolator) {
      if (fPerc2 == NULL) exportPercolator(fPerc, vRes);
      else {
        vector<mResults> v1;
        vector<mResults> v2;
        for(size_t z=0;z<vRes.size();z++){
          if(vRes[z].openMod.empty()) v1.push_back(vRes[z]);
          else v2.push_back(vRes[z]);
        }
        exportPercolator(fPerc,v1);
        exportPercolator(fPerc2,v2);
      }
    }

  }

  fclose(fTXT);
  if(pepxml!=NULL) {
    pepxml->write(fXML.c_str(),true);
    delete pepxml;
  }
  if(fPerc!=NULL) fclose(fPerc);
  if(fDiag!=NULL) fclose(fDiag);

  return true;

}

//Reads in raw/mzXML/mzML files. Other formats supported in MSToolkit as well.
bool MData::readSpectra(){

  MSReader   msr;
  Spectrum*   s;
  Spectrum   c;
  MSpectrum pls(*params);
  mPrecursor pre;

  int totalScans=0;
  int totalPeaks=0;
  int collapsedPeaks=0;
  int finalPeaks=0;
  int iPercent=0;
  int iTmp;

  deque<mMS2struct*> dMS2;
  vMS1Buffer.reserve(2000);
  memoryAllocate();
  initHardklor();

  Threading::CreateMutex(&mutexLockMS1);

  ThreadPool<mMS2struct*>* threadPool = new ThreadPool<mMS2struct*>(processMS2, params->threads, params->threads, 1);

  for (size_t a = 0; a<spec.size(); a++) delete spec[a];
  spec.clear();

  msr.setFilter(MS1);
  msr.addFilter(MS2);

  //Set progress meter
  printf("%2d%%", iPercent);
  fflush(stdout);

  s = new Spectrum;

  if (!msr.readFile(params->msFile.c_str(), *s)) return false;

  //temporary
  int nextMS2 = 0;

  while (s->getScanNumber()>0){

    totalScans++;
    if (s->size()<1) {
      msr.readFile(NULL, *s);
      continue;
    }

    if (s->getMsLevel() == 1) {
      vMS1Buffer.emplace_back(s);
      if (vMS1Buffer.size() == 10){  //When buffer is full, transfer to MS1 memory pool
        for (int a = 0; a<params->threads; a++) Threading::LockMutex(mutexHardklor[a]);
        while (spec.size()>0 && dMS1.size()>0 && dMS1.front()->getRTime()<spec.back()->getRTime() - 1){ //clear old memory
          delete dMS1.front();
          dMS1.pop_front();
        }
        for (size_t a = 0; a<vMS1Buffer.size(); a++){
          dMS1.emplace_back(vMS1Buffer[a]);
          vMS1Buffer[a] = NULL;
        }
        vMS1Buffer.clear();
        for (int a = 0; a<params->threads; a++) Threading::UnlockMutex(mutexHardklor[a]);
      }

    } else {
      mMS2struct* ms = new mMS2struct(s, params);
      dMS2.emplace_back(ms);
    }

    //Update progress meter
    iTmp = msr.getPercent();
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

    while (dMS2.size()>0 && dMS2[0]->state >= 3){ //copy and/or clear finished MS2 spectra
      if (dMS2[0]->state == 3) spec.push_back(dMS2[0]->pls);
      else delete dMS2[0]->pls;
      dMS2.pop_front();
      nextMS2--;
    }

    //Launch next MS2
    while (nextMS2<dMS2.size()) {
      if (dMS1.size()>0 && (dMS2[nextMS2]->s->getRTime() + 2)<dMS1.back()->getRTime()){
        //only launch this if there are enough MS1 for precursor analysis
        dMS2[nextMS2]->thread = true;
        threadPool->Launch(dMS2[nextMS2]);
        nextMS2++;
      } else break; //we got here because we need more MS1 first.
    }

    s = new Spectrum;
    msr.readFile(NULL, *s);

  }

  //finish flushing buffer
  for (int a = 0; a<params->threads; a++) Threading::LockMutex(mutexHardklor[a]);
  while (spec.size()>0 && dMS1.size()>0 && dMS1.front()->getRTime()<spec.back()->getRTime() - 1){ //clear old memory
    delete dMS1.front();
    dMS1.pop_front();
  }
  for (size_t a = 0; a<vMS1Buffer.size(); a++){
    dMS1.emplace_back(vMS1Buffer[a]);
    vMS1Buffer[a] = NULL;
  }
  vMS1Buffer.clear();
  for (int a = 0; a<params->threads; a++) Threading::UnlockMutex(mutexHardklor[a]);

  while (nextMS2<dMS2.size()) {
    dMS2[nextMS2]->thread = true;
    threadPool->Launch(dMS2[nextMS2]);
    nextMS2++;
  }

  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //finish processing last MS2 scans
  while (dMS2.size()>0){
    while (dMS2.size()>0 && dMS2[0]->state >= 3){ //copy and/or clear finished MS2 spectra
      if (dMS2[0]->state == 3) spec.push_back(dMS2[0]->pls);
      else delete dMS2[0]->pls;
      dMS2.pop_front();
      nextMS2--;
    }
  }

  //Finalize progress meter
  if (iPercent<100) printf("\b\b\b100%%");
  cout << endl;

  //clean up remaining memory
  while (dMS1.size()>0){
    delete dMS1.front();
    dMS1.pop_front();
  }

  delete threadPool;
  threadPool = NULL;
  memoryFree();
  releaseHardklor();
  Threading::DestroyMutex(mutexLockMS1);


  cout << "  " << spec.size() << " total spectra have enough data points for searching." << endl;
  //cout << totalScans << " total scans were loaded." <<  endl;
  //cout << totalPeaks << " total peaks in original data." << endl;
  //cout << collapsedPeaks << " peaks after collapsing." << endl;
  //cout << finalPeaks << " peaks after top N." << endl;
	//return true;

 // cout << "  " << specCounts << " spectra with " << peakCounts << " peaks will be analyzed." << endl;

  //Build mass list - this orders all precursor masses, with an index pointing to the actual
  //array position for the spectrum. This is because all spectra will have more than 1
  //precursor mass
  mMass m;
  massList.clear();
  for (int i = 0; i<spec.size(); i++){
    m.index = i;
    for (int j = 0; j<spec[i]->sizePrecursor(); j++){
      m.mass = spec[i]->getPrecursor(j).monoMass;
      massList.push_back(m);
    }
  }

  //sort mass list from low to high
  qsort(&massList[0], massList.size(), sizeof(mMass), compareMassList);

  if (bScans != NULL) delete[] bScans;
  bScans = new bool[spec.size()];

  return true;
}

void MData::releaseHardklor(){
  for (int a = 0; a<params->threads; a++){
    delete h[a];
    Threading::DestroyMutex(mutexHardklor[a]);
    delete averagine[a];
    delete mercury[a];
  }
  delete[] h;
  delete[] mutexHardklor;
  delete[] bHardklor;
  delete[] averagine;
  delete[] mercury;
  delete models;
}

void MData::setAdductSites(string s){
  for(size_t i=0;i<s.size();i++) adductSite[s[i]]=true;
}

void MData::setLog(MLog* c){
  mlog = c;
}

void MData::setParams(MParams* p){
  parObj=p;
}

void MData::setVersion(const char* v){
  strcpy(version,v);
}

int MData::size(){
  return (int)spec.size();
}

/*============================
  Private Utilities
============================*/
//First derivative method, returns base peak intensity of the set
void MData::centroid(Spectrum* s, MSpectrum* out, double resolution, int instrument){
  int i,j;
  float maxIntensity;
  int bestPeak;
  bool bLastPos;

	int nextBest;
	double FWHM;
  double maxMZ = (*s)[s->size()-1].mz+1.0;
	mSpecPoint centroid;

	vector<double> x;
	vector<double> y;
	vector<double> c;
	int left, right;
	bool bPoly;
	float lastIntensity;

	out->clear();

  bLastPos=false;
	for(i=0;i<s->size()-1;i++){

    if((*s)[i].intensity<(*s)[i+1].intensity) {
      bLastPos=true;
      continue;
    } else {
      if(bLastPos){
				bLastPos=false;

				//find max and add peak
				maxIntensity=0;
				for(j=i;j<i+1;j++){
				  if ((*s)[j].intensity>maxIntensity){
				    maxIntensity=(*s)[j].intensity;
				    bestPeak = j;
				  }
				}

				//walk left and right to find bounds above half max
				left=right=bestPeak;
				lastIntensity=maxIntensity;
				for(left=bestPeak-1;left>0;left--){
					if((*s)[left].intensity<(maxIntensity/3) || (*s)[left].intensity>lastIntensity){
						left++;
						break;
					}
					lastIntensity=(*s)[left].intensity;
				}
				lastIntensity=maxIntensity;
				for(right=bestPeak+1;right<s->size()-1;right++){
					if((*s)[right].intensity<(maxIntensity/3) || (*s)[right].intensity>lastIntensity){
						right--;
						break;
					}
					lastIntensity=(*s)[right].intensity;
				}

				//if we have at least 5 data points, try polynomial fit
				double r2;
				bPoly=false;
				if((right-left+1)>4){
					x.clear();
					y.clear();
					for(j=left;j<=right;j++){
						x.push_back((*s)[j].mz);
						y.push_back(log((*s)[j].intensity));
					}
					r2=polynomialBestFit(x,y,c);
					if(r2>0.95){
						bPoly=true;
            centroid.mass = -c[1] / (2 * c[2]) + c[3];
						centroid.intensity=(float)exp(c[0]-c[2]*(c[1]/(2*c[2]))*(c[1]/(2*c[2])));
					} else {

					}
				}

				if(!bPoly){
					//Best estimate of Gaussian centroid
					//Get 2nd highest point of peak
					if(bestPeak==s->size()) nextBest=bestPeak-1;
					else if((*s)[bestPeak-1].intensity > (*s)[bestPeak+1].intensity) nextBest=bestPeak-1;
					else nextBest=bestPeak+1;

					//Get FWHM
					switch(instrument){
						case 0: FWHM = (*s)[bestPeak].mz*sqrt((*s)[bestPeak].mz)/(20*resolution); break;  //Orbitrap
						case 1: FWHM = (*s)[bestPeak].mz*(*s)[bestPeak].mz/(400*resolution); break;				//FTICR
						default: break;
					}

					//Calc centroid MZ (in three lines for easy reading)
					centroid.mass = pow(FWHM,2)*log((*s)[bestPeak].intensity/(*s)[nextBest].intensity);
          centroid.mass /= GAUSSCONST*((*s)[bestPeak].mz - (*s)[nextBest].mz);
          centroid.mass += ((*s)[bestPeak].mz + (*s)[nextBest].mz) / 2;

					//Calc centroid intensity
          centroid.intensity = (float)((*s)[bestPeak].intensity / exp(-pow(((*s)[bestPeak].mz - centroid.mass) / FWHM, 2)*GAUSSCONST));
				}

				//some peaks are funny shaped and have bad gaussian fit.
				//if error is more than 10%, keep existing intensity
				if( fabs(((*s)[bestPeak].intensity - centroid.intensity) / centroid.intensity * 100) > 10 ||
            //not a good check for infinity
            centroid.intensity>9999999999999.9 ||
            centroid.intensity < 0 ) {
					centroid.intensity=(*s)[bestPeak].intensity;
				}

				//Hack until I put in mass ranges
        if (centroid.mass<0 || centroid.mass>maxMZ) {
					//do nothing if invalid mz
				} else {
					out->addPoint(centroid);
				}
			
      }

    }
  }

}

//Function tries to remove isotopes of signals by stacking the intensities on the monoisotopic peak
//Also creates an equal n+1 peak in case wrong monoisotopic peak was identified.
void MData::collapseSpectrum(MSpectrum& s){
  int i,j,k,n;
  int charge,z;
  int maxIndex;
  float max;
  float cutoff;
  vector<int> dist;

  vector<mSpecPoint> s2;

  while(true){
    max=0.1f;
    for(i=0;i<s.size();i++){
      if(s[i].intensity>max){
        max=s[i].intensity;
        maxIndex=i;
      }
    }

    //finish and exit function
    if(max<1) break;

    dist.clear();
    dist.push_back(maxIndex);

    //check right
    j=maxIndex+1;
    while (j<s.size() && (s[j].mass - s[maxIndex].mass)<1.1){
      if(s[j].intensity<1) {
        j++;
        continue;
      }
      charge=getCharge(s,maxIndex,j);

      if(charge==0){
        j++;
        continue;
      }

      //try stepping along at same charge state here out
      //note that if this doesn't work, it doesn't go back and look for a different charge state
      dist.push_back(j);
      k=j;
      n=j+1;
      while (n<s.size() && (s[n].mass - s[k].mass)<1.1){
        if(s[n].intensity<1) {
          n++;
          continue;
        }
        z=getCharge(s,k,n);
        if(z>0 && z<charge) {
          break;
        } else if (z == charge && (s[n].mass - s[k].mass)>(0.98 / charge) && (s[n].mass - s[k].mass)<(1.02 / charge)) {
          dist.push_back(n);
          k=n;
          n++;
        } else {
          n++;
        }
      }
      break;
    }

    //if nothing found to the right, quit here?
    if(dist.size()==1){
      s2.push_back(s[dist[0]]);
      s[dist[0]].intensity=0;
      continue;
    }

    //step to the left
    j=maxIndex-1;
    while (j >= 0 && (s[maxIndex].mass - s[j].mass)<1.1){
      if(s[j].intensity<1) {
        j--;
        continue;
      }
      z=getCharge(s,j,maxIndex);
      if(z!=charge){
        j--;
        continue;
      }

      //try stepping along at same charge state here out
      dist.push_back(j);
      k=j;
      n=j-1;
      while (n >= 0 && (s[k].mass - s[n].mass)<1.1){
        if(s[n].intensity<1) {
          n--;
          continue;
        }
        z=getCharge(s,n,k);
        //printf("\tleft\t%.6lf\t%.6lf\t%d\n",s[n].mz,s[k].mz-s[n].mz,z);
        if(z>0 && z<charge) {
          break;
        } else if (z == charge && s[k].mass - s[n].mass > 0.98 / charge && s[k].mass - s[n].mass < 1.02 / charge) {
          dist.push_back(n);
          k=n;
          n--;
        } else {
          n--;
        }
      }
      break;
    }


    //Only accept size of 2 if charge is 1 or 2
    if(dist.size()==2){
      if(charge<3){
        max = s[dist[0]].intensity + s[dist[1]].intensity;
        mSpecPoint sp;
        sp.mass = s[dist[0]].mass;
        sp.intensity = max;
        s2.push_back(sp);
        s[dist[1]].intensity = 0;
       // s2.add(s[dist[1]].mz,max);
      } else {
        s2.push_back(s[dist[0]]);
       // s2.add(s[dist[1]]);
      }
      s[dist[0]].intensity=0;
      //s[dist[1]].intensity=0;
    } else {
      cutoff=max/20;
      max=0;
      j=dist[0];
      k=dist[1];
      for(i=0;i<(int)dist.size();i++) {
        if(dist[i]<j && s[dist[i]].intensity>cutoff){
          k=j;
          j=dist[i];
        }
        if(s[dist[i]].intensity>cutoff){
          max+=s[dist[i]].intensity;
          s[dist[i]].intensity=0;
        }
      }
      mSpecPoint sp;
      sp.mass = s[j].mass;
      sp.intensity = max;
      s2.push_back(sp);
      //s2.add(s[k].mz,max);
    }

  }

  sort(s2.begin(), s2.end(), compareSpecPoint);
  s.clear();
  for(i=0;i<s2.size();i++) {
    if(i<s2.size()-1 && s2[i].mass==s2[i+1].mass){
      if (s2[i].intensity>s2[i + 1].intensity) s.addPoint(s2[i]);
      else s.addPoint(s2[i + 1]);
      i++;
    } else {
      s.addPoint(s2[i]);
    }
  }
  
}

int MData::compareInt(const void *p1, const void *p2){
  int d1 = *(int *)p1;
  int d2 = *(int *)p2;
  if(d1<d2) {
		return -1;
	} else if(d1>d2) {
  	return 1;
  } else {
	  return 0;
  }
}

int MData::compareMassList(const void *p1, const void *p2){
  mMass d1 = *(mMass *)p1;
  mMass d2 = *(mMass *)p2;
  if(d1.mass<d2.mass) {
		return -1;
	} else if(d1.mass>d2.mass) {
  	return 1;
  } else {
	  return 0;
  }
}

int MData::getCharge(MSpectrum& s, int index, int next){
  double mass;

  mass=s[next].mass-s[index].mass;
  if(mass>0.98 && mass<1.02) return 1;
  else if(mass>0.49 && mass<0.51) return 2;
  else if(mass>0.3267 && mass<0.34) return 3;
  else if(mass>0.245 && mass<0.255) return 4;
  else if(mass>0.196 && mass<0.204) return 5;
  else if(mass>0.1633 && mass<0.17) return 6;
  else return 0;
  /*
  if(mass>0.99 && mass<1.007) return 1;
  else if(mass>0.495 && mass<0.5035) return 2;
  else if(mass>0.33 && mass<0.335667) return 3;
  else if(mass>0.2475 && mass<0.25175) return 4;
  else if(mass>0.198 && mass<0.2014) return 5;
  else if(mass>0.165 && mass<0.1678333) return 6;
  else return 0;
  */

}

double MData::polynomialBestFit(vector<double>& x, vector<double>& y, vector<double>& coeff, int degree){
	if(degree>3){
		cout << "High order polynomials not supported with this function. Max=3" << endl;
		exit(1);
	}

	if(degree<2){
		cout << "Polynomials need at least two degrees. Min=2" << endl;
		exit(1);
	}

	int i,j,a;
	int n=(int)x.size();
	degree++;

	double sFactor=x[n/2];

	//set X matrix
	double** X = new double* [n];
	for(i=0;i<n;i++){
		X[i] = new double [degree];
		X[i][0] = 1.0;
		for(j=1;j<degree;j++) X[i][j]=X[i][j-1]*(x[i]-sFactor);
	}

	//make transpose of X
	double** Xt = new double* [degree];
	for(j=0;j<degree;j++) Xt[j] = new double [n];
	for(i=0;i<n;i++){
		for(j=0;j<degree;j++){
			Xt[j][i] = X[i][j];
		}
	}

	//matrix multiplication
	double** XtX = new double* [degree];
	for(i=0;i<degree;i++){
		XtX[i] = new double [degree];
		for(j=0;j<degree;j++){
			XtX[i][j]=0;
			for(a=0;a<n;a++) XtX[i][j]+=(Xt[i][a]*X[a][j]);
		}
	}

	//inverse using Gauss-Jordan Elimination
	double** XtXi = new double* [degree];
	for(i=0;i<degree;i++){
		XtXi[i] = new double [degree*2];
		for(j=0;j<degree*2;j++){
			if(j<degree) XtXi[i][j]=XtX[i][j];
			else if(j-degree==i) XtXi[i][j]=1;
			else XtXi[i][j]=0;
		}
	}
	double d;
	for(j=0;j<degree;j++){
		for(i=0;i<degree;i++){
			if(i==j) continue;
			if(XtXi[i][j]==0) continue;
			d=-XtXi[i][j]/XtXi[j][j];
			for(a=0;a<degree*2;a++) XtXi[i][a]+=(d*XtXi[j][a]);
		}
	}
	for(i=0;i<degree;i++){
		d=1/XtXi[i][i];
		for(j=0;j<degree*2;j++) XtXi[i][j]*=d;
	}

	//matrix multiplication
	double* Xty = new double [degree];
	for(i=0;i<degree;i++){
		Xty[i]=0;
		for(j=0;j<n;j++) Xty[i]+=Xt[i][j]*y[j];	
	}

	//matrix multiplication
	double* c = new double [degree];
	for(i=0;i<degree;i++){
		c[i]=0;
		for(j=0;j<degree;j++) c[i]+=XtXi[i][j+degree]*Xty[j];	
	}

	coeff.clear();
	for(i=0;i<degree;i++) {
		coeff.push_back(c[i]);
	}
	coeff.push_back(sFactor);

	vector<double> z;
	for(i=0;i<n;i++) z.push_back((x[i]-sFactor)*(x[i]-sFactor)*c[2]+(x[i]-sFactor)*c[1]+c[0]);

	//clean up memory
	delete [] c;
	delete [] Xty;
	for(i=0;i<degree;i++){
		delete [] XtXi[i];
		delete [] XtX[i];
		delete [] Xt[i];
	}
	for(i=0;i<n;i++) delete [] X[i];
	delete [] X;
	delete [] Xt;
	delete [] XtX;
	delete [] XtXi;

	double sxy=0;
  double sxx=0;
  double syy=0;
	double xavg=0;
	double yavg=0;
  for(i=0;i<n;i++){
		xavg+=z[i];
		yavg+=y[i];
	}
	xavg/=n;
	yavg/=n;
  for(i=0;i<n;i++){
    sxy += ((z[i]-xavg)*(y[i]-yavg));
    sxx += ((z[i]-xavg)*(z[i]-xavg));
    syy += ((y[i]-yavg)*(y[i]-yavg));
  }
	double r2 = (sxy*sxy)/(sxx*syy);

	return r2;

}

void MData::processMS2(mMS2struct* s){

  int j;

  formatMS2(s->s, s->pls);

  if (s->pls->size() <= params->minPeaks){
    s->state = 4;
    s->thread = false;
    return;
  }

  bool bAddHardklor = false;
  bool bAddEstimate = false;

  if (params->preferPrecursor == 1){
    if (s->pls->sizePrecursor() == 0){
      if (s->pls->getCharge()>0) bAddEstimate = true;
      else bAddHardklor = true;
    }
  } else if (params->preferPrecursor == 0){
    s->pls->clearPrecursors();
    bAddHardklor = true;
    if (s->pls->getCharge()) bAddEstimate = true;
  } else {
    bAddHardklor = true;
    if (s->pls->getCharge()) bAddEstimate = true;
  }

  int ret;
  if (bAddHardklor && params->precursorRefinement){
    //only do Hardklor analysis if data contain precursor scans
    //ret = pre.getSpecRange(*spec[i]);
    int tIndex;
    Threading::LockMutex(mutexLockMS1);
    for (tIndex = 0; tIndex<params->threads; tIndex++){
      if (!bHardklor[tIndex]){
        bHardklor[tIndex] = true;
        break;
      }
    }
    if (tIndex == params->threads) cout << "Thread overload" << endl;
    Threading::UnlockMutex(mutexLockMS1);

    Threading::LockMutex(mutexHardklor[tIndex]);
    ret = processPrecursor(s, tIndex);
    bHardklor[tIndex] = false;
    Threading::UnlockMutex(mutexHardklor[tIndex]);
    //Threading::LockMutex(mutexLockMS1); //is this necessary?
    //bHardklor[tIndex] = false;
    //Threading::UnlockMutex(mutexLockMS1);
  }

  if (bAddEstimate){
    mPrecursor pr;
    pr.monoMass = s->pls->getMZ()*s->pls->getCharge() - 1.007276466*s->pls->getCharge();
    pr.charge = s->pls->getCharge();
    pr.corr = -5;
    s->pls->setCharge(pr.charge);
    s->pls->addPrecursor(pr, params->topCount);
    for (int px = 1; px <= params->isotopeError; px++){
      if (px == 4) break;
      pr.monoMass -= 1.00335483;
      pr.corr -= 0.1;
      s->pls->addPrecursor(pr, params->topCount);
    }
  }

  //Now clean up any duplicate precursors. They should already be in order of priority. Use 5ppm as tolerance
  for (int k = 0; k<s->pls->sizePrecursor(); k++){
    for (int n = k + 1; n<s->pls->sizePrecursor(); n++){
      double m1 = s->pls->getPrecursor(k).monoMass;
      double m2 = s->pls->getPrecursor(n).monoMass;
      double m = (m1 - m2) / m1*1e6;
      if (fabs(m)<5){
        s->pls->erasePrecursor(n);
        n--;
      }
    }
  }

  if (s->pls->sizePrecursor()>0){
    //build singletList
    s->pls->resetSingletList();
    s->pls->peakCounts = s->pls->size();

  } else {
    s->state = 4; //no precursors, so advance state past transform to delete.
    s->thread = false;
    return;
  }

  Threading::LockMutex(mutexMemoryPool);
  for (j = 0; j<params->threads; j++){
    if (!memoryPool[j]){
      memoryPool[j] = true;
      break;
    }
  }
  Threading::UnlockMutex(mutexMemoryPool);

  if (j == params->threads){
    cout << "Error in KData::processMS2::state==2" << endl;
    exit(-1);
  }
  s->pls->kojakXCorr(tempRawData[j], tmpFastXcorrData[j], fastXcorrData[j], preProcess[j]);
  memoryPool[j] = false;

  s->state = 3;
  s->thread = false;
}

//Takes relative path and finds absolute path
bool MData::processPath(const char* in_path, char* out_path){
  char cwd[1024];
  char* ret=getcwd(cwd,1024);

  //if windows or unix in_path, just copy it to out_path
  if (strlen(in_path) > 0 && in_path[0] == '/'){ //unix
    strcpy(out_path,in_path);
    return true;
  }
  if (strlen(in_path) > 1 && in_path[1] == ':'){ //windows
    strcpy(out_path, in_path);
    return true;
  }

  //tokenize cwd
  char* tok;
  char str[1024];
  strcpy(str,cwd);
  string s;
  vector<string> v;

  tok = strtok(str, "\\/\n\r");
  while (tok != NULL){
    s=tok;
    v.push_back(s);
    tok = strtok(NULL, "\\/\n\r");
  }

  //tokenize in_path
  strcpy(str,in_path);

  tok = strtok(str, "\\/\n\r");
  while (tok != NULL){
    if (strcmp(tok, "..") == 0) {
      v.pop_back();
    } else if (strcmp(tok, ".") == 0){
      //do nothing
    } else {
      s=tok;
      v.push_back(s);
    }
    tok = strtok(NULL, "\\/\n\r");
  }

  //build absolute path
#ifdef _MSC_VER
  s.clear();
#else
  s.clear();
  s+=slashdir;
#endif
  for (size_t i = 0; i < v.size(); i++){
    s += v[i];
    s += slashdir;
  }
  s[s.size() - 1] = '\0';
  strcpy(out_path, &s[0]);
  return true;

}

string MData::processPeptide(mPeptide& pep, vector<mPepMod>& mod, int site, double massA, MDatabase& db){
  char tmp[32];
  size_t j, k;
  string seq = "";
  string peptide;

  db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, peptide);

  if (pep.nTerm && aa.getFixedModMass('$') != 0) {
    sprintf(tmp, "n[%.0lf]", aa.getFixedModMass('$'));
    seq += tmp;
  }
  for (k = 0; k<mod.size(); k++){ //check for n-terminal peptide mod
    if (mod[k].pos == 0 && mod[k].term){
      sprintf(tmp, "n[%.0lf]", mod[k].mass);
      seq += tmp;
    }
  }
  //uncomment to display open, unlocalized mod on the n-terminus
  //if(massA!=0 && site==-99){
  //  sprintf(tmp, "o[%.0lf]", massA);
  //  seq += tmp;
  //}
  for (j = 0; j<peptide.size(); j++) {
    seq += peptide[j];
    for (k = 0; k<mod.size(); k++){
      if (j == (unsigned int)mod[k].pos && !mod[k].term){
        sprintf(tmp, "[%.0lf]", mod[k].mass);
        seq += tmp;
      }
    }
    //uncomment to display the open, localized mod in the peptide sequence
    //if (j == (size_t)site){
    //  sprintf(tmp, "[%.0lf]", massA);
    //  seq += tmp;
    //}
  }
  for (k = 0; k<mod.size(); k++){ //check for c-terminal peptide mod
    if (mod[k].pos > 0 && mod[k].term){
      sprintf(tmp, "c[%.0lf]", mod[k].mass);
      seq += tmp;
    }
  }
  if (pep.cTerm && aa.getFixedModMass('%') != 0) {
    sprintf(tmp, "c[%.0lf]", aa.getFixedModMass('%'));
    seq += tmp;
  }

  return seq;
}

int MData::processPrecursor(mMS2struct* s, int tIndex){

  int j;
  float rt = s->pls->getRTime();
  double mz = s->pls->getMZ();
  float maxIntensity = 0;
  float maxRT = 0;
  int best = 0;
  int precursor = -1;
  int ret = 0;

  //Find the MS1 scan that contains the precursor ion within 6 seconds of when it was acquired.
  for (int i = 0; i<dMS1.size(); i++){
    if (dMS1[i]->getRTime()<rt - 0.167) continue;
    if (dMS1[i]->getRTime()>rt + 0.167) break;
    int j = findPeak(dMS1[i], mz, 10);

    if (j>-1){
      if (dMS1[i]->at(j).intensity>maxIntensity){
        maxIntensity = dMS1[i]->at(j).intensity;
        maxRT = dMS1[i]->getRTime();
        best = dMS1[i]->getScanNumber();
        precursor = i;
      }
    }
  }

  if (precursor<0){
    //cout << "Warning: Precursor not found for " << scanNum << " " << mz << endl;
    return ret;
  }

  //Get up to +/-15 sec of spectra around the max precursor intensity
  //This is done by extending on both sides until a gap is found or time is reached.
  //Additionally, stop if 2 scans found flanking either side (maximum 5 scans per precursor).
  vector<Spectrum*> vs;
  float rtHigh;
  float rtLow;
  rtHigh = rtLow = dMS1[precursor]->getRTime();
  int k = 0;
  for (int i = precursor; i<dMS1.size(); i++){
    if (dMS1[i]->getRTime()>maxRT + 0.25) break;
    j = findPeak(dMS1[i], mz, 10);
    if (j<0) break;
    vs.push_back(dMS1[i]);
    rtHigh = dMS1[i]->getRTime();
    k++;
    if (k == 3) break;
  }

  k = 0;
  int i = precursor;
  while (i>0){
    i--;
    if (dMS1[i]->getRTime()<maxRT - 0.25) break;
    j = findPeak(dMS1[i], mz, 10);
    if (j<0) break;
    vs.push_back(dMS1[i]);
    rtLow = dMS1[i]->getRTime();
    k++;
    if (k == 2) break;
  }

  //Average points between mz-1.5 and mz+2
  Spectrum sp;
  averageScansCentroid(vs, sp, mz - 1.0, mz + 1.5);
  if (sp.size() == 0) {
    cout << "\n   WARNING: Unexpected precursor scan data!";
    if (!params->ms1Centroid) cout << " Params are set to MS1 profile mode, but are MS1 scans centroided?" << endl;
    return ret;
  }
  sp.setScanNumber(dMS1[precursor]->getScanNumber());

  //Obtain the possible precursor charge states of the selected ion.
  //Find the index of the closest peak to the selected m/z.
  vector<int> preCharges;
  double tmz = fabs(mz - sp[0].mz);
  for (j = 1; j<sp.size(); j++){
    if (fabs(mz - sp[j].mz)<tmz) tmz = fabs(mz - sp[j].mz);
    else break;
  }
  j = j - 1;
  h[tIndex]->QuickCharge(sp, j, preCharges);

  //Clear corr
  double corr = 0;
  double monoMass = 0;
  int charge = 0;
  mPrecursor pre;
  h[tIndex]->GoHardklor(hs, &sp);

  //If nothing was found, really narrow down the window and try again.
  if (h[tIndex]->Size() == 0){
    averageScansCentroid(vs, sp, mz - 0.6, mz + 1.2);
    sp.setScanNumber(dMS1[precursor]->getScanNumber());
    h[tIndex]->GoHardklor(hs, &sp);
  }

  float intensity = 0;
  for (j = 0; j<h[tIndex]->Size(); j++){

    //Must have highest intensity and intersect isolated peak.
    if (h[tIndex]->operator[](j).intensity<intensity) continue;
    tmz = (h[tIndex]->operator[](j).monoMass + 1.007276466*h[tIndex]->operator[](j).charge) / h[tIndex]->operator[](j).charge;
    while (tmz<(s->pls->getMZ() + 0.01)){
      if (fabs(tmz - s->pls->getMZ())<0.01){
        monoMass = h[tIndex]->operator[](j).monoMass;
        charge = h[tIndex]->operator[](j).charge;
        corr = h[tIndex]->operator[](j).corr;
        intensity = h[tIndex]->operator[](j).intensity;
        ret = 1;
        break;
      }
      tmz += (1.00335483 / h[tIndex]->operator[](j).charge);
    }
  }

  //failing to match precursor peak, keep most intense precursor in presumed isolation window
  if (corr == 0){
    for (j = 0; j<h[tIndex]->Size(); j++){
      if (h[tIndex]->operator[](j).intensity>intensity){
        monoMass = h[tIndex]->operator[](j).monoMass;
        charge = h[tIndex]->operator[](j).charge;
        corr = h[tIndex]->operator[](j).corr;
        intensity = h[tIndex]->operator[](j).intensity;
        ret = 2;
      }
    }
  }

  if (corr>0){
    pre.monoMass = monoMass;
    pre.charge = charge;
    pre.corr = corr;
    pre.label = 0;
    s->pls->addPrecursor(pre, params->topCount);
    //also add isotope error
    if (params->isotopeError>0){
      pre.monoMass -= 1.00335483;
      pre.corr = -1;
      s->pls->addPrecursor(pre, params->topCount);
    }
    if (params->isotopeError>1){
      pre.monoMass -= 1.00335483;
      pre.corr = -2;
      s->pls->addPrecursor(pre, params->topCount);
    }
    if (params->isotopeError>2){
      pre.monoMass -= 1.00335483;
      pre.corr = -3;
      s->pls->addPrecursor(pre, params->topCount);
    }
  }

  //Assume two precursors with nearly identical mass (within precursor tolerance) are the same.
  //This can occur when checking multiple enrichment states.
  //Keep only the higher correlated precursor.
  if (s->pls->sizePrecursor()>1){
    bool bCheck = true;
    while (bCheck){
      bCheck = false;
      for (k = 0; k<s->pls->sizePrecursor() - 1; k++){
        for (j = k + 1; j<s->pls->sizePrecursor(); j++){
          if (fabs(s->pls->getPrecursor(k).monoMass - s->pls->getPrecursor(j).monoMass) / s->pls->getPrecursor(k).monoMass*1e6 < params->ppmPrecursor){
            if (s->pls->getPrecursor(k).corr>s->pls->getPrecursor(j).corr) s->pls->erasePrecursor(j);
            else s->pls->erasePrecursor(k);
            bCheck = true;
            break;
          }
        }
        if (bCheck) break;
      }
    }
  }

  return ret;
}

void MData::processPSM(MSpectrum& s, mScoreCard3& sc, mResults& r){
  mPrecursor p=s.getPrecursor(sc.precursor);
  r.monoMass=p.monoMass;
  r.charge=p.charge;

  r.psmMass=sc.mass+sc.massA;
  r.ppm=(r.monoMass-r.psmMass)/r.psmMass*1e6;
  r.eValue=sc.eVal;
  r.scoreMagnum=sc.simpleScore;

  //temporary structure for holding mod info
  typedef struct m2{
    double mass;
    int count;
    vector<char> pos;
  } m2;

  //iterate over all peptide variants
  vector<m2> mods;
  for(size_t a=0; a<sc.mSet.size(); a++){
    //iterate over all mods in this variant
    for(size_t b=0;b<sc.mSet[a].mods.size();b++){
      
      //convert mod pos
      char pos;
      if(sc.mSet[a].mods[b].term){
        if(sc.mSet[a].mods[b].pos>0) pos=127;
        else pos=126;
      } else pos = sc.mSet[a].mods[b].pos;

      //iterate over processed mods to find alternate positions
      size_t c;
      for(c=0;c<mods.size();c++){
        if(mods[c].mass==sc.mSet[a].mods[b].mass){
          if(a==0) mods[c].count++;
          size_t d;
          for(d=0;d<mods[c].pos.size();d++){            
            if(mods[c].pos[d]==sc.mSet[a].mods[b].pos) break;
          }
          if(d==mods[c].pos.size()) mods[c].pos.push_back(sc.mSet[a].mods[b].pos);
          break;
        }
      }
      if(c==mods.size()){ //add new mod
        m2 m;
        m.mass=sc.mSet[a].mods[b].mass;
        m.count=1;
        m.pos.push_back(sc.mSet[a].mods[b].pos);
        mods.push_back(m);
      }
    }
  }

  //process mod string
  for(size_t a=0;a<mods.size();a++){
    if (a>0) r.mods += ";";
    char str[512];
    sprintf(str, "%.6lf[%d;", mods[a].mass,mods[a].count);
    string st = str;
    for (size_t b = 0; b<mods[a].pos.size(); b++){
      if (b>0) st += ",";
      if (mods[a].pos[b] == 126) st += "n";
      else if (mods[a].pos[b] == 127) st += "c";
      else {
        sprintf(str, "%d", (int)mods[a].pos[b]+1);
        st += str;
      }
    }
    st += "]";
    r.mods += st;
  }

  //if open modification
  if(sc.massA!=0){
    vector<char> vs;
    for (size_t a = 0; a<sc.aSites.size(); a++){
      size_t b;
      for(b=0;b<vs.size();b++){
        if(vs[b]==sc.aSites[a]) break;
      }
      if(b==vs.size()) vs.push_back(sc.aSites[a]);
    }

    char str[32];
    sprintf(str,"%.6lf[%d;",sc.massA,(int)vs.size());
    r.openMod=str;
    for(size_t a=0;a<vs.size();a++){
      if(a>0) r.openMod+=",";
      if(sc.aSites[a]==-99) r.openMod+="0";
      else {
        sprintf(str,"%d",(int)vs[a]+1);
        r.openMod+=str;
      }
    }
    r.openMod+="]";
  }

  //convert all mods to pepxml mods
  for (size_t a = 0; a<sc.mSet.size(); a++){
    rMods2 modset;
    for(size_t b=0;b<sc.mSet[a].mods.size();b++){
      rMods rm;
      rm.mass=sc.mSet[a].mods[b].mass;
      rm.variable=true;
      rm.adduct=false;
      rm.term=false;
      if(sc.mSet[a].mods[b].term){
        rm.term=true;
        if(sc.mSet[a].mods[b].pos==0) {
          rm.pos='n';
          if (aa.getFixedModMass('n') == rm.mass || aa.getFixedModMass('$') == rm.mass) rm.variable = false;
        } else {
          rm.pos='c';
          if (aa.getFixedModMass('c') == rm.mass || aa.getFixedModMass('%') == rm.mass) rm.variable = false;
        }
      } else {
        rm.pos=sc.mSet[a].mods[b].pos;
        if(aa.getFixedModMass(r.peptide[rm.pos])==rm.mass) rm.variable=false;
      }
      modset.mods.push_back(rm);
    }
    //add adduct (if any)
    if (sc.massA != 0){
      rMods rm;
      rm.mass = sc.massA;
      rm.variable = true;
      rm.adduct = true;
      rm.pos = sc.aSites[a];
      if(rm.pos<0) rm.term=true;
      else rm.term=false;
      modset.mods.push_back(rm);
    }
    r.vMods.push_back(modset);
  }

}

void MData::processSpectrumInfo(MSpectrum& s, mResults& r){
  r.scanNumber=s.getScanNumber();
  r.rTimeSec=s.getRTime()*60;
  r.selectedMZ=s.getMZ();
}

int MData::compareScanBinRev2(const void *p1, const void *p2){
  mScanBin d1 = *(mScanBin *)p1;
  mScanBin d2 = *(mScanBin *)p2;
  if (d1.intensity>d2.intensity) {
    return -1;
  } else if (d1.intensity<d2.intensity) {
    return 1;
  } else {
    return 0;
  }
}