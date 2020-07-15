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

/*============================
  Constructors
============================*/
MData::MData(){
  bScans=NULL;
  params=NULL;
  for(int i=0;i<128;i++) adductSite[i]=false;
}

MData::MData(mParams* p){
  bScans=NULL;
  params=p;
  size_t i;
  for(i=0;i<p->fMods->size();i++) aa.addFixedMod((char)p->fMods->at(i).index,p->fMods->at(i).mass);
  for (i = 0; i<128; i++) adductSite[i] = false;
}

MData::~MData(){
  params=NULL;
  if(bScans!=NULL) delete[] bScans;
}


/*============================
  Operators
============================*/
MSpectrum& MData::operator [](const int& i){
  return spec[i];
}


/*============================
  Functions
============================*/
bool* MData::getAdductSites(){
  return &adductSite[0];
}

MSpectrum* MData::getSpectrum(const int& i){
  return &spec[i];
}

MSpectrum& MData::at(const int& i){
  return spec[i];
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
    for (int q = 0; q<spec[b].sizePrecursor(); q++){
      if (spec[b].getTopPeps(q)->peptideCount == 0) continue;
      if (spec[b].getTopPeps(q)->peptideFirst->simpleScore>spec[b].getScoreCard(0).simpleScore){
        if (spec[b].getTopPeps(q)->peptideFirst->simpleScore>maxScore) {
          maxScore = spec[b].getTopPeps(q)->peptideFirst->simpleScore;
          unknownMass = spec[b].getPrecursor(q).monoMass - spec[b].getTopPeps(q)->peptideFirst->mass;
        }
      }
    }
    if (maxScore>0) oddCount++;
    if (maxScore>3) bigCount++;
    if (maxScore>3 && spec[b].getScoreCard(0).simpleScore>0 && maxScore>spec[b].getScoreCard(0).simpleScore * 2) twoCount++;
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

//This function tries to assign best possible 18O2 and 18O4 precursor ion mass values
//for all MS2 spectra
bool MData::mapPrecursors(){
  
  int iPercent=0;
  int iTmp;
  
  unsigned int i;
  int j,k;

  MPrecursor pre(params);
  mMass      m;

  int peakCounts=0;
  int specCounts=0;
  int ret;

  int prePre=0;
  int foundPre=0;
  int noPre=0;

  //Open the data file in the precursor mapping object
  //if(!pre.setFile(&p)) return false;

  //Print progress
  printf("  Mapping precursors ... %2d%%",iPercent);
  fflush(stdout);

  //Iterate all MS/MS spectra
  for(i=0;i<spec.size();i++){

    //Update progress
    iTmp=(int)(i*100.0/spec.size());
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }

    //If instrument determined precursors are preferred, only compute precursors if none supplied
    if(params->preferPrecursor==1 && spec[i].sizePrecursor()>0){
      prePre++;
      specCounts++;
      peakCounts+=spec[i].size();
      continue;
    }

    //if (spec[i].getScanNumber() == 25028) {
    //  cout << "Current precursors: " << spec[i].sizePrecursor() << endl;
    //  for (j = 0; j<spec[i].sizePrecursor(); j++) cout << spec[i].getPrecursor(j).monoMass << " " << spec[i].getPrecursor(j).charge << endl;
    //}

    //Find precursor using object function. Take results and copy them to spectra
    if (params->precursorRefinement){
      ret=pre.getSpecRange(spec[i]);
    } else {
      ret=0;
    }

    //if (spec[i].getScanNumber() == 1008) {
    //  cout << "Current precursors: " << spec[i].sizePrecursor() << "\tret = " << ret << endl;
    //  for (j = 0; j<spec[i].sizePrecursor(); j++) cout << spec[i].getPrecursor(j).monoMass << " " << spec[i].getPrecursor(j).charge << endl;
    //}

    if(ret>0){

      //if supplementing instrument predicted precursor, then chance for precursor
      //to be seen twice (first by instrument, then by Hardklor). Keep Hardklor result.
      if(spec[i].getInstrumentPrecursor() && spec[i].sizePrecursor()>1){
        if(spec[i].getPrecursor(0).charge==spec[i].getPrecursor(1).charge && 
           fabs((spec[i].getPrecursor(0).monoMass-spec[i].getPrecursor(1).monoMass)/spec[i].getPrecursor(0).monoMass*1e6)<10.0 ){
          spec[i].erasePrecursor(0);
          spec[i].setInstrumentPrecursor(false);
        }
      }

      //if precursor prediction doesn't overlap selected ion, predict additional
      //precursors using presumed charge states and the selected ion.
      if(ret==2) {
        pre.estimatePrecursor(spec[i]);

        //if (spec[i].getScanNumber() == 25028) {
        //  cout << "Current precursors: " << spec[i].sizePrecursor() << endl;
        //  for (j = 0; j<spec[i].sizePrecursor(); j++) cout << spec[i].getPrecursor(j).monoMass << " " << spec[i].getPrecursor(j).charge << endl;
        //}

        //if supplementing instrument predicted precursor, then chance for precursor
        //to be seen twice (first by instrument, then by charge prediction). Keep instrument.
        if(spec[i].getInstrumentPrecursor() && spec[i].sizePrecursor()>1){
          for(k=1;k<spec[i].sizePrecursor();k++){
            if(spec[i].getPrecursor(0).charge==spec[i].getPrecursor(k).charge && 
               fabs((spec[i].getPrecursor(0).monoMass-spec[i].getPrecursor(k).monoMass)/spec[i].getPrecursor(0).monoMass*1e6)<10.0 ){
              spec[i].erasePrecursor(k);
              k--;
            }
          }
        }

      }

    } else { 
      
      //If no precursors found, estimate using mercury and selected mass
      pre.estimatePrecursor(spec[i]);
      
      //if supplementing instrument predicted precursor, then chance for precursor
      //to be seen twice (first by instrument, then by charge prediction). Keep instrument.
      if(spec[i].getInstrumentPrecursor() && spec[i].sizePrecursor()>1){
        for(k=1;k<spec[i].sizePrecursor();k++){
          if(spec[i].getPrecursor(0).charge==spec[i].getPrecursor(k).charge && 
             fabs(spec[i].getPrecursor(0).monoMass-spec[i].getPrecursor(k).monoMass)<0.01 ){
            spec[i].erasePrecursor(k);
            k--;
          }
        }
      }

    }

    if(spec[i].sizePrecursor()>0){
      foundPre++;
      specCounts++;
      peakCounts+=spec[i].size();

      //build singletList
      spec[i].resetSingletList();

    }

  }
 

  //Finalize the progress
  printf("\b\b\b100%%");
  cout << endl;

  cout << "  " << specCounts << " spectra with " << peakCounts << " peaks will be analyzed." << endl;

  //Build mass list - this orders all precursor masses, with an index pointing to the actual
  //array position for the spectrum. This is because all spectra will have more than 1
  //precursor mass
  massList.clear();
  for(i=0;i<spec.size();i++){
    m.index=i;
    for(j=0;j<spec[i].sizePrecursor();j++){
      m.mass=spec[i].getPrecursor(j).monoMass;
      massList.push_back(m);
    }
  }

  //sort mass list from low to high
  qsort(&massList[0],massList.size(),sizeof(mMass),compareMassList);

  if(bScans!=NULL) delete[] bScans;
  bScans = new bool[spec.size()];

  return true;
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
      fprintf(f, "\" link_site=\"%d\" score=\"%.4lf\" mass=\"%.4lf\" matches=\"%d\" con_fragments=\"%d\" additional_mass=\"%.4lf\"/>\n", (int)sc->site+1, sc->simpleScore, sc->mass, sc->match, sc->conFrag, sc->massA);
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

int MData::outputPepXML(PXWSpectrumQuery& sq, MDatabase& db, kResults& r){
  unsigned int i;
  unsigned int j;

  char c;
  char n;
  char score[32];

  string peptide;
  string protein;
  string sequence;
  string tStr;

  mPeptide pep;
  mScoreCard sc;

  PXWSearchHit sh;

  int siteA;

  sh.hit_rank=1;
  sh.peptide=r.peptide1;
  sh.calc_neutral_pep_mass=r.psmMass;
  sh.massdiff=r.psmMass-r.obsMass;
  
  for(i=0;i<r.mods1.size();i++){
    if (r.mods1[i].pos == 0 && r.mods1[i].term) sh.modInfo.mod_nterm_mass = r.mods1[i].mass;
    else if (r.mods1[i].pos == r.peptide1.length() - 1 && r.mods1[i].term) sh.modInfo.mod_cterm_mass = r.mods1[i].mass;
    else sh.modInfo.addMod((int)r.mods1[i].pos + 1, r.mods1[i].mass + aa.getAAMass(r.peptide1[r.mods1[i].pos]), r.mods1[i].mass, true, "param");
  }
  if (r.nTerm1 && aa.getFixedModMass('$')!=0)sh.modInfo.mod_nterm_mass += aa.getFixedModMass('$');
  if (r.cTerm1 && aa.getFixedModMass('%')!=0)sh.modInfo.mod_cterm_mass += aa.getFixedModMass('%');
  sh.modInfo.mod_nterm_mass += aa.getFixedModMass('n');
  sh.modInfo.mod_cterm_mass += aa.getFixedModMass('c');
  if(sh.modInfo.mod_nterm_mass!=0) sh.modInfo.mod_nterm_mass+=1.00782503;
  if(sh.modInfo.mod_cterm_mass!=0) sh.modInfo.mod_cterm_mass+=17.00273963;
  for(i=0;i<r.peptide1.size();i++){
    if(aa.getFixedModMass(r.peptide1[i])>0) {
      sh.modInfo.addMod(i + 1, aa.getAAMass(r.peptide1[i]), aa.getFixedModMass(r.peptide1[i]),false,"param");
    }
  }
  if (r.link1>0) sh.modInfo.addMod(r.link1, r.massB + aa.getAAMass(r.peptide1[r.link1 - 1]), r.massB,true, "adduct");

  sprintf(score,"%.4lf",r.score);
  sh.addScore("magnum_score",score);
  sprintf(score,"%.4lf",r.scoreDelta);
  sh.addScore("delta_score",score);
  sprintf(score,"%.4lf",r.ppm);
  sh.addScore("ppm_error",score);
  sprintf(score, "%.3e", r.eVal);
  sh.addScore("e_value", score);
  sh.addScore("reporter_ions",r.rIon.c_str());

  //Get proteins
  pep = db.getPeptide(r.pep1);
  sh.num_tot_proteins=(int)pep.map->size();
  for(j=0;j<pep.map->size();j++){
    protein="";
    for(i=0;i<db[pep.map->at(j).index].name.size();i++){
      if(params->truncate>0 && i==params->truncate) break;
      protein+=db[pep.map->at(j).index].name[i];
    }
    if(pep.map->at(j).start<1) n='-';
    else n=db[pep.map->at(j).index].sequence[pep.map->at(j).start-1];
    if(pep.map->at(j).stop+1==db[pep.map->at(j).index].sequence.size()) c='-';
    else c=db[pep.map->at(j).index].sequence[pep.map->at(j).stop+1];
    siteA = pep.map->at(j).start+r.link1;
    sh.addProtein(protein, c, n, siteA);
  }
    
  sq.addSearchHit(&sh,NULL,NULL,NULL);
  return (int)sq.sizeSearchHits()-1;
}

void MData::outputPepXML2(PXWSpectrumQuery& sq, int shIndex, kResults& r){
  size_t i;
  PXWModInfo m;
  for (i = 0; i<r.mods1.size(); i++){
    if (r.mods1[i].pos == 0 && r.mods1[i].term) m.mod_nterm_mass = r.mods1[i].mass;
    else if (r.mods1[i].pos == r.peptide1.length() - 1 && r.mods1[i].term) m.mod_cterm_mass = r.mods1[i].mass;
    else m.addMod((int)r.mods1[i].pos + 1, r.mods1[i].mass + aa.getAAMass(r.peptide1[r.mods1[i].pos]), r.mods1[i].mass, true, "param");
  }
  if (r.nTerm1 && aa.getFixedModMass('$') != 0) m.mod_nterm_mass += aa.getFixedModMass('$');
  if (r.cTerm1 && aa.getFixedModMass('%') != 0) m.mod_cterm_mass += aa.getFixedModMass('%');
  m.mod_nterm_mass += aa.getFixedModMass('n');
  m.mod_cterm_mass += aa.getFixedModMass('c');
  if (m.mod_nterm_mass != 0) m.mod_nterm_mass += 1.00782503;
  if (m.mod_cterm_mass != 0) m.mod_cterm_mass += 17.00273963;
  for (i = 0; i<r.peptide1.size(); i++){
    if (aa.getFixedModMass(r.peptide1[i])>0) {
      m.addMod(i + 1, aa.getAAMass(r.peptide1[i]), aa.getFixedModMass(r.peptide1[i]), false, "param");
    }
  }
  if (r.link1>0) m.addMod(r.link1, r.massB + aa.getAAMass(r.peptide1[r.link1 - 1]), r.massB, true, "adduct");
  sq.getSearchHit(shIndex).a->addAltModInfo(m);
}

bool MData::outputPercolator(FILE* f, MDatabase& db, kResults& r, int count){

  unsigned int i;
  unsigned int j;

  string peptide;
  string protein;
  string sequence;
  string tStr;
  string p1,p2;

  mPeptide pep;
  mScoreCard sc;
  mScoreCard sc2;

  //Export Results:
  if(r.decoy) fprintf(f,"D-");
  else fprintf(f,"T-");
  fprintf(f,"%d-%.2f",r.scanNumber,r.rTime);
  if(count>1) fprintf(f,"-%d",count);
  if(r.decoy) fprintf(f,"\t-1");
  else fprintf(f,"\t1");
  if(params->percVersion>2.04) fprintf(f,"\t%d",r.scanNumber);
  fprintf(f,"\t%.4lf",r.score);
  fprintf(f,"\t%.4lf",r.scoreDelta);
  fprintf(f,"\t%.6lf",-log10(r.eVal));
  
  //export charge state
  for (int z = 1; z<8; z++){
    if (r.charge == (int)z) fprintf(f, "\t1");
    else fprintf(f, "\t0");
  }
  if (r.charge>7) fprintf(f, "\t1");
  else fprintf(f, "\t0");

  fprintf(f,"\t%.4lf",r.psmMass);
  fprintf(f,"\t%.4lf",r.ppm);
  p1=r.modPeptide1;
  fprintf(f,"\t%d\t-.%s",(int)r.peptide1.size(),&p1[0]);
  fprintf(f,".-");

  //export proteins
  pep = db.getPeptide(r.pep1);
  for(j=0;j<pep.map->size();j++){
    protein="";
    for(i=0;i<db[pep.map->at(j).index].name.size();i++){
      if(params->truncate>0 && i==params->truncate) break;
      if(db[pep.map->at(j).index].name[i]==' ') protein+='_';
      else protein+=db[pep.map->at(j).index].name[i];
    }
    fprintf(f,"\t%s",&protein[0]);
  }

  fprintf(f,"\n");

  return true;
}

//Function deprecated. Should be excised.
bool MData::outputResults(MDatabase& db, MParams& par){
  size_t i;
  int j,n,d;
  char fName[256];
  char outPath[1024];
  char peptide[256];
  char tmp[16];
  char specID[256];

  mPeptide pep;
  mPeptide pep2;
  mPrecursor precursor;
  mScoreCard tmpSC;
  mScoreCard tmpSC2;

  mEnzymeRules enzyme;
  kResults res;

  PepXMLWriter p;
  pxwMSMSRunSummary rs;
  pxwSampleEnzyme enz;
  PXWSearchSummary ss;
  PXWSpectrumQuery sq;

  bool bBadFiles;
  bool bTarget1;
  //bool bDupe;
  bool bDiag;

  int scoreIndex;
  //int iDupe;

  float topScore;

  string tmpPep1;
  string tmpPep2;
  string outFile;

  FILE* fOut    = NULL;
  FILE* fSingle = NULL;
  FILE* fDiag   = NULL;

  //Open all the required output files.
  bBadFiles=false;
  sprintf(fName,"%s.magnum.txt",params->outFile);
  fOut=fopen(fName,"wt");
  if(fOut==NULL) bBadFiles=true;
  if(params->exportPercolator) {
    sprintf(fName,"%s.perc.single.txt",params->outFile);
    fSingle=fopen(fName,"wt");
    if(fSingle==NULL) bBadFiles=true;
  }
  if(params->exportPepXML) {
    rs.base_name=params->inFile; 
    rs.base_name = rs.base_name.substr(0, rs.base_name.find_last_of('.'));
    outFile=params->outFile;
    if (rs.base_name[0] == '/'){ //unix
      outFile = outFile.substr(outFile.find_last_of("/") + 1, outFile.size());
    } else { //assuming windows
      outFile = outFile.substr(outFile.find_last_of("\\") + 1, outFile.size());
    }
    rs.raw_data=params->ext;
    rs.raw_data_type="raw";
    rs.search_engine="Magnum";
    ss.base_name=rs.base_name;
    ss.search_engine="Magnum";
    ss.search_engine_version=version;
    processPath(params->dbFile,outPath);
    ss.search_database=outPath;
    for(i=0;i<par.xmlParams.size();i++){
      ss.parameters->push_back(par.xmlParams[i]);
    }
    enz.name=params->enzymeName;
    enz.cut.clear();
    enz.no_cut.clear();
    enz.minNumTermini=2;
    enz.maxNumInternalCleavages=params->miscleave;
    enzyme = db.getEnzymeRules();
    for (i = 65; i < 90; i++){
      if (enzyme.cutC[i]) enz.cut += (char)i;
      if (enzyme.exceptN[i]) enz.no_cut += (char)i;
    }
    if (enz.cut.size()>0){
      enz.sense = "C";
    } else {
      enz.sense = "N";
      for (i = 65; i < 90; i++){
        if (enzyme.cutN[i]) enz.cut += (char)i;
        if (enzyme.exceptC[i]) enz.no_cut += (char)i;
      }
    }
    sprintf(fName, "%s.pep.xml", params->outFile);
    if(!p.createPepXML(fName,rs,&enz,&ss)) bBadFiles=true;
  }
  
  if (params->diag->size()>0){ //create diagnostic file if needed
    sprintf(fName, "%s.diag.xml", params->outFile);
    fDiag = fopen(fName, "wt");
    if (fDiag == NULL) bBadFiles = true;
  }

  //check that all output files are valid
  if(bBadFiles){
    cout << "Error exporting results. Please make sure drive is writable." << endl;
    if(fOut!=NULL)    fclose(fOut);
    if(fSingle!=NULL) fclose(fSingle);
    if(fDiag!=NULL)   fclose(fDiag);
    return false;
  }

  //Put the headers on all the files
  fprintf(fOut,"Magnum version %s\n",version);
  fprintf(fOut,"Scan Number\tRet Time\tObs Mass\tCharge\tPSM Mass\tPPM Error\tScore\tdScore\te-value\tReporter Ions\tPeptide #1 Score\tPeptide #1\tLinked AA #1\tProtein #1\tProtein #1 Site\tPeptide #2 Score\tPeptide #2\tLinked AA #2\tProtein #2\tProtein #2 Site\tLinker Mass\n");
  if(params->exportPercolator){
    if(params->percVersion>2.04) fprintf(fSingle,"SpecId\tLabel\tscannr\tScore\tdScore\tnegLog10eVal\t");
    else fprintf(fSingle,"SpecId\tLabel\tScore\tdScore\tnegLog10eVal\t");
    fprintf(fSingle,"Charge1\tCharge2\tCharge3\tCharge4\tCharge5\tCharge6\tCharge7\tCharge8Plus\t");
    fprintf(fSingle,"Mass\tPPM\tLen\tPeptide\tProteins\n");
  }
  if(fDiag!=NULL){
    fprintf(fDiag,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fDiag,"<magnum_analysis date=\"now\">\n");
  }

  //Output top score for each spectrum
  //Must iterate through all possible precursors for that spectrum
  for(i=0;i<spec.size();i++) {

    //Check if we need to output diagnostic information
    bDiag=false;
    if(params->diag->size()>0){
      if(params->diag->at(0)==-1) {
        bDiag=true;
      } else {
        for (d = 0; d<params->diag->size(); d++){
          if (spec[i].getScanNumber() == params->diag->at(d)){
            bDiag=true;
            break;
          }
        }
      }
    }
    if(bDiag) outputDiagnostics(fDiag,spec[i],db);

    //Iterate over the top hits to find subgroupings.
    //Subgroupings are defined as same sequence, but alternate interpretations
    //Each novel subgroup has the same top score, but a unique sequence
    typedef struct pepXMLGroup {
      int shIndex;
      int pepIndex;
      vector<int> indexes;
    } pepXMLGroup;
    pepXMLGroup px;
    px.shIndex=-1;
    vector<pepXMLGroup> pxg;
    size_t q;
    topScore=spec[i].getScoreCard(0).simpleScore;
    if(topScore>0){
      for(j=0;j<20;j++){
        tmpSC = spec[i].getScoreCard(j);
        if(tmpSC.simpleScore<topScore) break;
        for(q=0;q<pxg.size();q++){
          if(pxg[q].pepIndex==tmpSC.pep) break;
        }
        if(q==pxg.size()) { //add new entry
          px.pepIndex=tmpSC.pep;
          px.indexes.clear();
          px.indexes.push_back(j);
          pxg.push_back(px);
        } else {
          pxg[q].indexes.push_back(j);
        }
      }
    }

    scoreIndex=0;
    tmpSC=spec[i].getScoreCard(scoreIndex);
    res.scanNumber=spec[i].getScanNumber();
    res.rTime=spec[i].getRTime();

    //if there are no matches to the spectrum, return null result and continue
    if(tmpSC.simpleScore==0){
      fprintf(fOut,"%d\t%.4f\t0\t0\t0\t0\t0\t0\t999.0\t-\t0\t-\t-\t-\t-\t0\t-\t-\t-\t-\t0\n",res.scanNumber,res.rTime);
      continue;
    }

    if(params->exportPepXML){
      sq.clear();
      sq.end_scan=res.scanNumber;
      sq.retention_time_sec=res.rTime;
      sq.start_scan=res.scanNumber;
    }

    //Check for Reporter Ions
    res.rIon.clear();
    for(j=0;j<params->rIons->size();j++){
      if(spec[i].checkReporterIon(params->rIons->at(j),params)){
        sprintf(tmp, "%.3lf", params->rIons->at(j));
        if(res.rIon.size()>0) res.rIon+=",";
        res.rIon+=tmp;
      }
    }
    if(res.rIon.size()==0) res.rIon="-";

    //Export top scoring peptide, plus any ties that occur after it.
    topScore=tmpSC.simpleScore;
    int count=0;
    while(tmpSC.simpleScore==topScore){

      count++;

      //Get precursor ion for the PSM
      precursor=spec[i].getPrecursor((int)tmpSC.precursor);
      res.obsMass = precursor.monoMass;
      res.charge  = precursor.charge;
      res.ppm = (tmpSC.mass +tmpSC.massA- precursor.monoMass) / precursor.monoMass*1e6;
      res.psmMass = tmpSC.mass+tmpSC.massA;
      res.hk = precursor.corr;

      if(params->exportPepXML){
        sq.assumed_charge=res.charge;
        sq.precursor_neutral_mass=res.obsMass;
        sprintf(specID,"%s.%d.%d.%d",&outFile[0],res.scanNumber,res.scanNumber,res.charge);
        sq.spectrum=specID;
      }

      //grab the next highest score that matches to the same precursor ion for the delta score
      //do not count ties - look for the first difference
      //if no other match has the same precursor, just take the lowest score in the list
      n=scoreIndex+1;
      while(n<19){
        tmpSC2=spec[i].getScoreCard(n++);
        if(tmpSC2.simpleScore==0) break;
        if(tmpSC2.simpleScore==topScore) continue;
        if(tmpSC2.precursor!=tmpSC.precursor) continue;

        //if peptides and link sites are the same, go to the next one
        //this no longer applies to the top result. duplicates may occur among lower results
        //if(tmpSC.link>-1 && tmpSC2.link>-1 && tmpSC2.pep1==tmpSC.pep1 && tmpSC2.pep2==tmpSC.pep2 && tmpSC2.k1==tmpSC.k1 && tmpSC2.k2==tmpSC.k2){
        //  cout << "Oddity 1: " << spec[i].getScanNumber() << endl;
        //  continue;
        //}
        break;
      }
      res.score       = tmpSC.simpleScore;
      res.scoreDelta  = tmpSC.simpleScore-tmpSC2.simpleScore;
      res.scorePepDif = 0;
      res.eVal        = tmpSC.eVal;
      
      MTopPeps* tp = spec[i].getTopPeps((int)tmpSC.precursor);
      mScoreCard* grr;
      grr = tp->peptideFirst;
      int rank=1;
      res.rankA=0;
      res.rankB=0;
      res.rank        = res.rankA+res.rankB;
      res.scoreA      = 0;
      res.scoreB      = 0;
      res.massA       = tmpSC.mass;
      res.massB       = tmpSC.massA;

      //Get the peptide sequence(s)
      pep = db.getPeptide(tmpSC.pep);
      db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
      res.peptide1 = peptide;
      res.mods1.clear();
      res.cTerm1 = pep.cTerm;
      res.nTerm1 = pep.nTerm;
      for(j=0;j<tmpSC.mods->size();j++) res.mods1.push_back(tmpSC.mods->at(j));
      res.peptide2 = "";
      res.modPeptide1 = processPeptide(pep,tmpSC.mods,tmpSC.site,tmpSC.massA,db);
      res.modPeptide2 = "";

      //Get the link positions - relative to the peptide
      res.link1 = tmpSC.site;
      if(res.link1>=0) res.link1++;

      //set link type
      res.type=0;

      //Get the peptide indexes
      res.pep1 = tmpSC.pep;
      res.linkable1 = false;

      //Edge case where single peptide is shared between linked and non-linked peptide lists
      //This occurs when the peptide appears multiple times in a database: internally and on
      //the c-terminus for amine reactive cross-linkers, for example.
      //bDupe=false;
      /*
      if(res.type==0){
        n=scoreIndex+1;
        iDupe=1;
        while(n<19){
          iDupe++;
          tmpSC2=spec[i].getScoreCard(n++);
          if(tmpSC2.simpleScore==0) break;
          if(tmpSC2.simpleScore!=topScore) break;

          //if peptides are the same, but different lists (linked vs. non), use second peptide as location
          if(tmpSC2.linkable1!=tmpSC.linkable1) {
            pep = db.getPeptide(res.pep1);
            db.getPeptideSeq(pep,tmpPep1);
            pep2 = db.getPeptide(tmpSC2.pep1);
            db.getPeptideSeq(pep2,tmpPep2);
            if(tmpPep1.compare(tmpPep2)==0){
              res.pep2=tmpSC2.pep1;
              res.linkable2=tmpSC2.linkable1;
              bDupe=true;
              break;
            }
          }
        }
      }
      */

      //Determine if target or decoy
      bTarget1=false;
      pep = db.getPeptide(res.pep1);
      for(j=0;j<pep.map->size();j++) if(db[pep.map->at(j).index].name.find(params->decoy)==string::npos) bTarget1=true;
      if(bTarget1) res.decoy=false;
      else res.decoy=true;
      tmpPep1=res.peptide1;
      sprintf(tmp,"(%d)",res.link1);
      tmpPep1+=tmp;

      //Export Results:
      fprintf(fOut,"%d",res.scanNumber);
      //fprintf(fOut, "\t%.4lf",res.hk); //this was for diagnostics of hardklor correlation results (or lack of)
      fprintf(fOut,"\t%.4f",res.rTime);
      fprintf(fOut,"\t%.4lf",res.obsMass);
      fprintf(fOut,"\t%d",res.charge);
      fprintf(fOut,"\t%.4lf",res.psmMass);
      fprintf(fOut,"\t%.4lf",res.ppm);
      fprintf(fOut,"\t%.4lf",res.score);
      fprintf(fOut,"\t%.4lf",res.scoreDelta);
      fprintf(fOut,"\t%.3e",res.eVal);
      fprintf(fOut,"\t%s",res.rIon.c_str());
      //fprintf(fOut,"\t%.4lf",res.scorePepDif);
      if (res.scoreA == 0)fprintf(fOut, "\t%.4lf", res.score);
      else fprintf(fOut,"\t%.4lf",res.scoreA);
      fprintf(fOut,"\t%s",&res.modPeptide1[0]);
      if(res.n15Pep1) fprintf(fOut,"-15N");
      fprintf(fOut,"\t%d",res.link1);

      //export protein
      fprintf(fOut,"\t");
      tmpPep1.clear();
      pep = db.getPeptide(res.pep1);
      for(j=0;j<pep.map->size();j++){
        if(j>0) fprintf(fOut,";");
        fprintf(fOut,"%s",&db[pep.map->at(j).index].name[0]);
        if(res.link1>=0) {
          sprintf(tmp, "%d", pep.map->at(j).start + res.link1);
          if(tmpPep1.size()>0) tmpPep1+=";";
          tmpPep1+=tmp;
          //fprintf(fOut,"(%d);",pep.map->at(j).start+res.link1); //put position from start of protein
        }

      }
      //if(bDupe){
      //  pep = db.getPeptide(res.pep2);
      //  for(j=0;j<pep.map->size();j++){
      //    fprintf(fOut,"%s;",&db[pep.map->at(j).index].name[0]);
      //    //if(res.link1>=0) fprintf(fOut,"(%d);",pep.map->at(j).start+res.link1); //only non-linked peptides
      //  }
      //}
      if(res.link1>-1) fprintf(fOut,"\t%s",&tmpPep1[0]);
      else fprintf(fOut,"\t-");
      
      fprintf(fOut,"\n");
      
      if(params->exportPercolator) outputPercolator(fSingle,db,res,count);
      /*if(params->exportPepXML) {
        for(q=0;q<pxg.size();q++){
          size_t qq;
          for(qq=0;qq<pxg[q].indexes.size();qq++){
            if(pxg[q].indexes[qq]==scoreIndex) break;
          }
          if (qq == pxg[q].indexes.size()) continue;
          if (pxg[q].shIndex<0) pxg[q].shIndex = outputPepXML(sq, db, res);
          else outputPepXML2(sq, pxg[q].shIndex, res);
          break;
        }
      }*/
      if (params->exportPepXML) outputPepXML(sq, db, res);

      //Get the next entry - it must also be exported if it has the same score
      //if(bDupe) scoreIndex+=iDupe;
      //else scoreIndex++;
      scoreIndex++;
      if(scoreIndex>=20) break;
      tmpSC=spec[i].getScoreCard(scoreIndex);
    }

    if(params->exportPepXML)p.writeSpectrumQuery(sq);

  }

  fclose(fOut);
  if(params->exportPercolator) {
    fclose(fSingle);
  }
  if(params->exportPepXML){
    p.closePepXML();
  }
  if (fDiag != NULL){
    fprintf(fDiag, "</magnum_analysis>\n");
    fclose(fDiag);
  }

  return true;

}

//Reads in raw/mzXML/mzML files. Other formats supported in MSToolkit as well.
bool MData::readSpectra(){

  MSReader   msr;
  Spectrum   s;
  Spectrum   c;
  //MSpectrum  pls(params->topCount,params->binSize,params->binOffset,params->threads);
  MSpectrum pls(*params);
  mSpecPoint sp;
  float      max;
  mPrecursor pre;

  int totalScans=0;
  int totalPeaks=0;
  int collapsedPeaks=0;
  int finalPeaks=0;
  int iPercent=0;
  int iTmp;

  int i;
  int j;

  spec.clear();
  msr.setFilter(MS2);

  //Set progress meter
  printf("%2d%%", iPercent);
  fflush(stdout);

  if(!msr.readFile(params->msFile,s)) return false;
  while(s.getScanNumber()>0){

    totalScans++;
    if(s.size()<1) {
      msr.readFile(NULL,s);
      continue;
    }

    //This is for the methods used in old grad school data
    /*
    if(s.getRTime()<15){
      msr.readFile(NULL,s);
      continue;
    }
    */

    pls.clear();
    pls.setRTime(s.getRTime());
    pls.setScanNumber(s.getScanNumber());
    max=0;

    //Check whether scans are centroided or not (info supplied by user in params)
    if(!params->ms2Centroid) {

      //If not centroided, do so now.
      centroid(s,c,params->ms2Resolution,params->instrument);

      totalPeaks+=c.size();
      
      //Collapse the isotope peaks
      if(params->specProcess==1 && c.size()>1) {
        collapseSpectrum(c);
        collapsedPeaks+=c.size();
      }

      //If user limits number of peaks to analyze, sort by intensity and take top N
      if(params->maxPeaks>0){
        if(c.size()>1) c.sortIntensityRev();
        if(c.size()<params->maxPeaks) j=c.size();
        else j=params->maxPeaks;
      } else {
        j=c.size();
      }
      for(i=0;i<j;i++){
        sp.mass=c[i].mz;
        sp.intensity=c[i].intensity;
        pls.addPoint(sp);
        if(sp.intensity>max) max=sp.intensity;
      }
      pls.setMaxIntensity(max);

      //Sort again by MZ, if needed
      if(pls.size()>1 && params->maxPeaks>0) pls.sortMZ();

      finalPeaks+=pls.size();

    } else {

      //Collapse the isotope peaks
      if(params->specProcess==1 && s.size()>1) collapseSpectrum(s);

      //If user limits number of peaks to analyze, sort by intensity and take top N
      if(params->maxPeaks>0){
        if(s.size()>1) s.sortIntensityRev();
        if(s.size()<params->maxPeaks) j=s.size();
        else j=params->maxPeaks;
      } else {
        j=s.size();
      }
      for(i=0;i<j;i++){
        sp.mass=s[i].mz;
        sp.intensity=s[i].intensity;
        pls.addPoint(sp);
        if(sp.intensity>max) max=sp.intensity;
      }
      pls.setMaxIntensity(max);
      
      //Sort again by MZ, if needed
      if(pls.size()>1 && params->maxPeaks>0) pls.sortMZ();
    }

    //Get any additional information user requested
    pls.setCharge(s.getCharge());
    pls.setMZ(s.getMZ());
    if(params->preferPrecursor>0){
      if(s.getMonoMZ()>0 && s.getCharge()>0){
        pre.monoMass=s.getMonoMZ()*s.getCharge()-s.getCharge()*1.007276466;
        pre.charge=s.getCharge();
        pre.corr=0;
        pls.addPrecursor(pre,params->topCount);
        pls.setInstrumentPrecursor(true);
      } else if(s.sizeZ()>0){
        for(j=0; j<s.sizeZ(); j++){
          pre.monoMass=s.atZ(j).mh-1.007276466;
          pre.charge=s.atZ(j).z;
          pre.corr=-5;
          pls.setCharge(pre.charge);
          pls.addPrecursor(pre, params->topCount);
          pls.setInstrumentPrecursor(true);
        }
      }
    }

    //Add spectrum (if it has enough data points) to data object and read next file
    if(pls.size()>=params->minPeaks) spec.push_back(pls);

    for(unsigned int d=0;d<params->diag->size();d++){
      if(pls.getScanNumber()==params->diag->at(d)){
        char diagStr[256];
        sprintf(diagStr,"diagnostic_spectrum_%d.txt",params->diag->at(d));
        FILE* f=fopen(diagStr,"wt");
        fprintf(f,"Scan: %d\t%d\n",pls.getScanNumber(),pls.size());
        for(int k=0;k<pls.size();k++) fprintf(f,"%.6lf\t%.0f\n",pls[k].mass,pls[k].intensity);
        fclose(f);
        break;
      }
    }

    //Update progress meter
    iTmp = msr.getPercent();
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

    msr.readFile(NULL,s);
  }

  //Finalize progress meter
  if(iPercent<100) printf("\b\b\b100%%");
  cout << endl;

  cout << "  " << spec.size() << " total spectra have enough data points for searching." << endl;
  //cout << totalScans << " total scans were loaded." <<  endl;
  //cout << totalPeaks << " total peaks in original data." << endl;
  //cout << collapsedPeaks << " peaks after collapsing." << endl;
  //cout << finalPeaks << " peaks after top N." << endl;
	return true;
}

void MData::setAdductSites(char* s){
  for(size_t i=0;i<strlen(s);i++) adductSite[s[i]]=true;
}

void MData::setVersion(const char* v){
  strcpy(version,v);
}

int MData::size(){
  return (int)spec.size();
}

void MData::xCorr(){
  cout << "  Transforming spectra ... ";

  //Threading is a little faster, but ultimately memory allocation is a large portion of the time in these operations.
  ThreadPool<MSpectrum*>* threadPool = new ThreadPool<MSpectrum*>(xCorrProc, params->threads, params->threads, 1);

  int iTmp;
  int iPercent = 0;
  printf("%2d%%", iPercent);
  fflush(stdout);
  for(size_t i=0;i<spec.size();i++) {
    //spec[i].xCorrScore(b);

    threadPool->WaitForQueuedParams();

    MSpectrum* a = &spec[i];
    threadPool->Launch(a);

    //Update progress meter
    iTmp = (int)((double)i / spec.size() * 100);
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

  }

  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //Finalize progress meter
  if(iPercent<100) printf("\b\b\b100%%");
  cout << endl;

  //clean up memory & release pointers
  delete threadPool;
  threadPool = NULL;
}

void MData::xCorrProc(MSpectrum* s){
  s->xCorrScore();
  s = NULL;
}

/*============================
  Private Utilities
============================*/
//First derivative method, returns base peak intensity of the set
void MData::centroid(Spectrum& s, Spectrum& out, double resolution, int instrument){
  int i,j;
  float maxIntensity;
  int bestPeak;
  bool bLastPos;

	int nextBest;
	double FWHM;
  double maxMZ = s[s.size()-1].mz+1.0;
	Peak_T centroid;

	vector<double> x;
	vector<double> y;
	vector<double> c;
	int left, right;
	bool bPoly;
	float lastIntensity;

	out.clear();

  bLastPos=false;
	for(i=0;i<s.size()-1;i++){

    if(s[i].intensity<s[i+1].intensity) {
      bLastPos=true;
      continue;
    } else {
      if(bLastPos){
				bLastPos=false;

				//find max and add peak
				maxIntensity=0;
				for(j=i;j<i+1;j++){
				  if (s[j].intensity>maxIntensity){
				    maxIntensity=s[j].intensity;
				    bestPeak = j;
				  }
				}

				//walk left and right to find bounds above half max
				left=right=bestPeak;
				lastIntensity=maxIntensity;
				for(left=bestPeak-1;left>0;left--){
					if(s[left].intensity<(maxIntensity/3) || s[left].intensity>lastIntensity){
						left++;
						break;
					}
					lastIntensity=s[left].intensity;
				}
				lastIntensity=maxIntensity;
				for(right=bestPeak+1;right<s.size()-1;right++){
					if(s[right].intensity<(maxIntensity/3) || s[right].intensity>lastIntensity){
						right--;
						break;
					}
					lastIntensity=s[right].intensity;
				}

				//if we have at least 5 data points, try polynomial fit
				double r2;
				bPoly=false;
				if((right-left+1)>4){
					x.clear();
					y.clear();
					for(j=left;j<=right;j++){
						x.push_back(s[j].mz);
						y.push_back(log(s[j].intensity));
					}
					r2=polynomialBestFit(x,y,c);
					if(r2>0.95){
						bPoly=true;
						centroid.mz=-c[1]/(2*c[2])+c[3];
						centroid.intensity=(float)exp(c[0]-c[2]*(c[1]/(2*c[2]))*(c[1]/(2*c[2])));
					} else {

					}
				}

				if(!bPoly){
					//Best estimate of Gaussian centroid
					//Get 2nd highest point of peak
					if(bestPeak==s.size()) nextBest=bestPeak-1;
					else if(s[bestPeak-1].intensity > s[bestPeak+1].intensity) nextBest=bestPeak-1;
					else nextBest=bestPeak+1;

					//Get FWHM
					switch(instrument){
						case 0: FWHM = s[bestPeak].mz*sqrt(s[bestPeak].mz)/(20*resolution); break;  //Orbitrap
						case 1: FWHM = s[bestPeak].mz*s[bestPeak].mz/(400*resolution); break;				//FTICR
						default: break;
					}

					//Calc centroid MZ (in three lines for easy reading)
					centroid.mz = pow(FWHM,2)*log(s[bestPeak].intensity/s[nextBest].intensity);
					centroid.mz /= GAUSSCONST*(s[bestPeak].mz-s[nextBest].mz);
					centroid.mz += (s[bestPeak].mz+s[nextBest].mz)/2;

					//Calc centroid intensity
					centroid.intensity=(float)(s[bestPeak].intensity/exp(-pow((s[bestPeak].mz-centroid.mz)/FWHM,2)*GAUSSCONST));
				}

				//some peaks are funny shaped and have bad gaussian fit.
				//if error is more than 10%, keep existing intensity
				if( fabs((s[bestPeak].intensity - centroid.intensity) / centroid.intensity * 100) > 10 ||
            //not a good check for infinity
            centroid.intensity>9999999999999.9 ||
            centroid.intensity < 0 ) {
					centroid.intensity=s[bestPeak].intensity;
				}

				//Hack until I put in mass ranges
				if(centroid.mz<0 || centroid.mz>maxMZ) {
					//do nothing if invalid mz
				} else {
					out.add(centroid);
				}
			
      }

    }
  }

}

//Function tries to remove isotopes of signals by stacking the intensities on the monoisotopic peak
//Also creates an equal n+1 peak in case wrong monoisotopic peak was identified.
void MData::collapseSpectrum(Spectrum& s){
  int i,j,k,n;
  int charge,z;
  int maxIndex;
  float max;
  float cutoff;
  vector<int> dist;

  Spectrum s2;

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
    while(j<s.size() && (s[j].mz-s[maxIndex].mz)<1.1){
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
      while(n<s.size() && (s[n].mz-s[k].mz)<1.1){
        if(s[n].intensity<1) {
          n++;
          continue;
        }
        z=getCharge(s,k,n);
        if(z>0 && z<charge) {
          break;
        } else if(z==charge && (s[n].mz-s[k].mz)>(0.98/charge) && (s[n].mz-s[k].mz)<(1.02/charge)) {
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
      s2.add(s[dist[0]]);
      s[dist[0]].intensity=0;
      continue;
    }

    //step to the left
    j=maxIndex-1;
    while(j>=0 && (s[maxIndex].mz-s[j].mz)<1.1){
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
      while(n>=0 && (s[k].mz-s[n].mz)<1.1){
        if(s[n].intensity<1) {
          n--;
          continue;
        }
        z=getCharge(s,n,k);
        //printf("\tleft\t%.6lf\t%.6lf\t%d\n",s[n].mz,s[k].mz-s[n].mz,z);
        if(z>0 && z<charge) {
          break;
        } else if(z==charge && s[k].mz-s[n].mz > 0.98/charge && s[k].mz-s[n].mz < 1.02/charge) {
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
        max=s[dist[0]].intensity+s[dist[1]].intensity;
        s2.add(s[dist[0]].mz,max);
        s[dist[1]].intensity=0;
       // s2.add(s[dist[1]].mz,max);
      } else {
        s2.add(s[dist[0]]);
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
      s2.add(s[j].mz,max);
      //s2.add(s[k].mz,max);
    }

  }

  s2.sortMZ();
  s.clearPeaks();
  for(i=0;i<s2.size();i++) {
    if(i<s2.size()-1 && s2[i].mz==s2[i+1].mz){
      if(s2[i].intensity>s2[i+1].intensity) s.add(s2[i]);
      else s.add(s2[i+1]);
      i++;
    } else {
      s.add(s2[i]);
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

int MData::getCharge(Spectrum& s, int index, int next){
  double mass;

  mass=s[next].mz-s[index].mz;
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

string MData::processPeptide(mPeptide& pep, vector<mPepMod>* mod, int site, double massA, MDatabase& db){
  char tmp[32];
  size_t j,k;
  string seq = "";
  string peptide;

  db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, peptide);

  if (pep.nTerm && aa.getFixedModMass('$') != 0) {
    sprintf(tmp, "n[%.0lf]", aa.getFixedModMass('$'));
    seq += tmp;
  }
  for (k = 0; k<mod->size(); k++){ //check for n-terminal peptide mod
    if (mod->at(k).pos == 0 && mod->at(k).term){
      sprintf(tmp, "n[%.0lf]", mod->at(k).mass);
      seq += tmp;
    }
  }
  for (j = 0; j<peptide.size(); j++) {
    seq += peptide[j];
    for (k = 0; k<mod->size(); k++){
      if (j == (unsigned int)mod->at(k).pos && !mod->at(k).term){
        sprintf(tmp, "[%.0lf]", mod->at(k).mass);
        seq += tmp;
      }
    }
    if(j==(size_t)site){
      sprintf(tmp, "[%.0lf]", massA);
      seq += tmp;
    }
  }
  for (k = 0; k<mod->size(); k++){ //check for c-terminal peptide mod
    if (mod->at(k).pos > 0 && mod->at(k).term){
      sprintf(tmp, "c[%.0lf]", mod->at(k).mass);
      seq += tmp;
    }
  }
  if (pep.cTerm && aa.getFixedModMass('%') != 0) {
    sprintf(tmp, "c[%.0lf]", aa.getFixedModMass('%'));
    seq += tmp;
  }

  return seq;
}
