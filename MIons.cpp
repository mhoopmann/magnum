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

#include "MIons.h"

using namespace std;

MIons::MIons(){
  int i;
  for(i=0;i<128;i++){
    aaMass[i]=0;
    aaFixedModMass[i]=0;
    aaMod[i].count=0;
    site[i]=false;
  }

  aaMass['A']=71.0371103;
  aaMass['C']=103.0091803;
  aaMass['D']=115.0269385;
  aaMass['E']=129.0425877;
  aaMass['F']=147.0684087;
  aaMass['G']=57.0214611;
  aaMass['H']=137.0589059;
  aaMass['I']=113.0840579;
  aaMass['K']=128.0949557;
  aaMass['L']=113.0840579;
  aaMass['M']=131.0404787;
  aaMass['N']=114.0429222;
  aaMass['P']=97.0527595;
  aaMass['Q']=128.0585714;
  aaMass['R']=156.1011021;
  aaMass['S']=87.0320244;
  aaMass['T']=101.0476736;
  aaMass['V']=99.0684087;
  aaMass['U']=150.9536303;
  aaMass['W']=186.0793065;
  aaMass['Y']=163.0633228;
  aaMass['c']=0;
  aaMass['n']=0;
  aaMass['$']=0;
  aaMass['%']=0;

  modList=NULL;
  pep1=NULL;
  maxModCount=0;
  ionCount=0;
}

MIons::~MIons(){
  if(modList!=NULL) delete [] modList;
  pep1=NULL;
}

void MIons::addFixedMod(char mod, double mass){
  aaMass[mod]+=mass;
  aaFixedModMass[mod]=mass;
}

void MIons::addMod(char mod, bool xl, double mass){
  aaMod[mod].mod[aaMod[mod].count].xl=xl;
  aaMod[mod].mod[aaMod[mod].count++].mass=mass;
}

void MIons::buildModIons2(bool bAdduct) {
  //cout << "buildModIons2: " << (int)bAdduct << endl;
  double mMass;

  mPrecursor.clear();
  pepLinks.clear();
  pepMass.clear();
  pepMods.clear();
  maxLink = -99;

  vPeaks.clear();
  mPeaks.clear();
  vPeaksRev.clear();
  mPeaksRev.clear();

  //set up boundaries
  ionCount = pep1Len;

  //peptide number
  pepCount = 0;

  modMask=string(pep1Len+2,'0');
  mPrecursor.insert(pair<string, size_t>(modMask, pepCount));

  //b- & y-ions from first peptide
  pepMass.push_back(pep1Mass);
  pepMassMin = pepMassMax = pep1Mass;
  mMass = aaMass['n'];
  pepMass[0] += aaMass['n'];
  if (nPep1) {
    mMass += aaMass['$'];
    pepMass[0] += aaMass['$'];
  }
  pepMods.emplace_back();
  if(bAdduct) {
    //cout << "Adduct: " << vPeaks.size() << "\t" << vPeaksRev.size() << endl;
    //TODO: Combine modIonsMaskRec and modIonsNew into single function.
    //So that you don't have to clear and repeat as is done below.
    modIonsMaskRec(0,mMass,-99,pepCount,0,-1,modMask);
    vPeaks.clear();
    mPeaks.clear();

    //displayPrecursors();
    //if (pepseq.compare("NFPPSQDASGDLYTTSSQLTLPATQCLAGK") == 0) displayPrecursors();
    mMass = pep1Mass+aaMass['n'];
    if(nPep1) mMass += aaMass['$'];
    map<string, size_t>::iterator it = mPrecursor.begin();
    while (it != mPrecursor.end()) {
      modIonsNew(it->first, it->second, mMass);
      it++;
    }
  } else {
    //cout << "NoAdduct" << endl;
    modIonsMaskRecNoAdduct(0, mMass, pepCount, 0, -1, modMask);
    //displayPrecursors();
    //if (pepseq.compare("NFPPSQDASGDLYTTSSQLTLPATQCLAGK") == 0) displayPrecursors();
    mMass = pep1Mass + aaMass['n'];
    if (nPep1) mMass += aaMass['$'];
    map<string, size_t>::iterator it = mPrecursor.begin();
    while (it != mPrecursor.end()) {
      modIonsNew(it->first, it->second, mMass,3);
      it++;
    }
  }
  pepCount++; //Double-check this?

  //build y-ions here
  //use same mMass?
  //mMass = 18.0105633 + aaMass['c'];
  //if (cPep1) mMass += aaMass['$']; //'%' instead??
  if (bAdduct) {
    map<string, size_t>::iterator it = mPrecursor.begin();
    while (it != mPrecursor.end()) {
      modIonsRecNew(it->first, it->second,pepMass[it->second]);
      it++;
    }
  } else {
    map<string, size_t>::iterator it = mPrecursor.begin();
    while (it != mPrecursor.end()) {
      modIonsRecNew(it->first, it->second,pepMass[it->second],3);
      it++;
    }
  }

}

//For diagnostics. Never called otherwise.
void MIons::displayPrecursors() {
  map<string, size_t>::iterator it = mPrecursor.begin();
  while (it != mPrecursor.end()) {
    cout << it->first << "\t" << it->second << endl;
    cout << "\t" << pepMass[it->second] << endl;
    if (pepMods[it->second].mods.size() > 0) {
      for (size_t a = 0;a < pepMods[it->second].mods.size();a++) {
        cout << "\t" << (int)pepMods[it->second].mods[a].pos << "," << pepMods[it->second].mods[a].mass << "," << (int)pepMods[it->second].mods[a].term;
      }
      cout << endl;
    }
    it++;
  }
}

void MIons::modIonsNew(const string& mask, size_t pepIndex, double mass, size_t stop) {
  //cout << "modIonsNew: " << mask << " " << pepIndex << " " << mass << endl;
  double m = mass;
  bool bAdduct = false;
  int iAdduct = -1;

  //add all modification masses
  for (int a = 0; a < mask.size() - 2; a++) {
    if (mask[a] != '0') {
      //cout << mask << " has mod at: " << a << "  new mass=" << m << " val=" << mask[a] - 49 << " aa:" << pep1[a] << " wtf:" << aaMod[pep1[a]].mod[mask[a] - 49].mass << endl;
      m += aaMod[pep1[a]].mod[mask[a] - 49].mass;
    }
  }
  double trueMass = m;
  m = 0;

  //process n-term here

  //process sequence
  for (int a = 0; a < mask.size() - stop; a++) {
    if (mask[a] == 'x') bAdduct = true;
    if (mask[a] != '0') m += aaMod[pep1[a]].mod[mask[a] - 49].mass;
    m += aaMass[pep1[a]];
    if(bAdduct) addPeakNew(-m, trueMass, pepIndex);
    else addPeakNew(m, trueMass, pepIndex);
  }

  //process c-term here
}

void MIons::modIonsRecNew(const string& mask, size_t pepIndex, double mass, size_t stop){
  //cout << "modIonsRecNew: " << mask << " " << pepIndex << " " << mass << endl;
  double m=mass;
  int iAdduct=-1;

  //add all modification masses
  for (int a = 0; a < mask.size() - 2; a++) {
    if(mask[a]=='x') iAdduct=a;
    //else if(mask[a]!='0'){
    //  cout << mask << " has mod at: " << a << "  new mass=" << m << " val=" << mask[a] - 49 << " aa:" << pep1[a] << " wtf:" << aaMod[pep1[a]].mod[mask[a] - 49].mass << endl;
    //  m+=aaMod[pep1[a]].mod[mask[a]-49].mass;
    //}
  }
  double trueMass=m;
  //cout << "trueMass: " << m << endl;

  //process n-term here

  //process sequence
  for(int a=0;a<mask.size()-stop;a++){
    if(mask[a]!='0' && mask[a]!='x'){
      if (a<iAdduct)addPeakRevNew((aaMass[pep1[a]]+aaMod[pep1[a]].mod[mask[a] - 49].mass) - m, trueMass, pepIndex);
      else {
        //cout << "Add peak1: " << m << " - " << (aaMass[pep1[a]] + aaMod[pep1[a]].mod[mask[a] - 49].mass) << " " << trueMass << endl;
        addPeakRevNew(m - (aaMass[pep1[a]] + aaMod[pep1[a]].mod[mask[a] - 49].mass), trueMass, pepIndex);
      }
      m -= (aaMass[pep1[a]]+ aaMod[pep1[a]].mod[mask[a] - 49].mass);
    } else {
      if(a <iAdduct)addPeakRevNew(aaMass[pep1[a]]-m, trueMass, pepIndex);
      else {
        //cout << "Add peak2: " << m << " - " << aaMass[pep1[a]] << " " << trueMass << endl;
        addPeakRevNew(m - aaMass[pep1[a]], trueMass, pepIndex);
      }
      m-=aaMass[pep1[a]];
    }
  }

  //process c-term here
}

//Abstract this and combine with addPeakRevNew
void MIons::addPeakNew(double mass, double pepMass, size_t pepIndex) {
  //cout << "addPeakNew: " << mass << "\t" << pepMass << "\t" << pepIndex << endl;
  map<int, size_t>::iterator it;
  map<int, size_t>::iterator it2;
  int peak = (int)(mass * 10);
  it = mPeaks.find(peak);
  if (it == mPeaks.end()) {
    mPeaks.insert(pair<int, size_t>(peak, vPeaks.size()));
    vPeaks.emplace_back();
    vPeaks.back().mass = mass;
    vPeaks.back().index.emplace_back();
    if (mass < 0) {
      vPeaks.back().index.back().pepMass = pepMass;
    }
    vPeaks.back().index.back().pepIndex.push_back(pepIndex);
  } else {
    if (mass < 0) {
      size_t a;
      for (a = 0; a < vPeaks[it->second].index.size(); a++) {
        if (vPeaks[it->second].index[a].pepMass == pepMass) break;
      }
      if (a == vPeaks[it->second].index.size()) {
        vPeaks[it->second].index.emplace_back();
        vPeaks[it->second].index[a].pepMass = pepMass;
      }
      vPeaks[it->second].index[a].pepIndex.push_back(pepIndex);
    } else vPeaks[it->second].index[0].pepIndex.push_back(pepIndex);
  }
}

void MIons::addPeakRevNew(double mass, double pepMass, size_t pepIndex){
  map<int,size_t>::iterator it;
  map<int, size_t>::iterator it2;
  int peak=(int)(mass*10);
  it=mPeaksRev.find(peak);
  if(it==mPeaksRev.end()){
    mPeaksRev.insert(pair<int,size_t>(peak,vPeaksRev.size()));
    vPeaksRev.emplace_back();
    vPeaksRev.back().mass=mass;
    vPeaksRev.back().index.emplace_back();
    if(mass<0){
      vPeaksRev.back().index.back().pepMass=pepMass;
    }
    vPeaksRev.back().index.back().pepIndex.push_back(pepIndex);
  } else {
    if(mass<0){
      size_t a;
      for(a=0;a<vPeaksRev[it->second].index.size();a++){
        if(vPeaksRev[it->second].index[a].pepMass==pepMass) break;
      }
      if(a== vPeaksRev[it->second].index.size()) {
        vPeaksRev[it->second].index.emplace_back();
        vPeaksRev[it->second].index[a].pepMass=pepMass;
      }
      vPeaksRev[it->second].index[a].pepIndex.push_back(pepIndex);
    } else vPeaksRev[it->second].index[0].pepIndex.push_back(pepIndex);
  }
}

double MIons::getAAMass(char aa){
  return aaMass[aa];
}

double MIons::getFixedModMass(char aa){
  return aaFixedModMass[aa];
}

//Iterates over the peptide amino acid sequence and creates a map of all variations
//that are possible (variable modifications, adducts, etc.).
void MIons::modIonsMaskRec(int pos, double mMass, int oSite, size_t pepNum, int depth, int modSite, string mask) {

  for (int i = pos; i < ionCount; i++) { //recursive function continues from last position

    //see if adduct attaches here
    if (oSite == -99 && i == 0 && nPep1 && site['n']) {  //if n-terminus can be linked, proceed as if it is linked
      pepCount++;
      string mask2(mask);
      mask2[pep1Len] = 'x';
      //cout << mask2 << endl;
      mPrecursor.insert(pair<string, size_t>(mask2, pepCount));
      pepMass.push_back(pepMass[pepNum]);
      pepMods.push_back(pepMods[pepNum]);
      modIonsMaskRec(i, mMass, i, pepCount, depth, modSite, mask2);

    } else if (oSite == -99 && modSite < i && site[pep1[i]]) { //otherwise attach adduct to available site
      pepCount++;
      string mask2(mask);
      mask2[i] = 'x';
      //cout << mask2 << endl;
      mPrecursor.insert(pair<string, size_t>(mask2, pepCount));
      pepMass.push_back(pepMass[pepNum]);
      pepMods.push_back(pepMods[pepNum]);
      modIonsMaskRec(i, mMass, i, pepCount, depth, modSite, mask2);
    }

    //only check modifications if adduct is not attached and we are allowed more modifications
    if (i > modSite && i != oSite && depth < maxModCount) {
      //Check if amino acid is on the modification list
      for (int j = 0; j < aaMod[pep1[i]].count; j++) {

        //skip mods if it is xl mod on a cut site - is this relevant in Magnum?
        if (aaMod[pep1[i]].mod[j].xl && i == pep1Len - 1 && !cPep1) continue;

        //Add masses
        pepCount++;
        string mask2(mask);
        mask2[i] = 49 + j;
        //cout << mask2 << endl;
        mPrecursor.insert(pair<string, size_t>(mask2, pepCount));
        pepMass.push_back(pepMass[pepNum] + aaMod[pep1[i]].mod[j].mass);
        pepMods.push_back(pepMods[pepNum]);
        mPepMod pm;
        pm.mass = aaMod[pep1[i]].mod[j].mass;
        pm.pos = (char)i;
        pm.term = false;
        pepMods.back().mods.push_back(pm);
        modIonsMaskRec(i, mMass + aaMod[pep1[i]].mod[j].mass, oSite, pepCount, depth + 1, i, mask2);

      }
    }

    mMass += aaMass[pep1[i]]; //add the amino acid mass
    if (i == ionCount) { //special case for the end of a peptide
      mMass += aaMass['c']; //add static mods to last amino acid
      pepMass[pepNum] += aaMass['c'];
      if (cPep1) {
        mMass += aaMass['%']; //add static mods to end of protein
        pepMass[pepNum] += aaMass['%'];
      }
    }
    addPeakNew(mMass, pepMass[pepNum], pepNum);

  }

  if (oSite == -99 && cPep1 && site['c']) { //if c-terminus can be linked, proceed as if it is linked
    pepCount++;
    string mask2(mask);
    mask2[pep1Len + 1] = 'x';
    //cout << mask2 << endl;
    mPrecursor.insert(pair<string, size_t>(mask2, pepCount));
    pepMass.push_back(pepMass[pepNum]);
    pepMods.push_back(pepMods[pepNum]);
    modIonsMaskRec(ionCount, mMass, ionCount, pepCount, depth, modSite, mask2);
  }

  //if the c-terminus can be modified, proceed with the modifications
  if (pos == ionCount && cPep1 && oSite != pos && depth < maxModCount) {
    //Check if amino acid is on the modification list
    for (int j = 0; j < aaMod[pep1[pos]].count; j++) {
      //Add masses
      pepCount++;
      string mask2(mask);
      mask2[pos] = 49 + j;
      //cout << mask2 << endl;
      mPrecursor.insert(pair<string, size_t>(mask2, pepCount));
      pepMass.push_back(pepMass[pepNum] + aaMod[pep1[pos]].mod[j].mass);
      pepMods.push_back(pepMods[pepNum]);
      mPepMod pm;
      pm.mass = aaMod[pep1[pos]].mod[j].mass;
      pm.pos = (char)pos;
      pm.term = false;
      pepMods.back().mods.push_back(pm);
      modIonsMaskRec(ionCount, mMass + aaMod[pep1[pos]].mod[j].mass, oSite, pepCount, depth + 1, pos, mask2);

    }
  }

  //check if mass boundaries have changed
  if (pepMass[pepNum] < pepMassMin)pepMassMin = pepMass[pepNum];
  if (pepMass[pepNum] > pepMassMax)pepMassMax = pepMass[pepNum];

  //Mark position of adduct
  while (pepLinks.size() <= pepNum) pepLinks.push_back(0);
  pepLinks[pepNum] = oSite;
  if (oSite > maxLink) maxLink = oSite;
}

void MIons::modIonsMaskRecNoAdduct(int pos, double mMass, size_t pepNum, int depth, int modSite, string mask) {

  for (int i = pos; i < ionCount - 1; i++) { //recursive function continues from last position

    //only check modifications if we are allowed more modifications
    if (i > modSite && depth < maxModCount) {
      //Check if amino acid is on the modification list
      for (int j = 0; j < aaMod[pep1[i]].count; j++) {

        //skip mods if it is xl mod on a cut site - is this relevant in Magnum?
        if (aaMod[pep1[i]].mod[j].xl && i == pep1Len - 1 && !cPep1) continue;

        //Add masses
        pepCount++;
        string mask2(mask);
        mask2[i] = 49 + j;
        //cout << mask2 << endl;
        mPrecursor.insert(pair<string, size_t>(mask2, pepCount));
        pepMass.push_back(pepMass[pepNum] + aaMod[pep1[i]].mod[j].mass);
        pepMods.push_back(pepMods[pepNum]);
        mPepMod pm;
        pm.mass = aaMod[pep1[i]].mod[j].mass;
        pm.pos = (char)i;
        pm.term = false;
        pepMods.back().mods.push_back(pm);
        modIonsMaskRecNoAdduct(i, mMass + aaMod[pep1[i]].mod[j].mass, pepCount, depth + 1, i, mask2);

      }
    }

    mMass += aaMass[pep1[i]]; //add the amino acid mass
    if (i == ionCount - 1) { //special case for the end of a peptide
      mMass += aaMass['c']; //add static mods to last amino acid
      pepMass[pepNum] += aaMass['c'];
      if (cPep1) {
        mMass += aaMass['%']; //add static mods to end of protein
        pepMass[pepNum] += aaMass['%'];
      }
    }
  }

  //check if mass boundaries have changed
  if (pepMass[pepNum] < pepMassMin)pepMassMin = pepMass[pepNum];
  if (pepMass[pepNum] > pepMassMax)pepMassMax = pepMass[pepNum];
}

double MIons::getModMass(int index){
  return modMassArray[index];
}

int MIons::getModMassSize(){
  return (int)modMassArray.size();
}

double* MIons::getMods(){
  return modList;
}

int MIons::getIonCount (){
  return ionCount;
}

void MIons::getPeptide(char *seq){
  strncpy(seq,pep1,pep1Len);
  seq[pep1Len]='\0';
}

int MIons::getPeptideLen(){
  return pep1Len;
}

void MIons::getPeptideMods(vector<mPepMod>& v){
  unsigned int i;
  mPepMod m;
  v.clear();
  for(i=0;i<modQueue.size();i++){
    if(modQueue[i].pos==-1) {
      m.pos=0; //convert n-terminus to first AA
      m.mass=modList[m.pos];
    } else {
      m.pos=(char)modQueue[i].pos;
      m.mass=modList[m.pos+1];
    }
    v.push_back(m);
  }
}

void MIons::setAAMass(char aa, double mass){
  aaMass[aa]=mass;
}

void MIons::setMaxModCount(int i){
  maxModCount=i;
}

void MIons::setPeptide(char* seq, int len, double mass, bool nTerm, bool cTerm){
  pep1=seq;
  pep1Len=len;
  pep1Mass=mass;
  nPep1=nTerm;
  cPep1=cTerm;

  char pp[256];
  strncpy(pp, seq, len);
  pp[len] = '\0';
  pepseq = pp;

}

/*============================
  Utilities
============================*/
int MIons::compareD(const void *p1, const void *p2){
  const double d1 = *(double *)p1;
  const double d2 = *(double *)p2;
  if(d1<d2) return -1;
  else if(d1>d2) return 1;
  else return 0;
}

