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
  peaks=NULL;
}

MIons::~MIons(){
  if(modList!=NULL) delete [] modList;
  if(peaks!=NULL) delete peaks;
  pep1=NULL;
}

MIonSet& MIons::operator[ ](const int& i){
  return sets[i];
}

MIonSet* MIons::at(const int& i){
  return &sets[i];
}


void MIons::addFixedMod(char mod, double mass){
  aaMass[mod]+=mass;
  aaFixedModMass[mod]=mass;
}

void MIons::addMod(char mod, bool xl, double mass){
  aaMod[mod].mod[aaMod[mod].count].xl=xl;
  aaMod[mod].mod[aaMod[mod].count++].mass=mass;
}

void MIons::addPeak(double mass, bool adduct, size_t& node, bool& ar, size_t& index){
  double pMass = mass;
  if (adduct) pMass = -pMass;
  size_t target;
  map<double, size_t>::iterator it = mP.find(pMass);
  if (it == mP.end()) {
    target = peaks->size();
    mP.insert(pair<double, size_t>(pMass, target));
    if (ar) {
      peaks->at(node).start[index].nextNode = target;
      peaks->at(node).start[index].nextIndex = 0;
    } else {
      peaks->at(node).next[index].nextNode = target;
      peaks->at(node).next[index].nextIndex = 0;
    }
    peaks->emplace_back();
    peaks->back().mass = pMass;
    index = peaks->back().next.size();
    peaks->back().next.emplace_back();

  } else {
    target = it->second;
    if (ar) {
      peaks->at(node).start[index].nextNode = target;
      peaks->at(node).start[index].nextIndex = peaks->at(target).next.size();
    } else {
      peaks->at(node).next[index].nextNode = target;
      peaks->at(node).next[index].nextIndex = peaks->at(target).next.size();
    }
    index = peaks->at(target).next.size();
    peaks->at(target).next.emplace_back();

  }
  node = target;
  ar = false;
}

void MIons::buildIons(){
	
  int b;
	int i;
  int j;
  int y;
  double fragMass;

  //set up boundaries
  ionCount=sets[0].len;
  b=0;
  y=ionCount-1;

	//b- & y-ions from first peptide
  fragMass=aaMass['n'];
  if (nPep1) fragMass+=aaMass['$'];
	for(i=0;i<ionCount;i++){
    fragMass+=aaMass[pep1[i]];
    if (i == ionCount - 1) {
      fragMass += aaMass['c'];
      if (cPep1) fragMass += aaMass['%'];
    }
    sets[0].aIons[0][b] = fragMass - 27.9949141;
    sets[0].bIons[0][b] = fragMass;
    sets[0].cIons[0][b] = fragMass + 17.026547;
    sets[0].xIons[0][y] = (pep1Mass-fragMass) + 25.9792649;
    sets[0].yIons[0][y] = (pep1Mass-fragMass);
    sets[0].zIons[0][y] = (pep1Mass-fragMass) - 16.0187224; //z. not z
    b++;
    y--;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<4;j++){
      sets[0].aIons[j][i] = ((sets[0].aIons[0][i]+1.007276466*j)/j);
      sets[0].bIons[j][i] = ((sets[0].bIons[0][i]+1.007276466*j)/j);
      sets[0].cIons[j][i] = ((sets[0].cIons[0][i]+1.007276466*j)/j);
      sets[0].xIons[j][i] = ((sets[0].xIons[0][i]+1.007276466*j)/j);
      sets[0].yIons[j][i] = ((sets[0].yIons[0][i]+1.007276466*j)/j);
      sets[0].zIons[j][i] = ((sets[0].zIons[0][i]+1.007276466*j)/j);
    }
  }

}

void MIons::buildModIons(int modSite){

  int b;
	int i;
  int j;
  int y;
	double mMass;

  //set up boundaries
  ionCount=sets[0].len; //pep1Len-1;
  b=0;
  y=ionCount-1;

	//b- & y-ions from first peptide
  mMass=aaMass['n'];
  if (nPep1) mMass += aaMass['$'];
	for(i=0;i<ionCount;i++){
    mMass+=aaMass[pep1[i]];
    if (i == ionCount - 1) {
      mMass += aaMass['n'];
      if (cPep1) mMass += aaMass['%'];
    }
    if(i>=modSite) { //negative mass indicates that the ion needs a modmass as well
      sets[0].xIons[0][y] = pep1Mass-mMass + 25.9792649;
      sets[0].yIons[0][y] = pep1Mass-mMass;
      sets[0].zIons[0][y] = pep1Mass-mMass - 16.0187224;   //z.
      sets[0].aIons[0][b] = -(mMass - 27.9949141);
      sets[0].bIons[0][b] = -mMass;
      sets[0].cIons[0][b] = -(mMass + 17.026547);
    } else {
      sets[0].xIons[0][y] = -(pep1Mass-mMass + 25.9792649);
      sets[0].yIons[0][y] = -(pep1Mass-mMass);
      sets[0].zIons[0][y] = -(pep1Mass-mMass - 16.0187224);
      sets[0].aIons[0][b] = mMass - 27.9949141;
      sets[0].bIons[0][b] = mMass;
      sets[0].cIons[0][b] = mMass + 17.026547;
    }
    y--;
    b++;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<4;j++){
      if(sets[0].aIons[0][i]<0) sets[0].aIons[j][i] = ((sets[0].aIons[0][i]-1.007276466*j)/j);
      else sets[0].aIons[j][i] = ((sets[0].aIons[0][i]+1.007276466*j)/j);
      if(sets[0].bIons[0][i]<0) sets[0].bIons[j][i] = ((sets[0].bIons[0][i]-1.007276466*j)/j);
      else sets[0].bIons[j][i] = ((sets[0].bIons[0][i]+1.007276466*j)/j);
      if(sets[0].cIons[0][i]<0) sets[0].cIons[j][i] = ((sets[0].cIons[0][i]-1.007276466*j)/j);
      else sets[0].cIons[j][i] = ((sets[0].cIons[0][i]+1.007276466*j)/j);
      if(sets[0].xIons[0][i]<0) sets[0].xIons[j][i] = ((sets[0].xIons[0][i]-1.007276466*j)/j);
      else sets[0].xIons[j][i] = ((sets[0].xIons[0][i]+1.007276466*j)/j);
      if(sets[0].yIons[0][i]<0) sets[0].yIons[j][i] = ((sets[0].yIons[0][i]-1.007276466*j)/j);
      else sets[0].yIons[j][i] = ((sets[0].yIons[0][i]+1.007276466*j)/j);
      if(sets[0].zIons[0][i]<0) sets[0].zIons[j][i] = ((sets[0].zIons[0][i]-1.007276466*j)/j);
      else sets[0].zIons[j][i] = ((sets[0].zIons[0][i]+1.007276466*j)/j);
    }
  }

}

void MIons::buildModIons2() {
  double mMass;

  if(peaks!=NULL) delete peaks;
  peaks = new vector<sNode2>;
  mP.clear();
  pepLinks.clear();
  pepMass.clear();
  pepMods.clear();
  maxLink = -99;

  peaks->emplace_back();
  peaks->at(0).start.emplace_back();
  peaks->at(0).start.back().pepNum = 0;

  //set up boundaries
  ionCount = pep1Len;

  //peptide number
  pepCount = 0;

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
  modIonsRec4(0, mMass, -99, pepCount, 0, -1, 0, true, 0);
  pepCount++; //Double-check this?

}

double MIons::getAAMass(char aa){
  return aaMass[aa];
}

double MIons::getFixedModMass(char aa){
  return aaFixedModMass[aa];
}

void MIons::addModIonSet(int index, char aa, int pos, int modIndex, int loopPos){
  int k,n;

  //Add masses
  MIonSet s = sets[index];

  for (k = pos; k<ionCount; k++){
    for (n = 1; n<4; n++){
      if (s.aIons[0][k]<0) s.aIons[n][k] -= (aaMod[aa].mod[modIndex].mass / n);
      else s.aIons[n][k] += (aaMod[aa].mod[modIndex].mass / n);
      if (s.bIons[0][k]<0) s.bIons[n][k] -= (aaMod[aa].mod[modIndex].mass / n);
      else s.bIons[n][k] += (aaMod[aa].mod[modIndex].mass / n);
      if (s.cIons[0][k]<0) s.cIons[n][k] -= (aaMod[aa].mod[modIndex].mass / n);
      else s.cIons[n][k] += (aaMod[aa].mod[modIndex].mass / n);
    }
  }
  for (k = ionCount - pos; k<ionCount; k++){
    for (n = 1; n<4; n++){
      if (s.xIons[0][k]<0) s.xIons[n][k] -= (aaMod[aa].mod[modIndex].mass / n);
      else s.xIons[n][k] += (aaMod[aa].mod[modIndex].mass / n);
      if (s.yIons[0][k]<0) s.yIons[n][k] -= (aaMod[aa].mod[modIndex].mass / n);
      else s.yIons[n][k] += (aaMod[aa].mod[modIndex].mass / n);
      if (s.zIons[0][k]<0) s.zIons[n][k] -= (aaMod[aa].mod[modIndex].mass / n);
      else s.zIons[n][k] += (aaMod[aa].mod[modIndex].mass / n);
    }
  }
  if (aa == 'n' || aa=='$') s.modNTerm=true;
  if (aa == 'c' || aa=='%') s.modCTerm=true;
  if (loopPos>-1) s.mods[loopPos] = aaMod[aa].mod[modIndex].mass;
  else s.mods[pos] = aaMod[aa].mod[modIndex].mass;
  s.mass += aaMod[aa].mod[modIndex].mass;
  s.difMass += aaMod[aa].mod[modIndex].mass;

  //Add to list
  sets.push_back(s);

}

void MIons::modIonsRec(int start, int link, int index, int depth, bool xl){
  int i,j;

  for(i=start;i<pep1Len;i++){

    //don't modify site where the cross-linker is bound.
    if (i == link) continue;

    //Check if amino acid is on the modification list
    for(j=0;j<aaMod[pep1[i]].count;j++){

      //skip mods if it is xl mod on a cut site
      if (aaMod[pep1[i]].mod[j].xl && i == pep1Len - 1 && !cPep1) continue;

      //Add masses
      addModIonSet(index,pep1[i],i,j);

      //solve another one
      if(depth+1<maxModCount) modIonsRec(i+1,link,(int)(sets.size())-1,depth+1,xl);
    }
  
  }

}

void MIons::modIonsRec2(int start, int link, int index, int depth, bool xl){
  int j;

  //if n-terminus can be linked, proceed as if it is linked
  if (link == 0 && nPep1 && site['n']){
    modIonsRec(start, -1, index, depth, xl);
  } else if (link == pep1Len - 1 && cPep1 && site['c']){ //if c-terminus can be linked, proceed as if it is linked
    modIonsRec(start, -1, index, depth, xl);
  }

  //now proceed with linking on an amino acid
  modIonsRec(start, link, index, depth, xl);

  //if at n-terminus and aa is not linkable, stop now; n-terminus cannot be modified because it is holding the linker
  if(link==0 && !site[pep1[0]]) return; 
  else if(link==pep1Len-1 && !site[pep1[pep1Len-1]]) return; //ditto for c-terminus

  //From here, proceed as if the link is on the amino acid with a terminal modification
  //Special case for protein n-terminus
  if (nPep1) {
    for (j = 0; j<aaMod['$'].count; j++){
      //Add masses
      addModIonSet(index, '$', 0, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  } else {    //Check peptide n-terminus mods
    for (j = 0; j<aaMod['n'].count; j++){
      //Add masses
      addModIonSet(index, 'n', 0, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  }

  //Special case for protein c-terminus
  if (cPep1){
    for (j = 0; j<aaMod['%'].count; j++){
      //Add masses
      addModIonSet(index, '%', pep1Len-1, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  } else { //Special case for peptide c-terminus
    for (j = 0; j<aaMod['c'].count; j++){
      //Add masses
      addModIonSet(index, 'c', pep1Len-1, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  }

}

void MIons::modIonsRec4(int pos, double mMass, int oSite, size_t pepNum, int depth, int modSite, size_t linkNode, bool linkAr, size_t linkIndex) {

  for (int i = pos; i < ionCount - 1; i++) { //recursive function continues from last position

    //see if adduct attaches here
    if (oSite == -99 && i == 0 && nPep1 && site['n']) {  //if n-terminus can be linked, proceed as if it is linked
      pepCount++;
      pepMass.push_back(pepMass[pepNum]);
      pepMods.push_back(pepMods[pepNum]);
      peaks->at(linkNode).start.emplace_back();
      peaks->at(linkNode).start.back().pepNum = pepCount;
      modIonsRec4(i, mMass, i, pepCount, depth, modSite, linkNode, true, peaks->at(linkNode).start.size() - 1);

    } else if (oSite == -99 && modSite < i && site[pep1[i]]) { //otherwise attach adduct to available site
      pepCount++;
      pepMass.push_back(pepMass[pepNum]);
      pepMods.push_back(pepMods[pepNum]);
      peaks->at(linkNode).start.emplace_back();
      peaks->at(linkNode).start.back().pepNum = pepCount;
      modIonsRec4(i, mMass, i, pepCount, depth, modSite, linkNode, true, peaks->at(linkNode).start.size() - 1);
    }

    //only check modifications if adduct is not attached and we are allowed more modifications
    if (i > modSite && i != oSite && depth < maxModCount) {
      //Check if amino acid is on the modification list
      for (int j = 0; j < aaMod[pep1[i]].count; j++) {

        //skip mods if it is xl mod on a cut site - is this relevant in Magnum?
        if (aaMod[pep1[i]].mod[j].xl && i == pep1Len - 1 && !cPep1) continue;

        //Add masses
        pepCount++;
        pepMass.push_back(pepMass[pepNum] + aaMod[pep1[i]].mod[j].mass);
        pepMods.push_back(pepMods[pepNum]);
        mPepMod pm;
        pm.mass = aaMod[pep1[i]].mod[j].mass;
        pm.pos = (char)i;
        pm.term = false;
        pepMods.back().mods.push_back(pm);
        peaks->at(linkNode).start.emplace_back();
        peaks->at(linkNode).start.back().pepNum = pepCount;
        modIonsRec4(i, mMass + aaMod[pep1[i]].mod[j].mass, oSite, pepCount, depth, i, linkNode, true, peaks->at(linkNode).start.size() - 1);

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

    if (oSite >= 0 && i >= oSite) addPeak(mMass, true, linkNode, linkAr, linkIndex);
    else addPeak(mMass, false, linkNode, linkAr, linkIndex);
    peaks->at(linkNode).next[linkIndex].pepNum = pepNum;
  }

  if (oSite == -99 && cPep1 && site['c']) { //if c-terminus can be linked, proceed as if it is linked
    pepCount++;
    pepMass.push_back(pepMass[pepNum]);
    pepMods.push_back(pepMods[pepNum]);
    peaks->at(linkNode).start.emplace_back();
    peaks->at(linkNode).start.back().pepNum = pepCount;
    modIonsRec4(ionCount, mMass, ionCount, pepCount, depth, modSite, linkNode, true, peaks->at(linkNode).start.size() - 1);

  }

  //check if mass boundaries have changed
  if (pepMass[pepNum] < pepMassMin)pepMassMin = pepMass[pepNum]; 
  if (pepMass[pepNum] > pepMassMax)pepMassMax = pepMass[pepNum];

  //Mark position of adduct
  while (pepLinks.size() <= pepNum) pepLinks.push_back(0);
  pepLinks[pepNum] = oSite;
  if (oSite>maxLink) maxLink = oSite;
}

void MIons::reset(){
  sets.clear();
  MIonSet k(pep1Len,pep1Mass);
  sets.push_back(k);
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

  sets.clear();
  MIonSet k(pep1Len,pep1Mass);
  sets.push_back(k);

}

int MIons::size(){
  return (int)sets.size();
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

