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

#ifndef _MIONS_H
#define _MIONS_H

#include "MStructs.h"
#include "MIonSet.h"
#include <vector>

typedef struct mModPos{
  int     pos;        //position of aa to add
  int     mod;        //index of mods at this position
  int     modCount;   //total mods prior to this position
  double  mass;       //total mass prior to this position
  double  modMass;    //total mod mass prior to this position
} mModPos;

typedef struct mModType{
  bool    xl;
  double  mass;
} mModType;

//max 10 mods per amino acid
typedef struct mMod{
  int       count;
  mModType  mod[10];
} mMod;

class MIons {
public:
  MIons();
  ~MIons();

  //Functions
  void      addFixedMod       (char mod, double mass);
  void      addMod            (char mod, bool xl, double mass);
  void      buildIons         ();
  void      buildModIons      (int modSite);
  double    getAAMass         (char aa);
  double    getFixedModMass   (char aa);
  void      modIonsRec        (int start, int link, int index, int depth, bool xl);
  void      modIonsRec2       (int start, int link, int index, int depth, bool xl);
  void      reset             ();

  //Accessors
  MIonSet&  operator[ ]   (const int& i);
  MIonSet*  at            (const int& i);
  int       getIonCount   ();
  double    getModMass    (int index);
  int       getModMassSize();
  double*   getMods       ();
  void      getPeptide    (char* seq);
  int       getPeptideLen ();
  void      getPeptideMods(std::vector<mPepMod>& v);
  int       size          ();

  //Modifiers
  void  setAAMass       (char aa, double mass);
  void  setMaxModCount  (int i);
  void  setPeptide      (char* seq, int len, double mass, bool nTerm, bool cTerm);

  //Data Members
  double* modList;
  bool site[128]; //possible sites of linkage based on parameters

private:

  void addModIonSet(int index, char aa, int pos, int modIndex, int loopPos=-1);
  void buildSeries(int setNum);
  void clearSeries();
  
  double  aaMass[128];
  double  aaFixedModMass[128];
  mMod    aaMod[128];   //inefficient memory usage, but not by much in the grand scheme.
  double  protFixedModMassC;
  double  protFixedModMassN;
  mMod  protModMassC;
  mMod  protModMassN;

  double  pep1Mass;
  int     modIndex;

  int ionCount;
  int maxModCount;
  int pep1Len;

  char* pep1;

  bool nPep1; //peptide has protein n-terminus
  bool cPep1; //peptide has protein c-terminus

  std::vector<mModPos> modQueue;
  std::vector<double>  modMassArray;
  std::vector<MIonSet> sets;

  //Utilities
  static int compareD(const void *p1,const void *p2);

};


#endif
