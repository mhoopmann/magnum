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

#ifndef _MTOPPEPS_H
#define _MTOPPEPS_H

#include "MStructs.h"
#include <list>

class MTopPeps{
public:

  MTopPeps();
  MTopPeps(const MTopPeps& c);
  ~MTopPeps();

  MTopPeps& operator=(const MTopPeps& c);

  //int singletBins;
  int peptideCount;
  int peptideMax;
  
  mScoreCard*    peptideFirst;   //pointer to start of linked list
  mScoreCard*    peptideLast;    //pointer to end of linked list
  //list<mSingletScoreCard*>**  singletList;

  void  checkPeptideScore(mScoreCard& s);
  //void  resetSingletList(double mass);

};

#endif
