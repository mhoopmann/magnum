#include "MTopPeps.h"

MTopPeps::MTopPeps(){
  peptideCount = 0;
  peptideFirst = NULL;
  peptideLast = NULL;
  peptideMax = 0;

  //singletList = NULL;
  //singletBins = 0;
}

MTopPeps::MTopPeps(const MTopPeps& c){
  peptideCount = c.peptideCount;
  peptideMax = c.peptideMax;
  peptideFirst = NULL;
  peptideLast = NULL;
  mScoreCard* sc = NULL;
  mScoreCard* tmp = c.peptideFirst;
  if (tmp != NULL) {
    peptideFirst = new mScoreCard(*tmp);
    sc = peptideFirst;
    tmp = tmp->next;
    while (tmp != NULL){
      sc->next = new mScoreCard(*tmp);
      sc->next->prev = sc;
      sc = sc->next;
      tmp = tmp->next;
    }
    peptideLast = sc;
  }

  /*
  singletBins = c.singletBins;
  if (singletBins == 0) singletList = NULL;
  else {
    singletList = new list<mSingletScoreCard*>*[singletBins];
    for (size_t j = 0; j<singletBins; j++){
      if (c.singletList[j] == NULL) singletList[j] = NULL;
      else {
        singletList[j] = new list<mSingletScoreCard*>;
        list<mSingletScoreCard*>::iterator it = c.singletList[j]->begin();
        while (it != c.singletList[j]->end()){
          singletList[j]->emplace_back(*it);
          it++;
        }
      }
    }
  }
  */
}

MTopPeps::~MTopPeps(){
  while (peptideFirst != NULL){
    mScoreCard* tmp = peptideFirst;
    peptideFirst = peptideFirst->next;
    delete tmp;
  }
  peptideLast = NULL;
}

MTopPeps& MTopPeps::operator=(const MTopPeps& c){
  if(this!=&c){
    peptideCount = c.peptideCount;
    peptideMax = c.peptideMax;
    
    while (peptideFirst != NULL){
      mScoreCard* tmp = peptideFirst;
      peptideFirst = peptideFirst->next;
      delete tmp;
    }
    peptideLast = NULL;

    mScoreCard* sc = NULL;
    mScoreCard* tmp = c.peptideFirst;
    if (tmp != NULL) {
      peptideFirst = new mScoreCard(*tmp);
      sc = peptideFirst;
      tmp = tmp->next;
      while (tmp != NULL){
        sc->next = new mScoreCard(*tmp);
        sc->next->prev = sc;
        sc = sc->next;
        tmp = tmp->next;
      }
      peptideLast = sc;
    }
  }
  return *this;
}

void MTopPeps::checkPeptideScore(mScoreCard& s){

  mScoreCard* sc;
  mScoreCard* cur;

  //If list is empty, add the score card
  if (peptideCount == 0){
    peptideFirst = new mScoreCard(s);
    peptideLast = peptideFirst;
    peptideCount++;
    return;
  }

  //check if we can just add to the end
  if (s.simpleScore<peptideLast->simpleScore){
    //check if we need to store the singlet
    if (peptideCount == peptideMax) return;

    peptideLast->next = new mScoreCard(s);
    peptideLast->next->prev = peptideLast;
    peptideLast = peptideLast->next;
    peptideCount++;
    return;
  }

  //check if it goes in the front
  if (s.simpleScore >= peptideFirst->simpleScore){
    peptideFirst->prev = new mScoreCard(s);
    peptideFirst->prev->next = peptideFirst;
    peptideFirst = peptideFirst->prev;

    if (peptideCount<peptideMax) {
      peptideCount++;
    } else {
      cur = peptideLast;
      peptideLast = peptideLast->prev;
      peptideLast->next = NULL;
      delete cur;
    }
    return;
  }


  //scan to find insertion point
  cur = peptideFirst->next;
  int i = 1;
  while (s.simpleScore < cur->simpleScore){
    i++;
    cur = cur->next;
  }

  sc = new mScoreCard(s);
  sc->prev = cur->prev;
  sc->next = cur;
  cur->prev->next = sc;
  cur->prev = sc;
  if (sc->prev == NULL) peptideFirst = sc;

  if (peptideCount<peptideMax) {
    peptideCount++;
  } else {
    cur = peptideLast;
    peptideLast = peptideLast->prev;
    peptideLast->next = NULL;
    delete cur;
  }

}

//void MTopPeps::resetSingletList(double mass){
  /*
  size_t j;
  if (singletList != NULL){
    for (j = 0; j<singletBins; j++){
      if (singletList[j] != NULL) delete singletList[j];
    }
    delete[] singletList;
  }
  singletBins = (int)(mass / 10 + 1);
  singletList = new list<mSingletScoreCard*>*[singletBins];
  for (j = 0; j<singletBins; j++) singletList[j] = NULL;
  */
//}
