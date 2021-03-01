#include "MLog.h"

using namespace std;

MLog::MLog(){
  clear();
}

void MLog::addDBWarning(std::string msg){
  string myMsg = "WARNING: " + msg;
  cout << myMsg << endl;
  vDBWarnings.push_back(myMsg);
}

void MLog::addError(std::string msg){
  strError = "\n";
  time_t timeNow;
  time(&timeNow);
  strError += ctime(&timeNow);
  strError.pop_back();
  strError += "\tERROR: " + msg;
  cout << strError << endl;
  exportLog();
  exit(-4004);
}

void MLog::addMessage(string msg, bool silent){
  if (!silent) cout << msg << endl;

  time_t timeNow;
  time(&timeNow);
  string myMsg = ctime(&timeNow);
  myMsg.pop_back();
  myMsg += "\t" + msg;

  vMsg.push_back(myMsg);
}

void MLog::addParameter(std::string msg){
  vParams.push_back(msg);
}

void MLog::addParameterWarning(std::string msg){
  string myMsg = "WARNING: " + msg;
  cout << myMsg << endl;
  vParamWarnings.push_back(myMsg);
}

void MLog::addWarning(size_t id, std::string msg){
  if (idIndex[id] != SIZE_MAX){
    vWarnings[idIndex[id]].count++;
  } else {
    string myMsg = "WARNING: " + msg;
    mWarning w;
    w.count = 1;
    w.msg = myMsg;
    idIndex[id] = vWarnings.size();
    vWarnings.push_back(w);
    cout << myMsg << endl;
  }
}

void MLog::clear(){
  for (size_t i = 0; i<LOGSZ; i++) idIndex[i] = SIZE_MAX;
  logFile.clear();
  vMsg.clear();
  vWarnings.clear();
  strError.clear();
}

void MLog::exportLog(){
  if (logFile.size() == 0) return;
  size_t i;
  FILE* f = fopen(logFile.c_str(), "wt");

  fprintf(f, "***** PARAMETERS *****\n");
  for (i = 0; i<vParams.size(); i++) fprintf(f, "%s\n", vParams[i].c_str());
  fprintf(f, "\n");
  for (i = 0; i<vParamWarnings.size(); i++) fprintf(f, "%s\n", vParamWarnings[i].c_str());

  fprintf(f, "\n\n***** DATABASE *****\n");
  fprintf(f, "%s", dbInfo.c_str());
  fprintf(f, "\n");
  for (i = 0; i<vDBWarnings.size(); i++) fprintf(f, "%s\n", vDBWarnings[i].c_str());

  fprintf(f, "\n\n***** MAGNUM LOG *****\n");
  for (i = 0; i<vMsg.size(); i++) fprintf(f, "%s\n", vMsg[i].c_str());

  if (strError.size()>0) {
    fprintf(f, "\n\n***** MAGNUM EXECUTION HALTED DUE TO ERROR!!! *****");
    fprintf(f, "%s\n", strError.c_str());
  }

  if (vWarnings.size()>0 || vParamWarnings.size()>0){
    fprintf(f, "\n\n***** WARNINGS *****\n");
    for (i = 0; i<vWarnings.size(); i++) fprintf(f, "(%d instances) %s\n", vWarnings[i].count, vWarnings[i].msg.c_str());
    if (vParamWarnings.size()>0) fprintf(f, "Check Parameter log above for warnings originating in configuration file.\n");
    if (vDBWarnings.size()>0) fprintf(f, "Check Database log above for warnings originating in the FASTA database file.\n");
  }

  fclose(f);
}

void MLog::setDBinfo(std::string fn, int prot, int pep, int adduct){
  dbInfo = "FASTA database: ";
  dbInfo += fn;
  char tStr[512];
  sprintf(tStr, "\nTotal proteins: %d\nTotal peptides to search: %d (%d with binding sites)", prot, pep, adduct);
  dbInfo += tStr;
}

void MLog::setLog(char* fn){
  string f = fn;
  setLog(f);
}

void MLog::setLog(std::string fn){
  logFile = fn;
}