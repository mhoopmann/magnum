#ifndef _MAGNUMMANAGER_H
#define _MAGNUMMANAGER_H

#include "MLog.h"
#include "MParams.h"

#define VERSION "1.2.0"
#define BDATE "September 14 2022"

class MagnumManager {
public:
  MagnumManager();

  void clearFiles();

  int setFile(const char* fn);
  int setFile(std::string& s);
  void setParam(const char* p);
  void setParam(std::string& s);
  bool setParams(const char* fn);
  bool setParams(std::string& s);

  bool getBaseFileName(std::string& base, std::string& fName, std::string& extP);
  int run();


private:
  std::vector<mFile> files;
  MLog log;
  std::string paramFile;
  MParams param_obj;
  mParams params;

};

#endif
