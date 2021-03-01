/*
Copyright 2021, Michael R. Hoopmann, Institute for Systems Biology

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

#ifndef _MPARAMS_H
#define _MPARAMS_H

#include "MLog.h"
#include "MStructs.h"
#include "pepXMLWriter.h"

#ifdef _MSC_VER
#include <direct.h>
#include <Windows.h>
#define getcwd _getcwd
#define slashdir '\\'
#else
#include <unistd.h>
#define slashdir '/'
#endif

class MParams {
public:
  MParams();
  MParams(mParams* p);
  ~MParams();

  std::vector<pxwBasicXMLTag> xmlParams;

  bool buildOutput(std::string in, std::string base, std::string ext);
  void parse(const char* cmd);
  bool parseConfig(const char* fname);
  void setLog(MLog* c);
  void setParams(mParams* p);

  std::string logFile;
 
private:

  MLog* log;
  mParams* params;
  
  bool checkMod(mMass m);
  void logParam(std::string name, std::string value);
  void logParam(pxwBasicXMLTag& t);
  bool processPath(const char* cwd, const char* in_path, char* out_path);
  void trimPath(std::string in, std::string& path, std::string& fname);
  void warn(const char* c, int i);
  void warn(std::string c, int i);

};

#endif
