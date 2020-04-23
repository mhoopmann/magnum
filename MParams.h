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

#ifndef _MPARAMS_H
#define _MPARAMS_H

#include "MStructs.h"
#include "pepXMLWriter.h"

#ifdef _MSC_VER
#include <direct.h>
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

  bool buildOutput(char* in, char* base, char* ext);
  bool parseConfig(char* fname);
 
private:

  mParams* params;
  
  bool checkMod(mMass m);
  void parse(char* cmd);
  bool processPath(const char* cwd, const char* in_path, char* out_path);
  void warn(const char* c, int i);

};

#endif
