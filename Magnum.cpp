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

#include "MAnalysis.h"
#include "MData.h"
#include "MDB.h"
#include "MIons.h"
#include "MParams.h"
#include "MagnumManager.h"

int main(int argc, char* argv[]){

  cout << "\nMagnum version " << VERSION << ", " << BDATE << endl;
  cout << "Copyright Michael Hoopmann, Institute for Systems Biology" << endl;
  cout << "Visit http://magnum-ms.org for full documentation." << endl;
  if(argc<2){
    cout << "Usage: Magnum <Config File> [<Data File>...]" << endl;
    cout << "\nNote: To create a default configuration file for Magnum, run the following command:" << endl;
    cout << "        Magnum --config" << endl;
    return 1;
  }

  if (strcmp(argv[1],"--config") == 0) {
    MParams p;
    p.exportDefault(VERSION);
    cout << "\nmagnum_default_params.conf file created." << endl;
    return 2;
  }

  cout << "\n****** Begin Magnum Analysis ******" << endl;

  time_t timeNow;
  time(&timeNow);
  cout << " Time at start: " << ctime(&timeNow) << endl;

  int i;
  int fc=1;
  MagnumManager manager;
  if(!manager.setParams(argv[1])) return -3;
  if(argc>2){
    manager.clearFiles();
    for (i = 2; i<argc; i++) fc = manager.setFile(argv[i]);
  }
  manager.run();

  time(&timeNow);
  cout << " Time at finish: " << ctime(&timeNow) << endl;
  cout << "\n****** Finished Magnum Analysis ******" << endl;
  return 0;
}

