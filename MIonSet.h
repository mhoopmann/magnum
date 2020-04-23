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

#ifndef _MIONSET_H
#define _MIONSET_H

class MIonSet{
public:
  MIonSet(int sz, double m);
  MIonSet(const MIonSet& k);
  ~MIonSet();

  MIonSet& operator=(const MIonSet& k);

  double** aIons;
  double** bIons;
  double** cIons;
  double** xIons;
  double** yIons;
  double** zIons;
  double* mods;
  double  mass;
  double  difMass;
  int     len;
  bool modNTerm;
  bool modCTerm;
};

#endif