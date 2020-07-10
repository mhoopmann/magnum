#ifndef _COMETDECOYS_H
#define _COMETDECOYS_H

#define MAX_DECOY_PEP_LEN 70
#define DECOY_SIZE        5000

typedef struct DecoysStruct
{
   char *szPeptide;
   double pdIonsN[MAX_DECOY_PEP_LEN];
   double pdIonsC[MAX_DECOY_PEP_LEN];
} DecoysStruct;

class MDecoys{
public:
  static DecoysStruct decoyIons[];
};

#endif

