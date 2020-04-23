#ifndef _COMETDECOYS_H
#define _COMETDECOYS_H

#define MAX_DECOY_PEP_LEN 40
#define DECOY_SIZE        3000

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

