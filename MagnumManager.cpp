#include "MagnumManager.h"
#include "MAnalysis.h"
#include "MData.h"
#include "MDB.h"
#include "MIons.h"

using namespace std;

MagnumManager::MagnumManager(){
  param_obj.setParams(&params);
  param_obj.setLog(&log);
}

void MagnumManager::clearFiles(){
  files.clear();
}

int MagnumManager::setFile(const char* fn){
  kFile f;
  f.input = fn;
  if (!getBaseFileName(f.base, f.input, f.ext)){
    cout << "  Error with batch file parameter:  " << fn << " - unknown file or extension." << endl;
    return -1;
  }
  files.push_back(f);
  return (int)files.size();
}

int MagnumManager::setFile(string& s){
  return setFile(s.c_str());
}

void MagnumManager::setParam(const char* p){
  param_obj.parse(p);
}

void MagnumManager::setParam(string& s){
  setParam(s.c_str());
}

bool MagnumManager::setParams(const char* fn){
  kFile f;
  paramFile = fn; //hang onto this

  cout << " Parameter file: " << fn << endl;
  param_obj.parseConfig(fn);
  f.input = params.msFile;
  if (!getBaseFileName(f.base, params.msFile, f.ext)){
    cout << "  Error with input file parameter: " << params.msFile << " - unknown file or extension." << endl;
    return false;
  }
  files.push_back(f);
  return true;

}

bool MagnumManager::setParams(string& s){
  return setParams(s.c_str());
}

int MagnumManager::run(){
  time_t timeNow;
  size_t i;

  //Step #1: Prepare from settings
  MData spec(&params);
  spec.setLog(&log);
  spec.setVersion(VERSION);
  spec.setAdductSites(params.adductSites);
  spec.setParams(&param_obj);

  //Step #2: Read in database and generate peptide lists
  MDatabase db;
  db.setLog(&log);
  for (i = 0; i<params.fMods.size(); i++) db.addFixedMod(params.fMods[i].index, params.fMods[i].mass);
  for (i = 0; i<params.aaMass.size(); i++) db.setAAMass((char)params.aaMass[i].index, params.aaMass[i].mass);
  if (!db.setEnzyme(params.enzyme.c_str())) exit(-3);
  db.setAdductSites(spec.getAdductSites());
  cout << "\n Reading FASTA database: " << params.dbFile << endl;
  if (!db.buildDB(params.dbFile.c_str(),params.decoy)){
    cout << "  Error opening database file: " << params.dbFile << endl;
    return -1;
  }
  if(params.buildDecoy) db.buildDecoy(params.decoy);
  db.buildPeptides(params.minPepMass, params.maxPepMass, params.miscleave, params.minPepLen, params.maxPepLen);
  log.setDBinfo(string(params.dbFile),db.getProteinDBSize(),db.getPeptideListSize(),db.adductPepCount);

  //Step #3: Read in spectra and map precursors
  //Iterate over all input files
  for (i = 0; i<files.size(); i++){

    //set up our log
    log.clear();
    param_obj.buildOutput(files[i].input, files[i].base, files[i].ext);
    log.setLog(param_obj.logFile);
    log.addMessage("Magnum version: " + string(VERSION), true);
    log.addMessage("Parameter file: " + paramFile, true);

    if (params.ext.compare(".mgf") == 0 && params.precursorRefinement){
      log.addError("Cannot perform precursor refinement using MGF files. Please disable by setting precursor_refinement=0");
      return -10;
    }

    //new file reading pipelines several steps to speed loading and transforming spectra
    log.addMessage("Reading and processing spectra data file: " + files[i].input, true);
    cout << " Reading and processing spectra data file: " << files[i].input.c_str() << " ... ";
    if (!spec.readSpectra()){
      log.addError("Error reading MS_data_file: " + files[i].input);
      return -2;
    }

    //for (size_t a = 0;a < spec.size();a++) {
    //  if (spec[a].getScanNumber() == 130224) spec[a].getPrecursor(1).monoMass = 1818.8927;
    //}

    //log.addMessage("Reading spectra data file: " + files[i].input, true);
    //cout << "\n Reading spectra data file: " << files[i].input << " ... ";
    //if (!spec.readSpectra()){
    //  log.addError("Error reading MS_data_file: " + files[i].input);
    //  return -2;
    //}
    //spec.mapPrecursors();

    //log.addMessage("Start transformation.", true);
    //time(&timeNow);
    //cout << "\n Start transformation: " << ctime(&timeNow);
    //spec.xCorr();
    //time(&timeNow);
    //cout << " Finished transformation: " << ctime(&timeNow) << endl;

    //Step #4: Analyze single peptides with open mods
    MAnalysis anal(params, &db, &spec);
    time(&timeNow);
    log.addMessage("Precompute expectation value histograms.", true);
    cout << " Precompute expectation value histograms: " << ctime(&timeNow);
    cout << "  Iterating spectra ... ";
    anal.doEValuePrecalc();
    time(&timeNow);
    cout << " Finished precompute expectation value histograms: " << ctime(&timeNow) << endl;

    log.addMessage("Start spectral search.", true);
    time(&timeNow);
    cout << " Start spectral search: " << ctime(&timeNow);
    log.addMessage("Scoring peptides.", true);
    cout << "  Scoring peptides ... ";
    anal.doPeptideAnalysis();

    log.addMessage("Finish spectral search.", true);
    time(&timeNow);
    cout << " Finished spectral search: " << ctime(&timeNow) << endl;

    //Step #5: Output results
    log.addMessage("Exporting results.", true);
    cout << " Exporting Results." << endl;
    spec.outputResults(db);

    log.addMessage("Finished Magnum analysis.", true);
    log.exportLog();

  }

  return 0;

}

bool MagnumManager::getBaseFileName(string& base, string& fName, string& extP) {
  char file[256];
  char ext[256];
  char *tok;
  char preExt[256];
  unsigned int i;

  strcpy(ext, "");

  strcpy(file, fName.c_str());
  tok = strtok(file, ".\n");
  while (tok != NULL){
    strcpy(preExt, ext);
    strcpy(ext, tok);
    tok = strtok(NULL, ".\n");
  }

  for (i = 0; i<strlen(ext); i++) ext[i] = toupper(ext[i]);
  for (i = 0; i<strlen(preExt); i++) preExt[i] = toupper(preExt[i]);

  base = fName;
  if (strcmp(ext, "MZML") == 0) {
    base=fName.substr(0,fName.size() - 5);
    extP = ".mzML";
    return true;
  }
  if (strcmp(ext, "MZXML") == 0) {
    base = fName.substr(0, fName.size() - 6);
    extP = ".mzXML";
    return true;
  }
  if (strcmp(ext, "GZ") == 0) {
    if (strcmp(preExt, "MZML") == 0){
      base = fName.substr(0, fName.size() - 8);
      extP = ".mzML.gz";
      return true;
    }
    if (strcmp(preExt, "MZXML") == 0) {
      base = fName.substr(0, fName.size() - 9);
      extP = ".mzXML.gz";
      return true;
    }
  }
  if (strcmp(ext, "RAW") == 0) {
    base = fName.substr(0, fName.size() - 4);
    extP = ".raw";
    return true;
  }
  if (strcmp(ext, "MGF") == 0) {
    base = fName.substr(0, fName.size() - 4);
    extP = ".mgf";
    return true;
  }
  if (strcmp(ext, "MS2") == 0) {
    base = fName.substr(0, fName.size() - 4);
    extP = ".ms2";
    return true;
  }
  return false;
}