#include "Events.h"

void Events::Loop(){}

void makeEventsClass (const char* filedir) {
  TFile *f = TFile::Open(filedir);
  TTree *t = (TTree*) f->Get("Events");
  t->MakeClass(); // this will generate Events.h file.
}

