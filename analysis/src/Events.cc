// Include the Events() definitions as part of this C++ file
#define Events_cxx
#include "nano/analysis/interface/Events.h"

void Events::Loop(){}

void makeEventsClass (const char* filedir) {
  TFile *f = TFile::Open(filedir);
  TTree *t = (TTree*) f->Get("Events");
  t->MakeClass(); // this will generate Events.h file.
}
