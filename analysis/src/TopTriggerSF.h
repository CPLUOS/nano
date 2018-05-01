#ifndef TopTriggerSF_H
#define TopTriggerSF_H
#include "TParticle.h"

double computeTrigSF(const TParticle& lep1, const TParticle& lep2, int direction=0);
double computeTrigSFInclusive(const TParticle& lep1, const TParticle& lep2, int direction=0);

#endif
