#include "nano/external/interface/TopTriggerSF.h"
#include "TParticle.h"
#include "nano/external/interface/TTbarModeDefs.h"
#include <cmath>

double computeTrigSF(const TParticle& lep1, const  TParticle& lep2, int direction)
{
  const int sumId = std::abs(lep1.GetPdgCode()) + std::abs(lep2.GetPdgCode());
  const int channel = sumId == 11+11 ? CH_ELEL : sumId == 13+13 ? CH_MUMU : CH_MUEL;

  auto getEta = [](const TParticle& lep)->double {
    if ( std::abs(lep.GetPdgCode()) == 11 ) return dynamic_cast<const TParticle&>(lep).GetWeight();
    return lep.Eta();
  };
  const double aeta1 = std::abs(getEta(lep1));
  const double aeta2 = std::abs(getEta(lep2));

  if ( channel == CH_ELEL ) {
    if ( aeta1 < 0.3 ) {
      if      ( aeta2 < 0.3 ) return 0.965 + direction*0.010;
      else if ( aeta2 < 0.6 ) return 0.989 + 0.5*direction*0.007;
      else if ( aeta2 < 1.2 ) return 1.003 + 0.5*direction*0.004;
      else if ( aeta2 < 1.7 ) return 0.990 + 0.5*direction*0.010;
      else if ( aeta2 < 2.4 ) return 0.997 + 0.5*direction*0.011;
    }
    else if ( aeta1 < 0.6 ) {
      if      ( aeta2 < 0.3 ) return 0.981 + direction*0.007;
      else if ( aeta2 < 0.6 ) return 0.972 + 0.5*direction*0.009;
      else if ( aeta2 < 1.2 ) return 0.990 + 0.5*direction*0.005;
      else if ( aeta2 < 1.7 ) return 1.000 + 0.5*direction*0.007;
      else if ( aeta2 < 2.4 ) return 0.997 + 0.5*direction*0.011;
    }
    else if ( aeta1 < 1.2 ) {
      if      ( aeta2 < 0.3 ) return 0.996 + direction*0.004;
      else if ( aeta2 < 0.6 ) return 0.989 + 0.5*direction*0.006;
      else if ( aeta2 < 1.2 ) return 0.985 + 0.5*direction*0.005;
      else if ( aeta2 < 1.7 ) return 0.982 + 0.5*direction*0.008;
      else if ( aeta2 < 2.4 ) return 0.993 + 0.5*direction*0.008;
    }
    else if ( aeta1 < 1.7 ) {
      if      ( aeta2 < 0.3 ) return 0.988 + direction*0.012;
      else if ( aeta2 < 0.6 ) return 0.997 + 0.5*direction*0.009;
      else if ( aeta2 < 1.2 ) return 0.995 + 0.5*direction*0.007;
      else if ( aeta2 < 1.7 ) return 0.984 + 0.5*direction*0.015;
      else if ( aeta2 < 2.4 ) return 0.997 + 0.5*direction*0.012;
    }
    else if ( aeta1 < 2.4 ) {
      if      ( aeta2 < 0.3 ) return 1.005 + direction*0.013;
      else if ( aeta2 < 0.6 ) return 0.995 + 0.5*direction*0.016;
      else if ( aeta2 < 1.2 ) return 0.989 + 0.5*direction*0.011;
      else if ( aeta2 < 1.7 ) return 0.985 + 0.5*direction*0.018;
      else if ( aeta2 < 2.4 ) return 0.996 + 0.5*direction*0.014;
    }
  }
  else if ( channel == CH_MUMU ) {
    if ( aeta1 < 0.3 ) {
      if      ( aeta2 < 0.3 ) return 0.990 + direction*0.006;
      else if ( aeta2 < 0.6 ) return 0.992 + 0.5*direction*0.005;
      else if ( aeta2 < 1.2 ) return 0.989 + 0.5*direction*0.005;
      else if ( aeta2 < 1.7 ) return 0.972 + 0.5*direction*0.014;
      else if ( aeta2 < 2.4 ) return 0.988 + 0.5*direction*0.016;
    }
    else if ( aeta1 < 0.6 ) {
      if      ( aeta2 < 0.3 ) return 1.001 + direction*0.004;
      else if ( aeta2 < 0.6 ) return 0.996 + 0.5*direction*0.005;
      else if ( aeta2 < 1.2 ) return 0.993 + 0.5*direction*0.004;
      else if ( aeta2 < 1.7 ) return 1.003 + 0.5*direction*0.006;
      else if ( aeta2 < 2.4 ) return 0.981 + 0.5*direction*0.016;
    }
    else if ( aeta1 < 1.2 ) {
      if      ( aeta2 < 0.3 ) return 1.000 + direction*0.004;
      else if ( aeta2 < 0.6 ) return 0.999 + 0.5*direction*0.003;
      else if ( aeta2 < 1.2 ) return 0.995 + 0.5*direction*0.003;
      else if ( aeta2 < 1.7 ) return 1.003 + 0.5*direction*0.003;
      else if ( aeta2 < 2.4 ) return 0.993 + 0.5*direction*0.007;
    }
    else if ( aeta1 < 1.7 ) {
      if      ( aeta2 < 0.3 ) return 0.985 + direction*0.013;
      else if ( aeta2 < 0.6 ) return 1.002 + 0.5*direction*0.006;
      else if ( aeta2 < 1.2 ) return 1.000 + 0.5*direction*0.003;
      else if ( aeta2 < 1.7 ) return 0.995 + 0.5*direction*0.004;
      else if ( aeta2 < 2.4 ) return 0.995 + 0.5*direction*0.006;
    }
    else if ( aeta1 < 2.4 ) {
      if      ( aeta2 < 0.3 ) return 0.998 + direction*0.020;
      else if ( aeta2 < 0.6 ) return 0.993 + 0.5*direction*0.017;
      else if ( aeta2 < 1.2 ) return 0.984 + 0.5*direction*0.009;
      else if ( aeta2 < 1.7 ) return 0.985 + 0.5*direction*0.006;
      else if ( aeta2 < 2.4 ) return 0.979 + 0.5*direction*0.007;
    }
  }
  else if ( channel == CH_MUEL ) {
    if ( aeta1 < 0.3 ) {
      if      ( aeta2 < 0.3 ) return 0.989 + direction*0.008;
      else if ( aeta2 < 0.6 ) return 0.997 + 0.5*direction*0.006;
      else if ( aeta2 < 1.2 ) return 1.000 + 0.5*direction*0.004;
      else if ( aeta2 < 1.7 ) return 0.998 + 0.5*direction*0.007;
      else if ( aeta2 < 2.4 ) return 0.981 + 0.5*direction*0.012;
    }
    else if ( aeta1 < 0.6 ) {
      if      ( aeta2 < 0.3 ) return 0.995 + direction*0.006;
      else if ( aeta2 < 0.6 ) return 0.992 + 0.5*direction*0.006;
      else if ( aeta2 < 1.2 ) return 1.000 + 0.5*direction*0.004;
      else if ( aeta2 < 1.7 ) return 1.003 + 0.5*direction*0.006;
      else if ( aeta2 < 2.4 ) return 1.000 + 0.5*direction*0.009;
    }
    else if ( aeta1 < 1.2 ) {
      if      ( aeta2 < 0.3 ) return 0.990 + direction*0.005;
      else if ( aeta2 < 0.6 ) return 0.999 + 0.5*direction*0.004;
      else if ( aeta2 < 1.2 ) return 0.997 + 0.5*direction*0.003;
      else if ( aeta2 < 1.7 ) return 1.001 + 0.5*direction*0.005;
      else if ( aeta2 < 2.4 ) return 1.004 + 0.5*direction*0.006;
    }
    else if ( aeta1 < 1.7 ) {
      if      ( aeta2 < 0.3 ) return 0.996 + direction*0.008;
      else if ( aeta2 < 0.6 ) return 1.002 + 0.5*direction*0.006;
      else if ( aeta2 < 1.2 ) return 0.997 + 0.5*direction*0.005;
      else if ( aeta2 < 1.7 ) return 0.990 + 0.5*direction*0.011;
      else if ( aeta2 < 2.4 ) return 0.975 + 0.5*direction*0.016;
    }
    else if ( aeta1 < 2.4 ) {
      if      ( aeta2 < 0.3 ) return 0.992 + direction*0.011;
      else if ( aeta2 < 0.6 ) return 0.986 + 0.5*direction*0.012;
      else if ( aeta2 < 1.2 ) return 0.989 + 0.5*direction*0.008;
      else if ( aeta2 < 1.7 ) return 1.006 + 0.5*direction*0.010;
      else if ( aeta2 < 2.4 ) return 1.002 + 0.5*direction*0.016;
    }
  }

  return 1;
};

double computeTrigSFInclusive(const TParticle& lep1, const TParticle& lep2, int direction)
{
  const int sumId = std::abs(lep1.GetPdgCode()) + std::abs(lep2.GetPdgCode());
  const int channel = sumId == 11+11 ? CH_ELEL : sumId == 13+13 ? CH_MUMU : CH_MUEL;

  if      ( channel == CH_ELEL ) return 0.953 + direction*0.009;
  else if ( channel == CH_MUMU ) return 0.948 + direction*0.002;
  else if ( channel == CH_MUEL ) return 0.975 + direction*0.004;

  return 1;
}
