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


double computeTrigSFForSL(const TParticle& lep, int direction)
{
  double dPT, dAbsEta;
  double dVal, dErr;
  
  dPT = lep.Pt();
  dAbsEta = std::abs(lep.Eta());
  
  dVal = 1.0;
  dErr = 0.0;
  
  if ( std::abs(lep.GetPdgCode()) == 11 ) {
    // From https://github.com/vallot/CATTools/blob/cat80x/CatAnalyzer/python/leptonSF_cff.py
    // Using trigSF_Ele27_WPTight_Gsf
    if ( 34 <= dPT && dPT <= 36 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.894309;
        dErr = 0.012535;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.894022;
        dErr = 0.010714;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.953297;
        dErr = 0.016973;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.926017;
        dErr = 0.006186;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.985899;
        dErr = 0.007009;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.003567;
        dErr = 0.005407;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.980676;
        dErr = 0.006028;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.958228;
        dErr = 0.006130;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.946115;
        dErr = 0.006253;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.981862;
        dErr = 0.005547;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.987135;
        dErr = 0.005887;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.988235;
        dErr = 0.005385;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.921569;
        dErr = 0.006106;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.942657;
        dErr = 0.017563;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.879945;
        dErr = 0.009160;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.853403;
        dErr = 0.009419;
      }
    } else if ( 36 <= dPT && dPT <= 38 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.916021;
        dErr = 0.019659;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.914810;
        dErr = 0.008468;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.956636;
        dErr = 0.015245;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.948565;
        dErr = 0.004623;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.991963;
        dErr = 0.006736;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.004662;
        dErr = 0.006284;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.000000;
        dErr = 0.004451;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.967862;
        dErr = 0.005075;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.957213;
        dErr = 0.004585;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.988166;
        dErr = 0.005560;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.998836;
        dErr = 0.005463;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.993135;
        dErr = 0.003624;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.944048;
        dErr = 0.005374;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.990528;
        dErr = 0.018001;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.915900;
        dErr = 0.007914;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.903797;
        dErr = 0.007511;
      }
    } else if ( 38 <= dPT && dPT <= 40 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.930556;
        dErr = 0.011661;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.910579;
        dErr = 0.007633;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.980519;
        dErr = 0.012045;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.961176;
        dErr = 0.007175;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.995475;
        dErr = 0.003748;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.010369;
        dErr = 0.004306;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.004711;
        dErr = 0.004399;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.976857;
        dErr = 0.004701;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.976914;
        dErr = 0.005267;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.996491;
        dErr = 0.004062;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.006896;
        dErr = 0.004310;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.000000;
        dErr = 0.003752;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.951163;
        dErr = 0.005359;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.950633;
        dErr = 0.012276;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.912658;
        dErr = 0.007570;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.923655;
        dErr = 0.010033;
      }
    } else if ( 40 <= dPT && dPT <= 42 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.951432;
        dErr = 0.009091;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.917589;
        dErr = 0.006555;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.971322;
        dErr = 0.010822;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.971231;
        dErr = 0.006230;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 1.000000;
        dErr = 0.002751;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.011312;
        dErr = 0.007732;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.006912;
        dErr = 0.002821;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.990465;
        dErr = 0.002675;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.992797;
        dErr = 0.004455;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.997698;
        dErr = 0.002821;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.013652;
        dErr = 0.003587;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.005618;
        dErr = 0.002975;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.966667;
        dErr = 0.003426;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.986163;
        dErr = 0.025177;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.940000;
        dErr = 0.008232;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.929009;
        dErr = 0.005620;
      }
    } else if ( 42 <= dPT && dPT <= 44 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.953086;
        dErr = 0.006033;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.916067;
        dErr = 0.005683;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.984108;
        dErr = 0.009589;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.976163;
        dErr = 0.006155;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 1.003341;
        dErr = 0.003356;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.019166;
        dErr = 0.004770;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.012586;
        dErr = 0.002557;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.997636;
        dErr = 0.003540;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.991774;
        dErr = 0.003714;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 1.004571;
        dErr = 0.005238;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.007821;
        dErr = 0.003225;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.011173;
        dErr = 0.003369;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.966292;
        dErr = 0.003528;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.968254;
        dErr = 0.015217;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.937956;
        dErr = 0.006401;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.931408;
        dErr = 0.007668;
      }
    } else if ( 44 <= dPT && dPT <= 46 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.943237;
        dErr = 0.005500;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.942377;
        dErr = 0.005579;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.978261;
        dErr = 0.010890;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.975446;
        dErr = 0.003312;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 1.002210;
        dErr = 0.002465;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.009956;
        dErr = 0.003330;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.010158;
        dErr = 0.002771;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.997674;
        dErr = 0.003071;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.991879;
        dErr = 0.004482;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 1.002260;
        dErr = 0.004242;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.008869;
        dErr = 0.003324;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.005501;
        dErr = 0.002918;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.974416;
        dErr = 0.005017;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.996337;
        dErr = 0.012212;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.942584;
        dErr = 0.007848;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.945783;
        dErr = 0.007168;
      }
    } else if ( 46 <= dPT && dPT <= 48 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.950898;
        dErr = 0.009839;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.937129;
        dErr = 0.006839;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.963700;
        dErr = 0.021348;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.983278;
        dErr = 0.004324;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.996721;
        dErr = 0.002675;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.010977;
        dErr = 0.010712;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.010090;
        dErr = 0.002961;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.997698;
        dErr = 0.005875;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.994253;
        dErr = 0.006768;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.994426;
        dErr = 0.003350;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.006586;
        dErr = 0.004112;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.002186;
        dErr = 0.006138;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.981111;
        dErr = 0.004123;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.979567;
        dErr = 0.014771;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.946429;
        dErr = 0.006778;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.949881;
        dErr = 0.007014;
      }
    } else if ( 48 <= dPT && dPT <= 50 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.944976;
        dErr = 0.019899;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.946619;
        dErr = 0.016582;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.956624;
        dErr = 0.029827;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.986681;
        dErr = 0.004705;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.995662;
        dErr = 0.003599;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.015419;
        dErr = 0.006038;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.000000;
        dErr = 0.006597;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 1.002288;
        dErr = 0.008074;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.990878;
        dErr = 0.005335;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.993363;
        dErr = 0.003818;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.005453;
        dErr = 0.005019;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.998915;
        dErr = 0.004206;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.979121;
        dErr = 0.005220;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.976275;
        dErr = 0.026695;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.944575;
        dErr = 0.009395;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.936842;
        dErr = 0.009524;
      }
    } else if ( 50 <= dPT && dPT <= 52 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.943396;
        dErr = 0.023797;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.942690;
        dErr = 0.013936;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.984468;
        dErr = 0.035073;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.982533;
        dErr = 0.005092;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.993528;
        dErr = 0.006190;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.004338;
        dErr = 0.006053;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.004420;
        dErr = 0.005522;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.996602;
        dErr = 0.007696;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 1.003421;
        dErr = 0.005208;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 1.001101;
        dErr = 0.005404;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.003268;
        dErr = 0.009057;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.997845;
        dErr = 0.009837;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.980435;
        dErr = 0.004974;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.949825;
        dErr = 0.043373;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.953846;
        dErr = 0.015126;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.936121;
        dErr = 0.010360;
      }
    } else if ( 52 <= dPT && dPT <= 54 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.949173;
        dErr = 0.012329;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.952663;
        dErr = 0.013893;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.959064;
        dErr = 0.030272;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.978166;
        dErr = 0.006795;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.998914;
        dErr = 0.007013;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.008715;
        dErr = 0.010630;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.007684;
        dErr = 0.007366;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.997753;
        dErr = 0.015731;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.996606;
        dErr = 0.007975;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.991257;
        dErr = 0.006650;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.994612;
        dErr = 0.008513;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.000000;
        dErr = 0.007559;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.970716;
        dErr = 0.006769;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.971326;
        dErr = 0.021586;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.947491;
        dErr = 0.022557;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.939323;
        dErr = 0.021704;
      }
    } else if ( 54 <= dPT && dPT <= 56 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.965721;
        dErr = 0.026023;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.938462;
        dErr = 0.016546;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.968458;
        dErr = 0.034503;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.962527;
        dErr = 0.006614;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.995704;
        dErr = 0.012514;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.009761;
        dErr = 0.010541;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.992432;
        dErr = 0.009689;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.998874;
        dErr = 0.008844;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.997760;
        dErr = 0.009904;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.994536;
        dErr = 0.008284;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.990364;
        dErr = 0.008343;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.000000;
        dErr = 0.009088;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.988004;
        dErr = 0.007510;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.993865;
        dErr = 0.050939;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.946512;
        dErr = 0.019651;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.925926;
        dErr = 0.019355;
      }
    } else if ( 56 <= dPT && dPT <= 60 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.949471;
        dErr = 0.027667;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.935260;
        dErr = 0.018492;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.963051;
        dErr = 0.036252;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.970842;
        dErr = 0.006360;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.983015;
        dErr = 0.007835;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 0.998927;
        dErr = 0.017050;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.996750;
        dErr = 0.008286;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.997768;
        dErr = 0.007660;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.998881;
        dErr = 0.007578;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.990302;
        dErr = 0.009099;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.004306;
        dErr = 0.011791;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.996798;
        dErr = 0.010801;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.981720;
        dErr = 0.006372;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.969988;
        dErr = 0.042970;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.949590;
        dErr = 0.014835;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.910734;
        dErr = 0.014390;
      }
    } else if ( 60 <= dPT && dPT <= 65 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.958968;
        dErr = 0.017420;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.923777;
        dErr = 0.016708;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.958824;
        dErr = 0.032011;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.968085;
        dErr = 0.007879;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.996795;
        dErr = 0.012760;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 1.010718;
        dErr = 0.010087;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.992497;
        dErr = 0.010088;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.990196;
        dErr = 0.009188;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.991189;
        dErr = 0.010874;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.982961;
        dErr = 0.008439;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.995767;
        dErr = 0.010092;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.986243;
        dErr = 0.008919;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.978655;
        dErr = 0.007446;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.942012;
        dErr = 0.042423;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.975000;
        dErr = 0.017746;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.974087;
        dErr = 0.019242;
      }
    } else if ( 65 <= dPT && dPT <= 70 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.915743;
        dErr = 0.026703;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.924294;
        dErr = 0.023382;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.993873;
        dErr = 0.041495;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.973517;
        dErr = 0.009945;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.996809;
        dErr = 0.011161;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 0.997881;
        dErr = 0.011081;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 1.004310;
        dErr = 0.010520;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.992400;
        dErr = 0.011552;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.998904;
        dErr = 0.015207;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.985043;
        dErr = 0.010942;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.988494;
        dErr = 0.013259;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.001057;
        dErr = 0.009807;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.988248;
        dErr = 0.009801;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.916760;
        dErr = 0.043291;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.949192;
        dErr = 0.022012;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.926257;
        dErr = 0.020738;
      }
    } else if ( 70 <= dPT && dPT <= 80 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.924276;
        dErr = 0.018035;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.934389;
        dErr = 0.019687;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 1.021871;
        dErr = 0.068847;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.978033;
        dErr = 0.007662;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.983385;
        dErr = 0.007332;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 0.989572;
        dErr = 0.008309;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.996825;
        dErr = 0.008737;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 1.004334;
        dErr = 0.010752;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.997833;
        dErr = 0.011999;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.988310;
        dErr = 0.009338;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.983454;
        dErr = 0.009590;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.977226;
        dErr = 0.007828;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.994686;
        dErr = 0.011303;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.941657;
        dErr = 0.071727;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.924528;
        dErr = 0.018424;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.931461;
        dErr = 0.017974;
      }
    } else if ( 80 <= dPT && dPT <= 90 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.916484;
        dErr = 0.024069;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.924009;
        dErr = 0.020784;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 1.005821;
        dErr = 0.042667;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.974147;
        dErr = 0.010190;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.989627;
        dErr = 0.008856;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 0.994775;
        dErr = 0.012357;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.993743;
        dErr = 0.009329;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 0.991570;
        dErr = 0.009523;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.979014;
        dErr = 0.009957;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.978102;
        dErr = 0.008929;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.984342;
        dErr = 0.011569;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.993782;
        dErr = 0.009092;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.991570;
        dErr = 0.010438;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.936291;
        dErr = 0.042373;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.914595;
        dErr = 0.026190;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.905681;
        dErr = 0.020353;
      }
    } else if ( 90 <= dPT && dPT <= 100 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.930997;
        dErr = 0.028849;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.911379;
        dErr = 0.029934;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.978771;
        dErr = 0.054512;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.976166;
        dErr = 0.016696;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 1.005241;
        dErr = 0.010903;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 0.987590;
        dErr = 0.012710;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.993782;
        dErr = 0.010408;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 1.004206;
        dErr = 0.012323;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 1.001062;
        dErr = 0.015582;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.978056;
        dErr = 0.013987;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 1.000000;
        dErr = 0.016677;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.984584;
        dErr = 0.011471;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 1.029252;
        dErr = 0.014994;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.990773;
        dErr = 0.066321;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.913943;
        dErr = 0.026488;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.906149;
        dErr = 0.026050;
      }
    } else if ( 100 <= dPT && dPT <= 120 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.935275;
        dErr = 0.024665;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.927095;
        dErr = 0.025083;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.921286;
        dErr = 0.049000;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.969419;
        dErr = 0.007704;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.995855;
        dErr = 0.011704;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 0.992754;
        dErr = 0.012193;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.996865;
        dErr = 0.013515;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 1.005297;
        dErr = 0.011345;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.983264;
        dErr = 0.010946;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.969072;
        dErr = 0.009824;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.993750;
        dErr = 0.012298;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 1.000000;
        dErr = 0.010051;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.976483;
        dErr = 0.008099;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 0.953725;
        dErr = 0.046474;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.937431;
        dErr = 0.027022;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.938710;
        dErr = 0.023348;
      }
    } else if ( 120 <= dPT && dPT <= 200 ) {
      if ( -2.100 <= dAbsEta && dAbsEta <= -1.800 ) {
        dVal = 0.942541;
        dErr = 0.024346;
      } else if ( -1.800 <= dAbsEta && dAbsEta <= -1.566 ) {
        dVal = 0.925357;
        dErr = 0.024083;
      } else if ( -1.566 <= dAbsEta && dAbsEta <= -1.442 ) {
        dVal = 0.959206;
        dErr = 0.042233;
      } else if ( -1.442 <= dAbsEta && dAbsEta <= -1.100 ) {
        dVal = 0.984536;
        dErr = 0.008789;
      } else if ( -1.100 <= dAbsEta && dAbsEta <= -0.800 ) {
        dVal = 0.977574;
        dErr = 0.007045;
      } else if ( -0.800 <= dAbsEta && dAbsEta <= -0.600 ) {
        dVal = 0.977688;
        dErr = 0.007926;
      } else if ( -0.600 <= dAbsEta && dAbsEta <= -0.300 ) {
        dVal = 0.979633;
        dErr = 0.006177;
      } else if ( -0.300 <= dAbsEta && dAbsEta <= 0.000 ) {
        dVal = 1.004175;
        dErr = 0.009824;
      } else if ( 0.000 <= dAbsEta && dAbsEta <= 0.300 ) {
        dVal = 0.986570;
        dErr = 0.008635;
      } else if ( 0.300 <= dAbsEta && dAbsEta <= 0.600 ) {
        dVal = 0.978484;
        dErr = 0.007162;
      } else if ( 0.600 <= dAbsEta && dAbsEta <= 0.800 ) {
        dVal = 0.995876;
        dErr = 0.008555;
      } else if ( 0.800 <= dAbsEta && dAbsEta <= 1.100 ) {
        dVal = 0.976791;
        dErr = 0.006268;
      } else if ( 1.100 <= dAbsEta && dAbsEta <= 1.442 ) {
        dVal = 0.985611;
        dErr = 0.009068;
      } else if ( 1.442 <= dAbsEta && dAbsEta <= 1.566 ) {
        dVal = 1.005747;
        dErr = 0.045165;
      } else if ( 1.566 <= dAbsEta && dAbsEta <= 1.800 ) {
        dVal = 0.927711;
        dErr = 0.024474;
      } else if ( 1.800 <= dAbsEta && dAbsEta <= 2.100 ) {
        dVal = 0.939197;
        dErr = 0.022935;
      }
    }
  } else if ( std::abs(lep.GetPdgCode()) == 13 ) {
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonWorkInProgressAndPagResults
    // https://calderon.web.cern.ch/calderon/MuonPOG/2016dataReRecoEfficiencies/trigger/
    // EfficienciesAndSF_RunBtoF.json
    // [ IsoMu24_OR_IsoTkMu24_PtEtaBins ][ pt_abseta_ratio ] is used
    if ( 26.0 <= dPT && dPT <= 30.0 ) {
      if ( 0.0 <= dAbsEta && dAbsEta <= 0.9 ) {
        dVal = 0.980355;
        dErr = 0.000357;
      } else if ( 1.2 <= dAbsEta && dAbsEta <= 2.1 ) {
        dVal = 0.982503;
        dErr = 0.000655;
      } else if ( 2.1 <= dAbsEta && dAbsEta <= 2.4 ) {
        dVal = 0.906213;
        dErr = 0.001369;
      } else if ( 0.9 <= dAbsEta && dAbsEta <= 1.2 ) {
        dVal = 0.956334;
        dErr = 0.000794;
      }
    } else if ( 30.0 <= dPT && dPT <= 40.0 ) {
      if ( 0.0 <= dAbsEta && dAbsEta <= 0.9 ) {
        dVal = 0.984093;
        dErr = 0.000116;
      } else if ( 1.2 <= dAbsEta && dAbsEta <= 2.1 ) {
        dVal = 0.995480;
        dErr = 0.000219;
      } else if ( 2.1 <= dAbsEta && dAbsEta <= 2.4 ) {
        dVal = 0.944554;
        dErr = 0.000492;
      } else if ( 0.9 <= dAbsEta && dAbsEta <= 1.2 ) {
        dVal = 0.965957;
        dErr = 0.000203;
      }
    } else if ( 40.0 <= dPT && dPT <= 50.0 ) {
      if ( 0.0 <= dAbsEta && dAbsEta <= 0.9 ) {
        dVal = 0.984901;
        dErr = 0.000090;
      } else if ( 1.2 <= dAbsEta && dAbsEta <= 2.1 ) {
        dVal = 0.999199;
        dErr = 0.000152;
      } else if ( 2.1 <= dAbsEta && dAbsEta <= 2.4 ) {
        dVal = 0.957530;
        dErr = 0.000391;
      } else if ( 0.9 <= dAbsEta && dAbsEta <= 1.2 ) {
        dVal = 0.967856;
        dErr = 0.000134;
      }
    } else if ( 50.0 <= dPT && dPT <= 60.0 ) {
      if ( 0.0 <= dAbsEta && dAbsEta <= 0.9 ) {
        dVal = 0.985119;
        dErr = 0.000197;
      } else if ( 1.2 <= dAbsEta && dAbsEta <= 2.1 ) {
        dVal = 0.999119;
        dErr = 0.000318;
      } else if ( 2.1 <= dAbsEta && dAbsEta <= 2.4 ) {
        dVal = 0.960395;
        dErr = 0.000833;
      } else if ( 0.9 <= dAbsEta && dAbsEta <= 1.2 ) {
        dVal = 0.968172;
        dErr = 0.000279;
      }
    } else if ( 60.0 <= dPT && dPT <= 120.0 ) {
      if ( 0.0 <= dAbsEta && dAbsEta <= 0.9 ) {
        dVal = 0.983893;
        dErr = 0.000329;
      } else if ( 1.2 <= dAbsEta && dAbsEta <= 2.1 ) {
        dVal = 0.998935;
        dErr = 0.000508;
      } else if ( 2.1 <= dAbsEta && dAbsEta <= 2.4 ) {
        dVal = 0.952406;
        dErr = 0.001345;
      } else if ( 0.9 <= dAbsEta && dAbsEta <= 1.2 ) {
        dVal = 0.964896;
        dErr = 0.000496;
      }
    } else if ( 120.0 <= dPT && dPT <= 200.0 ) {
      if ( 0.0 <= dAbsEta && dAbsEta <= 0.9 ) {
        dVal = 0.976029;
        dErr = 0.001831;
      } else if ( 1.2 <= dAbsEta && dAbsEta <= 2.1 ) {
        dVal = 1.005521;
        dErr = 0.002897;
      } else if ( 2.1 <= dAbsEta && dAbsEta <= 2.4 ) {
        dVal = 0.977007;
        dErr = 0.011693;
      } else if ( 0.9 <= dAbsEta && dAbsEta <= 1.2 ) {
        dVal = 0.946721;
        dErr = 0.002094;
      }
    } else if ( 200.0 <= dPT && dPT <= 500.0 ) {
      if ( 0.0 <= dAbsEta && dAbsEta <= 0.9 ) {
        dVal = 0.983280;
        dErr = 0.003238;
      } else if ( 1.2 <= dAbsEta && dAbsEta <= 2.1 ) {
        dVal = 0.982797;
        dErr = 0.009401;
      } else if ( 2.1 <= dAbsEta && dAbsEta <= 2.4 ) {
        dVal = 0.918096;
        dErr = 0.077093;
      } else if ( 0.9 <= dAbsEta && dAbsEta <= 1.2 ) {
        dVal = 0.949595;
        dErr = 0.006390;
      }
    }
  }
  
  return dVal + direction * dErr;
}


double computeTrigSFInclusive(const TParticle& lep1, const TParticle& lep2, int direction)
{
  const int sumId = std::abs(lep1.GetPdgCode()) + std::abs(lep2.GetPdgCode());
  const int channel = sumId == 11+11 ? CH_ELEL : sumId == 13+13 ? CH_MUMU : CH_MUEL;

  if      ( channel == CH_ELEL ) return 0.953 + direction*0.009;
  else if ( channel == CH_MUMU ) return 0.948 + direction*0.002;
  else if ( channel == CH_MUEL ) return 0.975 + direction*0.004;

  return 1;
}
