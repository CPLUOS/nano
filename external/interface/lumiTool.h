#ifndef LUMITOOL_H
#define LUMITOOL_H

#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <map>
#include <array>
#include <string>

class lumiTool {
public:
  lumiTool(std::string fileName);
  ~lumiTool();
  Bool_t LumiCheck(UInt_t run, UInt_t luminosityBlock);
private:
  std::map<UInt_t, std::vector<std::array<UInt_t, 2>>> lumiMap;
};

#endif
