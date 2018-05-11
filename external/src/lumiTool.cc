#include "lumiTool.h"

using namespace std;

lumiTool::lumiTool(string fileName)
{
  ifstream file(fileName.c_str());
  string temp;
  Size_t found1;
  //Size_t found2;
  vector<array<UInt_t, 2>> lumiVector;
  array<UInt_t, 2> lumiArray;
  UInt_t run;
  while (!file.eof()) {
    file >> temp;
    if((found1 = temp.find_first_of('"')) != std::string::npos) {
      temp = temp.substr(found1+1);
      found1 = temp.find_last_of('"');
      temp = temp.substr(0, found1);
      run = stoi(temp);
    } else if((found1 = temp.find("[[")) != std::string::npos) {
      temp = temp.substr(found1+2);
      found1 = temp.find_last_of(',');
      temp = temp.substr(0, found1);
      lumiArray[0] = stoi(temp);
    } else if((found1 = temp.find("]]")) != std::string::npos) {
      found1 = temp.find_first_of(']');
      temp = temp.substr(0, found1);
      lumiArray[1] = stoi(temp);
      lumiVector.push_back(lumiArray);
      lumiMap[run] = lumiVector;
      lumiVector.clear();
    } else if((found1 = temp.find_first_of('[')) != std::string::npos) {
      temp = temp.substr(found1+1);
      found1 = temp.find_last_of(',');
      temp = temp.substr(0, found1);
      lumiArray[0] = stoi(temp);
    } else if((found1 = temp.find_last_of(']')) != std::string::npos) {
      found1 = temp.find_first_of(']');
      temp = temp.substr(0, found1);
      lumiArray[1] = stoi(temp);
      lumiVector.push_back(lumiArray);
    }
  }
  file.close();
}

Bool_t lumiTool::LumiCheck(UInt_t run, UInt_t luminosityBlock){
  if (lumiMap.find(run) == lumiMap.end()) {
    return false;
  } else {
    for (UInt_t i = 0; i < lumiMap[run].size(); i++) {
      if (lumiMap[run][i][0] <= luminosityBlock && lumiMap[run][i][1] >= luminosityBlock) return true;
    }
    return false;
  }
}
