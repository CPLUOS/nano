#ifndef computePtEtaTable_H
#define computePtEtaTable_H


#include <TFile.h>
#include <TH2F.h>


class computePtEtaTable {
public:
  bool m_bValid;
  
  std::vector<std::vector<Float_t>> m_listVal;
  std::vector<std::vector<Float_t>> m_listErr;
  
  std::vector<Float_t> m_listPtBin;
  std::vector<Float_t> m_listEtaBin;
  
  bool m_bUseAbsEta;
  
public:
  computePtEtaTable(std::string strPath, std::string strHistName);
  ~computePtEtaTable() {};
  
  bool isValid() {return m_bValid;};
  
  double getFactor(Float_t fPt, Float_t fEta, int direction=0);
};


#endif // computePtEtaTable_H


