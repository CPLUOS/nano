#include "nano/external/interface/computePtEtaTable.h"

#include <iostream>


computePtEtaTable::computePtEtaTable(std::string strPath, std::string strHistName) {
  Int_t i, j;
  
  TFile *fData;
  TH2 *hData;
  
  bool bFlipped;
  
  TH1D *hEta, *hPt, *hTmp;
  Int_t nNEta, nNPt, nTmp;
  
  m_bValid = false;
  
  // Opening!
  fData = TFile::Open(strPath.c_str());
  if ( fData == NULL ) return;
  
  hData = (TH2 *)fData->Get(strHistName.c_str());
  m_bValid = true;
  
  // To extract the binning info
  hEta = hData->ProjectionX();
  hPt  = hData->ProjectionY();
  
  nNEta = hEta->GetNbinsX();
  nNPt  = hPt->GetNbinsX();
  
  bFlipped = false;
  
  // Maybe... some of data is flipped; eta along y and pt along x
  if (std::abs(hEta->GetBinLowEdge(nNEta+1)) > std::abs(hPt->GetBinLowEdge(nNPt+1))) {
    hTmp = hEta; hEta = hPt; hPt = hTmp;
    nTmp = nNEta; nNEta = nNPt; nNPt = nTmp;
    bFlipped = true;
  }
  
  // Some of them has eta region [-2.4, 2.4], while others have [0, 2.4]
  m_bUseAbsEta = hEta->GetBinLowEdge(1) > -0.1;
  
  // Keeping the binning info
  for (i = 1; i <= nNEta+1; i++) m_listEtaBin.push_back(hEta->GetBinLowEdge(i));
  for (i = 1; i <= nNPt +1; i++) m_listPtBin.push_back(hPt->GetBinLowEdge(i));
  
  // Making the store of weight and its uncertainty
  for (i = 1; i <= nNEta; i++) {
    m_listVal.emplace_back();
    m_listErr.emplace_back();
    
    for (j = 1; j <= nNPt; j++) {
      if (!bFlipped) {
        // Strange, the actual one starts at index 3
        m_listVal[i-1].push_back(hData->GetBinContent(i, j+2));
        m_listErr[i-1].push_back(hData->GetBinError(i, j+2));
      } else {
        m_listVal[i-1].push_back(hData->GetBinContent(j, i+2));
        m_listErr[i-1].push_back(hData->GetBinError(j, i+2));
      }
    }
  }
  
  fData->Close();
}

double computePtEtaTable::getFactor(Float_t fPt, Float_t fEta, int direction) {
  // Seeking the bin in which the value is contained
  auto GetIdxRange = [](Float_t fX, std::vector<Float_t> listEdges) {
    if (fX < listEdges[0]) return 0;
    else if (fX >= listEdges[listEdges.size()-1]) return ((Int_t)listEdges.size())-1;
    else return (int)(std::lower_bound(listEdges.begin(), listEdges.end(), fX)-listEdges.begin())-1;
  };
  
  Int_t nIdxPt, nIdxEta;
  
  if (m_bUseAbsEta) fEta = std::abs(fEta);
  
  nIdxPt  = GetIdxRange(fPt,  m_listPtBin);
  nIdxEta = GetIdxRange(fEta, m_listEtaBin);
  
  std::cout << fPt << ", " << fEta << ": (" 
    << nIdxPt << ", " << nIdxEta << ") - "
    << m_listVal[nIdxEta][nIdxPt] << ", " << m_listErr[nIdxEta][nIdxPt] << std::endl;
  return m_listVal[nIdxEta][nIdxPt]+direction*m_listErr[nIdxEta][nIdxPt];
}
