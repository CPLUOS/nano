#include "nano/external/interface/computePtEtaTable.h"

#include <iostream>


int computePtEtaTable::LoadData(std::string strPath, std::string strHistName) {
  Int_t i, j;
  
  TFile *fData;
  TH2 *hData;
  
  bool bFlipped;
  
  TH1D *hEta, *hPt, *hTmp;
  Int_t nNEta, nNPt, nTmp;
  
  Float_t fBinCX, fBinCY;
  Int_t nIdxX, nIdxY, nIdxZ;
  
  m_bValid = false;
  
  // Opening!
  fData = TFile::Open(strPath.c_str());
  if ( fData == NULL ) return -1;
  
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
        fBinCX = 0.5*(m_listEtaBin[i-1]+m_listEtaBin[i]);
        fBinCY = 0.5*(m_listPtBin[j-1]+m_listPtBin[j]);
      } else {
        fBinCX = 0.5*(m_listPtBin[j-1]+m_listPtBin[j]);
        fBinCY = 0.5*(m_listEtaBin[i-1]+m_listEtaBin[i]);
      }
      
      // There is a strange issue that the bin indices don't correspond to the displayed bins, 
      // which means, e.g., if I call hSrc.GetBinContent(1, 3) 
      // (see HLT_Ele32_eta2p1_WPTight_Gsf_FullRunRange.root)
      // then this method returns the content in (1, 1)-bin.
      // Even this phenomonon does not occur always; only some of histograms has this issue.
      // This complicated way to get the bin indices is for avoiding this problem.
      hData->GetBinXYZ(hData->FindBin(fBinCX, fBinCY), nIdxX, nIdxY, nIdxZ);
      
      m_listVal[i-1].push_back(hData->GetBinContent(nIdxX, nIdxY));
      m_listErr[i-1].push_back(hData->GetBinError(nIdxX, nIdxY));
    }
  }
  
  fData->Close();
  
  return 0;
}

double computePtEtaTable::getFactor(Float_t fPt, Float_t fEta, int direction) {
  if (!m_bValid) return -1;
  
  // Seeking the bin in which the value is contained
  auto GetIdxRange = [](Float_t fX, std::vector<Float_t> listEdges) {
    if (fX < listEdges[0]) return 0;
    else if (fX >= listEdges[listEdges.size()-1]) return ((Int_t)listEdges.size())-2;
    else return (int)(std::lower_bound(listEdges.begin(), listEdges.end(), fX)-listEdges.begin())-1;
  };
  
  Int_t nIdxPt, nIdxEta;
  
  if (m_bUseAbsEta) fEta = std::abs(fEta);
  
  nIdxPt  = GetIdxRange(fPt,  m_listPtBin);
  nIdxEta = GetIdxRange(fEta, m_listEtaBin);
  
  return m_listVal[nIdxEta][nIdxPt]+direction*m_listErr[nIdxEta][nIdxPt];
}
