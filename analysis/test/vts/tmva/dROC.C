#include <string>
#include <sstream>

void dROC (const char *a)
{
	//TString fname( a );
	TFile *f = TFile::Open(a, "READ");
	//f->cd("Method_BDT/GradBDT1Tree5000");

	//TH1 *h;
	f->cd("Method_BDT/GradBDT7Bagged0.7");
	//f->cd("Method_BDT/GradTree1K_Dep4");

	//f->cd("Method_BDT/BDT");
	//gDirectory->GetObject("MVA_GradBDT1Tree5000_effS",h);
	//TH1 *h = (TH1*)gDirectory->Get("MVA_GradBDT7Bagged0.7_effS");
	//TH1 *h = (TH1*)gDirectory->Get("MVA_GradBDT7Bagged0.7_effB");

	//TH1 *h = (TH1*)gDirectory->Get("MVA_BDT_rejBvsS");

	TH1 *h = (TH1*)gDirectory->Get("MVA_GradBDT7Bagged0.7_rejBvsS");
	//TH1 *h = (TH1*)gDirectory->Get("MVA_GradTree1K_Dep4_rejBvsS");
/*
	Float_t sum = 0.;

	for (Int_t k = 0; k < h3->GetNbinsX(); k++){
		sum += h3->GetBinContent(k) * 0.01;
	}
*/
	Int_t n = h->GetNbinsX();

	Double_t auc = 0.;

	for (Int_t i = 0; i<n; i++){
		auc += h->GetBinContent(i)*0.001;
	}

//	std::ostringstream auc_str;
//	auc_str << auc;
//	std::string auc_str = auc_str.str();

	TString train_type;
        TString test_type;
        TString var_type;

	auto auc_str = std::to_string(auc);
	std::cout << " data trained on (dijet, z+jet)" << endl;
	std::cin >> train_type;
	std::cout << " data tested on (dijet, z+jet)" << endl;
	std::cin >> test_type;
	std::cout << " variable set is (All, CMS, Moment, PT, etc.)" << endl;
	std::cin >> var_type;
	TString filename = "BDT_delphes_train_" + train_type + "_test_" + test_type + "_" + "var_" + var_type + "_AUC_" + auc_str + ".csv";

/*	TString file    = a;
	TString rocint = sum;

	TString outputfile = file + TString("_Delphes_AUC_") + rocint + TString(".csv");// + rocint;
*/
	ofstream myfile;
	myfile.open( filename );
	//myfile.open ( outputfile );
	myfile << "# x(Sig Efficiency), y(Bkg Rejection)" << endl;
	for (Int_t i=0; i<n; ++i){
       // for (Int_t i=0; i<n; i++){
	//	myfile << printf("%g %g\n", h->GetBinLowEdge(i)+h->GetBinWidth(i)/2, h->GetBinContent(i));
		myfile << h->GetBinLowEdge(i)+h->GetBinWidth(i)/2 << ", " <<  h->GetBinContent(i) << endl;
	//	myfile << h1->GetBinContent(i) << ", " << 1 - h2->GetBinContent(i) << endl;
	}
	myfile.close();
}
