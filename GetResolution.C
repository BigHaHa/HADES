{
	#include "ReadHistFunc.C"
	#include "FuncPlot.C"
	#include "FuncPlot.C"

	gROOT->Reset(); 

    gStyle->SetOptStat(0); 

	gSystem->Load("getCorr_C.so"); 

	FILE* output = fopen("Resolution.cc","wt");

	FileHist fh;
	TH1F* hPhiAB[11];
	Float_t chi;
	Float_t Rv[4][11];
	for(Int_t iC=0; iC<11; iC++){
		hPhiAB[iC] = fh.ReadTH1F("hPhiAB",iC);
	}

	for(Int_t iC=0; iC<11; iC++){
        if (hPhiAB[iC].Integral(90.,180.) == 0 || hPhiAB[iC].Integral(0.,180.)== 0 ){ chi = 0.; }
        else {
            chi = sqrt(-2.*TMath::Log(2.*(hPhiAB[iC].Integral(90.,180.)/hPhiAB[iC].Integral(0.,180.)))); 
            for(Int_t i=0; i<4; i++){ 
            	Rv[i][iC] = getCorr(i+1,chi); 
            	fprintf(output,"cRes[%i][%i]=%+5.4f;\n",i,iC,Rv[i][iC]);
            }
        }
	}



}