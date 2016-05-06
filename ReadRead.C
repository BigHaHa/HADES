{
	#include "ReadHistFunc.C"
	#include "FuncPlot.C"

	FileHist f;

	TH2F* hFW[11];
	TH1F* hPhiPR[3];
	for(Int_t iC=0; iC<11; iC++){
		hFW[iC] = f.ReadTH2F("hFWxyCCsmearM",iC);
	}
    hPhiPR[0] = f.ReadTH1F("hPhiPRpCent",4,0,3);
    hPhiPR[1] = f.ReadTH1F("hPhiPRRCent",4,0,3);
    hPhiPR[2] = f.ReadTH1F("hPhiPRFCent",4,0,3);

    for (Int_t i=0;i<3;i++){ hPhiPR[i]->Rebin(10); }

    TCanvas* canv;

    canv = makePlot(hPhiPR[0],"Raw",hPhiPR[1],"Recentered","#varphi-#Psi_{EP}, deg","#frac{dN}{d(#varphi-#Psi_{EP})}","R");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/phiPRR%i%i%i.png",4,0,3),"recreate");

    canv = makePlot(hPhiPR[0],"Raw",hPhiPR[2],"Flattened","#varphi-#Psi_{EP}, deg","#frac{dN}{d(#varphi-#Psi_{EP})}","R");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/phiPRF%i%i%i.png",4,0,3),"recreate");

}