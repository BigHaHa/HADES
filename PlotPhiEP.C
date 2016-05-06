#include <TMath.h>
#include <stdio.h>
  using namespace std;
 

void  PlotPhiEP(){

     TFile* infile1 = new TFile("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/result.root");
     TFile* infile2 = new TFile("/home/peter/Documents/WorkLocal/WorkFiles/HADES/GitFiles/result.root");

     TH1F* hPhiEP[   4];
     TH1F* hPhiEPrec[4];

     for (Int_t n =1;n<5;n++){
     	hPhiEP[   n] = (TH1F*)infile1->Get(Form("hPhiEP%i",n-1));
     	hPhiEPrec[n] = (TH1F*)infile2->Get(Form("hPhiEP%i",n-1));
          hPhiEP[   n]-> GetXaxis()->SetTitle(Form("%i #Psi_{%i}(EP), deg",n,n));
          hPhiEPrec[n]-> GetXaxis()->SetTitle(Form("%i #Psi_{%i}(EP), deg",n,n));
     }

     TCanvas* c = new TCanvas("c","",1000,700);

     c->Divide(2,2);

     c->cd();

     for (Int_t n=1;n<5;n++){
     	c->cd(n);
     	hPhiEP[   n]->DrawNormalized();
     	//hPhiEPrec[n]->DrawNormalized("same");
     }

     c->SaveAs("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PhiEP.png","recreate");

}
