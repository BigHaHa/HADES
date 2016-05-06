#include <TMath.h>
#include <stdio.h>
  using namespace std;
 

void  PrintPlots(){
    
  TFile* infile = new TFile("Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/result.root");

  const Int_t npt  = 10;
  const Int_t nrpdt= 10;
  
  Int_t       div  = 4;
  const Int_t dim  = 25; //(Int_t) npt*nrpdt/div;
  
  Int_t       k    = 1;
  Int_t       l    = 0;
    
  const Double_t binpt  [npt+1]  = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
  
  const Double_t binrpdt[nrpdt+1]= {0.01,0.21,0.41,0.61,0.81,1.01,1.21,1.41,1.61,1.81,2.01};
  
  TH1F* hMAll [npt][nrpdt];
  TH1F* hMPip [npt][nrpdt];
  TH1F* hMKap [npt][nrpdt];
  TH1F* hMPro [npt][nrpdt];
  TH1F* hMPim [npt][nrpdt];
    
  TCanvas *cPip[dim];
  TCanvas *cPim[dim];
  TCanvas *cPro[dim];
  TCanvas *cKap[dim];
  
  for (Int_t i=0;i<dim;i++){
    cPip[i] = new TCanvas(Form("cPip%1.0i",i),Form("Pip %1.0i",i),1000,1000);
    cPip[i]->Divide(2,2);
    cPim[i] = new TCanvas(Form("cPim%1.0i",i),Form("Pim %1.0i",i),1000,1000);
    cPim[i]->Divide(2,2);
    cPro[i] = new TCanvas(Form("cPro%1.0i",i),Form("Pro %1.0i",i),1000,1000);
    cPro[i]->Divide(2,2);
    cKap[i] = new TCanvas(Form("cKap%1.0i",i),Form("Kap %1.0i",i),1000,1000);
    cKap[i]->Divide(2,2);
  }
    
  for (Int_t i=0;i<npt;i++){
    
    for(Int_t j=0;j<nrpdt;j++){
      hMAll [i][j] = (TH1F*)infile->Get(Form("hM2pt%dy%d" , i,j));
      
      hMPip [i][j] = (TH1F*)infile->Get(Form("hM2Pip%dy%d", i,j));
      hMPim [i][j] = (TH1F*)infile->Get(Form("hM2Pim%dy%d", i,j));
      hMKap [i][j] = (TH1F*)infile->Get(Form("hM2Kap%dy%d", i,j));
      hMPro [i][j] = (TH1F*)infile->Get(Form("hM2P%dy%d"  , i,j));
                 
    }
  }
  //std::cout << "one done, go on"<< std::endl;
  
    for (Int_t i=0;i<npt;i++){
      for(Int_t j=0;j<nrpdt;j++){

       if(l > dim){ break;    }

       hMPip[i][j]->SetLineColor(50);
       hMPip[i][j]->SetFillColor(50);
       hMPro[i][j]->SetLineColor(51);
       hMPro[i][j]->SetFillColor(51);
       hMPim[i][j]->SetLineColor(4);
       hMPim[i][j]->SetFillColor(4);
       hMKap[i][j]->SetLineColor(5);
       hMKap[i][j]->SetFillColor(5);
              
       hMPip[i][j]->SetFillStyle(3001);
       hMKap[i][j]->SetFillStyle(3001);
       hMPim[i][j]->SetFillStyle(3001);
       hMPro[i][j]->SetFillStyle(3001);
       
       cPip[l]->cd(k);
       //cPip[i]->cd(j+1)->SetLogy();
       hMPip [i][j]->GetXaxis()->SetRangeUser(-0.1,0.4);
       //hMPip [i][j]->GetYaxis()->SetRangeUser(0.5,1e6);
       hMPip [i][j]->Draw();
       //hMPro [i][j]->Draw("same");
       //hMPim [i][j]->Draw("same");
       //hMKap [i][j]->Draw("same");
       hMAll [i][j]->Draw("same");
       cPip[l]->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PID/Pip/Pip%2.0d.png",l+1),"recreate");
       
       cPim[l]->cd(k);
       //cPip[i]->cd(j+1)->SetLogy();
       hMPim [i][j]->GetXaxis()->SetRangeUser(-0.4,0.1);
       //hMPip [i][j]->GetYaxis()->SetRangeUser(0.5,1e6);
       hMPim [i][j]->Draw();
       //hMPro [i][j]->Draw("same");
       //hMPip [i][j]->Draw("same");
       //hMKap [i][j]->Draw("same");
       hMAll [i][j]->Draw("same");
       cPim[l]->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PID/Pim/Pim%2.0d.png",l+1),"recreate");
       
       cPro[l]->cd(k);
       //cPro[i]->cd(j+1)->SetLogy();
       hMPro [i][j]->GetXaxis()->SetRangeUser(0.4,1.7);
       //hMPro [i][j]->GetYaxis()->SetRangeUser(0.5,1e6);
       hMPro [i][j]->Draw();
       //hMPip [i][j]->Draw("same");
       //hMPim [i][j]->Draw("same");
       //hMKap [i][j]->Draw("same");
       hMAll [i][j]->Draw("same");
       cPro[l]->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PID/Pro/Pro%2.0d.png",l+1),"recreate");
       
       cKap[l]->cd(k);
       //cKap[i]->cd(j+1)->SetLogy();
       hMKap [i][j]->GetXaxis()->SetRangeUser(0.,1.6);
       //hMKap [i][j]->GetYaxis()->SetRangeUser(0.5,1e6);
       hMKap [i][j]->Draw();
       //hMPip [i][j]->Draw("same");
       //hMPim [i][j]->Draw("same");
       //hMPro [i][j]->Draw("same");
       hMAll [i][j]->Draw("same");
       cKap[l]->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PID/Kap/Kap%2.0d.png",l+1),"recreate");
       
       if(k >= div){ k=0; l++; }
       
       k++;
       
       //std::cout << "i= "<< i << "j= "<< j << std::endl;
       
      }
    }
  }//eof