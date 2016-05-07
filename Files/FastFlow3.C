{
    gROOT->Reset(); 

    gStyle->SetOptStat(0); 

    const Int_t NPt   = 18; //-if-needed-this-limits-plotted-number-of-Pt-bins--
    const Int_t NCENT = 11;
    const Int_t NYo   = 9;
    const Int_t color[3]={2,8,9};
    Float_t vsum[3][6][NCENT][NYo][NPt];
    Float_t esum[3][6][NCENT][NYo][NPt];
    Float_t vMult[6][NCENT][NYo][NPt];
    Float_t vNtr[NCENT][NYo][NPt];
    Float_t cRes[4][NCENT];

    //#include "outflowEvCnt.cc"
    #include "outflowEvCntIII.cc"
    #include "Resolution.cc"

    for(Int_t iF=0;iF<3;iF++){
        for(Int_t n=0;n<4;n++){
            for(Int_t iC=0;iC<NCENT;iC++){
                for(Int_t iY=0;iY<NYo;iY++){
                    for(Int_t iP=0;iP<NPt-5;iP++){
                        vsum[iF][n][iC][iY][iP]/=cRes[n][iC];
                        esum[iF][n][iC][iY][iP]/=cRes[n][iC];
                    }
                }
            }
        }
    }

    const Int_t PrintCent = 4;
    const Int_t PrintYo   = 5;
    const Int_t PrintPt   = 17;

    Float_t PI=3.1415926535;
    Int_t MinSTATrequired=1600000; /*1600000 400000/*/ //below-this-number-of-events--errors-are-too-large--

    //----------------------------------------------------------------------------------------------------//  
    gStyle->SetPalette(1);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetPadRightMargin(0.16);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.08);
    gStyle->SetLabelOffset(0.015,"xyz");
    gStyle->SetLabelSize(0.042,"xyz");
    gStyle->SetLabelFont(62,"xyz") ;
    //-----------------------//------------------//
    //gStyle->SetTitleX(0.0) //title X location -//
    gStyle->SetTitleW(0.9);  //width-of-title ---//
    gStyle->SetTitleFontSize(0.05); //shifts-title-a-bit-down-
    gStyle->SetTitleFont(42,"t") ;
    gStyle->SetTitleOffset(1.15);
    gROOT->ForceStyle(); //-needed-to-apply-global-settings-(otherwice-histograms-keep-their-own-settings)--
    //----------------------------------------------------------------------------------------------------//

    //-alternative-LHCb-style--
    gROOT->LoadMacro("lhcbStyle.C");
    //lhcbStyle();

    Int_t doSubtr=0;

    TCanvas* canv = new TCanvas("canv" ,"",900,700);
    TCanvas* canv1= new TCanvas("canv1","",1350,750);
    TCanvas* canv2= new TCanvas("canv2","",1350,600);

    TH1F* hPhiPR[3][NCENT][NYo][NPt];
    TH1F* hPhiAB[3][NCENT][NYo][NPt];
    TH1F* hPhiEP[3][NCENT][NYo][NPt];
    TH1F* hMulTr[3][NCENT][NYo][NPt];
    TH1F* hNumTr[NCENT][NYo][NPt];
    TH1F* RatioPt;
    TH1F* RatMlPt;
    TH1F* RatioYo1;
    TH1F* RatioYo2;

    Char_t ch1[1024], ch2[1024], ch3[1024];
    Char_t ch[100];
    Int_t iC, iY, iPt;

    TH1F* hv[3][6][NCENT][NPt]; Char_t chn[100], cht[100]; 
    TH1F* hs[3][6][NCENT][NPt]; //-hv-subtracted
    for(Int_t iP=0; iP<NPt; iP++){
        for(Int_t iC=0; iC<NCENT; iC++){
            for(Int_t iv=0; iv<6 /*4*/; iv++){
                for(Int_t iF=0;iF<3;iF++){
                    sprintf(chn, "v%i_y0_centr_%i_Pt_%i | %i", iv+1, iC, iP,iF); printf("%s\n", chn);
                    sprintf(cht, "v%i:y0 centr:%i Pt:%i | %i", iv+1, iC, iP,iF); printf("%s\n", cht);
                    hv[iF][iv][iC][iP] = new TH1F(chn,cht,11, -1.1,1.1);

                    sprintf(chn, "vs%i_y0_centr_%i_Pt_%i | %i", iv+1, iC, iP,iF); printf("%s\n", chn);
                    sprintf(cht, "vs%i:y0 centr:%i Pt:%i | %i", iv+1, iC, iP,iF); printf("%s\n", cht);
                    hs[iF][iv][iC][iP] = new TH1F(chn,cht,11, -1.1,1.1); //-a-full-copy-of-hv--later-will-be-used-for-subtraction--
                }
            }
        }
    }

    TH1F* hvnYo[3][4][NCENT][NPt];
    TH1F* hvnPt[3][4][NCENT][NYo];
    for(Int_t iF=0;iF<3;iF++){
        for(Int_t n=0;n<4;n++){
            for(Int_t iC=0;iC<NCENT;iC++){
                for(Int_t iP=0;iP<NPt;iP++){
                    sprintf(ch, "hvnYo %i harm %i centr_bin %i Pt_bin %i",iF, n+1, iC, iP);
                    hvnYo[iF][n][iC][iP] = new TH1F(Form("hvnYo%i%i%i%i" ,iF,n,iC,iP),ch,NYo+2,-1.1,1.1);
                }
            }
        }
    }
    for(Int_t iF=0;iF<3;iF++){
        for(Int_t n=0;n<4;n++){
            for(Int_t iC=0;iC<NCENT;iC++){
                for(Int_t iY=0;iY<NYo;iY++){
                    sprintf(ch, "hvnPt %i harm %i centr_bin %i Yo_bin %i",iF, n+1, iC, iY);
                    hvnPt[iF][n][iC][iY] = new TH1F(Form("hvnPt%i%i%i%i",iF,n,iC,iY),ch,NPt,0.2,2.);
                }
            }
        }
    }
    TH1F* hMlYo[4][NCENT][NPt];
    TH1F* hMlPt[4][NCENT][NYo];
    for(Int_t n=0;n<4;n++){
        for(Int_t iC=0;iC<NCENT;iC++){
            for(Int_t iP=0;iP<NPt;iP++){
                sprintf(ch, "hvnYo %i harm %i centr_bin %i Pt_bin %i",iF, n+1, iC, iP);
                hMlYo[n][iC][iP] = new TH1F(Form("hMlYo%i%i%i" ,n,iC,iP),ch,NYo+2,-1.1,1.1);
            }
        }
    }
    for(Int_t n=0;n<4;n++){
        for(Int_t iC=0;iC<NCENT;iC++){
            for(Int_t iY=0;iY<NYo;iY++){
                sprintf(ch, "hvnPt %i harm %i centr_bin %i Yo_bin %i",iF, n+1, iC, iY);
                hMlPt[n][iC][iY] = new TH1F(Form("hMlPt%i%i%i" ,n,iC,iY),ch,NPt,0.2,2.);
            }
        }
    }
    canv1->Divide(3,2);
    canv2->Divide(3,1);

    for(Int_t iF=0;iF<3;iF++){
        for(Int_t i=0; i<2; i++){
            for(Int_t iC=0; iC<NCENT-4; iC++){
                for(Int_t iY=0; iY<NYo; iY++){
                    for(Int_t iP=0; iP<NPt-5; iP++){
                        if (iY>(NYo-1)/2 || i!=0) hvnPt[iF][i][iC][iY]->SetBinContent(iP+1, vsum[iF][i][iC][iY][iP]);
                        if (iY>(NYo-1)/2 || i!=0) hvnPt[iF][i][iC][iY]->SetBinError(  iP+1, esum[iF][i][iC][iY][iP]);
                        if (iY<(NYo-1)/2 && i==0) hvnPt[iF][i][iC][iY]->SetBinContent(iP+1,-vsum[iF][i][iC][iY][iP]);
                        if (iY<(NYo-1)/2 && i==0) hvnPt[iF][i][iC][iY]->SetBinError(  iP+1, esum[iF][i][iC][iY][iP]);
                    }
                    hvnPt[iF][i][iC][iY]->GetYaxis()->SetTitle(Form("v_{%i}",i+1));
                    hvnPt[iF][i][iC][iY]->GetXaxis()->SetTitle(Form("p_{T}, GeV"));
                }
            }
        }
    }

    for(Int_t iF=0;iF<3;iF++){
        for(Int_t i=0; i<2; i++){
            for(Int_t iC=0; iC<NCENT-4; iC++){
                for(Int_t iP=0; iP<NPt-5; iP++){
                    for(Int_t iY=0; iY<NYo; iY++){
                        hvnYo[iF][i][iC][iP]->SetBinContent(iY+2,vsum[iF][i][iC][iY][iP]);
                        hvnYo[iF][i][iC][iP]->SetBinError(  iY+2,esum[iF][i][iC][iY][iP]);
                    }
                    hvnYo[iF][i][iC][iP]->GetYaxis()->SetTitle(Form("v_{%i}",i+1));
                    hvnYo[iF][i][iC][iP]->GetXaxis()->SetTitle(Form("Y_{0}"));
                }
            }
        }
    }

    for(Int_t i=0; i<2; i++){
        for(Int_t iC=0; iC<NCENT-4; iC++){
            for(Int_t iY=0; iY<NYo; iY++){
                for(Int_t iP=0; iP<NPt-5; iP++){
                    hMlPt[i][iC][iY]->SetBinContent(iP+1, vNtr[iC][iY][iP]);
                }
                hMlPt[i][iC][iY]->GetYaxis()->SetTitle(Form("N_{tr}#timesN_{evts}"));
                hMlPt[i][iC][iY]->GetXaxis()->SetTitle(Form("p_{T}, GeV"));
            }
        }
    }

    for(Int_t i=0; i<2; i++){
        for(Int_t iC=0; iC<NCENT-4; iC++){
            for(Int_t iP=0; iP<NPt-5; iP++){
                for(Int_t iY=0; iY<NYo; iY++){
                    hMlYo[i][iC][iP]->SetBinContent(iY+2,vNtr[iC][iY][iP]);
                    if (iP ==0) hMlYo[i][iC][iP]->GetYaxis()->SetRangeUser(0.,2e6);
                }
                hMlYo[i][iC][iP]->GetYaxis()->SetTitle(Form("N_{tr}#timesN_{evts}"));
                hMlYo[i][iC][iP]->GetXaxis()->SetTitle(Form("Y_{0}"));
            }
        }
    }

    TLine* line1 = new TLine();
    TLine* line2 = new TLine();
    TLine* line3 = new TLine();
    line1->SetLineWidth(1);
    line1->SetLineStyle(1);
    line2->SetLineWidth(1);
    line2->SetLineStyle(2);
    line3->SetLineWidth(3);
    line3->SetLineStyle(2);

    for(Int_t iC=0; iC<NCENT-4; iC++){
        for(Int_t iY=0; iY<(NYo-1)/2; iY++){
            for(Int_t i=0;i<2;i++){
                count = 1;
                for(Int_t iF=0;iF<3;iF++){
                    canv1->cd(i+1);
                    hvnPt[iF][i][iC][iY]->SetLineColor(color[iF]);
                    hvnPt[iF][i][iC][iY]->SetMarkerColor(color[iF]);
                    hvnPt[iF][i][iC][iY]->SetMarkerStyle(25);
                    hvnPt[iF][i][iC][NYo-1-iY]->SetMarkerStyle(21);
                    hvnPt[iF][i][iC][NYo-1-iY]->SetLineColor(color[iF]);
                    hvnPt[iF][i][iC][NYo-1-iY]->SetMarkerColor(color[iF]);
                    if(iF==0) hvnPt[iF][i][iC][iY]->Draw("e1");
                    if(iF!=0) hvnPt[iF][i][iC][iY]->Draw("e1same");
                    hvnPt[iF][i][iC][NYo-1-iY]->Draw("e1same");

                    canv1->cd(i+4);
                    RatioPt = (TH1F*)hvnPt[iF][i][iC][NYo-1-iY]->Clone();
                    RatioPt->Divide(hvnPt[iF][i][iC][iY]);
                    RatioPt->GetYaxis()->SetRangeUser(-0.1,3.);
                    if(iF==0) RatioPt->Draw("e1");
                    if(iF!=0) RatioPt->Draw("e1same");
                    line1->DrawLine(0.2,1.,2.,1.);
                    line2->DrawLine(0.2,1.1,2.,1.1);
                    line2->DrawLine(0.2,0.9,2.,0.9);
                }
                canv1->cd(3);
                hMlPt[i][iC][iY]->SetLineColor(color[2]);
                hMlPt[i][iC][iY]->SetMarkerColor(color[2]);
                hMlPt[i][iC][iY]->SetMarkerStyle(25);
                hMlPt[i][iC][NYo-1-iY]->SetMarkerStyle(21);
                hMlPt[i][iC][NYo-1-iY]->SetLineColor(color[2]);
                hMlPt[i][iC][NYo-1-iY]->SetMarkerColor(color[2]);
                hMlPt[i][iC][iY]->Draw("p");
                hMlPt[i][iC][NYo-1-iY]->Draw("psame");

                canv1->cd(6);
                RatMlPt = (TH1F*)hMlPt[i][iC][NYo-1-iY]->Clone();
                RatMlPt->Divide(hMlPt[i][iC][iY]);
                RatMlPt->GetYaxis()->SetRangeUser(-0.1,3.);
                RatMlPt->Draw("p");
                line1->DrawLine(0.2,1.,2.,1.);
                line2->DrawLine(0.2,1.1,2.,1.1);
                line2->DrawLine(0.2,0.9,2.,0.9);
            }
            canv1->SetLogy()
            canv1->Update();
            canv1->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/FlowDetailed/VnvsPt/vnPt%iC%iY.png",iC,iY),"recreate");
        }
        for(Int_t iP=0;iP<NPt-5;iP++){
            for(Int_t i=0;i<2;i++){
                for(Int_t iF=0;iF<3;iF++){
                    canv2->cd(i+1);
                    hvnYo[iF][i][iC][iP]->SetLineColor(color[iF]);
                    //hvnYo[iF][i][iC][iP]->SetMarkerColor(color[iF]);
                    if(iF==0) hvnYo[iF][i][iC][iP]->Draw("e1");
                    if(iF!=0) hvnYo[iF][i][iC][iP]->Draw("e1same");
                    canv2->Update();
                    line3->DrawLine(0.,gPad->GetUymin(),0.,gPad->GetUymax());
                }
                canv2->cd(3);
                hMlYo[i][iC][iP]->SetLineColor(color[2]);
                hMlYo[i][iC][iP]->SetMarkerColor(color[2]);
                hMlYo[i][iC][iP]->SetMarkerStyle(21);
                hMlYo[i][iC][iP]->Draw("p");
                if (iP ==0) line3->DrawLine(0.,0.,0.,2e6);
                line3->DrawLine(0.,0.,0.,gPad->GetUymax());
            }
            canv2->SetLogy();
            canv2->Update();
            canv2->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/FlowDetailed/VnvsYo/vnYo%iC%iP.png",iC,iP),"recreate");
        }
    }
}
