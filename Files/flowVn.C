//flowVn(TString infile)
{
    gROOT->Reset(); 

    const Int_t PrintCent = 4;
    const Int_t PrintYo   = 5;
    const Int_t PrintPt   = 17;

    gSystem->Load("getCorr_C.so"); //this loads function [Float_t getCorr(Int_t i, Float_t chi)] from external file

    Float_t PI=3.1415926535;
    Int_t MinSTATrequired=1600000; /*1600000 400000/*/ //below-this-number-of-events--errors-are-too-large--

    TString currFName = gROOT->GetFile()->GetName();
    cout <<  gROOT->GetFile()->GetName() << endl;
    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    currFile->Print();

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

    const Int_t NPt   = 18; //-if-needed-this-limits-plotted-number-of-Pt-bins--
    const Int_t NCENT = 11;
    const Int_t NYo   = 9;

    Int_t doSubtr=0;

    FILE* output = fopen("outflowEvCnt.txt","wt");

    TCanvas c("Azimuthal flow of protons [Yo]x[Pt]","Azimuthal flow of protons [Yo]x[Pt]",1290,400);
    c.Divide(4,1);
    c.cd();
    TH1* h;
    TCanvas* canv = new TCanvas("canv" ,"",900,700);
    TCanvas* canv1= new TCanvas("canv1","",900,500);
    TCanvas* canv2= new TCanvas("canv2","",900,500);

    TH1F* hPhiPR[NCENT][NYo][NPt];
    TH1F* hPhiAB[NCENT][NYo][NPt];
    TH1F* hPhiEP[NCENT][NYo][NPt];
    TH1F* hMulTr[NCENT][NYo][NPt];

    TH1F* hPhiAB_sum[NCENT][NYo];
    TH1F* hPhiPR_sum[NCENT][NYo];

    Float_t vsum[6][NCENT][NYo][NPt];
    Float_t esum[6][NCENT][NYo][NPt];
    Float_t vnvsY0[6][NCENT][NYo];
    Float_t envsY0[6][NCENT][NYo];
    Float_t vnvsPtnegY0[6][NCENT][NPt];
    Float_t envsPtnegY0[6][NCENT][NPt];
    Float_t vnvsPtposY0[6][NCENT][NPt];
    Float_t envsPtposY0[6][NCENT][NPt];
    Float_t Rv[6][NCENT];

    Char_t ch1[1024], ch2[1024], ch3[1024];
    Int_t iC, iY, iPt;

    TH1F* hv[6][NCENT][NPt]; Char_t chn[100], cht[100]; 
    TH1F* hs[6][NCENT][NPt]; //-hv-subtracted
    for(Int_t iP=0; iP<NPt; iP++){
        for(Int_t iC=0; iC<NCENT; iC++){
            for(Int_t iv=0; iv<6 /*4*/; iv++){
                for(Int_t iF=0;iF<3;iF++){
                    sprintf(chn, "v%i_y0_centr_%i_Pt_%i", iv+1, iC, iP); printf("%s\n", chn);
                    sprintf(cht, "v%i:y0 centr:%i Pt:%i", iv+1, iC, iP); printf("%s\n", cht);
                    hv[iv][iC][iP] = new TH1F(chn,cht,11, -1.1,1.1);

                    sprintf(chn, "vs%i_y0_centr_%i_Pt_%i", iv+1, iC, iP); printf("%s\n", chn);
                    sprintf(cht, "vs%i:y0 centr:%i Pt:%i", iv+1, iC, iP); printf("%s\n", cht);
                    hs[iv][iC][iP] = new TH1F(chn,cht,11, -1.1,1.1); //-a-full-copy-of-hv--later-will-be-used-for-subtraction--
                }
            }
        }
    }

    //-rapidity-projection-(Pt-integrated)--
    TH1F* hvYo[4][NCENT]; Char_t ch[100];
    for(Int_t iC=0; iC<NCENT; iC++){
       for(Int_t iv=0; iv<4; iv++){
         sprintf(ch, "v%i vs y0 centr bin %i", iv+1, iC); printf("%s\n", ch);
         hvYo[iv][iC] = new TH1F(ch,ch,11, -1.1,1.1);
       }
    }
    TH1F* hvPtposYo[4][NCENT];
    TH1F* hvPtnegYo[4][NCENT];
    Char_t ch[100];
    for(Int_t iC=0; iC<NCENT; iC++){
       for(Int_t iv=0; iv<4; iv++){
         sprintf(ch, "v%ivsPt centr bin %i Y_{0}>0", iv+1, iC); printf("%s\n", ch);
         hvPtposYo[iv][iC] = new TH1F(ch,ch,NPt, 0.2,2.);
         sprintf(ch, "v%ivsPt centr bin %i Y_{0}<0", iv+1, iC); printf("%s\n", ch);
         hvPtnegYo[iv][iC] = new TH1F(ch,ch,NPt, 0.2,2.);
        }
    }
    TH1F* hvnYo[4][NCENT][NPt];
    TH1F* hvnPt[4][NCENT][NYo];
        for(Int_t n=0;n<4;n++){
            for(Int_t iC=0;iC<NCENT;iC++){
                for(Int_t iP=0;iP<NPt;iP++){
                    sprintf(ch, "hvnYo harm %i centr_bin %i Pt_bin %i", n+1, iC, iP);
                    hvnYo[n][iC][iP] = new TH1F(Form("hvnYo%i%i%i" ,n,iC,iP),ch,NYo+2,-1.1,1.1);
                }
                for(Int_t iY=0;iY<NYo;iY++){
                    sprintf(ch, "hvnPt harm %i centr_bin %i Yo_bin %i", n+1, iC, iY);
                    hvnPt[n][iC][iY] = new TH1F(Form("hvnPt%i%i%i",n,iC,iP),ch,NPt,0.2,2.);
                }
            }
        }

    //gSystem->Exec("sleep 1000");



    TF1 flows("flows", "( 1. +  2.*[0]*TMath::Cos(x[0]*3.1415926535/180.) + 2.*[1]*TMath::Cos(2.*x[0]*3.1415926535/180.) + 2.*[2]*TMath::Cos(3.*x[0]*3.1415926535/180.) + 2.*[3]*TMath::Cos(4.*x[0]*3.1415926535/180.) + 2.*[4]*TMath::Cos(5.*x[0]*3.1415926535/180.) + 2.*[5]*TMath::Cos(6.*x[0]*3.1415926535/180.) )", -180., 180.);

    TF1 v1line(    "v1line"    , "[0] + [1]*x");
    TF1 v3line(    "v3line"    , "[0] + [1]*x");
    TF1 v2parabola("v2parabola", "[0] + sqrt([1]*[1])*x*x");
    TF1 v4parabola("v4parabola", "[0] - sqrt([1]*[1])*x*x");

    TIter next( gDirectory->GetListOfKeys() ); //-dir-content--

    Int_t oncePRsum=0;
    //--reading-flow-histograms-from-ROOT-file-------------(begin)----
    TKey *key; 
    while( (key = (TKey*)next()) ){ 
        //cout << " Classname " << key->GetClassName() << endl;
        //cout << " Name "      << key->GetName()      << endl; 
        //cout << " Title "     << key->GetTitle()     << endl;
        printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());

        TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
            TH1 *h = (TH1*)obj;
            //h = (TH1*)obj;
            //h->Draw();
            //c.Update();
            sscanf(key->GetTitle(), "%s %i %s %i %s %i"    , ch1, &iC, ch2, &iY, ch3, &iPt);
            //printf("%s - %i - %s - %i\n", ch1,  iC, ch2,  iY);
            TString s(key->GetName());
            if( s.Contains("hPhiPRp") ){ //hPhiPRp
                printf("[PR]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiPR[iC][iY][iPt] = new TH1F;
                hPhiPR[iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiABp") ){ //hPhiABp
                printf("[AB]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiAB[iC][iY][iPt] = new TH1F;
                hPhiAB[iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiEPp") ){ //hPhiEPp
                printf("[EP]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiEP[iC][iY][iPt] = new TH1F;
                hPhiEP[iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hMulTrp") ){ //hMulTrp
                printf("[Tr]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hMulTr[iC][iY][iPt] = new TH1F;
                hMulTr[iC][iY][iPt] = (TH1F*)obj;
            }
        }
    }
    //--reading-flow-histograms-from-ROOT-file-------------( end )----

    //-now-do-flow-analysis---------------------//
    Float_t ent;   //--Number-of-entries--------//
    Float_t chi;   //--Ollitrault's--chi--------//
    Float_t ci[6]; //--correction-factors-------//
    Float_t vi[6]; //--Value-of-flow-v_i--------//
    Float_t ei[6]; //--Error-of-v_i-------------//
    for(Int_t iC=4; iC<NCENT; iC++){
        for(Int_t iY=0; iY<NYo; iY++){
            for(Int_t iP=0; iP<NPt /*20 16 18*/; iP++){ //18
                chi=0.;

                c.cd(1); 
                hPhiEP[iC][iY][iP].Draw();
                //ent=hPhiEP[iC][iY]->GetEntries(); hPhiEP[iC][iY]->Scale(90./ent); hPhiEP[iC][iY]->SetMinimum(0.);

                //-------reaction-plane-resolution-for-each-selection---
                c.cd(2); 
                hPhiAB[iC][iY][iP].Draw();
                if (hPhiAB[iC][iY][iP].Integral(90.,180.) == 0 || hPhiAB[iC][iY][iP].Integral(0.,180.)== 0 ){ chi = 0.; }
                else { chi = sqrt(-2.*TMath::Log(2.*(hPhiAB[iC][iY][iP].Integral(90.,180.)/hPhiAB[iC][iY][iP].Integral(0.,180.)))); for(Int_t i=0; i<4; i++){ Rv[i][iC] = getCorr(i+1,chi); }}
                //for(Int_t i=0; i<4; i++){ ci[i] = getCorr(i+1,chi[i]); }

                //for(Int_t i=0; i<4; i++){ ci[i] = getCorr(i+1,chi); }

                for(Int_t i=0; i<4; i++){ printf("ci[%i]=%6.4f  ", i, ci[i]); }
                printf("\n");

                //-global-evet-resolution--
                for(Int_t i=0; i<4; i++){ ci[i] = Rv[i][iC];        printf("Ri[%i]=%6.4f  ", i, ci[i]); }
                printf("\n");

                ci[4]=1.0; //-temporary-set-to-unity--because-it-is--not-yet-calculated--
                ci[5]=1.0; //-temporary-set-to-unity--because-it-is--not-yet-calculated--

                c.cd(3); 
                Float_t summ  = hPhiPR[iC][iY][iP]->GetSumOfWeights(); 
                Int_t phibins = hPhiPR[iC][iY][iP]->GetXaxis()->GetNbins();
                printf("*********************************[iC=%i][iY=%i]=Entries=%i\n",iC,iY,summ);

                if (phibins>0 && summ>0 ){ hPhiPR[iC][iY][iP]->Scale(phibins/summ); }
                hPhiPR[iC][iY][iP]->SetMinimum(0.);
                hPhiPR[iC][iY][iP]->Draw(); 
                if( summ > 0 ){ 
                    hPhiPR[iC][iY][iP]->Fit("flows");
                    Char_t ds[100];
                    TText txt;
                    for(Int_t i=0; i<6 /*4*/; i++){
                        vi[i] = flows->GetParameter(i);
                        ei[i] = flows->GetParError( i);
                        if(i==0) fprintf(output,"[*]--[%i][%i][%i]\n", iC, iY, iP);
                        fprintf(output,"Fitted: v%i=%+5.4f +/- %5.4f [corr=%6.4f] Corrected: v%i=%+5.4f +/- %5.4f\n", i+1, vi[i], ei[i], ci[i], i, vi[i]/ci[i], ei[i]/ci[i]);
                        sprintf(ds, "v%i=%+5.4f +/- %5.4f",i+1, vi[i]/ci[i], ei[i]/ci[i]);
                        txt.DrawTextNDC(0.25, 0.36-0.05*i, ds);

                        if(summ>12*100 /*12*1000*/  /*4800000*/  /*12.*100*/ ){
                            if(ci[i]>0.000001){
                                hv[i][iC][iP]->SetBinContent(iY+2, vi[i]/ci[i]);
                                hv[i][iC][iP]->SetBinError(  iY+2, ei[i]/ci[i]);
                                vsum[i][iC][iY][iP] = vi[i]/ci[i];
                                esum[i][iC][iY][iP] = ei[i]/ci[i];
     
                                //hs[i][iC][iP]->SetBinContent(iY+2, vi[i]/ci[i]);
                                //hs[i][iC][iP]->SetBinError(  iY+2, ei[i]/ci[i]);
                            }
                        }
                    }
                }
                c.cd(4); 
                hMulTr[iC][iY][iP].Draw();

                c.Update();
                canv->cd();
                hPhiPR[iC][iY][iP].Draw();
                for(Int_t i=0; i<4 /*4*/; i++){
                    sprintf(ds, "v%i=%+5.4f +/- %5.4f",i+1, vi[i]/ci[i], ei[i]/ci[i]);
                    txt.DrawTextNDC(0.25, 0.36-0.05*i, ds);
                }

                /*if(iC == PrintCent && iY == PrintYo && iP == PrintPt)*/ canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Flow/Info/hPhiPR%iC%iY%iP.png",iC,iY,iP),"recreate");

                //Int_t iTmp; scanf("%i",&iTmp);
            }
        }
    }

    canv1->Divide(2,1);
    canv2->Divide(2,1);

    for(Int_t iC=4; iC<NCENT; iC++){
        for(Int_t iY=0; iY<NYo; iY++){
            for(Int_t iP=0; iP<NPt; iP++){
                for(Int_t i=0; i<4; i++){
                    hvnYo[i][iC][iP]->SetBinContent(iY+2,vsum[i][iC][iY][iP]);
                    hvnPt[i][iC][iY]->SetBinContent(iP+1,vsum[i][iC][iY][iP]);
                    hvnYo[i][iC][iP]->SetBinError(  iY+2,esum[i][iC][iY][iP]);
                    hvnPt[i][iC][iY]->SetBinError(  iP+1,esum[i][iC][iY][iP]);
                }
            }
        }
    }
    for(Int_t iC=4; iC<NCENT; iC++){
        for(Int_t iY=0; iY<NYo; iY++){
            for(Int_t i=0;i<2;i++){
                canv1->cd(i+1);
                hvnPt[i][iC][iY]->Draw("e1");
            }
            canv1->Update();
            canv1->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Flow/VnvsPt/vnPt%iC%iY.png",iC,iY),"recreate");
        }
        for(Int_t iP=0;iP<NPt;iP++){
            for(Int_t i=0;i<2;i++){
                canv2->cd(i+1);
                hvnYo[i][iC][iP]->Draw("e1");
            }
            canv2->Update();
            canv2->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Flow/VnvsYo/vnYo%iC%iP.png",iC,iP),"recreate");
        }
    }


    /*
    for(Int_t iC=4; iC<NCENT; iC++){
        for(Int_t iY=0; iY<NYo; iY++){
            for(Int_t iP=0; iP<NPt; iP++){
                for(Int_t i=0; i<6; i++){
                    vnvsY0[i][iC][iY] += vsum[i][iC][iY][iP]/NPt;
                    if(iY>=0 && iY< 4) vnvsPtnegY0[i][iC][iP] += vsum[i][iC][iY][iP]/4;
                    if(iY> 4 && iY< 8) vnvsPtposY0[i][iC][iP] += vsum[i][iC][iY][iP]/4;

                    envsY0[i][iC][iY] += esum[i][iC][iY][iP]/NPt;
                    if(iY>=0 && iY< 4) envsPtnegY0[i][iC][iP] += esum[i][iC][iY][iP]/4;
                    if(iY> 4 && iY< 8) envsPtposY0[i][iC][iP] += esum[i][iC][iY][iP]/4;
                }
            }
        }
    }
    //Integrated over Pt. hvYo vs Y0 for 3 Centrality bins
    TCanvas* cVnvsYo = new TCanvas("cVnvsYo","vn vs Yo",1200,700);
    cVnvsYo->Divide(3,2);
    TLine* zero = new TLine();
    Int_t k = 1;
    for(Int_t i=0; i<2; i++){
        for(Int_t iC=4; iC<NCENT; iC++){
            cVnvsYo->cd(k);
            for(Int_t iY=0; iY<9; iY++){
                hvYo[i][iC]->SetBinContent(iY+2,vnvsY0[i][iC][iY]);
                hvYo[i][iC]->SetBinError(  iY+2,envsY0[i][iC][iY]);
            }
            hvYo[i][iC]->SetMarkerStyle(7);
            //hvYo[i][iC]->GetYaxis()->SetTitle(Form("v_{%i}",i));
            //hvYo[i][iC]->GetXaxis()->SetTitle("Y_{0}");
            hvYo[i][iC]->Draw("e1");
            cVnvsYo->Update();
            zero->SetLineWidth(2);
            zero->SetLineStyle(2);
            zero->DrawLine(0.,gPad->GetUymin(),0.,gPad->GetUymax());
            k++;
        }
    }

    //Integrated over Y0. hvPt vs Pt for 3 Centrality bins
    TCanvas* cVnvsPt = new TCanvas("cVnvsPt","vn vs Pt",1200,700);
    cVnvsPt->Divide(3,2);
    k = 1;
    for(Int_t i=0; i<2; i++){
        for(Int_t iC=4; iC<NCENT; iC++){
            cVnvsPt->cd(k);
            for(Int_t iP=0; iP<NPt; iP++){
                hvPtposYo[i][iC]->SetBinContent(iP+1,vnvsPtposY0[i][iC][iP]);
                hvPtposYo[i][iC]->SetBinError(  iP+1,envsPtposY0[i][iC][iP]);
                hvPtnegYo[0][iC]->SetBinContent(iP+1,fabs(vnvsPtnegY0[0][iC][iP]));
                hvPtnegYo[1][iC]->SetBinContent(iP+1,vnvsPtnegY0[1][iC][iP]);
                hvPtnegYo[i][iC]->SetBinError(  iP+1,envsPtnegY0[i][iC][iP]);
            }
            hvPtnegYo[i][iC]->SetMarkerStyle(25);
            hvPtnegYo[0][iC]->SetMinimum(-0.01);
            hvPtnegYo[i][iC]->Draw("e");
            hvPtposYo[i][iC]->SetMarkerStyle(21);
            hvPtposYo[i][iC]->Draw("esame");
            k++;
        }
    }

    cVnvsYo->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Flow/VnvsYo.png"),"recreate");
    cVnvsPt->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Flow/VnvsPt.png"),"recreate");
    */
}