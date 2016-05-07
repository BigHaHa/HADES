//flowVn(TString infile)
{
    gROOT->Reset(); 

    gStyle->SetOptStat(0); 

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
    //gStyle->SetTitleW(0.9);  //width-of-title ---//
    //gStyle->SetTitleFontSize(0.05); //shifts-title-a-bit-down-
    //gStyle->SetTitleFont(42,"t") ;
    gStyle->SetTitleSize(0.04);
    //gStyle->SetTitleOffset(1.15);
    gROOT->ForceStyle(); //-needed-to-apply-global-settings-(otherwice-histograms-keep-their-own-settings)--
    //----------------------------------------------------------------------------------------------------//

    //-alternative-LHCb-style--
    gROOT->LoadMacro("lhcbStyle.C");
    //lhcbStyle();

    const Int_t NPt   = 18; //-if-needed-this-limits-plotted-number-of-Pt-bins--
    const Int_t NCENT = 11;
    const Int_t NYo   = 9;
    const Int_t color[3]={2,8,9};

    Int_t doSubtr=0;

    FILE* output = fopen("outflowEvCntIII.cc","wt");

    TCanvas c("Azimuthal flow of protons [Yo]x[Pt]","Azimuthal flow of protons [Yo]x[Pt]",1290,400);
    c.Divide(4,1);
    c.cd();
    TH1* h;
    TCanvas* canv = new TCanvas("canv" ,"",900,700);
    TCanvas* canv1= new TCanvas("canv1","",900,500);
    TCanvas* canv2= new TCanvas("canv2","",900,500);

    TH1F* hPhiPR[3][NCENT][NYo][NPt];
    TH1F* hPhiAB[3][NCENT][NYo][NPt];
    TH1F* hPhiEP[3][NCENT][NYo][NPt];
    TH1F* hMulTr[3][NCENT][NYo][NPt];
    TH1F* hM2[   2][NCENT][NYo][NPt];
    TH1F* hNtr[     NCENT][NYo][NPt];

    TH1F* hPhiAB_sum[NCENT][NYo];
    TH1F* hPhiPR_sum[NCENT][NYo];

    Float_t vsum[3][6][NCENT][NYo][NPt];
    Float_t esum[3][6][NCENT][NYo][NPt];
    Float_t vnvsY0[3][6][NCENT][NYo];
    Float_t envsY0[3][6][NCENT][NYo];
    Float_t vnvsPtnegY0[6][NCENT][NPt];
    Float_t envsPtnegY0[6][NCENT][NPt];
    Float_t vnvsPtposY0[6][NCENT][NPt];
    Float_t envsPtposY0[6][NCENT][NPt];
    Float_t Rv[3][6][NCENT];

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
                    sprintf(ch, "hvnPt %i harm %i centr_bin %i Yo_bin %i", n+1, iC, iY);
                    hvnPt[iF][n][iC][iY] = new TH1F(Form("hvnPt%i%i%i%i",iF,n,iC,iY),ch,NPt,0.2,2.);
                }
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
                hPhiPR[0][iC][iY][iPt] = new TH1F;
                hPhiPR[0][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiABp") ){ //hPhiABp
                printf("[AB]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiAB[0][iC][iY][iPt] = new TH1F;
                hPhiAB[0][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiEPp") ){ //hPhiEPp
                printf("[EP]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiEP[0][iC][iY][iPt] = new TH1F;
                hPhiEP[0][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hMulTrp") ){ //hMulTrp
                printf("[Tr]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hMulTr[0][iC][iY][iPt] = new TH1F;
                hMulTr[0][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiPRR") ){ //hPhiPRp
                printf("[PR]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiPR[1][iC][iY][iPt] = new TH1F;
                hPhiPR[1][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiABR") ){ //hPhiABp
                printf("[AB]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiAB[1][iC][iY][iPt] = new TH1F;
                hPhiAB[1][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiEPR") ){ //hPhiEPp
                printf("[EP]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiEP[1][iC][iY][iPt] = new TH1F;
                hPhiEP[1][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hMulTrR") ){ //hMulTrp
                printf("[Tr]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hMulTr[1][iC][iY][iPt] = new TH1F;
                hMulTr[1][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiPRF") ){ //hPhiPRp
                printf("[PR]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiPR[2][iC][iY][iPt] = new TH1F;
                hPhiPR[2][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiABF") ){ //hPhiABp
                printf("[AB]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiAB[2][iC][iY][iPt] = new TH1F;
                hPhiAB[2][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hPhiEPF") ){ //hPhiEPp
                printf("[EP]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hPhiEP[2][iC][iY][iPt] = new TH1F;
                hPhiEP[2][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hMulTrF") ){ //hMulTrp
                printf("[Tr]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hMulTr[2][iC][iY][iPt] = new TH1F;
                hMulTr[2][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hRawM2") ){ //hPhiPRp
                printf("[PR]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hM2[0][iC][iY][iPt] = new TH1F;
                hM2[0][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hProM2") ){ //hPhiPRp
                printf("[PR]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hM2[1][iC][iY][iPt] = new TH1F;
                hM2[1][iC][iY][iPt] = (TH1F*)obj;
            }
            if( s.Contains("hNumTr") ){ //hPhiPRp
                printf("[PR]--%s[%i]%s[%i]%s[%i]--->\n",ch1,iC, ch2,iY, ch3,iPt);
                hNtr[iC][iY][iPt] = new TH1F;
                hNtr[iC][iY][iPt] = (TH1F*)obj;
            }
        }
    }
    //--reading-flow-histograms-from-ROOT-file-------------( end )----

    TCanvas* cMass2[NCENT-4];
    for (Int_t i=0; i<NCENT-4;i++){
    	cMass2[i] = new TCanvas(Form("cMass2%i",i+4),Form("Mass2 %i iC",i),1350,750);
    	cMass2[i]->Divide(2,2);
    }
    for(Int_t iC=0;iC<NCENT-4;iC++){
    	cMass2[iC]->cd(1);
    	cMass2[iC]->cd(1)->SetLogy();
    	hM2[0][iC][1][4]->SetTitle(Form("m^{2}, %i cent bin, -0.7<Y_{0}<-0.5, 0.5<p_{T}<0.6 GeV/c",iC+4));
    	hM2[0][iC][1][4]->GetYaxis()->SetTitle("N_{count}");
    	hM2[0][iC][1][4]->GetXaxis()->SetTitle("m^{2}, GeV^{2}/c^{4}");
    	hM2[0][iC][1][4]->Draw();
    	hM2[1][iC][1][4]->SetLineColor(51);
    	hM2[1][iC][1][4]->SetFillColor(51);
    	hM2[1][iC][1][4]->SetFillStyle(3001);
    	hM2[1][iC][1][4]->Draw("same");
    	cMass2[iC]->cd(2);
    	cMass2[iC]->cd(2)->SetLogy();
    	hM2[0][iC][2][4]->SetTitle(Form("m^{2}, %i cent bin, -0.5<Y_{0}<-0.3, 0.5<p_{T}<0.6 GeV/c",iC+4));
    	hM2[0][iC][2][4]->GetYaxis()->SetTitle("N_{count}");
    	hM2[0][iC][2][4]->GetXaxis()->SetTitle("m^{2}, GeV^{2}/c^{4}");
    	hM2[0][iC][2][4]->Draw();
    	hM2[1][iC][2][4]->SetLineColor(51);
    	hM2[1][iC][2][4]->SetFillColor(51);
    	hM2[1][iC][2][4]->SetFillStyle(3001);
    	hM2[1][iC][2][4]->Draw("same");
    	cMass2[iC]->cd(3);
    	cMass2[iC]->cd(3)->SetLogy();
    	hM2[0][iC][7][4]->SetTitle(Form("m^{2}, %i cent bin, 0.5<Y_{0}<0.7, 0.5<p_{T}<0.6 GeV/c",iC+4));
    	hM2[0][iC][7][4]->GetYaxis()->SetTitle("N_{count}");
    	hM2[0][iC][7][4]->GetXaxis()->SetTitle("m^{2}, GeV^{2}/c^{4}");
    	hM2[0][iC][7][4]->Draw();
    	hM2[1][iC][7][4]->SetLineColor(51);
    	hM2[1][iC][7][4]->SetFillColor(51);
    	hM2[1][iC][7][4]->SetFillStyle(3001);
    	hM2[1][iC][7][4]->Draw("same");
    	cMass2[iC]->cd(4);
    	cMass2[iC]->cd(4)->SetLogy();
    	hM2[0][iC][6][4]->SetTitle(Form("m^{2}, %i cent bin, 0.3<Y_{0}<0.5, 0.5<p_{T}<0.6 GeV/c",iC+4));
    	hM2[0][iC][7][4]->GetYaxis()->SetTitle("N_{count}");
    	hM2[0][iC][7][4]->GetXaxis()->SetTitle("m^{2}, GeV^{2}/c^{4}");
    	hM2[0][iC][6][4]->Draw();
    	hM2[1][iC][6][4]->SetLineColor(51);
    	hM2[1][iC][6][4]->SetFillColor(51);
    	hM2[1][iC][6][4]->SetFillStyle(3001);
    	hM2[1][iC][6][4]->Draw("same");

    	cMass2[iC]->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PID/Proton/hMass2%iC.png",iC),"recreate");
    }

    //-now-do-flow-analysis---------------------//
    Float_t ent;   //--Number-of-entries--------//
    Float_t chi;   //--Ollitrault's--chi--------//
    Float_t ci[6]; //--correction-factors-------//
    Float_t vi[6]; //--Value-of-flow-v_i--------//
    Float_t ei[6]; //--Error-of-v_i-------------//
    for(Int_t iC=0; iC<NCENT-4; iC++){
        for(Int_t iY=0; iY<NYo; iY++){
            for(Int_t iP=0; iP<NPt-5 /*20 16 18*/; iP++){ //18
                for(Int_t iF=0;iF<3;iF++){
                    chi=0.;

                    c.cd(1); 
                    hPhiEP[iF][iC][iY][iP].SetLineColor(color[iF]);
                    if(iF==0)hPhiEP[iF][iC][iY][iP].Draw();
                    if(iF!=0)hPhiEP[iF][iC][iY][iP].Draw("same");
                    //ent=hPhiEP[iC][iY]->GetEntries(); hPhiEP[iC][iY]->Scale(90./ent); hPhiEP[iC][iY]->SetMinimum(0.);

                    //-------reaction-plane-resolution-for-each-selection---
                    c.cd(2); 
                    hPhiAB[iF][iC][iY][iP].SetLineColor(color[iF]);
                    if(iF==0) hPhiAB[iF][iC][iY][iP].Draw();
                    if(iF!=0) hPhiAB[iF][iC][iY][iP].Draw("same");
                    if (hPhiAB[iF][iC][iY][iP].Integral(90.,180.) == 0 || hPhiAB[iF][iC][iY][iP].Integral(0.,180.)== 0 ){ chi = 0.; }
                    else { chi = sqrt(-2.*TMath::Log(2.*(hPhiAB[iF][iC][iY][iP].Integral(90.,180.)/hPhiAB[iF][iC][iY][iP].Integral(0.,180.)))); for(Int_t i=0; i<4; i++){ Rv[iF][i][iC] = getCorr(i+1,chi); }}
                    //for(Int_t i=0; i<4; i++){ ci[i] = getCorr(i+1,chi[i]); }

                    //for(Int_t i=0; i<4; i++){ ci[i] = getCorr(i+1,chi); }

                    for(Int_t i=0; i<4; i++){ printf("ci[%i]=%6.4f  ", i, ci[i]); }
                    printf("\n");

                    //-global-evet-resolution--
                    for(Int_t i=0; i<4; i++){ ci[i] = Rv[iF][i][iC];        printf("Ri[%i]=%6.4f  ", i, ci[i]); }
                    printf("\n");

                    ci[4]=1.0; //-temporary-set-to-unity--because-it-is--not-yet-calculated--
                    ci[5]=1.0; //-temporary-set-to-unity--because-it-is--not-yet-calculated--

                    c.cd(3); 
                    Float_t summ  = hPhiPR[iF][iC][iY][iP]->GetSumOfWeights(); 
                    Int_t phibins = hPhiPR[iF][iC][iY][iP]->GetXaxis()->GetNbins();
                    printf("*************************[iF=%i][iC=%i][iY=%i][iP=%i]=Entries=%i\n",iF,iC,iY,iP,summ);

                    if (phibins>0 && summ>0 ){ hPhiPR[iF][iC][iY][iP]->Scale(phibins/summ); }
                    hPhiPR[iF][iC][iY][iP]->SetMinimum(0.);
                    hPhiPR[iF][iC][iY][iP]->SetLineColor(color[iF]);
                    if(iF==0) hPhiPR[iF][iC][iY][iP]->Draw(); 
                    if(iF!=0) hPhiPR[iF][iC][iY][iP]->Draw("same");
                    if( summ > 0 ){ 
                        flows.SetLineColor(color[iF]);
                        hPhiPR[iF][iC][iY][iP]->Fit("flows","Q");
                        Char_t ds[100];
                        TText txt;
                        for(Int_t i=0; i<6 /*4*/; i++){
                            vi[i] = flows->GetParameter(i);
                            ei[i] = flows->GetParError( i);
                            if(vi[i] < -10. && vi[i] >10. || ci[i]<=0.000001 ){ vi[i] = 0.;}
                            if(ei[i] < -10. && ei[i] >10. || ci[i]<=0.000001 ){ ei[i] = 0.;}
                            if(i==0) fprintf(output,"//[*]--[%i][%i][%i][%i];\n", iF, iC, iY, iP);
                            if(iF==0 )fprintf(output,"vMult[%i][%i][%i][%i]=%+10.0f;//----------------------------------------\n",i,iC,iY,iP,hMulTr[iF][iC][iY][iP].Integral());
                            if(iF==0 && i==0 )fprintf(output,"vNtr[%i][%i][%i]=%+5.4f;//------------------------------------------\n",iC,iY,iP,hNtr[iC][iY][iP].Integral());
                            fprintf(output,"vsum[%i][%i][%i][%i][%i]=%+5.4f; esum[%i][%i][%i][%i][%i]=%+5.4f;\n",iF,i,iC,iY,iP,vi[i], iF,i,iC,iY,iP,ei[i]);
                            sprintf(ds, "v%i=%+5.4f +/- %5.4f",i+1, vi[i]/ci[i], ei[i]/ci[i]);
                            txt.DrawTextNDC(0.25, 0.36-0.05*i, ds);

                            if(summ>12*100 /*12*1000*/  /*4800000*/  /*12.*100*/ ){
                                if(ci[i]>0.000001){
                                    hv[iF][i][iC][iP]->SetBinContent(iY+2, vi[i]/ci[i]);
                                    hv[iF][i][iC][iP]->SetBinError(  iY+2, ei[i]/ci[i]);
                                    vsum[iF][i][iC][iY][iP] = vi[i]/ci[i];
                                    esum[iF][i][iC][iY][iP] = ei[i]/ci[i];
     
                                    //hs[i][iC][iP]->SetBinContent(iY+2, vi[i]/ci[i]);
                                    //hs[i][iC][iP]->SetBinError(  iY+2, ei[i]/ci[i]);
                                }
                            }
                        }
                    }
                    c.cd(4); 
                    hMulTr[iF][iC][iY][iP].SetLineColor(color[iF]);
                    if(iF==0) hMulTr[iF][iC][iY][iP].Draw();
                    if(iF!=0) hMulTr[iF][iC][iY][iP].Draw("same");

                    c.Update();
                    canv->cd();
                    if(iF==0) hPhiPR[iF][iC][iY][iP].Draw();
                    if(iF!=0) hPhiPR[iF][iC][iY][iP].Draw("same");
                    for(Int_t i=0; i<4 /*4*/; i++){
                        sprintf(ds, "v%i=%+5.4f +/- %5.4f",i+1, vi[i]/ci[i], ei[i]/ci[i]);
                        txt.SetTextSize(0.025);
                        txt.DrawTextNDC(0.15+0.05  , 0.36+0.05  , "Orig");
                        txt.DrawTextNDC(0.15+0.25  , 0.36+0.05  , "Recnt");
                        txt.DrawTextNDC(0.15+0.45  , 0.36+0.05  , "Flatt");
                        txt.DrawTextNDC(0.15+0.2*iF, 0.36-0.05*i, ds);
                    }
                }

                /*if(iC == PrintCent && iY == PrintYo && iP == PrintPt)*/ 
                canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/FlowDetailed/Info/hPhiPR%iC%iY%iP.png",iC,iY,iP),"recreate");

                //Int_t iTmp; scanf("%i",&iTmp);
            }
        }
    }

    /*canv1->Divide(2,1);
    canv2->Divide(2,1);

    for(Int_t iF=0;iF<3;iF++){
        for(Int_t i=0; i<2; i++){
            for(Int_t iC=4; iC<NCENT; iC++){
                for(Int_t iY=0; iY<NYo; iY++){
                    for(Int_t iP=0; iP<NPt; iP++){
                        if (iY>(NYo-1)/2 || i!=0) hvnPt[iF][i][iC][iY]->SetBinContent(iP+1, vsum[iF][i][iC][iY][iP]);
                        if (iY>(NYo-1)/2 || i!=0) hvnPt[iF][i][iC][iY]->SetBinError(  iP+1, esum[iF][i][iC][iY][iP]);
                        if (iY<(NYo-1)/2 && i==0) hvnPt[iF][i][iC][iY]->SetBinContent(iP+1,-vsum[iF][i][iC][iY][iP]);
                        if (iY<(NYo-1)/2 && i==0) hvnPt[iF][i][iC][iY]->SetBinError(  iP+1, esum[iF][i][iC][iY][iP]);
                    }
                }
            }
        }
    }

    for(Int_t iF=0;iF<3;iF++){
        for(Int_t i=0; i<2; i++){
            for(Int_t iC=4; iC<NCENT; iC++){
                for(Int_t iP=0; iP<NPt; iP++){
                    for(Int_t iY=0; iY<NYo; iY++){
                        hvnYo[iF][i][iC][iP]->SetBinContent(iY+2,vsum[iF][i][iC][iY][iP]);
                        hvnYo[iF][i][iC][iP]->SetBinError(  iY+2,esum[iF][i][iC][iY][iP]);
                    }
                }
            }
        }
    }

    

    for(Int_t iC=4; iC<NCENT; iC++){
        for(Int_t iY=0; iY<(NYo-1)/2; iY++){
            for(Int_t i=0;i<2;i++){
                for(Int_t iF=0;iF<3;iF++){
                    canv1->cd(i+1);
                    hvnPt[iF][i][iC][iY]->SetLineColor(color[iF]);
                    hvnPt[iF][i][iC][iY]->SetMarkerStyle(25);
                    hvnPt[iF][i][iC][NYo-1-iY]->SetMarkerStyle(21);
                    hvnPt[iF][i][iC][NYo-1-iY]->SetLineColor(color[iF]);
                    if(iF==0) hvnPt[iF][i][iC][iY]->Draw("e1");
                    if(iF!=0) hvnPt[iF][i][iC][iY]->Draw("e1same");
                              hvnPt[iF][i][iC][NYo-1-iY]->Draw("e1same");
                }
            }
            canv1->Update();
            canv1->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/FlowDetailed/VnvsPt/vnPt%iC%iY.png",iC,iY),"recreate");
        }
        for(Int_t iP=0;iP<NPt;iP++){
            for(Int_t i=0;i<2;i++){
                for(Int_t iF=0;iF<3;iF++){
                    canv2->cd(i+1);
                    hvnYo[iF][i][iC][iP]->SetLineColor(color[iF]);
                    if(iF==0) hvnYo[iF][i][iC][iP]->Draw("e1");
                    if(iF!=0) hvnYo[iF][i][iC][iP]->Draw("e1same");
                }
            }
            canv2->Update();
            canv2->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/FlowDetailed/VnvsYo/vnYo%iC%iP.png",iC,iP),"recreate");
        }
    }*/
}