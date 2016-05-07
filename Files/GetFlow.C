void InitgStyle(){
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
}
class flowhist{
 public:
 	flowhist(){
 		InitHists();
 	}

    void InitHists(){
	 //-hv-subtracted
        for(Int_t iP=0; iP<NPt; iP++){
            for(Int_t iC=0; iC<NCENT; iC++){
                for(Int_t iv=0; iv<6; iv++){
                    sprintf(chn, "v%i_y0_centr_%i_Pt_%i", iv+1, iC, iP);
                    sprintf(cht, "v%i:y0 centr:%i Pt:%i", iv+1, iC, iP);
                    hv[iv][iC][iP] = new TH1F(chn,cht,11, -1.1,1.1);

                    sprintf(chn, "vs%i_y0_centr_%i_Pt_%i", iv+1, iC, iP);
                    sprintf(cht, "vs%i:y0 centr:%i Pt:%i", iv+1, iC, iP);
                    hs[iv][iC][iP] = new TH1F(chn,cht,11, -1.1,1.1); //-a-full-copy-of-hv--later-will-be-used-for-subtraction--
                }
            }
        }
        for(Int_t n=0;n<4;n++){
            for(Int_t iC=0;iC<NCENT;iC++){
                for(Int_t iP=0;iP<NPt;iP++){
                    sprintf(ch, "hvnYo harm %i centr_bin %i Pt_bin %i", n+1, iC, iP);
                    hvnYo[n][iC][iP] = new TH1F(Form("hvnYo%i%i%i",n,iC,iP),ch,NYo+2,-1.1,1.1);
                }
                for(Int_t iY=0;iY<NYo;iY++){
                    sprintf(ch, "hvnPt harm %i centr_bin %i Yo_bin %i", n+1, iC, iY);
                    hvnPt[n][iC][iY] = new TH1F(Form("hvnPt%i%i%i",n,iC,iP),ch,NPt,0.2,2.);
                }
            }
        }
    }

    void ReadHist(TString currFName){
	    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    
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
    }

    void GetPhiPR(Int_t iC, Int_t iY, Int_t iPt, TH1F* h){
    	h = hPhiPR[iC][iY][iPt];
    }

    void GetPhiAB(Int_t iC, Int_t iY, Int_t iPt, TH1F* h){
    	h = hPhiAB[iC][iY][iPt];
    }

    void GetPhiEP(Int_t iC, Int_t iY, Int_t iPt, TH1F* h){
    	h = hPhiEP[iC][iY][iPt];
    }

    void GetMulTr(Int_t iC, Int_t iY, Int_t iPt, TH1F* h){
    	h = hMulTr[iC][iY][iPt];
    }

    void FitHist(TH1F* hPR, TH1F* hAB, TH1F* hEP){
    	for(Int_t iC=4; iC<NCENT; iC++){
            for(Int_t iY=0; iY<NYo; iY++){
                for(Int_t iP=0; iP<NPt; iP++){ //18
                    chi=0.;

                    if (hAB[iC][iY][iP].Integral(90.,180.) == 0 || hAB[iC][iY][iP].Integral(0.,180.)== 0 ){ chi = 0.; }
                    else { chi = sqrt(-2.*TMath::Log(2.*(hAB[iC][iY][iP].Integral(90.,180.)/hAB[iC][iY][iP].Integral(0.,180.)))); for(Int_t i=0; i<4; i++){ Rv[i][iC] = getCorr(i+1,chi); }}

                    //-global-evet-resolution--
                    for(Int_t i=0; i<4; i++){ ci[i] = Rv[i][iC]; }
                    ci[4]=1.0; //-temporary-set-to-unity--because-it-is--not-yet-calculated--
                    ci[5]=1.0;
                    Float_t summ  = hPR[iC][iY][iP]->GetSumOfWeights(); 
                    Int_t phibins = hPR[iC][iY][iP]->GetXaxis()->GetNbins();
                    if( summ > 0 ){ 
                        hPR[iC][iY][iP]->Fit("flows");
                        Char_t ds[100];
                        TText txt;
                        for(Int_t i=0; i<6; i++){
                            vi[i] = flows->GetParameter(i);
                            ei[i] = flows->GetParError( i);
                            sprintf(ds, "v%i=%+5.4f +/- %5.4f",i+1, vi[i]/ci[i], ei[i]/ci[i]);
                            txt.DrawTextNDC(0.25, 0.36-0.05*i, ds);

                            if(summ>12*100 ){
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
                }
            }
        }
    }

    void BuildVn(){
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
    }

    void GetVnYo(TH1F* h, Int_t i, Int_t iC, Int_t iP){
    	h = hvnYo[i][iC][iP];
    }

    void GetVnPt(TH1F* h, Int_t i, Int_t iC, Int_t iY){
    	h = hvnPt[i][iC][iY];
    }
 
    const Int_t NPt   = 18; //-if-needed-this-limits-plotted-number-of-Pt-bins--
    const Int_t NCENT = 11;
    const Int_t NYo   = 9;
    Float_t ent;   //--Number-of-entries--------//
    Float_t chi;   //--Ollitrault's--chi--------//
    Float_t ci[6]; //--correction-factors-------//
    Float_t vi[6]; //--Value-of-flow-v_i--------//
    Float_t ei[6]; //--Error-of-v_i-------------//
    TH1F* hPhiPR[NCENT][NYo][NPt];
    TH1F* hPhiAB[NCENT][NYo][NPt];
    TH1F* hPhiEP[NCENT][NYo][NPt];
    TH1F* hMulTr[NCENT][NYo][NPt];

    TH1F* hv[6][NCENT][NPt]; Char_t chn[100], cht[100]; 
    TH1F* hs[6][NCENT][NPt];

    TH1F* hvnYo[4][NCENT][NPt];
    TH1F* hvnPt[4][NCENT][NYo];
    //TF1 flows("flows", "( 1. +  2.*[0]*TMath::Cos(x[0]*3.1415926535/180.) + 2.*[1]*TMath::Cos(2.*x[0]*3.1415926535/180.) + 2.*[2]*TMath::Cos(3.*x[0]*3.1415926535/180.) + 2.*[3]*TMath::Cos(4.*x[0]*3.1415926535/180.) + 2.*[4]*TMath::Cos(5.*x[0]*3.1415926535/180.) + 2.*[5]*TMath::Cos(6.*x[0]*3.1415926535/180.) )", -180., 180.);

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
    Float_t summ ;
    Int_t phibins;
};

void GetFlow(){
	gROOT->Reset();
	gSystem->Load("getCorr_C.so");
	Float_t PI=3.1415926535;
	const Int_t NPt   = 18;
    const Int_t NCENT = 11;
    const Int_t NYo   = 9;

	TH1F* hPsiPR[1][NCENT][NYo][NPt];
	TH1F* hPsiAB[1][NCENT][NYo][NPt];
	TH1F* hPsiEP[1][NCENT][NYo][NPt];
	TH1F* hMulTr[1][NCENT][NYo][NPt];

	InitgStyle();

	flowhist flow;

	//flow.ReadHist("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/12apr108.root");
	/*for(Int_t iF=0;iF<1;iF++){
		for(Int_t iC=0;iC<NCENT;iC++){
			for(Int_t iY=0;iY<NYo;iY++){
				for(Int_t iP=0;iP<NPt;iP++){
					flow.GetPhiPR(iC,iY,iP,hPsiPR[iF][iC][iY][iP]);
					flow.GetPhiAB(iC,iY,iP,hPsiAB[iF][iC][iY][iP]);
					flow.GetPhiEP(iC,iY,iP,hPsiEP[iF][iC][iY][iP]);
					flow.GetMulTr(iC,iY,iP,hMulTr[iF][iC][iY][iP]);
				}
			}
		}
	}*/

}