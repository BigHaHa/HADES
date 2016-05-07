{
	gROOT->Reset(); 
	TString currFName = gROOT->GetFile()->GetName();
    cout <<  gROOT->GetFile()->GetName() << endl;
    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    currFile->Print();
    gSystem->Load("getCorr_C.so"); 
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
    gStyle->SetTitleFontSize(0.1); //shifts-title-a-bit-down-
    gStyle->SetTitleFont(42,"t") ;
    gROOT->ForceStyle(); //-needed-to-apply-global-settings-(otherwice-histograms-keep-their-own-settings)--
    //----------------------------------------------------------------------------------------------------//
    TIter next( gDirectory->GetListOfKeys() ); //-dir-content--

    Int_t oncePRsum=0;

    TH1F *hEPyc[  3];
    TH1F *hOEPyc[ 3];
    TH1F *hCnCent[6];
    for(Int_t n=0;n<6;n++){
        hCnCent[n] = new TH1F(Form("hCnCent%i",n),Form("C_{%i} vs Centrality",n+1),3,1.,4.);
    }

    TF1 flows("flows", "[6]*( 1.+2.*[0]*TMath::Cos(x[0]*3.1415926535/180.)+2.*[1]*TMath::Cos(2.*x[0]*3.1415926535/180.)+2.*[2]*TMath::Cos(3.*x[0]*3.1415926535/180.)+2.*[3]*TMath::Cos(4.*x[0]*3.1415926535/180.) + [4]*TMath::Cos(5.*x[0]*3.1415/180.) + [5]*TMath::Cos(6.*x[0]*3.1415/180.) + [7]*TMath::Sin(x[0]*3.1415/180.) + [8]*TMath::Sin(2.*x[0]*3.1415/180.) )", -180., 180.);
    
    FILE* output = fopen("Fourier.txt","wt");
    //--reading-flow-histograms-from-ROOT-file-------------(begin)----
    TKey *key; 
    while( (key = (TKey*)next()) ){ 
    	//printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());
    	TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
        	TH1 *h = (TH1*)obj;
        	Char_t ch1[1024], ch2[1024];
            Int_t iC;
            sscanf(key->GetTitle(), "%s %i %s "    , ch1, &iC, ch2);
            TString s(key->GetName());
            if( s.Contains("hFOPIWEP") ){
            	printf("[EP]-META-%s[%i]--->\n",ch1,iC);
                hEPyc[iC] = new TH1F;
                hEPyc[iC] = (TH1F*)obj;
            }
        }
    }

    cout << "reading complete" << endl;

    Float_t vi[6];
    Float_t ei[6];
    TCanvas* canv   = new TCanvas("canv","PhiEP All",1000,400);
    Char_t ds[100];
    TText txt;

    canv->Divide(3,1);
    for (Int_t im=0;im<3;im++){
        canv->cd(im+1);
        Float_t summ  = hEPyc[im]->GetSumOfWeights(); 
        Int_t phibins = hEPyc[im]->GetXaxis()->GetNbins();
        hEPyc[im]->SetMinimum(0.);
        //hOEPyc[im]->SetMinimum(0.);
        //if (im == 1) hEPyc[im]->GetYaxis()->SetRangeUser(0.,15e3);
        //if (im == 3) hEPyc[im]->GetYaxis()->SetRangeUser(0.,6e4);
        //if (im  > 3) hEPyc[im]->GetYaxis()->SetRangeUser(0.,1e5);
        hEPyc[im]->Draw();
        //hOEPyc[im]->Draw("same");
        fprintf(output,"-------------META--------------\n");
        if( summ > 0 ){ 
            flows.SetLineColor(2);
            hEPyc[im]->Fit("flows");
            for(Int_t i=0; i<6; i++){
                vi[i] = flows->GetParameter(i);
                ei[i] = flows->GetParError( i);
                if(i==0) { fprintf(output,"[*]-N-Cent-[%i]\n", im); }
                fprintf(output,"Fitted: v%i=%+5.6f +/- %5.6f;\n",i+1,vi[i],ei[i]);
                sprintf(ds, "v%i=%+5.6f +/- %5.6f",i+1, vi[i], ei[i]);
                //txt.SetTextSize(12);
                txt.DrawTextNDC(0.25, 0.36+0.05  , "META");
                txt.DrawTextNDC(0.15, 0.36-0.05*i, ds);
                hCnCent[i]->SetBinContent(im+1,vi[i]);
                hCnCent[i]->SetBinError(  im+1,ei[i]);
            }
        }
        canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PhiEP/PhiEPFOPI.pdf"),"recreate");
    }
    TCanvas* canvCent = new TCanvas("canvCent","Coeff. vs Centrality bins",1000,700);
    TLatex* Ltxt = new TLatex();
    canvCent->Divide(3,2);
    for(Int_t n=0;n<6;n++){
        canvCent->cd(n+1);
        hCnCent[n]->Draw("e");
        //Ltxt->SetTextSize(0.02);
        //Ltxt->DrawLatexNDC(0.1,0.15,"#frac{dN}{d#Psi}=1+2c_{1}cos(#Psi)+2c_{2}cos(2#Psi)+2c_{3}cos(3#Psi)+2c_{4}cos(4#Psi)+c_{5}sin(#Psi)+c_{6}sin(2#Psi)");
    }
    canvCent->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PhiEP/CoeffFOPI.png"),"recreate");


    TCanvas cResol("cResol","",1200,400);
    Float_t chi; //-Ollitrault-chi--
    Float_t Rv[4][3]; //resolution-for-v[1,2,3,4]-in-[3 centrality classes]--
    cResol.Divide(3,1);
    cResol.cd(1); hPhiAB_redd.Draw(); chi = sqrt(-2.*log(2.*(hPhiAB_redd.Integral(90.,180.)/hPhiAB_redd.Integral(0.,180.)))); sprintf(ds, "#chi=%+5.6f",chi); Ltxt.DrawLatexNDC(0.5, 0.8, ds); for(Int_t i=0; i<4; i++){ Rv[i][0] = getCorr(i+1,chi); /*sprintf(ds, "Rv_{%d0}=%+5.6f",i+1,Rv[i][0]); Ltxt.DrawLatexNDC(0.5, 0.78-0.05*(i+1), ds);*/ }
    cResol.cd(2); hPhiAB_gren.Draw(); chi = sqrt(-2.*log(2.*(hPhiAB_gren.Integral(90.,180.)/hPhiAB_gren.Integral(0.,180.)))); sprintf(ds, "#chi=%+5.6f",chi); Ltxt.DrawLatexNDC(0.5, 0.8, ds); for(Int_t i=0; i<4; i++){ Rv[i][1] = getCorr(i+1,chi); /*sprintf(ds, "Rv_{%d1}=%+5.6f",i+1,Rv[i][1]); Ltxt.DrawLatexNDC(0.5, 0.78-0.05*(i+1), ds);*/}
    cResol.cd(3); hPhiAB_blue.Draw(); chi = sqrt(-2.*log(2.*(hPhiAB_blue.Integral(90.,180.)/hPhiAB_blue.Integral(0.,180.)))); sprintf(ds, "#chi=%+5.6f",chi); Ltxt.DrawLatexNDC(0.5, 0.8, ds); for(Int_t i=0; i<4; i++){ Rv[i][2] = getCorr(i+1,chi); /*sprintf(ds, "Rv_{%d2}=%+5.6f",i+1,Rv[i][2]); Ltxt.DrawLatexNDC(0.5, 0.78-0.05*(i+1), ds);*/}
    cResol.Update();
    cResol.SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PhiEP/ResolFOPI.png"),"recreate");
}