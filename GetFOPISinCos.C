{
	gROOT->Reset(); 
	TString currFName = gROOT->GetFile()->GetName();
    cout <<  gROOT->GetFile()->GetName() << endl;
    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    TIter next( gDirectory->GetListOfKeys() ); //-dir-content--
    const Int_t NHARM = 6;
    const Int_t NMULT = 3;
    const Int_t NDAY  = 13;

    FILE* output = fopen("FOPICorPar.cc","wt");
    TKey *key; 
    TH2F* hCos[NHARM][NMULT];
    TH2F* hSin[NHARM][NMULT];
    while( (key = (TKey*)next()) ){ 
    	//printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());
    	TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TH2::Class() ) ) {
        	TH2 *h = (TH2*)obj;
        	Char_t ch1[1024], ch2[1024];
            Int_t iC1,iC2; // nHarmonic, nMultiplicity, nDay;
            sscanf(key->GetTitle(), "%s %i %i" , ch1, &iC1, &iC2);
            TString s(key->GetName());
            if( s.Contains("hCosFOPI") ){
            	printf("[EP]-%s[%i][%i]--->\n",ch1,iC1,iC2);
                hCos[iC1][iC2] = new TH2F;
                hCos[iC1][iC2] = (TH2F*)obj;
            }
            if( s.Contains("hSinFOPI") ){
            	printf("[EP]-%s[%i][%i]--->\n",ch1,iC1,iC2);
                hSin[iC1][iC2] = new TH2F;
                hSin[iC1][iC2] = (TH2F*)obj;
            }
        }
    }

    cout << "reading complete" << endl;
    TH1D* hCosPar[13];
    TH1D* hSinPar[13];
    fprintf(output,"//Averaged values of Sin[n*Phi] and Cos[n*Psi]\n");
    for(Int_t n=0;n<NHARM;n++){
    	for(Int_t imult=0;imult<NMULT;imult++){
            for(Int_t iday=0;iday<NDAY;iday++){
            	hCosPar[iday] = hCos[n][imult]->ProjectionY(Form(""),iday+2,iday+3-0.01);
            	hSinPar[iday] = hSin[n][imult]->ProjectionY(Form(""),iday+2,iday+3-0.01);
                fprintf(output,"SinFOPI[%i][%2i][%2i] = %2.6f; CosFOPI[%i][%2i][%2i] = %2.6f;\n",n,imult,iday,hSinPar[iday]->GetMean(),n,imult,iday,hCosPar[iday]->GetMean());
    	    }
    	}
    }
    fclose(output);
}