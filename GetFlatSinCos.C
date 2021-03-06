{
	gROOT->Reset(); 
	TString currFName = gROOT->GetFile()->GetName();
    cout <<  gROOT->GetFile()->GetName() << endl;
    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    TIter next( gDirectory->GetListOfKeys() ); //-dir-content--
    const Int_t NHARM = 6;
    const Int_t NMULT = 11;
    const Int_t NDAY  = 13;

    FILE* output = fopen("FlatSinCos.cc","wt");
    TKey *key; 
    TProfile* hCos[3][NHARM][NMULT];
    TProfile* hSin[3][NHARM][NMULT];
    while( (key = (TKey*)next()) ){ 
    	//printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());
    	TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TProfile::Class() ) ) {
        	TProfile *h = (TProfile*)obj;
        	Char_t ch1[1024], ch2[1024];
            Int_t iC1,iC2,iC3; // nHarmonic, nMultiplicity, nDay;
            sscanf(key->GetTitle(), "%s %i %i %i" , ch1, &iC1, &iC2, &iC3);
            TString s(key->GetName());
            if( s.Contains("hCosPsi") ){
            	printf("[EP]-%s[%i][%i][%i]--->\n",ch1,iC1,iC2,iC3);
                hCos[iC1][iC2][iC3] = new TProfile;
                hCos[iC1][iC2][iC3] = (TProfile*)obj;
            }
            if( s.Contains("hSinPsi") ){
            	printf("[EP]-%s[%i][%i][%i]--->\n",ch1,iC1,iC2,iC3);
                hSin[iC1][iC2][iC3] = new TProfile;
                hSin[iC1][iC2][iC3] = (TProfile*)obj;
            }
        }
    }

    cout << "reading complete" << endl;
    fprintf(output,"//Averaged values of Sin[n*Phi] and Cos[n*Psi]\n");
    for(Int_t iT=0;iT<3;iT++){
        for(Int_t n=0;n<NHARM;n++){
    	    for(Int_t imult=0;imult<NMULT;imult++){
                for(Int_t iday=0;iday<NDAY;iday++){
                    fprintf(output,"FlatCos[%i][%i][%2i][%2i] = %2.6f; FlatSin[%i][%i][%2i][%2i] = %2.6f;\n",iT,n,imult,iday,hCos[iT][n][imult]->GetBinContent(iday+2),iT,n,imult,iday,hSin[iT][n][imult]->GetBinContent(iday+2));
    	        }
    	    }
        }
    }
    fclose(output);
}