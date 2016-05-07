{
	gROOT->Reset(); 
	TString currFName = gROOT->GetFile()->GetName();
    cout <<  gROOT->GetFile()->GetName() << endl;
    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    TIter next( gDirectory->GetListOfKeys() ); //-dir-content--
    const Int_t NHARM = 6;
    const Int_t NMULT = 11;
    const Int_t NDAY  = 13;

    FILE* output = fopen("Recenter.cc","wt");
    TKey *key; 
    TProfile* hsumX[NMULT];
    TProfile* hsumY[NMULT];
    while( (key = (TKey*)next()) ){ 
    	//printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());
    	TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TProfile::Class() ) ) {
        	TProfile *h = (TProfile*)obj;
        	Char_t ch1[1024], ch2[1024];
            Int_t iC1; // nMultiplicity;
            sscanf(key->GetTitle(), "%s %i" , ch1, &iC1);
            TString s(key->GetName());
            if( s.Contains("hsumXmean") ){
            	printf("[EP]-%s[%i]--->\n",ch1,iC1);
                hsumX[iC1] = new TProfile;
                hsumX[iC1] = (TProfile*)obj;
            }
            if( s.Contains("hsumYmean") ){
            	printf("[EP]-%s[%i]--->\n",ch1,iC1);
                hsumY[iC1] = new TProfile;
                hsumY[iC1] = (TProfile*)obj;
            }
        }
    }

    cout << "reading complete" << endl;
    fprintf(output,"//X and Y shift in terms of normalized TVector2;\n");
    for(Int_t imult=0;imult<NMULT;imult++){
        for(Int_t iday=0;iday<NDAY;iday++){
            fprintf(output,"sumXmean[%2i][%2i] = %2.10f; sumYmean[%2i][%2i] = %2.10f; sumXsigma[%2i][%2i] = %2.10f; sumYsigma[%2i][%2i] = %2.10f;\n",imult,iday,hsumX[imult]->GetBinContent(iday+2),imult,iday,hsumY[imult]->GetBinContent(iday+2),imult,iday,hsumX[imult]->GetBinError(iday+2)*sqrt(hsumX[imult]->GetEntries()),imult,iday,hsumY[imult]->GetBinError(iday+2)*sqrt(hsumY[imult]->GetEntries()));
    	}
    }
    fclose(output);
}