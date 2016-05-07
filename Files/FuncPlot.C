TCanvas* makePlot(TH1F* h1,TString legh1, TH1F* h2, TString legh2, TString xtitle, TString ytitle, TString style = "R"){
	//Style: N - normal, R- drawn with ratio,
	if (strcmp("R",style.Data())==0 || strcmp("r",style.Data())==0){
	    
	    gStyle->SetOptTitle(0); 
	    gStyle->SetOptStat(0); 
	    //gStyle->SetLabelSize(0.04,"xy");
	    gStyle->SetLabelOffset(0.005,"xy");
	    //gStyle->SetTitleSize(0.02,"xy");
	    gStyle->SetTitleOffset(1.,"xy");
	    gStyle->SetTextFont(43);
	    //gStyle->SetTextSize(16);

	    TCanvas* c = new TCanvas("c","",200,10,700,500);
	    c->SetTicks(1,1);
	    c->cd();
	    TPad* p[2];

	    //Set margin size for pads
        Float_t padmarginX = 0.08;
        Float_t padmarginY = 0.1;
        Float_t textsize[2];

        TLine* line = new TLine();
        line->SetLineWidth(1);
        line->SetLineStyle(2);

        h1->GetYaxis()->SetTitle(ytitle.Data());
        h2->GetYaxis()->SetTitle(ytitle.Data());
        h1->GetXaxis()->SetTitle(xtitle.Data());
        h2->GetXaxis()->SetTitle(xtitle.Data());

	    p[0] = new TPad("p1", "p1", 0, 0.29, 1, 1, 0, 0, 0);
	    p[0]->SetBottomMargin(0);
	    p[0]->SetTopMargin(padmarginY);
	    p[0]->SetLeftMargin(padmarginX*1.5);
	    p[0]->SetRightMargin(padmarginX);
	    p[0]->SetTicks(1,1);
	    //p[0]->SetLogy();
	    p[0]->SetFixedAspectRatio();
	    textsize[0] = 18/(p[0]->GetWh()*p[0]->GetAbsHNDC());
	    p[0]->Draw();

	    p[1] = new TPad("p2", "p2", 0, 0, 1, 0.29, 0, 0, 0);
	    p[1]->SetBottomMargin(padmarginY*3);
	    p[1]->SetTopMargin(0);
	    p[1]->SetLeftMargin(padmarginX*1.5);
	    p[1]->SetRightMargin(padmarginX);
	    p[1]->SetTicks(1,1);
	    p[1]->SetFixedAspectRatio();
	    textsize[1] = 18/(p[1]->GetWh()*p[1]->GetAbsHNDC());
	    p[1]->Draw();

        p[0]->cd();
        h1->GetYaxis()->SetNdivisions(10);
        h2->GetYaxis()->SetNdivisions(10);
        h1->GetXaxis()->SetLabelSize(textsize[0]);
        h1->GetYaxis()->SetLabelSize(textsize[0]);
        h2->GetXaxis()->SetLabelSize(textsize[0]);
        h2->GetYaxis()->SetLabelSize(textsize[0]);
        h1->GetXaxis()->SetTitleSize(textsize[0]);
        h1->GetYaxis()->SetTitleSize(textsize[0]);
        h2->GetXaxis()->SetTitleSize(textsize[0]);
        h2->GetYaxis()->SetTitleSize(textsize[0]);
        h1->SetLineColor(2);
        h1->SetLineColor(8);
        h1->Draw("e1");
        h2->Draw("e1same");

        p[1]->cd();
        TH1F* hratio = (TH1F*)h1->Clone();
        hratio->Divide(h2);
        hratio->GetYaxis()->SetTitle("Ratio");
        hratio->GetXaxis()->SetTitle(xtitle.Data());
        hratio->GetYaxis()->SetNdivisions(5);
        hratio->GetYaxis()->SetRangeUser(0.89,1.12);
        hratio->GetXaxis()->SetLabelSize(textsize[1]);
        hratio->GetYaxis()->SetLabelSize(textsize[1]);
        hratio->GetXaxis()->SetTitleSize(textsize[1]);
        hratio->GetYaxis()->SetTitleSize(textsize[1]);
        hratio->GetYaxis()->CenterTitle();
        hratio->GetYaxis()->SetTitleOffset(0.7*textsize[0]/textsize[1]);
        hratio->SetTickLength(h1->GetTickLength()*textsize[1]/textsize[0]);
        hratio->Draw("e1");
        line->DrawLine(-180.,1.,180.,1.);
        
        p[0]->cd();
        TLegend *leg = new TLegend(0.65,0.09,0.84,0.25);
        leg->AddEntry(h1,legh1.Data(),"l");
		leg->AddEntry(h2,legh2.Data(),"l");
		leg->SetBorderSize(0);
		leg->Draw();

        return c;
	}
}