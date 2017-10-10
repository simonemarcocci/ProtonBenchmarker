void makeStackedCleanlinessHisto(){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

  TString trackorshower("track");

  TString yLabel("# "+trackorshower+"s");

  TH1D *h1 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessPrimaryProton");
  TH1D *h2 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessPrimaryMuonOrPion");
  TH1D *h3 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessPrimaryElectron");
  TH1D *h4 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessConv");
  TH1D *h5 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessInelastic");
  TH1D *h6 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessMuIoni");
  TH1D *h7 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessOther");

  h1->SetFillColor(TColor::GetColor(  136, 13, 30 ));
  h2->SetFillColor(TColor::GetColor(  249,   87,   56 ));
  h3->SetFillColor(TColor::GetColor( 238, 150,  75 ));
  h4->SetFillColor(TColor::GetColor(  244, 211, 94 ));
  h5->SetFillColor(TColor::GetColor(  203, 238, 243 ));
  h6->SetFillColor(TColor::GetColor(  93, 183, 222 ));
  h7->SetFillColor(TColor::GetColor(  13,  59,  102 ));
  
  h1->SetLineColor(kBlack);
  h2->SetLineColor(kBlack);
  h3->SetLineColor(kBlack);
  h4->SetLineColor(kBlack);
  h5->SetLineColor(kBlack);
  h6->SetLineColor(kBlack);
  h7->SetLineColor(kBlack);

  THStack *hs = new THStack("hs", "");
  hs->Add(h1);
  hs->Add(h2);
  hs->Add(h3);
  hs->Add(h4);
  hs->Add(h5);
  hs->Add(h6);
  hs->Add(h7);

  hs->Draw();
  hs->GetXaxis()->SetTitle("Reco/Truth matching cleanliness");
  hs->GetYaxis()->SetTitle(yLabel);
  hs->GetYaxis()->SetTitleOffset(1.4);
  hs->SetMaximum(hs->GetMaximum()*1.1);

  TLegend *leg = new TLegend(0.15, 0.5, 0.5, 0.85);
  leg->AddEntry(h1, "Primary Proton", "f");
  leg->AddEntry(h2, "Primary Mu/Pi", "f");
  leg->AddEntry(h3, "Primary Electron", "f");
  leg->AddEntry(h4, "Conv", "f");
  leg->AddEntry(h5, "Inelastic", "f");
  leg->AddEntry(h6, "MuIoni", "f");
  leg->AddEntry(h7, "Other", "f");

  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->Draw();

  c1->Modified();
}
