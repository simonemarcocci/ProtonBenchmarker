void makeStackedCleanlinessHisto(){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

  TH1D *h1 = (TH1D*)_file0->Get("recobenchmarker/showerCleanlinessPrimaryProton");
  TH1D *h2 = (TH1D*)_file0->Get("recobenchmarker/showerCleanlinessPrimaryMuonOrPion");
  TH1D *h3 = (TH1D*)_file0->Get("recobenchmarker/showerCleanlinessPrimaryElectron");
  TH1D *h4 = (TH1D*)_file0->Get("recobenchmarker/showerCleanlinessConv");
  TH1D *h5 = (TH1D*)_file0->Get("recobenchmarker/showerCleanlinessInelastic");
  TH1D *h6 = (TH1D*)_file0->Get("recobenchmarker/showerCleanlinessMuIoni");
  TH1D *h7 = (TH1D*)_file0->Get("recobenchmarker/showerCleanlinessOther");

  h1->SetFillColor(kGreen+3);
  h2->SetFillColor(kGreen-3);
  h3->SetFillColor(kGreen-9);
  h4->SetFillColor(kAzure-3);
  h5->SetFillColor(kViolet-5);
  h6->SetFillColor(kPink+6);
  h7->SetFillColor(kOrange+7);

  THStack *hs = new THStack("hs", "");
  hs->Add(h1);
  hs->Add(h2);
  hs->Add(h3);
  hs->Add(h4);
  hs->Add(h5);
  hs->Add(h6);
  hs->Add(h7);

  hs->Draw();

  hs->GetXaxis()->SetTitle("Reco/Truth matching completeness");
  c1->Modified();
}
