void makeStackedCleanlinessHisto(){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

  TString trackorshower("track");

  TH1D *h1 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessPrimaryProton");
  TH1D *h2 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessPrimaryMuonOrPion");
  TH1D *h3 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessPrimaryElectron");
  TH1D *h4 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessConv");
  TH1D *h5 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessInelastic");
  TH1D *h6 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessMuIoni");
  TH1D *h7 = (TH1D*)_file0->Get("recobenchmarker/"+trackorshower+"CleanlinessOther");

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
  hs->GetXaxis()->SetTitle("Reco/Truth matching cleanliness");
  hs->SetMaximum(hs->GetMaximum()*1.1);

  TLegend *leg = new TLegend(0.15, 0.5, 0.5, 0.85);
  leg->AddEntry(h1, "Primary Proton");
  leg->AddEntry(h2, "Primary Mu/Pi");
  leg->AddEntry(h3, "Primary Electron");
  leg->AddEntry(h4, "Conv");
  leg->AddEntry(h5, "Inelastic");
  leg->AddEntry(h6, "MuIoni");
  leg->AddEntry(h7, "Other");

  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->Draw();
  
  c1->Modified();
}
