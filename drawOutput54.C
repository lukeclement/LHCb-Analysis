void drawOutput() {
    //////////////////////
    // Example Root Macro for third year B->hhh Lab
    ////////////////////// 

    // Open the root file that was produced by running the example program
    TFile *f = new TFile("outputDataAll.root");

    // Get pointers to the example histograms that were made 
    TH1F *hx = (TH1F*)f->Get("h_PX");
    TH1F *hy = (TH1F*)f->Get("h_PY");
    TH1F *hz = (TH1F*)f->Get("h_PZ");
    TH1F *hPP = (TH1F*)f->Get("h_PP");
    TH1F *hPK= (TH1F*)f->Get("h_PK");
    TH1F *hm = (TH1F*)f->Get("h_M");
    TH1F *hrm= (TH1F*)f->Get("h_RM");
    TH1F *hrmp= (TH1F*)f->Get("h_RMP");

    TH2F *hxy= (TH2F*)f->Get("h_PKP");
    // Create a canvas onto which the histograms are plotted and which can be saved
    //TCanvas *c1 = new TCanvas("c1","",1200,800);
    // Draw the first histogram with a blue line, and an x-axis title
    //hPP->SetLineColor(kBlue);
    //hPP->GetXaxis()->SetTitle("Probability");
    //hPP->Draw();
    //hPK->SetLineColor(kRed);
    //hPK->Draw("Same");
    // Save the canvas as pdf file. Other possible formats include root (for later editing) eps, png.
    //c1->SaveAs("CutProbs.pdf");
    //c1->SaveAs("CutProbs.root");
    
    TCanvas *c2 = new TCanvas("c2","",1200,800);
    hrm->SetLineColor(kBlue);
    hrm->GetXaxis()->SetTitle("Mass/MeV/c^2");
    hrm->Draw();
    hrmp->SetLineColor(kRed);
    hrmp->Draw("same");

    c2->SaveAs("ResonantMassPostRemove.pdf");
    c2->SaveAs("ResonantMassPostRemove.root");

    TCanvas *c4 = new TCanvas("c4","",1200,800);
    hm->SetLineColor(kBlue);
    hm->GetXaxis()->SetTitle("Mass/MeV/c^2");
    hm->Draw();

    c4->SaveAs("BMassAR.pdf");
    c4->SaveAs("BMassAR.root");

    // TCanvas *c3 = new TCanvas("c3","",1200,800);
    //hxy->SetStats(0);
    //hxy->GetXaxis()->SetTitle("Probability of Kaon");
    // hxy->GetYaxis()->SetTitle("Probability of Pion");
    //gPad->SetLogz(1);
    //hxy->Draw("colz");
    // TLine *line = new TLine(-0.1,0.8,1.1,0.8);
    //line->SetLineColor(kBlack);
    //line->Draw("Same");
    //TLine *line = new TLine(0.8,-0.1,0.8,1.1);
    //line->SetLineColor(kBlack);
    //line->Draw("Same");
    //TLine *line = new TLine(-0.1,1,1.1,1);
    //line->SetLineColor(kBlack);
    //line->Draw("Same");
    //TLine *line = new TLine(1,-0.1,1,1.1);
    //line->SetLineColor(kBlack);
    //line->Draw("Same");

    //c3->SaveAs("ProbabilityCheckcut.pdf");
    //c3->SaveAs("ProbabilityCheckcut.root");
    //c3->SaveAs("ProbabilityCheckcut.png");
}
