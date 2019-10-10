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
    //TH2F *hxy = (TH2F*)f->Get("h_TXTY");
    TH1F *hp = (TH1F*)f->Get("h_pP");
    TH1F *hk = (TH1F*)f->Get("h_pK");

    // Create a canvas onto which the histograms are plotted and which can be saved
    TCanvas *c1 = new TCanvas("c1","",600,400);
    // Draw the first histogram with a blue line, and an x-axis title
    hx->SetLineColor(kBlue);
    hx->GetXaxis()->SetTitle("Momentum [MeV/c^{2}]");
    hx->Draw();

    // Draw the second histogram on the same plot with a red line
    hy->SetLineColor(kRed);
    hy->Draw("same");

    hz->SetLineColor(kRed);
    hz->Draw("same");

    // Save the canvas as pdf file. Other possible formats include root (for later editing) eps, png.
    c1->SaveAs("hx_hy_hz.pdf");

    // Repeat the above for the 2D histogram
    //TCanvas *c2 = new TCanvas("c2","",600,400);
    //hxy->SetStats(0);                           // remove the statistics box
    //hxy->GetXaxis()->SetTitle("Slope in x");    // add axis titles
    //hxy->GetYaxis()->SetTitle("Slope in y");
    //hxy->Draw("colz");                          // draw with a colour scale
    //c2->SaveAs("hxy.pdf");
    
    TCanvas *c3 = new TCanvas("c3","",600,400);
    hp->SetLineColor(kBlue);
    hp->SetStats(0);                           // remove the statistics box
    hp->GetXaxis()->SetTitle("Probability of Pion(Blue)/Kaon(Red)");    // add axis titles
    // hp->GetYaxis()->SetTitle("Counts");
    hp->Draw();
    hk->SetLineColor(kRed);
    hk->Draw("same");
    c3->SaveAs("PionKaonProbs.pdf");
    c3->SaveAs("PionKaonProbs.root");
}
