void drawOutput() {
    //////////////////////
    // Example Root Macro for third year B->hhh Lab
    ////////////////////// 

    // Open the root file that was produced by running the example program
    TFile *f = new TFile("outputPhaseSpace.root");

    // Get pointers to the example histograms that were made 
    TH1F *hx = (TH1F*)f->Get("h_PX");
    TH1F *hy = (TH1F*)f->Get("h_PY");
    TH1F *hz = (TH1F*)f->Get("h_PZ");
    TH1F *hm = (TH1F*)f->Get("h_M");

    // Create a canvas onto which the histograms are plotted and which can be saved
    TCanvas *c1 = new TCanvas("c1","",600,400);
    // Draw the first histogram with a blue line, and an x-axis title
    hm->SetLineColor(kBlue);
    hm->GetXaxis()->SetTitle("Mass [MeV/c^{2}]");
    hm->Draw();
    // Save the canvas as pdf file. Other possible formats include root (for later editing) eps, png.
    c1->SaveAs("Mass.pdf");
    c1->SaveAs("Mass.root");

    // Repeat the above for the 2D histogram
    //TCanvas *c2 = new TCanvas("c2","",600,400);
    //hxy->SetStats(0);                           // remove the statistics box
    //hxy->GetXaxis()->SetTitle("Slope in x");    // add axis titles
    //hxy->GetYaxis()->SetTitle("Slope in y");
    //hxy->Draw("colz");                          // draw with a colour scale
    //c2->SaveAs("hxy.pdf");
}
