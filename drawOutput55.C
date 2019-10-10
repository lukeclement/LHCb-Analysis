// Exponential background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0]*exp(-par[1]*x[0])+par[2];
}
// Gauss function
Double_t gauss(Double_t *x, Double_t *par) {
  // return (0.5*par[0]*par[1]/TMath::Pi()) /
  //TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
  // + .25*par[1]*par[1]);

 return par[2]/(pow(2*TMath::Pi()*par[0]*par[0],0.5)) *exp(-pow(x[0]-par[1],2)/(2*par[0]*par[0]));
}

Double_t gauss2(Double_t *x, Double_t *par){
  return par[2]/(pow(2*TMath::Pi()*par[0]*par[0],0.5)) *exp(-pow(x[0]-par[1],2)/(2*par[0]*par[0]));
  }

Double_t crystal(Double_t *x, Double_t *par){
  Double_t alpha=par[0];
      
  Double_t n=par[1];
  Double_t sigma=par[2];
  Double_t mean=par[3];
  Double_t A=pow(n/alpha,n)*exp(-(alpha*alpha)/2);
  Double_t B=n/alpha - alpha;
  Double_t C=n/alpha * 1/(n-1) * exp(-alpha*alpha/2);
  Double_t D=pow(TMath::Pi()/2,0.5)*(1+TMath::Erf(alpha/pow(2,0.5)));
  Double_t N=1/(sigma*(C+D));

  if(((x[0]-mean)/sigma)>-alpha){
    return par[4]*N*exp(-pow(x[0]-mean,2)/(2*sigma*sigma));
  }else{
    return par[4]*N*A*pow(B-((x[0]-mean)/sigma),-n);
  }

}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x,par)+crystal(x,&par[3])+gauss2(x,&par[8]);
}

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
    TH1F *hm = (TH1F*)f->Get("h_M1");
    TH1F *hrm= (TH1F*)f->Get("h_M2");
    TH1F *hrmp= (TH1F*)f->Get("h_RMP");

    TH2F *hxy= (TH2F*)f->Get("h_PKP");



    TCanvas *c4 = new TCanvas("c4","",1200,800);

    hm->GetYaxis()->SetRangeUser(0,1200);
    hm->SetLineColor(kBlue);
    hm->GetXaxis()->SetTitle("Mass/MeV/c^2");
    hm->Draw();
    //hrm->SetLineColor(kRed);
    //hrm->Draw("same");
    
 
    TF1 *bestFit = new TF1("bestFit",fitFunction,5100,5500,11);
    bestFit->SetNpx(500);
    bestFit->SetLineWidth(4);
    bestFit->SetLineColor(kMagenta);
    //0;Exp 1st term, 1;exp mean, 2;exp offset, 
    //3;alpha,4;n,5;sigma,6;mean,7;magnitude,
    //8;sigma,9;mean,10;magnitude
    bestFit->SetParameters(22,1,1,10,10,10,10,200,40,5000,1);
    //bestFit->SetParameters(1,0.5,1,1,1,1);
    // Set start values for some parameters

    bestFit->SetParameter(3,0.5);
    bestFit->SetParameter(4,0.5);
    bestFit->SetParameter(5,40);
    bestFit->SetParameter(6,5280);

    bestFit->SetParameter(9,5050);
    bestFit->SetParameter(8,92);

    bestFit->SetParLimits(1,0,30);
    bestFit->SetParLimits(0,0,30);

    //bestFit->SetParLimits(3,0.8,1.1);
    //bestFit->SetParLimits(4,3.8,4.2);
    bestFit->SetParLimits(3,0,2);
    bestFit->SetParLimits(4,1,5);
    bestFit->SetParLimits(5,0,50);
    bestFit->SetParLimits(6,5200,5350);

    bestFit->SetParLimits(8,10,70);
    bestFit->SetParLimits(9,4900,5200);
    bestFit->SetParLimits(10,0,1e8);


    hm->Fit("bestFit","V+", "ep");
    
    TF1 *backFcn = new TF1("backFcn",background,5100,5500,3);
    backFcn->SetLineColor(kGreen);

    //THIS IS THE TARGET
    TF1 *signalFcn = new TF1("signalFcn",crystal,5100,5500,5);
    signalFcn->SetLineColor(kBlue);
    signalFcn->SetNpx(500);
    //TARGET IS ABOVE

    TF1 *signalFcn2 = new TF1("signalFcn2",gauss2,5100,5500,3);
    signalFcn2->SetLineColor(kRed);
    signalFcn2->SetNpx(500);
    Double_t par[10];

    bestFit->GetParameters(par);
    
    backFcn->SetParameters(par);
    backFcn->Draw("same");

    signalFcn->SetParameters(&par[3]);
    Double_t value=signalFcn->Integral(5100,5500)/ hm->GetBinWidth(1);
    //Double_t error=signalFcn->IntegralError(5100,5500,signalFcn->GetParams(), signalFcn->GetCovarianceMatrix()->GetMatrixArray());
    //error=error/2;
    signalFcn->Draw("same");


    signalFcn2->SetParameters(&par[8]);
    signalFcn2->Draw("same");

    c4->SaveAs("BPlusFixedT.pdf");
    c4->SaveAs("BPlusFixedT.root");

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
    cout<<value<<endl;
    // cout<<error<<endl;
}
