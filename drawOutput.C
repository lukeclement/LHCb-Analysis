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
//Cruijff function (https://journals.aps.org/prd/pdf/10.1103/PhysRevD.82.051101)
Double_t cru(Double_t *x, Double_t *par){
  Double_t alpha=par[0];
  
  Double_t sigma=par[2];
  Double_t mean=par[3];
  Double_t mag=par[4];
  
  return mag*exp(-pow(x-mean,2)/(2*sigma*sigma + alpha*pow(x-mean,2) ) );
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
    TH1F *hPP = (TH1F*)f->Get("h_PP");//Pion probability
    TH1F *hPK = (TH1F*)f->Get("h_PK");//Kaon probability
    TH1F *hm1 = (TH1F*)f->Get("h_M1");//Mass of B+
    TH1F *hm2 = (TH1F*)f->Get("h_M2");//Mass of B-
    TH1F *hrmp= (TH1F*)f->Get("h_RMP");//Resonance masses of Pion/Pion
    TH1F *hrmk= (TH1F*)f->Get("h_RMK");//Resonance masses of Kaon/Pion

    TH2F *hxy  = (TH2F*)f->Get("h_DPP");//Dalitz plot of B+
    TH2F *hxyM = (TH2F*)f->Get("h_DPN");//Dalitz plot of B-
    TH2F *hxyE = (TH2F*)f->Get("h_DPPE");//Dalitz plot of B+
    TH2F *hxyME= (TH2F*)f->Get("h_DPNE");//Dalitz plot of B-
    TH2F *hkp  = (TH2F*)f->Get("h_prb");//Kaon/Pion 2D probability

    //Dalitz plot
    //B+ Plot
    TCanvas *c3 = new TCanvas("c3","",1200,800);
    hxy->SetStats(0);
    hxy->GetXaxis()->SetTitle("M^2 P+P-");
    hxy->GetYaxis()->SetTitle("M^2 K+P-");
    hxy->Draw("COLZ");

    c3->SaveAs("Dalitz+TestMK1.pdf");
    c3->SaveAs("Dalitz+TestMK1.root");
    c3->SaveAs("Dalitz+TestMK1.png");
    
    //B- Plot
    TCanvas *c1 = new TCanvas("c1","",1200,800);
    hxyM->SetStats(0);
    hxyM->GetXaxis()->SetTitle("M^2 P+P-");
    hxyM->GetYaxis()->SetTitle("M^2 K-P+");
    hxyM->Draw("COLZ");

    c1->SaveAs("Dalitz-TestMK1.pdf");
    c1->SaveAs("Dalitz-TestMK1.root");
    c1->SaveAs("Dalitz-TestMK1.png");
    //B+ Plot
    TCanvas *cE3 = new TCanvas("cE3","",1200,800);
    hxyE->SetStats(0);
    hxyE->GetXaxis()->SetTitle("M^2 P+P-");
    hxyE->GetYaxis()->SetTitle("M^2 K+P-");
    hxyE->Draw("COLZ");

    cE3->SaveAs("Dalitz+TestMK1E.pdf");
    cE3->SaveAs("Dalitz+TestMK1E.root");
    cE3->SaveAs("Dalitz+TestMK1E.png");
    
    //B- Plot
    TCanvas *cE1 = new TCanvas("cE1","",1200,800);
    hxyME->SetStats(0);
    hxyME->GetXaxis()->SetTitle("M^2 P+P-");
    hxyME->GetYaxis()->SetTitle("M^2 K-P+");
    hxyME->Draw("COLZ");

    cE1 ->SaveAs("Dalitz-TestMK1E.pdf");
    cE1->SaveAs("Dalitz-TestMK1E.root");
    cE1->SaveAs("Dalitz-TestMK1E.png");

    //DIFFERENCES
    //B+ Plot
    TCanvas *c10 = new TCanvas("c10","",1200,800);
    hxy->Add(hxyE,-1);
    TH2F *hxyA = new TH2F("hxyA","Dalitz plot for B+",30,0,30,30,0,30);
    for(int i=0;i<30;i++){
      for(int j=0;j<30;j++){
	Double_t v=hxy->GetBinContent(hxy->GetBin(i,j));
	if(v>=0){
	  hxyA->SetBinContent(i,j,v);
   	}
      }
    }

    hxyA->SetStats(0)
    hxyA->GetXaxis()->SetTitle("M^2 P+P-");
    hxyA->GetYaxis()->SetTitle("M^2 K+P-");
    hxyA->Draw("COLZ");

    c10->SaveAs("Dalitz+TestMK2.pdf");
    c10->SaveAs("Dalitz+TestMK2.root");
    c10->SaveAs("Dalitz+TestMK2.png");
    
    //B- Plot
    TCanvas *c11 = new TCanvas("c11","",1200,800);
    hxyM->Add(hxyME,-1);
    TH2F *hxyB = new TH2F("hxyB","Dalitz plot for B-",30,0,30,30,0,30);

    for(int i=0;i<30;i++){
      for(int j=0;j<30;j++){
	Double_t v=hxyM->GetBinContent(hxyM->GetBin(i,j));
	if(v>=0){
	  hxyB->SetBinContent(i,j,v);
   	}
      }
    }

    hxyB->SetStats(0);
    hxyB->GetXaxis()->SetTitle("M^2 P+P-");
    hxyB->GetYaxis()->SetTitle("M^2 K-P+");
    hxyB->Draw("COLZ");

    c11 ->SaveAs("Dalitz-TestMK2.pdf");
    c11->SaveAs("Dalitz-TestMK2.root");
    c11->SaveAs("Dalitz-TestMK2.png");



















    //Mass of B+ with fitting
    TCanvas *c2 = new TCanvas("c2","",1200,800);
    hm1->GetYaxis()->SetRangeUser(0,1200);
    hm1->SetLineColor(kBlue);
    hm1->GetXaxis()->SetTitle("Mass/MeV/c^2");
    hm1->Draw();
    
    //Fitting lines to data
    TF1 *bestFit = new TF1("bestFit",fitFunction,5100,5500,11);
    bestFit->SetNpx(500);
    bestFit->SetLineWidth(4);
    bestFit->SetLineColor(kMagenta);
    //0;Exp 1st term, 1;exp mean, 2;exp offset, 
    //3;alpha,4;n,5;sigma,6;mean,7;magnitude,
    //8;sigma,9;mean,10;magnitude
    bestFit->SetParameters(22,1,1,10,10,10,10,200,40,5000,1);
    //Set start values for some parameters
    
    //Crystal ball function
    bestFit->SetParameter(3,0.5);
    bestFit->SetParameter(4,0.5);
    bestFit->SetParameter(5,40);
    bestFit->SetParameter(6,5280);
    
    //Gauss for 4-body masses
    bestFit->SetParameter(9,5050);
    bestFit->SetParameter(8,92);
    
    //Setting limits
    //Exponential background fit limits
    bestFit->SetParLimits(1,0,30);
    bestFit->SetParLimits(0,0,30);
 
    //Crystal ball limits
    bestFit->SetParLimits(3,0,2);
    bestFit->SetParLimits(4,1,5);
    bestFit->SetParLimits(5,0,50);
    bestFit->SetParLimits(6,5200,5350);
    bestFit->SetParLimits(7,0,1e8);
    
    //Gauss limits
    bestFit->SetParLimits(8,10,70);
    bestFit->SetParLimits(9,4900,5100);
    bestFit->SetParLimits(10,0,1e8);

    //Fitting!
    hm1->Fit("bestFit","V+", "ep");
    
    //Background function
    TF1 *backFcn = new TF1("backFcn",background,5100,5500,3);
    backFcn->SetLineColor(kGreen);

    //THIS IS THE TARGET
    TF1 *signalFcn = new TF1("signalFcn",crystal,5100,5500,5);
    signalFcn->SetLineColor(kBlue);
    signalFcn->SetNpx(500);
    //TARGET IS ABOVE

    //Gauss signal function
    TF1 *signalFcn2 = new TF1("signalFcn2",gauss2,5100,5500,3);
    signalFcn2->SetLineColor(kRed);
    signalFcn2->SetNpx(500);
    Double_t parx[10];

    //Getting and sorting parameters
    bestFit->GetParameters(parx);
    
    backFcn->SetParameters(parx);
    backFcn->Draw("same");

    signalFcn->SetParameters(&parx[3]);
    Double_t valueP=signalFcn->Integral(5100,5500)/ hm1->GetBinWidth(1);
    //Double_t error=signalFcn->IntegralError(5100,5500,signalFcn->GetParams(), signalFcn->GetCovarianceMatrix()->GetMatrixArray());
    //error=error/2;
    signalFcn->Draw("same");

    signalFcn2->SetParameters(&parx[8]);
    signalFcn2->Draw("same");

    c2->SaveAs("B+Mass.pdf");
    c2->SaveAs("B+Mass.root");
    c2->SaveAs("B+Mass.png");

    //
    //

    //Mass of B- with fitting
    TCanvas *c4 = new TCanvas("c2","",1200,800);
    hm2->GetYaxis()->SetRangeUser(0,1200);
    hm2->SetLineColor(kBlue);
    hm2->GetXaxis()->SetTitle("Mass/MeV/c^2");
    hm2->Draw();
    
    //Fitting lines to data
    TF1 *bestFit2 = new TF1("bestFit2",fitFunction,5100,5500,11);
    bestFit2->SetNpx(500);
    bestFit2->SetLineWidth(4);
    bestFit2->SetLineColor(kMagenta);
    //0;Exp 1st term, 1;exp mean, 2;exp offset, 
    //3;alpha,4;n,5;sigma,6;mean,7;magnitude,
    //8;sigma,9;mean,10;magnitude
    bestFit2->SetParameters(22,1,1,10,10,10,10,200,40,5000,1);
    //Set start values for some parameters
    
    //Crystal ball function
    bestFit2->SetParameter(3,0.5);
    bestFit2->SetParameter(4,0.5);
    bestFit2->SetParameter(5,40);
    bestFit2->SetParameter(6,5280);
    
    //Gauss for 4-body masses
    bestFit2->SetParameter(9,5050);
    bestFit2->SetParameter(8,92);
    
    //Setting limits
    //Exponential background fit limits
    bestFit2->SetParLimits(1,0,30);
    bestFit2->SetParLimits(0,0,30);
 
    //Crystal ball limits
    bestFit2->SetParLimits(3,0,2);
    bestFit2->SetParLimits(4,1,5);
    bestFit2->SetParLimits(5,0,50);
    bestFit2->SetParLimits(6,5200,5350);
    bestFit2->SetParLimits(7,0,1e8);
    
    //Gauss limits
    bestFit2->SetParLimits(8,10,70);
    bestFit2->SetParLimits(9,4900,5100);
    bestFit2->SetParLimits(10,0,1e8);

    //Fitting!
    hm2->Fit("bestFit2","V+", "ep");
    
    //Background function
    TF1 *backFcn2 = new TF1("backFcn2",background,5100,5500,3);
    backFcn2->SetLineColor(kGreen);

    //THIS IS THE TARGET
    TF1 *signalFcn9 = new TF1("signalFcn2",crystal,5100,5500,5);
    signalFcn9->SetLineColor(kBlue);
    signalFcn9->SetNpx(500);
    //TARGET IS ABOVE

    //Gauss signal function
    TF1 *signalFcn22 = new TF1("signalFcn22",gauss2,5100,5500,3);
    signalFcn22->SetLineColor(kRed);
    signalFcn22->SetNpx(500);
    Double_t par2[10];

    //Getting and sorting parameters
    bestFit2->GetParameters(par2);
    
    backFcn2->SetParameters(par2);
    backFcn2->Draw("same");

    signalFcn9->SetParameters(&par2[3]);
    Double_t valueM=signalFcn9->Integral(5100,5500)/ hm2->GetBinWidth(1);
    //Double_t error=signalFcn->IntegralError(5100,5500,signalFcn->GetParams(), signalFcn->GetCovarianceMatrix()->GetMatrixArray());
    //error=error/2;
    signalFcn9->Draw("same");

    signalFcn22->SetParameters(&par2[8]);
    signalFcn22->Draw("same");

    c4->SaveAs("B-Mass.pdf");
    c4->SaveAs("B-Mass.root");
    c4->SaveAs("B-Mass.png");

    //Kaon/Pion Probabilities
    TCanvas *c5 = new TCanvas("c5","",1200,800);
    hkp->SetStats(0);
    hkp->GetXaxis()->SetTitle("Pion Probability");
    hkp->GetYaxis()->SetTitle("Kaon Probability");
    gPad->SetLogz(1);
    hkp->Draw("colz");

    c5->SaveAs("KaonvsPionP.pdf");
    c5->SaveAs("KaonvsPionP.root");
    c5->SaveAs("KaonvsPionP.png");
    

    
    
    TCanvas *c6 = new TCanvas("c6","",1200,800);
    hA = (TH2F*) hxyA -> GetAsymmetry(hxyB);


    TH2F *hAE = new TH2F("hxyA","Asymmetry error",30,0,30,30,0,30);
    for(int i=0;i<30;i++){
      for(int j=0;j<30;j++){
	Double_t A=hA->GetBinContent(hA->GetBin(i,j));
	Double_t N1=hxyA->GetBinContent(hxyA->GetBin(i,j));
	Double_t N2=hxyB->GetBinContent(hxyB->GetBin(i,j));
	hAE->SetBinContent(i,j,pow((1-A*A)/(N1+N2),0.5));
      }
    }


    for(int i=0;i<30;i++){
      for(int j=0;j<30;j++){
	Double_t A=hA->GetBinContent(hA->GetBin(i,j));
	Double_t AErr=hAE->GetBinContent(hAE->GetBin(i,j));
	bool sign=A>0;
	bool inverseSign=false;
	if(sign){
	  inverseSign=A-0*AErr < 0;
	}else{
	  inverseSign=A+0*AErr > 0;
	}
	hA->SetBinContent(i,j,pow(A*A,0.5));
	if(inverseSign||AErr==0){
	  hA->SetBinContent(i,j,0);
	}
      }
    }
    hA->SetStats(0);
    hA->GetXaxis()->SetTitle("M^2_P+P-");
    hA->GetYaxis()->SetTitle("M^2_K+-P-+");
    //gPad->SetLogz(1);
    hA->Draw("colz");

    
    c6->SaveAs("Asymmetry0sigma.pdf");
    c6->SaveAs("Asymmetry0sigma.root");
    c6->SaveAs("Asymmetry0sigma.png");

    
    TCanvas *c20 = new TCanvas("c20","",1200,800);
    hAE->SetStats(0);
    hAE->GetXaxis()->SetTitle("M^2_P+P-");
    hAE->GetYaxis()->SetTitle("M^2_K+-P-+");
    hAE->Draw("colz");

    c20->SaveAs("AsymmetryError.pdf");
    c20->SaveAs("AsymmetryError.root");
    c20->SaveAs("AsymmetryError.png");


    //Declaring N+ and N-
    /*Int_t bin=hA->GetBin(5,5);
    TH2F *nA = new TH2F("nA","Asymmetry without stuff",100,0,30,100,0,30);

    for(int i=0;i<100;i++){
      for(int j=0;j<100;j++){
	Double_t v=hA->GetBinContent(hA->GetBin(i,j));
	Double_t ve=hA->GetBinError(hA->GetBin(i,j));
	//if(v!=0 && ve!=0){
	if(v!=0 && v!=1 && v!=-1 && v / ve >= 1){
	  //nA->Fill((double)i*30/100,(double)j*30/100);
	  nA->SetBinContent(i,j,v);
	  cout<<i<<endl;
	  cout<<v<<endl;
	  cout<<hA->GetBinError(hA->GetBin(i,j))<<endl;
	}
      }
    }
    TCanvas *c7 = new TCanvas("c7","",1200,800);
    nA->SetStats(0);
    nA->GetXaxis()->SetTitle("M^2_P+P-");
    nA->GetYaxis()->SetTitle("M^2_K+-P-+");
    //gPad->SetLogz(1);
    nA->Draw("colz");

    
    c7->SaveAs("Asymmetrynew.pdf");
    c7->SaveAs("Asymmetrynew.root");
    c7->SaveAs("Asymmetrynew.png");
    cout<<bin<<endl;
    */
    cout<<valueP<<"->N+"<<endl;
    cout<<valueM<<"->N-"<<endl;
    Double_t A = (valueM-valueP)/(valueM+valueP);
    Double_t AE= pow((1-A*A)/(valueM+valueP),0.5);
    cout<<A<<"+-"<<AE<<"->A"<<endl;
}
