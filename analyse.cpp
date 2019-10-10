#include "Analysis.hpp"

// This is the analysis class, which realises the generic Analysis
// from Analysis.hpp
//
// Look in Analysis.hpp for the event variables available.
class MyAnalysis : public Analysis {
public:
    // Define your histograms here
    TH1F           *h_PX;
    TH1F           *h_PY;
    TH1F           *h_PZ;
    TH1F           *h_M1;
    TH1F           *h_M2;
    TH1F           *h_PK;   
    TH1F           *h_PP;
    TH1F           *h_RMK;
    TH1F           *h_RMP;   

    TH2F           *h_DPNE;
    TH2F           *h_DPPE;
    TH2F           *h_DPN;
    TH2F           *h_DPP;
    TH2F           *h_prb;

    double PX;
    double PY;
    double PZ;
    double E1;
    double E2;
    double E3;
    double background;
    
    int BPlus;
    int BMinus;

    std::vector<double> masses;

    void     BookHistos();

    Bool_t   Cut();
    void     Execute();

    int Identity(double,double,double,bool,bool,bool);
    double FindMass(double,double,double,double,double,double,double,double);
    int FindCharge(bool,double,double,double);
};

void MyAnalysis::BookHistos()
{
    // This function is only called once at the start of the program.
    // Book your histograms here. The format is object_name,
    // histogram_name, number_of_bins, minimum, maximum For a 2D
    // histogram, use TH2F with first the number of bins and limits
    // for the x axis and then for the y axis
    //
    // push_back() adds the histograms to a vector v_Histos.  This
    // will take care of writing out histograms in
    // Analysis::SaveHistos
    v_Histos.push_back( h_PX   = new TH1F("h_PX",  "Momentum in x", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PY   = new TH1F("h_PY",  "Momentum in y", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PZ   = new TH1F("h_PZ",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_M1   = new TH1F("h_M1",  "Mass of B+ Meson", 200, 5100, 5500) );
    v_Histos.push_back( h_M2   = new TH1F("h_M2",  "Mass of B- Meson", 200, 5100, 5500) );
    v_Histos.push_back( h_RMK  = new TH1F("h_RMK",  "Mass of Resonances(PK)", 1000, 0, 5000) );
    v_Histos.push_back( h_RMP  = new TH1F("h_RMP",  "Mass of Resonances(PP)", 1000, 0, 5000) );
    v_Histos.push_back( h_PP   = new TH1F("h_PP",  "Pion prob", 100, -0.1, 1.1) );
    v_Histos.push_back( h_PK   = new TH1F("h_PK",  "Kaon prob", 100, -0.1, 1.1) );
    v_Histos.push_back( h_DPN  = new TH2F("h_DPN", "Dalitz plot -", 30, 0, 30, 30,0,30) );
    v_Histos.push_back( h_DPP  = new TH2F("h_DPP", "Dalitz plot +", 30, 0, 30, 30,0,30) );
    v_Histos.push_back( h_DPNE = new TH2F("h_DPNE", "Dalitz plot -EXTRA", 30, 0, 30, 30,0,30) );
    v_Histos.push_back( h_DPPE = new TH2F("h_DPPE", "Dalitz plot +EXTRA", 30, 0, 30, 30,0,30) );
    v_Histos.push_back( h_prb  = new TH2F("h_prb", "Kaon and Pion probabilities", 200, -0.1, 1.1, 200,-0.1,1.1) );
    h_DPP->Sumw2();

}

int MyAnalysis::Identity(double probK, double probPi, double Pcharge,  bool K, bool Pm, bool Pp){
    //Function to perform particle identification
    double Kmin = 0.8;
    double Pimin= 0.2;
    //Override
    //return 1;

    //See if K
    if(probK>Kmin){
      if(!K && Pcharge!=0){
	return 1;
      }else{
	return 0;
      }
      //See if pi
    }else if(probPi>Pimin && probK<0.4){
      if(Pcharge==-1){
	if(!Pm){
	  return 2;
	}else{
	  return 0;
	}
      }else{
	if(!Pp){
	  return 3;
	}else{
	  return 0;
	}
      }
    }else{
      return 0;
    }
}

int MyAnalysis::FindCharge(bool target, double h1c, double h2c, double h3c){
  int tCharge=0;
  if(target){
    tCharge=1;
  }else{
    tCharge=-1;
  }
  if(h1c==tCharge){
    return 1;
  }else if(h2c==tCharge){
    return 2;
  }else if(h3c==tCharge){
    return 3;
  }

  return 0;
}

Bool_t MyAnalysis::Cut()
{
    // This function is called for every event from the Execute
    // function to define whether or not to accept this event.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.

    //Some useful constants
    double Kpm = 493.677;
    double K0  = 497.648;
    double Ppm = 139.571;

    //Selecting ony KPP events
  bool kaon      = false;
  bool pionMinus = false;
  bool pionPlus  = false;
  double muonMin = 0.8;
  masses.clear();
  if(H1_isMuon>muonMin||H2_isMuon>muonMin||H3_isMuon>muonMin){
    return false;
  }
  switch(Identity(H1_ProbK, H1_ProbPi, H1_Charge,kaon,pionMinus,pionPlus)){
  case 0:
    return false;break;
  case 1:
    kaon = true;
    masses.push_back(Kpm);break;
  case 2:
    pionMinus=true;
    masses.push_back(Ppm);break;
  case 3:
    pionPlus=true;
    masses.push_back(Ppm);break;
  }switch(Identity(H2_ProbK, H2_ProbPi, H2_Charge,kaon,pionMinus,pionPlus)){
  case 0:
    return false;break;
  case 1:
    kaon = true;
    masses.push_back(Kpm);break;
  case 2:
    pionMinus=true;
    masses.push_back(Ppm);break;
  case 3:
    pionPlus=true;
    masses.push_back(Ppm);break;
  }switch(Identity(H3_ProbK, H3_ProbPi, H3_Charge,kaon,pionMinus,pionPlus)){
  case 0:
    return false;break;
  case 1:
    kaon = true;
    masses.push_back(Kpm);break;
  case 2:
    pionMinus=true;
    masses.push_back(Ppm);break;
  case 3:
    pionPlus=true;
    masses.push_back(Ppm);break;
  }

  //Removing D mesons
  //Charges must sum to 0!!!
    double totalCharge=0;
    PX=0;
    PY=0;
    PZ=0;
    E1=0;
    E2=0;
    E3=0;
    bool KHit=false;
    int K=4;

    double M12 = pow(FindMass(H1_PX, H2_PX, H1_PY, H2_PY, H1_PZ, H2_PZ, masses[0],masses[1]),0.5)*1000;
    double M13 = pow(FindMass(H1_PX, H3_PX, H1_PY, H3_PY, H1_PZ, H3_PZ, masses[0],masses[2]),0.5)*1000;
    double M23 = pow(FindMass(H2_PX, H3_PX, H2_PY, H3_PY, H2_PZ, H3_PZ, masses[1],masses[2]),0.5)*1000;
    for(int i=0;i<3;i++){
      if(masses[i]==Kpm){
	K=i;
      }
    }

    //Limits to remove D
    double maxDMass=1895;
    double minDMass=1835;
    if(M12<maxDMass && M12>minDMass){
      return false;
    } if(M13<maxDMass && M13>minDMass){
      return false;
    } if(M23<maxDMass && M23>minDMass){
      return false;
    }
    //Add more limits to remove other mesons //NOT IMPLIMENTED
    double f0Min = 920;
    double f0Max = 1028;
    
    double rhoMin= 693;
    double rhoMax= 812;
    
    double PK1Min= 1465;
    double PK1Max= 1738;
    
    double PK2Min= 818;
    double PK2Max= 977;
    /*if(KHit){
      if(resonantMass<PK1Max && resonantMass>PK1Min){
	//return false;
      }
      if(resonantMass<PK2Max && resonantMass>PK2Min){
	//return false;
      }
    }else{
      if(resonantMass<f0Max && resonantMass>f0Min){
	//return false;
      }
      if(resonantMass<rhoMax && resonantMass>rhoMin){
	//return false;
      }
      }*/
    PX = H1_PX + H2_PX + H3_PX;
    PY = H1_PY + H2_PY + H3_PY;
    PZ = H1_PZ + H2_PZ + H3_PZ;
    
    E1 = pow(pow(masses[0],2) + pow(H1_PX,2) + pow(H1_PY,2) + pow(H1_PZ,2) ,0.5);
    E2 = pow(pow(masses[1],2) + pow(H2_PX,2) + pow(H2_PY,2) + pow(H2_PZ,2) ,0.5);
    E3 = pow(pow(masses[2],2) + pow(H3_PX,2) + pow(H3_PY,2) + pow(H3_PZ,2) ,0.5);
    
    double BMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
    if(BMass<5100 ||  BMass>5600){
      return false;
    }
    
    return true;
}

double MyAnalysis::FindMass(double PXA, double PXB, double PYA, double PYB, double PZA, double PZB, double MA, double MB){
  double px=PXA + PXB;
  double py=PYA + PYB;
  double pz=PZA + PZB;

  double ea=pow(pow(MA,2) + pow(PXA,2) + pow(PYA,2) + pow(PZA,2),0.5);
  double eb=pow(pow(MB,2) + pow(PXB,2) + pow(PYB,2) + pow(PZB,2),0.5);
  
  double M=pow(pow(ea+eb,2)-pow(px,2)-pow(py,2)-pow(pz,2),0.5)/1000;
  return M*M;
}

void MyAnalysis::Execute()
{
    BMinus = 0;
    BPlus  = 0;
    // This method gets called on every event.
    // In this example the momentum components are filled into histograms.

    // Call the Cut function to decide whether to plot this event or not
    // it returns if the cut function returns false
    if ( !Cut() )
	return;

    
    h_prb->Fill(H1_ProbPi, H1_ProbK);
    h_prb->Fill(H2_ProbPi, H2_ProbK);
    h_prb->Fill(H3_ProbPi, H3_ProbK);
    // Fill your histograms below
    PX = H1_PX + H2_PX + H3_PX;
    PY = H1_PY + H2_PY + H3_PY;
    PZ = H1_PZ + H2_PZ + H3_PZ;
    
    E1 = pow(pow(masses[0],2) + pow(H1_PX,2) + pow(H1_PY,2) + pow(H1_PZ,2) ,0.5);
    E2 = pow(pow(masses[1],2) + pow(H2_PX,2) + pow(H2_PY,2) + pow(H2_PZ,2) ,0.5);
    E3 = pow(pow(masses[2],2) + pow(H3_PX,2) + pow(H3_PY,2) + pow(H3_PZ,2) ,0.5);
    //std::cout << pow(pow(E1+E2+E3,2)+pow(P1+P2+P3,2),0.5) << std::endl;
    double superMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
    //std::cout<<"s1"<<std::endl;
    if(Identity(H1_ProbK, H1_ProbPi, H1_Charge,false,false,false)==1){
      if(H1_Charge==1){
	BPlus++;
	//superMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
	h_M1->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
      }else{
	BMinus++;
	//superMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
	h_M2->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
      }
    }else if(Identity(H2_ProbK, H2_ProbPi, H2_Charge,false,false,false)==1){
      if(H2_Charge==1){
	BPlus++;
	//superMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
	h_M1->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
      }else{
	BMinus++;
	//superMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
	h_M2->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
      }
    }else if(Identity(H3_ProbK, H3_ProbPi, H3_Charge,false,false,false)==1){
      if(H3_Charge==1){
	BPlus++;
	//superMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
	h_M1->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
      }else{
	BMinus++;
	//superMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
	h_M2->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
      }
    }

    bool onePi    = false;
    bool twoPi    = false;
    bool thrPi    = false;
    double piMass = 139.571;
    double kMass = 493.677;
    // std::cout<<"s2"<<std::endl;
    if(superMass>=5250 && superMass<=5350){
      if(BPlus==1){
	//Find the -
	switch(FindCharge(false,H1_Charge,H2_Charge,H3_Charge)){
	case 1:{
	  double M12=FindMass(H1_PX, H2_PX, H1_PY, H2_PY, H1_PZ, H2_PZ, masses[0],masses[1]);
	  double M13=FindMass(H1_PX, H3_PX, H1_PY, H3_PY, H1_PZ, H3_PZ, masses[0],masses[2]);
	  if(masses[0]==masses[1]){
	    h_DPP->Fill(M12,M13);
	  }else{
	    h_DPP->Fill(M13,M12);
	  }
	  break;}
	case 2:{
	  double M21=FindMass(H2_PX, H1_PX, H2_PY, H1_PY, H2_PZ, H1_PZ, masses[1],masses[0]);
	  double M23=FindMass(H2_PX, H3_PX, H2_PY, H3_PY, H2_PZ, H3_PZ, masses[1],masses[2]);
	  if(masses[0]==masses[1]){
	    h_DPP->Fill(M21,M23);
	  }else{
	    h_DPP->Fill(M23,M21);
	  }
	  break;}
	case 3:{
	  double M31=FindMass(H3_PX, H1_PX, H3_PY, H1_PY, H3_PZ, H1_PZ, masses[2],masses[0]);
	  double M32=FindMass(H3_PX, H2_PX, H3_PY, H2_PY, H3_PZ, H2_PZ, masses[2],masses[1]);
	  if(masses[2]==masses[0]){
	    h_DPP->Fill(M31,M32);
	  }else{
	    h_DPP->Fill(M32,M31);
	  }
	  break;}
	}
      }else{
	//Find the +
	switch(FindCharge(true,H1_Charge,H2_Charge,H3_Charge)){
	case 1:{
	  double M122=FindMass(H1_PX, H2_PX, H1_PY, H2_PY, H1_PZ, H2_PZ, masses[0],masses[1]);
	  double M132=FindMass(H1_PX, H3_PX, H1_PY, H3_PY, H1_PZ, H3_PZ, masses[0],masses[2]);
	  if(masses[0]==masses[1]){
	    h_DPN->Fill(M122,M132);
	  }else{
	    h_DPN->Fill(M132,M122);
	  }
	  break;}
	case 2:{
	  double M212=FindMass(H2_PX, H1_PX, H2_PY, H1_PY, H2_PZ, H1_PZ, masses[1],masses[0]);
	  double M232=FindMass(H2_PX, H3_PX, H2_PY, H3_PY, H2_PZ, H3_PZ, masses[1],masses[2]);
	  if(masses[0]==masses[1]){
	    h_DPN->Fill(M212,M232);
	  }else{
	    h_DPN->Fill(M232,M212);
	  }
	  break;}
	case 3:{
	  double M312=FindMass(H3_PX, H1_PX, H3_PY, H1_PY, H3_PZ, H1_PZ, masses[2],masses[0]);
	  double M322=FindMass(H3_PX, H2_PX, H3_PY, H2_PY, H3_PZ, H2_PZ, masses[2],masses[1]);
	  if(masses[2]==masses[0]){
	    h_DPN->Fill(M312,M322);
	  }else{
	    h_DPN->Fill(M322,M312);
	  }
	  break;}
	}
      }
    }
    
    //   std::cout<<"s3"<<std::endl;
    if(superMass>=5350&&superMass<=5450){
      //     std::cout<<"s3.1"<<std::endl;

      //     std::cout<<BPlus<<std::endl;
      if(BPlus==1){
	//Find the -
	switch(FindCharge(false,H1_Charge,H2_Charge,H3_Charge)){
	case 1:{
	  //	  std::cout<<"s3.31"<<std::endl;
	  double M12E=FindMass(H1_PX, H2_PX, H1_PY, H2_PY, H1_PZ, H2_PZ, masses[0],masses[1]);
	  double M13E=FindMass(H1_PX, H3_PX, H1_PY, H3_PY, H1_PZ, H3_PZ, masses[0],masses[2]);
	  //	  std::cout<<"s3.31.1"<<std::endl;
	  if(masses[0]==masses[1]){
	    //	    std::cout<<"s3.31.2T"<<std::endl;
	    h_DPPE->Fill(M12E,M13E);
	    //	    std::cout<<"s3.31.3T"<<std::endl;
	  }else{
	    //   std::cout<<"s3.31.2F"<<std::endl;
	    h_DPPE->Fill(M13E,M12E);
	    //   std::cout<<"s3.31.3F"<<std::endl;
	  }
	  // std::cout<<"s3.31.5"<<std::endl;
	  break;}
	case 2:{
	  // std::cout<<"s3.32"<<std::endl;
	  double M21E=FindMass(H2_PX, H1_PX, H2_PY, H1_PY, H2_PZ, H1_PZ, masses[1],masses[0]);
	  double M23E=FindMass(H2_PX, H3_PX, H2_PY, H3_PY, H2_PZ, H3_PZ, masses[1],masses[2]);
	  if(masses[0]==masses[1]){
	    h_DPPE->Fill(M21E,M23E);
	  }else{
	    h_DPPE->Fill(M23E,M21E);
	  }
	  break;}
	case 3:{
	  //std::cout<<"s3.33"<<std::endl;
	  double M31E=FindMass(H3_PX, H1_PX, H3_PY, H1_PY, H3_PZ, H1_PZ, masses[2],masses[0]);
	  double M32E=FindMass(H3_PX, H2_PX, H3_PY, H2_PY, H3_PZ, H2_PZ, masses[2],masses[1]);
	  if(masses[2]==masses[0]){
	    h_DPPE->Fill(M31E,M32E);
	  }else{
	    h_DPPE->Fill(M32E,M31E);
	  }
	  break;}
	}
      }else{
	//Find the +
	//std::cout<<"s3.2"<<std::endl;
	switch(FindCharge(true,H1_Charge,H2_Charge,H3_Charge)){
	case 1:{
	  double M122E=FindMass(H1_PX, H2_PX, H1_PY, H2_PY, H1_PZ, H2_PZ, masses[0],masses[1]);
	  double M132E=FindMass(H1_PX, H3_PX, H1_PY, H3_PY, H1_PZ, H3_PZ, masses[0],masses[2]);
	  if(masses[0]==masses[1]){
	    h_DPNE->Fill(M122E,M132E);
	  }else{
	    h_DPNE->Fill(M132E,M122E);
	  }
	  break;}
	case 2:{
	  double M212E=FindMass(H2_PX, H1_PX, H2_PY, H1_PY, H2_PZ, H1_PZ, masses[1],masses[0]);
	  double M232E=FindMass(H2_PX, H3_PX, H2_PY, H3_PY, H2_PZ, H3_PZ, masses[1],masses[2]);
	  if(masses[0]==masses[1]){
	    h_DPNE->Fill(M212E,M232E);
	  }else{
	    h_DPNE->Fill(M232E,M212E);
	  }
	  break;}
	case 3:{
	  double M312E=FindMass(H3_PX, H1_PX, H3_PY, H1_PY, H3_PZ, H1_PZ, masses[2],masses[0]);
	  double M322E=FindMass(H3_PX, H2_PX, H3_PY, H2_PY, H3_PZ, H2_PZ, masses[2],masses[1]);
	  if(masses[2]==masses[0]){
	    h_DPNE->Fill(M312E,M322E);
	  }else{
	    h_DPNE->Fill(M322E,M312E);
	  }
	  break;}
	}
      }
    }
    //std::cout<<"s4"<<std::endl;
   
    // 2D histogram of Dalitz plot
    //h_DP->Fill( H2_ProbK , H2_ProbPi);
    //h_DP->Fill( H3_ProbK , H3_ProbPi);
}


// The main function just calls the generic AnalysisMain function
// with the MyAnalysis class
//
// Normally you don't need to change this
int main(int argc, char* argv[])
{
    MyAnalysis* ana = new MyAnalysis();
    int res = ana->AnalysisMain(argc, argv);
    return res;
}
