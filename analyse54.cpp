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
    TH1F           *h_M;
    TH1F           *h_PK;   
    TH1F           *h_PP;
    TH1F           *h_RM;
    TH1F           *h_RMP;   

    TH2F           *h_PKP;

    double PX;
    double PY;
    double PZ;
    double E1;
    double E2;
    double E3;

    std::vector<double> masses;


    void     BookHistos();

    Bool_t   Cut();
    void     Execute();

  int Identity(double,double,double,bool,bool,bool);
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
    v_Histos.push_back( h_M   = new TH1F("h_M",  "Mass of B+- Meson", 1000, 4500, 6300) );
    v_Histos.push_back( h_RM   = new TH1F("h_RM",  "Mass of Resonances(PK)", 1000, 0, 5000) );
v_Histos.push_back( h_RMP   = new TH1F("h_RMP",  "Mass of Resonances(PP)", 1000, 0, 5000) );
    v_Histos.push_back( h_PP   = new TH1F("h_PP",  "Pion prob", 100, -0.1, 1.1) );
    v_Histos.push_back( h_PK   = new TH1F("h_PK",  "Kaon prob", 100, -0.1, 1.1) );
    v_Histos.push_back( h_PKP = new TH2F("h_PKP","ProbK vs ProbP", 500, -0.1,1.1, 500,-0.1, 1.1) );


}

int MyAnalysis::Identity(double probK, double probPi, double Pcharge,  bool K, bool Pm, bool Pp){
    //Function to perform particle identification
    double Kmin = 0.8;
    double Pimin= 0.8;
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
    }else if(probPi>Pimin){
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

    if(true){
      PX += H1_PX;
      E1= pow(pow(masses[0],2) + pow(H1_PX,2) + pow(H1_PY,2) + pow(H1_PZ,2) ,0.5);
      PY += H1_PY;
      PZ += H1_PZ;
      totalCharge+=H1_Charge;
      KHit=masses[0]==Kpm;
    }
    if(totalCharge+H2_Charge==0){
      PX += H2_PX;
      E2= pow(pow(masses[1],2) + pow(H2_PX,2) + pow(H2_PY,2) + pow(H2_PZ,2) ,0.5);
      PY += H2_PY;
      PZ += H2_PZ;
      totalCharge+=H2_Charge;
      KHit=KHit||masses[1]==Kpm;
    }
    if(totalCharge+H3_Charge==0){
      PX += H3_PX;
      E3= pow(pow(masses[2],2) + pow(H3_PX,2) + pow(H3_PY,2) + pow(H3_PZ,2) ,0.5);
      PY += H3_PY;
      PZ += H3_PZ;
      totalCharge+=H3_Charge;
      KHit=KHit||masses[2]==Kpm;
    } 
    double resonantMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
    //Limits to remove D
    double maxDMass=1895;
    double minDMass=1835;
    if(resonantMass<maxDMass && resonantMass>minDMass && KHit){
      return false;
    }
    //Add more limits to remove other mesons
    double f0Min = 920;
    double f0Max = 1028;
    
    double rhoMin= 693;
    double rhoMax= 812;
    
    double PK1Min= 1465;
    double PK1Max= 1738;
    
    double PK2Min= 818;
    double PK2Max= 977;
    if(KHit){
      if(resonantMass<PK1Max && resonantMass>PK1Min){
	return false;
      }
      if(resonantMass<PK2Max && resonantMass>PK2Min){
	return false;
      }
    }else{
      if(resonantMass<f0Max && resonantMass>f0Min){
	return false;
      }
      if(resonantMass<rhoMax && resonantMass>rhoMin){
	return false;
      }
    }
    PX = H1_PX + H2_PX + H3_PX;
    PY = H1_PY + H2_PY + H3_PY;
    PZ = H1_PZ + H2_PZ + H3_PZ;
    
    E1 = pow(pow(masses[0],2) + pow(H1_PX,2) + pow(H1_PY,2) + pow(H1_PZ,2) ,0.5);
    E2 = pow(pow(masses[1],2) + pow(H2_PX,2) + pow(H2_PY,2) + pow(H2_PZ,2) ,0.5);
    E3 = pow(pow(masses[2],2) + pow(H3_PX,2) + pow(H3_PY,2) + pow(H3_PZ,2) ,0.5);
    
    double BMass=pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5);
    if(BMass<5200){
      //return false;
    }
    
    return true;
}

void MyAnalysis::Execute()
{
    // This method gets called on every event.
    // In this example the momentum components are filled into histograms.

    // Call the Cut function to decide whether to plot this event or not
    // it returns if the cut function returns false
    if ( !Cut() )
	return;

    // Fill your histograms below.
    // fill the momentum of all three particles 
    h_PX->Fill( H1_PX );
    h_PX->Fill( H2_PX );
    h_PX->Fill( H3_PX );
    // the PY of all three particles
    h_PY->Fill( H1_PY );
    h_PY->Fill( H2_PY );
    h_PY->Fill( H3_PY );
    // the PZ of all three particles
    h_PZ->Fill( H1_PZ );
    h_PZ->Fill( H2_PZ );
    h_PZ->Fill( H3_PZ );
    
    h_PK->Fill( H1_ProbK );
    h_PK->Fill( H2_ProbK );
    h_PK->Fill( H3_ProbK );

    h_PP->Fill( H1_ProbPi );
    //h_PP->Fill( H2_ProbPi );
    //h_PP->Fill( H3_ProbPi );
    // 2D histogram of ProbP vs ProbK
    h_PKP->Fill( H1_ProbK , H1_ProbPi);
    //h_PKP->Fill( H2_ProbK , H2_ProbPi);
    //h_PKP->Fill( H3_ProbK , H3_ProbPi);
    PX = H1_PX + H2_PX + H3_PX;
    PY = H1_PY + H2_PY + H3_PY;
    PZ = H1_PZ + H2_PZ + H3_PZ;
    
    E1 = pow(pow(masses[0],2) + pow(H1_PX,2) + pow(H1_PY,2) + pow(H1_PZ,2) ,0.5);
    E2 = pow(pow(masses[1],2) + pow(H2_PX,2) + pow(H2_PY,2) + pow(H2_PZ,2) ,0.5);
    E3 = pow(pow(masses[2],2) + pow(H3_PX,2) + pow(H3_PY,2) + pow(H3_PZ,2) ,0.5);
    //std::cout << pow(pow(E1+E2+E3,2)+pow(P1+P2+P3,2),0.5) << std::endl;
    h_M->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));

    
    bool onePi    = false;
    bool twoPi    = false;
    bool thrPi    = false;
    double piMass = 139.571;
    double kMass = 493.677;
    PX=0;
    PY=0;
    PZ=0;
    E1=0;
    E2=0;
    E3=0;

    //Charges must sum to 0!!!
    double totalCharge=0;
    bool KHit =false;
    if(true){
      PX += H1_PX;
      E1= pow(pow(masses[0],2) + pow(H1_PX,2) + pow(H1_PY,2) + pow(H1_PZ,2) ,0.5);
      PY += H1_PY;
      PZ += H1_PZ;
      totalCharge+=H1_Charge;
      KHit=masses[0]==kMass;
    }
    if(totalCharge+H2_Charge==0){
      PX += H2_PX;
      E2= pow(pow(masses[1],2) + pow(H2_PX,2) + pow(H2_PY,2) + pow(H2_PZ,2) ,0.5);
      PY += H2_PY;
      PZ += H2_PZ;
      totalCharge+=H2_Charge;
      KHit=KHit || masses[1]==kMass;
    }
    if(totalCharge+H3_Charge==0){
      PX += H3_PX;
      E3= pow(pow(masses[2],2) + pow(H3_PX,2) + pow(H3_PY,2) + pow(H3_PZ,2) ,0.5);
      PY += H3_PY;
      PZ += H3_PZ;
      totalCharge+=H3_Charge;
      KHit=KHit||masses[2]==kMass;
    }
    if(KHit){
      h_RM->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
    }else{
      h_RMP->Fill(pow(pow(E1+E2+E3,2)-pow(PX,2)-pow(PY,2)-pow(PZ,2),0.5));
    }
    
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
