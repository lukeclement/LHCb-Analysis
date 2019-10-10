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
    TH2F           *h_TXTY;
    TH1F           *h_pK;
    TH1F           *h_pP;

    void     BookHistos();

    Bool_t   Cut();
    void     Execute();
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
    v_Histos.push_back( h_PX   = new TH1F("h_PX",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PY   = new TH1F("h_PY",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PZ   = new TH1F("h_PZ",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_TXTY = new TH2F("h_TXTY","", 100, -1,1, 100,-1, 1) );
    v_Histos.push_back( h_pK   = new TH1F("h_pK",  "", 100, -0.1,1.1));
    v_Histos.push_back( h_pP   = new TH1F("h_pP",  "", 100, -0.1,1.1));
}

Bool_t MyAnalysis::Cut()
{
    // This function is called for every event from the Execute
    // function to define whether or not to accept this event.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    // This example checks if the PZ component of particle 3 is greater than 0. 
    if ( false )
	return false;
    else
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
    // 2D histogram of PX/PZ vs PY/PZ
    //h_TXTY->Fill( H1_PX / H1_PZ, H1_PY / H1_PZ );
    
    // Probability histogram
    h_pK->Fill(H1_ProbK);
    h_pP->Fill(H1_ProbPi);
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
