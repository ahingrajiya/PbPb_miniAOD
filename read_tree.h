#include "call_libraries.h" // call libraries from ROOT and C++

// declare variables

// event quantities
float vertexz; // event z vertex
int hiBin;     // event centrality (used if use_centrality = true in input_variables.h)

// trigger quantities
int jet_trigger_bit; // jet HLT path trigger used for analysis (jet_trigger variable in input_variables.h)

// reco jets
int nref;           // number of jets
float jtpt[99999];  // jet pT after JEC (not always true, sometimes we need to correct rawpt a posterior, as done in pPb)
float jteta[99999]; // jet eta
float jtphi[99999]; // jet phi
float rawpt[99999]; // jet pT without JEC
float trackMax[9999];

// reco tracks
int ntrk;                            // number of track
std::vector<float> *trkpt = 0;       // track pT
std::vector<float> *trketa = 0;      // track eta
std::vector<float> *trkphi = 0;      // track phi
std::vector<float> *trkpterr = 0;    // track pT error (uncertainty)
std::vector<float> *trkdcaxy = 0;    // track dxy impact parameter (transverse distance between primary vertex and collision - distance of closest approuch - DCA)
std::vector<float> *trkdcaz = 0;     // track dz impact parameter (longitudinal distance between primary vertex and collision - distance of closest approuch - DCA)
std::vector<float> *trkdcaxyerr = 0; // track dxy error (uncertainty)
std::vector<float> *trkdcazerr = 0;  // track dxy error (uncertainty)
std::vector<float> *trkchi2 = 0;     // track reconstruction chi2 of the fitting
std::vector<float> *pfEcal = 0;      // particle flow energy deposit in ECAL
std::vector<float> *pfHcal = 0;      // particle flow energy deposit in HCAL
std::vector<float> *trkmva = 0;      // track mva for each step
std::vector<char> *trkcharge = 0;    // track charge
std::vector<char> *trknhits = 0;     // number of hits in the tracker
std::vector<char> *trknlayer = 0;    // number of layers with measurement in the tracker
std::vector<bool> *highpur;          // tracker steps MVA selection

// events quantities from gen
float weight; // event weight --> pthat weight
float pthat;  // pthat (initial parton pT)

// gen jets
int ngen;               // number of gen jets
float gen_jtpt[99999];  // gen jet pT
float gen_jteta[99999]; // gen jet eta
float gen_jtphi[99999]; // gen jet phi

// matched jets
float refpt[99999];              // jet pT matched with Gen pT
float refeta[99999];             // jet eta matched with Gen eta
float refphi[99999];             // jet phi matched with Gen phi
int refparton_flavor[99999];     // jet phi matched with Gen phi
int refparton_flavorForB[99999]; // jet phi matched with Gen phisour

// gen tracks
std::vector<float> *gen_trkpt = 0;  // gen particle pT
std::vector<float> *gen_trketa = 0; // gen particle eta
std::vector<float> *gen_trkphi = 0; // gen particle phi
std::vector<int> *gen_trkchg = 0;   // gen particle charge
std::vector<int> *gen_trkpid = 0;   // gen particle pid
std::vector<int> *gen_trksube = 0;  // gen particle pid

// All variables listed above are readed in the function bellow
/*
Function to read the Forest/Skim tree
Arguments ->  transfer quantities from trees to our variables
tree: input TChain from jet_analyzer.C file
all the arguments bellow are defined in input_variables.h
is_MC: true -> MC; false -> Data
use_WTA: true -> use WTA (winner-takes-all); false -> use E-Scheme
jet_trigger: string with trigger name
colliding_system: pp, pPb, PbPb, XeXe, ... (future)
colliding_energy: colliding energy in GeV -> 5020, 5440, 8160, 13000, ...
year_of_datataking: year of data taking
event_filterstr: string of event filters
event_filters: integer (0 or 1) from event filters
*/
void read_tree(TChain *tree, bool is_MC, bool use_WTA, TString jet_trigger, TString colliding_system, int colliding_energy, int year_of_datataking, std::vector<TString> event_filterstr, std::vector<int> event_filters)
{

    tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files

    // enable branches of interest -> see definition of each variables above

    // event quantities
    tree->SetBranchStatus(Form("%s", jet_trigger.Data()), 1);
    tree->SetBranchStatus("vz", 1);
    if (colliding_system == "PbPb" || colliding_system == "XeXe")
        tree->SetBranchStatus("hiBin", 1); // centrality only for PbPb and XeXe
    for (int i = 0; i < event_filterstr.size(); i++)
        tree->SetBranchStatus(Form("%s", event_filterstr[i].Data()), 1); // event filters

    tree->SetBranchAddress(Form("%s", jet_trigger.Data()), &jet_trigger_bit);
    tree->SetBranchAddress("vz", &vertexz);
    if (colliding_system == "PbPb" || colliding_system == "XeXe")
        tree->SetBranchAddress("hiBin", &hiBin); // centrality only for PbPb and XeXe
    for (int i = 0; i < event_filters.size(); i++)
        tree->SetBranchAddress(Form("%s", event_filterstr[i].Data()), &event_filters[i]);

    if (is_MC)
    {
        tree->SetBranchStatus("weight", 1);
        tree->SetBranchStatus("pthat", 1);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("pthat", &pthat);
    }

    // jet quantities
    tree->SetBranchStatus("nref", 1);
    //    tree->SetBranchStatus("jtpt", 1);
    tree->SetBranchStatus("rawpt", 1);
    tree->SetBranchStatus("trackMax", 1);
    if (use_WTA)
    {
        tree->SetBranchStatus("WTAeta", 1);
        tree->SetBranchStatus("WTAphi", 1);
    }
    else
    {
        tree->SetBranchStatus("jteta", 1);
        tree->SetBranchStatus("jtphi", 1);
    }

    tree->SetBranchAddress("nref", &nref);
    //    tree->SetBranchAddress("jtpt", &jtpt);
    tree->SetBranchAddress("rawpt", &rawpt);
    tree->SetBranchAddress("trackMax", &trackMax);
    if (use_WTA)
    {
        tree->SetBranchAddress("WTAeta", &jteta);
        tree->SetBranchAddress("WTAphi", &jtphi);
    }
    else
    {
        tree->SetBranchAddress("jteta", &jteta);
        tree->SetBranchAddress("jtphi", &jtphi);
    }

    // gen jet quantities
    if (is_MC)
    {
        tree->SetBranchStatus("ngen", 1);
        tree->SetBranchStatus("genpt", 1);
        if (use_WTA)
        {
            tree->SetBranchStatus("WTAgeneta", 1);
            tree->SetBranchStatus("WTAgenphi", 1);
        }
        else
        {
            tree->SetBranchStatus("geneta", 1);
            tree->SetBranchStatus("genphi", 1);
        }

        tree->SetBranchAddress("ngen", &ngen);
        tree->SetBranchAddress("genpt", &gen_jtpt);
        if (use_WTA)
        {
            tree->SetBranchAddress("WTAgeneta", &gen_jteta);
            tree->SetBranchAddress("WTAgenphi", &gen_jtphi);
        }
        else
        {
            tree->SetBranchAddress("geneta", &gen_jteta);
            tree->SetBranchAddress("genphi", &gen_jtphi);
        }
    }
    // matching quantities
    if (is_MC)
    {
        tree->SetBranchStatus("refpt", 1);
        tree->SetBranchAddress("refpt", &refpt);
        // tree->SetBranchStatus("refeta", 1);
        // tree->SetBranchAddress("refeta", &refeta);
        // tree->SetBranchStatus("refphi", 1);
        // tree->SetBranchAddress("refphi", &refphi);
        tree->SetBranchStatus("refparton_flavor", 1);
        tree->SetBranchAddress("refparton_flavor", &refparton_flavor);
        tree->SetBranchStatus("refparton_flavorForB", 1);
        tree->SetBranchAddress("refparton_flavorForB", &refparton_flavorForB);
    }

    // track quantities
    tree->SetBranchStatus("nTrk", 1);
    tree->SetBranchStatus("trkPt", 1);
    tree->SetBranchStatus("trkEta", 1);
    tree->SetBranchStatus("trkPhi", 1);
    tree->SetBranchStatus("trkPtError", 1);
    tree->SetBranchStatus("trkDxyFirstVtx", 1);
    tree->SetBranchStatus("trkDxyErrFirstVtx", 1);
    tree->SetBranchStatus("trkDzFirstVtx", 1);
    tree->SetBranchStatus("trkDzErrFirstVtx", 1);
    tree->SetBranchStatus("trkCharge", 1);
    tree->SetBranchStatus("highPurity", 1);
    tree->SetBranchStatus("pfEcal", 1);
    tree->SetBranchStatus("pfHcal", 1);

    tree->SetBranchAddress("nTrk", &ntrk);
    tree->SetBranchAddress("trkPt", &trkpt);
    tree->SetBranchAddress("trkEta", &trketa);
    tree->SetBranchAddress("trkPhi", &trkphi);
    tree->SetBranchAddress("trkPtError", &trkpterr);
    tree->SetBranchAddress("trkDxyFirstVtx", &trkdcaxy);
    tree->SetBranchAddress("trkDxyErrFirstVtx", &trkdcaxyerr);
    tree->SetBranchAddress("trkDzFirstVtx", &trkdcaz);
    tree->SetBranchAddress("trkDzErrFirstVtx", &trkdcazerr);
    tree->SetBranchAddress("trkCharge", &trkcharge);
    tree->SetBranchAddress("pfEcal", &pfEcal);
    tree->SetBranchAddress("pfHcal", &pfHcal);
    tree->SetBranchAddress("highPurity", &highpur);
    if (colliding_system == "PbPb" || colliding_system == "XeXe")
    {
        tree->SetBranchStatus("trkNormChi2", 1);
        tree->SetBranchStatus("trkNHits", 1);
        tree->SetBranchStatus("trkNLayers", 1);
        tree->SetBranchAddress("trkNormChi2", &trkchi2);
        tree->SetBranchAddress("trkNHits", &trknhits);
        tree->SetBranchAddress("trkNLayers", &trknlayer);
    }

    if (colliding_system == "PbPb" && colliding_energy == 5020 && year_of_datataking == 2018)
    { // special for 2018 PbPb MC
      // tree->SetBranchStatus("trkMVA", 1);
      // tree->SetBranchStatus("trkAlgo", 1);
      // tree->SetBranchAddress("trkMVA", &trkmva);
      // tree->SetBranchAddress("trkAlgo", &trkalgo);
    }

    // gen particle quantities
    if (is_MC)
    {
        tree->SetBranchStatus("pt", 1);
        tree->SetBranchStatus("eta", 1);
        tree->SetBranchStatus("phi", 1);
        tree->SetBranchStatus("chg", 1);
        tree->SetBranchStatus("pdg", 1);
        tree->SetBranchStatus("sube", 1);

        tree->SetBranchAddress("pt", &gen_trkpt);
        tree->SetBranchAddress("eta", &gen_trketa);
        tree->SetBranchAddress("phi", &gen_trkphi);
        tree->SetBranchAddress("chg", &gen_trkchg);
        tree->SetBranchAddress("pdg", &gen_trkpid);
        tree->SetBranchAddress("sube", &gen_trksube);
    }
}
