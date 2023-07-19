#include "call_libraries.h"       // call libraries from ROOT and C++
#include "input_variables_PbPb.h" // call inputs

int trkbinsize = (int)trk_pt_bins.size(); // track bins for jet-track correlation
int multbinsize =
    (int)multiplicity_centrality_bins.size(); // multiplicity or centrality bins for jet-track correlation

// -------------------------------- QA plots --------------------------------
// Event quantities
TH1I *Nevents = new TH1I("Nevents", "Nevents", 10, 0, 10);
TH1I *Nev_recoreco = new TH1I("Nev_recoreco", "Nev_recoreco", 1, 0, 1);
TH1I *Nev_recoreco_lead = new TH1I("Nev_recoreco_lead", "Nev_recoreco_lead", 1, 0, 1);
TH1I *Nev_recoreco_subl = new TH1I("Nev_recoreco_subl", "Nev_recoreco_subl", 1, 0, 1);
TH1I *Nev_recogen = new TH1I("Nev_recogen", "Nev_recogen", 1, 0, 1);
TH1I *Nev_recogen_lead = new TH1I("Nev_recogen_lead", "Nev_recogen_lead", 1, 0, 1);
TH1I *Nev_recogen_subl = new TH1I("Nev_recogen_subl", "Nev_recogen_subl", 1, 0, 1);
TH1I *Nev_genreco = new TH1I("Nev_genreco", "Nev_genreco", 1, 0, 1);
TH1I *Nev_genreco_lead = new TH1I("Nev_genreco_lead", "Nev_genreco_lead", 1, 0, 1);
TH1I *Nev_genreco_subl = new TH1I("Nev_genreco_subl", "Nev_genreco_subl", 1, 0, 1);
TH1I *Nev_gengen = new TH1I("Nev_gengen", "Nev_gengen", 1, 0, 1);
TH1I *Nev_gengen_lead = new TH1I("Nev_gengen_lead", "Nev_gengen_lead", 1, 0, 1);
TH1I *Nev_gengen_subl = new TH1I("Nev_gengen_subl", "Nev_gengen_subl", 1, 0, 1);
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 500, 0.0, 500.0);
TH1D *multiplicity_weighted = new TH1D("multiplicity_weighted", "multiplicity_weighted", 500, 0.0, 500.0);
TH1D *vzhist = new TH1D("vzhist", "vzhist", 80, -20, 20);
TH1D *vzhist_weighted = new TH1D("vzhist_weighted", "vzhist_weighted", 80, -20, 20);
TH1D *vzhist185_weighted = new TH1D("vzhist185_weighted", "vzhist185_weighted", 80, -20, 20);
TH1D *vzhist250_weighted = new TH1D("vzhist250_weighted", "vzhist250_weighted", 80, -20, 20);

TH1D *pthathist = new TH1D("pthathist", "pthathist", 230, 0, 460);
TH1D *pthathist_weighted = new TH1D("pthathist_weighted", "pthathist_weighted", 230, 0, 460);
TH1D *centrality = new TH1D("centrality", "centrality", 200, 0.0, 200.0);
TH1D *centrality_weighted = new TH1D("centrality_weight", "centrality_weight", 200, 0.0, 200.0);

// Track/Particle histograms
int bins4D_trk[4] = {500, 50, 64, multbinsize - 1};
double xmin4D_trk[4] = {0.0, -2.5, -TMath::Pi(), 0};
double xmax4D_trk[4] = {20.0, 2.5, TMath::Pi(), (double)multbinsize - 1};

// --> Reco
THnSparseD *hist_reco_trk = new THnSparseD("hist_reco_trk", "hist_reco_trk", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_reco_trk_corr =
    new THnSparseD("hist_reco_trk_corr", "hist_reco_trk_corr", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_reco_trk_weighted =
    new THnSparseD("hist_reco_trk_weighted", "hist_reco_trk_weighted", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// --> Gen
THnSparseD *hist_gen_trk = new THnSparseD("hist_gen_trk", "hist_gen_trk", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_gen_trk_weighted =
    new THnSparseD("hist_gen_trk_weighted", "hist_gen_trk_weighted", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// Tracks from Correlations
// Inclusive
THnSparseD *hist_trk_from_reco_reco_sig =
    new THnSparseD("hist_trk_from_reco_reco_sig", "hist_trk_from_reco_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_reco_gen_sig =
    new THnSparseD("hist_trk_from_reco_gen_sig", "hist_trk_from_reco_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_reco_sig =
    new THnSparseD("hist_trk_from_gen_reco_sig", "hist_trk_from_gen_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_gen_sig =
    new THnSparseD("hist_trk_from_gen_gen_sig", "hist_trk_from_gen_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_reco_reco_mix =
    new THnSparseD("hist_trk_from_reco_reco_mix", "hist_trk_from_reco_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_reco_gen_mix =
    new THnSparseD("hist_trk_from_reco_gen_mix", "hist_trk_from_reco_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_reco_mix =
    new THnSparseD("hist_trk_from_gen_reco_mix", "hist_trk_from_gen_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_gen_mix =
    new THnSparseD("hist_trk_from_gen_gen_mix", "hist_trk_from_gen_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
// Leading
THnSparseD *hist_LJ_trk_from_reco_reco_sig = new THnSparseD(
    "hist_LJ_trk_from_reco_reco_sig", "hist_LJ_trk_from_reco_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_reco_gen_sig = new THnSparseD(
    "hist_LJ_trk_from_reco_gen_sig", "hist_LJ_trk_from_reco_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_reco_sig = new THnSparseD(
    "hist_LJ_trk_from_gen_reco_sig", "hist_LJ_trk_from_gen_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_gen_sig = new THnSparseD(
    "hist_LJ_trk_from_gen_gen_sig", "hist_LJ_trk_from_gen_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_reco_reco_mix = new THnSparseD(
    "hist_LJ_trk_from_reco_reco_mix", "hist_LJ_trk_from_reco_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_reco_gen_mix = new THnSparseD(
    "hist_LJ_trk_from_reco_gen_mix", "hist_LJ_trk_from_reco_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_reco_mix = new THnSparseD(
    "hist_LJ_trk_from_gen_reco_mix", "hist_LJ_trk_from_gen_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_gen_mix = new THnSparseD(
    "hist_LJ_trk_from_gen_gen_mix", "hist_LJ_trk_from_gen_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
// Subleading
THnSparseD *hist_SLJ_trk_from_reco_reco_sig = new THnSparseD(
    "hist_SLJ_trk_from_reco_reco_sig", "hist_SLJ_trk_from_reco_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_reco_gen_sig = new THnSparseD(
    "hist_SLJ_trk_from_reco_gen_sig", "hist_SLJ_trk_from_reco_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_reco_sig = new THnSparseD(
    "hist_SLJ_trk_from_gen_reco_sig", "hist_SLJ_trk_from_gen_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_gen_sig = new THnSparseD(
    "hist_SLJ_trk_from_gen_gen_sig", "hist_SLJ_trk_from_gen_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_reco_reco_mix = new THnSparseD(
    "hist_SLJ_trk_from_reco_reco_mix", "hist_SLJ_trk_from_reco_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_reco_gen_mix = new THnSparseD(
    "hist_SLJ_trk_from_reco_gen_mix", "hist_SLJ_trk_from_reco_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_reco_mix = new THnSparseD(
    "hist_SLJ_trk_from_gen_reco_mix", "hist_SLJ_trk_from_gen_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_gen_mix = new THnSparseD(
    "hist_SLJ_trk_from_gen_gen_mix", "hist_SLJ_trk_from_gen_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// Jet histograms
int bins4D_jet[4] = {100, 80, 64, multbinsize - 1};
double xmin4D_jet[4] = {0.0, -4.0, -TMath::Pi(), 0};
double xmax4D_jet[4] = {500.0, 4.0, TMath::Pi(), (double)multbinsize - 1};

// --> Reco
TH1D *hist_reco_jet_weighted_nocut =
    new TH1D("hist_reco_jet_weighted_nocut", "hist_reco_jet_weighted_nocut", 100, 0.0, 500.0);
THnSparseD *hist_reco_jet = new THnSparseD("hist_reco_jet", "hist_reco_jet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_jet_corr =
    new THnSparseD("hist_reco_jet_corr", "hist_reco_jet_corr", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_jet_weighted =
    new THnSparseD("hist_reco_jet_weighted", "hist_reco_jet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_jet_corr_weighted =
    new THnSparseD("hist_reco_jet_corr_weighted", "hist_reco_jet_corr_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
// --> Gen
TH1D *hist_gen_jet_weighted_nocut =
    new TH1D("hist_gen_jet_weighted_nocut", "hist_gen_jet_weighted_nocut", 100, 0.0, 500.0);
THnSparseD *hist_gen_jet = new THnSparseD("hist_gen_jet", "hist_gen_jet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_jet_weighted =
    new THnSparseD("hist_gen_jet_weighted", "hist_gen_jet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// Leading and Subleading Jets
// --> Reco
TH1D *hist_reco_leadjet_pt_nocut =
    new TH1D("hist_reco_leadjet_pt_nocut", "hist_reco_leadjet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_reco_leadjet_pt_nocut_weighted =
    new TH1D("hist_reco_leadjet_pt_nocut_weighted", "hist_reco_leadjet_pt_nocut_weighted", 100, 0.0, 500.0);
TH1D *hist_reco_subljet_pt_nocut =
    new TH1D("hist_reco_subljet_pt_nocut", "hist_reco_subljet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_reco_subljet_pt_nocut_weighted =
    new TH1D("hist_reco_subljet_pt_nocut_weighted", "hist_reco_subljet_pt_nocut_weighted", 100, 0.0, 500.0);

THnSparseD *hist_reco_leadjet =
    new THnSparseD("hist_reco_leadjet", "hist_reco_leadjet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_leadjet_weighted =
    new THnSparseD("hist_reco_leadjet_weighted", "hist_reco_leadjet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_subljet =
    new THnSparseD("hist_reco_subljet", "hist_reco_subljet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_subljet_weighted =
    new THnSparseD("hist_reco_subljet_weighted", "hist_reco_subljet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// --> Gen
TH1D *hist_gen_leadjet_pt_nocut = new TH1D("hist_gen_leadjet_pt_nocut", "hist_gen_leadjet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_gen_leadjet_pt_nocut_weighted =
    new TH1D("hist_gen_leadjet_pt_nocut_weighted", "hist_gen_leadjet_pt_nocut_weighted", 100, 0.0, 500.0);
TH1D *hist_gen_subljet_pt_nocut = new TH1D("hist_gen_subljet_pt_nocut", "hist_gen_subljet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_gen_subljet_pt_nocut_weighted =
    new TH1D("hist_gen_subljet_pt_nocut_weighted", "hist_gen_subljet_pt_nocut_weighted", 100, 0.0, 500.0);

THnSparseD *hist_gen_leadjet =
    new THnSparseD("hist_gen_leadjet", "hist_gen_leadjet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_leadjet_weighted =
    new THnSparseD("hist_gen_leadjet_weighted", "hist_gen_leadjet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_subljet =
    new THnSparseD("hist_gen_subljet", "hist_gen_subljet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_subljet_weighted =
    new THnSparseD("hist_gen_subljet_weighted", "hist_gen_subljet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// Jets from Correlations
// Inclusive
THnSparseD *hist_jet_from_reco_reco_sig =
    new THnSparseD("hist_jet_from_reco_reco_sig", "hist_jet_from_reco_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_reco_gen_sig =
    new THnSparseD("hist_jet_from_reco_gen_sig", "hist_jet_from_reco_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_reco_sig =
    new THnSparseD("hist_jet_from_gen_reco_sig", "hist_jet_from_gen_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_gen_sig =
    new THnSparseD("hist_jet_from_gen_gen_sig", "hist_jet_from_gen_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_reco_reco_mix =
    new THnSparseD("hist_jet_from_reco_reco_mix", "hist_jet_from_reco_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_reco_gen_mix =
    new THnSparseD("hist_jet_from_reco_gen_mix", "hist_jet_from_reco_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_reco_mix =
    new THnSparseD("hist_jet_from_gen_reco_mix", "hist_jet_from_gen_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_gen_mix =
    new THnSparseD("hist_jet_from_gen_gen_mix", "hist_jet_from_gen_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
// Leading
THnSparseD *hist_lead_jet_from_reco_reco_sig = new THnSparseD(
    "hist_lead_jet_from_reco_reco_sig", "hist_lead_jet_from_reco_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_reco_gen_sig = new THnSparseD(
    "hist_lead_jet_from_reco_gen_sig", "hist_lead_jet_from_reco_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_reco_sig = new THnSparseD(
    "hist_lead_jet_from_gen_reco_sig", "hist_lead_jet_from_gen_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_gen_sig = new THnSparseD(
    "hist_lead_jet_from_gen_gen_sig", "hist_lead_jet_from_gen_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_reco_reco_mix = new THnSparseD(
    "hist_lead_jet_from_reco_reco_mix", "hist_lead_jet_from_reco_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_reco_gen_mix = new THnSparseD(
    "hist_lead_jet_from_reco_gen_mix", "hist_lead_jet_from_reco_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_reco_mix = new THnSparseD(
    "hist_lead_jet_from_gen_reco_mix", "hist_lead_jet_from_gen_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_gen_mix = new THnSparseD(
    "hist_lead_jet_from_gen_gen_mix", "hist_lead_jet_from_gen_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
// Subleading
THnSparseD *hist_subl_jet_from_reco_reco_sig = new THnSparseD(
    "hist_subl_jet_from_reco_reco_sig", "hist_subl_jet_from_reco_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_reco_gen_sig = new THnSparseD(
    "hist_subl_jet_from_reco_gen_sig", "hist_subl_jet_from_reco_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_reco_sig = new THnSparseD(
    "hist_subl_jet_from_gen_reco_sig", "hist_subl_jet_from_gen_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_gen_sig = new THnSparseD(
    "hist_subl_jet_from_gen_gen_sig", "hist_subl_jet_from_gen_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_reco_reco_mix = new THnSparseD(
    "hist_subl_jet_from_reco_reco_mix", "hist_subl_jet_from_reco_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_reco_gen_mix = new THnSparseD(
    "hist_subl_jet_from_reco_gen_mix", "hist_subl_jet_from_reco_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_reco_mix = new THnSparseD(
    "hist_subl_jet_from_gen_reco_mix", "hist_subl_jet_from_gen_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_gen_mix = new THnSparseD(
    "hist_subl_jet_from_gen_gen_mix", "hist_subl_jet_from_gen_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// --------------------------------------------------------------------------------------------------------
// Quenching studies
// Axis : 0 -> Aj, 1 -> Xj, 2 -> delta phi, 3 -> multiplicity
int bins4D_quenc[4] = {20, 20, 40, multbinsize - 1};
double xmin4D_quenc[4] = {0.0, 0.0, 0.0, 0};
double xmax4D_quenc[4] = {1.0, 1.0, TMath::Pi(), (double)multbinsize - 1};
THnSparseD *hist_reco_lead_reco_subl_quench = new THnSparseD(
    "hist_reco_lead_reco_subl_quench", "hist_reco_lead_reco_subl_quench", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench = new THnSparseD(
    "hist_gen_lead_gen_subl_quench", "hist_gen_lead_gen_subl_quench", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);

// Axis : 0 -> Aj, 1 -> Xj, 2 -> delta phi 2PC, 3 -> multiplicity
int bins4D_quenc_2pc[4] = {20, 20, 40, multbinsize - 1};
double xmin4D_quenc_2pc[4] = {0.0, 0.0, -TMath::Pi() / 2.0, 0};
double xmax4D_quenc_2pc[4] = {1.0, 1.0, 3.0 * TMath::Pi() / 2.0, (double)multbinsize - 1};
THnSparseD *hist_reco_lead_reco_subl_quench2pc = new THnSparseD("hist_reco_lead_reco_subl_quench2pc",
                                                                "hist_reco_lead_reco_subl_quench2pc",
                                                                4,
                                                                bins4D_quenc_2pc,
                                                                xmin4D_quenc_2pc,
                                                                xmax4D_quenc_2pc);
THnSparseD *hist_gen_lead_gen_subl_quench2pc = new THnSparseD("hist_gen_lead_gen_subl_quench2pc",
                                                              "hist_gen_lead_gen_subl_quench2pc",
                                                              4,
                                                              bins4D_quenc_2pc,
                                                              xmin4D_quenc_2pc,
                                                              xmax4D_quenc_2pc);

// Correlation studies
// Axis : 0 -> delta phi, 1 -> delta eta, 2 -> track pT, 3 -> multiplicity
int bins4D_jettrk[4] = {200, 500, trkbinsize - 1, multbinsize - 1};
double xmin4D_jettrk[4] = {-TMath::Pi() / 2.0, -5.0, 0, 0};
double xmax4D_jettrk[4] = {3.0 * TMath::Pi() / 2.0, 5.0, (double)trkbinsize - 1, (double)multbinsize - 1};

// Correlation: Reco Jet + Reco Track
THnSparseD *hist_correlation_signal_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_jet_reco_track_reco",
                                                                         "hist_correlation_signal_jet_reco_track_reco",
                                                                         4,
                                                                         bins4D_jettrk,
                                                                         xmin4D_jettrk,
                                                                         xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_reco_track_reco =
    new THnSparseD("hist_correlation_rotation_jet_reco_track_reco",
                   "hist_correlation_rotation_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_reco_track_reco = new THnSparseD("hist_correlation_mixing_jet_reco_track_reco",
                                                                         "hist_correlation_mixing_jet_reco_track_reco",
                                                                         4,
                                                                         bins4D_jettrk,
                                                                         xmin4D_jettrk,
                                                                         xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_reco_track_reco =
    new THnSparseD("hist_correlation_signal_lead_jet_reco_track_reco",
                   "hist_correlation_signal_lead_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_reco_track_reco =
    new THnSparseD("hist_correlation_rotation_lead_jet_reco_track_reco",
                   "hist_correlation_rotation_lead_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_reco_track_reco =
    new THnSparseD("hist_correlation_mixing_lead_jet_reco_track_reco",
                   "hist_correlation_mixing_lead_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_reco_track_reco =
    new THnSparseD("hist_correlation_signal_subl_jet_reco_track_reco",
                   "hist_correlation_signal_subl_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_reco_track_reco =
    new THnSparseD("hist_correlation_rotation_subl_jet_reco_track_reco",
                   "hist_correlation_rotation_subl_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_reco_track_reco =
    new THnSparseD("hist_correlation_mixing_subl_jet_reco_track_reco",
                   "hist_correlation_mixing_subl_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);

// Correlation: Reco Jet + Gen Track
THnSparseD *hist_correlation_signal_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_jet_reco_track_gen",
                                                                        "hist_correlation_signal_jet_reco_track_gen",
                                                                        4,
                                                                        bins4D_jettrk,
                                                                        xmin4D_jettrk,
                                                                        xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_reco_track_gen =
    new THnSparseD("hist_correlation_rotation_jet_reco_track_gen",
                   "hist_correlation_rotation_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_reco_track_gen = new THnSparseD("hist_correlation_mixing_jet_reco_track_gen",
                                                                        "hist_correlation_mixing_jet_reco_track_gen",
                                                                        4,
                                                                        bins4D_jettrk,
                                                                        xmin4D_jettrk,
                                                                        xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_reco_track_gen =
    new THnSparseD("hist_correlation_signal_lead_jet_reco_track_gen",
                   "hist_correlation_signal_lead_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_reco_track_gen =
    new THnSparseD("hist_correlation_rotation_lead_jet_reco_track_gen",
                   "hist_correlation_rotation_lead_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_reco_track_gen =
    new THnSparseD("hist_correlation_mixing_lead_jet_reco_track_gen",
                   "hist_correlation_mixing_lead_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_reco_track_gen =
    new THnSparseD("hist_correlation_signal_subl_jet_reco_track_gen",
                   "hist_correlation_signal_subl_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_reco_track_gen =
    new THnSparseD("hist_correlation_rotation_subl_jet_reco_track_gen",
                   "hist_correlation_rotation_subl_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_reco_track_gen =
    new THnSparseD("hist_correlation_mixing_subl_jet_reco_track_gen",
                   "hist_correlation_mixing_subl_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);

// Correlation: Gen Jet + Reco Track
THnSparseD *hist_correlation_signal_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_jet_gen_track_reco",
                                                                        "hist_correlation_signal_jet_gen_track_reco",
                                                                        4,
                                                                        bins4D_jettrk,
                                                                        xmin4D_jettrk,
                                                                        xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_gen_track_reco =
    new THnSparseD("hist_correlation_rotation_jet_gen_track_reco",
                   "hist_correlation_rotation_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_gen_track_reco = new THnSparseD("hist_correlation_mixing_jet_gen_track_reco",
                                                                        "hist_correlation_mixing_jet_gen_track_reco",
                                                                        4,
                                                                        bins4D_jettrk,
                                                                        xmin4D_jettrk,
                                                                        xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_gen_track_reco =
    new THnSparseD("hist_correlation_signal_lead_jet_gen_track_reco",
                   "hist_correlation_signal_lead_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_gen_track_reco =
    new THnSparseD("hist_correlation_rotation_lead_jet_gen_track_reco",
                   "hist_correlation_rotation_lead_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_gen_track_reco =
    new THnSparseD("hist_correlation_mixing_lead_jet_gen_track_reco",
                   "hist_correlation_mixing_lead_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_gen_track_reco =
    new THnSparseD("hist_correlation_signal_subl_jet_gen_track_reco",
                   "hist_correlation_signal_subl_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_gen_track_reco =
    new THnSparseD("hist_correlation_rotation_subl_jet_gen_track_reco",
                   "hist_correlation_rotation_subl_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_gen_track_reco =
    new THnSparseD("hist_correlation_mixing_subl_jet_gen_track_reco",
                   "hist_correlation_mixing_subl_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);

// Correlation: Gen Jet + Gen Track
THnSparseD *hist_correlation_signal_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_jet_gen_track_gen",
                                                                       "hist_correlation_signal_jet_gen_track_gen",
                                                                       4,
                                                                       bins4D_jettrk,
                                                                       xmin4D_jettrk,
                                                                       xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_gen_track_gen = new THnSparseD("hist_correlation_rotation_jet_gen_track_gen",
                                                                         "hist_correlation_rotation_jet_gen_track_gen",
                                                                         4,
                                                                         bins4D_jettrk,
                                                                         xmin4D_jettrk,
                                                                         xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_gen_track_gen = new THnSparseD("hist_correlation_mixing_jet_gen_track_gen",
                                                                       "hist_correlation_mixing_jet_gen_track_gen",
                                                                       4,
                                                                       bins4D_jettrk,
                                                                       xmin4D_jettrk,
                                                                       xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_gen_track_gen =
    new THnSparseD("hist_correlation_signal_lead_jet_gen_track_gen",
                   "hist_correlation_signal_lead_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_gen_track_gen =
    new THnSparseD("hist_correlation_rotation_lead_jet_gen_track_gen",
                   "hist_correlation_rotation_lead_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_gen_track_gen =
    new THnSparseD("hist_correlation_mixing_lead_jet_gen_track_gen",
                   "hist_correlation_mixing_lead_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_gen_track_gen =
    new THnSparseD("hist_correlation_signal_subl_jet_gen_track_gen",
                   "hist_correlation_signal_subl_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_gen_track_gen =
    new THnSparseD("hist_correlation_rotation_subl_jet_gen_track_gen",
                   "hist_correlation_rotation_subl_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_gen_track_gen =
    new THnSparseD("hist_correlation_mixing_subl_jet_gen_track_gen",
                   "hist_correlation_mixing_subl_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);

// Correlation: Include sube > 0
THnSparseD *hist_correlation_signal_subg0_jet_reco_track_reco =
    new THnSparseD("hist_correlation_signal_subg0_jet_reco_track_reco",
                   "hist_correlation_signal_subg0_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_reco_track_gen =
    new THnSparseD("hist_correlation_signal_subg0_jet_reco_track_gen",
                   "hist_correlation_signal_subg0_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_gen_track_reco =
    new THnSparseD("hist_correlation_signal_subg0_jet_gen_track_reco",
                   "hist_correlation_signal_subg0_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_gen_track_gen =
    new THnSparseD("hist_correlation_signal_subg0_jet_gen_track_gen",
                   "hist_correlation_signal_subg0_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);

THnSparseD *hist_correlation_signal_subg0_lead_jet_reco_track_reco =
    new THnSparseD("hist_correlation_signal_subg0_lead_jet_reco_track_reco",
                   "hist_correlation_signal_subg0_lead_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_reco_track_gen =
    new THnSparseD("hist_correlation_signal_subg0_lead_jet_reco_track_gen",
                   "hist_correlation_signal_subg0_lead_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_gen_track_reco =
    new THnSparseD("hist_correlation_signal_subg0_lead_jet_gen_track_reco",
                   "hist_correlation_signal_subg0_lead_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_gen_track_gen =
    new THnSparseD("hist_correlation_signal_subg0_lead_jet_gen_track_gen",
                   "hist_correlation_signal_subg0_lead_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);

THnSparseD *hist_correlation_signal_subg0_subl_jet_reco_track_reco =
    new THnSparseD("hist_correlation_signal_subg0_subl_jet_reco_track_reco",
                   "hist_correlation_signal_subg0_subl_jet_reco_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_reco_track_gen =
    new THnSparseD("hist_correlation_signal_subg0_subl_jet_reco_track_gen",
                   "hist_correlation_signal_subg0_subl_jet_reco_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_gen_track_reco =
    new THnSparseD("hist_correlation_signal_subg0_subl_jet_gen_track_reco",
                   "hist_correlation_signal_subg0_subl_jet_gen_track_reco",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_gen_track_gen =
    new THnSparseD("hist_correlation_signal_subg0_subl_jet_gen_track_gen",
                   "hist_correlation_signal_subg0_subl_jet_gen_track_gen",
                   4,
                   bins4D_jettrk,
                   xmin4D_jettrk,
                   xmax4D_jettrk);

// histograms for matched jets and parton flavor studies
TH1D *hist_matched_jet_weighted_nocut =
    new TH1D("hist_matched_jet_weighted_nocut", "hist_matched_jet_weighted_nocut", 100, 0.0, 500.0);
THnSparseD *hist_matched_jet =
    new THnSparseD("hist_matched_jet", "hist_matched_jet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_matched_jet_weighted =
    new THnSparseD("hist_matched_jet_weighted", "hist_matched_jet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

TH1D *hist_matched_jet_pt_parton_from_u =
    new TH1D("hist_matched_jet_pt_parton_from_u", "hist_matched_jet_pt_parton_from_u", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_from_d =
    new TH1D("hist_matched_jet_pt_parton_from_d", "hist_matched_jet_pt_parton_from_d", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_from_s =
    new TH1D("hist_matched_jet_pt_parton_from_s", "hist_matched_jet_pt_parton_from_s", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_from_c =
    new TH1D("hist_matched_jet_pt_parton_from_c", "hist_matched_jet_pt_parton_from_c", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_from_b =
    new TH1D("hist_matched_jet_pt_parton_from_b", "hist_matched_jet_pt_parton_from_b", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_from_t =
    new TH1D("hist_matched_jet_pt_parton_from_t", "hist_matched_jet_pt_parton_from_t", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_from_g =
    new TH1D("hist_matched_jet_pt_parton_from_g", "hist_matched_jet_pt_parton_from_g", 100, 0.0, 500.0);

TH1D *hist_matched_jet_pt_parton_B_from_u =
    new TH1D("hist_matched_jet_pt_parton_B_from_u", "hist_matched_jet_pt_parton_B_from_u", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_B_from_d =
    new TH1D("hist_matched_jet_pt_parton_B_from_d", "hist_matched_jet_pt_parton_B_from_d", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_B_from_s =
    new TH1D("hist_matched_jet_pt_parton_B_from_s", "hist_matched_jet_pt_parton_B_from_s", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_B_from_c =
    new TH1D("hist_matched_jet_pt_parton_B_from_c", "hist_matched_jet_pt_parton_B_from_c", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_B_from_b =
    new TH1D("hist_matched_jet_pt_parton_B_from_b", "hist_matched_jet_pt_parton_B_from_b", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_B_from_t =
    new TH1D("hist_matched_jet_pt_parton_B_from_t", "hist_matched_jet_pt_parton_B_from_t", 100, 0.0, 500.0);
TH1D *hist_matched_jet_pt_parton_B_from_g =
    new TH1D("hist_matched_jet_pt_parton_B_from_g", "hist_matched_jet_pt_parton_B_from_g", 100, 0.0, 500.0);

// histograms for Jet Energy Scale (JES)
// double jetptbin[18] = {30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,140.0,160.0,180.0,210.0,250.0,300.0,380.0,500.0}; //just to remember

int bins3D_jetrecovsgen[3] = {50, 50, multbinsize - 1};
double xmin3D_jetrecovsgen[3] = {0.0, 0, 0};
double xmax3D_jetrecovsgen[3] = {500.0, 500, (double)multbinsize - 1};

// JES JER 2D Histogram
THnSparseD *hist_matchcorrectedpt_vs_refpt = new THnSparseD("hist_matchcorrectedpt_vs_refpt",
                                                            "hist_matchcorrectedpt_vs_refpt",
                                                            3,
                                                            bins3D_jetrecovsgen,
                                                            xmin3D_jetrecovsgen,
                                                            xmax3D_jetrecovsgen);
THnSparseD *hist_matchcorrectedpt_vs_refpt_weighted = new THnSparseD("hist_matchcorrectedpt_vs_refpt_weighted",
                                                                     "hist_matchcorrectedpt_vs_refpt_weighted",
                                                                     3,
                                                                     bins3D_jetrecovsgen,
                                                                     xmin3D_jetrecovsgen,
                                                                     xmax3D_jetrecovsgen);

// JES JER Manual Matching 2D histograms
THnSparseD *hist_matchcorrectedpt_vs_genpt = new THnSparseD("hist_matchcorrectedpt_vs_genpt",
                                                            "hist_matchcorrectedpt_vs_genpt",
                                                            3,
                                                            bins3D_jetrecovsgen,
                                                            xmin3D_jetrecovsgen,
                                                            xmax3D_jetrecovsgen);
THnSparseD *hist_matchcorrectedpt_vs_genpt_weighted = new THnSparseD("hist_matchcorrectedpt_vs_genpt_weighted",
                                                                     "hist_matchcorrectedpt_vs_genpt_weighted",
                                                                     3,
                                                                     bins3D_jetrecovsgen,
                                                                     xmin3D_jetrecovsgen,
                                                                     xmax3D_jetrecovsgen);

int bins3D_jes[4] = {500, 50, multbinsize - 1, 80};
double xmin3D_jes[4] = {0.0, 0, 0, -4.0};
double xmax3D_jes[4] = {5.0, 500, (double)multbinsize - 1, 4.0};

int bins3D_jer[3] = {500, 50, multbinsize - 1};
double xmin3D_jer[3] = {-5.0, 0, 0};
double xmax3D_jer[3] = {5.0, 500, (double)multbinsize - 1};
// JES Ratio Histogram
THnSparseD *hist_jes_ratio_matchedcorrectedpt_refpt = new THnSparseD("hist_jes_ratio_matchedcorrectedpt_refpt",
                                                                     "hist_jes_ratio_matchedcorrectedpt_refpt",
                                                                     4,
                                                                     bins3D_jes,
                                                                     xmin3D_jes,
                                                                     xmax3D_jes);
THnSparseD *hist_jes_ratio_matchedcorrectedpt_refpt_weighted =
    new THnSparseD("hist_jes_ratio_matchedcorrectedpt_refpt_weighted",
                   "hist_jes_ratio_matchedcorrectedpt_refpt_weighted",
                   4,
                   bins3D_jes,
                   xmin3D_jes,
                   xmax3D_jes);

// JER ratio histograms
THnSparseD *hist_jer_ratio_matchedcorrectedpt_refpt = new THnSparseD("hist_jer_ratio_matchedcorrectedpt_refpt",
                                                                     "hist_jer_ratio_matchedcorrectedpt_refpt",
                                                                     3,
                                                                     bins3D_jer,
                                                                     xmin3D_jer,
                                                                     xmax3D_jer);
THnSparseD *hist_jer_ratio_matchedcorrectedpt_refpt_weighted =
    new THnSparseD("hist_jer_ratio_matchedcorrectedpt_refpt_weighted",
                   "hist_jer_ratio_matchedcorrectedpt_refpt_weighted",
                   3,
                   bins3D_jer,
                   xmin3D_jer,
                   xmax3D_jer);

// JES Manula Matching Ratio Histogram
THnSparseD *hist_jes_ratio_matchedcorrectedpt_genpt = new THnSparseD("hist_jes_ratio_matchedcorrectedpt_genpt",
                                                                     "hist_jes_ratio_matchedcorrectedpt_genpt",
                                                                     3,
                                                                     bins3D_jes,
                                                                     xmin3D_jes,
                                                                     xmax3D_jes);
THnSparseD *hist_jes_ratio_matchedcorrectedpt_genpt_weighted =
    new THnSparseD("hist_jes_ratio_matchedcorrectedpt_genpt_weighted",
                   "hist_jes_ratio_matchedcorrectedpt_genpt_weighted",
                   3,
                   bins3D_jes,
                   xmin3D_jes,
                   xmax3D_jes);

// JER Manula Matching ratio histograms
THnSparseD *hist_jer_ratio_matchedcorrectedpt_genpt = new THnSparseD("hist_jer_ratio_matchedcorrectedpt_genpt",
                                                                     "hist_jer_ratio_matchedcorrectedpt_genpt",
                                                                     3,
                                                                     bins3D_jer,
                                                                     xmin3D_jer,
                                                                     xmax3D_jer);
THnSparseD *hist_jer_ratio_matchedcorrectedpt_genpt_weighted =
    new THnSparseD("hist_jer_ratio_matchedcorrectedpt_genpt_weighted",
                   "hist_jer_ratio_matchedcorrectedpt_genpt_weighted",
                   3,
                   bins3D_jer,
                   xmin3D_jer,
                   xmax3D_jer);

TH1D *dR = new TH1D("dR", "dR", 100, 0.0, 10.0);
TH1D *dR_Manual = new TH1D("dR_Manual", "dR_Manual", 100, 0.0, 10.0);

int bins3D_jtrk[3] = {100, 100, 8};
double xmin3D_jtrk[3] = {0.0, 0, 0};
double xmax3D_jtrk[3] = {500.0, 100, 8};

THnSparseD *hist_jet_multiplicity =
    new THnSparseD("hist_jet_multiplicity", "hist_jet_multiplicity", 3, bins3D_jtrk, xmin3D_jtrk, xmax3D_jtrk);
THnSparseD *hist_jet_multiplicity_weighted = new THnSparseD(
    "hist_jet_multiplicity_weighted", "hist_jet_multiplicity_weighted", 3, bins3D_jtrk, xmin3D_jtrk, xmax3D_jtrk);

int bins4D_jtparton[4] = {100, 50, 64, 8};
double xmin4D_jtparton[4] = {0.0, -2.5, -TMath::Pi(), 0};
double xmax4D_jtparton[4] = {500.0, 2.5, TMath::Pi(), 8};
THnSparseD *hist_flavor_jets =
    new THnSparseD("hist_flavor_jets", "hist_flavor_jets", 4, bins4D_jtparton, xmin4D_jtparton, xmax4D_jtparton);
THnSparseD *hist_flavor_jets_weighted = new THnSparseD(
    "hist_flavor_jets_weighted", "hist_flavor_jets_weighted", 4, bins4D_jtparton, xmin4D_jtparton, xmax4D_jtparton);

TH2D *hist_hibin_mult = new TH2D("hist_hibin_mult", "hist_hibin_mult", 200, 0.0, 200.0, 3000, 0.0, 3000.0);
TH2D *hist_hibin_mult_weighted =
    new TH2D("hist_hibin_mult_weighted", "hist_hibin_mult_weighted", 200, 0.0, 200.0, 3000, 0.0, 3000.0);

TH1D *hist_hibin = new TH1D("hist_hibin", "hist_hibin", 200, 0., 200.);
TH1D *hist_hibin_weighted = new TH1D("hist_hibin_weighted", "hist_hibin_weighted", 200, 0., 200.);
TH1D *hist_ntrkoff = new TH1D("hist_ntrkoff", "hist_ntrkoff", 1000, 0.0, 5000.0);
TH1D *hist_ntrkoff_weighted = new TH1D("hist_ntrkoff_weighted", "hist_ntrkoff_weighted", 1000, 0.0, 5000.0);

TH1D *hist_ntrkoff_10 = new TH1D("hist_ntrkoff_10", "hist_ntrkoff_10", 500, 0.0, 500.0);
TH1D *hist_ntrkoff_10_weighted = new TH1D("hist_ntrkoff_10_weighted", "hist_ntrkoff_10_weighted", 500, 0.0, 500.0);

TH1D *hist_ntrkoff_60 = new TH1D("hist_ntrkoff_60", "hist_ntrkoff_60", 500, 0.0, 500.0);
TH1D *hist_ntrkoff_60_weighted = new TH1D("hist_ntrkoff_60_weighted", "hist_ntrkoff_60_weighted", 500, 0.0, 500.0);

TH1D *hist_ntrkoff_120 = new TH1D("hist_ntrkoff_120", "hist_ntrkoff_120", 500, 0.0, 500.0);
TH1D *hist_ntrkoff_120_weighted = new TH1D("hist_ntrkoff_120_weighted", "hist_ntrkoff_120_weighted", 500, 0.0, 500.0);

TH1D *hist_gen_multiplicity = new TH1D("hist_gen_multiplicity", "gen_multiplicity", 1000, 0.0, 5000.0);
TH1D *hist_gen_multiplicity_weighted =
    new TH1D("hist_gen_multiplicity_weighted", "gen_multiplicity_weighted", 1000, 0.0, 5000.0);

TH1D *hist_gen_ntrkoff_10 = new TH1D("hist_gen_ntrkoff_10", "hist_gen_ntrkoff_10", 500, 0.0, 500.0);
TH1D *hist_gen_ntrkoff_10_weighted =
    new TH1D("hist_gen_ntrkoff_10_weighted", "hist_gen_ntrkoff_10_weighted", 500, 0.0, 500.0);

TH1D *hist_gen_ntrkoff_185 = new TH1D("hist_gen_ntrkoff_185", "hist_gen_ntrkoff_185", 500, 0.0, 500.0);
TH1D *hist_gen_ntrkoff_185_weighted =
    new TH1D("hist_gen_ntrkoff_185_weighted", "hist_gen_ntrkoff_185_weighted", 500, 0.0, 500.0);

TH1D *hist_gen_ntrkoff_250 = new TH1D("hist_gen_ntrkoff_250", "hist_gen_ntrkoff_250", 500, 0.0, 500.0);
TH1D *hist_gen_ntrkoff_250_weighted =
    new TH1D("hist_gen_ntrkoff_250_weighted", "hist_gen_ntrkoff_250_weighted", 500, 0.0, 500.0);

TH2D *hist_ntrkoff_eff = new TH2D("hist_ntrkoff_eff", "hist_ntrkoff_eff", 1000, 0.0, 5000.0, 1000, 0.0, 5000.0);
TH2D *hist_ntrkoff_eff_weighted =
    new TH2D("hist_ntrkoff_eff_weighted", "hist_ntrkoff_eff_weighted", 1000, 0.0, 5000.0, 1000, 0.0, 5000.0);

float pt_edges[44] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45,                       // 6
                      0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, // 10
                      1.0, 1.05, 1.1, 1.15, 1.2,                             // 5
                      1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,                // 8
                      2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0,             // 8
                      20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0};            // 7
float eta_edges[19] = {
    -2.4, -2.0, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4};
float cent_edges[18] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0};

TH2D *hist_eta_pt_efficiency_185 = new TH2D("hist_eta_pt_efficiency_185", "hist_eta_pt_efficiency_185", 18, eta_edges, 43, pt_edges);
TH2D *hist_eta_pt_efficiency_gen_185 = new TH2D("hist_eta_pt_efficiency_gen_185", "hist_eta_pt_efficiency_gen_185", 18, eta_edges, 43, pt_edges);

TH2D *hist_eta_pt_efficiency_185_weighted = new TH2D("hist_eta_pt_efficiency_185_weighted", "hist_eta_pt_efficiency_185_weighted", 18, eta_edges, 43, pt_edges);
TH2D *hist_eta_pt_efficiency_gen_185_weighted = new TH2D("hist_eta_pt_efficiency_gen_185_weighted", "hist_eta_pt_efficiency_gen_185_weighted", 18, eta_edges, 43, pt_edges);

TH2D *hist_eta_pt_efficiency_250 = new TH2D("hist_eta_pt_efficiency_250", "hist_eta_pt_efficiency_250", 18, eta_edges, 43, pt_edges);
TH2D *hist_eta_pt_efficiency_gen_250 = new TH2D("hist_eta_pt_efficiency_gen_250", "hist_eta_pt_efficiency_gen_250", 18, eta_edges, 43, pt_edges);

TH2D *hist_eta_pt_efficiency_250_weighted = new TH2D("hist_eta_pt_efficiency_250_weighted", "hist_eta_pt_efficiency_250_weighted", 18, eta_edges, 43, pt_edges);
TH2D *hist_eta_pt_efficiency_gen_250_weighted = new TH2D("hist_eta_pt_efficiency_gen_250_weighted", "hist_eta_pt_efficiency_gen_250_weighted", 18, eta_edges, 43, pt_edges);

TH3D *hist_reco_matched_eta_pt = new TH3D("hist_reco_matched_eta_pt", "hist_reco_matched_eta_pt", 18, eta_edges, 43, pt_edges, 17, cent_edges);
TH3D *hist_gen_matched_eta_pt = new TH3D("hist_gen_matched_eta_pt", "hist_gen_matched_eta_pt", 18, eta_edges, 43, pt_edges, 17, cent_edges);
TH3D *hist_fake_matched_eta_pt = new TH3D("hist_fake_matched_eta_pt", "hist_fake_matched_eta_pt", 18, eta_edges, 43, pt_edges, 17, cent_edges);

TH2D *hist_mult_10_60 = new TH2D("hist_mult_10_60", "hist_mult_10_60", 400, 0, 400, 400, 0, 400);
TH2D *hist_mult_60_120 = new TH2D("hist_mult_60_120", "hist_mult_60_120", 400, 0, 400, 400, 0, 400);
TH2D *hist_mult_120_185 = new TH2D("hist_mult_120_185", "hist_mult_120_185", 400, 0, 400, 400, 0, 400);
TH2D *hist_mult_185_250 = new TH2D("hist_mult_185_250", "hist_mult_185_250", 400, 0, 400, 400, 0, 400);
TH2D *hist_mult_250_400 = new TH2D("hist_mult_250_400", "hist_mult_250_400", 400, 0, 400, 400, 0, 400);
// Evaluate uncertainties correctly at ROOT
void sw2()
{
    Nevents->Sumw2();
    Nev_recoreco->Sumw2();
    Nev_recoreco_lead->Sumw2();
    Nev_recoreco_subl->Sumw2();
    Nev_recogen->Sumw2();
    Nev_recogen_lead->Sumw2();
    Nev_recogen_subl->Sumw2();
    Nev_genreco->Sumw2();
    Nev_genreco_lead->Sumw2();
    Nev_genreco_subl->Sumw2();
    Nev_gengen->Sumw2();
    Nev_gengen_lead->Sumw2();
    Nev_gengen_subl->Sumw2();
    multiplicity->Sumw2();
    multiplicity_weighted->Sumw2();
    vzhist->Sumw2();
    vzhist_weighted->Sumw2();
    pthathist->Sumw2();
    pthathist_weighted->Sumw2();
    hist_reco_trk->Sumw2();
    hist_reco_trk_corr->Sumw2();
    hist_reco_trk_weighted->Sumw2();
    hist_gen_trk->Sumw2();
    hist_gen_trk_weighted->Sumw2();
    hist_reco_jet_weighted_nocut->Sumw2();
    hist_reco_jet->Sumw2();
    hist_reco_jet_corr->Sumw2();
    hist_reco_jet_weighted->Sumw2();
    hist_reco_jet_corr_weighted->Sumw2();
    hist_reco_leadjet_pt_nocut->Sumw2();
    hist_reco_leadjet_pt_nocut_weighted->Sumw2();
    hist_reco_subljet_pt_nocut->Sumw2();
    hist_reco_subljet_pt_nocut_weighted->Sumw2();
    hist_reco_leadjet->Sumw2();
    hist_reco_leadjet_weighted->Sumw2();
    hist_reco_subljet->Sumw2();
    hist_reco_subljet_weighted->Sumw2();
    hist_gen_jet_weighted_nocut->Sumw2();
    hist_gen_jet->Sumw2();
    hist_gen_jet_weighted->Sumw2();
    hist_gen_leadjet_pt_nocut->Sumw2();
    hist_gen_leadjet_pt_nocut_weighted->Sumw2();
    hist_gen_subljet_pt_nocut->Sumw2();
    hist_gen_subljet_pt_nocut_weighted->Sumw2();
    hist_gen_leadjet->Sumw2();
    hist_gen_leadjet_weighted->Sumw2();
    hist_gen_subljet->Sumw2();
    hist_gen_subljet_weighted->Sumw2();
    hist_correlation_signal_jet_reco_track_reco->Sumw2();
    hist_correlation_rotation_jet_reco_track_reco->Sumw2();
    hist_correlation_mixing_jet_reco_track_reco->Sumw2();
    hist_correlation_signal_lead_jet_reco_track_reco->Sumw2();
    hist_correlation_rotation_lead_jet_reco_track_reco->Sumw2();
    hist_correlation_mixing_lead_jet_reco_track_reco->Sumw2();
    hist_correlation_signal_subl_jet_reco_track_reco->Sumw2();
    hist_correlation_rotation_subl_jet_reco_track_reco->Sumw2();
    hist_correlation_mixing_subl_jet_reco_track_reco->Sumw2();
    hist_correlation_signal_jet_reco_track_gen->Sumw2();
    hist_correlation_rotation_jet_reco_track_gen->Sumw2();
    hist_correlation_mixing_jet_reco_track_gen->Sumw2();
    hist_correlation_signal_lead_jet_reco_track_gen->Sumw2();
    hist_correlation_rotation_lead_jet_reco_track_gen->Sumw2();
    hist_correlation_mixing_lead_jet_reco_track_gen->Sumw2();
    hist_correlation_signal_subl_jet_reco_track_gen->Sumw2();
    hist_correlation_rotation_subl_jet_reco_track_gen->Sumw2();
    hist_correlation_mixing_subl_jet_reco_track_gen->Sumw2();
    hist_correlation_signal_jet_gen_track_reco->Sumw2();
    hist_correlation_rotation_jet_gen_track_reco->Sumw2();
    hist_correlation_mixing_jet_gen_track_reco->Sumw2();
    hist_correlation_signal_lead_jet_gen_track_reco->Sumw2();
    hist_correlation_rotation_lead_jet_gen_track_reco->Sumw2();
    hist_correlation_mixing_lead_jet_gen_track_reco->Sumw2();
    hist_correlation_signal_subl_jet_gen_track_reco->Sumw2();
    hist_correlation_rotation_subl_jet_gen_track_reco->Sumw2();
    hist_correlation_mixing_subl_jet_gen_track_reco->Sumw2();
    hist_correlation_signal_jet_gen_track_gen->Sumw2();
    hist_correlation_rotation_jet_gen_track_gen->Sumw2();
    hist_correlation_mixing_jet_gen_track_gen->Sumw2();
    hist_correlation_signal_lead_jet_gen_track_gen->Sumw2();
    hist_correlation_rotation_lead_jet_gen_track_gen->Sumw2();
    hist_correlation_mixing_lead_jet_gen_track_gen->Sumw2();
    hist_correlation_signal_subl_jet_gen_track_gen->Sumw2();
    hist_correlation_rotation_subl_jet_gen_track_gen->Sumw2();
    hist_correlation_mixing_subl_jet_gen_track_gen->Sumw2();
    hist_correlation_signal_subg0_jet_reco_track_reco->Sumw2();
    hist_correlation_signal_subg0_jet_reco_track_gen->Sumw2();
    hist_correlation_signal_subg0_jet_gen_track_reco->Sumw2();
    hist_correlation_signal_subg0_jet_gen_track_gen->Sumw2();
    hist_correlation_signal_subg0_lead_jet_reco_track_reco->Sumw2();
    hist_correlation_signal_subg0_lead_jet_reco_track_gen->Sumw2();
    hist_correlation_signal_subg0_lead_jet_gen_track_reco->Sumw2();
    hist_correlation_signal_subg0_lead_jet_gen_track_gen->Sumw2();
    hist_correlation_signal_subg0_subl_jet_reco_track_reco->Sumw2();
    hist_correlation_signal_subg0_subl_jet_reco_track_gen->Sumw2();
    hist_correlation_signal_subg0_subl_jet_gen_track_reco->Sumw2();
    hist_correlation_signal_subg0_subl_jet_gen_track_gen->Sumw2();
    hist_jet_from_reco_reco_sig->Sumw2();
    hist_jet_from_reco_gen_sig->Sumw2();
    hist_jet_from_gen_reco_sig->Sumw2();
    hist_jet_from_gen_gen_sig->Sumw2();
    hist_jet_from_reco_reco_mix->Sumw2();
    hist_jet_from_reco_gen_mix->Sumw2();
    hist_jet_from_gen_reco_mix->Sumw2();
    hist_jet_from_gen_gen_mix->Sumw2();
    hist_lead_jet_from_reco_reco_sig->Sumw2();
    hist_lead_jet_from_reco_gen_sig->Sumw2();
    hist_lead_jet_from_gen_reco_sig->Sumw2();
    hist_lead_jet_from_gen_gen_sig->Sumw2();
    hist_lead_jet_from_reco_reco_mix->Sumw2();
    hist_lead_jet_from_reco_gen_mix->Sumw2();
    hist_lead_jet_from_gen_reco_mix->Sumw2();
    hist_lead_jet_from_gen_gen_mix->Sumw2();
    hist_subl_jet_from_reco_reco_sig->Sumw2();
    hist_subl_jet_from_reco_gen_sig->Sumw2();
    hist_subl_jet_from_gen_reco_sig->Sumw2();
    hist_subl_jet_from_gen_gen_sig->Sumw2();
    hist_subl_jet_from_reco_reco_mix->Sumw2();
    hist_subl_jet_from_reco_gen_mix->Sumw2();
    hist_subl_jet_from_gen_reco_mix->Sumw2();
    hist_subl_jet_from_gen_gen_mix->Sumw2();
    hist_trk_from_reco_reco_sig->Sumw2();
    hist_trk_from_reco_gen_sig->Sumw2();
    hist_trk_from_gen_reco_sig->Sumw2();
    hist_trk_from_gen_gen_sig->Sumw2();
    hist_trk_from_reco_reco_mix->Sumw2();
    hist_trk_from_reco_gen_mix->Sumw2();
    hist_trk_from_gen_reco_mix->Sumw2();
    hist_trk_from_gen_gen_mix->Sumw2();
    hist_LJ_trk_from_reco_reco_sig->Sumw2();
    hist_LJ_trk_from_reco_gen_sig->Sumw2();
    hist_LJ_trk_from_gen_reco_sig->Sumw2();
    hist_LJ_trk_from_gen_gen_sig->Sumw2();
    hist_LJ_trk_from_reco_reco_mix->Sumw2();
    hist_LJ_trk_from_reco_gen_mix->Sumw2();
    hist_LJ_trk_from_gen_reco_mix->Sumw2();
    hist_LJ_trk_from_gen_gen_mix->Sumw2();
    hist_SLJ_trk_from_reco_reco_sig->Sumw2();
    hist_SLJ_trk_from_reco_gen_sig->Sumw2();
    hist_SLJ_trk_from_gen_reco_sig->Sumw2();
    hist_SLJ_trk_from_gen_gen_sig->Sumw2();
    hist_SLJ_trk_from_reco_reco_mix->Sumw2();
    hist_SLJ_trk_from_reco_gen_mix->Sumw2();
    hist_SLJ_trk_from_gen_reco_mix->Sumw2();
    hist_SLJ_trk_from_gen_gen_mix->Sumw2();
    hist_reco_lead_reco_subl_quench->Sumw2();
    hist_reco_lead_reco_subl_quench2pc->Sumw2();
    hist_gen_lead_gen_subl_quench->Sumw2();
    hist_gen_lead_gen_subl_quench2pc->Sumw2();
    hist_matched_jet_weighted_nocut->Sumw2();
    hist_matched_jet->Sumw2();
    hist_matched_jet_weighted->Sumw2();
    hist_matched_jet_pt_parton_from_u->Sumw2();
    hist_matched_jet_pt_parton_from_d->Sumw2();
    hist_matched_jet_pt_parton_from_s->Sumw2();
    hist_matched_jet_pt_parton_from_c->Sumw2();
    hist_matched_jet_pt_parton_from_b->Sumw2();
    hist_matched_jet_pt_parton_from_t->Sumw2();
    hist_matched_jet_pt_parton_from_g->Sumw2();
    hist_matched_jet_pt_parton_B_from_u->Sumw2();
    hist_matched_jet_pt_parton_B_from_d->Sumw2();
    hist_matched_jet_pt_parton_B_from_s->Sumw2();
    hist_matched_jet_pt_parton_B_from_c->Sumw2();
    hist_matched_jet_pt_parton_B_from_b->Sumw2();
    hist_matched_jet_pt_parton_B_from_t->Sumw2();
    hist_matched_jet_pt_parton_B_from_g->Sumw2();
    hist_matchcorrectedpt_vs_refpt->Sumw2();
    hist_matchcorrectedpt_vs_refpt_weighted->Sumw2();
    hist_matchcorrectedpt_vs_genpt->Sumw2();
    hist_matchcorrectedpt_vs_genpt_weighted->Sumw2();

    hist_jes_ratio_matchedcorrectedpt_refpt->Sumw2();
    hist_jes_ratio_matchedcorrectedpt_refpt_weighted->Sumw2();
    hist_jer_ratio_matchedcorrectedpt_refpt->Sumw2();
    hist_jer_ratio_matchedcorrectedpt_refpt_weighted->Sumw2();
    hist_jes_ratio_matchedcorrectedpt_genpt->Sumw2();
    hist_jes_ratio_matchedcorrectedpt_genpt_weighted->Sumw2();
    hist_jer_ratio_matchedcorrectedpt_genpt->Sumw2();
    hist_jer_ratio_matchedcorrectedpt_genpt_weighted->Sumw2();
    dR->Sumw2();
    dR_Manual->Sumw2();
    hist_jet_multiplicity->Sumw2();
    hist_jet_multiplicity_weighted->Sumw2();
    hist_flavor_jets_weighted->Sumw2();
    hist_flavor_jets->Sumw2();
    hist_hibin_mult->Sumw2();
    hist_hibin_mult_weighted->Sumw2();
    hist_hibin->Sumw2();
    hist_hibin_weighted->Sumw2();
    hist_ntrkoff->Sumw2();
    hist_ntrkoff_weighted->Sumw2();
    hist_ntrkoff_10->Sumw2();
    hist_ntrkoff_60->Sumw2();
    hist_ntrkoff_120->Sumw2();
    hist_ntrkoff_10_weighted->Sumw2();
    hist_ntrkoff_60_weighted->Sumw2();
    hist_ntrkoff_120_weighted->Sumw2();
    hist_gen_multiplicity->Sumw2();
    hist_gen_multiplicity_weighted->Sumw2();
    hist_gen_ntrkoff_10->Sumw2();
    hist_gen_ntrkoff_185->Sumw2();
    hist_gen_ntrkoff_250->Sumw2();
    hist_gen_ntrkoff_10_weighted->Sumw2();
    hist_gen_ntrkoff_185_weighted->Sumw2();
    hist_gen_ntrkoff_250_weighted->Sumw2();
    hist_ntrkoff_eff->Sumw2();
    hist_ntrkoff_eff_weighted->Sumw2();
    hist_eta_pt_efficiency_250->Sumw2();
    hist_eta_pt_efficiency_250_weighted->Sumw2();
    hist_eta_pt_efficiency_185->Sumw2();
    hist_eta_pt_efficiency_185_weighted->Sumw2();
    hist_eta_pt_efficiency_gen_185->Sumw2();
    hist_eta_pt_efficiency_gen_250->Sumw2();
    hist_eta_pt_efficiency_gen_185_weighted->Sumw2();
    hist_eta_pt_efficiency_gen_250_weighted->Sumw2();
    centrality_weighted->Sumw2();
    centrality->Sumw2();
    vzhist185_weighted->Sumw2();
    vzhist250_weighted->Sumw2();
    hist_reco_matched_eta_pt->Sumw2();
    hist_gen_matched_eta_pt->Sumw2();
    hist_fake_matched_eta_pt->Sumw2();
    hist_mult_250_400->Sumw2();
    hist_mult_185_250->Sumw2();
    hist_mult_120_185->Sumw2();
    hist_mult_60_120->Sumw2();
    hist_mult_10_60->Sumw2();
}

// write QA histograms
/*
--> Arguments
isMC: true for MC and false for Data
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_QA_hist(bool isMC, bool doleadsubl)
{
    Nevents->Write();
    Nev_recoreco->Write();
    if (doleadsubl)
    {
        Nev_recoreco_lead->Write();
        Nev_recoreco_subl->Write();
    }
    if (isMC)
    {
        Nev_recogen->Write();
        Nev_genreco->Write();
        Nev_gengen->Write();
        if (doleadsubl)
        {
            Nev_recogen_lead->Write();
            Nev_genreco_lead->Write();
            Nev_gengen_lead->Write();
            Nev_recogen_subl->Write();
            Nev_genreco_subl->Write();
            Nev_gengen_subl->Write();
        }
    }
    hist_mult_250_400->Write();
    hist_mult_185_250->Write();
    hist_mult_120_185->Write();
    hist_mult_60_120->Write();
    hist_mult_10_60->Write();
    multiplicity->Write();
    multiplicity_weighted->Write();
    hist_gen_multiplicity->Write();
    hist_gen_multiplicity_weighted->Write();
    vzhist->Write();
    vzhist_weighted->Write();
    vzhist185_weighted->Write();
    vzhist250_weighted->Write();
    hist_hibin_mult->Write();
    hist_hibin_mult_weighted->Write();
    hist_hibin->Write();
    hist_hibin_weighted->Write();
    hist_ntrkoff->Write();
    hist_ntrkoff_weighted->Write();
    hist_ntrkoff_10->Write();
    hist_ntrkoff_60->Write();
    hist_ntrkoff_120->Write();
    hist_ntrkoff_10_weighted->Write();
    hist_ntrkoff_60_weighted->Write();
    hist_ntrkoff_120_weighted->Write();
    hist_gen_ntrkoff_10->Write();
    hist_gen_ntrkoff_185->Write();
    hist_gen_ntrkoff_250->Write();
    hist_gen_ntrkoff_10_weighted->Write();
    hist_gen_ntrkoff_185_weighted->Write();
    hist_gen_ntrkoff_250_weighted->Write();
    hist_ntrkoff_eff->Write();
    hist_ntrkoff_eff_weighted->Write();

    hist_eta_pt_efficiency_250->Write();
    hist_eta_pt_efficiency_250_weighted->Write();
    hist_eta_pt_efficiency_185->Write();
    hist_eta_pt_efficiency_185_weighted->Write();
    hist_eta_pt_efficiency_gen_185->Write();
    hist_eta_pt_efficiency_gen_250->Write();
    hist_eta_pt_efficiency_gen_185_weighted->Write();
    hist_eta_pt_efficiency_gen_250_weighted->Write();

    centrality_weighted->Write();
    centrality->Write();

    hist_reco_matched_eta_pt->Write();
    hist_gen_matched_eta_pt->Write();
    hist_fake_matched_eta_pt->Write();

    if (isMC)
    {
        pthathist->Write();
        pthathist_weighted->Write();
    }
    // tracks
    // reco
    hist_reco_trk->Write();
    hist_reco_trk_corr->Write();
    hist_reco_trk_weighted->Write();

    // gen
    if (isMC)
    {
        hist_gen_trk->Write();
        hist_gen_trk_weighted->Write();
    }
    // jets
    // reco
    hist_reco_jet_weighted_nocut->Write();
    hist_reco_jet->Write();
    hist_reco_jet_corr->Write();
    hist_reco_jet_weighted->Write();
    hist_reco_jet_corr_weighted->Write();

    if (doleadsubl)
    {
        hist_reco_leadjet_pt_nocut->Write();
        hist_reco_leadjet_pt_nocut_weighted->Write();
        hist_reco_subljet_pt_nocut->Write();
        hist_reco_subljet_pt_nocut_weighted->Write();
        hist_reco_leadjet->Write();
        hist_reco_leadjet_weighted->Write();
        hist_reco_subljet->Write();
        hist_reco_subljet_weighted->Write();
    }
    if (isMC)
    {
        hist_gen_jet_weighted_nocut->Write();
        hist_gen_jet->Write();
        hist_gen_jet_weighted->Write();
        if (doleadsubl)
        {
            hist_gen_leadjet_pt_nocut->Write();
            hist_gen_leadjet_pt_nocut_weighted->Write();
            hist_gen_subljet_pt_nocut->Write();
            hist_gen_subljet_pt_nocut_weighted->Write();
            hist_gen_leadjet->Write();
            hist_gen_leadjet_weighted->Write();
            hist_gen_subljet->Write();
            hist_gen_subljet_weighted->Write();
        }
    }
}

// Reco-Reco correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_recoreco_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl)
{
    if (doinclusive)
    {
        hist_correlation_signal_jet_reco_track_reco->Write();
        hist_correlation_signal_subg0_jet_reco_track_reco->Write();
        hist_jet_from_reco_reco_sig->Write();
        hist_trk_from_reco_reco_sig->Write();
        if (rotation)
            hist_correlation_rotation_jet_reco_track_reco->Write();
        if (mixing)
        {
            hist_correlation_mixing_jet_reco_track_reco->Write();
            hist_jet_from_reco_reco_mix->Write();
            hist_trk_from_reco_reco_mix->Write();
        }
    }

    if (doleadsubl)
    {
        hist_correlation_signal_lead_jet_reco_track_reco->Write();
        hist_correlation_signal_subg0_lead_jet_reco_track_reco->Write();
        hist_lead_jet_from_reco_reco_sig->Write();
        hist_LJ_trk_from_reco_reco_sig->Write();
        if (rotation)
            hist_correlation_rotation_lead_jet_reco_track_reco->Write();
        if (mixing)
        {
            hist_correlation_mixing_lead_jet_reco_track_reco->Write();
            hist_lead_jet_from_reco_reco_mix->Write();
            hist_LJ_trk_from_reco_reco_mix->Write();
        }
        hist_correlation_signal_subl_jet_reco_track_reco->Write();
        hist_correlation_signal_subg0_subl_jet_reco_track_reco->Write();
        hist_subl_jet_from_reco_reco_sig->Write();
        hist_SLJ_trk_from_reco_reco_sig->Write();
        if (rotation)
            hist_correlation_rotation_subl_jet_reco_track_reco->Write();
        if (mixing)
        {
            hist_correlation_mixing_subl_jet_reco_track_reco->Write();
            hist_subl_jet_from_reco_reco_mix->Write();
            hist_SLJ_trk_from_reco_reco_mix->Write();
        }
    }
}

// Reco-Gen correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_recogen_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl)
{
    if (doinclusive)
    {
        hist_correlation_signal_jet_reco_track_gen->Write();
        hist_correlation_signal_subg0_jet_reco_track_gen->Write();
        hist_jet_from_reco_gen_sig->Write();
        hist_trk_from_reco_gen_sig->Write();
        if (rotation)
            hist_correlation_rotation_jet_reco_track_gen->Write();
        if (mixing)
        {
            hist_correlation_mixing_jet_reco_track_gen->Write();
            hist_jet_from_reco_gen_mix->Write();
            hist_trk_from_reco_gen_mix->Write();
        }
    }

    if (doleadsubl)
    {
        hist_correlation_signal_lead_jet_reco_track_gen->Write();
        hist_correlation_signal_subg0_lead_jet_reco_track_gen->Write();
        hist_lead_jet_from_reco_gen_sig->Write();
        hist_LJ_trk_from_reco_gen_sig->Write();
        if (rotation)
            hist_correlation_rotation_lead_jet_reco_track_gen->Write();
        if (mixing)
        {
            hist_correlation_mixing_lead_jet_reco_track_gen->Write();
            hist_lead_jet_from_reco_gen_mix->Write();
            hist_LJ_trk_from_reco_gen_mix->Write();
        }
        hist_correlation_signal_subl_jet_reco_track_gen->Write();
        hist_correlation_signal_subg0_subl_jet_reco_track_gen->Write();
        hist_subl_jet_from_reco_gen_sig->Write();
        hist_SLJ_trk_from_reco_gen_sig->Write();
        if (rotation)
            hist_correlation_rotation_subl_jet_reco_track_gen->Write();
        if (mixing)
        {
            hist_correlation_mixing_subl_jet_reco_track_gen->Write();
            hist_subl_jet_from_reco_gen_mix->Write();
            hist_SLJ_trk_from_reco_gen_mix->Write();
        }
    }
}

// Gen-Reco correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_genreco_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl)
{
    if (doinclusive)
    {
        hist_correlation_signal_jet_gen_track_reco->Write();
        hist_correlation_signal_subg0_jet_gen_track_reco->Write();
        hist_jet_from_gen_reco_sig->Write();
        hist_trk_from_gen_reco_sig->Write();
        if (rotation)
            hist_correlation_rotation_jet_gen_track_reco->Write();
        if (mixing)
        {
            hist_correlation_mixing_jet_gen_track_reco->Write();
            hist_jet_from_gen_reco_mix->Write();
            hist_trk_from_gen_reco_mix->Write();
        }
    }

    if (doleadsubl)
    {
        hist_correlation_signal_lead_jet_gen_track_reco->Write();
        hist_correlation_signal_subg0_lead_jet_gen_track_reco->Write();
        hist_lead_jet_from_gen_reco_sig->Write();
        hist_LJ_trk_from_gen_reco_sig->Write();
        if (rotation)
            hist_correlation_rotation_lead_jet_gen_track_reco->Write();
        if (mixing)
        {
            hist_correlation_mixing_lead_jet_gen_track_reco->Write();
            hist_lead_jet_from_gen_reco_mix->Write();
            hist_LJ_trk_from_gen_reco_mix->Write();
        }
        hist_correlation_signal_subl_jet_gen_track_reco->Write();
        hist_correlation_signal_subg0_subl_jet_gen_track_reco->Write();
        hist_subl_jet_from_gen_reco_sig->Write();
        hist_SLJ_trk_from_gen_reco_sig->Write();
        if (rotation)
            hist_correlation_rotation_subl_jet_gen_track_reco->Write();
        if (mixing)
        {
            hist_correlation_mixing_subl_jet_gen_track_reco->Write();
            hist_subl_jet_from_gen_reco_mix->Write();
            hist_SLJ_trk_from_gen_reco_mix->Write();
        }
    }
}

// Gen-Gen correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_gengen_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl)
{
    if (doinclusive)
    {
        hist_correlation_signal_jet_gen_track_gen->Write();
        hist_correlation_signal_subg0_jet_gen_track_gen->Write();
        hist_jet_from_gen_gen_sig->Write();
        hist_trk_from_gen_gen_sig->Write();
        if (rotation)
            hist_correlation_rotation_jet_gen_track_gen->Write();
        if (mixing)
        {
            hist_correlation_mixing_jet_gen_track_gen->Write();
            hist_jet_from_gen_gen_mix->Write();
            hist_trk_from_gen_gen_mix->Write();
        }
    }

    if (doleadsubl)
    {
        hist_correlation_signal_lead_jet_gen_track_gen->Write();
        hist_correlation_signal_subg0_lead_jet_gen_track_gen->Write();
        hist_lead_jet_from_gen_gen_sig->Write();
        hist_LJ_trk_from_gen_gen_sig->Write();
        if (rotation)
            hist_correlation_rotation_lead_jet_gen_track_gen->Write();
        if (mixing)
        {
            hist_correlation_mixing_lead_jet_gen_track_gen->Write();
            hist_lead_jet_from_gen_gen_mix->Write();
            hist_LJ_trk_from_gen_gen_mix->Write();
        }
        hist_correlation_signal_subl_jet_gen_track_gen->Write();
        hist_correlation_signal_subg0_subl_jet_gen_track_gen->Write();
        hist_subl_jet_from_gen_gen_sig->Write();
        hist_SLJ_trk_from_gen_gen_sig->Write();
        if (rotation)
            hist_correlation_rotation_subl_jet_gen_track_gen->Write();
        if (mixing)
        {
            hist_correlation_mixing_subl_jet_gen_track_gen->Write();
            hist_subl_jet_from_gen_gen_mix->Write();
            hist_SLJ_trk_from_gen_gen_mix->Write();
        }
    }
}

// Jet quenching histograms
/*
--> Arguments
isMC: true for MC and false for Data
*/
void w_jetquenching_hist(bool isMC)
{
    hist_reco_lead_reco_subl_quench->Write();
    hist_reco_lead_reco_subl_quench2pc->Write();
    if (isMC)
    {
        hist_gen_lead_gen_subl_quench->Write();
        hist_gen_lead_gen_subl_quench2pc->Write();
    }
}

// Matched and parton jet spectra histograms
void w_QA_parton_hist()
{
    hist_matched_jet_weighted_nocut->Write();
    hist_matched_jet->Write();
    hist_matched_jet_weighted->Write();
    hist_matched_jet_pt_parton_from_u->Write();
    hist_matched_jet_pt_parton_from_d->Write();
    hist_matched_jet_pt_parton_from_s->Write();
    hist_matched_jet_pt_parton_from_c->Write();
    hist_matched_jet_pt_parton_from_b->Write();
    hist_matched_jet_pt_parton_from_t->Write();
    hist_matched_jet_pt_parton_from_g->Write();
    hist_matched_jet_pt_parton_B_from_u->Write();
    hist_matched_jet_pt_parton_B_from_d->Write();
    hist_matched_jet_pt_parton_B_from_s->Write();
    hist_matched_jet_pt_parton_B_from_c->Write();
    hist_matched_jet_pt_parton_B_from_b->Write();
    hist_matched_jet_pt_parton_B_from_t->Write();
    hist_matched_jet_pt_parton_B_from_g->Write();
    hist_jet_multiplicity->Write();
    hist_jet_multiplicity_weighted->Write();
    hist_flavor_jets_weighted->Write();
    hist_flavor_jets->Write();
}

// JES histograms
void w_jes_hist()
{
    hist_matchcorrectedpt_vs_refpt->Write();
    hist_matchcorrectedpt_vs_refpt_weighted->Write();
    hist_matchcorrectedpt_vs_genpt->Write();
    hist_matchcorrectedpt_vs_genpt_weighted->Write();
    hist_jes_ratio_matchedcorrectedpt_refpt->Write();
    hist_jes_ratio_matchedcorrectedpt_refpt_weighted->Write();
    hist_jer_ratio_matchedcorrectedpt_refpt->Write();
    hist_jer_ratio_matchedcorrectedpt_refpt_weighted->Write();
    hist_jes_ratio_matchedcorrectedpt_genpt->Write();
    hist_jes_ratio_matchedcorrectedpt_genpt_weighted->Write();
    hist_jer_ratio_matchedcorrectedpt_genpt->Write();
    hist_jer_ratio_matchedcorrectedpt_genpt_weighted->Write();
    dR->Write();
    dR_Manual->Write();
}
