#include "call_libraries.h"            // call libraries from ROOT and C++
#include "trk_efficiency_correction.h" // track efficiency correction
#include "weights.h"                   // weights applied
// #include "JECorr.h"
/*
Find Ntrk offline -> updated for all systems (and easy to update for future systems)
The Ntrk offline is a definition with specific cuts (we should not change it). The track systematics must be applied using the input_variables.h!
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
hp: track high purity workflow
pterr: track pT uncertainty
dcaxy: track DCA in the transverse plane
dcaxyerr: track DCA in the transverse plane uncertainty
dcaz: track DCA in the longitudinal plane
dcazerr: track DCA in the longitudinal plane uncertainty
chi2: track chi2 of reconstruction
ndof: track number of degrees of freedom reconstruction
nlayer: track number of layers with measurements
nhits: track number of hits with measurements
algo: track MVA algorith step
mva: track MVA algorith value [-1,1]
*/
float multiplicity_weight_60(int ntrkoff_)
{
  TF1 *mult_weight_60 = new TF1("mult_weight_60", "0.574258 * TMath::Exp(-0.5 * TMath::Power(((x - 29.4429 ) / 5.55299 ), 2)) + 0.332276 * TMath::Exp(-0.5 * TMath::Power(((x - 22.9841 ) / 4.47453 ), 2)) + 0.647315 * TMath::Exp(-0.5 * TMath::Power(((x - 38.6441 ) /7.57712 ), 2))");
  float multweight_4 = mult_weight_60->Eval((float)ntrkoff_);
  return multweight_4;
}
float multiplicity_weight_120(int ntrkoff_)
{
  TF1 *mult_weight_120 = new TF1("mult_weight_120", "0.753348 * TMath::Exp(-0.5 * TMath::Power(((x - 50.1855 ) / 6.819 ), 2)) + 0.479375 * TMath::Exp(-0.5 * TMath::Power(((x - 61.4951 ) / 9.62468 ), 2)) ");
  float multweight_3 = mult_weight_120->Eval((float)ntrkoff_);
  return multweight_3;
}
float multiplicity_weight_185(int ntrkoff_)
{
  TF1 *mult_weight_185 = new TF1("mult_weight_185", "0.535919 * TMath::Exp(-0.5 * TMath::Power(((x - 75.9566 ) / 7.55684 ), 2)) + 0.0815166 * TMath::Exp(-0.5 * TMath::Power(((x - 91.8459 ) / 1.25429 ), 2)) + 0.470442 * TMath::Exp(-0.5 * TMath::Power(((x - 83.5405 ) /12.1473 ), 2))");
  float multweight_1 = mult_weight_185->Eval((float)ntrkoff_);
  return multweight_1;
}

float multiplicity_weight_250(int ntrkoff_)
{
  TF1 *mult_weight_250 = new TF1("mult_weight_250", "0.818452 * TMath::Exp(-0.5 * TMath::Power(((x - 105.573) / 11.3248 ), 2)) + 0.0729919 * TMath::Exp(-0.5 * TMath::Power(((x - 124.524 ) /  15.3399 ), 2))  -0.131896 * TMath::Exp(-0.5 * TMath::Power(((x - 114.532 ) /1.80402 ), 2))");
  float multweight_2 = mult_weight_250->Eval((float)ntrkoff_);
  return multweight_2;
}

float vz_weight_10(float vzbin)
{
  TF1 *Vz_weight_10 = new TF1("Vz_weight_10", "-6.43079e-06 * TMath::Power(x, 4)  + 8.28848e-06 * TMath::Power(x, 3) - 7.73088e-04 * TMath::Power(x, 2) + 6.30374e-03 * TMath::Power(x, 1) + 0.97476 ");
  float vzweight_10 = Vz_weight_185->Eval(vzbin);
  return vzweight_10;
}

float vz_weight_185(float vzbin)
{
  TF1 *Vz_weight_185 = new TF1("Vz_weight_185", "-1.23656e-05 * TMath::Power(x, 4)  + 2.04872e-05 * TMath::Power(x, 3)  -2.04554e-03 * TMath::Power(x, 2) + 6.71352e-03 * TMath::Power(x, 1) + 9.48790e-01");
  float vzweight_185 = Vz_weight_185->Eval(vzbin);
  return vzweight_185;
}

float vz_weight_250(float vzbin)
{
  TF1 *Vz_weight_250 = new TF1("Vz_weight_250", "-2.81671e-06 * TMath::Power(x, 4) + 2.28879e-05* TMath::Power(x, 3)  -4.04633e-03 * TMath::Power(x, 2) - 1.92669e-03 * TMath::Power(x, 1) +  9.66937e-01");
  float vzweight_250 = Vz_weight_250->Eval(vzbin);
  return vzweight_250;
}

int get_Ntrkoff_new(TString col_sys, int col_energy, int yearofdatataking, int size, std::vector<float> *eta, std::vector<float> *pt, std::vector<char> *charge, std::vector<bool> *hp, std::vector<float> *pterr, std::vector<float> *dcaxy, std::vector<float> *dcaxyerr, std::vector<float> *dcaz, std::vector<float> *dcazerr)
{
  int Ntrk_off = 0;
  for (int ii = 0; ii < size; ii++)
  {
    if (fabs(eta->at(ii)) >= 2.4)
      continue;
    if (fabs(charge->at(ii)) == 0)
      continue;
    if (hp->at(ii) == false)
      continue;
    if (fabs(pterr->at(ii) / pt->at(ii)) >= 0.1)
      continue;
    if (fabs(dcaxy->at(ii) / dcaxyerr->at(ii)) >= 3.0)
      continue;
    if (fabs(dcaz->at(ii) / dcazerr->at(ii)) >= 3.0)
      continue;
    if (col_sys == "PbPb" && col_energy == 5020 && yearofdatataking == 2018)
    {
      if (pt->at(ii) <= 1.0)
        continue;
    }
    Ntrk_off = Ntrk_off + 1;
  }
  return Ntrk_off;
}

int get_Ntrkoff_old(TString col_sys, int col_energy, int yearofdatataking, int size, std::vector<float> *eta, std::vector<float> *pt, std::vector<char> *charge, std::vector<bool> *hp, std::vector<float> *pterr, std::vector<float> *dcaxy, std::vector<float> *dcaxyerr, std::vector<float> *dcaz, std::vector<float> *dcazerr, std::vector<float> *chi2, std::vector<char> *nlayer, std::vector<char> *nhits, std::vector<float> *pfEcal, std::vector<float> *pfHcal)
{
  int Ntrk_off = 0;
  for (int ii = 0; ii < size; ii++)
  {
    if (fabs(eta->at(ii)) >= 2.4)
      continue;
    if (fabs(charge->at(ii)) == 0)
      continue;
    if (hp->at(ii) == false)
      continue;
    if (fabs(pterr->at(ii) / pt->at(ii)) >= 0.1)
      continue;
    if (fabs(dcaxy->at(ii) / dcaxyerr->at(ii)) >= 3.0)
      continue;
    if (fabs(dcaz->at(ii) / dcazerr->at(ii)) >= 3.0)
      continue;
    double calomatching = ((pfEcal->at(ii) + pfHcal->at(ii)) / cosh(eta->at(ii))) / pt->at(ii);
    if (col_sys == "PbPb" && col_energy == 5020 && yearofdatataking == 2018)
    {
      if (pt->at(ii) <= 0.5)
        continue;
      if ((chi2->at(ii) / nlayer->at(ii)) >= 0.18)
        continue;
      if (nhits->at(ii) < 11)
        continue;
      if (pt->at(ii) > 20.0)
      {
        if (calomatching <= 0.5)
          continue;
      } // is this applicable in pp or pPb?
    }
    Ntrk_off = Ntrk_off + 1;
  }
  return Ntrk_off;
}

int gen_ntrkoff_new(std::vector<float> *gentrkpt, std::vector<float> *gentrketa, std::vector<int> *gentrkcharge)
{
  float ntrkoff = 0;
  for (int i = 0; i < (int)gentrkpt->size(); i++)
  {
    if (gentrkpt->at(i) <= 0.4)
      continue;
    if (gentrkcharge->at(i) == 0)
      continue;
    if (fabs(gentrketa->at(i)) > 2.4)
      continue;

    ntrkoff = ntrkoff + 1;
  }
  return ntrkoff;
}
/*
Calculate jet Aj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float asymmetry(float pt_leading, float pt_subleading)
{
  float Avariable = (pt_leading - pt_subleading) / (pt_leading + pt_subleading);
  return Avariable;
}

/*
Calculate jet Xj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float xjvar(float pt_leading, float pt_subleading)
{
  float XJvariable = pt_subleading / pt_leading;
  return XJvariable;
}

/*
Calculate Delta Eta
--> Arguments
eta1: eta of first object
eta2: eta of second object
*/
float deltaeta(float eta1, float eta2)
{
  float deltaEta = (eta1 - eta2);
  return deltaEta;
}

/*
Calculate Delta Phi
--> Arguments
phi1: eta of first object
phi2: eta of second object
*/
float deltaphi(float phi1, float phi2)
{
  float deltaPhi = (phi1 - phi2);
  while (deltaPhi > (TMath::Pi()))
  {
    deltaPhi += -2 * TMath::Pi();
  }
  while (deltaPhi < (-1.0 * TMath::Pi()))
  {
    deltaPhi += 2 * TMath::Pi();
  }
  return deltaPhi;
}

/*
Calculate Delta Phi in the range [-pi/2 , 3/2 pi]
--> Arguments
phi1: eta of first object
phi2: eta of second object
*/
float deltaphi2PC(float phi1, float phi2)
{
  float deltaPhi = (phi1 - phi2);
  if (deltaPhi > TMath::Pi())
    deltaPhi = deltaPhi - 2. * TMath::Pi();
  if (deltaPhi <= -TMath::Pi() / 2.)
    deltaPhi = deltaPhi + 2. * TMath::Pi();
  return deltaPhi;
}

/*
Calculate delta R (distance)
--> Arguments
eta1: eta of first object
phi1: eta of first object
eta2: eta of second object
phi2: eta of second object
*/
float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deltaR = sqrt(pow(deltaeta(eta1, eta2), 2) + pow(deltaphi(phi1, phi2), 2));
  return deltaR;
}

/*
Find the leading and subleading jets, return jet leading and subleading pt, eta and phi
--> Arguments
pt: jet pT
eta: jet Eta
phi: jet Phi
leadpt: leading jet pT
leadeta: leading jet Eta
leadphi: leading jet Phi
sublpt: subleading jet pT
subleta: subleading jet Eta
sublphi: subleading jet Phi
*/
void find_leading_subleading(float pt,
                             float eta,
                             float phi,
                             float &leadpt,
                             float &leadeta,
                             float &leadphi,
                             float &sublpt,
                             float &subleta,
                             float &sublphi)
{
  if (pt > leadpt)
  {
    sublpt = leadpt;
    leadpt = pt;
    leadeta = eta;
    leadphi = phi;
  }
  else if (sublpt < pt)
  {
    sublpt = pt;
    subleta = eta;
    sublphi = phi;
  }
}

/*
Find bin dynamically
--> Arguments
quant_vec: vector with binning
quant: variable
*/
int find_my_bin(std::vector<double> quant_vec, double quant)
{
  int bin = -999;
  for (int ii = 0; ii < quant_vec.size() - 1; ii++)
  {
    if (quant >= quant_vec[ii] && quant < quant_vec[ii + 1])
    {
      bin = ii;
    }
  }
  return bin;
}

/*
Measure the correlation between objects
--> Arguments
jets: vector with jet informations
jet_w: vector with jet weight informations
tracks: vector with track informations
trk_w: vector with track weight informations
histo_corr: multidimentional histogram for correlations {Delta Phi, Delta Eta, track pT bin, multiplicity or centrality}
histjet: multidimentional histogram for jets in correlations {pT, Eta, Phi}
histtrk: multidimentional histogram for tracks in correlations {pT, Eta, Phi}
event_weight: event weight vector for each event
mult: multiplicity or centrality vector for each event
do_rotation: true means apply/fill rotation method, otherwise use false
N_rot: number of rotation (only use if do_rotation is true)
histo_rot: histogram using rotation method (only use if do_rotation is true)
sube_trk: vector with sube track (MC embedded samples) sube == 0 means PYTHIA embedded tracks while sube > 0 means the other MC (HYDJET, EPOS, ...)
histo_corr_subeg0: if sube > 0 save in this histogram
*/
void correlation(std::vector<TVector3> jets,
                 std::vector<double> jets_w,
                 std::vector<TVector3> tracks,
                 std::vector<double> tracks_w,
                 THnSparse *histo_corr,
                 THnSparse *histjet,
                 THnSparse *histtrk,
                 float event_weight,
                 int mult,
                 bool do_rotation,
                 int N_rot,
                 THnSparse *histo_rot,
                 std::vector<int> sube_trk,
                 THnSparse *histo_corr_subeg0)
{
  // get correlation histograms
  for (int a = 0; a < jets.size(); a++)
  { // start loop over jets
    double jet_weight = jets_w[a];
    for (int b = 0; b < tracks.size(); b++)
    { // start loop over tracks
      double trkpt = tracks[b].Pt();
      double trketa = tracks[b].Eta();
      int subetrk = sube_trk[b];
      // track efficiency correction for reco
      double trk_weight = tracks_w[b];
      // Find track and multiplicity bins
      int trkbin = (int)find_my_bin(trk_pt_bins, trkpt);
      int multcentbin = (int)find_my_bin(multiplicity_centrality_bins, (float)mult);
      // Fill jet and track quantities
      double x4D_jet[4] = {jets[a].Pt(), jets[a].Eta(), jets[a].Phi(), (double)multcentbin};
      histjet->Fill(x4D_jet, jet_weight * event_weight);
      double x4D_trk[4] = {tracks[b].Pt(), tracks[b].Eta(), tracks[b].Phi(), (double)multcentbin};
      histtrk->Fill(x4D_trk, trk_weight * event_weight);
      // Fill correlation histograms
      double del_phi = deltaphi2PC(jets[a].Phi(), tracks[b].Phi());
      double del_eta = deltaeta(jets[a].Eta(), tracks[b].Eta());
      double x4D[4] = {del_phi, del_eta, (double)trkbin, (double)multcentbin};
      if (subetrk == 0)
      {
        histo_corr->Fill(x4D, jet_weight * trk_weight * event_weight * trkpt);
      }
      else
      {
        histo_corr_subeg0->Fill(x4D, jet_weight * trk_weight * event_weight * trkpt);
      }
    }
    // get rotation histograms
    if (do_rotation)
    {
      for (int c = 0; c < N_rot; c++)
      {
        TRandom2 *r = new TRandom2();                // make a random number
        float alpha = r->Uniform(2.0 * TMath::Pi()); // make a random number between 0 and 2pi
        float newphi = jets[a].Phi() + alpha;        // add to the jet phi and reorder the axis
        if (newphi > 3.0 * TMath::Pi() / 2.)
          newphi -= 2.0 * TMath::Pi();
        if (newphi < -TMath::Pi() / 2.)
          newphi += 2.0 * TMath::Pi();
        float neweta = -1.0 * jets[a].Eta(); // invert the jet eta
        for (int d = 0; d < tracks.size();
             d++)
        { // start loop over tracks using the new jet eta and phi (similar as above)
          double trkpt = tracks[d].Pt();
          double trketa = tracks[d].Eta();
          // track efficiency correction for reco
          double trk_weight = tracks_w[d];
          // Find track and multiplicity bins
          int trkbin = (int)find_my_bin(trk_pt_bins, trkpt);
          int multcentbin = (int)find_my_bin(multiplicity_centrality_bins, (float)mult);
          // Fill correlation histograms
          double del_phi_rot = deltaphi2PC(newphi, tracks[d].Phi());
          double del_eta_rot = deltaeta(neweta, tracks[d].Eta());
          double x4D_rot[4] = {del_phi_rot, del_eta_rot, (double)trkbin, (double)multcentbin};
          histo_rot->Fill(x4D_rot, jet_weight * trk_weight * event_weight * trkpt);
        }
      }
    } // end rotation
  }   // end jet track loop
}

/*
Function to fill vectors to be used during the mixing (the ones with &)
--> Arguments
similar_events: true for using tracks only if event has one jet within the jet cut or falt for all tracks
nev: histogram with number of events stored to be used in the mixing
jets: vector with jet informations
jet_weight: vector with jet weight informations
tracks: vector with track informations
trk_weight: vector with track weight informations
mult: multiplicity or centrality
vertexz: Z vertex in centimeters
weight: event weight
ev_jet_vector: vector to be used in the mixing with jet information for each event
ev_jet_weight_vector: vector to be used in the mixing with jet weight information for each event
ev_track_vector: vector to be used in the mixing with track information for each event
ev_trk_weight_vector: vector to be used in the mixing with track weight information for each event
multvec: vector to be used in the mixing with event multiplicity or centrality information
vzvec: vector to be used in the mixing with event Z vertex position information
weightvec: vector to be used in the mixing with event weight information
*/
void fillvectors(bool similar_events,
                 TH1 *nev,
                 std::vector<TVector3> jets,
                 std::vector<double> jet_weight,
                 std::vector<TVector3> tracks,
                 std::vector<double> trk_weight,
                 int mult,
                 double vertexz,
                 double weight,
                 std::vector<std::vector<TVector3>> &ev_jet_vector,
                 std::vector<std::vector<double>> &ev_jet_weight_vector,
                 std::vector<std::vector<TVector3>> &ev_track_vector,
                 std::vector<std::vector<double>> &ev_trk_weight_vector,
                 std::vector<int> &multvec,
                 std::vector<double> &vzvec,
                 std::vector<double> &weightvec)
{
  if (similar_events)
  {
    if (jets.size() > 0 && tracks.size() > 0)
    {
      nev->Fill(0);
      ev_jet_vector.push_back(jets);
      ev_track_vector.push_back(tracks);
      ev_jet_weight_vector.push_back(jet_weight);
      ev_trk_weight_vector.push_back(trk_weight);
      multvec.push_back(mult);
      vzvec.push_back(vertexz);
      weightvec.push_back(weight);
    }
  }
  else
  {
    nev->Fill(0);
    ev_jet_vector.push_back(jets);
    ev_track_vector.push_back(tracks);
    ev_jet_weight_vector.push_back(jet_weight);
    ev_trk_weight_vector.push_back(trk_weight);
    multvec.push_back(mult);
    vzvec.push_back(vertexz);
    weightvec.push_back(weight);
  }
}

/*
  JetEnergyScale
  This funtion matches given recojet with appropriate gen jet based on the least distance.
  Multidimensional histogram of matched recopt, gen pt and multiplicity/centrality
  Arguments:
  hist_matched_vs_gen: THnSpase of matched pt on X-axis and gen pt on Y-axis and multiplicity on Z-axis
  hist_matched_vs_gen_weight : THnSpase of matched pt on X-axis and gen pt on Y-axis and multiplicity on Z-axis with weights
  hist_jes : Mathchedpt/genpt on X-axis and Genpt on Y-axis and multiplicity obn Z-axis
  hist_jes_weighted : Mathchedpt/genpt on X-axis and Genpt on Y-axis and multiplicity obn Z-axis with weights
  dR : 1D histogram for distance between jets
  Array jetpt = rawpt[nref]
  Array jeteta = jteta[nref]
  Array jetphi = jtphi[nref]
  nref = size of reco level jet variable array
  Array genpt = gen_jtp[ngen]
  Array geneta = gen_jteta[ngen]
  Array genphi = gen_jtphi[ngen]
  ngen = size of gen level jet variable array
  weight : event level weight
  mult = multiplicity of the event
  jet_pt_min = minimum cut of jetpt
  jet_pt_max = maximum cut of jetpt
  jet_eta_min = minimum cut of jeteta
  jet_eta_max = maximum cut of jeteta
  isMC = for monte carlo only
  system: colliding system
  year: year of data-taking
  energy: colliding energy
*/
void JetEnergyScaleWithMatching(THnSparseD *hist_matched_vs_gen,
                                THnSparseD *hist_matched_vs_gen_weighted,
                                THnSparseD *hist_jes,
                                THnSparseD *hist_jes_weighted,
                                THnSparseD *hist_jer,
                                THnSparseD *hist_jer_weighted,
                                TH1D *dR,
                                std::vector<TVector3> reco_jets,
                                std::vector<TVector3> gen_jets,
                                float min_dR,
                                float jet_pt_min,
                                float jet_pt_max,
                                float jet_eta_min,
                                float jet_eta_max,
                                float weight,
                                int mult,
                                bool isMC,
                                string system,
                                int year,
                                int energy)
{
  int multcentbin = (int)find_my_bin(
      multiplicity_centrality_bins,
      (float)
          mult); // Finding multiplicity bin for histogram													   // Minimum dR that will be filled. Initialized with values that will be larger than the largest possibel; distance betweeen two jets.

  for (int i = 0; i < reco_jets.size(); i++) // Recolevel jet loop
  {
    float new_dR = min_dR;
    float ratio_jes = 0.0;
    float ratio_jer =
        0.0;                                                                  // Ratio of matchedpt/genpt																  // Ratio of matchedpt/genpt
    float matched_pt = 0.0;                                                   // Matched genpt that will be filled in histogram
    if (reco_jets[i].Eta() < jet_eta_min || reco_jets[i].Eta() > jet_eta_max) // Jeteta cut at reco level
      continue;
    if (reco_jets[i].Pt() > jet_pt_min && reco_jets[i].Pt() < jet_pt_max) // Jetpt cut at recolevel
    {
      for (int j = 0; j < gen_jets.size(); j++) // Gen level jet loop
      {
        float delta_R = deltaR(reco_jets[i].Eta(),
                               reco_jets[i].Phi(),
                               gen_jets[j].Eta(),
                               gen_jets[j].Phi()); // Calculation of distance between gen and reco jet
        if (delta_R < new_dR)                      // Smallest distance between reco and gen jet
        {
          ratio_jes = reco_jets[i].Pt() / gen_jets[j].Pt();
          ratio_jer = (reco_jets[i].Pt() - gen_jets[j].Pt()) / gen_jets[j].Pt();
          matched_pt = gen_jets[j].Pt();
          new_dR = delta_R;
        }
      }
      if (ratio_jes == 0.0)
        continue;
      dR->Fill(new_dR);
      double xjetw = get_jetpT_weight(isMC, system, year, energy, reco_jets[i].Pt());
      double yjetw = get_jetpT_weight(isMC, system, year, energy, matched_pt);
      double x3D_xjetsvsyjets[3] = {reco_jets[i].Pt(), matched_pt, (double)multcentbin};
      hist_matched_vs_gen_weighted->Fill(x3D_xjetsvsyjets, xjetw * yjetw * weight);
      hist_matched_vs_gen->Fill(x3D_xjetsvsyjets);
      double x3D_jes[3] = {ratio_jes, matched_pt, (double)multcentbin};
      hist_jes_weighted->Fill(x3D_jes, xjetw * yjetw * weight);
      hist_jes->Fill(x3D_jes);
      double x3D_jer[3] = {ratio_jer, matched_pt, (double)multcentbin};
      hist_jer_weighted->Fill(x3D_jer, xjetw * yjetw * weight);
      hist_jer->Fill(x3D_jer);
      // Once matched gen jet is found for Reco jet i.e. closest gen jet to reco jet then we fill all the histograms with appropriate quantities
    }
  }
}

void JetEnergyScale(THnSparseD *hist_matched_vs_gen,
                    THnSparseD *hist_matched_vs_gen_weighted,
                    THnSparseD *hist_jes,
                    THnSparseD *hist_jes_weighted,
                    THnSparseD *hist_jer,
                    THnSparseD *hist_jer_weighted,
                    TH1D *dR,
                    std::vector<TVector3> reco_jets,
                    std::vector<TVector3> ref_jets,
                    float jet_pt_min,
                    float jet_pt_max,
                    float jet_eta_min,
                    float jet_eta_max,
                    float weight,
                    int mult,
                    bool isMC,
                    string system,
                    int year,
                    int energy)
{
  // int multcentbin = (int)find_my_bin(multiplicity_centrality_bins, (float)mult); // Finding multiplicity bin for histogram
  int multcentbin =
      (int)find_my_bin(multiplicity_centrality_bins, (float)mult); // Finding multiplicity bin for histogram

  // Minimum dR that will be filled. Initialized with values that will be larger than the largest possibel; distance betweeen two jets.
  for (int i = 0; i < reco_jets.size(); i++) // Recolevel jet loop
  {
    if (ref_jets[i].Pt() < 0.0)
      continue;
    if (reco_jets[i].Eta() > jet_eta_max)
      continue;
    if (reco_jets[i].Eta() < jet_eta_min)
      continue;
    if (reco_jets[i].Pt() > jet_pt_min_cut && reco_jets[i].Pt() < jet_pt_max_cut)
    {
      float delta_R = deltaR(
          reco_jets[i].Eta(),
          reco_jets[i].Phi(),
          ref_jets[i].Eta(),
          ref_jets[i]
              .Phi()); // Calculation of distance between gen and reco jet																		  // Smallest distance between reco and gen jet
      dR->Fill(delta_R);
      double xjetw = get_jetpT_weight(isMC, system, year, energy, reco_jets[i].Pt());
      double yjetw = get_jetpT_weight(isMC, system, year, energy, ref_jets[i].Pt());
      double x3D_xjetsvsyjets[3] = {reco_jets[i].Pt(), ref_jets[i].Pt(), (double)multcentbin};
      hist_matched_vs_gen_weighted->Fill(x3D_xjetsvsyjets, xjetw * yjetw * weight);
      hist_matched_vs_gen->Fill(x3D_xjetsvsyjets);
      double x3D_jes[4] = {reco_jets[i].Pt() / ref_jets[i].Pt(), ref_jets[i].Pt(), (double)multcentbin, ref_jets[i].Eta()};
      hist_jes_weighted->Fill(x3D_jes, xjetw * yjetw * weight);
      hist_jes->Fill(x3D_jes);
      double x3D_jer[3] = {
          (reco_jets[i].Pt() - ref_jets[i].Pt()) / ref_jets[i].Pt(), ref_jets[i].Pt(), (double)multcentbin};
      hist_jer_weighted->Fill(x3D_jer, xjetw * yjetw * weight);
      hist_jer->Fill(x3D_jer);
    }
  }
}
