#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/JetCorrector.h"              // reader for JEC
#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/JetUncertainty.h"            // reader for JEU
#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/call_libraries.h"            // call libraries from ROOT and C++
#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/histogram_definition_PbPb.h" // define histograms
#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/random_mixing.h"             // random mixing function
#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/read_tree.h"                 // read the TChains
#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/uiclogo.h"                   // print UIC jets and start/stop time
#include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/vector_definition.h"         // define the vectors for mixing
// #include "/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_8_0_28/src/myProcesses/jettrackcorrelation_analyzer/function_definition.h"
void PbPb_studies(TString input_file, TString outputfilenumber)
{
  clock_t sec_start, sec_end, sec_start_mix, sec_end_mix;
  sec_start = clock(); // start timing measurement

  TDatime *date = new TDatime();

  printwelcome(true); // welcome message

  print_start(); // start timing print

  if (!do_inclusejettrack_correlation && do_leading_subleading_jettrack_correlation)
    do_Xj_or_Ajcut =
        true; // if only leading and subleading jet+track correlation are measured we make sure that Xj and Aj are true

  if (!is_MC)
    do_pthatcut = false; // MC only
  if (!is_MC)
    do_pid = false; // MC only

  if (colliding_system != "pPb")
    do_CM_pPb = false; // Only do center-of-mass for pPb

  // print important informations in the output file
  TString data_or_mc;
  if (!is_MC)
  {
    data_or_mc = "Data";
  }
  else
  {
    data_or_mc = "MC";
  }
  TString simev;
  if (similar_events)
  {
    simev = "simevs";
  }
  else
  {
    simev = "";
  }
  TString ref_sample = "norefsample";
  if (do_mixing && !do_rotation)
  {
    ref_sample = Form("mix%ievsMult%iDVz%.1f%s", N_ev_mix, Mult_or_Cent_range, DVz_range, simev.Data());
  }
  else if (!do_mixing && do_rotation)
  {
    ref_sample = Form("rot%ievs", N_of_rot);
  }
  else if (do_mixing && do_rotation)
  {
    ref_sample =
        Form("mix%ievsMult%iDVz%.1f%s_rot%ievs", N_ev_mix, Mult_or_Cent_range, DVz_range, simev.Data(), N_of_rot);
  }
  TString jet_axis;
  if (use_WTA)
  {
    jet_axis = "WTA";
  }
  else
  {
    jet_axis = "ESC";
  }
  TString smear;
  if (do_jet_smearing)
  {
    smear = "_smearing_";
  }
  else
  {
    smear = "";
  }
  if (!do_pid)
    particles = "CH";
  TString jet_type;
  if (do_inclusejettrack_correlation && !do_leading_subleading_jettrack_correlation)
    jet_type = "Incl";
  if (!do_inclusejettrack_correlation && do_leading_subleading_jettrack_correlation)
    jet_type = "LeadSubl";
  if (do_inclusejettrack_correlation && do_leading_subleading_jettrack_correlation)
    jet_type = "InclLeadSubl";
  TString XjAj;
  if (do_Xj_or_Ajcut)
  {
    XjAj = Form("_Ajmin_%.1f_Ajmax_%.1f_Xjmin_%.1f_Xjmax_%.1f", Ajmin, Ajmax, xjmin, xjmax);
  }
  else
  {
    XjAj = "";
  }

  // In case of wrong input, printout error message and kill the job
  if (year_of_datataking != 2012 && year_of_datataking != 2016 && year_of_datataking != 2017 &&
      year_of_datataking != 2018)
  {
    cout << "Data and MC not supported: choose 2012 for pp at 8 TeV, 2016 for pPb at 8.16 TeV, 2017 for pp at 5.02 TeV "
            "or XeXe at 5.44 TeV and 2018 for PbPb at 5.02 TeV"
         << endl;
    return;
  }
  if (colliding_system != "pp" && colliding_system != "pPb" && colliding_system != "XeXe" &&
      colliding_system != "PbPb")
  {
    cout << "Data and MC not supported: choose pp for proton-proton, pPb for proton-lead, PbPb for lead-lead and XeXe "
            "for xenon-xenon"
         << endl;
    return;
  }
  if (sNN_energy_GeV != 5020 && sNN_energy_GeV != 5440 && sNN_energy_GeV != 8000 && sNN_energy_GeV != 8160 &&
      sNN_energy_GeV != 13000)
  {
    cout << "Data and MC not supported: 5020 for pp 2017 or PbPb 2018, 5440 for XeXe, 8000 for pp 2018, 8160 for pPb "
            "2016"
         << endl;
    return;
  }

  // Read JEC file
  vector<string> Files;
  Files.push_back(Form("/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/aux_files/%s_%i/%s", colliding_system.Data(), sNN_energy_GeV, JEC_file.Data()));
  JetCorrector JEC(Files);

  // Track or particle efficiency file
  TFile *fileeff = TFile::Open(Form("/afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/phd/aux_files/%s_%i/%s", colliding_system.Data(), sNN_energy_GeV, trk_eff_file.Data()));
  cout << endl;
  TH2 *reff2D = nullptr;
  TH2 *rsec2D = nullptr;
  TH2 *rfak2D = nullptr;
  TH2 *rmul2D = nullptr;
  TH3 *reff3D = nullptr;
  TH3 *rsec3D = nullptr;
  TH3 *rfak3D = nullptr;
  TH3 *rmul3D = nullptr;
  fileeff->GetObject("eff", reff2D);  // Absolute Efficiency
  fileeff->GetObject("fake", rfak2D); // Fake Reconstruction Fraction
  fileeff->GetObject("sec", rsec2D);  // Multiple Reconstruction Fraction
  fileeff->GetObject("mult", rmul2D); // Non-Primary Reconstruction Fraction
  vector<TH2 *> eff_histos = {reff2D, rfak2D, rsec2D, rmul2D};
  vector<TH3 *> eff_histos3D = {reff3D, rfak3D, rsec3D, rmul3D};

  // Print the input in the screen/log
  print_input(data_or_mc, fileeff, colliding_system);
  cout << endl;

  // Read the input file(s)
  fstream inputfile;
  inputfile.open(Form("%s", input_file.Data()), ios::in);
  if (!inputfile.is_open())
  {
    cout << "List of input files not founded!" << endl;
    return;
  }
  {
    cout << "List of input files founded! --> " << input_file.Data() << endl;
  }

  // Make a chain and a vector of file names
  std::vector<TString> file_name_vector;
  string file_chain;
  while (getline(inputfile, file_chain))
  {
    file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov/%s", file_chain.c_str()));
  }
  inputfile.close();

  // Read the trees to be added in the Chain
  TChain *hlt_tree = new TChain("hltanalysis/HltTree");
  TChain *jet_tree = new TChain(Form("%s/t", jet_collection.Data()));
  TChain *trk_tree = new TChain("PbPbTracks/trackTree");
  TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *gen_tree;
  if (is_MC)
  {
    gen_tree = new TChain("HiGenParticleAna/hi");
  }
  TChain *ski_tree = new TChain("skimanalysis/HltTree");

  // add all the trees to the chain
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end();
       listIterator++)
  {
    cout << "Adding file " << *listIterator << " to the chains" << endl;
    hlt_tree->Add(*listIterator);
    trk_tree->Add(*listIterator);
    hea_tree->Add(*listIterator);
    jet_tree->Add(*listIterator);
    ski_tree->Add(*listIterator);
    if (is_MC)
    {
      gen_tree->Add(*listIterator);
    }
  }

  // Connect all chains
  hlt_tree->AddFriend(trk_tree);
  hlt_tree->AddFriend(hea_tree);
  hlt_tree->AddFriend(jet_tree);
  hlt_tree->AddFriend(ski_tree);
  if (is_MC)
  {
    hlt_tree->AddFriend(gen_tree);
  }

  // Read the desired branchs in the trees
  read_tree(hlt_tree, is_MC, use_WTA, jet_trigger.Data(), colliding_system.Data(), sNN_energy_GeV, year_of_datataking, event_filter_str, event_filter_bool); // access the tree informations

  // Use sumw2() to make sure about histogram uncertainties in ROOT
  sw2();

  int nevents = hlt_tree->GetEntries(); // number of events
  // int nevents = 100;
  cout << "Total number of events in those files: " << nevents << endl;
  cout << endl;
  cout << "-------------------------------------------------" << endl;

  // Start loop over events
  double nev = (double)nevents;
  // double nev = 1000;

  for (int i = 0; i < nevents; i++)
  {
    hlt_tree->GetEntry(i);

    if (i != 0 && (i % 10000) == 0)
    {
      double alpha = (double)i;
      cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;
    }

    //		if(i != 0 && i % 10000 == 0 ) break;

    Nevents->Fill(0); // filled after each event cut

    // Booleans to remove events which does not pass the Aj or Xj selection
    bool pass_Aj_or_Xj_reco_cut = false;
    bool pass_Aj_or_Xj_gen_cut = false;
    if (do_Xj_or_Ajcut)
    {
      pass_Aj_or_Xj_reco_cut = false;
      pass_Aj_or_Xj_gen_cut = false;
    }

    // Apply trigger
    // if(jet_trigger_bit != 1 ) continue;
    Nevents->Fill(1);

    // Apply event filters
    for (int ii = 0; ii < event_filter_bool.size(); ii++)
      if (event_filter_bool[ii] != 1)
        continue;
    Nevents->Fill(2);

    // Vectors used for objects

    // reco jets and tracks
    std::vector<TVector3> tracks_reco;
    std::vector<int> sube_tracks_reco;
    std::vector<TVector3> jets_reco;
    std::vector<TVector3> lead_jets_reco;
    std::vector<TVector3> subl_jets_reco;
    std::vector<double> track_w_reco;
    std::vector<double> jet_w_reco;
    std::vector<double> lead_jet_w_reco;
    std::vector<double> subl_jet_w_reco;

    // gen jets and tracks
    std::vector<TVector3> tracks_gen;
    std::vector<int> sube_tracks_gen;
    std::vector<TVector3> jets_gen;
    std::vector<TVector3> lead_jets_gen;
    std::vector<TVector3> subl_jets_gen;
    std::vector<double> track_w_gen;
    std::vector<double> jet_w_gen;
    std::vector<double> lead_jet_w_gen;
    std::vector<double> subl_jet_w_gen;

    // matched jets
    std::vector<TVector3> matched_jets_matched;
    std::vector<TVector3> matched_corrected_jets_matched;

    std::vector<TVector3> reco_jets_matched;
    std::vector<TVector3> gen_jets_matched;

    // apply event selections
    // Vz
    if (vertexz <= vz_cut_min || vertexz >= vz_cut_max)
      continue;
    Nevents->Fill(3);

    // pthat (MC only)
    if (do_pthatcut)
    {
      if (pthat <= pthatmin || pthat > pthatmax)
        continue;
    }
    Nevents->Fill(4);
    // cout << "here" << endl;

    // multiplicity or centrality

    int trksize = (int)ntrk;

    int mult = get_Ntrkoff_old(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trknlayer, trknhits, pfEcal, pfHcal);

    if (mult < multiplicity_centrality_bins[0] || mult > multiplicity_centrality_bins[multiplicity_centrality_bins.size() - 1])
      continue; // centrality of multiplicity range
    int multcentbin = (int)find_my_bin(multiplicity_centrality_bins, (float)mult);
    Nevents->Fill(5);

    // event weight(s), this must be applied in all histograms

    int mult_ntrkoff = get_Ntrkoff_new(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr);

    int old_ntrkoff = mult;
    // if (mult_ntrkoff < 30 || mult_ntrkoff > 149)
    // continue;

    if (old_ntrkoff > 10 && old_ntrkoff < 61)
    {
      hist_mult_10_60->Fill(old_ntrkoff, mult_ntrkoff, weight);
    }
    else if (old_ntrkoff > 60 && old_ntrkoff < 121)
    {
      hist_mult_60_120->Fill(old_ntrkoff, mult_ntrkoff, weight);
    }
    else if (old_ntrkoff > 120 && old_ntrkoff < 186)
    {
      hist_mult_120_185->Fill(old_ntrkoff, mult_ntrkoff, weight);
    }
    else if (old_ntrkoff > 185 && old_ntrkoff < 251)
    {
      hist_mult_185_250->Fill(old_ntrkoff, mult_ntrkoff, weight);
    }
    else if (old_ntrkoff > 250 && old_ntrkoff < 401)
    {
      hist_mult_250_400->Fill(old_ntrkoff, mult_ntrkoff, weight);
    }

    if (mult_ntrkoff < 35)
    {
      hist_ntrkoff_10->Fill(mult_ntrkoff);
      hist_ntrkoff_10_weighted->Fill(mult_ntrkoff, weight);
    }
    if (mult_ntrkoff > 10 && mult_ntrkoff < 56)
    {
      hist_ntrkoff_60->Fill(mult_ntrkoff);
      hist_ntrkoff_60_weighted->Fill(mult_ntrkoff, weight);
    }
    if (mult_ntrkoff > 29 && mult_ntrkoff < 86)
    {
      hist_ntrkoff_120->Fill(mult_ntrkoff);
      hist_ntrkoff_120_weighted->Fill(mult_ntrkoff, weight);
    }

    double event_weight = weight;

    hist_hibin_mult->Fill((double)hiBin, (double)mult_ntrkoff);
    hist_hibin_mult_weighted->Fill((double)hiBin, (double)mult_ntrkoff, event_weight);

    // Fill vertex, pthat and multiplicity/centrality histograms
    vzhist->Fill(vertexz);
    vzhist_weighted->Fill(vertexz, event_weight);
    pthathist->Fill(pthat);
    pthathist_weighted->Fill(pthat, event_weight);
    multiplicity->Fill(mult);
    multiplicity_weighted->Fill(mult, event_weight);

    // Reconstruction level	(Data and MC)
    // Start loop over reco tracks (trksize is number of reco tracks)
    /*for (int j = 0; j < trksize; j++)
    {
      // Define track/particle kinematics
      float trk_pt = trkpt->at(j);
      float trk_eta = trketa->at(j);
      float trk_phi = trkphi->at(j);
      // cout << trk_pt << endl;
      // In pPb case, for the center-of-mass correction if needed
      if (colliding_system == "pPb" && do_CM_pPb)
      {
        if (is_pgoing)
        {
          trk_eta = trk_eta - 0.465;
        }
        else
        {
          trk_eta = -trk_eta - 0.465;
        }
      }

      // Apply track selection (see read_tree.h to see what each variable means)
      if (fabs(trk_eta) > trk_eta_cut)
        continue;
      if (trk_pt <= trk_pt_min_cut)
        continue;
      if (highpur->at(j) == false)
        continue;
      if (fabs(trkpterr->at(j) / trkpt->at(j)) >= trk_pt_resolution_cut)
        continue;
      if (fabs(trkdcaxy->at(j) / trkdcaxyerr->at(j)) >= trk_dca_xy_cut)
        continue;
      if (fabs(trkdcaz->at(j) / trkdcazerr->at(j)) >= trk_dca_z_cut)
        continue;
      double calomatching = ((pfEcal->at(j) + pfHcal->at(j)) / cosh(trketa->at(j))) / trk_pt;
      if (colliding_system == "PbPb" || colliding_system == "XeXe")
      {
        if ((trkchi2->at(j) / trknlayer->at(j)) >= chi2_ndf_nlayer_cut)
          continue;
        if (trknhits->at(j) < nhits)
          continue;
        if (trk_pt > 20.0 && fabs(calomatching) <= calo_matching)
          continue;
      }
      // Track efficiency correction
      double trk_weight = 1.0;
      trk_weight = trk_weight * getTrkCorrWeight(eff_histos[0], eff_histos[1], eff_histos[2], eff_histos[3], trk_pt, trk_eta);

      // Track QA histogram filling
      double x4D_reco_trk[4] = {trk_pt, trk_eta, trk_phi, (double)multcentbin};
      hist_reco_trk->Fill(x4D_reco_trk);
      hist_reco_trk_corr->Fill(x4D_reco_trk, trk_weight);
      hist_reco_trk_weighted->Fill(x4D_reco_trk, trk_weight * event_weight);

      // Track vector filling
      TVector3 GoodTracks;
      GoodTracks.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
      tracks_reco.push_back(GoodTracks);
      sube_tracks_reco.push_back(0);                                                                                                      // set == 0 because this is only valid for gen
      double trk_etamix_weight = get_trketamix_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, trk_eta, true); // weight to deal with Seagull (test)
      track_w_reco.push_back(trk_weight * trk_etamix_weight);                                                                             // save weight to apply in the mixing

    } // End loop over tracks */

    // Start loop over jets
    float leadrecojet_pt = -999, leadrecojet_eta = -999, leadrecojet_phi = -999; // leading jet quantities
    float sublrecojet_pt = -999, sublrecojet_eta = -999, sublrecojet_phi = -999; // subleading jet quantities

    int jetsize = (int)nref; // number of jets in an event
    for (int j = 0; j < jetsize; j++)
    {
      // Define jet kinematics
      // float jet_pt = jtpt[j];
      float jet_rawpt = rawpt[j];
      float jet_eta = jteta[j];
      float jet_phi = jtphi[j];
      // cout << jet_rawpt << endl;

      // In pPb case, for the center-of-mass correction if needed
      if (colliding_system == "pPb" && do_CM_pPb)
      {
        if (is_pgoing)
        {
          jet_eta = jet_eta - 0.465;
        }
        else
        {
          jet_eta = -jet_eta - 0.465;
        }
      }

      // Apply JEC
      JEC.SetJetPT(jet_rawpt);
      JEC.SetJetEta(jet_eta);
      JEC.SetJetPhi(jet_phi);
      float jet_pt_corr = JEC.GetCorrectedPT();
      // Jet eta cut

      if (jet_eta < jet_eta_min_cut || jet_eta > jet_eta_max_cut)
        continue;
      // find_leading_subleading(jet_pt_corr, jet_eta, jet_phi, leadrecojet_pt, leadrecojet_eta, leadrecojet_phi, sublrecojet_pt, sublrecojet_eta, sublrecojet_phi); // Find leading and subleading jets                                                                                             // Fill histogram without any cut

      if (is_MC)
      {
        // Matched for JES studies and parton fraction
        hist_matched_jet_weighted_nocut->Fill(refpt[j], event_weight); // Fill histogram without any cut
        // Fill jet pT for parton fractions --> no cuts applied
        if (fabs(refparton_flavor[j]) == 1)
          hist_matched_jet_pt_parton_from_d->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavor[j]) == 2)
          hist_matched_jet_pt_parton_from_u->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavor[j]) == 3)
          hist_matched_jet_pt_parton_from_s->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavor[j]) == 4)
          hist_matched_jet_pt_parton_from_c->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavor[j]) == 5)
          hist_matched_jet_pt_parton_from_b->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavor[j]) == 6)
          hist_matched_jet_pt_parton_from_t->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavor[j]) == 21)
          hist_matched_jet_pt_parton_from_g->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavorForB[j]) == 1)
          hist_matched_jet_pt_parton_B_from_d->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavorForB[j]) == 2)
          hist_matched_jet_pt_parton_B_from_u->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavorForB[j]) == 3)
          hist_matched_jet_pt_parton_B_from_s->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavorForB[j]) == 4)
          hist_matched_jet_pt_parton_B_from_c->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavorForB[j]) == 5)
          hist_matched_jet_pt_parton_B_from_t->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavorForB[j]) == 6)
          hist_matched_jet_pt_parton_B_from_b->Fill(refpt[j], event_weight);
        if (fabs(refparton_flavorForB[j]) == 21)
          hist_matched_jet_pt_parton_B_from_g->Fill(refpt[j], event_weight);
        // Fill matched and reco vectors without any cut for JES study
        if (refpt[j] < 0.0)
          continue;
        TVector3 Good_M_Jets;
        Good_M_Jets.SetPtEtaPhi(refpt[j], refeta[j], refphi[j]);
        matched_jets_matched.push_back(Good_M_Jets);
        TVector3 Good_R_Jets;
        Good_R_Jets.SetPtEtaPhi(jet_pt_corr, jet_eta, jet_phi);
        reco_jets_matched.push_back(Good_R_Jets);
      }

      if (jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut)
      {                                                                                                                        // Jet pT cut
        double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr); // Jet weight (specially for MC)
        // resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
        jet_weight = jet_weight * get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr, do_jet_smearing, 0.663); // Jet smearing (For systematics)
        // Fill reco jet QA histograms
        double x4D_reco_jet[4] = {jet_rawpt, jet_eta, jet_phi, (double)multcentbin};
        hist_reco_jet->Fill(x4D_reco_jet);
        hist_reco_jet_weighted->Fill(x4D_reco_jet, event_weight * jet_weight);
        double x4D_reco_jet_corr[4] = {jet_pt_corr, jet_eta, jet_phi, (double)multcentbin};
        hist_reco_jet_corr->Fill(x4D_reco_jet_corr);
        hist_reco_jet_corr_weighted->Fill(x4D_reco_jet_corr, event_weight * jet_weight);
        // Fill reco jet vectors
        TVector3 GoodJets;
        GoodJets.SetPtEtaPhi(jet_pt_corr, jet_eta, jet_phi);
        jets_reco.push_back(GoodJets);
        jet_w_reco.push_back(jet_weight);
      }
      // Fill matched jet vectors with jet pT cut
      if (is_MC && refpt[j] > jet_pt_min_cut && refpt[j] < jet_pt_max_cut)
      {
        double x4D_match_jet[4] = {refpt[j], refeta[j], refphi[j], (double)multcentbin};
        hist_matched_jet->Fill(x4D_match_jet);
        hist_matched_jet_weighted->Fill(x4D_match_jet, event_weight);
      }

    } // End loop over jets
    if (do_Xj_or_Ajcut)
    {
      // leading/subleading jets
      if (jetsize > 1)
      {
        Nevents->Fill(6);

        double ljet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrecojet_pt); // Jet weight (specially for MC)
        // resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
        ljet_weight = ljet_weight * get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrecojet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

        double sljet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrecojet_pt); // Jet weight (specially for MC)
        // resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
        sljet_weight = sljet_weight * get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrecojet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

        // Fill leading and subleading pT histograms without cuts
        hist_reco_leadjet_pt_nocut->Fill(leadrecojet_pt);
        hist_reco_leadjet_pt_nocut_weighted->Fill(leadrecojet_pt, event_weight * ljet_weight);
        hist_reco_subljet_pt_nocut->Fill(sublrecojet_pt);
        hist_reco_subljet_pt_nocut_weighted->Fill(sublrecojet_pt, event_weight * sljet_weight);

        // leading/subleading pT cuts
        if (leadrecojet_pt > leading_pT_min && sublrecojet_pt > subleading_pT_min)
        {
          Nevents->Fill(7);

          // Fill leading/subleading jet quenching quantities
          double delta_phi_reco = deltaphi(leadrecojet_phi, sublrecojet_phi);
          double delta_phi_reco_2pc = deltaphi2PC(leadrecojet_phi, sublrecojet_phi);
          double Aj_reco = asymmetry(leadrecojet_pt, sublrecojet_pt);
          double Xj_reco = xjvar(leadrecojet_pt, sublrecojet_pt);
          double x4D_reco[4] = {Xj_reco, Aj_reco, delta_phi_reco, (double)multcentbin};
          hist_reco_lead_reco_subl_quench->Fill(x4D_reco, event_weight * ljet_weight * sljet_weight);
          double x4D_reco2pc[4] = {Xj_reco, Aj_reco, delta_phi_reco_2pc, (double)multcentbin};
          hist_reco_lead_reco_subl_quench2pc->Fill(x4D_reco2pc, event_weight * ljet_weight * sljet_weight);

          // leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
          if (fabs(leadrecojet_phi - sublrecojet_phi) > leading_subleading_deltaphi_min)
          {
            Nevents->Fill(8);

            // Fill leading and subleading jet QA histograms
            double x4D_lead[4] = {leadrecojet_pt, leadrecojet_eta, leadrecojet_phi, (double)multcentbin};
            hist_reco_leadjet->Fill(x4D_lead);
            hist_reco_leadjet_weighted->Fill(x4D_lead, event_weight * ljet_weight);
            double x4D_sublead[4] = {sublrecojet_pt, sublrecojet_eta, sublrecojet_phi, (double)multcentbin};
            hist_reco_subljet->Fill(x4D_sublead);
            hist_reco_subljet_weighted->Fill(x4D_sublead, event_weight * sljet_weight);

            // Fill leading and subleading jet vectors
            TVector3 GoodLeadingJets_reco;
            GoodLeadingJets_reco.SetPtEtaPhi(leadrecojet_pt, leadrecojet_eta, leadrecojet_phi);
            lead_jets_reco.push_back(GoodLeadingJets_reco);
            lead_jet_w_reco.push_back(ljet_weight);
            TVector3 GoodSubLeadingJets_reco;
            GoodSubLeadingJets_reco.SetPtEtaPhi(sublrecojet_pt, sublrecojet_eta, sublrecojet_phi);
            subl_jets_reco.push_back(GoodSubLeadingJets_reco);
            subl_jet_w_reco.push_back(sljet_weight);

            if ((Xj_reco >= xjmin && Xj_reco <= xjmax) && (Aj_reco >= Ajmin && Aj_reco <= Ajmax))
              pass_Aj_or_Xj_reco_cut = true; // if we apply Xj or Aj cuts
          }
        }
      }
    }

    // Measure correlations and filling mixing vectors
    // Reco-Reco
    // Inclusive jets

    // Generator level --> MC only
    int gentrksize; // number of gen tracks/particles
    if (is_MC)
      gentrksize = (int)gen_trkpt->size();
    int gen_jetsize; // number of gen jets
    if (is_MC)
      gen_jetsize = (int)ngen;

    if (is_MC)
    {
      // Start loop over gen particles

      /* for (int j = 0; j < gentrksize; j++)
       {
         // Define track/particle kinematics
         float gtrk_pt = gen_trkpt->at(j);
         float gtrk_eta = gen_trketa->at(j);
         float gtrk_phi = gen_trkphi->at(j);

         // In pPb case, for the center-of-mass correction if needed
         if (colliding_system == "pPb" && do_CM_pPb)
         {
           if (is_pgoing)
           {
             gtrk_eta = gtrk_eta - 0.465;
           }

           {
             gtrk_eta = -gtrk_eta - 0.465;
           }
         }

         // Kinematic and charge cuts
         if (fabs(gtrk_eta) > trk_eta_cut)
           continue;
         if (gen_trkpt->at(j) <= trk_pt_min_cut)
           continue;
         if (!do_pid)
           if (gen_trkchg->at(j) == 0)
             continue;
         if (do_pid)
         {
           if (fabs(gen_trkpid->at(j)) != particlepid)
             continue;
         }

         // Track/particle QA histogram filling
         double x4D_gen_trk[4] = {gtrk_pt, gtrk_eta, gtrk_phi, (double)multcentbin};
         hist_gen_trk->Fill(x4D_gen_trk);
         hist_gen_trk_weighted->Fill(x4D_gen_trk, event_weight);
         // Track/particle vector filling
         TVector3 GoodTracks_gen;
         GoodTracks_gen.SetPtEtaPhi(gtrk_pt, gtrk_eta, gtrk_phi);
         tracks_gen.push_back(GoodTracks_gen);
         sube_tracks_gen.push_back(gen_trksube->at(j));                                                                                        // get sube from the tree
         double trk_etamix_weight = get_trketamix_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gtrk_eta, false); // weight to deal with Seagull (test)
         track_w_gen.push_back(trk_etamix_weight);                                                                                             // save weight to apply in the mixing
       }
       */

      // Start loop over gen jets
      float leadgenjet_pt = -999, leadgenjet_eta = -999, leadgenjet_phi = -999; // leading jet quantities
      float sublgenjet_pt = -999, sublgenjet_eta = -999, sublgenjet_phi = -999; // subleading jet quantities

      for (int j = 0; j < gen_jetsize; j++)
      {
        // Define jet kinematics
        float gjet_pt = gen_jtpt[j];
        float gjet_eta = gen_jteta[j];
        float gjet_phi = gen_jtphi[j];

        // cout << gjet_pt << endl;

        // In pPb case, for the center-of-mass correction if needed
        if (colliding_system == "pPb" && do_CM_pPb)
        {
          if (is_pgoing)
          {
            gjet_eta = gjet_eta - 0.465;
          }
          else
          {
            gjet_eta = -gjet_eta - 0.465;
          }
        }

        if (gjet_eta < jet_eta_min_cut || gjet_eta > jet_eta_max_cut)
          continue; // jet eta cut

        // Fill vector for JES (no pT cuts)
        TVector3 Good_G_Jets;
        Good_G_Jets.SetPtEtaPhi(gjet_pt, gjet_eta, gjet_phi);
        gen_jets_matched.push_back(Good_G_Jets);
        // find_leading_subleading(gjet_pt, gjet_eta, gjet_phi, leadgenjet_pt, leadgenjet_eta, leadgenjet_phi, sublgenjet_pt, sublgenjet_eta, sublgenjet_phi); // Find leading and subleading jets
        hist_gen_jet_weighted_nocut->Fill(gjet_pt, event_weight); // Fill jet pT without cut

        if (gjet_pt > jet_pt_min_cut && gjet_pt < jet_pt_max_cut)
        { // Jet pT cut

          double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gjet_pt); // Jet weight (specially for MC)
          // resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
          jet_weight = jet_weight * get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gjet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

          // Fill gen jet QA histograms
          double x4D_gen_jet[4] = {gjet_pt, gjet_eta, gjet_phi, (double)multcentbin};
          hist_gen_jet->Fill(x4D_gen_jet);
          hist_gen_jet_weighted->Fill(x4D_gen_jet, event_weight * jet_weight);

          // Fill gen jet vectors
          TVector3 GoodJets_gen;
          GoodJets_gen.SetPtEtaPhi(gjet_pt, gjet_eta, gjet_phi);
          jets_gen.push_back(GoodJets_gen);
          jet_w_gen.push_back(jet_weight);
        }
      }

      JetEnergyScaleWithMatching(hist_matchcorrectedpt_vs_genpt, hist_matchcorrectedpt_vs_genpt_weighted, hist_jes_ratio_matchedcorrectedpt_genpt, hist_jes_ratio_matchedcorrectedpt_genpt_weighted, hist_jer_ratio_matchedcorrectedpt_genpt, hist_jer_ratio_matchedcorrectedpt_genpt_weighted, dR_Manual, reco_jets_matched, gen_jets_matched, 0.8, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, event_weight, mult, is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV);

      // leading/subleading jets
      if (do_Xj_or_Ajcut)
      {
        if (gen_jetsize > 1)
        {
          double ljet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadgenjet_pt); // Jet weight (specially for MC)
          // resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
          ljet_weight = ljet_weight * get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadgenjet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

          double sljet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublgenjet_pt); // Jet weight (specially for MC)
          // resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
          sljet_weight = sljet_weight * get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublgenjet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

          // Fill leading and subleading pT histograms without cuts
          hist_gen_leadjet_pt_nocut->Fill(leadgenjet_pt);
          hist_gen_leadjet_pt_nocut_weighted->Fill(leadgenjet_pt, event_weight * ljet_weight);
          hist_gen_subljet_pt_nocut->Fill(sublgenjet_pt);
          hist_gen_subljet_pt_nocut_weighted->Fill(sublgenjet_pt, event_weight * sljet_weight);

          // leading/subleading pT cuts
          if (leadgenjet_pt > leading_pT_min && sublgenjet_pt > subleading_pT_min)
          {
            // Fill leading/subleading jet quenching quantities
            double delta_phi_gen = deltaphi(leadgenjet_phi, sublgenjet_phi);
            double delta_phi_gen_2pc = deltaphi2PC(leadgenjet_phi, sublgenjet_phi);
            double Aj_gen = asymmetry(leadgenjet_pt, sublgenjet_pt);
            double Xj_gen = xjvar(leadgenjet_pt, sublgenjet_pt);
            double x4D_gen[4] = {Xj_gen, Aj_gen, delta_phi_gen, (double)multcentbin};
            hist_gen_lead_gen_subl_quench->Fill(x4D_gen, event_weight * ljet_weight * sljet_weight);
            double x4D_gen2pc[4] = {Xj_gen, Aj_gen, delta_phi_gen_2pc, (double)multcentbin};
            hist_gen_lead_gen_subl_quench2pc->Fill(x4D_gen2pc, event_weight * ljet_weight * sljet_weight);

            // leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
            if (fabs(leadgenjet_phi - sublgenjet_phi) > leading_subleading_deltaphi_min)
            {
              // Fill leading and subleading jet QA histograms
              double x4D_lead[4] = {leadgenjet_pt, leadgenjet_eta, leadgenjet_phi, (double)multcentbin};
              hist_gen_leadjet->Fill(x4D_lead);
              hist_gen_leadjet_weighted->Fill(x4D_lead, event_weight * ljet_weight);
              double x4D_sublead[4] = {sublgenjet_pt, sublgenjet_eta, sublgenjet_phi, (double)multcentbin};
              hist_gen_subljet->Fill(x4D_sublead);
              hist_gen_subljet_weighted->Fill(x4D_sublead, event_weight * sljet_weight);

              // Fill leading and subleading jet vectors
              TVector3 GoodLeadingJets_gen;
              GoodLeadingJets_gen.SetPtEtaPhi(leadgenjet_pt, leadgenjet_eta, leadgenjet_phi);
              lead_jets_gen.push_back(GoodLeadingJets_gen);
              lead_jet_w_gen.push_back(ljet_weight);
              TVector3 GoodSubLeadingJets_gen;
              GoodSubLeadingJets_gen.SetPtEtaPhi(sublgenjet_pt, sublgenjet_eta, sublgenjet_phi);
              subl_jets_gen.push_back(GoodSubLeadingJets_gen);
              subl_jet_w_gen.push_back(sljet_weight);

              if ((Xj_gen >= xjmin && Xj_gen < xjmax) && (Aj_gen >= Ajmin && Aj_gen < Ajmax))
                pass_Aj_or_Xj_gen_cut = true; // if we apply Xj or Aj cuts
            }
          }
        }
      }

      // Calculate JES
      // JetEnergyScale(hist_matchcorrectedpt_vs_genpt, hist_matchcorrectedpt_vs_genpt_weighted, hist_jes_matchedcorrectedgen, hist_jes_matchedcorrectedgen_weighted, dR, reco_jets_matched, gen_jets_matched, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, event_weight, mult, is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV);
      JetEnergyScale(hist_matchcorrectedpt_vs_refpt, hist_matchcorrectedpt_vs_refpt_weighted, hist_jes_ratio_matchedcorrectedpt_refpt, hist_jes_ratio_matchedcorrectedpt_refpt_weighted, hist_jer_ratio_matchedcorrectedpt_refpt, hist_jer_ratio_matchedcorrectedpt_refpt_weighted, dR, reco_jets_matched, matched_jets_matched, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, event_weight, mult, is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV);
    }
    /*
        if (is_MC)
        {
          int refparton;
          for (int j = 0; j < nref; j++)
          {
            if (refpt[j] < 40.0)
              continue;
            if (fabs(refeta[j]) > 1.6)
              continue;

            if (fabs(refparton_flavorForB[j]) == 21)
            {
              refparton = 7;
            }
            else if (fabs(refparton_flavorForB[j]) >= 1 && fabs(refparton_flavorForB[j]) <= 6)

            {
              refparton = fabs(refparton_flavorForB[j]);
            }
            else
            {
              refparton = 0;
            }

            double x4D_jtparton[4] = {refpt[j], refeta[j], refphi[j], (double)refparton};
            hist_flavor_jets_weighted->Fill(x4D_jtparton, event_weight);
            hist_flavor_jets->Fill(x4D_jtparton);

            int jetmult = 0;
            for (int k = 0; k < gentrksize; k++)
            {
              if (gen_trkchg->at(k) == 0)
                continue;
              if (fabs(gen_trketa->at(k)) > 2.4)
                continue;
              if (gen_trkpt->at(k) < 2.0)
                continue;
              float delR = deltaR(refeta[j], refphi[j], gen_trketa->at(k), gen_trkphi->at(k));

              if (delR < 0.4)
              {
                // cout << refparton << endl;
                jetmult = jetmult + 1;
              }
            }

            double x3D_jtrk[3] = {refpt[j], (double)jetmult, (double)refparton};
            hist_jet_multiplicity_weighted->Fill(x3D_jtrk, event_weight);
            hist_jet_multiplicity->Fill(x3D_jtrk);
          }
        }*/
  } // End loop over events

  // Output file name
  cout << endl;
  cout << "Writing histograms on ";
  cout << endl;

  // Make an output file
  // string file_output = Form("%s_%s_%iGeV_%s_%s_%s_Jptmin_%.1f_Jptmax_%.1f_Jetamin_%.1f_Jetamax_%.1f_%s%s%s_%s_%s_%i_f%i", colliding_system.Data(), data_or_mc.Data(), sNN_energy_GeV, jet_type.Data(), jet_collection.Data(), jet_trigger.Data(), jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, jet_axis.Data(), smear.Data(), XjAj.Data(), ref_sample.Data(), particles.Data(), date->GetDate(), ouputfilenumber); // output file
  // std::replace(file_output.begin(), file_output.end(), '.', 'p');                                                                                                                                                                                                                                                                                                                                                                       // replace . to p
  // std::replace(file_output.begin(), file_output.end(), '-', 'N');                                                                                                                                                                                                                                                                                                                                                                       // replace - to N for negative

  // Open, write and close the output file
  TFile *MyFile = new TFile(Form("/eos/user/a/ahingraj/outputs/%s.root", outputfilenumber.Data()), "RECREATE");
  if (MyFile->IsOpen())
    cout << "output file: " << outputfilenumber.Data() << ".root" << endl;
  MyFile->cd();

  // Write in different folders (see histogram_definition.h)
  MyFile->mkdir("QA_histograms");
  MyFile->cd("QA_histograms");
  w_QA_hist(is_MC, do_leading_subleading_jettrack_correlation);

  MyFile->mkdir("correlation_reco_reco_histograms");
  MyFile->cd("correlation_reco_reco_histograms");
  w_recoreco_hist(do_mixing, do_rotation, do_inclusejettrack_correlation, do_leading_subleading_jettrack_correlation);

  if (is_MC)
  {
    MyFile->mkdir("correlation_reco_gen_histograms");
    MyFile->cd("correlation_reco_gen_histograms");
    w_recogen_hist(do_mixing, do_rotation, do_inclusejettrack_correlation, do_leading_subleading_jettrack_correlation);

    MyFile->mkdir("correlation_gen_reco_histograms");
    MyFile->cd("correlation_gen_reco_histograms");
    w_genreco_hist(do_mixing, do_rotation, do_inclusejettrack_correlation, do_leading_subleading_jettrack_correlation);

    MyFile->mkdir("correlation_gen_gen_histograms");
    MyFile->cd("correlation_gen_gen_histograms");
    w_gengen_hist(do_mixing, do_rotation, do_inclusejettrack_correlation, do_leading_subleading_jettrack_correlation);
  }

  MyFile->mkdir("jetquenching_histograms");
  MyFile->cd("jetquenching_histograms");
  w_jetquenching_hist(is_MC);

  if (is_MC)
  {
    MyFile->mkdir("QA_matchedparton_histograms");
    MyFile->cd("QA_matchedparton_histograms");
    w_QA_parton_hist();

    MyFile->mkdir("JES");
    MyFile->cd("JES");
    w_jes_hist();
  }

  MyFile->Close();

  cout << endl;
  cout << "------------------------------------- DONE --------------------------------------" << endl;
  cout << endl;

  sec_end = clock(); // stop time counting
  cout << "========================================" << endl;
  cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
  cout << "========================================" << endl;
  print_stop(); // Print time, date and hour when it stops
}
