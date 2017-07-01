#ifndef CU_ttH_EDA_Ntuple_h
#define CU_ttH_EDA_Ntuple_h

#include <algorithm>
#include <map>

#include "TClass.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_event_vars.h"
#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

/*
 *
 * Ntuple class
 *
 */
//#ifdef __CINT__
//#pragma link C++ class CU_ttH_EDA_Ntuple+;
//#endif

class CU_ttH_EDA_Ntuple //: public TClass
{

private:
  // private member functions

  inline void fill_ntuple_electrons(const std::vector<pat::Electron> &,
                                    const MiniAODHelper &);
  inline void fill_ntuple_muons(const std::vector<pat::Muon> &, const std::vector<TLorentzVector> &,
                                const MiniAODHelper &);
  inline void fill_ntuple_jets(const std::vector<pat::Jet> &,
                               const MiniAODHelper &);
  inline void fill_ntuple_bjets(const std::vector<pat::Jet> &,
                                 const MiniAODHelper &);
  inline void fill_ntuple_gen(const CU_ttH_EDA_event_vars &);

  inline void fill_SF_SL(const CU_ttH_EDA_event_vars &);
  inline void fill_SF_DL(const CU_ttH_EDA_event_vars &);

public:
  /// function member
  CU_ttH_EDA_Ntuple();
  ~CU_ttH_EDA_Ntuple();

  void initialize();
  void set_up_branches(TTree *);
  void write_ntuple_SL(const CU_ttH_EDA_event_vars &, const MiniAODHelper &);
  void write_ntuple_DL(const CU_ttH_EDA_event_vars &, const MiniAODHelper &);

  /// variables

  // event variables
  int nEvent;
  int ls;  // luminosity section number
  int run; // run number

  int npv;
  int ttHf_cat;
  int ttHFGenFilter;
  int SL_tag;
  int DL_tag;
  int FH_tag;

  int n_ele;
  int n_mu;
  int n_lep;

  double mll;
  double ht;

  // muons
  std::vector<double> mu_pt;
  std::vector<double> mu_eta;
  std::vector<double> mu_phi;
  std::vector<double> mu_E;
  std::vector<int> mu_charge;
  std::vector<double> mu_iso;

  /*
  int    mu0_jetNDauChargedMVASel;
  double mu0_miniRelIso;
  double mu0_miniIsoCharged;
  double mu0_miniIsoNeutral;
  double mu0_jetPtRel;
  double mu0_jetPtRatio;
  double mu0_jetCSV;
  double mu0_sip3D;
  double mu0_dxy;
  double mu0_dz;
  double mu0_segmentCompatibility;
  double mu0_leptonMVA;
  double mu0_mediumID;
  double mu0_dpt_div_pt;
  int    mu0_iscutsel;
  int    mu0_ismvasel;
  int    mu0_isfakeablesel;
  int    mu1_jetNDauChargedMVASel;
  double mu1_miniRelIso;
  double mu1_miniIsoCharged;
  double mu1_miniIsoNeutral;
  double mu1_jetPtRel;
  double mu1_jetPtRatio;
  double mu1_jetCSV;
  double mu1_sip3D;
  double mu1_dxy;
  double mu1_dz;
  double mu1_segmentCompatibility;
  double mu1_leptonMVA;
  double mu1_mediumID;
  double mu1_dpt_div_pt;
  int    mu1_iscutsel;
  int    mu1_ismvasel;
  int    mu1_isfakeablesel;
  */

  // electrons
  std::vector<double> ele_pt;
  std::vector<double> ele_eta;
  std::vector<double> ele_phi;
  std::vector<double> ele_E;
  std::vector<int> ele_charge;
  std::vector<double> ele_iso;

  /*
  int    ele0_jetNDauChargedMVASel;
  double ele0_miniRelIso;
  double ele0_miniIsoCharged;
  double ele0_miniIsoNeutral;
  double ele0_jetPtRel;
  double ele0_jetPtRatio;
  double ele0_jetCSV;
  double ele0_sip3D;
  double ele0_dxy;
  double ele0_dz;
  double ele0_ntMVAeleID;
  double ele0_leptonMVA;
  int    ele0_isChargeConsistent;
  int    ele0_passesConversionVeto;
  int    ele0_nMissingHits;
  int    ele0_iscutsel;
  int    ele0_ismvasel;
  int    ele0_isfakeablesel;
  int    ele1_jetNDauChargedMVASel;
  double ele1_miniRelIso;
  double ele1_miniIsoCharged;
  double ele1_miniIsoNeutral;
  double ele1_jetPtRel;
  double ele1_jetPtRatio;
  double ele1_jetCSV;
  double ele1_sip3D;
  double ele1_dxy;
  double ele1_dz;
  double ele1_ntMVAeleID;
  double ele1_leptonMVA;
  int    ele1_isChargeConsistent;
  int    ele1_passesConversionVeto;
  int    ele1_nMissingHits;
  int    ele1_iscutsel;
  int    ele1_ismvasel;
  int    ele1_isfakeablesel;
  */

  // jets
  int n_jets;
  int n_btags;

  std::vector<double> jet_pt;
  std::vector<double> jet_eta;
  std::vector<double> jet_phi;
  std::vector<double> jet_E;
  std::vector<double> jet_CSV;

  std::vector<double> bjet_pt;
  std::vector<double> bjet_eta;
  std::vector<double> bjet_phi;
  std::vector<double> bjet_E;
  std::vector<double> bjet_CSV;

  // MET
  double PFMETpt;
  double PFMETphi;

  // Gen-Level Info

  std::vector<double> genmu_pt;
  std::vector<double> genmu_eta;
  std::vector<double> genmu_phi;
  std::vector<double> genmu_E;
  std::vector<int> genmu_charge;

  std::vector<double> genele_pt;
  std::vector<double> genele_eta;
  std::vector<double> genele_phi;
  std::vector<double> genele_E;
  std::vector<int> genele_charge;

  std::vector<double> genjet_pt;
  std::vector<double> genjet_eta;
  std::vector<double> genjet_phi;
  std::vector<double> genjet_E;


  // SF and event weight

  std::vector<double> lep_sf_id;
  std::vector<double> lep_sf_iso;
  double lep_sf_trig;

  double b_weight;
  double gen_weight;
  double PU_weight;
  double pdf_weight_up;
  double pdf_weight_down;
  double q2_weight_up;
  double q2_weight_down;

};

#endif
