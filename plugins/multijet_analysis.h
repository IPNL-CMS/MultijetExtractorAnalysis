#ifndef MULTIJET_ANALYSIS_H
#define MULTIJET_ANALYSIS_H

#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include <Extractors/PatExtractor/interface/MuonExtractor.h>
#include <Extractors/PatExtractor/interface/ElectronExtractor.h>
#include <Extractors/PatExtractor/interface/JetMETExtractor.h>
#include <Extractors/PatExtractor/interface/PhotonExtractor.h>
#include <Extractors/PatExtractor/interface/VertexExtractor.h>
#include <Extractors/PatExtractor/interface/PatExtractor.h>

class MuonExtractor;
class ElectronExtractor;
class JetMETExtractor;
class PhotonExtractor;
class VertexExtractor;
class PatExtractor;


class multijet_analysis: public patextractor::Plugin {
	public:
		multijet_analysis(const edm::ParameterSet& iConfig);
		~multijet_analysis();

		virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, PatExtractor& extractor);		
		virtual void analyze(const edm::EventSetup& iSetup, PatExtractor& extractor);
		void fillTree();
		void reset();
		float GetMuu();
		float GetMJB(float ptLeading, float ptRecoil);
		int CountJetsPt30();
		int CountJetsPuLoose();
		TLorentzVector getRecoilLorentzVector();
		int getNgoodJets();
		float computeAlpha(TLorentzVector recoil);
		float computeBeta(int n_jet);
	
	private:
		TTree*  m_tree_Multijet;
		int m_n_jets;
		int m_n_jets_puLoose;
		int m_n_goodJets;
		int m_n_muons;
		int m_n_photons;
		int m_n_electrons;
		int m_n_jets_pt30;
		int m_n_jets_recoil;
		TClonesArray* m_leadingjet_lorentzvector;
		TClonesArray* m_leadingjetgen_lorentzvector;
		TClonesArray* m_recoil_lorentzvector;
		TClonesArray* m_jets_recoil_lorentzvector;
		TClonesArray* m_jetsgen_recoil_lorentzvector;
		  
		std::shared_ptr<MuonExtractor> m_muon;
		  
		std::shared_ptr<ElectronExtractor> m_electron;
		  
		std::shared_ptr<JetMETExtractor> m_jetMet;
		
		std::shared_ptr<PhotonExtractor> m_photon;
		
		std::shared_ptr<VertexExtractor> m_vertex;
		
		float m_Muu;
		float m_MJB;
		
		//float m_met;
		TClonesArray* m_met_lorentzvector;
		float m_secondjetpt;
		float m_alpha;
		float m_beta;
		float m_A;
		
		//config
		float m_MET_Pt_min;
		
		float m_JET1_Pt_min;
		float m_JET1_Eta_max;
		
		float m_JETS_Number_min;
		float m_JETS_Pt_min;
		float m_JETS_Eta_max;
		
		//float m_VERTEX_Number_min;
		float m_VERTEX_Tracks_min;
		
		float m_JET2_Pt_max;
		float m_JET2_A_max;
		
		float m_ALPHA_max;
		float m_BETA_min;
		
		float m_ELE_Iso_max;
		
		float m_MU_Iso_max;
		
		float m_PHOT_Iso_max;
		
		static const int 	m_jets_MAX       = 200;
		int 	m_jet_isPFJetLoose[m_jets_MAX];
		int 	m_jet_puJetId[m_jets_MAX];
		
		static const int 	m_muons_MAX       = 100;
		int 	m_muon_isLooseMuon[m_muons_MAX];
		int 	m_muon_isSoftMuon[m_muons_MAX];
		int 	m_muon_isTightMuon[m_muons_MAX];
		int 	m_muon_isHighPtMuon[m_muons_MAX];
		int 	m_muon_isIsolatedMuon[m_muons_MAX];
		
		static const int 	m_photons_MAX       = 100;
		int 	m_photon_isLoosePhoton[m_photons_MAX];
		int 	m_photon_isMediumPhoton[m_photons_MAX];
		int 	m_photon_isTightPhoton[m_photons_MAX];
		
// 		static const int 	m_electrons_MAX       = 200;
// 		int 	m_electron_isGoodElectron[m_electrons_MAX];
// 		int 	m_electron_isIsolatedElectron[m_electrons_MAX];
		  
		  //Selections
		int isLooseMuon(int index);
		int isSoftMuon(int index);
		int isTightMuon(int index);
		int isHighPtMuon(int index);
		int isIsolatedMuon(int index) ;
		
		int isLoosePhoton(int index);
		int isMediumPhoton(int index);
		int isTightPhoton(int index);
		
		int isGoodIsolatedElectron(int index);
		
		int ElectronSel();
		int JetSel1();
		int JetSel2();
		int FirstJetSel();
		int VertexSel();
		int SecondJetSel(TLorentzVector recoil);
		int AlphaSel(TLorentzVector recoil);
		int BetaSel();
		  
		int m_multijet_isSel;
		  
		// Cut ; -1 event drop before arriving to this cut ; 0 cut failed, 1 cut passed
		int m_pass_electron_cut;
		int m_pass_vertex_cut;
		int m_pass_Jet_cut1;
		int m_pass_Jet_cut2;
		int m_pass_1stJet_cut;
		int m_pass_2ndJet_cut;
		int m_pass_alpha_cut;
		int m_pass_beta_cut;
};

#endif 
