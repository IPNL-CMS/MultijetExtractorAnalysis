#include <Extractors/MultijetExtractorAnalysis/plugins/multijetExtractorAnalysis.h>


using namespace std;

float computeDeltaPhi(float phi1, float phi2) {
	float deltaPhi = TMath::Abs((phi1) - (phi2));
	if(deltaPhi>TMath::Pi()){
		deltaPhi = 2*TMath::Pi()-deltaPhi;
	}
	return deltaPhi;
}

//namespace patextractor {

multijetExtractorAnalysis::multijetExtractorAnalysis(const edm::ParameterSet& cmsswSettings): Plugin(cmsswSettings)
{
  // Initialize the analysis parameters using the ParameterSet cmsswSettings
  
 
   // Set everything to 0
  m_met_lorentzvector = new TClonesArray("TLorentzVector");
  m_recoil_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjet_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjetgen_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjetraw_lorentzvector = new TClonesArray("TLorentzVector");
  m_jets_recoil_lorentzvector = new TClonesArray("TLorentzVector");
  m_jetsgen_recoil_lorentzvector = new TClonesArray("TLorentzVector");
  
  //reset();

  /// Tree definition
  m_tree_Multijet = new TTree("Multijet", "Analysis info");

  /// Branches definition
  
  m_tree_Multijet->Branch("n_totJets"         , &m_n_totJets             , "n_totJets/I");
  m_tree_Multijet->Branch("n_goodJets"         , &m_n_goodJets             , "n_goodJets/I");
  m_tree_Multijet->Branch("n_puLooseJets"         , &m_n_puLooseJets             , "n_puLooseJets/I");
  m_tree_Multijet->Branch("n_PFLooseJets"         , &m_n_PFLooseJets             , "n_PFLooseJets/I");
  m_tree_Multijet->Branch("jet_isPFJetLoose",  &m_jet_isPFJetLoose,   "jet_isPFJetLoose[n_totJets]/I"); 
  m_tree_Multijet->Branch("jet_puJetId",  &m_jet_puJetId,   "jet_puJetId[n_totJets]/I"); 
  m_tree_Multijet->Branch("goodJetsIndex",  &m_goodJetsIndex); 
  m_tree_Multijet->Branch("HT"         , &m_HT             , "HT/F");
   
  m_tree_Multijet->Branch("n_muons"         , &m_n_muons             , "n_muons/I");
  m_tree_Multijet->Branch("n_muons_loose"         , &m_n_muons_loose             , "n_muons_loose/I");
  m_tree_Multijet->Branch("n_muons_soft"         , &m_n_muons_soft             , "n_muons_soft/I");
  m_tree_Multijet->Branch("n_muons_tight"         , &m_n_muons_tight             , "n_muons_tight/I");
  m_tree_Multijet->Branch("n_muons_highPt"         , &m_n_muons_highPt             , "n_muons_highPt/I");
  m_tree_Multijet->Branch("muon_isLooseMuon",  &m_muon_isLooseMuon,   "muon_isLooseMuon[n_muons]/I");  
  m_tree_Multijet->Branch("muon_isSoftMuon",  &m_muon_isSoftMuon,   "muon_isSoftMuon[n_muons]/I");  
  m_tree_Multijet->Branch("muon_isTightMuon",  &m_muon_isTightMuon,   "muon_isTightMuon[n_muons]/I"); 
  m_tree_Multijet->Branch("muon_isHighPtMuon",  &m_muon_isHighPtMuon,    "muon_isHighPtMuon[n_muons]/I");   
  m_tree_Multijet->Branch("muon_isIsolatedMuon",  &m_muon_isIsolatedMuon,   "muon_isIsolatedMuon[n_muons]/I");  

  m_tree_Multijet->Branch("n_photons"         , &m_n_photons             , "n_photons/I");
  m_tree_Multijet->Branch("n_photons_loose"         , &m_n_photons_loose            , "n_photons_loose/I");
  m_tree_Multijet->Branch("n_photons_medium"         , &m_n_photons_medium             , "n_photons_medium/I");
  m_tree_Multijet->Branch("n_photons_tight"         , &m_n_photons_tight             , "n_photons_tight/I");  
  m_tree_Multijet->Branch("photon_isLoosePhoton",  &m_photon_isLoosePhoton,   "photon_isLoosePhoton[n_photons]/I"); 
  m_tree_Multijet->Branch("photon_isMediumPhoton",  &m_photon_isMediumPhoton,   "photon_isMediumPhoton[n_photons]/I");  
  m_tree_Multijet->Branch("photon_isTightPhoton",  &m_photon_isTightPhoton,   "photon_isTightPhoton[n_photons]/I");   
  m_tree_Multijet->Branch("photon_pt",  &m_photon_pt,   "photon_pt[n_photons]/I"); 

//   m_tree_Multijet->Branch("n_electrons"         , &m_n_electrons             , "n_electrons/I");
//   m_tree_Multijet->Branch("electron_isGoodElectron",  &m_electron_isGoodElectron,   "electron_isGoodElectron[n_electrons]/I");  
//   m_tree_Multijet->Branch("electron_isIsolatedElectron",  &m_electron_isIsolatedElectron,  
//   "electron_isIsolatedElectron[n_electrons]/I");  
  
  
  m_tree_Multijet->Branch("n_jets_pt30"         , &m_n_jets_pt30             , "n_jets_pt30/I");
  m_tree_Multijet->Branch("isSel"              , &m_multijet_isSel                 , "isSel/I");
  
// cuts
  m_tree_Multijet->Branch("pass_electron_cut", &m_pass_electron_cut, "pass_electron_cut/I");
  m_tree_Multijet->Branch("pass_vertex_cut", &m_pass_vertex_cut, "pass_vertex_cut/I");
  m_tree_Multijet->Branch("pass_Jet_cut1", &m_pass_Jet_cut1, "pass_Jet_cut1/I");
  m_tree_Multijet->Branch("pass_Jet_cut2", &m_pass_Jet_cut2, "pass_Jet_cut2/I");
  m_tree_Multijet->Branch("pass_recoil_cut", &m_pass_recoil_cut, "pass_recoil_cut/I");
  m_tree_Multijet->Branch("pass_1stJet_cut", &m_pass_1stJet_cut, "pass_1stJet_cut/I");
  m_tree_Multijet->Branch("pass_2ndJet_cut", &m_pass_2ndJet_cut, "pass_2ndJet_cut/I");
  m_tree_Multijet->Branch("pass_alpha_cut", &m_pass_alpha_cut, "pass_alpha_cut/I");
  m_tree_Multijet->Branch("pass_beta_cut", &m_pass_beta_cut, "pass_beta_cut/I");
  
  m_tree_Multijet->Branch("MJB"         , &m_MJB             , "MJB/F");
  
  //m_tree_Multijet->Branch("MET", &m_met, "MET/F");
  m_tree_Multijet->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
  m_tree_Multijet->Branch("leadingjet_4vector","TClonesArray",&m_leadingjet_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("leadingjetgen_4vector","TClonesArray",&m_leadingjetgen_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("leadingjetraw_4vector","TClonesArray",&m_leadingjetraw_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("recoil_4vector","TClonesArray",&m_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("n_jets_recoil", &m_n_jets_recoil, "n_jets_recoil/I");
  m_tree_Multijet->Branch("jets_recoil_4vector","TClonesArray",&m_jets_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("jetsgen_recoil_4vector","TClonesArray",&m_jetsgen_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("secondjetpt", &m_secondjetpt, "secondjetpt/F");
  m_tree_Multijet->Branch("alpha", &m_alpha, "alpha/F");
  m_tree_Multijet->Branch("beta", &m_beta, "beta/F");
  m_tree_Multijet->Branch("A", &m_A, "A/F");

  

  //removePUJets(); 
  m_removePUJets = cmsswSettings.getParameter<edm::ParameterSet>("PUJets").getParameter<bool>("removePUJets");  
 
  //FirstJetSel(); 
  m_JET1_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("firstJet").getParameter<double>("eta_max");
  
  //JetSel()	
  m_JETS_Number_min = cmsswSettings.getParameter<edm::ParameterSet>("jet").getParameter<double>("number_min");
  m_JETS_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("jet").getParameter<double>("pt_min");
  m_JETS_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("jet").getParameter<double>("eta_max");
  
  //Recoil jets	
  m_RECOILJETS_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("recoilJets").getParameter<double>("pt_min");
  m_RECOILJETS_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("recoilJets").getParameter<double>("eta_max");
  
  //VertexSel()	
  //m_VERTEX_Number_min = cmsswSettings.getParameter<edm::ParameterSet>("vertex").getParameter<double>("number_min");
  m_VERTEX_Tracks_min = cmsswSettings.getParameter<edm::ParameterSet>("vertex").getParameter<double>("tracks_min");
  
  //SecondJetSel()	
  m_JET2_Pt_max = cmsswSettings.getParameter<edm::ParameterSet>("secondJet").getParameter<double>("pt_max");
  m_JET2_A_max = cmsswSettings.getParameter<edm::ParameterSet>("secondJet").getParameter<double>("a_max");
	
  //AlphaSel()		
  m_ALPHA_max = cmsswSettings.getParameter<edm::ParameterSet>("alpha").getParameter<double>("alpha_max");
  
  //BetaSel()
  m_BETA_min = cmsswSettings.getParameter<edm::ParameterSet>("beta").getParameter<double>("beta_min");
 
  //ElectronSel 
  m_ELE_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("electrons").getParameter<double>("isolation_max");
		
  //MuonSel()
  m_MU_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("muons").getParameter<double>("isolation_max");

  //PhotonSel()		
  m_PHOT_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("photons").getParameter<double>("isolation_max");
}


multijetExtractorAnalysis::~multijetExtractorAnalysis()
{

}

std::vector<int> multijetExtractorAnalysis::getGoodJetsIndex() {
	std::vector<int> myVector;
	int n_jet = m_jetMet->getSize();
	
	if (n_jet > 0) {
		for(int i=0; i<n_jet; i++) {
			if(m_removePUJets) {
				if(m_jet_puJetId[i] == 7) {
					myVector.push_back(i);
				}
			}
			else {
				myVector.push_back(i);
			}
		}
	}
	return myVector;
}

float multijetExtractorAnalysis::computeHT()
{
	float HT = 0;
	int n_jet = m_goodJetsIndex.size();
	
	if (! n_jet)
		return 0;
	
	for(int i=0; i<n_jet; i++) {
		HT = HT + m_jetMet->getP4(m_goodJetsIndex.at(i))->Pt();
	}
	return HT;
}

int multijetExtractorAnalysis::isGoodIsolatedElectron(int index) 
{
	int isGood = 0;
	TLorentzVector *eP = m_electron->getEleLorentzVector(index);
	if(fabs(eP->Pt()) <= 20)
		return isGood;
		
	if (fabs(eP->Eta()) >= 2.5)
		return isGood;
		
	if (m_electron->getEleeidMVATrigV0(index) <= 0 ||
	m_electron->getEleeidMVATrigV0(index) >= 1)
		return isGood;	
		
	if (m_electron->getRhoCorrectedRelativeIsolation(index) >= m_ELE_Iso_max)
		return isGood;
		
	isGood = 1;
	return isGood;
}

int multijetExtractorAnalysis::getN_PFJets() 
{
	int n_jet = m_jetMet->getSize();
	int n_PFLooseJets = 0;
	
	if(n_jet == 0) n_PFLooseJets = 0;

	if(n_jet != 0) {
		for(int i=0; i<n_jet; i++) {
			if(m_jet_isPFJetLoose[i] == 1) {
				n_PFLooseJets = n_PFLooseJets + 1;
			}
		}
	}
	
	return n_PFLooseJets;
}


int multijetExtractorAnalysis::isLooseMuon(int index) 
{
	int isLoose = 0;
	if(m_muon->getMuisGlobal(index) == 1 || m_muon->getMuisTracker(index))
	{
		isLoose = 1;
	}
	return isLoose;
}

int multijetExtractorAnalysis::isSoftMuon(int index) 
{
	int isSoft = 0;
	if(! m_muon->getMuisGood(index))
		return isSoft;
		
	if (m_muon->getTrackerLayersWithMeasurements(index) <= 5)
		return isSoft;
		
	if (m_muon->getPixelLayerWithMeasurement(index) <= 1)
		return isSoft;	
		
	if (m_muon->getMunormChi2(index) >= 1.8)
		return isSoft;
		
	if (m_muon->getdZ(index) >= 30.)
		return isSoft;
	
	if (m_muon->getMudB(index) >= 3.)
		return isSoft;
	
	isSoft = 1;
	return isSoft;
}

int multijetExtractorAnalysis::isTightMuon(int index) 
{
	int isTight = 0;
	if (! m_muon->getMuisGlobal(index))
		return isTight;

	if (m_muon->getMunormChi2(index) >= 10.)
		return isTight;

	if (m_muon->getTrackerLayersWithMeasurements(index) <= 5)
		return isTight;

	if (m_muon->getGlobalTrackNumberOfValidMuonHits(index) <= 0)
		return isTight;

	if (m_muon->getNumberOfMatchedStations(index) <= 1)
		return isTight;

	if (m_muon->getMudB(index) >= 0.2)
		return isTight;

	if (m_muon->getdZ(index) >= 0.5)
		return isTight;

	if (m_muon->getMunValPixelHits(index) <= 0)
		return isTight;	
	
	isTight = 1;
	return isTight;


}


int multijetExtractorAnalysis::isHighPtMuon(int index) 
{
	int isHighPt = 0;
	if (m_muon->getMuisHighPt(index)) {
		isHighPt = 1;
	}
	return isHighPt;
}

int multijetExtractorAnalysis::isIsolatedMuon(int index) 
{
	int isIsolated = 0;
	if(m_muon->getDeltaBetaCorrectedRelativeIsolation(index) >= m_MU_Iso_max)
		return isIsolated;	
	isIsolated = 1;
	return isIsolated;
}

int multijetExtractorAnalysis::isLoosePhoton(int index) 
{
	int isLoose = 0;
	if(m_photon->getHadTowOverEm(index) >= 0.05)
	  return isLoose;
	  
	if(m_photon->getSigmaIetaIeta(index) >= 0.0012)
	  return isLoose;
	  
	if(m_photon->hasMatchedPromptElectron(index))
	  return isLoose;
	  
	 if(m_photon->getChargedHadronsIsolation(index) >= 2.6)
	  return isLoose;
	  
	if(m_photon->getNeutralHadronsIsolation(index) >= (3.5 + 0.04 * m_photon->getP4(index)->Pt()))
	  return isLoose;
	  
	if(m_photon->getPhotonIsolation(index) >= (1.3 + 0.005 * m_photon->getP4(index)->Pt()))
	  return isLoose;
	  
	isLoose = 1;

	return isLoose;
}

int multijetExtractorAnalysis::isMediumPhoton(int index) 
{
	int isMedium = 0;
	if(m_photon->getHadTowOverEm(index) >= 0.05)
	  return isMedium;
	  
	if(m_photon->getSigmaIetaIeta(index) >= 0.0011)
	  return isMedium;
	  
	if(m_photon->hasMatchedPromptElectron(index))
	  return isMedium;
	  
	 if(m_photon->getChargedHadronsIsolation(index) >= 1.5)
	  return isMedium;
	  
	if(m_photon->getNeutralHadronsIsolation(index) >= (1.0 + 0.04 * m_photon->getP4(index)->Pt()))
	  return isMedium;
	  
	if(m_photon->getPhotonIsolation(index) >= (0.7 + 0.005 * m_photon->getP4(index)->Pt()))
	  return isMedium;
	  
	 isMedium = 1;

	return isMedium;
}

int multijetExtractorAnalysis::isTightPhoton(int index) 
{
	int isTight = 0;
	if(m_photon->getHadTowOverEm(index) >= 0.05)
	  return isTight;
	  
	if(m_photon->getSigmaIetaIeta(index) >= 0.0011)
	  return isTight;
	  
	if(m_photon->hasMatchedPromptElectron(index))
	  return isTight;
	  
	 if(m_photon->getChargedHadronsIsolation(index) >= 0.7)
	  return isTight;
	  
	if(m_photon->getNeutralHadronsIsolation(index) >= (0.4 + 0.04 * m_photon->getP4(index)->Pt()))
	  return isTight;
	  
	if(m_photon->getPhotonIsolation(index) >= (0.5 + 0.005 * m_photon->getP4(index)->Pt()))
	  return isTight;
	  
	 isTight = 1;

	return isTight;
}

int multijetExtractorAnalysis::ElectronSel()
{
	int n_ele = m_electron->getSize();
	int pass_ele_sel = 0;

	if (!n_ele)
		pass_ele_sel = 1;
		
	int IsGoodIsolated = 0;
	int n_isol_ele = 0;

	for (int i = 0; i < n_ele; i++)
	{
		IsGoodIsolated = isGoodIsolatedElectron(i) ;
		if(IsGoodIsolated == 1) 
			n_isol_ele = n_isol_ele + 1;
	}
	if(n_isol_ele == 0) pass_ele_sel = 1;
	else if(n_isol_ele >= 1) pass_ele_sel = 0;
	return pass_ele_sel;
}

int multijetExtractorAnalysis::RecoilSel(TLorentzVector recoil_p4)
{
	int n_jet = m_goodJetsIndex.size();
	
	if (! n_jet)
		return 0;

	if(recoil_p4.Pt() <= 210.)
		return 0;
	
	return 1;
}

int multijetExtractorAnalysis::FirstJetSel()
{
	if (! m_goodJetsIndex.size())
		return 0;

	if(fabs(m_jetMet->getP4(m_goodJetsIndex.at(0))->Eta()) >= m_JET1_Eta_max)
		return 0;
	
	return 1;
}

int multijetExtractorAnalysis::JetSel1()
{
	int n_jet = m_goodJetsIndex.size();
	int isOK = 0;
	int n_validgoodjet = 0;

	if(n_jet != 0) {
		for(int i=0; i<n_jet; i++) {//at least m_JETS_Number_min good jets with pt > m_JETS_Pt_min GeV and |eta| < m_JETS_Eta_max
			TLorentzVector *jetP = m_jetMet->getP4(m_goodJetsIndex.at(i));
			if(fabs(jetP->Pt()) > m_JETS_Pt_min && fabs(jetP->Eta()) < m_JETS_Eta_max && m_jet_isPFJetLoose[m_goodJetsIndex.at(i)] == 1) {	
				n_validgoodjet = n_validgoodjet + 1;
			}
		}
		if(n_validgoodjet >= m_JETS_Number_min) isOK = 1;
		else isOK = 0;
	}
	else isOK=0;
	
	return isOK;
}

int multijetExtractorAnalysis::JetSel2()
{
	int n_jet = m_goodJetsIndex.size();
	int isOK = 0;
	int n_fakejet = 0;
	
	if(n_jet != 0) {
		for(int i=0; i<n_jet; i++) {
			if(fabs(m_jetMet->getP4(m_goodJetsIndex.at(i))->Pt()) > 20.) {
				if(m_jet_isPFJetLoose[m_goodJetsIndex.at(i)] == 0) {
					n_fakejet = n_fakejet + 1;
				}			
			}
		}
		if(n_fakejet < 1) {
			isOK = 1;
		}
		else isOK=0;
	}
	else isOK=0;
	
	return isOK;
}

int multijetExtractorAnalysis::CountJetsPt30()
{
	int n_jet = m_jetMet->getSize();
	int n_jet_pt30 = 0;
		
	if(n_jet != 0) {	
		for (int i = 0; i < n_jet; i++)
		{
			TLorentzVector *jetP = m_jetMet->getP4(i);
			if(fabs(jetP->Pt())>30) 
			{
				n_jet_pt30 = n_jet_pt30 + 1;
			}		
		}
	}
	
	return n_jet_pt30;
}

int multijetExtractorAnalysis::CountJetsPuLoose()
{
	int n_jet = m_jetMet->getSize();
	int n_jet_puLoose = 0;
		
	if(n_jet != 0) {	
		for (int i = 0; i < n_jet; i++)
		{
			if(m_jetMet->getPuJetFullId(i) < 6) 
			{
				n_jet_puLoose = n_jet_puLoose + 1;
			}		
		}
	}
	
	return n_jet_puLoose;
}

int multijetExtractorAnalysis::CountPhotonsLoose(int n_photons)
{
	int n_photons_loose = 0;
		
	if(n_photons != 0) {	
		for (int i = 0; i < n_photons; i++)
		{
			if(m_photon_isLoosePhoton[i] == 1) 
			{
				n_photons_loose = n_photons_loose + 1;
			}		
		}
	}
	
	return n_photons_loose;
}

int multijetExtractorAnalysis::CountPhotonsMedium(int n_photons)
{
	int n_photons_medium = 0;
		
	if(n_photons != 0) {	
		for (int i = 0; i < n_photons; i++)
		{
			if(m_photon_isMediumPhoton[i] == 1) 
			{
				n_photons_medium = n_photons_medium + 1;
			}		
		}
	}
	
	return n_photons_medium;
}

int multijetExtractorAnalysis::CountPhotonsTight(int n_photons)
{
	int n_photons_tight = 0;
		
	if(n_photons != 0) {	
		for (int i = 0; i < n_photons; i++)
		{
			if(m_photon_isTightPhoton[i] == 1) 
			{
				n_photons_tight = n_photons_tight + 1;
			}		
		}
	}
	
	return n_photons_tight;
}

int multijetExtractorAnalysis::CountMuonsLoose(int n_muons)
{
	int n_muons_loose = 0;
		
	if(n_muons != 0) {	
		for (int i = 0; i < n_muons; i++)
		{
			if(m_muon_isLooseMuon[i] == 1) 
			{
				n_muons_loose = n_muons_loose + 1;
			}		
		}
	}
	
	return n_muons_loose;
}

int multijetExtractorAnalysis::CountMuonsSoft(int n_muons)
{
	int n_muons_soft = 0;
		
	if(n_muons != 0) {	
		for (int i = 0; i < n_muons; i++)
		{
			if(m_muon_isSoftMuon[i] == 1) 
			{
				n_muons_soft = n_muons_soft + 1;
			}		
		}
	}
	
	return n_muons_soft;
}

int multijetExtractorAnalysis::CountMuonsTight(int n_muons)
{
	int n_muons_tight = 0;
		
	if(n_muons != 0) {	
		for (int i = 0; i < n_muons; i++)
		{
			if(m_muon_isTightMuon[i] == 1) 
			{
				n_muons_tight = n_muons_tight + 1;
			}		
		}
	}
	
	return n_muons_tight;
}

int multijetExtractorAnalysis::CountMuonsHighPt(int n_muons)
{
	int n_muons_highPt = 0;
		
	if(n_muons != 0) {	
		for (int i = 0; i < n_muons; i++)
		{
			if(m_muon_isHighPtMuon[i] == 1) 
			{
				n_muons_highPt = n_muons_highPt + 1;
			}		
		}
	}
	
	return n_muons_highPt;
}


int multijetExtractorAnalysis::VertexSel()
{
	int n_vtx = m_vertex->getSize();

	if (!n_vtx)
		return 0;
	else if(m_vertex->getNtracks(0)>=m_VERTEX_Tracks_min)
		return 1;
	else return 0;
}

TLorentzVector multijetExtractorAnalysis::getRecoilLorentzVector()
{
	TLorentzVector recoil;
	int n_jet = m_goodJetsIndex.size();
	
	if (! n_jet)
		recoil = TLorentzVector(0.,0.,0.,0.);
	else {
		for (int i = 1; i < n_jet; i++)
		{
			TLorentzVector *jetP = m_jetMet->getP4(m_goodJetsIndex.at(i));
			if(jetP->Pt() > m_RECOILJETS_Pt_min && fabs(jetP->Eta()) < m_RECOILJETS_Eta_max) {				
				recoil += *jetP;
			}

		}
	}
	
	return recoil;
	
}



int multijetExtractorAnalysis::SecondJetSel(TLorentzVector recoil)
{
	int isOK = 0;
	int n_jet = m_goodJetsIndex.size();
		
	if(n_jet != 0) {
		float ptRecoil = recoil.Pt();
		float secondjetpt = m_jetMet->getP4(m_goodJetsIndex.at(1))->Pt();
		float A = fabs(secondjetpt)/fabs(ptRecoil);
		if(secondjetpt<m_JET2_Pt_max && A < m_JET2_A_max)
			isOK = 1;
	}
	return isOK;
}

float multijetExtractorAnalysis::computeAlpha(TLorentzVector recoil)
{
	float phiRecoil = recoil.Phi();
	float phiJet1 = m_jetMet->getP4(m_goodJetsIndex.at(0))->Phi();
	float deltaPhi = computeDeltaPhi(phiRecoil, phiJet1);
		
	float alpha =  TMath::Abs(deltaPhi - TMath::Pi());
		
	return alpha;
}

int multijetExtractorAnalysis::AlphaSel(TLorentzVector recoil)
{
	int isOK = 0;
	int n_jet = m_goodJetsIndex.size();
		
	if (n_jet != 0) {
		float alpha =  computeAlpha(recoil);
		
		if(alpha < m_ALPHA_max) {
			isOK = 1;
		}
	}
	return isOK;
}

float multijetExtractorAnalysis::computeBeta(int n_jet)
{
	float beta = -1;
	if(n_jet < 2) {
		return beta;
	}
	//if(n_jet != 0) {
	if(n_jet > 1) {//we need at least 2 jets to compute a deltaPhi
		//int indexClosestJet = 1;
		float phiJet1 = m_jetMet->getP4(m_goodJetsIndex.at(0))->Phi();
		float phiJet2 = m_jetMet->getP4(m_goodJetsIndex.at(1))->Phi();
		beta = computeDeltaPhi(phiJet1, phiJet2);
	
		for (int i = 1; i < n_jet; i++) {
			TLorentzVector *jetP = m_jetMet->getP4(m_goodJetsIndex.at(i));
			float phiJet_tmp = jetP->Phi();
			float deltaPhi_tmp = computeDeltaPhi(phiJet1, phiJet_tmp);
			if(deltaPhi_tmp<beta) {
				//indexClosestJet = i;
				beta = deltaPhi_tmp;
			}			
		}
	}
	return beta;
}

int multijetExtractorAnalysis::BetaSel()
{
	int isOK = 0;
	int n_jet = m_goodJetsIndex.size();

	//if(n_jet != 0) {
	if(n_jet > 1) {//we need at least 2 jets to compute a deltaPhi

		float minDeltaPhi = computeBeta(n_jet);
	

		if(minDeltaPhi>m_BETA_min) {
			isOK = 1;		
		}
	}
	return isOK;
}

float multijetExtractorAnalysis::GetMJB(float ptLeading, float ptRecoil)
{
	float MJB = -1;
		
	MJB = fabs(ptLeading)/fabs(ptRecoil);
	
	return MJB;

}
	


#define CHECK_RES_AND_RETURN(res, var) \
  if (res != 1) { \
    var = 0; \
    m_multijet_isSel = res; \
    fillTree(); \
    return; \
  } else { \
    var = 1; \
  }

void multijetExtractorAnalysis::analyze(const edm::Event& event, const edm::EventSetup& iSetup, PatExtractor& extractor) {
  analyze(iSetup, extractor);
}

void multijetExtractorAnalysis::analyze(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  // Do the analysis
  
	reset();
	
	//int pass_all_cuts = 0;
  
	m_muon     = std::static_pointer_cast<MuonExtractor>(extractor.getExtractor("muons"));
	m_electron = std::static_pointer_cast<ElectronExtractor>(extractor.getExtractor("electrons"));
	m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("JetMET"));
	m_photon   = std::static_pointer_cast<PhotonExtractor>(extractor.getExtractor("photons"));
	m_vertex   = std::static_pointer_cast<VertexExtractor>(extractor.getExtractor("vertex"));

	new((*m_met_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getMETLorentzVector(0)));		

	m_n_totJets = m_jetMet->getSize();
	
	if(m_n_totJets == 0) {
		fillTree();
    		return; 
	}
	
	if(m_n_totJets != 0) {
	  for(int i=0; i<m_n_totJets; i++) {
	    m_jet_isPFJetLoose[i] =  m_jetMet->isPFJetLoose(i);
	    m_jet_puJetId[i]      =  m_jetMet->getPuJetFullId(i);
	  }
	}
	
	m_goodJetsIndex = getGoodJetsIndex();
	m_n_goodJets = m_goodJetsIndex.size();	
	
	if(m_n_goodJets == 0) {
		fillTree();
    		return; 
	}
	
	
	if(m_removePUJets) {
		if(m_jet_puJetId[0] < 7) {
			fillTree();
    			return; 
		}	
	}
	
	m_HT = computeHT();


	
	m_n_puLooseJets = CountJetsPuLoose();
	m_n_jets_pt30 = CountJetsPt30();
	m_n_PFLooseJets = getN_PFJets() ;	
	
	m_n_muons = m_muon->getSize();	
	m_n_photons = m_photon->getSize();
	
	if(m_n_muons != 0) {
	  for(int i=0; i<m_n_muons; i++) {
	    m_muon_isLooseMuon[i] = isLooseMuon(i);
	    m_muon_isSoftMuon[i] = isSoftMuon(i);	
	    m_muon_isTightMuon[i] = isTightMuon(i);
	    m_muon_isHighPtMuon[i] = isHighPtMuon(i);
	    m_muon_isIsolatedMuon[i] = isIsolatedMuon(i);
	  }
	  m_n_muons_loose = CountMuonsLoose(m_n_muons);
	  m_n_muons_soft = CountMuonsSoft(m_n_muons);
	  m_n_muons_tight = CountMuonsTight(m_n_muons);
	  m_n_muons_highPt = CountMuonsHighPt(m_n_muons);
	}

	if(m_n_photons != 0) {
	  for(int i=0; i<m_n_photons; i++) {
	    m_photon_pt[i] = m_photon->getP4(i)->Pt();
	    m_photon_isLoosePhoton[i] = isLoosePhoton(i);
	    m_photon_isMediumPhoton[i] = isMediumPhoton(i);	
	    m_photon_isTightPhoton[i] = isTightPhoton(i);
	  }
	  m_n_photons_loose = CountPhotonsLoose(m_n_photons);
	  m_n_photons_medium = CountPhotonsMedium(m_n_photons);
	  m_n_photons_tight = CountPhotonsTight(m_n_photons);
	}

	TLorentzVector recoil = getRecoilLorentzVector();
	
	new((*m_leadingjet_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getP4(m_goodJetsIndex.at(0))));
	new((*m_leadingjetgen_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getGenP4(m_goodJetsIndex.at(0))));
    new((*m_leadingjetraw_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getRawP4(m_goodJetsIndex.at(0))));
	
	m_n_jets_recoil = 0;
	for (int i = 1; i < m_n_goodJets; i++)
	{
		TLorentzVector *jetP = m_jetMet->getP4(m_goodJetsIndex.at(i));
		TLorentzVector *jetgenP = m_jetMet->getGenP4(m_goodJetsIndex.at(i));
		if(jetP->Pt() > m_RECOILJETS_Pt_min && fabs(jetP->Eta()) < m_RECOILJETS_Eta_max) {				
			new((*m_jets_recoil_lorentzvector)[m_n_jets_recoil]) TLorentzVector(*jetP);
			new((*m_jetsgen_recoil_lorentzvector)[m_n_jets_recoil]) TLorentzVector(*jetgenP);
			m_n_jets_recoil ++;
		}
	}
	new((*m_recoil_lorentzvector)[0]) TLorentzVector(recoil);
	
	if (m_jets_recoil_lorentzvector->At(0)) {
		m_secondjetpt = ((TLorentzVector*)m_jets_recoil_lorentzvector->At(0))->Pt();
	}
	m_alpha = computeAlpha(recoil);
	m_beta = computeBeta(m_n_goodJets);	
	float ptrecoil = ((TLorentzVector*) (m_recoil_lorentzvector->At(0)))->Pt();
	float leadingjetpt = ((TLorentzVector*) (m_leadingjet_lorentzvector->At(0)))->Pt();
	m_A = fabs(m_secondjetpt)/fabs(ptrecoil);	

	
	int res = ElectronSel();
	CHECK_RES_AND_RETURN(res, m_pass_electron_cut);

	res = VertexSel();
	CHECK_RES_AND_RETURN(res, m_pass_vertex_cut);

	res = JetSel1();
	CHECK_RES_AND_RETURN(res, m_pass_Jet_cut1);
	
	res = JetSel2();
	CHECK_RES_AND_RETURN(res, m_pass_Jet_cut2);	
	
	res = RecoilSel(recoil);
	CHECK_RES_AND_RETURN(res, m_pass_recoil_cut);

	res = FirstJetSel();
	CHECK_RES_AND_RETURN(res, m_pass_1stJet_cut);
	
	res = SecondJetSel(recoil);
	CHECK_RES_AND_RETURN(res, m_pass_2ndJet_cut);

	res = AlphaSel(recoil);
	CHECK_RES_AND_RETURN(res, m_pass_alpha_cut);

	res = BetaSel();
	CHECK_RES_AND_RETURN(res, m_pass_beta_cut);


	m_multijet_isSel = 1;	

	
	m_MJB = GetMJB(leadingjetpt,ptrecoil);

	fillTree();
}



void multijetExtractorAnalysis::fillTree()
{
  m_tree_Multijet->Fill();
}

void multijetExtractorAnalysis::reset()
{
  m_goodJetsIndex.clear();
  m_n_totJets               	 = -1;
  m_n_goodJets               	 = -1;
  m_n_puLooseJets             	 = -1;
  m_n_PFLooseJets               	 = -1;
  m_n_muons               	 = -1;
  m_n_muons_loose = -1;
  m_n_muons_soft = -1;
  m_n_muons_tight = -1;
  m_n_muons_highPt = -1;
  m_n_electrons               	 = -1;
  m_n_jets_pt30          	 = -1;
  m_n_jets_recoil                = -1;
  for (int i=0;i<m_muons_MAX;++i) 
  {
    m_muon_isLooseMuon[i]    = -1;
    m_muon_isSoftMuon[i]     = -1;	
    m_muon_isTightMuon[i]    = -1;
    m_muon_isHighPtMuon[i]   = -1;
    m_muon_isIsolatedMuon[i] = -1;
  }
  m_n_photons		       = -1;
  m_n_photons_loose = -1;
  m_n_photons_medium = -1;
  m_n_photons_tight = -1;
  for (int i=0;i<m_photons_MAX;++i) 
  {
    m_photon_pt[i]               = -1;
    m_photon_isLoosePhoton[i]    = -1;
    m_photon_isMediumPhoton[i]   = -1;	
    m_photon_isTightPhoton[i]    = -1;
  }
    
  m_secondjetpt = -99999.;
  m_A = -1;
  m_alpha = -1;
  m_beta = -1;
		
//   for (int i=0;i<m_electrons_MAX;++i) 
//   {
//     m_electron_isGoodElectron[i]     = -1;
//     m_electron_isIsolatedElectron[i] = -1;
//   }
  for(int i=0; i<m_jets_MAX; i++) {
	m_jet_isPFJetLoose[i] =  -1;
	m_jet_puJetId[i]      =  -1;
  }
  m_HT                  	 = -1;
  m_MJB                  	 = -1;
  //m_met                  	 = -1;
  m_pass_electron_cut	         = -1; 
  m_pass_Jet_cut1         	 = -1;
  m_pass_Jet_cut2         	 = -1;
  m_pass_recoil_cut      	 = -1;
  m_pass_1stJet_cut      	 = -1;
  m_pass_2ndJet_cut      	 = -1;
  m_pass_vertex_cut      	 = -1;
  m_pass_alpha_cut       	 = -1;
  m_pass_beta_cut        	 = -1;
 
  m_multijet_isSel       	 = -1;
  
  if (m_leadingjet_lorentzvector)
    m_leadingjet_lorentzvector->Clear();
    
  if (m_leadingjetgen_lorentzvector)
    m_leadingjetgen_lorentzvector->Clear();

  if (m_leadingjetraw_lorentzvector)
    m_leadingjetraw_lorentzvector->Clear();
    
  if (m_recoil_lorentzvector)
    m_recoil_lorentzvector->Clear();
    
  if (m_jets_recoil_lorentzvector)
    m_jets_recoil_lorentzvector->Clear();
    
  if (m_jetsgen_recoil_lorentzvector)
    m_jetsgen_recoil_lorentzvector->Clear();

  if (m_met_lorentzvector)
    m_met_lorentzvector->Clear();
}

//}

DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, multijetExtractorAnalysis, "multijetExtractorAnalysis");
