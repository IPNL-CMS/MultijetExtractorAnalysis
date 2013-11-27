#include <Extractors/MultijetExtractorAnalysis/plugins/multijet_analysis.h>


using namespace std;

float computeDeltaPhi(float phi1, float phi2) {
	float deltaPhi = TMath::Abs((phi1) - (phi2));
	if(deltaPhi>TMath::Pi()){
		deltaPhi = 2*TMath::Pi()-deltaPhi;
	}
	return deltaPhi;
}

//namespace patextractor {

multijet_analysis::multijet_analysis(const edm::ParameterSet& cmsswSettings): Plugin(cmsswSettings)
{
  // Initialize the analysis parameters using the ParameterSet cmsswSettings
  
 
   // Set everything to 0
  m_met_lorentzvector = new TClonesArray("TLorentzVector");
  m_recoil_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjet_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjetgen_lorentzvector = new TClonesArray("TLorentzVector");
  m_jets_recoil_lorentzvector = new TClonesArray("TLorentzVector");
  
  //reset();

  /// Tree definition
  m_tree_Multijet = new TTree("Multijet", "Analysis info");

  /// Branches definition
  m_tree_Multijet->Branch("n_jets"         , &m_n_jets             , "n_jets/I");
  m_tree_Multijet->Branch("n_jets_puLoose"         , &m_n_jets_puLoose             , "n_jets_puLoose/I");
  m_tree_Multijet->Branch("n_goodJets"         , &m_n_goodJets             , "n_goodJets/I");
  m_tree_Multijet->Branch("jet_isPFJetLoose",  &m_jet_isPFJetLoose,   "jet_isPFJetLoose[n_jets]/I"); 
  m_tree_Multijet->Branch("jet_puJetId",  &m_jet_puJetId,   "jet_puJetId[n_jets]/I"); 
   
  m_tree_Multijet->Branch("n_muons"         , &m_n_muons             , "n_muons/I");
  m_tree_Multijet->Branch("muon_isLooseMuon",  &m_muon_isLooseMuon,   "muon_isLooseMuon[n_muons]/I");  
  m_tree_Multijet->Branch("muon_isSoftMuon",  &m_muon_isSoftMuon,   "muon_isSoftMuon[n_muons]/I");  
  m_tree_Multijet->Branch("muon_isTightMuon",  &m_muon_isTightMuon,   "muon_isTightMuon[n_muons]/I"); 
  m_tree_Multijet->Branch("muon_isHighPtMuon",  &m_muon_isHighPtMuon,    "muon_isHighPtMuon[n_muons]/I");   
  m_tree_Multijet->Branch("muon_isIsolatedMuon",  &m_muon_isIsolatedMuon,   "muon_isIsolatedMuon[n_muons]/I");  

  m_tree_Multijet->Branch("n_photons"         , &m_n_photons             , "n_photons/I");
  m_tree_Multijet->Branch("photon_isLoosePhoton",  &m_photon_isLoosePhoton,   "photon_isLoosePhoton[n_photons]/I"); 
  m_tree_Multijet->Branch("photon_isMediumPhoton",  &m_photon_isMediumPhoton,   "photon_isMediumPhoton[n_photons]/I");  
  m_tree_Multijet->Branch("photon_isTightPhoton",  &m_photon_isTightPhoton,   "photon_isTightPhoton[n_photons]/I");   

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
  m_tree_Multijet->Branch("pass_1stJet_cut", &m_pass_1stJet_cut, "pass_1stJet_cut/I");
  m_tree_Multijet->Branch("pass_2ndJet_cut", &m_pass_2ndJet_cut, "pass_2ndJet_cut/I");
  m_tree_Multijet->Branch("pass_alpha_cut", &m_pass_alpha_cut, "pass_alpha_cut/I");
  m_tree_Multijet->Branch("pass_beta_cut", &m_pass_beta_cut, "pass_beta_cut/I");
  
  m_tree_Multijet->Branch("MJB"         , &m_MJB             , "MJB/F");
  
  //m_tree_Multijet->Branch("MET", &m_met, "MET/F");
  m_tree_Multijet->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
  m_tree_Multijet->Branch("leadingjet_4vector","TClonesArray",&m_leadingjet_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("leadingjetgen_4vector","TClonesArray",&m_leadingjetgen_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("recoil_4vector","TClonesArray",&m_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("n_jets_recoil"         , &m_n_jets_recoil            
  , "n_jets_recoil/I");
  m_tree_Multijet->Branch("jets_recoil_4vector","TClonesArray",&m_jets_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("secondjetpt", &m_secondjetpt, "secondjetpt/F");
  m_tree_Multijet->Branch("alpha", &m_alpha, "alpha/F");
  m_tree_Multijet->Branch("beta", &m_beta, "beta/F");
  m_tree_Multijet->Branch("A", &m_A, "A/F");
  
  m_tree_Multijet->Branch("Muu"         , &m_Muu             , "Muu/F");
  
  
 
  //FirstJetSel(); 
  m_JET1_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("firstJet").getParameter<double>("pt_min");
  m_JET1_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("firstJet").getParameter<double>("eta_max");
  
  //JetSel()	
  m_JETS_Number_min = cmsswSettings.getParameter<edm::ParameterSet>("jet").getParameter<double>("number_min");
  m_JETS_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("jet").getParameter<double>("pt_min");
  m_JETS_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("jet").getParameter<double>("eta_max");
  
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


multijet_analysis::~multijet_analysis()
{

}

int multijet_analysis::isGoodIsolatedElectron(int index) 
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

int multijet_analysis::getNgoodJets() 
{
	int n_jet = m_jetMet->getSize();
	int n_goodjet = 0;
	
	if(n_jet == 0) n_goodjet = 0;

	if(n_jet != 0) {
		for(int i=0; i<n_jet; i++) {
			if(m_jet_isPFJetLoose[i] == 1) {
				n_goodjet = n_goodjet + 1;
			}
		}
	}
	
	return n_goodjet;
}


int multijet_analysis::isLooseMuon(int index) 
{
	int isLoose = 0;
	if(m_muon->getMuisGlobal(index) == 1 || m_muon->getMuisTracker(index))
	{
		isLoose = 1;
	}
	return isLoose;
}

int multijet_analysis::isSoftMuon(int index) 
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

int multijet_analysis::isTightMuon(int index) 
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


int multijet_analysis::isHighPtMuon(int index) 
{
	int isHighPt = 0;
	if (m_muon->getMuisHighPt(index)) {
		isHighPt = 1;
	}
	return isHighPt;
}

int multijet_analysis::isIsolatedMuon(int index) 
{
	int isIsolated = 0;
	if(m_muon->getDeltaBetaCorrectedRelativeIsolation(index) >= m_MU_Iso_max)
		return isIsolated;	
	isIsolated = 1;
	return isIsolated;
}

int multijet_analysis::isLoosePhoton(int index) 
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

int multijet_analysis::isMediumPhoton(int index) 
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

int multijet_analysis::isTightPhoton(int index) 
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

int multijet_analysis::ElectronSel()
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


int multijet_analysis::FirstJetSel()
{
	int n_jet = m_jetMet->getSize();
	
	if (! n_jet)
		return 0;
	
	int indexMaxPt = 0;
	float maxPt = fabs(m_jetMet->getP4(0)->Pt());
	
	for (int i = 0; i < n_jet; i++)
	{
		TLorentzVector *jetP = m_jetMet->getP4(i);
		if(fabs(jetP->Pt())>maxPt) 
		{
			indexMaxPt = i;
			maxPt = fabs(jetP->Pt());
		}		
	}

	if(maxPt <= m_JET1_Pt_min && fabs(m_jetMet->getP4(indexMaxPt)->Eta())>= m_JET1_Eta_max)
		return 0;
	
	return 1;
}

int multijet_analysis::JetSel1()
{
	int n_jet = m_jetMet->getSize();
	int isOK = 0;
	int n_goodjet = 0;
	int n_validgoodjet = 0;

	if(n_jet != 0) {
		n_goodjet = getNgoodJets() ;
		if(n_goodjet >= m_JETS_Number_min) {
			for(int i=0; i<m_JETS_Number_min; i++) {//at least m_JETS_Number_min good jets with pt > m_JETS_Pt_min GeV and |eta| < m_JETS_Eta_max
			TLorentzVector *jetP = m_jetMet->getP4(i);
				if(fabs(jetP->Pt()) > m_JETS_Pt_min && fabs(jetP->Eta()) < m_JETS_Eta_max) {	
					n_validgoodjet = n_validgoodjet + 1;
				}
			}
			if(n_validgoodjet >= m_JETS_Number_min) isOK = 1;
			else isOK = 0;
		}
		else isOK = 0;
	}
	else isOK=0;
	
	return isOK;
}

int multijet_analysis::JetSel2()
{
	int n_jet = m_jetMet->getSize();
	int isOK = 0;
	int n_fakejet = 0;
	
// 	if(n_jet != 0) {
// 		float n_jet_pt30 = CountJetsPt30();
// 		if(n_jet_pt30 == n_jet)	
// 			isOK = 1;
// 	}

	if(n_jet != 0) {
		for(int i=0; i<n_jet; i++) {
			if(fabs(m_jetMet->getP4(i)->Pt()) > 20.) {
				if(m_jet_isPFJetLoose[i] == 0) {
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

int multijet_analysis::CountJetsPt30()
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

int multijet_analysis::CountJetsPuLoose()
{
	int n_jet = m_jetMet->getSize();
	int n_jet_puLoose = 0;
		
	if(n_jet != 0) {	
		for (int i = 0; i < n_jet; i++)
		{
			if(m_jetMet->getPuJetId(i) == 1) 
			{
				n_jet_puLoose = n_jet_puLoose + 1;
			}		
		}
	}
	
	return n_jet_puLoose;
}


int multijet_analysis::VertexSel()
{
	int n_vtx = m_vertex->getSize();

	if (!n_vtx)
		return 0;
	else if(m_vertex->getNtracks(0)>=m_VERTEX_Tracks_min)
		return 1;
	else return 0;
}

TLorentzVector multijet_analysis::getRecoilLorentzVector()
{
	TLorentzVector recoil;
	int n_jet = m_jetMet->getSize();
	
	if (! n_jet)
		recoil = TLorentzVector(0.,0.,0.,0.);
	else {
		for (int i = 1; i < n_jet; i++)
		{
			TLorentzVector *jetP = m_jetMet->getP4(i);
			if(1) { //put here the condition on the PU jet ID
				recoil += *jetP;
			}
		}
	}
	
	return recoil;
	
}



int multijet_analysis::SecondJetSel(TLorentzVector recoil)
{
	int isOK = 0;
	int n_jet = m_jetMet->getSize();
		
	if(n_jet != 0) {
		float ptRecoil = recoil.Pt();
		float secondjetpt = m_jetMet->getP4(1)->Pt();
		float A = fabs(secondjetpt)/fabs(ptRecoil);
		if(secondjetpt<m_JET2_Pt_max && A < m_JET2_A_max)
			isOK = 1;
	}
	return isOK;
}

float multijet_analysis::computeAlpha(TLorentzVector recoil)
{
	float phiRecoil = recoil.Phi();
	float phiJet1 = m_jetMet->getP4(0)->Phi();
	float deltaPhi = computeDeltaPhi(phiRecoil, phiJet1);
		
	float alpha =  TMath::Abs(deltaPhi - TMath::Pi());
		
	return alpha;
}

int multijet_analysis::AlphaSel(TLorentzVector recoil)
{
	int isOK = 0;
	int n_jet = m_jetMet->getSize();
		
	if (n_jet != 0) {
		float alpha =  computeAlpha(recoil);
		
		if(alpha < m_ALPHA_max) {
			isOK = 1;
		}
	}
	return isOK;
}

float multijet_analysis::computeBeta(int n_jet)
{
	float beta = -1;
	if(n_jet < 2) {
		return beta;
	}
	//if(n_jet != 0) {
	if(n_jet > 1) {//we need at least 2 jets to compute a deltaPhi
		//int indexClosestJet = 1;
		float phiJet1 = m_jetMet->getP4(0)->Phi();
		float phiJet2 = m_jetMet->getP4(1)->Phi();
		beta = computeDeltaPhi(phiJet1, phiJet2);
	
		for (int i = 1; i < n_jet; i++) {
			TLorentzVector *jetP = m_jetMet->getP4(i);
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

int multijet_analysis::BetaSel()
{
	int isOK = 0;
	int n_jet = m_jetMet->getSize();

	//if(n_jet != 0) {
	if(n_jet > 1) {//we need at least 2 jets to compute a deltaPhi

		float minDeltaPhi = computeBeta(n_jet);
	

		if(minDeltaPhi>m_BETA_min) {
			isOK = 1;		
		}
	}
	return isOK;
}

float multijet_analysis::GetMJB(float ptLeading, float ptRecoil)
{
	float MJB = -1;
		
	MJB = fabs(ptLeading)/fabs(ptRecoil);
	
	return MJB;

}
	

float multijet_analysis::GetMuu()
{
	int n_mu = m_muon->getSize();
	
	int n_goodMu = 0;
	
	float Muu = -1;
	
	float E1, E2, px1, px2, py1, py2, pz1, pz2;
	
	if (n_mu == 2) {
		for (int i = 0; i < 2; i++) {
			if (! m_muon->getMuisGlobal(i))
				continue;

			TLorentzVector *muP = m_muon->getMuLorentzVector(i);

			if (fabs(muP->Pt()) <= 26)
				continue;

			if (fabs(muP->Eta()) >= 2.1)
				continue;

			if (m_muon->getMunormChi2(i) >= 10.)
				continue;

			if (m_muon->getTrackerLayersWithMeasurements(i) <= 5)
				continue;

			if (m_muon->getGlobalTrackNumberOfValidMuonHits(i) <= 0)
				continue;

			if (m_muon->getNumberOfMatchedStations(i) <= 1)
				continue;

			if (m_muon->getMudB(i) >= 0.2)
				continue;

			if (m_muon->getdZ(i) >= 0.5)
				continue;

			if (m_muon->getMunValPixelHits(i) <= 0)
				continue;
			
			if (i == 0) {
				E1  = muP->E() ;
				px1 = muP->Px() ;
				py1 = muP->Py() ;
				pz1 = muP->Pz() ; 
			}
			else {
				E2  = muP->E() ;
				px2 = muP->Px() ;
				py2 = muP->Py() ;
				pz2 = muP->Pz() ; 			
			}
			
			n_goodMu = n_goodMu + 1;
		}
		
		if(n_goodMu == 2) {
			Muu = sqrt(pow(E1+E2,2) - pow(px1+px2,2) - pow(py1+py2,2) - pow(pz1+pz2,2));		
		}		
	}

	
	return Muu;
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

void multijet_analysis::analyze(const edm::Event& event, const edm::EventSetup& iSetup, PatExtractor& extractor) {
  analyze(iSetup, extractor);
}

void multijet_analysis::analyze(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  // Do the analysis
  
	reset();
	
	//int pass_all_cuts = 0;
  
	m_muon     = std::static_pointer_cast<MuonExtractor>(extractor.getExtractor("muons"));
	m_electron = std::static_pointer_cast<ElectronExtractor>(extractor.getExtractor("electrons"));
	m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("JetMET"));
	m_photon   = std::static_pointer_cast<PhotonExtractor>(extractor.getExtractor("photons"));
	m_vertex   = std::static_pointer_cast<VertexExtractor>(extractor.getExtractor("vertex"));
	
	m_Muu = GetMuu();
	
	//m_met = m_jetMet->getMETLorentzVector(0)->Pt();
	new((*m_met_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getMETLorentzVector(0)));

	

	m_n_jets = m_jetMet->getSize();
	m_n_jets_puLoose = CountJetsPuLoose();
	
	if(m_n_jets == 0) {
		fillTree();
    		return; 
	}
	
	m_n_jets_pt30 = CountJetsPt30();
	m_n_muons = m_muon->getSize();
	m_n_goodJets = getNgoodJets() ;
	m_n_photons = m_photon->getSize();
	
	if(m_n_jets != 0) {
	  for(int i=0; i<m_n_jets; i++) {
	    m_jet_isPFJetLoose[i] =  m_jetMet->isPFJetLoose(i);
	    m_jet_puJetId[i]      =  m_jetMet->getPuJetId(i);
	  }
	}
	
	if(m_n_muons != 0) {
	  for(int i=0; i<m_n_muons; i++) {
	    m_muon_isLooseMuon[i] = isLooseMuon(i);
	    m_muon_isSoftMuon[i] = isSoftMuon(i);	
	    m_muon_isTightMuon[i] = isTightMuon(i);
	    m_muon_isHighPtMuon[i] = isHighPtMuon(i);
	    m_muon_isIsolatedMuon[i] = isIsolatedMuon(i);
	  }
	}

	if(m_n_photons != 0) {
	  for(int i=0; i<m_n_photons; i++) {
	    m_photon_isLoosePhoton[i] = isLoosePhoton(i);
	    m_photon_isMediumPhoton[i] = isMediumPhoton(i);	
	    m_photon_isTightPhoton[i] = isTightPhoton(i);
	  }
	}
	
	int res = ElectronSel();
	CHECK_RES_AND_RETURN(res, m_pass_electron_cut);

	res = VertexSel();
	CHECK_RES_AND_RETURN(res, m_pass_vertex_cut);

	res = JetSel1();
	CHECK_RES_AND_RETURN(res, m_pass_Jet_cut1);
	
	res = JetSel2();
	CHECK_RES_AND_RETURN(res, m_pass_Jet_cut2);

	res = FirstJetSel();
	CHECK_RES_AND_RETURN(res, m_pass_1stJet_cut);
	
	TLorentzVector recoil = getRecoilLorentzVector();

	res = SecondJetSel(recoil);
	CHECK_RES_AND_RETURN(res, m_pass_2ndJet_cut);

	res = AlphaSel(recoil);
	CHECK_RES_AND_RETURN(res, m_pass_alpha_cut);

	res = BetaSel();
	CHECK_RES_AND_RETURN(res, m_pass_beta_cut);


	
// 	m_pass_electron_cut = ElectronSel();
// 	m_pass_Jet_cut1 = JetSel1();
// 	m_pass_1stJet_cut = FirstJetSel();
// 	m_pass_2ndJet_cut = SecondJetSel();
// 	m_pass_photon_cut = PhotonSel();
// 	m_pass_vertex_cut = VertexSel();
// 	m_pass_alpha_cut = AlphaSel();
// 	m_pass_beta_cut = BetaSel();
	
	//pass_all_cuts = m_pass_electron_cut + m_pass_muon_cut + m_pass_1stJet_cut + m_pass_beta_cut + m_pass_Jet_cut	                + m_pass_2ndJet_cut + m_pass_photon_cut + m_pass_vertex_cut + m_pass_MET_cut + m_pass_alpha_cut;
	
	
// 	if(pass_all_cuts != 10)
// 		m_multijet_isSel = 0;
// 	else {
// 		m_multijet_isSel = 1;
// 		m_MJB = GetMJB();
// 	}

	m_multijet_isSel = 1;
	
	new((*m_leadingjet_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getP4(0)));
	new((*m_leadingjetgen_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getGenP4(0)));
	
	m_n_jets_recoil = 0;
	for (int i = 1; i < m_n_jets; i++)
	{
		TLorentzVector *jetP = m_jetMet->getP4(i);
		if(1) { //put here the condition on the PU jet ID
			new((*m_jets_recoil_lorentzvector)[m_n_jets_recoil]) TLorentzVector(*jetP);
			m_n_jets_recoil ++;
		}
	}
	new((*m_recoil_lorentzvector)[0]) TLorentzVector(recoil);
	
	m_secondjetpt = m_jetMet->getP4(1)->Pt();
	m_alpha = computeAlpha(recoil);
	m_beta = computeBeta(m_n_jets);
	
	float ptrecoil = ((TLorentzVector*) (m_recoil_lorentzvector->At(0)))->Pt();
	float leadingjetpt = ((TLorentzVector*) (m_leadingjet_lorentzvector->At(0)))->Pt();
	m_A = fabs(m_secondjetpt)/fabs(ptrecoil);
	
	m_MJB = GetMJB(leadingjetpt,ptrecoil);

	fillTree();
}



void multijet_analysis::fillTree()
{
  m_tree_Multijet->Fill();
}

void multijet_analysis::reset()
{
  m_n_jets               	 = -1;
  m_n_jets_puLoose             	 = -1;
  m_n_goodJets               	 = -1;
  m_n_muons               	 = -1;
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
  
  for (int i=0;i<m_photons_MAX;++i) 
  {
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
  m_MJB                  	 = -1;
  m_Muu                  	 = -1;
  //m_met                  	 = -1;
  m_pass_electron_cut	         = -1; 
  m_pass_Jet_cut1         	 = -1;
  m_pass_Jet_cut2         	 = -1;
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
    
  if (m_recoil_lorentzvector)
    m_recoil_lorentzvector->Clear();
    
  if (m_jets_recoil_lorentzvector)
    m_jets_recoil_lorentzvector->Clear();

  if (m_met_lorentzvector)
    m_met_lorentzvector->Clear();
}

//}

DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, multijet_analysis, "multijet_analysis");
