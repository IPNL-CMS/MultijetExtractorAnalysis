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
  m_genmet_lorentzvector = new TClonesArray("TLorentzVector");
  m_met_puSubstract_lorentzvector = new TClonesArray("TLorentzVector");
  m_recoil_lorentzvector = new TClonesArray("TLorentzVector");
  m_genrecoil_lorentzvector = new TClonesArray("TLorentzVector");
  m_pu_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjet_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjetgen_lorentzvector = new TClonesArray("TLorentzVector");
  m_leadingjetraw_lorentzvector = new TClonesArray("TLorentzVector");
  m_jets_recoil_lorentzvector = new TClonesArray("TLorentzVector");
  m_jets_pu_lorentzvector = new TClonesArray("TLorentzVector");
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
  
  m_tree_Multijet->Branch("MJB"         , &m_MJB             , "MJB/F");
  
  //m_tree_Multijet->Branch("MET", &m_met, "MET/F");
  m_tree_Multijet->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
  m_tree_Multijet->Branch("genmet_4vector","TClonesArray",&m_genmet_lorentzvector, 1000, 0);
  m_tree_Multijet->Branch("met_puSubstract_4vector","TClonesArray",&m_met_puSubstract_lorentzvector, 1000, 0);
  m_tree_Multijet->Branch("leadingjet_4vector","TClonesArray",&m_leadingjet_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("leadingjetgen_4vector","TClonesArray",&m_leadingjetgen_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("leadingjetraw_4vector","TClonesArray",&m_leadingjetraw_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("recoil_4vector","TClonesArray",&m_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("genrecoil_4vector","TClonesArray",&m_genrecoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("pu_4vector","TClonesArray",&m_pu_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("n_jets_recoil", &m_n_jets_recoil, "n_jets_recoil/I");
  m_tree_Multijet->Branch("jets_recoil_4vector","TClonesArray",&m_jets_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("jets_pu_4vector","TClonesArray",&m_jets_pu_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("jetsgen_recoil_4vector","TClonesArray",&m_jetsgen_recoil_lorentzvector, 5000, 0);
  m_tree_Multijet->Branch("secondjetpt", &m_secondjetpt, "secondjetpt/F");
  m_tree_Multijet->Branch("alpha", &m_alpha, "alpha/F");
  m_tree_Multijet->Branch("beta", &m_beta, "beta/F");
  m_tree_Multijet->Branch("A", &m_A, "A/F");

  

  //removePUJets(); 
  m_removePUJets = cmsswSettings.getParameter<edm::ParameterSet>("PUJets").getParameter<bool>("removePUJets");  
  m_PUJets_Id_min = cmsswSettings.getParameter<edm::ParameterSet>("PUJets").getParameter<double>("id_min");
 
  //Recoil jets	
  m_RECOILJETS_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("recoilJets").getParameter<double>("pt_min");
  m_RECOILJETS_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("recoilJets").getParameter<double>("eta_max");
  
  //VertexSel()	
  //m_VERTEX_Number_min = cmsswSettings.getParameter<edm::ParameterSet>("vertex").getParameter<double>("number_min");
  m_VERTEX_Tracks_min = cmsswSettings.getParameter<edm::ParameterSet>("vertex").getParameter<double>("tracks_min");

  //MuonSel()
  m_MU_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("muons").getParameter<double>("isolation_max");

}


multijetExtractorAnalysis::~multijetExtractorAnalysis()
{

}

std::vector<int> multijetExtractorAnalysis::getGoodJetsIndex() {
	std::vector<int> myVector;
	int n_jet = m_jetMet->getSize();


  TLorentzVector pu;
	
	if (n_jet > 0) {
		for(int i=0; i<n_jet; i++) {
			if(m_removePUJets) {
				if(m_jet_puJetId[i] >= m_PUJets_Id_min) {
					myVector.push_back(i);
				}
        else {
          TLorentzVector *jetP = m_jetMet->getP4(i);	
          pu += *jetP;
		    	new((*m_jets_pu_lorentzvector)[i]) TLorentzVector(*jetP);       
        }
			}
			else {
				myVector.push_back(i);
			}
		}
	}
  new((*m_pu_lorentzvector)[0]) TLorentzVector(pu);
	return myVector;
}

float multijetExtractorAnalysis::computeHT()
{
	float HT = 0;
	int n_jet = m_goodJetsIndex.size();
	
	if (! n_jet)
		return 0;
	
	for(int i=0; i<n_jet; i++) {
   	TLorentzVector *jetP = m_jetMet->getP4(m_goodJetsIndex.at(i));
    if(jetP->Pt() > m_RECOILJETS_Pt_min && fabs(jetP->Eta()) < m_RECOILJETS_Eta_max) {	
      HT = HT + m_jetMet->getP4(m_goodJetsIndex.at(i))->Pt();
    }
	}
	return HT;
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
  // see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Loose_Muon
  // function isLooseMuon if isPFMuon and (isGlobalMuon or isTrackerMuon)
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

int multijetExtractorAnalysis::ElectronSel() // no isolated electron (with loose ID)
{
	int n_ele = m_electron->getSize();
	int pass_ele_sel = 0;

  if (!n_ele)
    pass_ele_sel = 1;
 
  int n_isol_ele = 0;

  for (int i = 0; i < n_ele; i++)	{
    if (m_electron->passLooseId(i)) { // in runII, a loose electron is an isolated electron: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      n_isol_ele = n_isol_ele + 1;
    }
  }

  if(n_isol_ele == 0) pass_ele_sel = 1;
  else if(n_isol_ele >= 1) pass_ele_sel = 0;
	return pass_ele_sel;
}

int multijetExtractorAnalysis::MuonSel() // no isolated muon (with loose ID)
{
	int n_muon = m_muon->getSize();
	int pass_muon_sel = 0;

	if (!n_muon)
		pass_muon_sel = 1;
	
	int n_isol_muon = 0;

	for (int i = 0; i < n_muon; i++) {
	  if (isIsolatedMuon(i) && m_muon->isLooseMuon(i)) {
			n_isol_muon = n_isol_muon + 1;
    }
	}

	if(n_isol_muon == 0) pass_muon_sel = 1;
	else if(n_isol_muon >= 1) pass_muon_sel = 0;
	return pass_muon_sel;
}

int multijetExtractorAnalysis::PhotonSel()  // no isolated photon (with loose ID)
{
	int n_pho = m_photon->getSize();
	int pass_pho_sel = 0;

	if (!n_pho)
		pass_pho_sel = 1;

  int n_isol_pho = 0;

	for (int i = 0; i < n_pho; i++)	{
	  if (m_photon->passLooseId(i)) { // in runII, a loose photon is an isolated photon: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
			n_isol_pho = n_isol_pho + 1;
	  }
  }

	if(n_isol_pho == 0) pass_pho_sel = 1;
	else if(n_isol_pho >= 1) pass_pho_sel = 0;
	return pass_pho_sel;
}



int multijetExtractorAnalysis::JetSel1()
{
	int n_jet = m_goodJetsIndex.size();
	int isOK = 0;
	int n_validgoodjet = 0;

	if(n_jet != 0) {
		for(int i=0; i<n_jet; i++) {//at least 3 good jets with pt > 25 GeV and |eta| < 2.8
			TLorentzVector *jetP = m_jetMet->getP4(m_goodJetsIndex.at(i));
			if(fabs(jetP->Pt()) > 25. && fabs(jetP->Eta()) < 2.8 && m_jet_isPFJetLoose[m_goodJetsIndex.at(i)] == 1) {	
				n_validgoodjet = n_validgoodjet + 1;
			}
		}
		if(n_validgoodjet >= 3) isOK = 1;
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
	else return 1;
 /* else if(m_vertex->getNtracks(0)>=m_VERTEX_Tracks_min)*/
		//return 1;
	/*else return 0;*/
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



float multijetExtractorAnalysis::computeAlpha(TLorentzVector recoil)
{
	float phiRecoil = recoil.Phi();
	float phiJet1 = m_jetMet->getP4(m_goodJetsIndex.at(0))->Phi();
	float deltaPhi = computeDeltaPhi(phiRecoil, phiJet1);
		
	float alpha =  TMath::Abs(deltaPhi - TMath::Pi());
		
	return alpha;
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
  
	m_muon     = std::static_pointer_cast<MuonExtractor>(extractor.getExtractor("muon_PF"));
	m_electron = std::static_pointer_cast<ElectronExtractor>(extractor.getExtractor("electron_PF"));
	m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("jetmet"));
	m_photon   = std::static_pointer_cast<PhotonExtractor>(extractor.getExtractor("photon"));
	m_vertex   = std::static_pointer_cast<VertexExtractor>(extractor.getExtractor("Vertices"));

	new((*m_met_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getMETLorentzVector(0)));		
	new((*m_genmet_lorentzvector)[0]) TLorentzVector(*(m_jetMet->getGenMETLorentzVector(0)));		

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
		if(m_jet_puJetId[0] < m_PUJets_Id_min) {
			fillTree();
    			return; 
		}	
	}
	
	m_HT = computeHT();

  TLorentzVector metCorr = *(m_jetMet->getMETLorentzVector(0));
  TLorentzVector pu = *((TLorentzVector*)m_pu_lorentzvector->At(0));
  metCorr.SetPx(metCorr.Px() + pu.Px());
  metCorr.SetPy(metCorr.Py() + pu.Py());
  metCorr.SetE(sqrt(metCorr.Py()*metCorr.Py() + metCorr.Px()*metCorr.Px()));
  new((*m_met_puSubstract_lorentzvector)[0]) TLorentzVector(metCorr);

	
	m_n_puLooseJets = CountJetsPuLoose();
	m_n_jets_pt30 = CountJetsPt30();
	m_n_PFLooseJets = getN_PFJets() ;	
	
	m_n_muons = m_muon->getSize();	
	m_n_photons = m_photon->getSize();
	
	if(m_n_muons != 0) {
	  for(int i=0; i<m_n_muons; i++) {
	    m_muon_isLooseMuon[i] = m_muon->isLooseMuon(i);
	    m_muon_isSoftMuon[i] = m_muon->isSoftMuon(i);
	    m_muon_isTightMuon[i] = m_muon->isTightMuon(i);
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
	    m_photon_isLoosePhoton[i] = m_photon->passLooseId(i); // in runII, a loose photon is an isolated photon: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
	    m_photon_isMediumPhoton[i] = m_photon->passMediumId(i);	
	    m_photon_isTightPhoton[i] = m_photon->passTightId(i);
	  }
	  m_n_photons_loose = CountPhotonsLoose(m_n_photons);
	  m_n_photons_medium = CountPhotonsMedium(m_n_photons);
	  m_n_photons_tight = CountPhotonsTight(m_n_photons);
	}

	TLorentzVector recoil = getRecoilLorentzVector();
  TLorentzVector genrecoil;
	
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
      genrecoil += *jetgenP;
		}
	}

  if(recoil.Pt())
    new((*m_recoil_lorentzvector)[0]) TLorentzVector(recoil);

  if(genrecoil.Pt())
    new((*m_genrecoil_lorentzvector)[0]) TLorentzVector(genrecoil);
	
	if (m_jets_recoil_lorentzvector->At(0)) {
		m_secondjetpt = ((TLorentzVector*)m_jets_recoil_lorentzvector->At(0))->Pt();
	}
	m_alpha = computeAlpha(recoil);
	m_beta = computeBeta(m_n_goodJets);	
  float ptrecoil = -1.;
  if(m_recoil_lorentzvector->GetEntriesFast()) {
	  ptrecoil = ((TLorentzVector*) (m_recoil_lorentzvector->At(0)))->Pt();
  }
	float leadingjetpt = ((TLorentzVector*) (m_leadingjet_lorentzvector->At(0)))->Pt();
	m_A = fabs(m_secondjetpt)/fabs(ptrecoil);	

	
	int res = ElectronSel();
	CHECK_RES_AND_RETURN(res, m_pass_electron_cut);

	res = MuonSel();
	CHECK_RES_AND_RETURN(res, m_pass_muon_cut);

	res = PhotonSel();
	CHECK_RES_AND_RETURN(res, m_pass_photon_cut);

	res = VertexSel();
	CHECK_RES_AND_RETURN(res, m_pass_vertex_cut);

	res = JetSel1();
	CHECK_RES_AND_RETURN(res, m_pass_Jet_cut1);
	
	res = JetSel2();
	CHECK_RES_AND_RETURN(res, m_pass_Jet_cut2);	

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
  m_pass_vertex_cut      	 = -1;
 
  m_multijet_isSel       	 = -1;
  
  if (m_leadingjet_lorentzvector)
    m_leadingjet_lorentzvector->Clear();
    
  if (m_leadingjetgen_lorentzvector)
    m_leadingjetgen_lorentzvector->Clear();

  if (m_leadingjetraw_lorentzvector)
    m_leadingjetraw_lorentzvector->Clear();
    
  if (m_recoil_lorentzvector)
    m_recoil_lorentzvector->Clear();

  if (m_genrecoil_lorentzvector)
    m_genrecoil_lorentzvector->Clear();

  if (m_pu_lorentzvector)
    m_pu_lorentzvector->Clear();
    
  if (m_jets_recoil_lorentzvector)
    m_jets_recoil_lorentzvector->Clear();

  if (m_jets_pu_lorentzvector)
    m_jets_pu_lorentzvector->Clear();
    
  if (m_jetsgen_recoil_lorentzvector)
    m_jetsgen_recoil_lorentzvector->Clear();

  if (m_met_lorentzvector)
    m_met_lorentzvector->Clear();

  if (m_genmet_lorentzvector)
    m_genmet_lorentzvector->Clear();

  if (m_met_puSubstract_lorentzvector)
    m_met_puSubstract_lorentzvector->Clear();
}

//}

DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, multijetExtractorAnalysis, "multijetExtractorAnalysis");
