#include "SampleAnalyzer/User/Analyzer/darkmatter.h"
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool darkmatter::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  Manager()->AddRegionSelection("Z -> l l");
  Manager()->AddCut("le 2");
  Manager()->AddCut("opp charge");
  Manager()->AddCut("lep pt" );
  Manager()->AddCut("z veto");
  Manager()->AddCut("dilep pt");
  Manager()->AddCut("jet veto");
  Manager()->AddCut("bjet veto");
  Manager()->AddCut("pu veto");
  Manager()->AddCut("tau veto");
  Manager()->AddCut("met");
  Manager()->AddCut("dphi llmet");
  Manager()->AddCut("met pt");
  Manager()->AddCut("dphi jetmet");
  Manager()->AddCut("dr");
  cout << "END   Initialization" << endl;
  nevent =0;
  ee =0;
  mm =0;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void darkmatter::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  cout << "END   Finalization" << endl;
  cout << nevent <<endl;
  cout << "ee" << ee << "  mm" << mm <<endl;
}

double darkmatter::reliso04(const EventFormat& event, const RecLeptonFormat& lep){
      return PHYSICS->Isol->eflow->relIsolation(lep,event.rec(),0.4,0.,IsolationEFlow::ALL_COMPONENTS);
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool darkmatter::Execute(SampleFormat& sample, const EventFormat& event)
{

  double myEventWeight;
  if(Configuration().IsNoEventWeight()) myEventWeight=1.;
  else if(event.mc()->weight()!=0.) myEventWeight=event.mc()->weight();
  else
  {
    WARNING << "Found one event with a zero weight. Skipping..." << endmsg;
    return false;
  }
  Manager()->InitializeForNewEvent(myEventWeight);

  // ***************************************************************************
  // Example of analysis with reconstructed objects
  // Concerned samples : 
  //   - LHCO samples
  //   - LHE/STDHEP/HEPMC samples after applying jet-clustering algorithm
  // ***************************************************************************
  if (event.rec()!=0)
  {
    cout << "---------------NEW EVENT-------------------" << endl;
    selectedels.clear();
    selectedmus.clear();
    selectedjets.clear();

    // lepton selection
    for (unsigned int i=0;i<event.rec()->electrons().size();i++)
    {
      const RecLeptonFormat& elec = event.rec()->electrons()[i];
      if (elec.pt()<20) continue;
      if (fabs(elec.eta())>2.4) continue;
      if ( fabs(elec.eta()) < 1.479 && reliso04(event, elec) > 0.077 ) continue;
      if ( fabs(elec.eta()) > 1.479 && reliso04(event, elec) > 0.068 ) continue;
      selectedels.push_back(&elec);
    }

    for (unsigned int i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& mu = event.rec()->muons()[i];
      if (mu.pt()<20) continue;
      if (fabs(mu.eta())>2.4) continue;
      if (reliso04(event,mu) >0.15) continue;
      selectedmus.push_back(&mu);
    }

    // evnet selection
    if (!Manager()->ApplyCut((selectedels.size() >= 2) || (selectedmus.size() >= 2),"le 2")) return true;

    std::vector<const RecLeptonFormat*> leps;
    string channel;

    SORTER->sort(selectedels, PTordering);
    SORTER->sort(selectedmus, PTordering);

    // get dilep
    if ( selectedels.size()>=2 && selectedmus.size()<2 ) {
	channel = "ee";
	leps.push_back(selectedels[0]);
	leps.push_back(selectedels[1]);
	selectedels.erase(selectedels.begin(),selectedels.begin()+2);
	ee +=1;
    }
    else if ( selectedmus.size()>=2 && selectedels.size()<2 ) {
	channel = "mm";
	leps.push_back(selectedmus[0]);
	leps.push_back(selectedmus[1]);
	selectedmus.erase(selectedmus.begin(),selectedmus.begin()+2);
	mm +=1;
    }
    else if ( selectedmus.size()>=2 && selectedels.size()>=2 ) {
	cout << "there're more than 4 leptons!" << endl; return true;
    }
    else return true;
    auto dilep = *leps[0]+*leps[1];

    for (unsigned int i=0;i<event.rec()->jets().size();i++) {
      const RecJetFormat& jet = event.rec()->jets()[i];
      if (jet.pt()<20) continue;
      if (fabs(jet.eta())>5) continue;
      if (jet.dr(leps[0])<0.4) continue;
      if (jet.dr(leps[1])<0.4) continue;
      selectedjets.push_back(&jet); 
    }

    if (!Manager()->ApplyCut(leps[0]->charge()*leps[1]->charge() < 0,"opp charge")) return true;

    bool passElectron = channel == "ee" && (leps[0]->pt()>25 && leps[1]->pt()>20);
    bool passMuon = channel == "mm" && (leps[0]->pt()>20 && leps[1]->pt()>20);
    //if (!Manager()->ApplyCut(passElectron,"lep pt" )) return true;
    if (!Manager()->ApplyCut(passElectron||passMuon, "lep pt" )) return true;

    if (!Manager()->ApplyCut(fabs(dilep.m()-Zmass) < 15, "z veto" )) return true;
    if (!Manager()->ApplyCut(dilep.pt()>60, "dilep pt" )) return true;

    int jetveto = 0;
    bool bjetveto = false;
    for (unsigned int i=0;i<event.rec()->jets().size();i++){
        const RecJetFormat& jet = event.rec()->jets()[i];
	if (jet.pt() >30) jetveto += 1;
        if (jet.btag() && jet.pt()>20 && fabs(jet.eta())<2.4 ) bjetveto = true;
    }
    if (!Manager()->ApplyCut(jetveto>1, "jet veto" )) return true;
    if (!Manager()->ApplyCut(!bjetveto, "bjet veto" )) return true;

    bool puveto = false;
    for (auto& e : selectedels) { cout << e->pt() << endl; if (e->pt() > 10) puveto=true; } 
    for (auto& m : selectedmus) { cout << m->pt() << endl; if (m->pt() > 10) puveto=true; } 
    if (!Manager()->ApplyCut(!puveto,"pu veto" )) return true;

    bool tauveto = false;
    for (unsigned int i=0;i<event.rec()->taus().size();i++)
    {
      const RecTauFormat& tau = event.rec()->taus()[i];
      if (tau.pt() > 18) tauveto = true;
    }
    if (!Manager()->ApplyCut(!tauveto, "tau veto" )) return true;

    if (!Manager()->ApplyCut(event.rec()->MET().pt() > 100, "met" )) return true;

    MALorentzVector pTmiss = event.rec()->MET().momentum();
    auto dphiLepsMet = dilep.dphi_0_pi(pTmiss);
    if (!Manager()->ApplyCut(dphiLepsMet > 2.6, "dphi llmet")) return true;
    if (!Manager()->ApplyCut(fabs(event.rec()->MET().pt()-dilep.pt())/dilep.pt() < 0.4, "met pt")) return true;

    if (selectedjets.size() >0){
      auto dphiJetMet = selectedjets[0]->dphi_0_pi(pTmiss);
      if (!Manager()->ApplyCut(dphiJetMet > 0.5, "dphi jetmet")) return true;
      if (!Manager()->ApplyCut(leps[0]->dr(leps[1]) < 1.8, "dr")) return true;
    }

    nevent += 1;
  }
  return true;
}

