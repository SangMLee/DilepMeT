#include "SampleAnalyzer/User/Analyzer/DilepMet.h"
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool DilepMet::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  myEvent = 0;  
  cout << "BEGIN Initialization" << endl; 
  Manager()->AddRegionSelection("Z -> l l");
  
  Manager()->AddCut("le 2");
  Manager()->AddCut("l*l < 0");
  Manager()->AddCut("diel pt cut");
  Manager()->AddCut("dilep mass");
  Manager()->AddCut("dilep pt > 60");
  Manager()->AddCut("jet veto");
  Manager()->AddCut("B-tagged Jet Veto");
  Manager()->AddCut("Met Cut");
  Manager()->AddCut("Dphi");
  Manager()->AddCut("Met Pt");
  Manager()->AddCut("Jet Met");
  Manager()->AddCut("dR(ll) < 1.8");
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void DilepMet::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  cout << "Events Passed  :  " << myEvent << endl; 
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool DilepMet::Execute(SampleFormat& sample, const EventFormat& event)
{
  //Event Weight  
  double myEventWeight;
  if(Configuration().IsNoEventWeight()) myEventWeight=1.;
  else if(event.mc()->weight()!=0.) myEventWeight=event.mc()->weight();
  else
  {
    WARNING << "Found one event with a zero weight. Skipping..." << endmsg;
    return false;
  }
  Manager()->InitializeForNewEvent(myEventWeight);

  if (event.rec()!=0)
  {
    //Containers
    vector<const RecJetFormat*>vetoJets,vetoBJets;
    vector<const RecLeptonFormat*>selectMu,selectElec,vetoMu,vetoElec;
    vector<const RecTauFormat*>vetoTau;
    // Looking through the reconstructed electron collection
    for (unsigned int i=0;i<event.rec()->electrons().size();i++)
    {
      const RecLeptonFormat& elec = event.rec()->electrons()[i];
// ~~~~    // double pfIso = PHYSICS->Isol->tracker->relIsolation(elec,event.rec(),0.4,0.);
    //  if pfIso 
      double pt = elec.pt();
      double eta = fabs(elec.eta());

      if ((pt > 20) && (eta < 2.4))
          selectElec.push_back(&elec);
      if (pt > 10)
          vetoElec.push_back(&elec);
             
    }

    // Looking through the reconstructed muon collection
    for (unsigned int i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& mu = event.rec()->muons()[i];
      //double pfIso = PHYSICS->Isol->tracker->relIsolation(mu,event.rec(),0.4,0.);
      //if (pfIso > 0.15) continue;
      double pt = mu.pt();
      double eta = fabs(mu.eta());
      
      if ((pt > 20) && (eta < 2.4)) 
        selectMu.push_back(&mu);
      if (pt > 5) 
        vetoMu.push_back(&mu); 
    }
    
    // Looking through the reconstructed hadronic tau collection
    for (unsigned int i=0;i<event.rec()->taus().size();i++)
    {
      const RecTauFormat& tau = event.rec()->taus()[i];
      double pt = tau.pt();
      double eta = fabs(tau.eta());

      if ((pt > 18.) && (eta < 2.3))
          vetoTau.push_back(&tau);
    }

    // Looking through the reconstructed jet collection
    for (unsigned int i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];
      double pt = jet.pt();
      double eta = fabs(jet.eta()); 

      if ((pt > 30.) && (eta < 5))
        vetoJets.push_back(&jet); 
      if (jet.btag()){
          if ((pt > 20.) && (eta < 2.4))
            vetoBJets.push_back(&jet);
      }
    }
    
    // MET 
    MALorentzVector pTmiss = event.rec()->MET().momentum();
    double MET = pTmiss.Pt();

    // Sorting Muons and Electrons in Pt order 
    SORTER->sort(selectElec, PTordering);
    SORTER->sort(selectMu, PTordering);
    // Multiple Muon 
    bool isElec;
    if (!Manager()->ApplyCut((selectElec.size() >= 2) || (selectMu.size() >= 2),"le 2")) return true;
    if (selectElec.size() >= 2)  
        isElec = true;
    else
        isElec = false;  

    // Opposite Charge
    const RecLeptonFormat &lep1 = isElec ? *selectElec[0] : *selectMu[0];
    const RecLeptonFormat &lep2 = isElec ? *selectElec[1] : *selectMu[1];
    if (!Manager()->ApplyCut(lep1.charge()*lep2.charge() < 0., "l*l < 0")) return true;
    
    // Leading elec 
    if (isElec){
       if(!Manager()->ApplyCut(lep1.pt() > 25. && lep2.pt() > 20., "diel pt cut")) return true;
    }

    // Dilep Cuts 
    auto dilep = lep1 + lep2;
    if (!Manager()->ApplyCut(fabs(dilep.m() - 91.1876) < 15.,"dilep mass")) return true; 
    
    if (!Manager()->ApplyCut(dilep.pt() > 60., "dilep pt > 60")) return true;            
    
    // Jet Veto
    if (!Manager()->ApplyCut(vetoJets.size() <= 1,"jet veto")) return true;
    
    // B-Jet Veto 
    if (!Manager()->ApplyCut(vetoBJets.size() == 0.,"B-tagged Jet Veto")) return true;
    
    //MET Cuts
   if (!Manager()->ApplyCut(MET > 100.,"Met Cut")) return true;

    //Delta phi cut 
    auto DiPTphi = dilep.dphi_0_pi(pTmiss);
    if (!Manager()->ApplyCut(DiPTphi > 2.6,"Dphi")) return true;
    
    //Met - dilep.Pt 
    if (!Manager()->ApplyCut(fabs((MET - dilep.pt())/ dilep.pt()) < 0.4,"Met Pt")) return true;
    
   //Single jet
    if (vetoJets.size() == 1){
        auto jetPT = vetoJets[0];
        auto Jetphi = jetPT->dphi_0_pi(pTmiss);
        auto DRll = lep1.dr(lep2) ;
        if (!Manager()->ApplyCut(Jetphi > 0.5,"Jet Met")) return true;
        if (!Manager()->ApplyCut(DRll < 1.8, "dR(ll) < 1.8")) return true;
    }
    myEvent ++;
  }
  
  return true;
}

