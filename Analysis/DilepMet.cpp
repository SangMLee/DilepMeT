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
  Manager()->AddCut("dimu pt cut");
  Manager()->AddCut("d pt cut");
  Manager()->AddCut("dilep pt > 60");
  Manager()->AddCut("many jet");
  Manager()->AddCut("b-jet veto");
  Manager()->AddCut("tau veto");
  Manager()->AddCut("mu > 5");

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
    cout << "---------------NEW EVENT-------------------" << endl;
   
    //Containers
    vector<const RecJetFormat*>selectJets,selectBJets;
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

      if ((pt > 18) && (eta < 2.3))
          vetoTau.push_back(&tau);
    }

    // Looking through the reconstructed jet collection
    for (unsigned int i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];
      if (jet.pt() < 20) continue;
      if (fabs(jet.eta()) > 2.4) continue; 
      cout << "----------------------------------" << endl;
      cout << "Jet" << endl;
      cout << "----------------------------------" << endl;
      cout << "jet: index=" << i+1 
           << " charge=" << jet.charge() << endl;
      cout << "px=" << jet.px()
           << " py=" << jet.py()
           << " pz=" << jet.pz()
           << " e="  << jet.e()
           << " m="  << jet.m() << endl;
      cout << "pt=" << jet.pt() 
           << " eta=" << jet.eta() 
           << " phi=" << jet.phi() << endl;
      cout << "b-tag=" << jet.btag()
           << " true b-tag (before eventual efficiency)=" 
           << jet.true_btag() << endl;
      cout << "EE/HE=" << jet.EEoverHE()
           << " ntracks=" << jet.ntracks() << endl;
      cout << endl;
    }

    // Sorting Muons and Electrons in Pt order 
    SORTER->sort(selectElec, PTordering);
    SORTER->sort(selectMu, PTordering);
  //  cout << "No. electons" << selectElec.size()<<endl;
  //  bool multiel;
  //  bool multimu;
  //  if (selectElec.size() >= 2)
  //      multiel = true;
  //  else 
  //      multiel = false; 
  //
  //  if (selectMu.size() >= 2)
  //      multimu = true;
  //      
  //  else
  //      multimu = false; 
  //  cout <<"Multiel: " << multiel << endl;
  //  cout <<"Multimu: " << multimu << endl;
  //
  //  cout << "No. muons" << selectMu.size() << endl;
   // if (!Manager()->ApplyCut(multiel|| multimu,"le 2")) return true;

    // Multiple Muon 
    bool isElec;
    if (!Manager()->ApplyCut((selectElec.size() >= 2) || (selectMu.size() >= 2),"le 2")) return true;
    if (selectElec.size() >= 2)  
        isElec = true;
    else
        isElec = false;  
    cout << "pass" << endl;
    cout << "is electron" << isElec << endl;
    //Charge cut  
    for ( auto& el : selectElec) {
        cout << "electron pt:  " << el->pt() ;
        cout << "electron charge:  " << el->charge() << endl;
    }
    for ( auto& el : event.rec()->electrons()){
        cout << "rec elec pt:  " << el.pt() << "  rec elec charge:  "<< el.charge()<<endl;
    }
    
    const RecLeptonFormat &lep1 = isElec ? *selectElec[0] : *selectMu[0];
    const RecLeptonFormat &lep2 = isElec ? *selectElec[1] : *selectMu[1];
    cout << "selected electron 1 pt:  " << lep1.pt() << "  selec elect charge 1:  " << lep1.charge()<<endl; 
    cout << "selected electron 2 pt:  " << lep2.pt() << "  selec elect charge 2:  " << lep2.charge()<<endl; 
    //const RecLeptonFormat &elec1 = selectElec[0], &elec2 = selectElec[1];
    //const RecLeptonFormat &mu1 = selectMu[0], &mu2 = selectMu[1];
    if (!Manager()->ApplyCut(lep1.charge()*lep2.charge() < 0., "l*l < 0")) return true; 

    if(!Manager()->ApplyCut(lep1.pt() > 25. && lep2.pt() > 20., "  ")) return true; 
    //Leading lepton pt cut
   // if (selectElec[0].pt() < 25) return false;
   // if (selectMu[0].pt() < 20) return false;

   // const RecLeptonFormat &mu1 = selectMu[0], &mu2 = selectMu[1];
   // if (mu1.charge()*mu2.charge() > 0) return false;
    
   // MALorentzVector pTmiss = event.rec()->MET().momentum();
   // double MET = pTmiss.Pt();
    // Transverse missing energy (MET)
    cout << "MET pt=" << event.rec()->MET().pt()
         << " phi=" << event.rec()->MET().phi() << endl;
    cout << endl;
   // if (MET <= 100){
      //  return false;
   // }
    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.rec()->MHT().pt()
              << " phi=" << event.rec()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.rec()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.rec()->THT() << endl;
    cout << endl;
    myEvent ++; 
  }
  
  return true;
}

