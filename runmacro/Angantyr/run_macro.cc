#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include <vector>
#include <fstream>
#include <map>
#include <iostream>
#include "headers.h"
#include "TMath.h"
      

bool removeparticle(double pt, int cent, TF1* func)
{
  if(pt<5) return false;
  static TRandom3 rndm(0);  
  float surviveprob = func->Eval(cent);
  return (rndm.Uniform(0,1) > surviveprob);
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <job_index>\n";
    return 1;
  }

  int jobIndex = atoi(argv[1]); // e.g., 0, 1, 2, ...

  Long64_t nEvent = atoi(argv[2]);
  std::string acctype= argv[3];
  std::string vartype = argv[4];
  std::string collisionenergy;
  float etalow = -1; float etahigh =-1;

  float jetPtMinCut = -1;
  float jetEtaMaxCut = 1.3;
  float hadronEtaMaxCut = 0.5;
  float hadronPtMinCut = 0.2;
  float hadronJet_hadronEtaMaxCut = 0.5;
  float hadronJet_jetEtaMaxCut = 1.3;
  float hadronJetPtMinCut = 7;
  float hadronJetPtMaxCut = 30;
  float hadronJet_deltaPhiCut  = 1./4*M_PI;

  bool use_quenching{false};
  use_quenching = atoi(argv[5]);

  TF1 *f_quench = new TF1("f_quench",quenchFunct,0,100,1);

  if(acctype == "LHC"){
    collisionenergy = "5360.";
    etalow = 3.2; etahigh = 4.9;
    jetEtaMaxCut = 2.1;
    jetPtMinCut = 10;
    hadronPtMinCut = 0.2;
    hadronEtaMaxCut=2.5;
    hadronJet_hadronEtaMaxCut = 0.9;
    hadronJet_jetEtaMaxCut = 0.7;
    hadronJetPtMinCut = 12;
    hadronJetPtMaxCut = 50;
    hadronJet_deltaPhiCut  = 0.6;
  }
  else if(acctype == "RHIC"){
    collisionenergy = "200.";
    etalow = 2.1; etahigh = 5.1;
    jetEtaMaxCut = 1.3;
    jetPtMinCut = 5;
    hadronPtMinCut = 0.2;
    hadronJet_hadronEtaMaxCut = 1.5;
    hadronJet_jetEtaMaxCut = 1.3;
    hadronJetPtMinCut = 7;
    hadronJetPtMaxCut = 30;
    hadronJet_deltaPhiCut  = 1./4*M_PI;
  }
  else{
    std::cerr << "Collision type error. User either LHC or RHIC" << std::endl;
    return 1;
  }

  if(vartype != "nch" && vartype !="sumet"){
    std::cerr << "centrality calibration variable not set" << std::endl;
    return 1;
  }
  std::string outputfilename = Form("out_%s_%s_job%d.root",acctype.c_str(), vartype.c_str(), jobIndex);
  std::string outputdirectory = (use_quenching) ? Form("output%s_quench",acctype.c_str()) : Form("output%s",acctype.c_str());
  TFile *wf = new TFile(Form("%s/%s",outputdirectory.c_str(), outputfilename.c_str()),"RECREATE");
  int count=0;
  const int nCentBins=200;
  const int nCent = 10;
  int centbin[nCent+1] = {0,10,20,30,40,50,60,70,80,90,100};
 
 /* 
  double ptjet_binning[] = {
    // 0–20 GeV in 1.0
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
    6.0, 7.0, 8.0, 9.0, 10.0,
    11.0, 12.0, 13.0, 14.0, 15.0,
    16.0, 17.0, 18.0, 19.0, 20.0,
    // 20–30 GeV in 2.0
    22.0, 24.0, 26.0, 28.0, 30.0,
    // 30–42 GeV in 4.0
    34.0, 38.0, 42.0,
    // 42–50 GeV in 8.0
    50.0,
    // 50–60 GeV in 10.0
    60.0
  };
  const int nPtJetBins = sizeof(ptjet_binning)/sizeof(ptjet_binning[0]) - 1;  

  double ptcharge_binning[] = {
    // 0–10 GeV in 0.5
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
    5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    // 10–16 GeV in 1.0
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
    // 16–24 GeV in 2.0
    18.0, 20.0, 22.0, 24.0,
    // 24–40 GeV in 4.0
    28.0, 32.0, 36.0, 40.0,
    // 40–60 GeV in 10.0
    50.0, 60.0
  };
  const int nPtChargeBins = sizeof(ptcharge_binning)/sizeof(ptcharge_binning[0]) - 1;
*/
  
  float ptjetmin = 0;
  float ptjetmax = 300;
  float ptchargemin = 0;
  float ptchargemax = 100;
  const int nPtJetBins = 300;
  const int nPtChargeBins = 200;
  const int nBinsNcoll = 300;

  TH1D* hsumet = new TH1D("hsumet",";Sum E_{T};",1000,0,200);
  TH1D* hcent = new TH1D("hcent",";Centrality (%);",100,0,100);
  TH1D* hncollall = new TH1D("hncollall",";ncoll;",nBinsNcoll, 0, nBinsNcoll);
  TH1D* hncollallcharge = new TH1D("hncollallcharge",";ncoll;",nBinsNcoll, 0, nBinsNcoll);
  TH1D* hncoll[nCent];
  TH1D* hncollcharge[nCent];
  TH1D* hjetptall = new TH1D("hjetptall",";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin, ptjetmax);
  TH1D* hjetptall_charged = new TH1D("hjetptall_charged",";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin, ptjetmax);
  TH1D* hchargeptall = new TH1D("hchargeptall",";p_{T}^{charge} [GeV];",nPtChargeBins,ptchargemin,ptchargemax);
  TH1D* hchargejetptall = new TH1D("hchargejetptall",";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin, ptjetmax);
  TH1D* hchargejetptall_charged = new TH1D("hchargejetptall_charged",";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin, ptjetmax);
  TH1D* hjetpt[nCent];
  TH1D* hjetpt_charged[nCent];
  TH2D* h2jetptncollref =  new TH2D("h2jetncollref",";N_{coll};p_{T}^{jet} [GeV];",nBinsNcoll, 0, nBinsNcoll,nPtJetBins,ptjetmin,ptjetmax);
  TH2D* h2jetptncollref_charged =  new TH2D("h2jetncollref_charged",";N_{coll};p_{T}^{jet} [GeV];",nBinsNcoll, 0, nBinsNcoll,nPtJetBins,ptjetmin,ptjetmax);
  TH2D* h2chargeptncollref =  new TH2D("h2chargencollref",";N_{coll};p_{T}^{charge} [GeV];",nBinsNcoll, 0, nBinsNcoll,nPtChargeBins,ptchargemin,ptchargemax);
  TH2D* h2chargejetptncollref =  new TH2D("h2chargejetncollref",";N_{coll};p_{T}^{jet} [GeV];",nBinsNcoll, 0, nBinsNcoll,nPtJetBins,ptjetmin,ptjetmax);
  TH2D* h2chargejetptncollref_charged =  new TH2D("h2chargejetncollref_charged",";N_{coll};p_{T}^{jet} [GeV];",nBinsNcoll, 0, nBinsNcoll,nPtJetBins,ptjetmin,ptjetmax);

  TH1D* hchargept[nCent];
  TH1D* hchargejetpt[nCent];
  TH1D* hchargejetpt_charged[nCent];

  for(int ic=0;ic<nCent;ic++){
    hncoll[ic] = new TH1D(Form("hncoll_%d",ic),";ncoll;",nBinsNcoll, 0, nBinsNcoll);
    hncollcharge[ic] = new TH1D(Form("hncollcharge_%d",ic),";ncoll;",nBinsNcoll, 0, nBinsNcoll);
    hjetpt[ic] = new TH1D(Form("hjetpt_%d",ic),";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin,ptjetmax);
    hjetpt_charged[ic] = new TH1D(Form("hjetpt_charged_%d",ic),";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin,ptjetmax);
    hchargept[ic] = new TH1D(Form("hchargept_%d",ic),";p_{T}^{charge} [GeV];",nPtChargeBins,ptchargemin,ptchargemax);
    hchargejetpt[ic] = new TH1D(Form("hchargejetpt_%d",ic),";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin,ptjetmax);
    hchargejetpt_charged[ic] = new TH1D(Form("hchargejetpt_charged_%d",ic),";p_{T}^{jet} [GeV];",nPtJetBins,ptjetmin,ptjetmax);
  }

  //centrality boundary loading
  auto boundaries = load_boundaries(Form("Centrality_calib_%s_%s.txt",acctype.c_str(), vartype.c_str()));

  bool DEBUG = false;
  
  ROOT::EnableImplicitMT(4); // Enable Multi-threading
  long int nCountDebug = 500;
  Pythia pythia;
  pythia.readString("Beams:idA = 1000080160");
  pythia.readString("Beams:idB = 1000080160");

//  pythia.particleData.addParticle(1000080160, "16O", 6, 30, 0, 15.994915);
//  pythia.particleData.addParticle(1000791970, "197Au", 6, 158, 0, 196.96657);
//pythia.readString("Beams:idA = 1000791970");
//pythia.readString("Beams:idB = 1000791970");

  pythia.readString(Form("Beams:eCM = %s",collisionenergy.c_str()));     

  pythia.readString("Beams:frameType = 1");

  pythia.readString("Angantyr:CollisionModel = 2");
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  pythia.readString("HeavyIon:SigFitDefPar = 2.15,17.24,0.33");
  pythia.readString("HeavyIon:SigFitNGen = 20");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d", 100 + jobIndex));
  
  pythia.readString(Form("Next:numberCount = %ld",nCountDebug));
  if (!pythia.init()) return 1;

  long int countmultMB[10] = {0};
  long int countmultMBcharge[10] = {0};

  SlowJet slowJetCharge( -1, 0.2, jetPtMinCut, jetEtaMaxCut, 3);
  SlowJet slowJetFull( -1, 0.2, jetPtMinCut, jetEtaMaxCut, 2);

  TRandom3 *rnd = new TRandom3();
  rnd->SetSeed(jobIndex);
  const double nparres=0.1;

  for(Long64_t iEvent=0;iEvent<nEvent;iEvent++){
    if (!pythia.next()) continue;

    int code = pythia.info.code();
    if ( ( iEvent % nCountDebug ) == 0 ) std::cout << " run " << count << " events, out of " << iEvent << " attempted / " << (double) iEvent/nEvent * 100 << "(%) ... current code " << code << std::endl; 
    const double weight = pythia.info.weight();
    
    auto hiPtr = pythia.info.hiInfo;
    auto scPtr = hiPtr->subCollisionsPtr();

    int nCollG=0;
    for (auto sc : *scPtr) {
      if (sc.type == SubCollision::ABS) ++nCollG;
    }

    int nColl = hiPtr->nCollND();
    if(nColl<=0) continue;
    
    
    int particle_n = 0;
    double sum_et=0;
  
    for (int j = 0; j < pythia.event.size(); ++j) {
      const auto& p = pythia.event[j];
      if (!p.isFinal()) continue;
      float eta = p.eta();
      float pt = p.pT();
      if(pt<0.05) continue;
      float e = p.e();
      float et = calcEt(pt,eta,e);
      if(fabs(eta)< etalow || fabs(eta)>etahigh) continue;
      if(p.isCharged()) particle_n++;

      sum_et += et; 
    }

    int cent =-1;
    if(vartype == "nch"){
      double nparsmear = rnd->Gaus(particle_n, sqrt(particle_n) * nparres);
      cent = (nCentBins -1 -get_quantile_bin(nparsmear, boundaries))/2;
    }
    else if(vartype == "sumet"){
     cent = (nCentBins -1 -get_quantile_bin(sum_et, boundaries))/2;
    }
    int centfillbin=-1;
    for(int ic=0; ic<nCent;ic++){
      if(cent >= centbin[ic] && cent<centbin[ic+1]){
        centfillbin = ic;
        break;
      }
      if(cent==100) centfillbin = nCent-1; 
    }
    if (centfillbin < 0 || centfillbin > 10) continue;
   
    Event filteredEvent; 
    if(use_quenching)
    {
      filteredEvent.init();
      filteredEvent.append(pythia.event[0]);
      filteredEvent.append(pythia.event[1]);
      for (int j = 0; j < pythia.event.size(); ++j) {
        const auto& p = pythia.event[j];
        if (!p.isFinal()) continue;
        float eta = p.eta();
        float pt = p.pT();
        if(removeparticle(pt, cent, f_quench)) continue;
        filteredEvent.append(p);
      }
      slowJetFull.analyze( filteredEvent );
      slowJetCharge.analyze( filteredEvent );
    }
    else{
      slowJetFull.analyze( pythia.event );
      slowJetCharge.analyze( pythia.event );
    }
    
    hsumet->Fill(sum_et);
    hncollall->Fill(nColl);
    hncoll[centfillbin]->Fill(nColl);
    hcent->Fill(cent);
    countmultMB[centfillbin]++;

    
    int eventsize = (use_quenching) ? filteredEvent.size() : pythia.event.size();
    for (int j = 0; j < eventsize; ++j) {
      const auto& p = (use_quenching) ? filteredEvent[j] : pythia.event[j];
      if (!p.isFinal()) continue;
      if (!p.isCharged()) continue;
      float eta = p.eta();
      float pt = p.pT();
      if(pt<hadronPtMinCut) continue;
      if(fabs(eta)>hadronEtaMaxCut) continue;
      hchargept[centfillbin]->Fill(pt);
      hchargeptall->Fill(pt);
      h2chargeptncollref->Fill(nColl,pt);

      //For I_CP
      if (pt>hadronJetPtMinCut && pt<hadronJetPtMaxCut && fabs(eta)< hadronJet_hadronEtaMaxCut) {
        countmultMBcharge[centfillbin]++;
        hncollallcharge->Fill(nColl);
        hncollcharge[centfillbin]->Fill(nColl);
        if ( slowJetFull.sizeJet() == 0  && slowJetCharge.sizeJet() == 0) continue;
        for(int ij=0; ij<slowJetFull.sizeJet(); ij++){
          float jetpt = slowJetFull.pT(ij);
          float jeteta = slowJetFull.p(ij).eta();
          if(fabs(jeteta)>hadronJet_jetEtaMaxCut) continue;
          if(fabs(jetpt)<jetPtMinCut) continue;
          float jetphi = slowJetFull.phi(ij);
          float chargephi = p.phi();
          float deltaphi = fabs(chargephi - jetphi);
          if(deltaphi > M_PI) deltaphi = 2*M_PI - deltaphi;
          if(fabs(deltaphi - M_PI) < hadronJet_deltaPhiCut) continue;
          hchargejetpt[centfillbin]->Fill(jetpt);
          hchargejetptall->Fill(jetpt);
          h2chargejetptncollref->Fill(nColl,jetpt);
        }
        for(int ij=0; ij<slowJetCharge.sizeJet(); ij++){
          float jetpt = slowJetCharge.pT(ij);
          float jeteta = slowJetCharge.p(ij).eta();
          if(fabs(jeteta)>hadronJet_jetEtaMaxCut) continue;
          if(fabs(jetpt)<jetPtMinCut) continue;
          float jetphi = slowJetCharge.phi(ij);
          float chargephi = p.phi();
          float deltaphi = fabs(chargephi - jetphi);
          if(deltaphi > M_PI) deltaphi = 2*M_PI - deltaphi;
          if(fabs(deltaphi - M_PI) < hadronJet_deltaPhiCut) continue;
          hchargejetpt_charged[centfillbin]->Fill(jetpt);
          hchargejetptall_charged->Fill(jetpt);
          h2chargejetptncollref_charged->Fill(nColl,jetpt);
        }
      }
    }
    
    //Fill general jets
    if ( slowJetFull.sizeJet() == 0  && slowJetCharge.sizeJet() == 0) continue;
    double jetmaxb = 0;
    int jetmaxbidx = -1;
    for(int ij=0; ij<slowJetFull.sizeJet(); ij++){
      float jetpt = slowJetFull.pT(ij);
      float jeteta = slowJetFull.p(ij).eta();
      if(fabs(jeteta)<jetEtaMaxCut && jetpt > jetPtMinCut){
        if(jetmaxb < jetpt){
          jetmaxb = jetpt;
          jetmaxbidx = ij;
        }
        hjetpt[centfillbin]->Fill(jetpt);
        hjetptall->Fill(jetpt);
        h2jetptncollref->Fill(nColl,jetpt);
      }
    }

    //charged jet
    for(int ij=0; ij<slowJetCharge.sizeJet(); ij++){
      float jetpt = slowJetCharge.pT(ij);
      float jeteta = slowJetCharge.p(ij).eta();
      if(fabs(jeteta)>jetEtaMaxCut) continue;
      if(fabs(jetpt)<jetPtMinCut) continue;
      hjetpt_charged[centfillbin]->Fill(jetpt);
      hjetptall_charged->Fill(jetpt);
      h2jetptncollref_charged->Fill(nColl,jetpt);
    }
    count++;
  }

  pythia.stat();
  std::cout << "total events : " << count << std::endl;
  
  TH1D* hcountmb = new TH1D("hcountmb",";;",10,0,10);
  TH1D* hcountmbcharge = new TH1D("hcountmbcharge",";;",10,0,10);
  wf->cd();
  hjetptall->Write();
  hjetptall_charged->Write();
  hchargeptall->Write();
  hchargejetptall->Write();
  hcent->Write();
  hncollall->Write();
  hncollallcharge->Write();
  hsumet->Write();
  for(int i=0; i<nCent;i++){
    hncoll[i]->Write();
    hncollcharge[i]->Write();
    hjetpt[i]->Write();
    hjetpt_charged[i]->Write();
    hchargept[i]->Write();
    
    hchargejetpt[i]->Write();
    hchargejetpt_charged[i]->Write();
    hcountmb->SetBinContent(i+1, countmultMB[i]); // to be used for counting events in each centrality 
    hcountmbcharge->SetBinContent(i+1, countmultMBcharge[i]);
  }
  hcountmb->Write();
  hcountmbcharge->Write();
  h2jetptncollref->Write();
  h2chargeptncollref->Write();
  h2chargejetptncollref->Write();
  h2jetptncollref_charged->Write();
  h2chargejetptncollref_charged->Write();

  wf->Close();
  return 0;
}
