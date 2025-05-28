#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1D.h"
#include <vector>
#include <fstream>
#include <map>
#include <iostream>
#include "headers.h"

int main(int argc, char* argv[])
{
  //gSystem->Load("./MyDict_cxx.so");
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <job_index>\n";
    return 1;
  }

  int jobIndex = atoi(argv[1]); // e.g., 0, 1, 2, ...

  Long64_t nEvent = atoi(argv[2]);

  bool DEBUG = false;
  
  ROOT::EnableImplicitMT(4); // Enable Multi-threading
  long int nCountDebug = 500;
  Pythia pythia;
//  pythia.particleData.addParticle(1000080160, "16O", 6, 30, 0, 15.994915);
  pythia.readString("Beams:idA = 1000080160");
  pythia.readString("Beams:idB = 1000080160");
//  pythia.particleData.addParticle(1000791970, "197Au", 6, 158, 0, 196.96657);
//pythia.readString("Beams:idA = 1000791970");
//pythia.readString("Beams:idB = 1000791970");

  pythia.readString("Beams:eCM = 200.");     

  pythia.readString("Beams:frameType = 1");

  pythia.readString("Angantyr:CollisionModel = 2");
  pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
  pythia.readString("HeavyIon:SigFitNGen = 20");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d", 100 + jobIndex));
  
  pythia.readString(Form("Next:numberCount = %ld",nCountDebug));
  if (!pythia.init()) return 1;

  TFile *wf = new TFile(Form("out_job%d.root", jobIndex), "RECREATE");
  double sumW = 0;
  int count=0;
  const int nCentBins=200;
  const int nCent = 10;
  int centbin[nCent+1] = {0,10,20,30,40,50,60,70,80,90,100};
  TH1D* hsumet = new TH1D("hsumet",";Sum E_{T};",1000,0,200);
  TH1D* hcent = new TH1D("hcent",";Centrality (%);",100,0,100);
  TH1D* hncollall = new TH1D("hncollall",";ncoll;",100,0,100);
  TH1D* hncoll[nCent];
  TH1D* hjetptall = new TH1D("hjetptall",";p_{T}^{jet} [GeV];",20,0,30);
  TH1D* hjetpt[nCent];
  TH1D* hjetptallncoll = new TH1D("hjetptallncoll",";p_{T}^{jet} [GeV];",20,0,30);
  TH1D* hjetptncoll[nCent];
  TH1D* hjetptncollg[nCent];
  for(int ic=0;ic<nCent;ic++){
    hncoll[ic] = new TH1D(Form("hncoll_%d",ic),";ncoll;",100,0,100);
    hjetpt[ic] = new TH1D(Form("hjetpt_%d",ic),";p_{T}^{jet} [GeV];",20,0,30);
    hjetptncoll[ic] = new TH1D(Form("hjetptncoll_%d",ic),";p_{T}^{jet} [GeV];",20,0,30);
    hjetptncollg[ic] = new TH1D(Form("hjetptncollg_%d",ic),";p_{T}^{jet} [GeV];",20,0,30);
  }
  auto boundaries = load_boundaries("Centrality_calib.txt");


  SlowJet slowJet( -1, 0.2, 5, 1.3, 2 );

  for(Long64_t iEvent=0;iEvent<nEvent;iEvent++){
    if (!pythia.next()) continue;

    int code = pythia.info.code();
    if ( ( iEvent % nCountDebug ) == 0 ) std::cout << " run " << count << " events, out of " << iEvent << " attempted / " << (double) iEvent/nEvent * 100 << "(%) ... current code " << code << std::endl; 
    slowJet. analyze( pythia.event );
    const double weight = pythia.info.weight();
    
    auto hiPtr = pythia.info.hiInfo;
    auto scPtr = hiPtr->subCollisionsPtr();

    int nCollG=0;
    for (auto sc : *scPtr) {
      if (sc.type == SubCollision::ABS) ++nCollG;
    }

    int nColl = hiPtr->nCollND();
    if(nColl==0) continue;
    
    int particle_n = 0;
    double sum_et=0;
    bool isfired1 = false;
    bool isfired2 = false;
    for (int j = 0; j < pythia.event.size(); ++j) {
      const auto& p = pythia.event[j];
      if (!p.isFinal()) continue;
      float eta = p.eta();
      float pt = p.pT();
      if(pt<0.05) continue;
      float e = p.e();
      float et = calcEt(pt,eta,e);
      if(eta> 2.1 && eta<5.1){sum_et += et; particle_n++;isfired1=true;}
      if(eta< -2.1 && eta>-5.1){sum_et += et; particle_n++;isfired2=true;}
    }

    if(!(sum_et>0)) continue;

    int cent = (nCentBins -1 -get_quantile_bin(sum_et, boundaries))/2;
    std::cout << "cent " << cent << std::endl;
    int centfillbin=-1;
    for(int ic=0; ic<nCent;ic++){
      if(cent >= centbin[ic] && cent<centbin[ic+1]){
        centfillbin = ic;
        break;
      }
      if(cent==100) centfillbin = nCent-1; 
    }
    
    hsumet->Fill(sum_et);
    hncollall->Fill(nColl);
    hncoll[centfillbin]->Fill(nColl);
    hcent->Fill(cent);

    if ( slowJet.sizeJet() == 0 ) continue;
    for(int ij=0; ij<slowJet.sizeJet(); ij++){
      float jetpt = slowJet.pT(ij);
      float jeteta = slowJet.p(ij).eta();
      if(fabs(jeteta)>1.3) continue;
      if(fabs(jetpt)<5) continue;
      hjetpt[centfillbin]->Fill(jetpt);
      hjetptncoll[centfillbin]->Fill(jetpt,nColl);
      hjetptncollg[centfillbin]->Fill(jetpt,nCollG);
      hjetptall->Fill(jetpt);
      hjetptallncoll->Fill(jetpt,nColl);
    }
    count++;
  }

  pythia.stat();
  std::cout << "total events : " << count << std::endl;
  wf->cd();
  hjetptall->Write();
  hjetptallncoll->Write();
  hcent->Write();
  hncollall->Write();
  hsumet->Write();
  for(int i=0; i<nCent;i++){
    hncoll[i]->Write();
    hjetpt[i]->Write();
    hjetptncoll[i]->Write();
    hjetptncollg[i]->Write();
  }
  wf->Close();
  return 0;
}
