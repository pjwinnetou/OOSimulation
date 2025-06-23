#include "anaheaders.h"
void drawPlot(bool isRHIC = true, bool isAngantyr = false,std::string cent_type="nch")
{
  double xmaxcharge = (isRHIC) ? 30 : 80;
  double xmaxjet = (isRHIC) ? 50 : 250;
  double hadronetacut = (isRHIC) ? 0.5 : 2.5;
  double hadronetacutForICP = (isRHIC) ? 1.5 : 0.9;
  double jetetacut = (isRHIC) ? 1.3 : 2.1;
  double jetetacutForICP = (isRHIC) ? 1.3 : 0.7;
  float etalow_cent = (isRHIC) ? 2.1 : 3.2;
  float etahigh_cent = (isRHIC) ? 5.1 : 4.9;
  float hadronptlowForICP = (isRHIC) ? 7 : 12;
  float hadronpthighForICP = (isRHIC) ? 30 : 50;

  double xpos = 0.451;
  double ypos = 0.82;
  double ydiff = 0.057;
  double xshift = 0.26;
  
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  std::string collen = (isRHIC) ? "200GeV" : "5p36TeV";
  std::string colltxt = (isRHIC) ? "200 GeV" : "5.36 TeV";
  std::string simtype = (isAngantyr) ? "Angantyr" : "Hijing";
  std::string simtxt = (isAngantyr) ? "Pythia-Angantyr" : "Hijing";
  std::string cent_type_txt = (cent_type=="nch") ? "N_{ch}" : "Sum E_{T}";
    
  TFile *f1 = new TFile(Form("../histfiles/output_%s_OO_%s_%s_merged.root",simtype.c_str(),collen.c_str(),cent_type.c_str()),"read");
  TH1D* hjetpt[10];
  TH1D* hjetpt_charged[10];
  TH1D* hchargept[10];
  TH1D* hchargejetpt[10];
  TH1D* hchargejetpt_charged[10];
  TH1D* hncoll[10];
  TH1D* hncollcharge[10];

  TH1D* hsumet = (TH1D*) f1->Get("hsumet");
  TH1D* hncollall = (TH1D*) f1->Get("hncollall");
  TH1D* hncollallcharge = (TH1D*) f1->Get("hncollallcharge");
  TH1D* hcent = (TH1D*) f1->Get("hcent");
  TH1D* hcountmb = (TH1D*) f1->Get("hcountmb"); 
  TH1D* hcountmbcharge = (TH1D*) f1->Get("hcountmbcharge"); 
  for(int i=0;i<10;i++){
    hchargept[i] = (TH1D*)f1->Get(Form("hchargept_%d",i));
    hchargejetpt[i] = (TH1D*)f1->Get(Form("hchargejetpt_%d",i));
    hchargejetpt_charged[i] = (TH1D*)f1->Get(Form("hchargejetpt_charged_%d",i));
    hjetpt[i] = (TH1D*)f1->Get(Form("hjetpt_%d",i));
    hjetpt_charged[i] = (TH1D*)f1->Get(Form("hjetpt_charged_%d",i));
    hncoll[i] = (TH1D*) f1->Get(Form("hncoll_%d",i));
    hncollcharge[i] = (TH1D*) f1->Get(Form("hncollcharge_%d",i));
    ncollmean[i] = hncoll[i]->GetMean();
    nevent[i] = hcountmb->GetBinContent(i+1); 
    neventall += nevent[i];
    neventcharge[i] = hcountmbcharge->GetBinContent(i+1); 
  }

  //2D Normalization
  TH2D* h2jetncollref = (TH2D*) f1->Get("h2jetncollref"); 
  TH2D* h2jetncollref_charged = (TH2D*) f1->Get("h2jetncollref_charged"); 
  TH2D* h2chargencollref= (TH2D*) f1->Get("h2chargencollref"); 
  TH2D* h2chargejetncollref = (TH2D*) f1->Get("h2chargejetncollref"); 
  TH2D* h2chargejetncollref_charged = (TH2D*) f1->Get("h2chargejetncollref_charged"); 

  NormBias2D(h2jetncollref,hncollall);
  NormBias2D(h2jetncollref_charged,hncollall);
  NormBias2D(h2chargencollref,hncollall);
  NormBias2D(h2chargejetncollref,hncollallcharge);
  NormBias2D(h2chargejetncollref_charged,hncollallcharge);

  //Specific for 60-80%
  double nevent6080 = nevent[6] + nevent[7]; 
  double neventcharge6080 = neventcharge[6] + neventcharge[7]; 
  TH1D* hjetpt6080 = (TH1D*) hjetpt[6]->Clone("hjet6080");
  hjetpt6080->Add(hjetpt[7]);
  TH1D* hjetpt6080_charged = (TH1D*) hjetpt_charged[6]->Clone("hjet6080_charged");
  hjetpt6080_charged->Add(hjetpt_charged[7]);

  TH1D* hchargept6080 = (TH1D*) hchargept[6]->Clone("hcharge6080");
  hchargept6080->Add(hchargept[7]);

  TH1D* hchargejetpt6080 = (TH1D*) hchargejetpt[6]->Clone("hchargejet6080");
  hchargejetpt6080->Add(hchargejetpt[7]);
  TH1D* hchargejetpt6080_charged = (TH1D*) hchargejetpt_charged[6]->Clone("hchargejet6080_charged");
  hchargejetpt6080_charged->Add(hchargejetpt_charged[7]);

  TH1D* hncoll6080 = (TH1D*) hncoll[6]->Clone("hncoll6080");
  hncoll6080->Add(hncoll[7]);
  double ncoll6080 = hncoll6080->GetMean();

  TH1D* hncollcharge6080 = (TH1D*) hncollcharge[6]->Clone("hncollcharge6080");
  hncollcharge6080->Add(hncollcharge[7]);
  double ncollcharge6080 = hncollcharge6080->GetMean();

  //Make Truth Hist
  TH1D* hjetpt6080truth = (TH1D*) hjetpt6080->Clone("hjetpt6080truth");
  TH1D* hjetpt6080truth_charged = (TH1D*) hjetpt6080_charged->Clone("hjetpt6080truth_charged");
  TH1D* hchargept6080truth = (TH1D*) hchargept6080->Clone("hchargept6080truth");
  TH1D* hchargejetpt6080truth = (TH1D*) hchargejetpt6080->Clone("hchargejetpt6080truth");
  TH1D* hchargejetpt6080truth_charged = (TH1D*) hchargejetpt6080->Clone("hchargejetpt6080truth_charged");
  MakeTruthHist(hjetpt6080truth, h2jetncollref, hncoll6080);
  MakeTruthHist(hjetpt6080truth_charged, h2jetncollref_charged, hncoll6080);
  MakeTruthHist(hchargept6080truth, h2chargencollref, hncoll6080);
  MakeTruthHist(hchargejetpt6080truth, h2chargejetncollref, hncollcharge6080);
  MakeTruthHist(hchargejetpt6080truth_charged, h2chargejetncollref_charged, hncollcharge6080);

  TH1D* hjetpt_truth[10];
  TH1D* hjetpt_truth_charged[10];
  TH1D* hchargept_truth[10];
  TH1D* hchargejetpt_truth[10];
  TH1D* hchargejetpt_truth_charged[10];

  TH1D* hchargept_truth_all = (TH1D*) hchargept[0]->Clone("hchargept_truth_all");
  TH1D* hchargept_truth_all_fill = (TH1D*) hchargept[0]->Clone("hchargept_truth_all");

  for(int i=0; i<10; i++){
    hjetpt_truth[i] = (TH1D*) hjetpt[i]->Clone(Form("hjetpt_truth%d",i));
    hjetpt_truth_charged[i] = (TH1D*) hjetpt_charged[i]->Clone(Form("hjetpt_truth_charged_%d",i));
    hchargept_truth[i] = (TH1D*) hchargept[i]->Clone(Form("hchargept_truth%d",i));
    hchargejetpt_truth[i] = (TH1D*) hchargejetpt[i]->Clone(Form("hchargejetpt_truth%d",i));
    hchargejetpt_truth_charged[i] = (TH1D*) hchargejetpt_charged[i]->Clone(Form("hchargejetpt_truth_charged_%d",i));
    
    MakeTruthHist(hjetpt_truth[i], h2jetncollref, hncoll[i]);
    MakeTruthHist(hjetpt_truth_charged[i], h2jetncollref_charged, hncoll[i]);
    MakeTruthHist(hchargept_truth[i], h2chargencollref, hncoll[i]);
    MakeTruthHist(hchargejetpt_truth[i], h2chargejetncollref, hncollcharge[i]);
    MakeTruthHist(hchargejetpt_truth_charged[i], h2chargejetncollref_charged, hncollcharge[i]);
    
    //Rebin after threshold
    if(isRHIC){
      PartialRebinCountsFill(hjetpt[i], 10,2);
      PartialRebinCountsFill(hjetpt_charged[i], 10,2);
      PartialRebinCountsFill(hchargejetpt[i], 10,2);
      PartialRebinCountsFill(hchargejetpt_charged[i], 10,2);

      PartialRebinCountsFill(hjetpt_truth[i], 10,2);
      PartialRebinCountsFill(hjetpt_truth_charged[i], 10,2);
      PartialRebinCountsFill(hchargejetpt_truth[i], 10,2);
      PartialRebinCountsFill(hchargejetpt_truth_charged[i], 10,2);
      
      PartialRebinCountsFill(hjetpt[i], 30,2);
      PartialRebinCountsFill(hjetpt_charged[i], 30,2);
      PartialRebinCountsFill(hchargejetpt[i], 30,2);
      PartialRebinCountsFill(hchargejetpt_charged[i], 30,2);

      PartialRebinCountsFill(hjetpt_truth[i], 30,2);
      PartialRebinCountsFill(hjetpt_truth_charged[i], 30,2);
      PartialRebinCountsFill(hchargejetpt_truth[i], 30,2);
      PartialRebinCountsFill(hchargejetpt_truth_charged[i], 30,2);

      //charged
      PartialRebinCountsFill(hchargept[i], 10,2);
      PartialRebinCountsFill(hchargept_truth[i], 10,2);
      PartialRebinCountsFill(hchargept[i], 20,2);
      PartialRebinCountsFill(hchargept_truth[i], 20,2);

      if(i==0){
        PartialRebinCountsFill(hjetpt6080, 10,2);
        PartialRebinCountsFill(hjetpt6080_charged, 10,2);
        PartialRebinCountsFill(hjetpt6080truth, 10,2);
        PartialRebinCountsFill(hjetpt6080truth_charged, 10,2);
        PartialRebinCountsFill(hchargejetpt6080, 10,2);
        PartialRebinCountsFill(hchargejetpt6080_charged, 10,2);
        PartialRebinCountsFill(hchargejetpt6080truth, 10,2);
        PartialRebinCountsFill(hchargejetpt6080truth_charged, 10,2);
        
        PartialRebinCountsFill(hjetpt6080, 30,2);
        PartialRebinCountsFill(hjetpt6080_charged, 30,2);
        PartialRebinCountsFill(hjetpt6080truth, 30,2);
        PartialRebinCountsFill(hjetpt6080truth_charged, 30,2);
        PartialRebinCountsFill(hchargejetpt6080, 30,2);
        PartialRebinCountsFill(hchargejetpt6080_charged, 30,2);
        PartialRebinCountsFill(hchargejetpt6080truth, 30,2);
        PartialRebinCountsFill(hchargejetpt6080truth_charged, 30,2);
        
        PartialRebinCountsFill(hchargept6080, 10,2);
        PartialRebinCountsFill(hchargept6080truth, 10,2);
        PartialRebinCountsFill(hchargept6080, 20,2);
        PartialRebinCountsFill(hchargept6080truth, 20,2);
      }
    }
  }

  MakeTruthHist(hchargept_truth_all,h2chargencollref,hncollall);
  MakeTruthHistFill(hchargept_truth_all_fill,h2chargencollref,hncollall);
  TH1D* hchargeptall = (TH1D*) f1->Get("hchargeptall");
  std::cout << "neventall " << neventall << " vs " << hncollall->Integral() << std::endl;

  hchargeptall->Scale(1./neventall);
  std::cout << "cross check reco vs truth for 0-100% --> bin 10 center reco / truth : " << hchargeptall->GetBinCenter(10) << " / " << hchargept_truth_all->GetBinCenter(10) << " -- content : " << hchargeptall->GetBinContent(10) << " / " << hchargept_truth_all->GetBinContent(10) << " / " << hchargept_truth_all_fill->GetBinContent(10) << std::endl;


  //-------------------------------------------------
  //Make hist for R_CP / I_CP
  //-------------------------------------------------
  TH1D* hRCP_jetptcent =(TH1D*) hjetpt[0]->Clone("hjetptcent");
  TH1D* hRCP_jetptcent_charged =(TH1D*) hjetpt_charged[0]->Clone("hjetptcent_charged");
  TH1D* hRCP_chargeptcent =(TH1D*) hchargept[0]->Clone("hchargeptcent");
  TH1D* hICP_chargejetptcent =(TH1D*) hchargejetpt[0]->Clone("hchargejetptcent");
  TH1D* hICP_chargejetptcent_charged =(TH1D*) hchargejetpt_charged[0]->Clone("hchargejetptcent_charged");
  
  TH1D* hRCP_jetptcenttruth =(TH1D*) hjetpt_truth[0]->Clone("hjetptcent_truth");
  TH1D* hRCP_jetptcenttruth_charged =(TH1D*) hjetpt_truth_charged[0]->Clone("hjetptcent_truth_charged");
  TH1D* hRCP_chargeptcenttruth =(TH1D*) hchargept_truth[0]->Clone("hchargeptcent_truth");
  TH1D* hICP_chargejetptcenttruth =(TH1D*) hchargejetpt_truth[0]->Clone("hchargejetptcent_truth");
  TH1D* hICP_chargejetptcenttruth_charged =(TH1D*) hchargejetpt_truth_charged[0]->Clone("hchargejetptcent_truth_charged");

  //for 0-10 / 60-80% Reco / truth
  CalcRCP(hRCP_jetptcent, hjetpt6080, nevent[0], nevent6080, ncollmean[0], ncoll6080);
  CalcRCP(hRCP_jetptcent_charged, hjetpt6080_charged, nevent[0], nevent6080, ncollmean[0], ncoll6080);
  CalcRCP(hRCP_chargeptcent, hchargept6080, nevent[0], nevent6080, ncollmean[0], ncoll6080);
  CalcICP(hICP_chargejetptcent, hchargejetpt6080, neventcharge[0], neventcharge6080);
  CalcICP(hICP_chargejetptcent_charged, hchargejetpt6080_charged, neventcharge[0], neventcharge6080);
  
  CalcRCP(hRCP_jetptcenttruth, hjetpt6080truth, 1, 1, ncollmean[0], ncoll6080);
  CalcRCP(hRCP_jetptcenttruth_charged, hjetpt6080truth_charged, 1, 1, ncollmean[0], ncoll6080);
  CalcRCP(hRCP_chargeptcenttruth, hchargept6080truth, 1, 1, ncollmean[0], ncoll6080);
  CalcICP(hICP_chargejetptcenttruth, hchargejetpt6080truth, 1, 1);
  CalcICP(hICP_chargejetptcenttruth_charged, hchargejetpt6080truth_charged, 1, 1);

  TwoColorHistRedBlue(hRCP_jetptcent,hRCP_jetptcenttruth);
  TwoColorHistRedBlue(hRCP_jetptcent_charged,hRCP_jetptcenttruth_charged);
  TwoColorHistRedBlue(hRCP_chargeptcent,hRCP_chargeptcenttruth);
  TwoColorHistRedBlue(hICP_chargejetptcent,hICP_chargejetptcenttruth);
  TwoColorHistRedBlue(hICP_chargejetptcent_charged,hICP_chargejetptcenttruth_charged);

  //-------------------------------------------------
  //Make bias factor
  //-------------------------------------------------
  
  TH1D* hbias_jetpt6080 = (TH1D*) hjetpt6080->Clone("hbias_jet6080");
  TH1D* hbias_chargept6080 = (TH1D*) hchargept6080->Clone("hbias_charge6080");
  TH1D* hbias_chargejetpt6080 = (TH1D*) hchargejetpt6080->Clone("hbias_chargejet6080");
  MakeBiasHist(hbias_jetpt6080, nevent6080, hjetpt6080truth);
  MakeBiasHist(hbias_chargept6080, nevent6080, hchargept6080truth);
  MakeBiasHist(hbias_chargejetpt6080, neventcharge6080, hchargejetpt6080truth);
  hbias_jetpt6080->SetLineColor(kBlue+1);
  hbias_chargept6080->SetLineColor(kBlue+1);
  hbias_chargejetpt6080->SetLineColor(kBlue+1);

  TH1D* hbias_jetpt[10];
  TH1D* hbias_jetpt_charged[10];
  TH1D* hbias_chargept[10];
  TH1D* hbias_chargejetpt[10];
  TH1D* hbias_chargejetpt_charged[10];
  for(int i=0;i<10; i++){
    hbias_jetpt[i] = (TH1D*) hjetpt[i]->Clone(Form("hbias_jetpt%d",i));
    hbias_jetpt_charged[i] = (TH1D*) hjetpt_charged[i]->Clone(Form("hbias_jetpt_charged%d",i));
    hbias_chargept[i] = (TH1D*) hchargept[i]->Clone(Form("hbias_chargept%d",i));
    hbias_chargejetpt[i] = (TH1D*) hchargejetpt[i]->Clone(Form("hbias_chargejetpt%d",i));
    hbias_chargejetpt_charged[i] = (TH1D*) hchargejetpt_charged[i]->Clone(Form("hbias_chargejetpt_charged%d",i));

    MakeBiasHist(hbias_jetpt[i], nevent[i], hjetpt_truth[i]);
    MakeBiasHist(hbias_jetpt_charged[i], nevent[i], hjetpt_truth_charged[i]);
    MakeBiasHist(hbias_chargept[i], nevent[i], hchargept_truth[i]);
    MakeBiasHist(hbias_chargejetpt[i], neventcharge[i], hchargejetpt_truth[i]);
    MakeBiasHist(hbias_chargejetpt_charged[i], neventcharge[i], hchargejetpt_truth_charged[i]);
  }



  //-------------------------------------------------
  // Drawing
  //-------------------------------------------------
  
  TCanvas* c_jet_full_spectra= new TCanvas("c_jet_full_spectra","",4,45,550,520);
  c_jet_full_spectra->SetLeftMargin(0.13);
  c_jet_full_spectra->SetRightMargin(0.05);
  c_jet_full_spectra->SetBottomMargin(0.13);
  c_jet_full_spectra->SetTopMargin(0.08);
  c_jet_full_spectra->SetTicks(1,1);
  c_jet_full_spectra->SetLogy();
  c_jet_full_spectra->cd();

  TH1D *hjetfullspectracent = (TH1D*) hjetpt[0]->Clone("hjetfullspectracent");
  TH1D *hjetfullspectraperi = (TH1D*) hjetpt[6]->Clone("hjetfullspectraperi");
  hjetfullspectraperi->Add(hjetpt[7]);
  
  TwoColorHistRedBlue(hjetfullspectracent,hjetfullspectraperi);

  double totalentries = hncollall->Integral();
  TH1ScaleByWidth(hjetfullspectracent);
  TH1ScaleByWidth(hjetfullspectraperi);
  hjetfullspectracent->Scale(1./totalentries);
  hjetfullspectraperi->Scale(1./totalentries);
  hjetfullspectracent->SetFillStyle(4000);
  hjetfullspectracent->GetYaxis()->SetRangeUser(1e-12,hjetfullspectracent->GetMaximum()*5);
  //hjetfullspectracent->GetXaxis()->SetLimits(0, xmax);
  hjetfullspectracent->GetXaxis()->SetRangeUser(0,xmaxjet);
  hjetfullspectracent->GetYaxis()->CenterTitle();
  hjetfullspectracent->SetTitle("");
  hjetfullspectracent->GetYaxis()->SetTitle("dN/dp_{T}");
  hjetfullspectracent->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hjetfullspectracent->GetYaxis()->SetTitleSize(0.048);
  hjetfullspectracent->GetYaxis()->SetLabelSize(0.045);
  hjetfullspectracent->GetXaxis()->SetLabelSize(0.045);
  hjetfullspectracent->GetYaxis()->SetTitleOffset(1.2) ;
  hjetfullspectracent->GetXaxis()->SetTitleOffset(1.0) ;
  hjetfullspectracent->GetXaxis()->SetTitleSize(0.05) ;
  hjetfullspectracent->GetYaxis()->SetMaxDigits(3);
  hjetfullspectracent->GetYaxis()->SetNoExponent(false);
  hjetfullspectracent->GetXaxis()->CenterTitle();
  hjetfullspectracent->Draw();
  hjetfullspectraperi->Draw("same");

  TLegend* l1= new TLegend(0.67,0.47,0.84,0.603);
  SetLegendStyle(l1);
  l1->SetTextFont(43);
  l1->SetTextSize(17);
  l1->SetBorderSize(0);
  l1->AddEntry(hjetfullspectracent,"0-10%","l");
  l1->AddEntry(hjetfullspectraperi,"60-80%","l");
  l1->Draw("same");
  
  drawText("Anti-k_{T} R=0.2, Full jet",xpos,ypos-ydiff,1,16);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacut),xpos,ypos-2*ydiff,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-3*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_jet_full_spectra->SaveAs(Form("plots/c_spectra_jet_full_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  TCanvas* c_jet_charged_spectra= new TCanvas("c_jet_charged_spectra","",4,45,550,520);
  c_jet_charged_spectra->SetLeftMargin(0.13);
  c_jet_charged_spectra->SetRightMargin(0.05);
  c_jet_charged_spectra->SetBottomMargin(0.13);
  c_jet_charged_spectra->SetTopMargin(0.08);
  c_jet_charged_spectra->SetTicks(1,1);
  c_jet_charged_spectra->SetLogy();
  c_jet_charged_spectra->cd();

  TH1D *hjetchargedspectracent = (TH1D*) hjetpt_charged[0]->Clone("hjetchargedspectracent");
  TH1D *hjetchargedspectraperi = (TH1D*) hjetpt_charged[6]->Clone("hjetchargedspectraperi");
  hjetchargedspectraperi->Add(hjetpt_charged[7]);
  TwoColorHistRedBlue(hjetchargedspectracent,hjetchargedspectraperi);

  TH1ScaleByWidth(hjetchargedspectracent);
  TH1ScaleByWidth(hjetchargedspectraperi);
  hjetchargedspectracent->Scale(1./totalentries);
  hjetchargedspectraperi->Scale(1./totalentries);
  hjetchargedspectracent->SetFillStyle(4000);
  hjetchargedspectracent->GetYaxis()->SetRangeUser(1e-12,hjetchargedspectracent->GetMaximum()*5);
  //hjetchargedspectracent->GetXaxis()->SetLimits(0, xmax);
  hjetchargedspectracent->GetXaxis()->SetRangeUser(0,xmaxjet);
  hjetchargedspectracent->GetYaxis()->CenterTitle();
  hjetchargedspectracent->SetTitle("");
  hjetchargedspectracent->GetYaxis()->SetTitle("dN/dp_{T}");
  hjetchargedspectracent->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hjetchargedspectracent->GetYaxis()->SetTitleSize(0.048);
  hjetchargedspectracent->GetYaxis()->SetLabelSize(0.045);
  hjetchargedspectracent->GetXaxis()->SetLabelSize(0.045);
  hjetchargedspectracent->GetYaxis()->SetTitleOffset(1.2) ;
  hjetchargedspectracent->GetXaxis()->SetTitleOffset(1.0) ;
  hjetchargedspectracent->GetXaxis()->SetTitleSize(0.05) ;
  hjetchargedspectracent->GetYaxis()->SetMaxDigits(3);
  hjetchargedspectracent->GetYaxis()->SetNoExponent(false);
  hjetchargedspectracent->GetXaxis()->CenterTitle();
  hjetchargedspectracent->Draw();
  hjetchargedspectraperi->Draw("same");
  l1->Draw("same");
  
  drawText("Anti-k_{T} R=0.2, Charged jet",xpos,ypos-ydiff,1,16);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacut),xpos,ypos-2*ydiff,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-3*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_jet_charged_spectra->SaveAs(Form("plots/c_spectra_jet_charged_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));
  
  TCanvas* c_spectra_charged= new TCanvas("c_spectra_charged","",4,45,550,520);
  c_spectra_charged->SetLeftMargin(0.13);
  c_spectra_charged->SetRightMargin(0.05);
  c_spectra_charged->SetBottomMargin(0.13);
  c_spectra_charged->SetTopMargin(0.08);
  c_spectra_charged->SetTicks(1,1);
  c_spectra_charged->SetLogy();
  c_spectra_charged->cd();

  TH1D *hspectrachargecent = (TH1D*) hchargept[0]->Clone("hspectrachargecent");
  TH1D *hspectrachargeperi = (TH1D*) hchargept[6]->Clone("hspectrachargeperi");
  hspectrachargeperi->Add(hchargept[7]);

  TwoColorHistRedBlue(hspectrachargecent,hspectrachargeperi);

  TH1ScaleByWidth(hspectrachargecent);
  TH1ScaleByWidth(hspectrachargeperi);
  hspectrachargecent->SetLineWidth(2);
  hspectrachargeperi->SetLineWidth(2);
  hspectrachargecent->SetLineColor(kRed+1);
  hspectrachargeperi->SetLineColor(kBlue+1);
  hspectrachargecent->Scale(1./totalentries);
  hspectrachargeperi->Scale(1./totalentries);
  hspectrachargecent->SetFillStyle(4000);
  hspectrachargecent->GetYaxis()->SetRangeUser(1e-12,hspectrachargecent->GetMaximum()*5);
  //hspectrachargecent->GetXaxis()->SetLimits(0, xmax);
  hspectrachargecent->GetXaxis()->SetRangeUser(0,xmaxcharge);
  hspectrachargecent->GetYaxis()->CenterTitle();
  hspectrachargecent->SetTitle("");
  hspectrachargecent->GetYaxis()->SetTitle("dN/dp_{T}");
  hspectrachargecent->GetXaxis()->SetTitle("p_{T}^{h^{#pm}} [GeV]");
  hspectrachargecent->GetYaxis()->SetTitleSize(0.048);
  hspectrachargecent->GetYaxis()->SetLabelSize(0.045);
  hspectrachargecent->GetXaxis()->SetLabelSize(0.045);
  hspectrachargecent->GetYaxis()->SetTitleOffset(1.2) ;
  hspectrachargecent->GetXaxis()->SetTitleOffset(1.0) ;
  hspectrachargecent->GetXaxis()->SetTitleSize(0.05) ;
  hspectrachargecent->GetYaxis()->SetMaxDigits(3);
  hspectrachargecent->GetYaxis()->SetNoExponent(false);
  hspectrachargecent->GetXaxis()->CenterTitle();
  hspectrachargecent->Draw();
  hspectrachargeperi->Draw("same");

  l1->Draw("same");
  
  drawText(Form("|#eta^{h^{#pm}}| < %.1f",hadronetacut),xpos,ypos-ydiff,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-2*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_spectra_charged->SaveAs(Form("plots/c_spectra_charged_hadrons_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  //R_CP
  TCanvas* c_RCP_jet_full= new TCanvas("c_RCP_jet_full","",4,45,550,520);
  c_RCP_jet_full->SetLeftMargin(0.13);
  c_RCP_jet_full->SetRightMargin(0.05);
  c_RCP_jet_full->SetBottomMargin(0.13);
  c_RCP_jet_full->SetTopMargin(0.08);
  c_RCP_jet_full->SetTicks(1,1);
  c_RCP_jet_full->cd();

  hRCP_jetptcent->SetFillStyle(4000);
  //hRCP_jetptcent->GetXaxis()->SetLimits(0, xmax);
  hRCP_jetptcent->GetXaxis()->SetRangeUser(0,xmaxjet);
  hRCP_jetptcent->GetYaxis()->SetRangeUser(0,4);
  hRCP_jetptcent->GetYaxis()->CenterTitle();
  hRCP_jetptcent->SetTitle("");
  hRCP_jetptcent->GetYaxis()->SetTitle("R_{CP}");
  hRCP_jetptcent->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hRCP_jetptcent->GetYaxis()->SetTitleSize(0.048);
  hRCP_jetptcent->GetYaxis()->SetLabelSize(0.045);
  hRCP_jetptcent->GetXaxis()->SetLabelSize(0.045);
  hRCP_jetptcent->GetYaxis()->SetTitleOffset(1.2) ;
  hRCP_jetptcent->GetXaxis()->SetTitleOffset(1.0) ;
  hRCP_jetptcent->GetXaxis()->SetTitleSize(0.05) ;
  hRCP_jetptcent->GetYaxis()->SetMaxDigits(3);
  hRCP_jetptcent->GetYaxis()->SetNoExponent(false);
  hRCP_jetptcent->GetXaxis()->CenterTitle();
  hRCP_jetptcent->Draw();
  hRCP_jetptcenttruth->Draw("same");
  
  TLegend* l_RCP_fulljet= new TLegend(0.31,0.63,0.48,0.753);
  SetLegendStyle(l_RCP_fulljet);
  l_RCP_fulljet->SetTextFont(43);
  l_RCP_fulljet->SetTextSize(17);
  l_RCP_fulljet->SetBorderSize(0);
  l_RCP_fulljet->AddEntry(hRCP_jetptcent,"Reco","l");
  l_RCP_fulljet->AddEntry(hRCP_jetptcenttruth,"Truth","l");
  l_RCP_fulljet->Draw("same");
  dashedLine(0,1,xmaxjet,1,1,1);
  
  drawText("0-10%/60-80%",xpos,ypos-ydiff,1,16);
  drawText("Anti-k_{T} R=0.2, Full jet",xpos,ypos-2*ydiff,1,16);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacut),xpos,ypos-3*ydiff,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-4*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_RCP_jet_full->SaveAs(Form("plots/c_RCP_jet_full_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  TCanvas* c_RCP_jet_charged= new TCanvas("c_RCP_jet_charged","",4,45,550,520);
  c_RCP_jet_charged->SetLeftMargin(0.13);
  c_RCP_jet_charged->SetRightMargin(0.05);
  c_RCP_jet_charged->SetBottomMargin(0.13);
  c_RCP_jet_charged->SetTopMargin(0.08);
  c_RCP_jet_charged->SetTicks(1,1);
  c_RCP_jet_charged->cd();

  hRCP_jetptcent->SetFillStyle(4000);
  //hRCP_jetptcent->GetXaxis()->SetLimits(0, xmax);
  hRCP_jetptcent_charged->GetXaxis()->SetRangeUser(0,xmaxjet);
  hRCP_jetptcent_charged->GetYaxis()->SetRangeUser(0,4);
  hRCP_jetptcent_charged->GetYaxis()->CenterTitle();
  hRCP_jetptcent_charged->SetTitle("");
  hRCP_jetptcent_charged->GetYaxis()->SetTitle("R_{CP}");
  hRCP_jetptcent_charged->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hRCP_jetptcent_charged->GetYaxis()->SetTitleSize(0.048);
  hRCP_jetptcent_charged->GetYaxis()->SetLabelSize(0.045);
  hRCP_jetptcent_charged->GetXaxis()->SetLabelSize(0.045);
  hRCP_jetptcent_charged->GetYaxis()->SetTitleOffset(1.2) ;
  hRCP_jetptcent_charged->GetXaxis()->SetTitleOffset(1.0) ;
  hRCP_jetptcent_charged->GetXaxis()->SetTitleSize(0.05) ;
  hRCP_jetptcent_charged->GetYaxis()->SetMaxDigits(3);
  hRCP_jetptcent_charged->GetYaxis()->SetNoExponent(false);
  hRCP_jetptcent_charged->GetXaxis()->CenterTitle();
  hRCP_jetptcent_charged->Draw();
  hRCP_jetptcenttruth_charged->Draw("same");
  
  l_RCP_fulljet->Draw("same");
  dashedLine(0,1,xmaxjet,1,1,1);
  
  drawText("0-10%/60-80%",xpos,ypos-ydiff,1,16);
  drawText("Anti-k_{T} R=0.2, Charged jet",xpos,ypos-2*ydiff,1,16);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacut),xpos,ypos-3*ydiff,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-4*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_RCP_jet_charged->SaveAs(Form("plots/c_RCP_jet_charged_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  //Charged R_CP
  TCanvas* c_RCP_charged= new TCanvas("c_RCP_charged","",4,45,550,520);
  c_RCP_charged->SetLeftMargin(0.13);
  c_RCP_charged->SetRightMargin(0.05);
  c_RCP_charged->SetBottomMargin(0.13);
  c_RCP_charged->SetTopMargin(0.08);
  c_RCP_charged->SetTicks(1,1);
  c_RCP_charged->cd();

  hRCP_chargeptcent->SetFillStyle(4000);
  //hRCP_chargeptcent->GetXaxis()->SetLimits(0, xmax);
  hRCP_chargeptcent->GetXaxis()->SetRangeUser(0,xmaxcharge);
  hRCP_chargeptcent->GetYaxis()->SetRangeUser(0,4);
  hRCP_chargeptcent->GetYaxis()->CenterTitle();
  hRCP_chargeptcent->SetTitle("");
  hRCP_chargeptcent->GetYaxis()->SetTitle("R_{CP}");
  hRCP_chargeptcent->GetXaxis()->SetTitle("p_{T}^{h^{#pm}} [GeV]");
  hRCP_chargeptcent->GetYaxis()->SetTitleSize(0.048);
  hRCP_chargeptcent->GetYaxis()->SetLabelSize(0.045);
  hRCP_chargeptcent->GetXaxis()->SetLabelSize(0.045);
  hRCP_chargeptcent->GetYaxis()->SetTitleOffset(1.2) ;
  hRCP_chargeptcent->GetXaxis()->SetTitleOffset(1.0) ;
  hRCP_chargeptcent->GetXaxis()->SetTitleSize(0.05) ;
  hRCP_chargeptcent->GetYaxis()->SetMaxDigits(3);
  hRCP_chargeptcent->GetYaxis()->SetNoExponent(false);
  hRCP_chargeptcent->GetXaxis()->CenterTitle();
  hRCP_chargeptcent->Draw();
  hRCP_chargeptcenttruth->SetLineColor(kBlue);
  hRCP_chargeptcenttruth->SetLineWidth(2);
  hRCP_chargeptcenttruth->Draw("same");
  
  TLegend* l_RCP_charge= new TLegend(0.31,0.63,0.48,0.753);
  SetLegendStyle(l_RCP_charge);
  l_RCP_charge->SetTextFont(43);
  l_RCP_charge->SetTextSize(17);
  l_RCP_charge->SetBorderSize(0);
  l_RCP_charge->AddEntry(hRCP_chargeptcent,"Reco","l");
  l_RCP_charge->AddEntry(hRCP_chargeptcenttruth,"Truth","l");
  l_RCP_charge->Draw("same");
  dashedLine(0,1,xmaxcharge,1,1,1);
  
  //drawText("Anti-k_{T} R=0.2, Full jet",xpos,ypos-ydiff,1,16);
  drawText("0-10%/60-80%",xpos,ypos-ydiff,1,16);
  drawText(Form("|#eta^{h^{#pm}}| < %.1f",hadronetacut),xpos,ypos-2*ydiff,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-3*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_RCP_charged->SaveAs(Form("plots/c_RCP_charged_hadrons_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  //I_CP
  TCanvas* c_ICP_jet_full= new TCanvas("c_ICP_jet_full","",4,45,550,520);
  c_ICP_jet_full->SetLeftMargin(0.13);
  c_ICP_jet_full->SetRightMargin(0.05);
  c_ICP_jet_full->SetBottomMargin(0.13);
  c_ICP_jet_full->SetTopMargin(0.08);
  c_ICP_jet_full->SetTicks(1,1);
  c_ICP_jet_full->cd();

  hICP_chargejetptcent->SetFillStyle(4000);
  //hICP_chargejetptcent->GetXaxis()->SetLimits(0, xmax);
  hICP_chargejetptcent->GetXaxis()->SetRangeUser(0,xmaxjet);
  hICP_chargejetptcent->GetYaxis()->SetRangeUser(0,4);
  hICP_chargejetptcent->GetYaxis()->CenterTitle();
  hICP_chargejetptcent->SetTitle("");
  hICP_chargejetptcent->GetYaxis()->SetTitle("I_{CP}");
  hICP_chargejetptcent->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hICP_chargejetptcent->GetYaxis()->SetTitleSize(0.048);
  hICP_chargejetptcent->GetYaxis()->SetLabelSize(0.045);
  hICP_chargejetptcent->GetXaxis()->SetLabelSize(0.045);
  hICP_chargejetptcent->GetYaxis()->SetTitleOffset(1.2) ;
  hICP_chargejetptcent->GetXaxis()->SetTitleOffset(1.0) ;
  hICP_chargejetptcent->GetXaxis()->SetTitleSize(0.05) ;
  hICP_chargejetptcent->GetYaxis()->SetMaxDigits(3);
  hICP_chargejetptcent->GetYaxis()->SetNoExponent(false);
  hICP_chargejetptcent->GetXaxis()->CenterTitle();
  hICP_chargejetptcent->Draw();
  hICP_chargejetptcenttruth->Draw("same");
  
  TLegend* l_ICP_full = new TLegend(0.28,0.65,0.45,0.783);
  SetLegendStyle(l_ICP_full);
  l_ICP_full->SetTextFont(43);
  l_ICP_full->SetTextSize(15);
  l_ICP_full->SetBorderSize(0);
  l_ICP_full->AddEntry(hICP_chargejetptcent,"Reco","l");
  l_ICP_full->AddEntry(hICP_chargejetptcenttruth,"Truth","l");
  l_ICP_full->Draw("same");
  dashedLine(0,1,xmaxjet,1,1,1);
  
  double xposicp = 0.41; double yposicp = 0.78;
  double ydifficp = 0.055;
  drawText("0-10%/60-80%",xposicp,yposicp,1,15);
  drawText("Anti-k_{T} R=0.2, Full jet",xposicp,yposicp-ydifficp,1,15);
  drawText(Form("|#eta^{h^{#pm}}| < %.1f, %.f < p_{T}^{h^{#pm}} < %.f",hadronetacutForICP,hadronptlowForICP,hadronpthighForICP),xposicp,yposicp-2*ydifficp,1,15);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacutForICP),xposicp,yposicp-3*ydifficp,1,15);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xposicp,yposicp-4*ydifficp,1,15);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_ICP_jet_full->SaveAs(Form("plots/c_ICP_full_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));
  
  
  TCanvas* c_ICP_jet_charged= new TCanvas("c_ICP_jet_charged","",4,45,550,520);
  c_ICP_jet_charged->SetLeftMargin(0.13);
  c_ICP_jet_charged->SetRightMargin(0.05);
  c_ICP_jet_charged->SetBottomMargin(0.13);
  c_ICP_jet_charged->SetTopMargin(0.08);
  c_ICP_jet_charged->SetTicks(1,1);
  c_ICP_jet_charged->cd();

  hICP_chargejetptcent_charged->SetFillStyle(4000);
  hICP_chargejetptcent_charged->GetXaxis()->SetRangeUser(0,xmaxjet);
  hICP_chargejetptcent_charged->GetYaxis()->SetRangeUser(0,4);
  hICP_chargejetptcent_charged->GetYaxis()->CenterTitle();
  hICP_chargejetptcent_charged->SetTitle("");
  hICP_chargejetptcent_charged->GetYaxis()->SetTitle("I_{CP}");
  hICP_chargejetptcent_charged->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hICP_chargejetptcent_charged->GetYaxis()->SetTitleSize(0.048);
  hICP_chargejetptcent_charged->GetYaxis()->SetLabelSize(0.045);
  hICP_chargejetptcent_charged->GetXaxis()->SetLabelSize(0.045);
  hICP_chargejetptcent_charged->GetYaxis()->SetTitleOffset(1.2) ;
  hICP_chargejetptcent_charged->GetXaxis()->SetTitleOffset(1.0) ;
  hICP_chargejetptcent_charged->GetXaxis()->SetTitleSize(0.05) ;
  hICP_chargejetptcent_charged->GetYaxis()->SetMaxDigits(3);
  hICP_chargejetptcent_charged->GetYaxis()->SetNoExponent(false);
  hICP_chargejetptcent_charged->GetXaxis()->CenterTitle();
  hICP_chargejetptcent_charged->Draw();
  hICP_chargejetptcenttruth_charged->Draw("same");
  
  l_ICP_full->Draw("same");
  dashedLine(0,1,xmaxjet,1,1,1);
  
  drawText("0-10%/60-80%",xposicp,yposicp,1,15);
  drawText("Anti-k_{T} R=0.2, Charged jet",xposicp,yposicp-ydifficp,1,15);
  drawText(Form("|#eta^{h^{#pm}}| < %.1f, %.f < p_{T}^{h^{#pm}} < %.f",hadronetacutForICP,hadronptlowForICP,hadronpthighForICP),xposicp,yposicp-2*ydifficp,1,15);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacutForICP),xposicp,yposicp-3*ydifficp,1,15);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xposicp,yposicp-4*ydifficp,1,15);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_ICP_jet_charged->SaveAs(Form("plots/c_ICP_charged_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));


  //Cent, Ncoll
  TCanvas* c_cent= new TCanvas("c_cent","",4,45,550,520);
  c_cent->SetLeftMargin(0.13);
  c_cent->SetRightMargin(0.05);
  c_cent->SetBottomMargin(0.13);
  c_cent->SetTopMargin(0.08);
  c_cent->SetTicks(1,1);
  c_cent->cd();

  hcent->Scale(1./hcent->GetEntries()*100);
  hcent->SetLineColor(kBlack);
  hcent->SetLineWidth(2);
  hcent->SetFillStyle(4000);
  hcent->GetYaxis()->CenterTitle();
  hcent->SetTitle("");
  hcent->GetXaxis()->SetRangeUser(0,100);
  hcent->GetYaxis()->SetRangeUser(0,2);
  hcent->GetYaxis()->SetTitle("");
  hcent->GetXaxis()->SetTitle("Centrality (%)");
  hcent->GetYaxis()->SetTitleSize(0.048);
  hcent->GetYaxis()->SetLabelSize(0.045);
  hcent->GetXaxis()->SetLabelSize(0.045);
  hcent->GetYaxis()->SetTitleOffset(1.2) ;
  hcent->GetXaxis()->SetTitleOffset(1.0) ;
  hcent->GetXaxis()->SetTitleSize(0.05) ;
  hcent->GetYaxis()->SetMaxDigits(3);
  hcent->GetYaxis()->SetNoExponent(false);
  hcent->GetXaxis()->CenterTitle();
  hcent->Draw();

  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-3*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_cent->SaveAs(Form("plots/c_centrality_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  TCanvas* c_ncollall= new TCanvas("c_ncollall","",4,45,550,520);
  c_ncollall->SetLeftMargin(0.13);
  c_ncollall->SetRightMargin(0.05);
  c_ncollall->SetBottomMargin(0.13);
  c_ncollall->SetTopMargin(0.08);
  c_ncollall->SetTicks(1,1);
  c_ncollall->SetLogy();
  c_ncollall->cd();

  hncollall->SetFillStyle(4000);
  hncollall->GetYaxis()->CenterTitle();
  hncollall->SetTitle("");
  hncollall->GetYaxis()->SetTitle("");
  hncollall->GetXaxis()->SetTitle("N_{coll}");
  hncollall->GetYaxis()->SetTitleSize(0.048);
  hncollall->GetYaxis()->SetLabelSize(0.045);
  hncollall->GetXaxis()->SetLabelSize(0.045);
  hncollall->GetYaxis()->SetTitleOffset(1.2) ;
  hncollall->GetXaxis()->SetTitleOffset(1.0) ;
  hncollall->GetXaxis()->SetTitleSize(0.05) ;
  hncollall->GetYaxis()->SetMaxDigits(3);
  hncollall->GetYaxis()->SetNoExponent(false);
  hncollall->GetXaxis()->CenterTitle();
  hncollall->Draw();

  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xpos,ypos-3*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  c_ncollall->SaveAs(Form("plots/c_ncollall_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  TCanvas* c_ncolldiv= new TCanvas("c_ncolldiv","",4,45,580,520);
  c_ncolldiv->SetLeftMargin(0.11);
  c_ncolldiv->SetRightMargin(0.20);
  c_ncolldiv->SetBottomMargin(0.13);
  c_ncolldiv->SetTopMargin(0.08);
  c_ncolldiv->SetTicks(1,1);
  c_ncolldiv->SetLogy();
  c_ncolldiv->cd();

  TLegend* lncoll= new TLegend(0.82,0.35,1.10,0.793);
  SetLegendStyle(lncoll);
  lncoll->SetTextFont(43);
  lncoll->SetTextSize(17);
  lncoll->SetBorderSize(0);
  
  std::string perc="%";
  double max = hncoll[0]->GetMaximum();
  hncoll[0]->GetYaxis()->SetRangeUser(1e-1,max*100);
  const Int_t ncolors = 10;
  Int_t colorIdx[ncolors];

  const int stops = 7;
  double stops_arr[stops] = {0.0, 0.17, 0.33, 0.5, 0.67, 0.83, 1.0};
  double red[stops]   = {1.0, 1.0, 1.0, 0.2, 0.0, 0.5, 1.0};
  double green[stops] = {0.0, 0.5, 1.0, 1.0, 0.6, 0.0, 0.6};
  double blue[stops]  = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  Int_t startColor = TColor::CreateGradientColorTable(stops, stops_arr, red, green, blue, ncolors);
  for (int i = 0; i < ncolors; ++i)
    colorIdx[i] = startColor + i;

  for(int i=0;i<10; i++){
    hncoll[i]->SetLineColor(colorIdx[i]);
    hncoll[i]->SetLineWidth(2);
    hncoll[i]->SetFillStyle(4000);
    hncoll[i]->GetYaxis()->CenterTitle();
    hncoll[i]->SetTitle("");
    hncoll[i]->GetYaxis()->SetTitle("");
    hncoll[i]->GetXaxis()->SetTitle("N_{coll}");
    hncoll[i]->GetYaxis()->SetTitleSize(0.048);
    hncoll[i]->GetYaxis()->SetLabelSize(0.045);
    hncoll[i]->GetXaxis()->SetLabelSize(0.045);
    hncoll[i]->GetYaxis()->SetTitleOffset(1.2) ;
    hncoll[i]->GetXaxis()->SetTitleOffset(1.0) ;
    hncoll[i]->GetXaxis()->SetTitleSize(0.05) ;
    hncoll[i]->GetYaxis()->SetMaxDigits(3);
    hncoll[i]->GetYaxis()->SetNoExponent(false);
    hncoll[i]->GetXaxis()->CenterTitle();
    hncoll[i]->Draw("hist same");
    lncoll->AddEntry(hncoll[i],Form("%d-%d%s",i*10,(i+1)*10,perc.c_str()),"l");
    
    //set color also for bias
    hbias_jetpt[i]->SetLineWidth(2);
    hbias_jetpt[i]->SetLineColor(colorIdx[i]);
    hbias_jetpt[i]->SetMarkerColor(colorIdx[i]);
    hbias_jetpt_charged[i]->SetLineWidth(2);
    hbias_jetpt_charged[i]->SetLineColor(colorIdx[i]);
    hbias_jetpt_charged[i]->SetMarkerColor(colorIdx[i]);
    hbias_chargept[i]->SetLineWidth(2);
    hbias_chargept[i]->SetLineColor(colorIdx[i]);
    hbias_chargept[i]->SetMarkerColor(colorIdx[i]);
    hbias_chargejetpt[i]->SetLineWidth(2);
    hbias_chargejetpt[i]->SetLineColor(colorIdx[i]);
    hbias_chargejetpt_charged[i]->SetLineColor(colorIdx[i]);
  }

  lncoll->Draw("same");
  drawText(Form("Centrality calib. w/ %s",cent_type_txt.c_str()),xpos+0.03,ypos-1.7*ydiff,1,16);
  drawText(Form("(%.1f < |#eta| < %.1f)",etalow_cent,etahigh_cent),xpos+0.03,ypos-2.8*ydiff,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.26,0.85,1,18);
  c_ncolldiv->SaveAs(Form("plots/c_ncoll_div_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  //Bias factor jet 
  TH1D* hdummy = new TH1D("hdummy",";;",300,0,300);
  hdummy->SetFillStyle(4000);
  hdummy->GetYaxis()->CenterTitle();
  hdummy->SetTitle("");
  hdummy->GetYaxis()->SetTitle("Bias factor");
  hdummy->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hdummy->GetYaxis()->SetTitleSize(0.048);
  hdummy->GetYaxis()->SetLabelSize(0.045);
  hdummy->GetXaxis()->SetLabelSize(0.045);
  hdummy->GetYaxis()->SetTitleOffset(1.2) ;
  hdummy->GetXaxis()->SetTitleOffset(1.0) ;
  hdummy->GetXaxis()->SetTitleSize(0.05) ;
  hdummy->GetYaxis()->SetRangeUser(0,2);
  hdummy->GetXaxis()->SetRangeUser(0,xmaxjet);

  TCanvas *cbias_jet_full = new TCanvas("cbias_jet_full","",700,600);
  cbias_jet_full->SetLeftMargin(0.13);
  cbias_jet_full->SetRightMargin(0.20);
  cbias_jet_full->SetBottomMargin(0.13);
  cbias_jet_full->SetTopMargin(0.08);
  cbias_jet_full->SetTicks(1,1);
  cbias_jet_full->cd();
  hdummy->Draw();

  TLegend* lbias_jet= new TLegend(0.83,0.35,1.13,0.793);
  SetLegendStyle(lbias_jet);
  lbias_jet->SetTextFont(43);
  lbias_jet->SetTextSize(17);
  lbias_jet->SetBorderSize(0);
  
  for(int i=0; i<10; i++){
    hbias_jetpt[i]->Draw("same");
    lbias_jet->AddEntry(hbias_jetpt[i],Form("%d-%d%s",i*10,(i+1)*10,perc.c_str()),"l");
  }

  double xposbias = 0.22; double yposbias=0.88; double ydiffbias = 0.05;
  drawText("Anti-k_{T} R=0.2, Full jet",xposbias,yposbias-2*ydiffbias,1,16);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacut),xposbias,yposbias-3*ydiffbias,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xposbias,yposbias-4*ydiffbias,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.30,0.86,1,18);
  dashedLine(0,1,xmaxjet,1,1,1);
  lbias_jet->Draw("same");
  cbias_jet_full->SaveAs(Form("plots/c_bias_jet_full_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  TCanvas *cbias_jet_charged = new TCanvas("cbias_jet_charged","",700,600);
  cbias_jet_charged->SetLeftMargin(0.13);
  cbias_jet_charged->SetRightMargin(0.20);
  cbias_jet_charged->SetBottomMargin(0.13);
  cbias_jet_charged->SetTopMargin(0.08);
  cbias_jet_charged->SetTicks(1,1);
  cbias_jet_charged->cd();
  hdummy->Draw();
  for(int i=0; i<10; i++)
    hbias_jetpt_charged[i]->Draw("same");

  drawText("Anti-k_{T} R=0.2, Charged jet",xposbias,yposbias-2*ydiffbias,1,16);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacut),xposbias,yposbias-3*ydiffbias,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xposbias,yposbias-4*ydiffbias,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.30,0.86,1,18);
  dashedLine(0,1,xmaxjet,1,1,1);
  lbias_jet->Draw("same");
  cbias_jet_charged->SaveAs(Form("plots/c_bias_jet_charged_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));
  
  
  //Bias factor charged 
  TCanvas *cbias_charged = new TCanvas("cbias_charged","",700,600);
  cbias_charged->SetLeftMargin(0.13);
  cbias_charged->SetRightMargin(0.20);
  cbias_charged->SetBottomMargin(0.13);
  cbias_charged->SetTopMargin(0.08);
  cbias_charged->SetTicks(1,1);
  cbias_charged->cd();

  hdummy->GetXaxis()->SetRangeUser(0,xmaxcharge);
  hdummy->GetXaxis()->SetTitle("p_{T}^{charged} [GeV]");
  hdummy->Draw();

  for(int i=0; i<10; i++)
    hbias_chargept[i]->Draw("same");

  drawText(Form("|#eta^{h^{#pm}}| < %.1f",hadronetacut),xposbias,yposbias-2*ydiffbias,1,16);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xposbias,yposbias-3*ydiffbias,1,16);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.30,0.86,1,18);
  dashedLine(0,1,xmaxcharge,1,1,1);
  lbias_jet->Draw("same");
  cbias_charged->SaveAs(Form("plots/c_bias_charged_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));
  
  //Bias factor I_CP
  xposicp = 0.2; yposicp=0.45;
  TCanvas *cbias_ICP_jet_full = new TCanvas("cbias_ICP_jet_full","",700,600);
  cbias_ICP_jet_full->SetLeftMargin(0.13);
  cbias_ICP_jet_full->SetRightMargin(0.20);
  cbias_ICP_jet_full->SetBottomMargin(0.13);
  cbias_ICP_jet_full->SetTopMargin(0.08);
  cbias_ICP_jet_full->SetTicks(1,1);
  cbias_ICP_jet_full->cd();

  hdummy->GetXaxis()->SetRangeUser(0,xmaxjet);
  hdummy->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hdummy->Draw();

  for(int i=0; i<10; i++)
    hbias_chargejetpt[i]->Draw("same");
  
  drawText("Anti-k_{T} R=0.2, Full jet",xposicp,yposicp-ydifficp,1,15);
  drawText(Form("|#eta^{h^{#pm}}| < %.1f, %.f < p_{T}^{h^{#pm}} < %.f",hadronetacutForICP,hadronptlowForICP,hadronpthighForICP),xposicp,yposicp-2*ydifficp,1,15);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacutForICP),xposicp,yposicp-3*ydifficp,1,15);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xposicp,yposicp-4*ydifficp,1,15);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  dashedLine(0,1,xmaxjet,1,1,1);
  lbias_jet->Draw("same");
  cbias_ICP_jet_full->SaveAs(Form("plots/c_bias_ICP_fulljet_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));

  TCanvas *cbias_ICP_jet_charged = new TCanvas("cbias_ICP_jet_charged","",700,600);
  cbias_ICP_jet_charged->SetLeftMargin(0.13);
  cbias_ICP_jet_charged->SetRightMargin(0.20);
  cbias_ICP_jet_charged->SetBottomMargin(0.13);
  cbias_ICP_jet_charged->SetTopMargin(0.08);
  cbias_ICP_jet_charged->SetTicks(1,1);
  cbias_ICP_jet_charged->cd();
  hdummy->Draw();

  for(int i=0; i<10; i++)
    hbias_chargejetpt_charged[i]->Draw("same");
  
  drawText("Anti-k_{T} R=0.2, Charged jet",xposicp,yposicp-ydifficp,1,15);
  drawText(Form("|#eta^{h^{#pm}}| < %.1f, %.f < p_{T}^{h^{#pm}} < %.f",hadronetacutForICP,hadronptlowForICP,hadronpthighForICP),xposicp,yposicp-2*ydifficp,1,15);
  drawText(Form("|#eta^{jet}| < %.1f",jetetacutForICP),xposicp,yposicp-3*ydifficp,1,15);
  drawText(Form("Centrality calib. w/ %s (%.1f < |#eta| < %.1f)",cent_type_txt.c_str(),etalow_cent,etahigh_cent),xposicp,yposicp-4*ydifficp,1,15);
  drawText(Form("%s O+O  #sqrt{s_{NN}} = %s",simtxt.c_str(),colltxt.c_str()),0.38,0.85,1,18);
  dashedLine(0,1,xmaxjet,1,1,1);
  lbias_jet->Draw("same");
  cbias_ICP_jet_charged->SaveAs(Form("plots/c_bias_ICP_chargedjet_%s_%s_%s.pdf",simtype.c_str(),collen.c_str(),cent_type.c_str()));
} 
