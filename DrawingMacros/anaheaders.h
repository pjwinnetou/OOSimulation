#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include "commonUtility.h"
#include "Style_jaebeom.h"
using namespace std;
double nevent[10];
double neventall=0;
double neventcharge[10];
double ncollmean[10];

std::string formatNumber(double num) {
  if (num >= 1e9)      return Form("%.1fB", num / 1e9);  
  else if (num >= 1e6) return Form("%.1fM", num / 1e6); 
  else if (num >= 1e3) return Form("%.1fk", num / 1e3); 
  else                 return Form("%.1f", num);       
}

void NormBias2D(TH2D* h2d, TH1D* hncoll){
  int nbinsx = h2d->GetNbinsX();
  int nbinsy = h2d->GetNbinsY();

  h2d->Sumw2();
  for(int i=1; i<=nbinsx;++i){
    double ncollval = h2d->GetXaxis()->GetBinCenter(i);
    int ncollbin = hncoll->FindFixBin(ncollval);
    double ncollentries = hncoll->GetBinContent(ncollbin);
    if(ncollentries<=0){
      continue;
    }
    for(int j=1; j<=nbinsy;++j){
      double con = h2d->GetBinContent(i,j);
      con = con/ncollentries;
      double err = h2d->GetBinError(i,j) / ncollentries;
      h2d->SetBinContent(i,j, con);
      h2d->SetBinError(i,j, err);
    }
  }
}

void MakeTruthHist(TH1D* hist, TH2D* h2d, TH1D* hncollinput){
  int nbinsx = h2d->GetNbinsX();
  int nbinsy = h2d->GetNbinsY();

  hist->Reset();
  hist->Sumw2();

  TH1D* hncoll_norm = (TH1D*) hncollinput->Clone("hncoll_norm");
  hncoll_norm->Scale(1./hncoll_norm->Integral()); // make as probability hist.
  std::cout << "hncoll norm integral " << hncoll_norm->Integral() << std::endl;

  for(int i=1; i<=nbinsx;++i){
    double prob = hncoll_norm->GetBinContent(i); 
    for(int j=1; j<=nbinsy;++j){
      double ybincenter = h2d->GetYaxis()->GetBinCenter(j);
      double con = h2d->GetBinContent(i,j) * prob;
      double err = h2d->GetBinError(i,j) * prob;

      int ybin = hist->FindFixBin(ybincenter);
      double oldval = hist->GetBinContent(ybin);
      if(ybin != j) std::cout << "smth wrong...." << std::endl;
      double olderr = hist->GetBinError(ybin);

      hist->SetBinContent(ybin, oldval + con);
      hist->SetBinError(ybin, std::sqrt(olderr*olderr + err*err));
    }
  }
}

void MakeTruthHistFill(TH1D* hist, TH2D* h2d, TH1D* hncollinput){
  int nbinsx = h2d->GetNbinsX();
  int nbinsy = h2d->GetNbinsY();

  hist->Reset();
  hist->Sumw2();

  TH1D* hncoll_norm = (TH1D*) hncollinput->Clone("hncoll_norm");
  hncoll_norm->Scale(1./hncoll_norm->Integral()); // make as probability hist.
  std::cout << "hncoll norm integral " << hncoll_norm->Integral() << std::endl;

  for(int i=1; i<=nbinsx;++i){
    double prob = hncoll_norm->GetBinContent(i); 
    for(int j=1; j<=nbinsy;++j){
      double ybincenter = h2d->GetYaxis()->GetBinCenter(j);
      double con = h2d->GetBinContent(i,j) * prob;
      hist->Fill(ybincenter, con);
    }
  }
}

void CalcRCP(TH1D* hcent, TH1D* hperi, double ncent, double nperi, int ncollcent, int ncollperi){
  
  hcent->Sumw2();
  hperi->Sumw2();
  hcent->Divide(hperi);
  hcent->Scale(nperi/ncent);
  hcent->Scale((double)ncollperi/ncollcent);
}  

void CalcICP(TH1D* hcent, TH1D* hperi, double ncent, double nperi){
  
  hcent->Sumw2();
  hperi->Sumw2();
  hcent->Divide(hperi);
  hcent->Scale(nperi/ncent);
}  

void MakeBiasHist(TH1D* hreco, double nevent, TH1D* htruth){
  hreco->Divide(htruth);
  hreco->Scale(1./nevent);
  hreco->SetLineWidth(2);
  hreco->SetLineColor(kBlack);
}  

void PartialRebin_Weighted(TH1D*& h, double threshold, int groupSize) {
    int nbins = h->GetNbinsX();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    int bin_thresh = h->FindBin(threshold);

    std::vector<double> edges;
    for (int i = 1; i < bin_thresh; ++i)
        edges.push_back(h->GetBinLowEdge(i));

    int bin = bin_thresh;
    while (bin <= nbins) {
        edges.push_back(h->GetBinLowEdge(bin));
        bin += groupSize;
    }
    edges.push_back(xmax); 

    TH1D* hnew = (TH1D*)h->Clone(Form("%s_rebinned", h->GetName()));
    hnew->Reset();  
    hnew->SetBins(edges.size() - 1, edges.data());


    for (int i = 1; i <= hnew->GetNbinsX(); ++i) {
        double low = hnew->GetBinLowEdge(i);
        double high = hnew->GetBinLowEdge(i + 1);

        int binLow = h->FindBin(low + 1e-6);
        int binHigh = h->FindBin(high - 1e-6);

        double weightedSum = 0;
        double sumWeights = 0;

        for (int j = binLow; j <= binHigh; ++j) {
            double y = h->GetBinContent(j);
            double err = h->GetBinError(j);
            if (err <= 0) continue;
            double w = 1.0 / (err * err);
            weightedSum += w * y;
            sumWeights += w;
        }

        if (sumWeights > 0) {
            double y_avg = weightedSum / sumWeights;
            double err_avg = std::sqrt(1.0 / sumWeights);
            hnew->SetBinContent(i, y_avg);
            hnew->SetBinError(i, err_avg);
        }
    }

    delete h;
    h = hnew;
}

void PartialRebinCountsFill(TH1D*& h, double threshold, int groupSize) {
    int nbins = h->GetNbinsX();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    int bin_thresh = h->FindBin(threshold);

    std::vector<double> edges;
    for (int i = 1; i < bin_thresh; ++i)
        edges.push_back(h->GetBinLowEdge(i));

    int bin = bin_thresh;
    while (bin <= nbins) {
        edges.push_back(h->GetBinLowEdge(bin));
        bin += groupSize;
    }
    edges.push_back(xmax);  

    TH1D* hnew = (TH1D*)h->Clone(Form("%s_rebinned", h->GetName()));
    hnew->Reset();  
    hnew->SetBins(edges.size() - 1, edges.data());

    for (int i = 1; i <= hnew->GetNbinsX(); ++i) {
        double low = hnew->GetBinLowEdge(i);
        double high = hnew->GetBinLowEdge(i + 1);

        int binLow = h->FindBin(low + 1e-6);
        int binHigh = h->FindBin(high - 1e-6);

        double sum = 0;
        double err2 = 0;

        for (int j = binLow; j <= binHigh; ++j) {
            sum += h->GetBinContent(j);
            err2 += std::pow(h->GetBinError(j), 2);
        }

        hnew->SetBinContent(i, sum);
        hnew->SetBinError(i, std::sqrt(err2));
    }

    delete h;
    h = hnew;
}
