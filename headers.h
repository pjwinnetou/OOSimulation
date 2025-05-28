#include <iostream>
#include <cmath>
#include <TROOT.h>

std::vector<double> load_boundaries(const std::string& filename) {
  std::vector<double> bins;
  std::ifstream fin(filename);
  double val;
  while (fin >> val) {
    bins.push_back(val);
  }   
  return bins;
}   

int get_quantile_bin(double val, const std::vector<double>& boundaries) {
  auto it = std::upper_bound(boundaries.begin(), boundaries.end(), val);
  int bin = std::max(0, std::min((int)(it - boundaries.begin()) - 1, (int)boundaries.size() - 2));
  return bin;
}   

Double_t calcEt(double pt, double eta, double e){ 
  double massSquared = e*e - pt*pt * std::cosh(eta) * std::cosh(eta);
  if (massSquared < 0) massSquared = 0;
  double mass = std::sqrt(massSquared);
  double Et = std::sqrt(pt*pt + mass*mass) / std::cosh(eta);
  return Et; 
}

