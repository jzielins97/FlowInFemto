#ifndef _FEMTOFLOWDATABASE_H_
#define _FEMTOFLOWDATABASE_H_

#include <math.h>
#include <iostream>

//ROOT includes
#include <TSQLServer.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TString.h>
#include <TSQLStatement.h>
#include <TRandom.h>

class FemtoFlowDataBase {
public:
  //constructor with specified particle pdg
  FemtoFlowDataBase( int pdg,
                     const char* tableName = "PbPb",
                     double energy = 2760,
                     const char* centrality = "0-5%",
                     const char* experiment = "ALICE");
  ~FemtoFlowDataBase();
  //public methods
  double getPhi(double pT); //returns a random phi with a distribution from spherical harmonics
  TGraphErrors* getGraph(int vParam); //returns a graph of vn(pT) from the database
  TH2D* getStatsHisto(); // returns 2D histogram of all pairs pT and phi(pT)
  void showParams(); // prints parametrs used while connecting with the database
  //setters
  void setTableName(const char* tableName);
  void setExperiment(const char* experiment);
  void setEnergy(double energy);
  void setCentrality(const char* centrality);
  int downloadGraphs();


private:
  void getVms(double pT);


  int kPdg; // pdg of the particle we want to get vm for
  TSQLServer* kServer; // pointer to the connected MySQL server with the database
  const char* kTableName; // name of the table storing vm in the database
  const char* kExperiment; //experiment from which measurements were made
  Double_t kEnergy; //collision sqrt(sNN) in GeV
  const char* kCentrality; //centrality of the collision
  TF1 *kFlow; //spherical harmonics for the given pT
  TF1 *kDist; //distribution for drawing random phi
  double* kPT; //table with pT range (stores what's the lowest and highest pT for which v parametrs were found in the database)
  double* kVm; //table with vm parameters at given time

  TGraphErrors* gVm[2]; // graph to store values of v2 and v3 from the database
  TH2D* hStats; // 2D histogram filled with pairs (pT,phi(pT)) where phi(pT) are returns from eval from flow harmonics
};

double flowHarmonics(double *x, double *par); // function to generate angle distribution according to the flow (v2 and v3)
double flowDistribution(double *x, double *par);  // function to get phi with distribution from flowHarmonics

#endif
