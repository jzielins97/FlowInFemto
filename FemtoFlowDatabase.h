/*   Class for connecting to the flow coefficients database
//   and generating random angles with flow distribtions.
//   Author: Jakub Zielinski
//   Warsaw University of Technology
//   email: jakub.stanislaw.zielinski@cern.ch
*/

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

class FemtoFlowDatabase {
private:
  //Data memebers
  Int_t fPDG; // pdg of the particle we want to get vm for
  TSQLServer* kServer; // pointer to the connected MySQL server with the database
  const char* fTableName; // name of the table storing vm in the database
  const char* fExperiment; // experiment from which measurements were made
  Double_t fEnergy; // collision sqrt(sNN) in GeV
  const char* fCentrality; //centrality of the collision
  const char* fEta; // absolute pseudorapidity range
  
  TF1 *fFlow; // flow harmonics for the given pT
  TF1 *fIntegral; // distribution for drawing random phi (based on the flow harmonics)
  
  Double_t* fPT; // table with pT range (stores what's the lowest and highest pT for which v parameters were found in the database)
  Int_t fGrPos;
  Double_t* fVn; // table with vn parameters at given time
  Double_t fNVn; // number of Vn found in the database

  TGraphErrors* fVnGraph[2]; // graph to store values of v2 and v3 from the database
  TH2D* fStats; // 2D histogram filled with pairs (pT,phi(pT)) where phi(pT) are returns from eval from flow harmonics

  //Methods
  void GetVns(Double_t pT);

public:
  enum EVnName {kV2, kV3};
  enum EGraphPosition {kCentral, kLow, kHigh};
  //constructor with specified particle pdg
  FemtoFlowDatabase( Int_t pdg,
                     const char* tableName = "PbPb",
                     Double_t energy = 2760,
                     const char* centrality = "0-5%",
		     const char* eta = "<0.8",
                     const char* experiment = "ALICE");
  ~FemtoFlowDatabase();
  //public methods
  Int_t DownloadGraphs(); //filling graphs with v parameters from the database
  Double_t GetPhi(Double_t pT, Double_t psi3=0.0); //returns a random phi with a distribution from spherical harmonics (pT - momentum, psi3 - Reaction Plane angle for v3)
  void ShowParams(); // prints parametrs used while connecting with the database
//setters
  void SetCentrality(const char* centrality){ fCentrality = centrality; };
  void SetEnergy(Double_t energy){ fEnergy = energy; };
  void SetEta(const char* eta);
  void SetExperimentName(const char* experiment){ fExperiment = experiment; };
  void SetGraphPosition(Int_t pos) { fGrPos = pos; };
  void SetTableName(const char* tableName){ fTableName = tableName; };
//getters
  const char* GetCentrality(){ return fCentrality; };
  Double_t GetEnergy(){ return fEnergy; };
  const char* GetEta(){ return fEta; };
  const char* GetExperimentName(){ return fExperiment; };
  Int_t GetPDG(){ return fPDG; };
  const char* GetTableName(){ return fTableName; };
  Int_t GetNVn(){ return fNVn; };
  Double_t GetV1(){ return fVn[0]; };
  Double_t GetV2(){ return fVn[1]; };
  Double_t GetV3(){ return fVn[2]; };
  // Double_t GetV4(){ return fVn[3];};
  /******/

  TF1* GetFlowHarmonics(Double_t psi3=0.0);
  TF1* GetFlowIntegral(Double_t psi3=0.0);

  /* Returns a pointer to TGraphErrors object with a graph of the vn parameter (n = vParam)
   which was created with the database */
  TGraphErrors* GetGraph(Int_t vParam){ return fVnGraph[vParam]; }; //vParam = 0 -> v2; vParam = 1 -> v3 
  TH2D* GetStatsHisto(){ return fStats; }; // returns 2D histogram of all pairs pT and phi(pT)
};

Double_t flowHarmonics(Double_t *x, Double_t *par); // function to generate angle distribution according to the flow (v2 and v3)
Double_t flowIntegral(Double_t *x, Double_t *par);  // function to get phi with distribution from flowHarmonics
#endif
