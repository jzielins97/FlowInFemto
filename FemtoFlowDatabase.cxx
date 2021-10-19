/*   Class for connecting to the flow coefficients database
//   and generating random angles with flow distribtions.
//   Author: Jakub Zielinski
//   Warsaw University of Technology
//   email: jakub.stanislaw.zielinski@cern.ch
*/

#include "FemtoFlowDatabase.h"

Double_t flowHarmonics(Double_t *x, Double_t *par)
{
    Double_t val = 1.;
    for (Int_t i = 0; i < 4; i++)
    {
        val += 2. * par[i] * TMath::Cos((i + 1.) * x[0]);
    }
    return val;
}

Double_t flowIntegral(Double_t *x, Double_t *par)
{
  Double_t val = x[0];
  for(Int_t i=0; i<4; i++)
  {
    val += 2. * par[i] / (i+1.) * TMath::Sin((i + 1.) * x[0]);
    //val += par[i] * TMath::ASin((i + 2.) * x[0] / 2.0) / (i + 2.);
  }
  return val;
}

/*
* Basic constructor for the class.
* Note: there have to be separate objects for particles with different pdg
*/
FemtoFlowDatabase::FemtoFlowDatabase( Int_t pdg, const char* tableName, Double_t energy, const char* centrality, const char* eta, const char* experiment)
{
  kServer = TSQLServer::Connect("mysql://localhost/flow","vm","Pass1234"); //connect to the database
  fFlow = new TF1(Form("flow_%d",pdg), flowHarmonics, -TMath::Pi(), TMath::Pi(), 4); //create function for flowharmonics
  fFlow->SetParameters(0.0,0.0,0.0,0.0);
  fIntegral = new TF1(Form("int_%d", pdg), flowIntegral, -TMath::Pi(), TMath::Pi(), 4);
  fIntegral->SetParameters(0.0,0.0,0.0,0.0);
  fVm = new double[4];
  fVm[0] = 0;
  fVm[1] = 0;
  fVm[2] = 0;
  fVm[3] = 0;
  fPT = new double[2];
  fPDG = pdg;
  fTableName = tableName;
  fCentrality = centrality;
  fExperiment = experiment;
  fEnergy = energy;
  fEta = eta;

  fVmGraph[0] = new TGraphErrors();
  fVmGraph[0]->SetTitle(Form("v_{2} for %d;p_{T} (GeV);v_{2}", pdg));

  fVmGraph[1] = new TGraphErrors();
  fVmGraph[1]->SetTitle(Form("v_{3} for %d;p_{T} (GeV);v_{3}", pdg));

  fStats = new TH2D(Form("hStats_%d",pdg), Form("#phi(pT) vs pT histogram (%d);pT (GeV); #phi(pT);entries",pdg),100,0,3,100,-TMath::Pi(),TMath::Pi());
}

/*
* Basic destructor of the class
*/
FemtoFlowDatabase::~FemtoFlowDatabase(){
  if(kServer != nullptr) delete kServer;
  if(fFlow != nullptr) delete fFlow;
  if(fVm != nullptr) delete fVm;
  if(fPT != nullptr) delete fPT;

  if(fVmGraph[0] != nullptr) delete fVmGraph[0];
  if(fVmGraph[1] != nullptr) delete fVmGraph[1];
  if(fVmGraph[1] != nullptr) delete fVmGraph[2];

  if(fStats != nullptr) delete fStats;
}

/*
* Function connects to the database and downloads all vn(pT), then stores them
* as TGraphErrors objects inside the class object.
* Parameters are selected according to the fTableName, fExperiment, fEnergy,
* fCentrality and fPDG
* It returns number of parameters vn for which data was found in the database
*/
Int_t FemtoFlowDatabase::DownloadGraphs(){
  Int_t parameters = 0;

  const char* vm_param[4] = {"v2","v3", "v4", "v5"}; //for now only v2 and v3 stored in the database. We work on that later
  Int_t n; // number of rows found in the database for the given vn
  Double_t p;
  Double_t p_low = 0;
  Double_t p_high = 0;
  Double_t vm=0;
  Double_t statM = 0;
  Double_t statP = 0;
  Double_t sysM = 0;
  Double_t sysP = 0;
  this->fPT[0] = 10.0; //begining of the pT range
  this->fPT[1] = 0.0; // end of the pT range
  for(Int_t i=0; i<1;i++){ //not looking for v3 anymore
    TString sql_statement = Form("SELECT pT, pT_LOW, pT_HIGH, %s, %s_statM, %s_statP, %s_sysM, %s_sysP FROM %s WHERE %s IS NOT NULL AND energy = %f AND experiment = \"%s\" AND centrality = \"%s\" AND pdg = %d;" /*AND eta = \"%s\";"*/,
				 vm_param[i],vm_param[i],vm_param[i],vm_param[i],vm_param[i],fTableName,vm_param[i],fEnergy,fExperiment,fCentrality,fPDG); //fEta
    TSQLStatement* stmt = kServer->Statement(sql_statement.Data(),2048);
    n=0;

    if(stmt->Process()){
      stmt->StoreResult();
      while(stmt->NextResultRow()){
        p = stmt->GetDouble(0);
        p_low = stmt->GetDouble(1);
        p_high = stmt->GetDouble(2);
        vm  = stmt->GetDouble(3);
        if(vm < 1e-5) vm = 0;
        statM = stmt->GetDouble(4);
        statP = stmt->GetDouble(5);
        sysM = stmt->GetDouble(6);
        sysP = stmt->GetDouble(7);

        if(p<this->fPT[0]) this->fPT[0] = p; //updating lower range limit
        if(p>this->fPT[1]) this->fPT[1] = p; //updating higher range limit

        fVmGraph[i]->SetPoint(n,p,vm);
        // setting error bars according to the database
        if(p == p_high) fVmGraph[i]->SetPointError(n, (p_high - p_low)/2.0,statP+sysP);
        else fVmGraph[i]->SetPointError(n, p_high - p, statP+statM);
        n++;
      }
      if(n>0) parameters++; // a parameter found
    }
    else{
      std::cout<<"ERROR: there was a problem receiving data from the database"<<std::endl;
    }
    stmt->Close();
  }
  return parameters;
}

/*
* Returns random phi from the distribution according to the flow harmonics
*/
Double_t FemtoFlowDatabase::GetPhi(Double_t pT, Double_t eta = 0.0){
  //std::cout<<"\tGetting vms"<<std::endl;
  Double_t phi = 0.0;
  this->GetVms(pT, eta);

  fIntegral->SetParameters(this->fVm[0],this->fVm[1],this->fVm[2],this->fVm[3]);
  phi = fIntegral->GetX(gRandom->Rndm() * TMath::TwoPi() - TMath::Pi());
  fStats->Fill(pT,phi);
  return phi;
}

/*******************************************
*getVms(Double_t)* - function to set Vm parameters
(v1, v2, v3, etc) for the given pT (saved
in class stracture). Parameters are taken
from the databse.

**inputs*:
  double pT - value of the tranverse momentum
  for which we calculate the spherical harmonics

*******************************************/
  void FemtoFlowDatabase::GetVms(Double_t pT, Double_t eta){
    
    fVm[0] = 0.0; //- 0.75/0.8 * 10e-3 * eta; (estimate for v1 for pion-kaon)
    fVm[1] = 0.0;  // v2
    fVm[2] = 0.0;  // v3 - not using anymore (19.10.2021 jzielins)
    
    for(int i=0; i<1;i++){
      fVm[i+1] = fVmGraph[i]->Eval(pT);
    }

    fVm[3] = 0.0;
    //std::cout<<fVm[0]<<" "<<fVm[1]<<" "<<fVm[2]<<" "<<fVm[3]<<std::endl;
  }

TF1* FemtoFlowDatabase::GetFlowHarmonics(){
  fFlow->SetParameters(this->fVm[0],this->fVm[1],this->fVm[2],this->fVm[3]);
  return fFlow;
}

TF1* FemtoFlowDatabase::GetFlowIntegral(){
  fIntegral->SetParameters(this->fVm[0],this->fVm[1],this->fVm[2],this->fVm[3]);
  return fIntegral;
}

/*
* Prints all the information, which it will use while retrieving data from the database
* Also says if it managed to connect to the database
*/

void FemtoFlowDatabase::ShowParams(){
  if(kServer->IsConnected()){
    std::cout<<"The flow database is connected"<<std::endl;
    std::cout<<"Your settings: "<<std::endl;
    std::cout<<"\ttable: "<<fTableName<<std::endl;
    std::cout<<"\texperiment: "<<fExperiment<<std::endl;
    std::cout<<"\tenergy: "<<fEnergy<<std::endl;
    std::cout<<"\tcentrality: "<<fCentrality<<std::endl;
    std::cout<<"\teta: "<<fEta<<std::endl;
  }else{
    std::cout<<"ERROR::database not connected"<<std::endl;
  }
}
