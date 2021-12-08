/*   Class for connecting to the flow coefficients database
//   and generating random angles with flow distribtions.
//   Author: Jakub Zielinski
//   Warsaw University of Technology
//   email: jakub.stanislaw.zielinski@cern.ch
*/

#include "FemtoFlowDatabase.h"
#define V_START 2 // not using 
#define V_STOP 3

Double_t flowHarmonics(Double_t *x, Double_t *par)
{
    Double_t val = 1.;
    for (Int_t i = V_START; i < V_STOP+1; i++)
    {
      if( i == 3){ // only for v3	
	val += 2. * par[i] * TMath::Cos(i * (x[0] - par[0]));
	continue;
      }
      val += 2. * par[i] * TMath::Cos(i * x[0]);
    }
    return val;
}

Double_t flowIntegral(Double_t *x, Double_t *par)
{
  Double_t val = x[0];
  for(Int_t i=V_START; i<V_STOP+1; i++)
  {
    if( i == 3){ // only for v3	
      val += 2. * par[i] / i * TMath::Sin(i * (x[0]-par[0]));
      continue;
    }  
    val += 2. * par[i] / i * TMath::Sin(i * x[0]);
  }
  return val;
}

/*
* Basic constructor for the class.
* Note: there have to be separate objects for particles with different pdg
*/
FemtoFlowDatabase::FemtoFlowDatabase( Int_t pdg, const char* tableName, Double_t energy, const char* centrality, const char* eta, const char* experiment)
{
  // Connect to the database server
  kServer = TSQLServer::Connect("mysql://localhost/flow","vm","Pass1234"); //connect to the database
  // Create functions for phi distribution and it's integral
  fFlow = new TF1(Form("flow_%d",pdg), flowHarmonics, -1.1*TMath::Pi(), 1.1*TMath::Pi(), 4); //create function for flowharmonics
  fFlow->SetParameters(0.0,0.0,0.0,0.0);
  fIntegral = new TF1(Form("int_%d", pdg), flowIntegral, -1.1*TMath::Pi(), 1.1*TMath::Pi(), 4);
  fIntegral->SetParameters(0.0,0.0,0.0,0.0);
  // Begin with 0 Vn found in the database (it is changed whenever DownloadGraphs is used"
  fNVn = 0;
  // Create an array for Vn 
  fVn = new double[3];
  fVn[0] = 0; //v1
  fVn[1] = 0; //v2
  fVn[2] = 0; //v3
  fPT = new double[2];
  // Setting up parameters for the database query
  fPDG = pdg;
  fTableName = tableName;
  fCentrality = centrality;
  fExperiment = experiment;
  fEnergy = energy;
  if(strcmp(eta,"*") == 0) fEta = "%%";
  else fEta = eta;

  // Graph position (points can be moved to low error and high error to calculate uncertainty)
  fGrPos = kCentral;
  // Graph to store V2 data
  fVnGraph[0] = new TGraphErrors();
  fVnGraph[0]->SetTitle(Form("v_{2} for %d;p_{T} (GeV);v_{2}", pdg));

  // Graph to store V3
  fVnGraph[1] = new TGraphErrors();
  fVnGraph[1]->SetTitle(Form("v_{3} for %d;p_{T} (GeV);v_{3}", pdg));

  // 2D Graph for statistics
  fStats = new TH2D(Form("hStats_%d",pdg), Form("#phi(pT) vs pT histogram (%d);pT (GeV); #phi(pT);entries",pdg),100,0,3,100,-TMath::Pi(),TMath::Pi());
}

/*
* Basic destructor of the class
*/
FemtoFlowDatabase::~FemtoFlowDatabase(){
  if(kServer != nullptr) delete kServer;
  if(fFlow != nullptr) delete fFlow;
  if(fVn != nullptr) delete fVn;
  if(fPT != nullptr) delete fPT;

  if(fVnGraph[0] != nullptr) delete fVnGraph[0];
  if(fVnGraph[1] != nullptr) delete fVnGraph[1];
  if(fVnGraph[1] != nullptr) delete fVnGraph[2];

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
  fNVn = 0;

  const char* vn_param[] = {"v2","v3"}; //for now only v2 and v3 stored in the database. We work on that later
  Int_t n; // number of rows found in the database for the given vn
  Double_t p;
  Double_t p_low = 0;
  Double_t p_high = 0;
  Double_t vn=0;
  Double_t statM = 0;
  Double_t statP = 0;
  Double_t sysM = 0;
  Double_t sysP = 0;
  this->fPT[0] = 10.0; //begining of the pT range
  this->fPT[1] = 0.0; // end of the pT range
  for(Int_t i=0; i<2;i++){
    TString sql_statement = Form("SELECT pT, pT_LOW, pT_HIGH, %s, %s_statM, %s_statP, %s_sysM, %s_sysP FROM %s WHERE %s IS NOT NULL AND energy = %f AND experiment = \"%s\" AND centrality = \"%s\" AND pdg = %d AND eta LIKE \"%s\";",
				 vn_param[i],vn_param[i],vn_param[i],vn_param[i],vn_param[i],fTableName,vn_param[i],fEnergy,fExperiment,fCentrality,fPDG, fEta);
    TSQLStatement* stmt = kServer->Statement(sql_statement.Data(),2048);
    n=0;

    if(stmt->Process()){
      stmt->StoreResult();
      while(stmt->NextResultRow()){
        p = stmt->GetDouble(0);
        p_low = stmt->GetDouble(1);
        p_high = stmt->GetDouble(2);
        vn  = stmt->GetDouble(3);
        if(vn < 1e-5) vn= 0;
        statM = stmt->GetDouble(4);
        statP = stmt->GetDouble(5);
        sysM = stmt->GetDouble(6);
        sysP = stmt->GetDouble(7);

        if(p<this->fPT[0]) this->fPT[0] = p; //updating lower range limit
        if(p>this->fPT[1]) this->fPT[1] = p; //updating higher range limit

	fVnGraph[i]->SetPoint(n,p,vn);
	// setting error bars according to the database
	fVnGraph[i]->SetPointError(n, (p_high - p_low) / 2., statP+sysP);
        n++;
      }
      for(int ip=0; ip<fVnGraph[i]->GetN(); ip++){
	if(ip==0) fVnGraph[i]->SetPointError(0, fVnGraph[i]->GetErrorX(1), fVnGraph[i]->GetErrorY(0));
	switch(fGrPos){
	case kCentral:
	  break;
	case kLow:
	  fVnGraph[i]->SetPoint(ip,fVnGraph[i]->GetPointX(ip)-fVnGraph[i]->GetErrorX(ip),fVnGraph[i]->GetPointY(ip));
	  break;
	case kHigh:
	  fVnGraph[i]->SetPoint(ip,fVnGraph[i]->GetPointX(ip)+fVnGraph[i]->GetErrorX(ip),fVnGraph[i]->GetPointY(ip));
	  break;
	default:
	  break;
	}
      }  
      if(n>0) fNVn++; // a parameter found
    }
    else{
      std::cout<<"ERROR: there was a problem receiving data from the database"<<std::endl;
    }
    stmt->Close();
  }
  return fNVn;
}

/*
* Returns random phi from the distribution according to the flow harmonics
*/
Double_t FemtoFlowDatabase::GetPhi(Double_t pT, Double_t psi3){
  //std::cout<<"\tGetting vns"<<std::endl;
  Double_t phi = 0.0;
  // update v2 and v3
  this->GetVns(pT);
  // calculate flow distributions
  // std::cout<<"\tv2:"<<this->fVn[1]<<"\tv3:"<<this->fVn[2]<<"\tPsi3:"<<fPsi3<<std::endl;
  fIntegral->SetParameters(psi3,0.0,this->fVn[1],this->fVn[2]);
  phi = fIntegral->GetX(gRandom->Rndm() * TMath::TwoPi() - TMath::Pi());
  // if(abs(phi)>=TMath::Pi()) std::cout<<"\tphi="<<phi<<std::endl;
  fStats->Fill(pT,phi);
  return phi;
}

/*******************************************
*GetVns(Double_t)* - function to set Vn parameters
(v1, v2, v3, etc) for the given pT (saved
in class stracture). Parameters are taken
from the databse.

**inputs*:
  double pT - value of the tranverse momentum
  for which we calculate the spherical harmonics

*******************************************/
  void FemtoFlowDatabase::GetVns(Double_t pT){
    
    fVn[0] = 0.0;  // - 0.75/0.8 * 10e-3 * eta; //(estimate for v1 for pion-kaon)
    fVn[1] = 0.0;  // v2
    fVn[2] = 0.0;  // v3
    
    for(int i=0; i<fNVn;i++){
      fVn[i+1] = fVnGraph[i]->Eval(pT);
    }

    //std::cout<<fVn[0]<<" "<<fVn[1]<<" "<<fVn[2]<<std::endl;
  }

TF1* FemtoFlowDatabase::GetFlowHarmonics(Double_t psi3){
  fFlow->SetParameters(psi3,this->fVn[0],this->fVn[1],this->fVn[2]);
  return fFlow;
}

TF1* FemtoFlowDatabase::GetFlowIntegral(Double_t psi3){
  fIntegral->SetParameters(psi3,this->fVn[0],this->fVn[1],this->fVn[2]);
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
    std::cout<<"\tgraph position:";
    switch(fGrPos){
    case kCentral:
      std::cout<<"kCentral"<<std::endl;
      break;
    case kLow:
      std::cout<<"kLow"<<std::endl;
      break;
    case kHigh:
      std::cout<<"kHigh"<<std::endl;
      break;
    default:
      std::cout<<"kCentral"<<std::endl;
      break;
    } 
  }else{
    std::cout<<"ERROR::database not connected"<<std::endl;
  }
}

void FemtoFlowDatabase::SetEta(const char* eta){
  if(strcmp(eta,"*") == 0) fEta = "%%";
  else fEta = eta;
}
