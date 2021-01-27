#include "FemtoFlowDataBase.h"

double flowHarmonics(double *x, double *par)
{
    double val = 1.;
    for (int i = 0; i < 4; i++)
    {
        val += 2. * par[i] * TMath::Cos((i + 2.) * x[0]);
    }
    return val;
}

double flowDistribution(double *x, double *par)
{
  double val = x[0];
  for(int i=0; i<4; i++)
  {
    val += 2. * par[i] / (i+2.) * TMath::Sin((i + 2.) * x[0]);
    //val += par[i] * TMath::ASin((i + 2.) * x[0] / 2.0) / (i + 2.);
  }
  return val;
}


/*
* Basic constructor for the class.
* Note: there have to be separate objects for particles with different pdg
*/
FemtoFlowDataBase::FemtoFlowDataBase( int pdg, const char* tableName, double energy, const char* centrality,const char* experiment)
{
  kServer = TSQLServer::Connect("mysql://localhost/flow","vm","P4$$w0rd"); //connect to the database
  // kFlow = new TF1(Form("flow_%d",pdg), flowHarmonics, 0, TMath::TwoPi(), 4); //create function for flowharmonics
  kFlow = new TF1(Form("flow_%d",pdg), flowHarmonics, 0, TMath::TwoPi(), 4); //create function for flowharmonics
  kDist = new TF1(Form("dist_%d", pdg), flowDistribution, 0, 1, 4);
  kFlow->SetParameters(0.0,0.0,0.0,0.0);
  kDist->SetParameters(0.0,0.0,0.0,0.0);
  kVm = new double[4];
  kPT = new double[2];
  kPdg = pdg;
  kTableName = tableName;
  kCentrality = centrality;
  kExperiment = experiment;
  kEnergy = energy;

  gVm[0] = new TGraphErrors();
  gVm[0]->SetTitle(Form("v_{2} for %d;p_{T} (GeV);v_{2}", pdg));

  gVm[1] = new TGraphErrors();
  gVm[1]->SetTitle(Form("v_{3} for %d;p_{T} (GeV);v_{3}", pdg));

  hStats = new TH2D(Form("hStats_%d",pdg), Form("#phi(pT) vs pT histogram (%d);pT (GeV); #phi(pT);entries",pdg),100,0,10,100,0.0,TMath::TwoPi());
}

/*
* Basic destructor of the class
*/
FemtoFlowDataBase::~FemtoFlowDataBase(){
  if(kServer != nullptr) delete kServer;
  if(kFlow != nullptr) delete kFlow;
  if(kVm != nullptr) delete kVm;
  if(kPT != nullptr) delete kPT;

  if(gVm[0] != nullptr) delete gVm[0];
  if(gVm[1] != nullptr) delete gVm[1];

  if(hStats != nullptr) delete hStats;
}


/*
* Function connects to the database and downloads all vn(pT), then stores them
* as TGraphErrors objects inside the class object.
* Parameters are selected according to the kTableName, kExperiment, kEnergy,
* kCentrality and pdg
* It returns number of parameters vn for which data was found in the database
*/
int FemtoFlowDataBase::downloadGraphs(){
  int parameters = 0;

  const char* vm_param[4] = {"v2","v3", "v4", "v5"}; //for now only v2 and v3 stored in the database. We work on that later
  int n; // number of rows found in the database for the given vn
  double p;
  double p_low = 0;
  double p_high = 0;
  double vm=0;
  double statM = 0;
  double statP = 0;
  double sysM = 0;
  double sysP = 0;
  this->kPT[0] = 10.0; //begining of the pT range
  this->kPT[1] = 0.0; // end of the pT range
  for(int i=0; i<2;i++){
    TString sql_statement = Form("SELECT pT, pT_LOW, pT_HIGH, %s, %s_statM, %s_statP, %s_sysM, %s_sysP FROM %s WHERE %s IS NOT NULL AND energy = %f AND experiment = \"%s\" AND centrality = \"%s\" AND pdg = %d;",
                                  vm_param[i],vm_param[i],vm_param[i],vm_param[i],vm_param[i],kTableName,vm_param[i],kEnergy,kExperiment,kCentrality,kPdg);
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

        if(p<this->kPT[0]) this->kPT[0] = p; //updating lower range limit
        if(p>this->kPT[1]) this->kPT[1] = p; //updating higher range limit

        gVm[i]->SetPoint(n,p,vm);
        // setting error bars according to the database
        if(p == p_high) gVm[i]->SetPointError(n, (p_high - p_low)/2.0,statP+sysP);
        else gVm[i]->SetPointError(n, p_high - p, statP+statM);
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
double FemtoFlowDataBase::getPhi(double pT){
  double phi = 0.0;
  this->getVms(pT);

  kFlow->SetParameters(this->kVm[0],this->kVm[1],this->kVm[2],this->kVm[3]);
  phi = kDist->Eval(gRandom->Rndm());
  hStats->Fill(pT,phi);
  return phi;
}

/*******************************************
*getVms(double)* - function to set Vm parameters
(v1, v2, v3, etc) for the given pT (saved
in class stracture). Parameters are taken
from the databse.

**inputs*:
  double pT - value of the tranverse momentum
  for which we calculate the spherical harmonics

*******************************************/
void FemtoFlowDataBase::getVms(double pT){

  kVm[0] = 0.05; // v2
  kVm[1] = 0.0;  // v3
  kVm[2] = 0.0;  // v4 is not in the database (yet)
  kVm[3] = 0.0;  // v5 is not in the database (yet)

  for(int i=0; i<2;i++){
    kVm[i] = gVm[i]->Eval(pT);
  }
}

/*
* Returns a pointer to TGraphErrors object with a graph of the vn parameter (n = vParam)
* which was created with the database
*/
TGraphErrors* FemtoFlowDataBase::getGraph(int vParam){
  return gVm[vParam]; //vParam = 0 -> v2; vParam = 1 -> v3
}

/*
* Prints all the information, which it will use while retrieving data from the database
* Also says if it managed to connect to the database
*/

void FemtoFlowDataBase::showParams(){
  if(kServer->IsConnected()){
    std::cout<<"The flow database is connected"<<std::endl;
    std::cout<<"Your settings: "<<std::endl;
    std::cout<<"\ttable: "<<kTableName<<std::endl;
    std::cout<<"\texperiment: "<<kExperiment<<std::endl;
    std::cout<<"\tenergy: "<<kEnergy<<std::endl;
    std::cout<<"\tcentrality: "<<kCentrality<<std::endl;
  }else{
    std::cout<<"ERROR::database not connected"<<std::endl;
  }
}

void FemtoFlowDataBase::setTableName(const char* tableName){
  kTableName = tableName;
}

void FemtoFlowDataBase::setExperiment(const char* experiment){
  kExperiment = experiment;
}

void FemtoFlowDataBase::setEnergy(double energy){
  kEnergy = energy;
}

void FemtoFlowDataBase::setCentrality(const char* centrality){
  kCentrality = centrality;
}

TH2D* FemtoFlowDataBase::getStatsHisto(){
  return hStats;
}
