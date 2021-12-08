#include "FemtoFlowDatabase.h"
#include "CorrFctnDirectYlm.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TList.h"
#include "TLorentzVector.h"
//#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "TRandom.h"
#include "TStopwatch.h"

Double_t GetMass(int pdg){
  // if (pdg == 211) return 0.139570;
  // if (pdg == 321) return 0.493677;
  // if (pdg == 2212) return 0.938272;
  // if (pdg == 2112) return 0.939565;
  // if (pdg == 123456) return 2*0.939565;
  // return 0;
  switch(pdg){
  case 211: //pion
    return 0.139570;
  case 321: //kaon
    return 0.493677;
  case 2212: //proton
    return 0.938272;
  case 2112: //neutron
    return 0.939565;
  case 1000010020: //deuteron
    return 2*0.939565;
  default:
    return 0;
  }
}

int main(int argc, char** argv){
  gRandom->SetSeed(0); //set new seed each time you run the program
  int nBins = 100; //50;
  int pdg1, pdg2;
  TFile* fin1;
  TFile* fin2;
  TList* fOutput = new TList();

  TH1F* h_pt1=0;
  TH1F* h_pt2=0;

  TH2D* h_scan_part1; // pT of particle one with k* < 0.04
  TH2D* h_scan_part2; // pT of particle one with k* < 0.04

  CorrFctnDirectYlm* fCylm;
  FemtoFlowDatabase* dFlowPart1;
  FemtoFlowDatabase* dFlowPart2;
  Int_t grPos = FemtoFlowDatabase::kCentral;

  Double_t pt1_min = 0.19;
  Double_t pt1_max = 1.5;
  Double_t pt2_min = 0.19;
  Double_t pt2_max = 1.5;

  // Double_t massPi = 0.139570;
  // Double_t massK = 0.493677;
  
  const char* centrality = "0-5%";
  const char* experiment = "ALICE";
  const char* eta = "<0.8";
  Double_t sNN = 2760; // energy for of the collision

  if(argc < 6){
    std::cout<<"Error: not enough arguments. Needs: 5, was given: "<<argc-1<<"."<<std::endl;
    std::cout<<"<particle #1 pdg> <pT distribution file #1> <particle #2 pdg> <pT distribution #2> <N> --optional <Graph Position> <centrality> <eta> <energy> <experiment>"<<std::endl;
    return -1;
  }
  if(argc > 6) centrality = argv[6];
  if(argc > 7) grPos = atoi(argv[7]);
  if(argc > 8) eta = argv[8];
  if(argc > 9) sNN = atof(argv[9]);
  if(argc > 10) experiment = argv[10];

//setting up first particle
  //database:
  pdg1 = atoi(argv[1]);
  Double_t massPart1 = GetMass(pdg1);
  dFlowPart1 = new FemtoFlowDatabase(TMath::Abs(pdg1));
  dFlowPart1->SetExperimentName(experiment);
  dFlowPart1->SetEnergy(sNN);
  dFlowPart1->SetEta(eta);
  dFlowPart1->SetTableName("PbPb");
  dFlowPart1->SetCentrality(centrality);
  dFlowPart1->SetGraphPosition(grPos);
  dFlowPart1->ShowParams();
  std::cout<<"\tm1="<<massPart1<<std::endl;
  int foundVnForPart1 = dFlowPart1->DownloadGraphs();
  std::cout<<"\tFor particle "<<pdg1<<" found "<<foundVnForPart1<<" v params"<<std::endl;
  
  //geting pT distribution histogram
  fin1 = new TFile(argv[2]);
  if(!fin1->GetListOfKeys()->Contains("hpt")){
    std::cout<<"Error: there is no object called \"hpt\" in "<<argv[2]<<". Make sure that momentum distribution histogram is called \"hpt\"."<<std::endl;
    return -1;
  }
  h_pt1 = (TH1F*)fin1->Get("hpt")->Clone("hpt1");
  std::cout<<"\tThere are "<<h_pt1->GetEntries()<<" entries in the pT histogram"<<std::endl;  

//setting up second particle
  //database:
  pdg2 = atoi(argv[3]);
  Double_t massPart2 = GetMass(pdg2);
  dFlowPart2 = new FemtoFlowDatabase(TMath::Abs(pdg2));
  dFlowPart2->SetExperimentName(experiment);
  dFlowPart2->SetEnergy(sNN);
  dFlowPart2->SetEta(eta);
  dFlowPart2->SetTableName("PbPb");
  dFlowPart2->SetCentrality(centrality);
  dFlowPart2->SetGraphPosition(grPos);
  dFlowPart2->ShowParams();
  std::cout<<"\tm2="<<massPart2<<std::endl;
  int foundVnForPart2 = dFlowPart2->DownloadGraphs();
  std::cout<<"\tFor particle "<<pdg2<<" found "<<foundVnForPart2<<" v params"<<std::endl;
  
  //geting pT distribution histogram  
  fin2 = new TFile(argv[4]);
  if(!fin2->GetListOfKeys()->Contains("hpt")){
    std::cout<<"Error: there is no object called \"hpt\" in "<<argv[2]<<". Make sure that momentum distribution histogram is called \"hpt\"."<<std::endl;
    return -1;
  }
  h_pt2 = (TH1F*)fin2->Get("hpt")->Clone("hpt2");
  std::cout<<"\tThere are "<<h_pt2->GetEntries()<<" entries in the pT histogram"<<std::endl;
/* end of setting up database ********************/

  h_scan_part1 = new TH2D(Form("h_scan1_%d",pdg1), Form("%d;p_{T} [GeV/c^{2}];k* [GeV/c^{2}]",pdg1), (int)((pt1_max-pt1_min)*100),pt1_min,pt1_max,16,0.0,0.04);

  h_scan_part2 = new TH2D(Form("h_scan2_%d",pdg2), Form("%d;p_{T} [GeV/c^{2}];k* [GeV/c^{2}]",pdg2), (int)((pt2_max-pt2_min)*100),pt2_min,pt2_max,16,0.0,0.04);
  
  int N = (int)atof(argv[5]);
  std::cout<<"N="<<N<<std::endl;

  fCylm = new CorrFctnDirectYlm("Cylm",1,nBins,0.0,1);
  fCylm->AddToList(fOutput);

  ROOT::Math::PtEtaPhiM4D<double> v1(1, 1, 1, 1);
  ROOT::Math::PtEtaPhiM4D<double> v2(1, 1, 1, 1);
  ROOT::Math::PtEtaPhiM4D<double> vsum(1, 1, 1, 1);

  TLorentzVector v1b;
  TLorentzVector v2b;
  TLorentzVector vsumb;
  TStopwatch* timer = new TStopwatch();
  timer->Start();
  std::cout << "Loop for filling the random" << std::endl;
  for (int i = 0; i < N; i++)
  {
    Double_t psi3 = gRandom->Rndm() * TMath::TwoPi() - TMath::Pi(); // Reaction Plane angle for v3
    
    Double_t pt1 = h_pt1->GetRandom();
    Double_t pt2 = h_pt2->GetRandom();
    while(pt1 < pt1_min || pt1 > pt1_max) pt1 = h_pt1->GetRandom();
    while(pt2 < pt2_min || pt2 > pt2_max) pt2 = h_pt2->GetRandom();

    Double_t eta1 = 2. * (gRandom->Rndm() - 0.5) * 0.8;
    Double_t eta2 = 2. * (gRandom->Rndm() - 0.5) * 0.8;
    
    Double_t phi1 = dFlowPart1->GetPhi(pt1, psi3);
    Double_t phi2 = dFlowPart2->GetPhi(pt2, psi3);


    v1.SetCoordinates(pt1, eta1, phi1, massPart1);
    v2.SetCoordinates(pt2, eta2, phi2, massPart2);

    v1b.SetPxPyPzE(v1.Px(), v1.Py(), v1.Pz(), sqrt(massPart1 * massPart1 + v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz()));
    v2b.SetPxPyPzE(v2.Px(), v2.Py(), v2.Pz(), sqrt(massPart2 * massPart2 + v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz()));
    
    // sum the 4-vector components of both particles
    //total px, py, pz and E
    Double_t tpx = v1b.Px() + v2b.Px();
    Double_t tpy = v1b.Py() + v2b.Py();
    Double_t tpz = v1b.Pz() + v2b.Pz();
    Double_t te  = v1b.E()  + v2b.E();
    //total pt and mt
    Double_t tpt = tpx * tpx + tpy * tpy;
    Double_t tmt = te * te - tpz * tpz;
    Double_t tm  = sqrt(tmt - tpt);
    tmt = sqrt(tmt);
    tpt = sqrt(tpt);
    // bt
    Double_t tbt = tpt / tmt;
    //calculation of ko, ks, kl
    Double_t tbm = tpz / te;
    Double_t tgm = te / tmt;
    Double_t mkl = tgm * (v1b.Pz() - tbm * v1b.E());
    Double_t met = tgm * (v1b.E() - tbm * v1b.Pz());

    Double_t mko = (v1b.Px() * tpx + v1b.Py() * tpy) / tpt;
    Double_t mks = (-v1b.Px() * tpy + v1b.Py() * tpx) / tpt;
    // Double_t mkol = mko;
    mko = (tmt / tm) * (mko - (tpt / tmt) * met);

    Double_t mkv = TMath::Sqrt(mko*mko + mks*mks + mkl*mkl);
    if(mkv<1) fCylm->AddRealPair(mko,mks,mkl,1.0);
    if(mkv<0.04){
      h_scan_part1->Fill(pt1, mkv);
      h_scan_part2->Fill(pt2, mkv);
    }

    /***** creating den histogram **************/
    phi1 = gRandom->Rndm()*TMath::TwoPi() - TMath::Pi();
    phi2 = gRandom->Rndm()*TMath::TwoPi() - TMath::Pi();

    v1.SetCoordinates(pt1, eta1, phi1, massPart1);
    v2.SetCoordinates(pt2, eta2, phi2, massPart2);

    v1b.SetPxPyPzE(v1.Px(), v1.Py(), v1.Pz(), sqrt(massPart1 * massPart1 + v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz()));
    v2b.SetPxPyPzE(v2.Px(), v2.Py(), v2.Pz(), sqrt(massPart2 * massPart2 + v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz()));

    // sum the 4-vector components of both particles
    //total px, py, pz and E
    tpx = v1b.Px() + v2b.Px();
    tpy = v1b.Py() + v2b.Py();
    tpz = v1b.Pz() + v2b.Pz();
    te  = v1b.E()  + v2b.E();
    //total pt and mt
    tpt = tpx * tpx + tpy * tpy;
    tmt = te * te - tpz * tpz;
    tm  = sqrt(tmt - tpt);
    tmt = sqrt(tmt);
    tpt = sqrt(tpt);
    // bt
    tbt = tpt / tmt;
    //calculation of ko, ks, kl
    tbm = tpz / te;
    tgm = te / tmt;
    mkl = tgm * (v1b.Pz() - tbm * v1b.E());
    met = tgm * (v1b.E() - tbm * v1b.Pz());

    mko = (v1b.Px() * tpx + v1b.Py() * tpy) / tpt;
    mks = (-v1b.Px() * tpy + v1b.Py() * tpx) / tpt;
    // mkol = mko;
    mko = (tmt / tm) * (mko - (tpt / tmt) * met);
          
    mkv = TMath::Sqrt(mko*mko + mks*mks + mkl*mkl);
    if(mkv < 1) fCylm->AddMixedPair(mko,mks,mkl,1.0);

    if((i+1)%(int)(0.01*N) == 0){
      printf("\r%6d/%d",i+1,N);
      fflush(stdout);
    }
    //printProgress(i+1,N);
  }
  timer->Stop();
  std::cout<<std::endl;
  std::cout<<"time: "<<timer->RealTime()<<"; CPU: "<<timer->CpuTime()<<std::endl;
  fCylm->FillCovariances();
  fCylm->Finish();
  // fCylm->AddToList(fOutput);

  TFile* fout = new TFile("output_cylm.root","RECREATE");
  fCylm->Write();
  fout->Close();
  
  std::cout<<"Saving my list"<<std::endl;
  fout = new TFile("output_list.root","RECREATE");
  fout->cd();
  fOutput->Write();
  // fout->ls();
  fout->Close();
  
  fOutput->Clear();

  std::cout<<"Saving database statistics"<<std::endl;
  fout = new TFile("database_stats.root","RECREATE");
  fout->cd();
  dFlowPart1->GetGraph(FemtoFlowDatabase::kV2)->Write("hV2_pi");
  dFlowPart1->GetGraph(FemtoFlowDatabase::kV3)->Write("hV3_pi");
  dFlowPart1->GetStatsHisto()->Write("hStats_pi");
  h_scan_part1->Write();
  dFlowPart2->GetGraph(FemtoFlowDatabase::kV2)->Write("hV2_K");
  dFlowPart2->GetGraph(FemtoFlowDatabase::kV3)->Write("hV3_K");
  dFlowPart2->GetStatsHisto()->Write("hStats_K");
  h_scan_part2->Write();
  fout->Close();

  //fCylm->Wrie();

  delete fout;
  delete h_pt1;
  delete h_pt2;
  delete fin1;
  delete fin2;
  std::cout<<"done"<<std::endl;



  return 0;
}
