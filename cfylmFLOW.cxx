#include "FemtoFlowDataBase.h"
#include "CorrFctnDirectYlm.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "TRandom.h"
#include "TStopwatch.h"

int main(int argc, char** argv){
  int pdg1, pdg2;
  TFile* fin1;
  TFile* fin2;
  TList* fOutput = new TList();

  TH1F* h_pt1=0;
  TH1F* h_pt2=0;

  CorrFctnDirectYlm* fCylm;
  FemtoFlowDataBase* dFlowPart1;
  FemtoFlowDataBase* dFlowPart2;

  Double_t massPi = 0.139570;
  Double_t massK = 0.493677;

  const char* centrality = "0-5%";

  if(argc < 6){
    std::cout<<"Error: not enough arguments. Needs: 4, was given: "<<argc-1<<"."<<std::endl;
    std::cout<<"<particle #1 pdg> <pT distribution file #1> <particle #2 pdg> <pT distribution #2> <N>"<<std::endl;
    return -1;
  }
  if(argc > 6) centrality = argv[6];

//setting up first particle
  //database:
  pdg1 = atoi(argv[1]);
  dFlowPart1 = new FemtoFlowDataBase(TMath::Abs(pdg1));
  dFlowPart1->SetExperimentName("ALICE");
  dFlowPart1->SetEnergy(2760);
  dFlowPart1->SetTableName("PbPb");
  dFlowPart1->SetCentrality(centrality);
  dFlowPart1->ShowParams();
  int foundVnForPart1 = dFlowPart1->DownloadGraphs();
  std::cout<<"\tFor particle "<<pdg1<<" found "<<foundVnForPart1<<" v params"<<std::endl;
  //geting pT distribution histogram
  fin1 = new TFile(argv[2]);
  fin1->GetObject("hpt1",h_pt1);
  std::cout<<"\tThere are "<<h_pt1->GetEntries()<<" entries in the histogram"<<std::endl;
  

//setting up second particle
  //database:
  pdg2 = atoi(argv[3]);
  dFlowPart2 = new FemtoFlowDataBase(TMath::Abs(pdg2));
  dFlowPart2->SetCentrality(centrality);
  dFlowPart2->ShowParams();
  int foundVnForPart2 = dFlowPart2->DownloadGraphs();
  std::cout<<"\tFor particle "<<pdg2<<" found "<<foundVnForPart2<<" v params"<<std::endl;

  fin2 = new TFile(argv[4]);
  fin2->GetObject("hpt2",h_pt2);
  std::cout<<"\tThere are "<<h_pt2->GetEntries()<<" entries in the histogram"<<std::endl;
  
  int N = atoi(argv[5]);
  std::cout<<"N="<<N<<std::endl;

  fCylm = new CorrFctnDirectYlm("Cylm",1,100,0.0,0.5);

  ROOT::Math::PtEtaPhiE4D<double> v1(1, 1, 1, 1);
  ROOT::Math::PtEtaPhiE4D<double> v2(1, 1, 1, 1);
  ROOT::Math::PtEtaPhiE4D<double> vsum(1, 1, 1, 1);

  TLorentzVector v1b;
  TLorentzVector v2b;
  TLorentzVector vsumb;
  TStopwatch* timer = new TStopwatch();
  timer->Start();
  std::cout << "Loop for filling the random" << std::endl;
  for (int i = 0; i < N; i++)
  {
    
    Double_t pt1 = h_pt1->GetRandom();
    Double_t pt2 = h_pt2->GetRandom();
    while(pt1 < 0.19 || pt1 > 1.5) pt1 = h_pt1->GetRandom();
    while(pt2 < 0.19 || pt2 > 1.5) pt2 = h_pt2->GetRandom();

    Double_t phi1 = dFlowPart1->GetPhi(pt1);
    Double_t phi2 = dFlowPart2->GetPhi(pt2);

    Double_t eta1 = 2. * (gRandom->Rndm() - 0.5) * 0.8;
    Double_t eta2 = 2. * (gRandom->Rndm() - 0.5) * 0.8;

    v1.SetCoordinates(pt1, eta1, phi1, massPi);
    v2.SetCoordinates(pt2, eta2, phi2, massK);

    v1b.SetPxPyPzE(v1.Px(), v1.Py(), v1.Pz(), sqrt(massPi * massPi + v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz()));
    v2b.SetPxPyPzE(v2.Px(), v2.Py(), v2.Pz(), sqrt(massK * massK + v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz()));

    // vsumb = v1b + v2b;

    // TVector3 boost;// = vsumb.BoostVector();
    // boost.SetXYZ(0.0,0.0,(v1b.Pz()+v2b.Pz())/(v1b.E()+v2b.E()));
    // v1b.Boost(-boost);
    // v2b.Boost(-boost);

    // std::cout<<"v1.Px="<<v1.Px()<<" v1.Py="<<v1.Py()<<" v2.Px="<<v2.Px()<<" v2.Py="<<v2.Py()<<std::endl;
    // std::cout<<boost.X()<<". "<<boost.Y()<<std::endl;
    // std::cout<<"v1b.Px="<<v1b.Px()<<" v1b.Py="<<v1b.Py()<<" v1b.Px="<<v2b.Px()<<" v2b.Py="<<v2b.Py()<<std::endl;
    // std::cout<<"\tv1b.Px+v2b.Px="<<v1b.Px()+v2b.Px()<<" v1b.Py+v2b.Py="<<v1b.Py()+v2b.Py()<<std::endl;

    // if(isnan(v1b.Px()) || isnan(v2b.Px())){
    //   std::cout<<"\tThere is a nan"<<std::endl;
    // }
    // Float_t kt = 0.5 * TMath::Hypot(v1b.Px() + v2b.Px(), v1b.Py() + v2b.Py());

    // Float_t pXsum = v1b.Px() + v2b.Px();
    // Float_t pYsum = v2b.Py() + v2b.Py();
    // Float_t pXdif = v1b.Px() - v2b.Px();
    // Float_t pYdif = v1b.Py() - v2b.Py();

    // Float_t mko = 0.5 * (pXsum * pXdif + pYsum * pYdif) / kt;
    // Float_t mks = (v1b.Px() * v2b.Py() - v1b.Py() * v2b.Px()) / kt;
    // Float_t mkl = v1b.Pz() - v2b.Pz();

    //std::cout<<mko<<","<<mks<<","<<mkl<<"("<<TMath::Sqrt(mko*mko + mks*mks + mkl*mkl)<<")"<<std::endl;

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
      
    fCylm->AddRealPair(mko,mks,mkl,1.0);

    /***** creating den histogram **************/
    phi1 = gRandom->Rndm()*TMath::TwoPi() - TMath::Pi();
    phi2 = gRandom->Rndm()*TMath::TwoPi() - TMath::Pi();

    v1.SetCoordinates(pt1, eta1, phi1, massPi);
    v2.SetCoordinates(pt2, eta2, phi2, massK);

    v1b.SetPxPyPzE(v1.Px(), v1.Py(), v1.Pz(), sqrt(massPi * massPi + v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz()));
    v2b.SetPxPyPzE(v2.Px(), v2.Py(), v2.Pz(), sqrt(massK * massK + v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz()));

    // vsumb = v1b + v2b;

    // boost = vsumb.BoostVector();
    // boost.SetXYZ(0.0,0.0,(v1b.Pz()+v2b.Pz())/(v1b.E()+v2b.E()));
    // v1b.Boost(-boost);
    // v2b.Boost(-boost);
    // kt = 0.5 * TMath::Hypot(v1b.X() + v2b.X(), v1b.Y() + v2b.Y());

    // pXsum = v1b.X() + v2b.X();
    // pYsum = v2b.Y() + v2b.Y();
    // pXdif = v1b.X() - v2b.X();
    // pYdif = v1b.Y() - v2b.Y();

    // mko = 0.5 * (pXsum * pXdif + pYsum * pYdif) / kt;
    // mks = (v1b.X() * v2b.Y() - v1b.Y() * v2b.X()) / kt;
    // mkl = v1b.Z() - v2b.Z();

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
          
    fCylm->AddMixedPair(mko,mks,mkl,1.0);

    if((i+1)%1000 == 0){
      printf("\r%6d/%d",i+1,N);
      fflush(stdout);
    }
    //printProgress(i+1,N);
  }
  timer->Stop();
  std::cout<<std::endl;
  std::cout<<"time: "<<timer->RealTime()<<"; CPU: "<<timer->CpuTime()<<std::endl;
  fCylm->Finish();
  fCylm->AddToList(fOutput);


  // std::cout<<"\nAdded fOutput to fCylm ("<<fOutput->GetEntries()<<","<<fOutput->GetSize()<<")"<<std::endl;
  // std::cout<<std::endl;
  TFile* fout = new TFile("output_cylm.root","RECREATE");
  fCylm->Write();
  fout->Close();
  
  std::cout<<"Saving my list"<<std::endl;
  fout = new TFile("output_list.root","RECREATE");
  fout->cd();
  fOutput->Write("fOutput",1);
  // fout->ls();
  fout->Close();

  //  std::cout<<"Saving cf (spherical harmonics)"<<std::endl;
  //  fout = new TFile("cf_hists.root","RECREATE");
  //  fout->cd();
  //  for(int il=0;il<4;il++){
  //    for(int im=0;im<=il;im++){
  //      fCylm->GetCfnRealHist(il,im)->Write(Form("CfReYlm%d%dCylm",il,im));
  //      fCylm->GetCfnImagHist(il,im)->Write(Form("CfImYlm%d%dCylm",il,im));
  //    }
  //  }      
  //  fout->Close();
  
  fOutput->Clear();

  std::cout<<"Saving database statistics"<<std::endl;
  fout = new TFile("database_stats.root","RECREATE");
  dFlowPart1->GetStatsHisto()->Write("hStats_pi");
  dFlowPart2->GetStatsHisto()->Write("hStats_K");
  fout->Close();

  //fCylm->Wrie();

  // std::cout<<"its time to delete"<<std::endl;
  // std::cout<<"\tfout"<<std::endl;
  delete fout;
  // std::cout<<"\th_pt1"<<std::endl;
  delete h_pt1;
  // std::cout<<"\th_pt2"<<std::endl;
  delete h_pt2;
  // std::cout<<"\tfin1"<<std::endl;
  delete fin1;
  // std::cout<<"\tfin2"<<std::endl;
  delete fin2;
  // std::cout<<"\tdFlowPart1"<<std::endl;
  //delete dFlowPart1;
  //std::cout<<"\tdFlowPart2"<<std::endl;
  //delete dFlowPart2;
  //std::cout<<"\tfCylm"<<std::endl;
  //delete fCylm;
  std::cout<<"done"<<std::endl;



  return 0;
}
