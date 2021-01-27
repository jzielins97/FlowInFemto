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

  TH1D* h_pt1=0;
  TH1D* h_pt2=0;

  CorrFctnDirectYlm* fCylm;
  FemtoFlowDataBase* dFlowPart1;
  FemtoFlowDataBase* dFlowPart2;

  Double_t massPi = 0.139570;
  Double_t massK = 0.493677;

  if(argc < 6){
    std::cout<<"Error: not enough arguments. Needs: 4, was given: "<<argc-1<<"."<<std::endl;
    std::cout<<"<particle #1 pdg> <pT distribution file #1> <particle #2 pdg> <pT distribution #2> <N>"<<std::endl;
    return -1;
  }

//setting up first particle
  //database:
  pdg1 = atoi(argv[1]);
  dFlowPart1 = new FemtoFlowDataBase(pdg1);
  int foundVnForPart1 = dFlowPart1->downloadGraphs();
  std::cout<<"\tFor particle "<<pdg1<<" found "<<foundVnForPart1<<" v params"<<std::endl;
  //geting pT distribution histogram
  fin1 = new TFile(argv[2]);
  fin1->GetObject("hpt1",h_pt1);
  std::cout<<"\tThere are "<<h_pt1->GetEntries()<<" entries in the histogram"<<std::endl;

//setting up second particle
  //database:
  pdg2 = atoi(argv[3]);
  dFlowPart2 = new FemtoFlowDataBase(pdg2);
  int foundVnForPart2 = dFlowPart2->downloadGraphs();
  std::cout<<"\tFor particle "<<pdg2<<" found "<<foundVnForPart2<<" v params"<<std::endl;

  fin2 = new TFile(argv[4]);
  fin2->GetObject("hpt2",h_pt2);
  std::cout<<"\tThere are "<<h_pt2->GetEntries()<<" entries in the histogram"<<std::endl;

  int N = atoi(argv[5]);
  std::cout<<"\tN="<<N<<std::endl;

  fCylm = new CorrFctnDirectYlm("Cylm",3,100,0.0,0.5);

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

      Double_t phi1 = dFlowPart1->getPhi(pt1);
      Double_t phi2 = dFlowPart2->getPhi(pt2);

      Double_t eta1 = 2. * (gRandom->Rndm() - 0.5) * 0.8;
      Double_t eta2 = 2. * (gRandom->Rndm() - 0.5) * 0.8;

      v1.SetCoordinates(pt1, eta1, phi1, massPi);
      v2.SetCoordinates(pt2, eta2, phi2, massK);

      v1b.SetPxPyPzE(v1.Px(), v1.Py(), v1.Pz(), sqrt(massPi * massPi + v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz()));
      v2b.SetPxPyPzE(v2.Px(), v2.Py(), v2.Pz(), sqrt(massK * massK + v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz()));

      vsumb = v1b + v2b;

      TVector3 boost;// = vsumb.BoostVector();
      boost.SetXYZ(0.0,0.0,(v1b.Pz()+v2b.Pz())/(v1b.E()+v2b.E()));
      v1b.Boost(-boost);
      v2b.Boost(-boost);

      // std::cout<<"v1.Px="<<v1.Px()<<" v1.Py="<<v1.Py()<<" v2.Px="<<v2.Px()<<" v2.Py="<<v2.Py()<<std::endl;
      // std::cout<<boost.X()<<". "<<boost.Y()<<std::endl;
      // std::cout<<"v1b.Px="<<v1b.Px()<<" v1b.Py="<<v1b.Py()<<" v1b.Px="<<v2b.Px()<<" v2b.Py="<<v2b.Py()<<std::endl;
      // std::cout<<"\tv1b.Px+v2b.Px="<<v1b.Px()+v2b.Px()<<" v1b.Py+v2b.Py="<<v1b.Py()+v2b.Py()<<std::endl;

      if(isnan(v1b.Px()) || isnan(v2b.Px())){
        std::cout<<"\tThere is a nan"<<std::endl;
      }
      Float_t kt = 0.5 * TMath::Hypot(v1b.Px() + v2b.Px(), v1b.Py() + v2b.Py());

      Float_t pXsum = v1b.Px() + v2b.Px();
      Float_t pYsum = v2b.Py() + v2b.Py();
      Float_t pXdif = v1b.Px() - v2b.Px();
      Float_t pYdif = v1b.Py() - v2b.Py();

      Float_t mko = 0.5 * (pXsum * pXdif + pYsum * pYdif) / kt;
      Float_t mks = (v1b.Px() * v2b.Py() - v1b.Py() * v2b.Px()) / kt;
      Float_t mkl = v1b.Pz() - v2b.Pz();

      //std::cout<<mko<<","<<mks<<","<<mkl<<"("<<TMath::Sqrt(mko*mko + mks*mks + mkl*mkl)<<")"<<std::endl;
      fCylm->AddRealPair(mko,mks,mkl,1.0);

/***** creating den histogram **************/
      phi1 = gRandom->Rndm()*TMath::TwoPi();
      phi2 = gRandom->Rndm()*TMath::TwoPi();


      v1.SetCoordinates(pt1, eta1, phi1, massPi);
      v2.SetCoordinates(pt2, eta2, phi2, massK);

      v1b.SetPxPyPzE(v1.Px(), v1.Py(), v1.Pz(), sqrt(massPi * massPi + v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz()));
      v2b.SetPxPyPzE(v2.Px(), v2.Py(), v2.Pz(), sqrt(massK * massK + v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz()));

      vsumb = v1b + v2b;

      boost.SetXYZ(0.0,0.0,(v1b.Pz()+v2b.Pz())/(v1b.E()+v2b.E()));
      v1b.Boost(-boost);
      v2b.Boost(-boost);
      kt = 0.5 * TMath::Hypot(v1b.X() + v2b.X(), v1b.Y() + v2b.Y());

      pXsum = v1b.X() + v2b.X();
      pYsum = v2b.Y() + v2b.Y();
      pXdif = v1b.X() - v2b.X();
      pYdif = v1b.Y() - v2b.Y();

      mko = 0.5 * (pXsum * pXdif + pYsum * pYdif) / kt;
      mks = (v1b.X() * v2b.Y() - v1b.Y() * v2b.X()) / kt;
      mkl = v1b.Z() - v2b.Z();

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

  std::cout<<"Saving my list"<<std::endl;
  TFile* fout = new TFile("output_list.root","RECREATE");
  fout->cd();
  fOutput->Write("fOutput",1);
  // fout->ls();

  fout->Close();
  fOutput->Clear();

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
