int drawCF{
  int drawnFunctions = 0;
  TFile* fin = new TFile("output_list.root");
  TList* list = (TList*) fin->FindObject("fOutput");

  TH1D* hNum[];
  TH1D* hDen[];
  TH1D* hCF[];
  int l=0;
  int m=0;

}
