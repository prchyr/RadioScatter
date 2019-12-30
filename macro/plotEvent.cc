void plotEvent(TString infile, int entry=0, int rxIndex=0, double noise=0., int geom=0){
  auto ff=TFile::Open(infile);
  auto tree=(TTree*)ff->Get("tree");
  auto rs=new RadioScatterEvent();
  tree->SetBranchAddress("event", &rs);
  tree->GetEntry(entry);
  rs->plotEvent(0,rxIndex, noise,geom);
}
