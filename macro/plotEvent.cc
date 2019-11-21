void plotEvent(TString infile, int entry, int rxIndex){
  auto ff=TFile::Open(infile);
  auto tree=(TTree*)ff->Get("tree");
  auto rs=new RadioScatterEvent();
  tree->SetBranchAddress("event", &rs);
  tree->GetEntry(entry);
  rs->plotEvent(0,rxIndex);
}
