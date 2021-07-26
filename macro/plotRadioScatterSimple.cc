/*
This macro makes a simple plot from a RadioScatterEvent. 

 */

void plotRadioScatterSimple(TString infile, int antenna=0, int entry=0){
  //open the file (provide absolute path, or relative to current directory)
  auto ff=TFile::Open(infile);
  //get the tree (always named "tree" in radioscatter files)
  //this tree contains the RadioScatterEvent object
  auto tree=(TTree*)ff->Get("tree");
  //create a new radioscatterevent object
  auto event=new RadioScatterEvent();
  //point the tree to the 
  tree->SetBranchAddress("event", &event);
  //get whatever entry you specified (first entry is default)
  tree->GetEntry(entry);
  //plot it.
  event->plotEvent(0,antenna, 0, 0, 120, 100);
}
