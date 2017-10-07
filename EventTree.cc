#include "EventTree.hh"


ClassImp(EventTree)

EventTree::EventTree(){
dummy=-1;
}

double EventTree::power(){
double val;
for(int i=0;i<eventHist->GetNbinsX();i++){
val+=eventHist->GetBinContent(i)*eventHist->GetBinContent(i);
}
return val;
}

