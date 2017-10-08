#include "EventTree.hh"


ClassImp(EventTree)

EventTree::EventTree(){
  dummy=-1;
}

double EventTree::power(){
  double val;
  int entries=eventHist->GetNbinsX();
  for(int i=0;i<entries;i++){
    val+=eventHist->GetBinContent(i)*eventHist->GetBinContent(i);
  }
  return val;
}

// double EventTree::duration(){
//   double val, lastval, avg;
//   int entries=eventHist->GetNbinsX();
//   for(int i=0;i<entries;i++){
//     val=eventHist->GetBinContent(i);
//     if(lastval=0&&avg=0&&val!=0){
      
//   }
  
//   return duration;
// }

TGraph * EventTree::getEnvelope(double cutoff){
  vector<double> xx, yy;
  int entries=reHist->GetNbinsX();
  for(int i=0;i<entries;i++){
    xx.push_back(reHist->GetBinCenter(i));
    yy.push_back(sqrt(pow(reHist->GetBinContent(i), 2)+pow(imHist->GetBinContent(i), 2)));
  }
  
  if(cutoff){
    vector<double> out;  
    double w = cutoff*2.*pi*1e6;
    //    double T = 1/(sampleRate*1e9);
    double T = 1/(1e10);
    double a, b, c, value;

    a = w*T;
    b = exp(-w*T);
	
    //cout<<setprecision(12);
    //cout<<"filter coefficients"<<endl<<a<<endl<<b<<endl<<c<<endl;
    int size = yy.size();
    for(int i=0;i<size;i++){
      if(i>0){
	
	value = a*yy[i]+b*out[i-1];
	
	out.push_back(value);
      }
      if(i==0){
	value = a*yy[i];
	out.push_back(value);
      } 
    }
  
  
  TGraph *og=new TGraph(xx.size(), &xx[0], &out[0]);
  return og;
  }
  else{
    TGraph *og=new TGraph(xx.size(), &xx[0], &yy[0]);
  }
}
