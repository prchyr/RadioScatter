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

double EventTree::rms(){
  double val;
  int entries = eventGraph->GetN();
  double * yy = eventGraph->GetY();
  for(int i=0;i<entries;i++){
    val+=(yy[i]*yy[i])/entries;
  }
  return sqrt(val);
}

double EventTree::startFreq(){
  double freq;
  return freq;
}

double EventTree::stopFreq(){
  double freq;
  return freq;
}

double EventTree::bandWidth(){
  double width;
  return width;
}

double EventTree::chirpSlope(){
  double slope;
  return slope;
}

double EventTree::duration(){

  TGraph *og = getComplexEnvelope(100);
  double * xx=og->GetX();
  double * yy=og->GetY();
  int n=og->GetN();
  double lastval=0, val, time1=0, time2=0, avg=0, last_avg=0;
  //  std::deque<double> avg_vec(10, 0);
  double thresh=.01;
    double highthresh=.01, lowthresh=.0001;
  for(int i=10;i<n-10;i++){
    // avg_vec.pop_front();
    // avg_vec.push_back(yy[i]);
    // for(int j=0;j<10;j++){
    //   avg+=avg_vec[j]/10.;
    // }
    // //    cout<<avg<<endl;
    // if(avg>thresh&&last_avg<thresh){
    //   time1=xx[i];
    // }
    // if(avg<thresh&&last_avg>thresh){
    //   time2=xx[i];
    // }
    // last_avg=avg;
    // avg=0;
    if(yy[i]>highthresh&&lastval<=highthresh){
      time1=xx[i];
    }
    if(yy[i]<lowthresh&&lastval>=lowthresh){
      time2=xx[i];
    }
    lastval=yy[i];
  }
  //  cout<<time2<<" "<<time1<<endl;
  return time2-time1;

}

TGraph * EventTree::getComplexEnvelope(double cutoff){
  vector<double> xx, yy;
  int entries=reHist->GetNbinsX();
  for(int i=0;i<entries;i++){
    xx.push_back(reHist->GetBinCenter(i));
    yy.push_back(sqrt(pow(reHist->GetBinContent(i), 2)+pow(imHist->GetBinContent(i), 2)));
  }
  
  if(cutoff>0){
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
