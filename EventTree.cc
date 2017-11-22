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

double EventTree::peakV(){
  return eventHist->GetMaximum();
}

int EventTree::triggered(double thresh){
  int trig=peakV()>thresh?1:0;
  return trig;
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
  TGraph og= getComplexEnvelope(100);
  double * xx=og.GetX();
  double * yy=og.GetY();
  int n=og.GetN();
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

TGraph EventTree::getComplexEnvelope(double cutoff){
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
  
  
  TGraph og(xx.size(), &xx[0], &out[0]);

  return og;
  }
  else{
    TGraph og(xx.size(), &xx[0], &yy[0]);

    return og;
  }
}

TH1F* EventTree::getSpectrum(bool dbflag){
  int size = eventGraph->GetN();
  int nbins = eventHist->GetNbinsX();
  double samprate = nbins/(eventHist->GetXaxis()->GetXmax()-eventHist->GetXaxis()->GetXmin());//in samples/ns
  Float_t band = samprate;//in ghz
  Float_t bandwidth = band/2.;//nyquist
  Float_t timebase = size*(1./samprate);//ns
  //inhist = eventHist;
  //cout<<nbins<<" "<<samprate<<" "<<bandwidth<<" "<<timebase<<endl;
  TH1 *out = 0;
  out = eventHist->FFT(out, "mag");


  TH1F *outhist = new TH1F("freq","",nbins/2,0,bandwidth);
  // vector<double> xx, yy;
  for (Int_t i=0;i<=nbins/2;i++) {
    Double_t y = out->GetBinContent(i);
     if(dbflag==true){
       y = (10.*log10(pow(y, 2.)/50.))+30.;//normalized mv->dbW->dbm
     }
     //y=10.*log10(y)+30;                                                
  //   xx.push_back(i*(timebase/size));
  //   yy.push_back(y);
     //outhist->SetBinContent(i, y-80-(10.*log10(band)));//dmb/hz with 80 db of gain             
     outhist->SetBinContent(i, y);                                  
  }
  // //  double rebinval = 256;
  // //  outhist->Rebin(size/rebinval);//rebin to make it easier to read       
  // //  outhist->Scale(rebinval/size);
  // //cout<<outhist->Integral()<<endl;
  // outhist->GetXaxis()->SetRangeUser(0,(band/200000));
  // outhist->GetYaxis()->SetTitle("dBm/Hz");
  // //   outhist->GetYaxis()->SetTitle("normalized amplitude");      
  // outhist->GetXaxis()->SetTitle("freq (MHz)");
  // outhist->GetYaxis()->SetTitleOffset(1.5);
  // out->Delete();
  // inhist->Delete();
  // //outhist->Draw();
  //TGraph gr(xx.size(), &xx[0], &yy[0]);
  //TH1F*outt=0;
  //TH1F* outt = (TH1F*)outhist->Clone();
  //  outhist->Delete();
  out->Delete();
  return outhist;
  //  outhist->Delete();

}

double EventTree::peakFreq(){
  TH1F * hist = getSpectrum();
  double val = hist->GetBinCenter(hist->GetMaximumBin());
  hist->Delete();
  return val;
}
