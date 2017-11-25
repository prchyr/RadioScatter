#include "RadioScatterEvent.hh"


ClassImp(RadioScatterEvent)

RadioScatterEvent::RadioScatterEvent(){
  dummy=-1;
}

double RadioScatterEvent::power(){
  double val;
  int entries=eventHist->GetNbinsX();
  for(int i=0;i<entries;i++){
    val+=eventHist->GetBinContent(i)*eventHist->GetBinContent(i);
  }
  return val;
}

double RadioScatterEvent::peakV(){
  return eventHist->GetMaximum();
}

int RadioScatterEvent::triggered(double thresh){
  int trig=peakV()>thresh?1:0;
  return trig;
}

double RadioScatterEvent::rms(){
  double val;
  int entries = eventGraph->GetN();
  double * yy = eventGraph->GetY();
  for(int i=0;i<entries;i++){
    val+=(yy[i]*yy[i])/entries;
  }
  return sqrt(val);
}

double RadioScatterEvent::startFreq(){
  double freq;
  return freq;
}

double RadioScatterEvent::stopFreq(){
  double freq;
  return freq;
}

double RadioScatterEvent::bandWidth(){
  double width;
  return width;
}

double RadioScatterEvent::chirpSlope(){
  double slope;
  return slope;
}

double RadioScatterEvent::duration(){
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

TGraph RadioScatterEvent::getComplexEnvelope(double cutoff){
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

TH1F* RadioScatterEvent::getSpectrum(bool dbflag){
  int size = eventGraph->GetN();
  int nbins = eventHist->GetNbinsX();
  double samprate = nbins/(eventHist->GetXaxis()->GetXmax()-eventHist->GetXaxis()->GetXmin());//in samples/ns
  Float_t band = samprate;//in ghz
  Float_t bandwidth = band/2.;//nyquist
  Float_t timebase = size*(1./samprate);//ns

  TH1 *out = 0;
  out = eventHist->FFT(out, "mag");

  spectrumHist->SetBins(nbins/2,0,bandwidth);

  for (Int_t i=0;i<=nbins/2;i++) {
    Double_t y = out->GetBinContent(i);
     if(dbflag==true){
       y = (10.*log10(pow(y, 2.)/50.))+30.;//normalized mv->dbW->dbm
     }
     //y=10.*log10(y)+30;                                                
  //   xx.push_back(i*(timebase/size));
  //   yy.push_back(y);
     //spectrumHist->SetBinContent(i, y-80-(10.*log10(band)));//dmb/hz with 80 db of gain             
     spectrumHist->SetBinContent(i, y);                                  
  }

  out->Delete();

  return spectrumHist;

}

double RadioScatterEvent::peakFreq(){
  TH1F * hist = getSpectrum();
  double val = hist->GetBinCenter(hist->GetMaximumBin());
  //  hist->Delete();
  return val;
}

void RadioScatterEvent::spectrogram(Int_t binsize, Int_t overlap){
  Int_t size = eventHist->GetNbinsX();
  Float_t xmax = eventHist->GetXaxis()->GetXmax();
  Float_t xmin = eventHist->GetXaxis()->GetXmin();

  Int_t nbins = size/overlap;
  char*timebuff;
  Float_t samplerate = size/(xmax-xmin);
  Float_t bandwidth = 1e9*samplerate;
  TH1F *in = new TH1F("inhi", "inhi", binsize, 0, binsize);  
  TH1*outt=0;
  //TH2F *outhist=new TH2F("outhist", "spectrogram", nbins, xmin, xmax, (binsize), 0, samplerate);
  cout<<binsize<<" "<<samplerate<<endl;
  spectrogramHist->SetBins(nbins, xmin, xmax, (binsize), 0, samplerate);
  Int_t start = 0;
  
  for(int i=0;i<=nbins;i++){
    for(int j = 0;j<=binsize;j++){
      if((j+start)>=size)break;
      in->SetBinContent(j, eventHist->GetBinContent(j+start));
   
    }

    outt=in->FFT(outt, "mag");
    outt->Scale(1./sqrt(binsize));
    for(int j = 0;j<=(binsize);j++){
      Double_t y = outt->GetBinContent(j);
      y = (10.*log10(pow(y, 2.)/50.));//mv->dbm
	//y=10.*log10(y)+30;
            spectrogramHist->SetBinContent(i,j,(y-(10.*log10(bandwidth/binsize))));//dmb/hz
	    //	spectrogramHist->SetBinContent(i,j,y);//dmb

    }
    start+=overlap-1;
  }
  //cout<<"here"<<endl;
  gPad->SetRightMargin(.15);
  spectrogramHist->GetYaxis()->SetRangeUser(0, spectrogramHist->GetYaxis()->GetXmax()/2);
  spectrogramHist->GetXaxis()->SetTitle("Time (ns)");
  spectrogramHist->GetYaxis()->SetTitleOffset(1.3);
  spectrogramHist->GetYaxis()->SetTitle("Frequency (GHz)");
  spectrogramHist->GetZaxis()->SetTitle("dBm/Hz");
  spectrogramHist->GetZaxis()->SetTitleOffset(1.5);
  spectrogramHist->SetStats(0);
  spectrogramHist->Draw("colz");
  //maxbins->Draw("same");
  //maxvals->SetLineColor(kRed);
  //maxvals->Draw("same");
  outt->Delete();
  in->Delete();

  // spectrogramHist->Delete();
 
  //return outdat;
}


int RadioScatterEvent::plotEvent(int bins, int overlap){
  TCanvas *c = new TCanvas("plotEvent","plotEvent", 800, 400);
  c->Divide(2, 0);
  c->cd(1)->SetGrid();
  Float_t bandwidth = 1e9*sampleRate;//bandwith in Hz
  Float_t kB = 1.831e-23;
  Float_t thermal_noise = sqrt(kB*300.*50.*bandwidth)*1000.;//thermal noise (mV)
  int nbins = eventHist->GetNbinsX();
  //  cout<<bandwidth<<thermal_noise<<endl<<setprecision(10)<<nbins;
  Float_t noise, r1, r2, val;
  int j=0;
  for(int i=0;i<nbins;i++){
    j=i-(sampleRate*20);
    ran->Rannor(r1, r2);
    //    if(noise_flag==1){
      noise = r1*thermal_noise;
      //ev->eventHist->Fill(i, noise);
      eventHist->AddBinContent(i, noise);
      //    val=ev->eventHist->GetBinContent(i);
      //}
  }
  eventHist->GetXaxis()->SetTitle("ns");
  eventHist->GetYaxis()->SetTitle("mV");
  eventHist->Draw("histl");
  eventHist->SetStats(0);
  //timme->Draw();
  c->cd(2);
  // fft = dofft(eventHist);
  // fft->GetXaxis()->SetRangeUser(0, fft->GetXaxis()->GetXmax()/2);
  // fft->SetTitle("psd");
  // fft->GetYaxis()->SetTitleOffset(1.3);
  // fft->GetXaxis()->SetTitle("freq (GHz)");
  // fft->GetYaxis()->SetTitle("power (dBm/Hz)");
  // fft->Draw();
  // c->cd(3);
  //  g.plot(vals, "with lines");
  spectrogram(bins, overlap);
  //  TImage *img =TImage::Create();
  //  img->FromPad(c);
  //img->WriteImage("54mhz_100TeV_proton.png");

}
