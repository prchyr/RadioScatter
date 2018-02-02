#include "RadioScatterEvent.hh"


ClassImp(RadioScatterEvent)

RadioScatterEvent::RadioScatterEvent(){
  dummy=-1;
}


double RadioScatterEvent::power(int txindex, int rxindex){
  double val;
  int entries=eventHist[txindex][rxindex]->GetNbinsX();
  for(int i=0;i<entries;i++){
    val+=eventHist[txindex][rxindex]->GetBinContent(i)*eventHist[txindex][rxindex]->GetBinContent(i);
  }
  return val;
}

double RadioScatterEvent::peakV(int txindex, int rxindex){
  return eventHist[txindex][rxindex]->GetMaximum();
}

double RadioScatterEvent::primaryParticleEnergy(){
  return primaryEnergy*nPrimaries*1000000.;//g4 primaries in MeV units
}
//simple thing-did any receiver trigger?
int RadioScatterEvent::triggered(double thresh){
  int trig=0;
  for(int i=0;i<ntx;i++){
    for(int j=0;j<nrx;j++){
      trig=peakV(i,j)>=thresh?1:0;
      if(trig==1)return 1;
    }
  }
  return trig;
}

double RadioScatterEvent::rms(int txindex, int rxindex){
  double val=0;
  int entries = eventHist[txindex][rxindex]->GetNbinsX();
  //  double * yy = eventHist[txindex][rxindex]->GetY();
  for(int i=0;i<entries;i++){
    double yy= eventHist[txindex][rxindex]->GetBinContent(i);
    val+=(yy*yy)/entries;
  }
  return sqrt(val);
}

double RadioScatterEvent::startFreq(){
  double freq=0;
  return freq;
}

double RadioScatterEvent::stopFreq(){
  double freq=0;
  return freq;
}

double RadioScatterEvent::bandWidth(){
  double width=0;
  return width;
}

double RadioScatterEvent::chirpSlope(){
  double slope=0;
  return slope;
}

double RadioScatterEvent::pathLength(int txindex, int rxindex){
  return((tx[txindex]-position).mag()+(rx[rxindex]-position).mag());
}
double RadioScatterEvent::duration(int txindex, int rxindex){
  TGraph og= getComplexEnvelope(txindex,rxindex,100);
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

TGraph RadioScatterEvent::getComplexEnvelope(int txindex, int rxindex,double cutoff){
  vector<double> xx, yy;
  int entries=reHist[txindex][rxindex]->GetNbinsX();
  for(int i=0;i<entries;i++){
    xx.push_back(reHist[txindex][rxindex]->GetBinCenter(i));
    yy.push_back(sqrt(pow(reHist[txindex][rxindex]->GetBinContent(i), 2)+pow(imHist[txindex][rxindex]->GetBinContent(i), 2)));
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

TH1F* RadioScatterEvent::getSpectrum(int txindex, int rxindex,bool dbflag){
  //  int size = eventGraph[txindex][rxindex]->GetN();

  int nbins = eventHist[txindex][rxindex]->GetNbinsX();
  int size=nbins;
  double samprate = nbins/(eventHist[txindex][rxindex]->GetXaxis()->GetXmax()-eventHist[txindex][rxindex]->GetXaxis()->GetXmin());//in samples/ns
  Float_t band = samprate;//in ghz
  Float_t bandwidth = band/2.;//nyquist
  Float_t timebase = size*(1./samprate);//ns

  TH1 *out = 0;
  out = eventHist[txindex][rxindex]->FFT(out, "mag");

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

double RadioScatterEvent::peakFreq(int txindex, int rxindex){
  TH1F * hist = getSpectrum(txindex, rxindex);
  double val = hist->GetBinCenter(hist->GetMaximumBin());
  //  hist->Delete();
  return val;
}

void RadioScatterEvent::spectrogram(int txindex, int rxindex,Int_t binsize, Int_t overlap){
  Int_t size = eventHist[txindex][rxindex]->GetNbinsX();
  Float_t xmax = eventHist[txindex][rxindex]->GetXaxis()->GetXmax();
  Float_t xmin = eventHist[txindex][rxindex]->GetXaxis()->GetXmin();

  Int_t nbins = size/overlap;
  char*timebuff;
  Float_t samplerate = size/(xmax-xmin);
  Float_t bandwidth = 1e9*samplerate;
  TH1F *in = new TH1F("inhi", "inhi", binsize, 0, binsize);  
  TH1*outt=0;
  //TH2F *outhist=new TH2F("outhist", "spectrogram", nbins, xmin, xmax, (binsize), 0, samplerate);
  //  cout<<binsize<<" "<<samplerate<<endl;
  spectrogramHist->SetBins(nbins, xmin, xmax, (binsize), 0, samplerate);
  Int_t start = 0;
  
  for(int i=0;i<=nbins;i++){
    for(int j = 0;j<=binsize;j++){
      if((j+start)>=size)break;
      in->SetBinContent(j, eventHist[txindex][rxindex]->GetBinContent(j+start));
   
    }

    outt=in->FFT(outt, "mag");
    outt->Scale(1./sqrt(binsize));
    for(int j = 0;j<=(binsize);j++){
      Double_t y = outt->GetBinContent(j);
      y = (10.*log10(pow(y, 2.)/50.));//mv->dbm
	//y=10.*log10(y)+30;
      spectrogramHist->SetBinContent(i,j,(y-(10.*log10(bandwidth))));//dmb/hz
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


int RadioScatterEvent::plotEvent(int txindex, int rxindex, int show_geom, int bins, int overlap){
  // TCanvas *c=0;
  // TCanvas *openc = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
  // if(!openc){
  //   c = new TCanvas("plotEvent","plotEvent", 800, 400);
  // }
  // else{
  //   cout<<openc->GetName()<<endl;
  //   c=openc;
  // }
  TCanvas *c = new TCanvas("plotEvent", "plotEvent", 800, 400);
  if(show_geom==0){
    c->Divide(1, 0);
    c->GetPad(1)->Divide(2, 0);
  }
  else{
    c->Divide(1, 2);
    c->GetPad(1)->Divide(2, 0);

    c->SetWindowSize(700, 700);
    
  }

  c->cd(1)->cd(1)->SetGrid();
  Float_t bandwidth = 1e9*sampleRate;//bandwith in Hz
  Float_t kB = 1.831e-23;
  Float_t thermal_noise = sqrt(kB*300.*50.*bandwidth)*1000.;//thermal noise (mV)
  int nbins = eventHist[txindex][rxindex]->GetNbinsX();
  //  cout<<bandwidth<<thermal_noise<<endl<<setprecision(10)<<nbins;
  Float_t noise, r1, r2, val;
  int j=0;
  for(int i=0;i<nbins;i++){
    j=i-(sampleRate*20);
    ran->Rannor(r1, r2);
    //    if(noise_flag==1){
      noise = r1*thermal_noise;
      //ev->eventHist->Fill(i, noise);
      eventHist[txindex][rxindex]->AddBinContent(i, noise);
      //    val=ev->eventHist->GetBinContent(i);
      //}
  }
  eventHist[txindex][rxindex]->GetXaxis()->SetTitle("ns");
  eventHist[txindex][rxindex]->GetYaxis()->SetTitle("mV");
  eventHist[txindex][rxindex]->Draw("histl");
  eventHist[txindex][rxindex]->SetStats(0);
  //timme->Draw();
  c->cd(1)->cd(2);
  // fft = dofft(eventHist);
  // fft->GetXaxis()->SetRangeUser(0, fft->GetXaxis()->GetXmax()/2);
  // fft->SetTitle("psd");
  // fft->GetYaxis()->SetTitleOffset(1.3);
  // fft->GetXaxis()->SetTitle("freq (GHz)");
  // fft->GetYaxis()->SetTitle("power (dBm/Hz)");
  // fft->Draw();
  // c->cd(3);
  //  g.plot(vals, "with lines");
  spectrogram(txindex, rxindex, bins, overlap);
  //  TImage *img =TImage::Create();
  //  img->FromPad(c);
  //img->WriteImage("54mhz_100TeV_proton.png");

  if(show_geom>0){
    c->cd(2);
    TH3F *rxhist = new TH3F("rxhist", "rxhist", 101, 1, -1, 101, 1, -1, 101, 1, -1);
    
    for(int i=0;i<nrx;i++){
      if(i==rxindex)continue;
      rxhist->Fill(1.+rx[i].z()/1000., 1.+rx[i].x()/1000., 1.+rx[i].y()/1000., 1.);
    }
    rxhist->BufferEmpty();

    TH3F *txhist = new TH3F("txhist", "txhist", 101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    TH3F *vertexhist = new TH3F("vertexhist", "vertexhist", 101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    TH3F *triggeredhist = new TH3F("triggeredhist", "triggeredhist", 101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    
    for(int i=0;i<ntx;i++){
      txhist->Fill(1.+tx[i].z()/1000., 1.+tx[i].x()/1000., 1.+tx[i].y()/1000., 1.);
      //      cout<<tx[i].z()<<" "<<tx[i].x()<<" "<<tx[i].y()<<endl;    
    }
    triggeredhist->Fill(1.+rx[rxindex].z()/1000., 1.+rx[rxindex].x()/1000., 1.+rx[rxindex].y()/1000., 1.);
    vertexhist->Fill(position.z()/1000., position.x()/1000., position.y()/1000., 1.);

    txhist->SetMarkerStyle(3);
    txhist->SetMarkerColor(kRed);
    rxhist->SetMarkerStyle(8);
    rxhist->SetMarkerColor(kBlack);
    rxhist->SetTitle("geometry");
    vertexhist->SetMarkerStyle(34);
    vertexhist->SetMarkerColor(kViolet);
    txhist->SetMarkerSize(2);
    rxhist->SetMarkerSize(2);
    vertexhist->SetMarkerSize(2);
    triggeredhist->SetMarkerColor(kGreen);
    triggeredhist->SetMarkerStyle(8);
    triggeredhist->SetMarkerSize(2.5);
    
    TPolyLine3D *line = new TPolyLine3D(2);
    double scale = rxhist->GetXaxis()->GetXmax()/4.;
    line->SetPoint(0,position.z()/1000., position.x()/1000., position.y()/1000.);
    line->SetPoint(1,(position.z()/1000.)+((direction.z())*scale), (position.x()/1000.)+((direction.x())*scale), (position.y()/1000.)+((direction.y())*scale));
    line->SetLineWidth(2);
    line->SetLineColor(kViolet);
    rxhist->Draw("p");
    txhist->Draw("p same");
    vertexhist->Draw("p same");
    triggeredhist->Draw("p same");
    line->Draw("same");

    TLegend *leg = new TLegend(.7,.7,.9,.9);
    leg->AddEntry(rxhist, "receivers", "p");
    leg->AddEntry(txhist, "transmitter", "p");
    leg->AddEntry(vertexhist, "vertex", "p");
    leg->AddEntry(triggeredhist, "this rx", "p");
    leg->Draw();
    //    c->Update();
  }
  
}
