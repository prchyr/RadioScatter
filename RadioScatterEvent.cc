/*
this is radioscatter. copyright s. prohira 

released under GPL3.
 
 */
#include "RadioScatterEvent.hh"


ClassImp(RadioScatterEvent)

RadioScatterEvent::RadioScatterEvent(){

  dummy=-1;
}

int RadioScatterEvent::reset(){
  rx.clear();
  tx.clear();
  nPrimaries=0;
  totNScatterers=0;
  return 1;
}
//not memory safe.
TGraph * RadioScatterEvent::getGraph(int txindex, int rxindex){
  TGraph *outgr=new TGraph(eventHist[txindex][rxindex]);
  return outgr;
}
double RadioScatterEvent::integratedPower(int txindex, int rxindex){
  double val;
  int entries=eventHist[txindex][rxindex]->GetNbinsX();
  for(int i=0;i<entries;i++){
    val+=eventHist[txindex][rxindex]->GetBinContent(i)*eventHist[txindex][rxindex]->GetBinContent(i);
  }
  return val;
}

double RadioScatterEvent::integratedPower(int txindex, int rxindex, double tlow, double thigh, double dcoffset){
  double val;
  //  int entries=eventHist[txindex][rxindex]->GetNbinsX();
  int binlow=eventHist[txindex][rxindex]->FindBin(tlow);
  int binhigh=eventHist[txindex][rxindex]->FindBin(thigh);
  for(int i=binlow;i<binhigh;i++){
    val+=(dcoffset+eventHist[txindex][rxindex]->GetBinContent(i))*(dcoffset+eventHist[txindex][rxindex]->GetBinContent(i));
  }
  return val;
}

double RadioScatterEvent::integratedPowerAroundPeak(int txindex, int rxindex, double window){
  double val;
  int halfwindow=(int)(window/2.*sampleRate);
  TH1F *h1=eventHist[txindex][rxindex];
  int binlow=h1->GetMaximumBin()-halfwindow;
  int binhigh=h1->GetMaximumBin()+halfwindow;
  for(int i=binlow;i<binhigh;i++){
    val+=(h1->GetBinContent(i))*(h1->GetBinContent(i));
  }
  return val;
}

double RadioScatterEvent::integratedVoltage(int txindex, int rxindex){
  double val;
  int entries=eventHist[txindex][rxindex]->GetNbinsX();
  for(int i=0;i<entries;i++){
    val+=eventHist[txindex][rxindex]->GetBinContent(i);
  }
  return val;
}

double RadioScatterEvent::integratedVoltage(int txindex, int rxindex, double tlow, double thigh, double dcoffset){
  double val;
  //  int entries=eventHist[txindex][rxindex]->GetNbinsX();
  int binlow=eventHist[txindex][rxindex]->FindBin(tlow);
  int binhigh=eventHist[txindex][rxindex]->FindBin(thigh);
  for(int i=binlow;i<binhigh;i++){
    val+=(dcoffset+eventHist[txindex][rxindex]->GetBinContent(i));
  }
  return val;
}

double RadioScatterEvent::peakPowerMW(int txindex, int rxindex){
  double val=pow((peakV(txindex, rxindex)*1000), 2)/50.;

  return val;
}
double RadioScatterEvent::peakPowerW(int txindex, int rxindex){
  double val=pow((peakV(txindex, rxindex)), 2)/50.;

  return val;
}

double RadioScatterEvent::effectiveCrossSection(int txindex, int rxindex){
  double r_tx = (tx[txindex].vect()-position).mag()/1000.;//m
  double r_rx = (rx[rxindex].vect()-position).mag()/1000.;
  double lambda = 2.997e8/(freq*1.e9);
  cout<<r_rx<<" "<<r_tx<<" "<<lambda<<" "<<peakPowerW(txindex, rxindex)<<" txGain: "<<txGain<<" txpower: "<<txPowerW<<endl;
  return (pow(4.*pi, 3)*pow(r_tx, 2)*pow(r_rx, 2)*peakPowerW(txindex, rxindex))/(txPowerW*txGain*rxGain*pow(lambda, 2));
}


double RadioScatterEvent::thermalNoiseRMS(){
  double bandwidth = 1e9*sampleRate/2;//bandwith in Hz
  double kB = 1.831e-23;
  return  sqrt(kB*300.*50.*bandwidth);//thermal noise RMS (mV)
}

double RadioScatterEvent::peakV(int txindex, int rxindex){
  return reHist[txindex][rxindex]->GetMaximum();
}

double RadioScatterEvent::primaryParticleEnergy(){
  return primaryEnergy*nPrimaries*1000000.;//g4 primaries in MeV units
}
//simple thing-did any receiver trigger?
int RadioScatterEvent::triggered(double thresh, int n_antennas){
  int trig=0,num=0;
  for(int i=0;i<ntx;i++){
    for(int j=0;j<nrx;j++){
      trig=peakV(i,j)>=thresh?1:0;
      if(trig==1)num++;
      if(num>=n_antennas)return 1;
    }
  }
  return trig;
}
//how many triggered?
int RadioScatterEvent::nTriggered(double thresh){
  int trig=0,num=0;
  for(int i=0;i<ntx;i++){
    for(int j=0;j<nrx;j++){
      trig=peakV(i,j)>=thresh?1:0;
      if(trig==1)num++;
    }
  }
  return num;
}

int RadioScatterEvent::trigSingle(double thresh, int ant){
 return peakV(0,ant)>=thresh?1:0;
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

double RadioScatterEvent::pathLengthMM(int txindex, int rxindex){
  return((tx[txindex]-position).mag()+(rx[rxindex]-position).mag());
}
double RadioScatterEvent::pathLengthM(int txindex, int rxindex){
  return((tx[txindex]-position).mag()+(rx[rxindex]-position).mag())/1000.;
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

TH1F * RadioScatterEvent::getComplexEnvelope(int txindex, int rxindex,double cutoff){
  vector<double> xx, yy;
  int entries=reHist[txindex][rxindex]->GetNbinsX();
  ceHist->SetBins(entries, reHist[txindex][rxindex]->GetXaxis()->GetXmin(), reHist[txindex][rxindex]->GetXaxis()->GetXmax());
  

  for(int i=0;i<entries;i++){
    xx.push_back(reHist[txindex][rxindex]->GetBinCenter(i));
    yy.push_back(sqrt(pow(reHist[txindex][rxindex]->GetBinContent(i), 2)+pow(imHist[txindex][rxindex]->GetBinContent(i), 2)));
    ceHist->SetBinContent(i, sqrt(pow(reHist[txindex][rxindex]->GetBinContent(i)+imHist[txindex][rxindex]->GetBinContent(i), 2)));
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
	ceHist->SetBinContent(i, value);
	out.push_back(value);
      }
      if(i==0){
	value = a*yy[i];
	out.push_back(value);
	ceHist->SetBinContent(i, value);
      } 
    }
  
  
    //TGraph og(xx.size(), &xx[0], &out[0]);

  return ceHist;
  }
  else{
    //     TGraph og(xx.size(), &xx[0], &yy[0]);
     return ceHist;
  //   return og;
   }
  
}

TH1F* RadioScatterEvent::getSpectrum(int txindex, int rxindex,bool dbflag){
  //  int size = eventGraph[txindex][rxindex]->GetN();

  int nbins = reHist[txindex][rxindex]->GetNbinsX();
  int size=nbins;
  double samprate = nbins/(reHist[txindex][rxindex]->GetXaxis()->GetXmax()-reHist[txindex][rxindex]->GetXaxis()->GetXmin());//in samples/ns
  Float_t band = samprate*1.e9;//in hz
  Float_t bandwidth = samprate/2.;//nyquist
  Float_t timebase = size*(1./samprate);//ns

  TH1 *out = 0;
  out = reHist[txindex][rxindex]->FFT(out, "mag");

  spectrumHist->SetBins(nbins/2,0,bandwidth);

  for (Int_t i=0;i<=nbins/2;i++) {
    Double_t y = out->GetBinContent(i);
     if(dbflag==true){
       y = (10.*log10(pow(y, 2.)/50.))+30.;//normalized mv->dbW->dbm
       spectrumHist->SetBinContent(i, y-(10.*log10(band)));//dmb/hz            
     }
     //y=10.*log10(y)+30;                                                
  //   xx.push_back(i*(timebase/size));
  //   yy.push_back(y);
     else{
       spectrumHist->SetBinContent(i, y);
     }
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
  //  gPad->SetRightMargin(.15);
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

double RadioScatterEvent::sineSubtract(int txindex, int rxindex, double rangestart, double rangeend, double p0, double p1, double p2){
  TF1 *fun = new TF1("fun", "[0]*sin(2.*pi*[1]*x+[2])", rangestart, rangeend);
  //  TF1 *fun = new TF1("fun", "[0]*sin(2.*pi*.8*x+[1])", rangestart, rangeend);
  fun->SetParameters(p0, p1, p2);

  
  fun->SetParLimits(0, p0-(.1*p0),p0+(.1*p0));
  fun->SetParLimits(1, p1-(.1*p1), p1+(.1*p1));
  fun->SetParLimits(2, 0, 100.);

  eventHist[txindex][rxindex]->Fit(fun, "R");
  TF1 *fun2 = fun;
  fun2->SetParameter(0, -fun->GetParameter(0));
  fun2->SetParameter(1, fun->GetParameter(1));
  fun2->SetParameter(2, fun->GetParameter(2));
  //  fun->Draw();
  eventHist[txindex][rxindex]->Eval(fun2, "A");
 
  return fun->GetParameter(0)/abs(fun->GetParameter(0))*fun->GetParameter(2)/(fun->GetParameter(1)*2.*pi);
  //  return (fun->GetParameter(0)/abs(fun->GetParameter(0)))*fun->GetParameter(1)/(p2*2.*pi);
}

int RadioScatterEvent::backgroundSubtract(int txindex, int rxindex, TH1F *bSubHist){
  eventHist[txindex][rxindex]->Add(bSubHist, -1);
  return 1;
}

TH1F* RadioScatterEvent::makeBackgroundSubtractHist(int txindex, int rxindex, TString bfile){
  
  TFile *ff=TFile::Open(bfile);
  TTree *tree = (TTree*)ff->Get("tree");
  RadioScatterEvent *rse = new RadioScatterEvent();
  tree->SetBranchAddress("event", &rse);
  TH1F *bSubHist = new TH1F("bsub", "bsub", eventHist[txindex][rxindex]->GetNbinsX(), 0, eventHist[txindex][rxindex]->GetNbinsX()/sampleRate);
  int entries = tree->GetEntries();
  for(int i=0;i<entries;i++){
    tree->GetEntry(i);
    bSubHist->Add(rse->eventHist[txindex][rxindex]);
  }
  bSubHist->Scale(1./entries);
  //ff->Close();
  //  delete(ff);
  //delete(tree);
  return bSubHist;
}


int RadioScatterEvent::plotEvent(int txindex, int rxindex, int noise_flag, int show_geom, int bins, int overlap){
  // TCanvas *c=0;
  TSeqCollection *canlist = gROOT->GetListOfCanvases();
  TCanvas *openc = (TCanvas*)canlist->At(canlist->GetEntries()-1);
  int redraw_canvas = 0;
  if(canlist->GetEntries()==0||openc->GetCanvasImp()==NULL){
    redraw_canvas = 1;
    ccc = new TCanvas("plotEvent","plotEvent", 800, 400);
  }
  else{
     cout<<openc->GetName()<<endl;
     ccc=openc;
   }
  ccc->SetName("plotEvent");
  ccc->SetTitle("plotEvent");
  ///  TCanvas *c = new TCanvas("plotEvent", "plotEvent", 800, 400);
  if(show_geom==0&&redraw_canvas==1){
    ccc->Divide(1, 0);
    ccc->GetPad(1)->Divide(2, 0);
    ccc->SetWindowSize(800,400);
    cout<<"single event view"<<endl;
  }
  else if (show_geom>0&&redraw_canvas==1){
    //vertical canvas
    //    ccc->Divide(1, 2);
    //    ccc->GetPad(1)->SetPad(.005, .6525, .995, .995);
    //ccc->GetPad(1)->Divide(2, 0);
    //ccc->GetPad(2)->SetPad(.005, .005, .995, .6475);
    //ccc->SetWindowSize(500, 700);
    //horizontal canvas
    ccc->Divide(2, 0);
    ccc->GetPad(1)->SetPad(.005, .005, .3475, .995);
    ccc->GetPad(1)->Divide(0, 2);
    ccc->GetPad(2)->SetPad(.3525, .005, .995, .995);
    ccc->SetWindowSize(1200, 700);
    
  }

  ccc->cd(1)->cd(1)->SetGrid();
  Float_t bandwidth = 1e9*sampleRate/20.;//bandwith in Hz/10
  Float_t kB = 1.831e-23;
  Float_t thermal_noise = sqrt(kB*300.*50.*bandwidth);//thermal noise (V)
  int nbins = eventHist[txindex][rxindex]->GetNbinsX();
  //  cout<<bandwidth<<thermal_noise<<endl<<setprecision(10)<<nbins;
  Float_t noise, r1, r2, val;
  int j=0;
  for(int i=0;i<nbins;i++){
    j=i-(sampleRate*20);
    ran->Rannor(r1, r2);
    if(noise_flag==1){
      noise = r1*thermal_noise;
      //ev->eventHist->Fill(i, noise);
      eventHist[txindex][rxindex]->AddBinContent(i, noise);
      //    val=ev->eventHist->GetBinContent(i);
      }
    if(noise_flag>1){
      noise=r1*(double)noise_flag;
      eventHist[txindex][rxindex]->AddBinContent(i, noise);
      //    val=ev->eventHist->GetBinContent(i);
      }
    if(SINE_SUBTRACT==1){
      sineSubtract(txindex, rxindex);
    }
  }
  eventHist[txindex][rxindex]->GetXaxis()->SetTitle("Time (ns)");
  eventHist[txindex][rxindex]->GetYaxis()->SetTitle("V");
  eventHist[txindex][rxindex]->Draw("histl");
  eventHist[txindex][rxindex]->SetStats(0);
  //timme->Draw();
  ccc->cd(1)->cd(2);
  // fft = dofft(eventHist);
  // fft->GetXaxis()->SetRangeUser(0, fft->GetXaxis()->GetXmax()/2);
  // fft->SetTitle("psd");
  // fft->GetYaxis()->SetTitleOffset(1.3);
  // fft->GetXaxis()->SetTitle("freq (GHz)");
  // fft->GetYaxis()->SetTitle("power (dBm/Hz)");
  // fft->Draw();
  // ccc->cd(3);
  //  g.plot(vals, "with lines");
  spectrogram(txindex, rxindex, bins, overlap);
  //  TImage *img =TImage::Create();
  //  img->FromPad(c);
  //img->WriteImage("54mhz_100TeV_proton.png");

  if(show_geom>0){
    rxhist->Reset();
    txhist->Reset();
    vertexhist->Reset();
    triggeredhist->Reset();
    pointingHist->Reset();
    ccc->cd(2);
    //TH3F *rxhist = new TH3F("rxhist", "rxhist", 101, 1, -1, 101, 1, -1, 101, 1, -1);
    
    for(int i=0;i<nrx;i++){
      if(i==rxindex)continue;
      rxhist->Fill(1.+rx[i].z()/1000., 1.+rx[i].x()/1000., 1.+rx[i].y()/1000., eventHist[txindex][rxindex]->GetMaximum());
    }
    rxhist->BufferEmpty();
    RXHIST_FILLED=1;
  

    txhist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    vertexhist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    triggeredhist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    pointingHist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    
    for(int i=0;i<ntx;i++){
      txhist->Fill(1.+tx[i].z()/1000., 1.+tx[i].x()/1000., 1.+tx[i].y()/1000., 1.);
      //      cout<<tx[i].z()<<" "<<tx[i].x()<<" "<<tx[i].y()<<endl;    
    }
    triggeredhist->Fill(1.+rx[rxindex].z()/1000., 1.+rx[rxindex].x()/1000., 1.+rx[rxindex].y()/1000., 1.);
    vertexhist->Fill(1.+position.z()/1000., 1.+position.x()/1000., 1.+position.y()/1000., 1.);

    //see if we can reconstruct the event
    if(POINTING_MAP_BUILT==0){
      buildMap();
    }
    HepLorentzVector vvv = pointUsingMap();

    //pointingHist->Fill(1.+vvv.z()/1000., 1.+vvv.x()/1000., rxhist->GetZaxis()->GetXmax()+(vvv.y()/1000.), 1.);
    //    HepLorentzVector vvv = pointingTest();
    
    //    pointingHist->Fill(1.+vvv.z(), 1.+vvv.x(), vvv.y(), 1.);
    pointingHist->Fill(1.+vvv.z()/1000., 1.+vvv.x()/1000., (vvv.y()/1000.), 1.);

    txhist->SetMarkerStyle(3);
    txhist->SetMarkerColor(kRed);
    rxhist->SetMarkerStyle(8);
    rxhist->SetMarkerColor(kBlack);
    rxhist->SetTitle("geometry");
    rxhist->GetXaxis()->SetTitle("x (meters)");
    rxhist->GetXaxis()->SetTitleOffset(1.5);
    rxhist->GetYaxis()->SetTitle("y (meters)");
    rxhist->GetYaxis()->SetTitleOffset(1.5);
    rxhist->GetZaxis()->SetTitle("z (meters)");
    rxhist->GetZaxis()->SetTitleOffset(1.5);
    vertexhist->SetMarkerStyle(34);
    vertexhist->SetMarkerColor(kViolet);
    pointingHist->SetMarkerStyle(21);
    pointingHist->SetMarkerColor(kBlue);
    pointingHist->SetMarkerSize(1.5);
    txhist->SetMarkerSize(3);
    rxhist->SetMarkerSize(2);
    vertexhist->SetMarkerSize(2);
    triggeredhist->SetMarkerColor(kGreen);
    triggeredhist->SetMarkerStyle(8);
    triggeredhist->SetMarkerSize(2.5);

    TPolyLine3D *shower_indicator_line = new TPolyLine3D(2);

    double scale = rxhist->GetXaxis()->GetXmax()/4.;
    //cout<<scale<<endl;
    shower_indicator_line->SetPoint(0,1.+position.z()/1000., 1.+position.x()/1000., 1.+position.y()/1000.);

    shower_indicator_line->SetPoint(1,(1.+position.z()/1000.)+((direction.z())*scale), (1.+position.x()/1000.)+((direction.x())*scale), (1.+position.y()/1000.)+((direction.y())*scale));
    //cout<<"nere"<<endl;
    shower_indicator_line->SetLineWidth(3);
    shower_indicator_line->SetLineColor(kViolet);

    //cout<<"nereee"<<endl;
    rxhist->Draw("p");
    txhist->Draw("p same");
    vertexhist->Draw("p same");
    triggeredhist->Draw("p same");
    shower_indicator_line->Draw("same");
    pointingHist->Draw("p same");
    TLegend *leg = new TLegend(.7,.7,.9,.9);

    leg->AddEntry(txhist, "transmitter", "p");
    leg->AddEntry(vertexhist, "vertex", "p");
    leg->AddEntry(shower_indicator_line, "shower", "l");
    leg->AddEntry(rxhist, "receivers", "p");
    leg->AddEntry(triggeredhist, "this receiver", "p");
    leg->AddEntry(pointingHist, "source approx", "p");
    leg->Draw();
    //    ccc->Update();
  }
  
}

int RadioScatterEvent::plotEventNotebook(int txindex, int rxindex, int noise_flag, int show_geom, int bins, int overlap){

  TCanvas * ccc = new TCanvas();

  ccc->SetName("plotEvent");
  ccc->SetTitle("plotEvent");
  ///  TCanvas *c = new TCanvas("plotEvent", "plotEvent", 800, 400);
  if(show_geom==0){
    ccc->Divide(1, 0);
    ccc->GetPad(1)->Divide(2, 0);
    ccc->SetWindowSize(800,400);
    cout<<"single event view"<<endl;
  }
  else if (show_geom>0){
    //vertical canvas
    //    ccc->Divide(1, 2);
    //    ccc->GetPad(1)->SetPad(.005, .6525, .995, .995);
    //ccc->GetPad(1)->Divide(2, 0);
    //ccc->GetPad(2)->SetPad(.005, .005, .995, .6475);
    //ccc->SetWindowSize(500, 700);
    //horizontal canvas
    ccc->Divide(2, 0);
    ccc->GetPad(1)->SetPad(.005, .005, .3475, .995);
    ccc->GetPad(1)->Divide(0, 2);
    ccc->GetPad(2)->SetPad(.3525, .005, .995, .995);
    ccc->SetWindowSize(1200, 700);
    
  }

  ccc->cd(1)->cd(1)->SetGrid();
  Float_t bandwidth = 1e9*sampleRate/20.;//bandwith in Hz/10
  Float_t kB = 1.831e-23;
  Float_t thermal_noise = sqrt(kB*300.*50.*bandwidth);//thermal noise (V)
  int nbins = eventHist[txindex][rxindex]->GetNbinsX();
  //  cout<<bandwidth<<thermal_noise<<endl<<setprecision(10)<<nbins;
  Float_t noise, r1, r2, val;
  int j=0;
  for(int i=0;i<nbins;i++){
    j=i-(sampleRate*20);
    ran->Rannor(r1, r2);
    if(noise_flag==1){
      noise = r1*thermal_noise;
      //ev->eventHist->Fill(i, noise);
      eventHist[txindex][rxindex]->AddBinContent(i, noise);
      //    val=ev->eventHist->GetBinContent(i);
      }
    if(noise_flag>1){
      noise=r1*(double)noise_flag;
      eventHist[txindex][rxindex]->AddBinContent(i, noise);
      //    val=ev->eventHist->GetBinContent(i);
      }
    if(SINE_SUBTRACT==1){
      sineSubtract(txindex, rxindex);
    }
  }
  eventHist[txindex][rxindex]->GetXaxis()->SetTitle("Time (ns)");
  eventHist[txindex][rxindex]->GetYaxis()->SetTitle("V");
  eventHist[txindex][rxindex]->Draw("histl");
  eventHist[txindex][rxindex]->SetStats(0);
  //timme->Draw();
  ccc->cd(1)->cd(2);
  // fft = dofft(eventHist);
  // fft->GetXaxis()->SetRangeUser(0, fft->GetXaxis()->GetXmax()/2);
  // fft->SetTitle("psd");
  // fft->GetYaxis()->SetTitleOffset(1.3);
  // fft->GetXaxis()->SetTitle("freq (GHz)");
  // fft->GetYaxis()->SetTitle("power (dBm/Hz)");
  // fft->Draw();
  // ccc->cd(3);
  //  g.plot(vals, "with lines");
  spectrogram(txindex, rxindex, bins, overlap);
  //  TImage *img =TImage::Create();
  //  img->FromPad(c);
  //img->WriteImage("54mhz_100TeV_proton.png");

  if(show_geom>0){
    rxhist->Reset();
    txhist->Reset();
    vertexhist->Reset();
    triggeredhist->Reset();
    pointingHist->Reset();
    ccc->cd(2);
    //TH3F *rxhist = new TH3F("rxhist", "rxhist", 101, 1, -1, 101, 1, -1, 101, 1, -1);
    
    for(int i=0;i<nrx;i++){
      if(i==rxindex)continue;
      rxhist->Fill(1.+rx[i].z()/1000., 1.+rx[i].x()/1000., 1.+rx[i].y()/1000., eventHist[txindex][rxindex]->GetMaximum());
    }
    rxhist->BufferEmpty();
    RXHIST_FILLED=1;
  

    txhist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    vertexhist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    triggeredhist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    pointingHist->SetBins(101, rxhist->GetXaxis()->GetXmin(), rxhist->GetXaxis()->GetXmax(), 101,rxhist->GetYaxis()->GetXmin(), rxhist->GetYaxis()->GetXmax(),101, rxhist->GetZaxis()->GetXmin(), rxhist->GetZaxis()->GetXmax());
    
    for(int i=0;i<ntx;i++){
      txhist->Fill(1.+tx[i].z()/1000., 1.+tx[i].x()/1000., 1.+tx[i].y()/1000., 1.);
      //      cout<<tx[i].z()<<" "<<tx[i].x()<<" "<<tx[i].y()<<endl;    
    }
    triggeredhist->Fill(1.+rx[rxindex].z()/1000., 1.+rx[rxindex].x()/1000., 1.+rx[rxindex].y()/1000., 1.);
    vertexhist->Fill(1.+position.z()/1000., 1.+position.x()/1000., 1.+position.y()/1000., 1.);

    //see if we can reconstruct the event
    if(POINTING_MAP_BUILT==0){
      buildMap();
    }
    HepLorentzVector vvv = pointUsingMap();

    //pointingHist->Fill(1.+vvv.z()/1000., 1.+vvv.x()/1000., rxhist->GetZaxis()->GetXmax()+(vvv.y()/1000.), 1.);
    //    HepLorentzVector vvv = pointingTest();
    
    //    pointingHist->Fill(1.+vvv.z(), 1.+vvv.x(), vvv.y(), 1.);
    pointingHist->Fill(1.+vvv.z()/1000., 1.+vvv.x()/1000., (vvv.y()/1000.), 1.);

    txhist->SetMarkerStyle(3);
    txhist->SetMarkerColor(kRed);
    rxhist->SetMarkerStyle(8);
    rxhist->SetMarkerColor(kBlack);
    rxhist->SetTitle("geometry");
    rxhist->GetXaxis()->SetTitle("x (meters)");
    rxhist->GetXaxis()->SetTitleOffset(1.5);
    rxhist->GetYaxis()->SetTitle("y (meters)");
    rxhist->GetYaxis()->SetTitleOffset(1.5);
    rxhist->GetZaxis()->SetTitle("z (meters)");
    rxhist->GetZaxis()->SetTitleOffset(1.5);
    vertexhist->SetMarkerStyle(34);
    vertexhist->SetMarkerColor(kViolet);
    pointingHist->SetMarkerStyle(21);
    pointingHist->SetMarkerColor(kBlue);
    pointingHist->SetMarkerSize(1.5);
    txhist->SetMarkerSize(3);
    rxhist->SetMarkerSize(2);
    vertexhist->SetMarkerSize(2);
    triggeredhist->SetMarkerColor(kGreen);
    triggeredhist->SetMarkerStyle(8);
    triggeredhist->SetMarkerSize(2.5);

    TPolyLine3D *shower_indicator_line = new TPolyLine3D(2);

    double scale = rxhist->GetXaxis()->GetXmax()/4.;
    //cout<<scale<<endl;
    shower_indicator_line->SetPoint(0,1.+position.z()/1000., 1.+position.x()/1000., 1.+position.y()/1000.);

    shower_indicator_line->SetPoint(1,(1.+position.z()/1000.)+((direction.z())*scale), (1.+position.x()/1000.)+((direction.x())*scale), (1.+position.y()/1000.)+((direction.y())*scale));
    //cout<<"nere"<<endl;
    shower_indicator_line->SetLineWidth(3);
    shower_indicator_line->SetLineColor(kViolet);

    //cout<<"nereee"<<endl;
    rxhist->Draw("p");
    txhist->Draw("p same");
    vertexhist->Draw("p same");
    triggeredhist->Draw("p same");
    shower_indicator_line->Draw("same");
    pointingHist->Draw("p same");
    TLegend *leg = new TLegend(.7,.7,.9,.9);

    leg->AddEntry(txhist, "transmitter", "p");
    leg->AddEntry(vertexhist, "vertex", "p");
    leg->AddEntry(shower_indicator_line, "shower", "l");
    leg->AddEntry(rxhist, "receivers", "p");
    leg->AddEntry(triggeredhist, "this receiver", "p");
    leg->AddEntry(pointingHist, "source approx", "p");
    leg->Draw();
    //    ccc->Update();
  }
  ccc->Draw();
}

//this is not the best. requires 5 antennas. but pointing is pretty OK...
//TODO: find another (better) method.

// HepLorentzVector RadioScatterEvent::findSource(int debug){
//   HepLorentzVector source;
//   source.setX(9999999);
//   source.setY(9999999);
//   source.setZ(9999999);
//   //if we don't have enough ants to point, return nonsense.
//   if(triggered(thermalNoiseRMS()*2., 5)==0){
    
//     return source;
//   }
//   else{ 
//     HepLorentzVector dr[nrx];
//     double aa[nrx], bb[nrx], cc[nrx], dd[nrx];
//     int num_ants=0;
//     vector<double>gsl_a_dat, gsl_b_dat;
//     double tmin=9999999.;
//     for(int i=0;i<nrx;i++){
//       getComplexEnvelope(0, i, 200);//puts results in ceHist
//       rx[i].setT(ceHist->GetXaxis()->GetBinCenter(ceHist->GetMaximumBin()));
//       tmin=rx[i].t()<tmin?rx[i].t():tmin;
//     }
//     //works well for the 15 antenna case
//     for(int i=1;i<nrx;i++){
//       if(!trigSingle(thermalNoiseRMS()*2., i))continue;
//       dr[i].setT((rx[i].t()-rx[0].t()));
//       dr[i].setX((rx[i].x())-(rx[0].x()));
//       dr[i].setY((rx[i].y())-(rx[0].y()));
//       dr[i].setZ((rx[i].z())-(rx[0].z()));
//       //      if(debug==1)cout<<dr[i].t()<<" "<<dr[i].x()<<endl;
//       if(i>1){
//     	aa[i]=(2.*dr[i].x()/dr[i].t())-(2.*dr[1].x()/dr[1].t());
//     	bb[i]=(2.*dr[i].y()/dr[i].t())-(2.*dr[1].y()/dr[1].t());
//     	cc[i]=(2.*dr[i].z()/dr[i].t())-(2.*dr[1].z()/dr[1].t());
//     	dd[i]=dr[i].t()-dr[1].t()-((pow(dr[i].x(), 2)+pow(dr[i].y(), 2)+pow(dr[i].z(), 2))/dr[i].t())+((pow(dr[1].x(), 2)+pow(dr[1].y(), 2)+pow(dr[1].z(), 2))/dr[1].t());
      
//     	gsl_a_dat.push_back(aa[i]);
//     	gsl_a_dat.push_back(bb[i]);
//     	gsl_a_dat.push_back(cc[i]);
      
//     	gsl_b_dat.push_back(-dd[i]);
//     	num_ants++;
//       }
//     }
//     if(num_ants<5)return source;
//     gsl_matrix_view amat = gsl_matrix_view_array(&gsl_a_dat[0], num_ants, 3);

//     gsl_vector_view bvec = gsl_vector_view_array(&gsl_b_dat[0], num_ants);

//     gsl_vector *tau = gsl_vector_alloc(3);

//     gsl_vector *xvec = gsl_vector_alloc(3);
//     gsl_vector *resid = gsl_vector_alloc(num_ants);


    

//     //gsl_linalg_QR_decomp(&amat.matrix, tau);
//     //gsl_linalg_QR_lssolve(&amat.matrix, tau, &bvec.vector, xvec, resid);

//     gsl_permutation *p = gsl_permutation_alloc(3);
//     int signum=1;
//     gsl_vector *norm = gsl_vector_alloc(3);
//     gsl_linalg_QRPT_decomp(&amat.matrix, tau, p, &signum, norm);
//     gsl_linalg_QRPT_lssolve(&amat.matrix, tau, p, &bvec.vector, xvec, resid);

    
//     source.setX(gsl_vector_get(xvec, 0));
//     source.setY(gsl_vector_get(xvec, 1));
//     source.setZ(gsl_vector_get(xvec, 2));
//     // gsl_matrix_free(&amat.matrix);
//     // gsl_vector_free(&bvec.vector);
//     // gsl_vector_free(xx);
//     // gsl_vector_free(tau);
//     // gsl_vector_free(resid);
//     if(debug==1){
//       //      gsl_vector_fprintf (stdout, xvec, "%g");
//       //cout<<endl<<aa[4]<<" "<<bb[4]<<" "<<cc[4]<<endl;
//       //cout<<"matrix: "<<endl;

//       //gsl_matrix_fprintf(stdout, &amat.matrix, "%g");
//       cout<<"pointed: "<<source.x()/1000.<<" "<<source.y()/1000.<<" "<<source.z()/1000.<<endl;
//       cout<<"true: "<<position.x()/1000.<<" "<<position.y()/1000.<<" "<<position.z()/1000.<<endl;
//     }
//     return source;
//   }
// }

int RadioScatterEvent::buildMap(){
  //decided on 40x40x40 points of resoultion. arbitrary.
  //dt is [64000][3] and source is  [64000] which is a linearized matrix.
  double ndiv=40;
  //this sets the volume where we'll calculate the dts
  //redundant but that's OK
  if(RXHIST_FILLED==0){
    for(int i=0;i<nrx;i++){
      rxhist->Fill(1.+rx[i].z()/1000., 1.+rx[i].x()/1000., 1.+rx[i].y()/1000., 1.);
    }
    rxhist->BufferEmpty();
    RXHIST_FILLED=1;
  }

  double xmin=rxhist->GetYaxis()->GetXmin()*m;
  double ymin=rxhist->GetZaxis()->GetXmin()*m;
  double zmin=rxhist->GetXaxis()->GetXmin()*m;

  double xmax=rxhist->GetYaxis()->GetXmax()*m;
  double ymax=rxhist->GetZaxis()->GetXmax()*m;
  double zmax=rxhist->GetXaxis()->GetXmax()*m;


  double xdiv = (xmax-xmin)/ndiv;
  double ydiv = (ymax-ymin)/ndiv;
  double zdiv = (zmax-zmin)/ndiv;

  //make a linearized 'matrix' of time delays
  //and another for the x,y,z coords for each calculated delay (store as hep3vector)
  for(int i=0;i<ndiv;i++){
    for(int j=0;j<ndiv;j++){
      for(int k=0;k<ndiv;k++){
	int index=(((i*ndiv)+j)*ndiv)+k;
	source[index].setX(xmin+(i*xdiv));
	source[index].setY(ymin+(j*ydiv));
	source[index].setZ(zmin+(k*zdiv));
	//antenna pairs hardcoded for now, but will eventually be a loop over all antennas (TODO)
	double dt0 = ((source[index].vect()-rx[1].vect()).mag()-(source[index].vect()-rx[0].vect()).mag())/c_light;
	double dt1 = ((source[index].vect()-rx[2].vect()).mag()-(source[index].vect()-rx[0].vect()).mag())/c_light;
	double dt2 = ((source[index].vect()-rx[4].vect()).mag()-(source[index].vect()-rx[0].vect()).mag())/c_light;
	dt[index][0]=dt0;
	dt[index][1]=dt1;
	dt[index][2]=dt2;

      }
    }
  }
  cout<<"pointing map built"<<endl;
  POINTING_MAP_BUILT=1;
  return 1;
}


HepLorentzVector RadioScatterEvent::pointUsingMap(){
  if(POINTING_MAP_BUILT==0){
    buildMap();
  }
  //find the measured dts.
  for(int i=0;i<nrx;i++){
    //find the peak of the event envelope (corresponds to shower max)
    getComplexEnvelope(0, i, 200);//puts results in ceHist
    rx[i].setT(ceHist->GetBinCenter(ceHist->GetMaximumBin()));
  }

  //calculate the dts for this event for the corresponding antenna pairs
  //hardcoded for now, will update later
  double dt0=rx[1].t()-rx[0].t();
  double dt1=rx[2].t()-rx[0].t();
  double dt2=rx[4].t()-rx[0].t();

  //make sure we have 4 antennas (requisite number to point 3 pairs)
  if(triggered(thermalNoiseRMS()*2, 4)==0){
    HepLorentzVector dummy(99999999,99999999,99999999,99999999);
    cout<<"can't point, not enough antennas"<<endl;
    return dummy;
  }
  
  int ind0, ind1, ind2, index;
  double min0=999999, min1=999999,min2=999999, min=9999;
//minimize time differences
  for(int j=0;j<64000;j++){
    //pair 0
    if(abs(dt[j][0]-dt0)<min0){
      ind0=j;
      min0=abs(dt[j][0]-dt0);
    }
    //pair 1
    if(abs(dt[j][1]-dt1)<min1){
      ind1=j;
      min1=abs(dt[j][1]-dt1);
    }
    //pair 2
    if(abs(dt[j][2]-dt2)<min2){
      ind2=j;
      min2=abs(dt[j][2]-dt2);
    }
    //minize the sum of all three, works OK.
    //don't need the others really, can do this with just 1 loop through
    //the array. can do interpolation on the values here to acheive
    //better pointing.
    if(abs(dt[j][0]-dt0)+abs(dt[j][1]-dt1)+abs(dt[j][2]-dt2)<min){
      index=j;
      min=abs(dt[j][0]-dt0)+abs(dt[j][1]-dt1)+abs(dt[j][2]-dt2);
    }
  }
  //for now, return the vector at that index. in future, can calculate a better
  //source vector using interpolation
  return source[index];
    
}
