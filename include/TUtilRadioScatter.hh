/*
  copyright s. prohira 2018
  released under the GNU General Public License version 3
*/


/*
The TUtilRadioScatter namespace provides numerous useful functions which act primarily on TGraphs. 

The FFT namespace has FFT functions built on ROOT's implementation of FFTW3. 

The SVD namespace has SVD functions built on the ROOT implementation of the GSL linear algebra tools.

The global system of units is: 

time: ns
length: mm
frequency: GHz
power spectral density: dBm/Hz
amplitude: V
charge: nC

 */

#ifndef TUTILRS_BASE_H
#define TUTILRS_BASE_H


#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TObject.h"
#include "TNamed.h"
#include "TLine.h"
#include "TFile.h"
#include "TSystemDirectory.h"
#include "TList.h"

//#include "cnpy.h"
#include "TSystemFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "THStack.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "TVectorD.h"
#include "TMatrixTBase.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TDecompSVD.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <complex>

//#include "TUtilGraph.hh"



using namespace std;







//class TUtilGraph;

namespace TUtilRadioScatter{

  /*****************************************************
  
  some useful units. inspired by the CLHEP global system of units. 

    use these to keep your numbers in the global system of units defined above:
    ns, GHz, mm, nC

    if you want to write something in terms of MHz, for example, just do 
    
    freq = 1200*MHz

    and then freq will have units of GHz, as it should.
  */
  /*lengths
 
    for example, if you wanted to calculate the time it took for a signal 
    to propagate 75 feet, you'd do:
    
    75*ft/c_light
    
    and it would return the correct time in nanoseconds. 
  */
static constexpr double mm = 1;
  static constexpr double m = 1000.*mm;
  static constexpr double ft = .3047*m;
  static constexpr double cm = .01*m;
//static constexpr double mm = .001*m;


  //time
  static constexpr double ns = 1.;
  static constexpr double us = ns*1e3;
  static constexpr double ms = ns*1e6;
  static constexpr double s = ns*1e9;

  //energy

  static constexpr double MeV = 1.;
  static constexpr double GeV = 1000.*MeV;
  static constexpr double KeV = .001*MeV;
  static constexpr double eV = 1e-6*MeV;
  
  //frequency
  static constexpr double GHz = 1.;
  static constexpr double MHz = .001*GHz;
  static constexpr double kHz = 1e-6*GHz;
  static constexpr double Hz = 1e-9*GHz;
  //constants
  static complex<double> I=sqrt(complex<double>(-1));//imaginary unit, used for complex stuff
  static constexpr double c_light = 2.9979246e8*m/s; //  mm/ns
  static constexpr double pi = 3.14159265358979323846; //radians
  static constexpr double twoPi = 2.*pi; //radians
  static constexpr double z_0=50; //ohms
  static constexpr double deg=pi/180.; //radians
  static constexpr double degree=pi/180.; //radians
  static constexpr double kB=8.617343e-11;//MeV/kelvin
  static constexpr double m_e=0.510998;//MeV/c^2
  static constexpr double kelvin=1;//kelvin
  static constexpr double classic_electr_radius = 2.8179403e-12;
  static constexpr double kBJoulesKelvin=1.38e-23;//J/kelvin
  static constexpr double rho=1.168e-3;//sea level density
  static constexpr double x_0=36.7;//radiation length in air
  static constexpr double e_0=.078;//ionization energy 


  
/*******************************************************
EXPERIMENTAL::NOT WORKING YET
  
  useful generic vector class (T can be anything, double, TGraph, TH1D, TH2D, etc) with one addition:

  ****can be simply written to a root file (inherit from TObject).*** 

can be extended to use other 'vector' features in the future. there is certainly a smarter way to do this, but i'm not a programmer.

  */

  
  //1D version
  template<class T>
  class TVec1D: public TNamed{
  private:
    std::vector<T> vec;
  public:
    TVec1D(){
      vec=std::vector<T>();
    }
    TVec1D(int N){
      vec=std::vector<T>(N);
    }
    TVec1D(int N, T const &  init){
      vec=std::vector<T>(N, init);
    }
    void push_back(T const & elem);
    T size();
    T  & operator[]( int index);
    void clear();
    
    ClassDefNV(TVec1D, 1);
    
  };
  template <class T> void TVec1D<T>::push_back(T const & elem) {
    vec.push_back(elem);
  }
  template <class T> T & TVec1D<T>::operator[]( int index){
    return vec[index];
  }
  template <class T> T TVec1D<T>::size( ){
    return vec.size();
  }
  template <class T> void TVec1D<T>::clear(){
    vec.clear();
  }

  //2D version

  template<class T>
  class TVec2D: public TNamed{
  private:
    std::vector<TUtilRadioScatter::TVec1D<T>> vec;
  public:
    TVec2D(){
      vec=std::vector<TUtilRadioScatter::TVec1D<T>>();
    }
     TVec2D(int N){
       vec=std::vector<TUtilRadioScatter::TVec1D<T>>(N);
    }
     TVec2D(int N, int M){
       vec=std::vector<TUtilRadioScatter::TVec1D<T>>(N, TUtilRadioScatter::TVec1D<T>(M));
     }
     TVec2D(int N, int M, T const & init){
       vec=std::vector<TUtilRadioScatter::TVec1D<T>>(N, TUtilRadioScatter::TVec1D<T>(M, init));
    }
     
    TUtilRadioScatter::TVec1D<T> & operator[](const int index);
    //    T & operator[][](int ind1, int ind2);
    void clear();

    ClassDef(TVec2D, 1);
  };

  template <class T> TUtilRadioScatter::TVec1D<T> & TVec2D<T>::operator[](const int index){
    return vec[index];
  }
  //template <class T> T & TVec1D<T>::operator[][]( int ind1, int ind2){
  // return vec[ind1][ind2];
  // }
  template <class T> void TVec2D<T>::clear(){
    for(int i=0;i<vec.size();i++){
      vec[i].clear();
    }
  }

  /**********************************************

utilities.

  ***********************************************/
  
  //volts to dbm/hz
  double vToDbmHz(double bandwidthGSs, double re, double im=0);
  //volts to dbm/Ghz
  double vToDbmGHz(double bandwidthGSs, double re, double im=0);
  //make an axis with linearly increasing values.
  double * makeIndices(int n, double step, double offset=0);
  //assign an array of offsets to the x axis of a graph. both must be the same length
  int assignXOffset(TGraph *inGr, double * offsets, double constant=1.);
  //the normalized sinc function: sin(pi x)/(pi x)
  double sinc(double x);
  //return a TGraph interpoltaed using simple sinc interpolation.(slow)
  TGraph * sincInterpolateGraph(TGraph *inGr, double interpGSs);
  //experimental
  TGraph * sincInterpolateGraphDev(TGraph *inGr, double interpGSs);
  //approximated (fast) sinc interpolation.(will not be very good for low N.) sacrifice accuracy for speed.
  TGraph * sincInterpolateGraphFast(TGraph *inGr, double interpGSs, int N=10);
  /*interpolate a tgraph.
    
    type 0 is a ROOT akima interpolation.(default)
    type 1 is a sinc interpolation
    type 2 is a fast sinc interpolation. not as accurate as the sinc, but faster. the parameter N is only used for this type.

    the sinc interpolation is best/most accurate for truly band-limited signals sampled near the nyquist frequency. the fast method works very well, and is usually recommended over the full sinc for most applications. the akima is better for oversampled signals where the frequency of interest is far below the nyquist frequency.
  */
  int getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph, double interpGSs, int type=2, int N=10);
  //return the interpolated graph (memory use) using ROOT's akima method
  TGraph * interpolateGraph(TGraph * inGraph, double interpGSs);
  //return the interpolated graph, evaluated at times. the input graph can be unevenly sampled.
  TGraph * interpolateGraph(TGraph * inGraph, vector<double> times);
  //same as above but different input types
  TGraph * interpolateGraph(TGraph * inGraph, int N, double* times);
  //normalize a graph
  TGraph * normalize(TGraph *inGr);
  //normalize to peak
  TGraph * normToPeak(TGraph *inGr);
  //normalize a 2d graph
  TGraph2D * normalize(TGraph2D * inGr);
  //normalize a 2d histo, happens in place! returns the norm.
  double normalize(TH2D * inGr, double ymin, double ymax);
  //return the norm of a th2d
  double norm(TH2D * inGr, double ymin, double ymax);
  //return a chunk of a graph, specified by x-axis values.
  //if shift_to_zero==1, the time axis is shifted such that it starts at 0.
  //  return the cumulative distribution function of a graph. if normed is 1, the graph is scaled such that x and y axes spread from 0 to 1.
  TGraph *CDF(TGraph *inGr, int normed=0);
  //get a chunk of a graph from start to end. delay to zero shifts the time so that the returned graph starts at t=0
  //uses TGraph->Eval() to get a very specific point in the graph. slow.
  TGraph * getChunkOfGraph(TGraph *ingr, double start, double end, int delay_to_zero=0);
//get a chunk of a graph from start to end. delay to zero shifts the time so that the returned graph starts at t=0
  TGraph * getChunkOfGraphFast(TGraph *ingr, double start, double end, int delay_to_zero=0);
  TGraph * getNSamplesFrom(TGraph *ingr, double start, int nSamples, int delay_to_zero);
  //get a specific range of samples
  TGraph * getTheseSamples(TGraph *ingr, int sampStart, int sampEnd, int delay_to_zero);
  //cross correlation of two graphs. returns the cross-correlation graph
  //maxDelay is the maximum starting offset between gr1 and gr2. defaults
  //to the full length of the graphs.
  //t_low(high) are the lowest(highest) times over which to
  //calculate the correlation. defaults to the full length of the graphs.
  TGraph * crossCorrelate(TGraph * gr1, TGraph * gr2, double max_delay=999999., double t_low=0., double t_high=999999.);
  //same as above but allows you to supply a window function, to only
  //correlate parts of the graph that you want.
  TGraph * crossCorrelateWindowed(TGraph * gr1, TGraph * gr2, TGraph *wingraph, double max_delay=999999., double t_low=0., double t_high=999999.);
  //the same as the crossCorrelate() function, but returns gr2 shifted in time
  //to the point of peak cross correlation with gr1.
  TGraph * align(TGraph * gr1, TGraph * gr2, double max_delay=999999., double t_low=0., double t_high=999999.);
  /*align gr2 to gr1, but returning othGr, delayed by the amount needed to
  align gr2 to gr1 at the point of peak correlation. 
  for example, gr1 and gr2 are events which are causal, and othGr
  is some other graph which you'd like to align with these, but can't for
  whatever reason (contaminated with CW, etc.). 
  */
  TGraph * alignToOther(TGraph * gr1, TGraph * gr2, TGraph* othGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //align a large number of graphs to the first graph in the set.
  vector<TGraph*> alignMultiple(vector<TGraph*> inGr, double max_delay=999999., double t_low=0., double t_high=999999.);

//return the average of a bunch of aligned graphs
  TGraph* alignMultipleAndAverage(vector<TGraph*> inGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //align a large number of graphs to the first, then truncate them all to t_min and t_max
  vector<TGraph*> alignMultipleAndTruncate(vector<TGraph*> inGr, double max_delay, double t_min, double t_max, double t_low=0., double t_high=999999.);
  //align a large set to a reference set.
  vector<TGraph*> alignMultipleToOther(vector<TGraph*> inGr, vector<TGraph*> othGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //delay a graph. 
  TGraph *delayGraph(TGraph *ingr, double delay);
  //same but with no mem usage.
  int delayGraph(TGraph * ingr, TGraph *outgr, double delay);
  //plot \Delta(gr1[i], gr2[i]) for each graph point i
  TH1F * plotResiduals(TGraph *gr1, TGraph *gr2, int nbins=40, double min=1, double max=-1);
  //average some graphs
  TGraph * avgGraph(vector<TGraph*> inGr);
  //get the absolute value of a graph
  TGraph * absGraph(TGraph *inGr);
  int absHist(TH2D *inGr);
  //add 2 TGraphs. if constant is -1, they are subtracted.
  TGraph * add(TGraph * g1, TGraph * g2, double constant=1.);
  TGraph2D * add(TGraph2D * g1, TGraph2D * g2, double constant=1.);
  //dot product of 2 graphs
  double dot(TGraph *g1, TGraph *g2);
  double dot(TGraph *g1, TGraph *g2, double tLow, double tHigh);
  //multiply two graphs: out(t)=g1(t)*consant*g2(t)
  TGraph * mult(TGraph *g1, TGraph *g2, double constant=1.);
  //divide 2 graphs: out(t)=g1(t)/constant*g2(t). if there is a divide
  //by zero, that entry will just be 0.
  TGraph * divide(TGraph *g1, TGraph *g2, double constant);
  //shift a graph along the y axis by the factor
  TGraph * shiftY(TGraph *g1, double factor);
  //shift a tprofile along the y axis by the factor
  TGraphErrors * shiftY(TProfile *g1, double factor);
  //shift a graph along the x axis by the factor
  TGraph * shiftX(TGraph *g1, double factor);
  //scale a TGraph by a constant factor
  TGraph * scale(TGraph *g1, double factor);
  //stretch a TGraph in time by a factor
  TGraph *stretch(TGraph *g1, double factor);
  //find the mean of a TGraph. range is optional
  double mean(TGraph *gr, double t_low=0., double t_high=999999.);
  //get the 'peakiness' of a graph
  double peakiness(TGraph *inGr);
  //get the max value (wrapper of TMath::max)
  double max(TGraph *gr);
  double maxInRange(TGraph *gr, double t_low=0., double t_high=999999.);
  //get the x-axis location of t max value (wraper of TMath::LocMax)
  double locMax(TGraph *gr);
  double locMaxInRange(TGraph *gr, double t_low=0., double t_high=999999.);
    double min(TGraph *gr);
  double minInRange(TGraph *gr, double t_low=0., double t_high=999999.);
  //get the x-axis location of t min value (wraper of TMath::LocMin)
  double locMin(TGraph *gr);
  double locMinInRange(TGraph *gr, double t_low=0., double t_high=999999.);
  //this snr is a ratio of the rms of the peak (e.g. the peak value/sqrt(2)) over the rms of the waveform up to the peak
  double snr(TGraph *gr, double noiseStart=0., double noiseEnd=0.);
  //sort a vector. returns a vector of pairs with (index, value) sorted by value in order where 0 is ascending, 1 is decending 
  vector<pair<int, double>> sort(vector<double> inVec, int order=0);
  //return the phase at a single point, degrees=0, rad=1
  /*DON"T USE*/
  double getInstPhase(TGraph *inGr, double t, int deg0rad1=1);
  //get the initial phase and amplitude of a graph
  int getInitPhaseAmp(TGraph *event, double freq, double * amp, double *phase);
  //get the index of a certain value. will always under-estimate
  int getIndex(TGraph *gr, double t);
  int getIndex(TGraph2D * gr, double t);
  vector<double> toVector(int *data, int n);
  vector<double> toVector(double *data, int n);
  TVec1D<double> toTVec1D(double *data, int n);

  //get the power (square all values of a graph and divide by 50 ohms)
  TGraph * power(TGraph *gr);
  //square all values of a graph
  TGraph * squared(TGraph *gr);
  int squareHist(TH2D *inGr);
  //remove the mean of a TGraph. range to compute the mean over is optional.
  //the mean computed within a sub-range will be removed from the full graph.
  TGraph * removeMean(TGraph *gr, double t_low=0., double t_high=999999.);
  //flip an array ([0]=[n-1]) (NOT WORKING)
  double * flip(int n, double * in);
  //flip a graph
  TGraph * flip(TGraph *inGr);
  //swap the x and y arrays of a TGraph
  TGraph * swap(TGraph * inGr);
  //same as removeMean but the original graph is changed
  int removeMeanInPlace(TGraph *gr, double t_low=0., double t_high=999999.);
  //make CW with given parameters.
  TGraph * makeCW(double freq,  double amp, double t_min=0., double t_max=1000., double GSs=20., double phase=0.);
  //sample CW at the given times
  TGraph * sampledCW(double freq,  double amp, int N, double * times, double phase);
  TGraph * sampledCW(double freq,  double amp, vector<double> times, double phase);
  //integrate a TGraph. lower and upper bounds are optional.
  double integrate(TGraph * gr, double t_low=0, double t_high=999999.); //get the sum
  double sumGraph(TGraph * gr, double t_low=0, double t_high=999999.);
  //get the RMS
  double rms(TGraph * gr, double t_low, double t_high);
  //get the power
  double avgPower(TGraph *gr, double t_low, double t_high);
  //get the amplitude
  double amplitude(TGraph * gr, double t_low, double t_high);
//take the derivative. if direction=-1, takes derivative along other direction of axis.
  TGraph * derivative(TGraph *gr, int direction=1);
  //get the observer graph from a retarded graph. tUnits is the multiplier for ns (eg for ms, tUinits=1e6);
  TGraph * gObs(TGraph *inGr, double thetaDeg, double tUnits=1.);
  //integrate a TGraph but square it first to get units of power.
  double integratePower(TGraph * gr, double t_low=0, double t_high=999999.);

  //integrate a histogram
  double integrate(TH2D *h, double xmin, double xmax, double ymin, double ymax);
  double integrateWithError(TH2D *h, double xmin, double xmax, double ymin, double ymax, double & err);
  //integrate a graph bin-by-bin, putting the result in a tgraph
  //binNS is the desired bin width in nanoseconds
  TGraph * integrateByBin(TGraph *gr, double binNS);
  //simple 2 pole lowpass filter
  TGraph * lowpassFilter(TGraph *ingr, double cutoff, int order=2);
  //highpass filter
  TGraph * highpassFilter(TGraph *ingr, double cutoff, int order=1);
  //a bandpass filter made of one of each of the above filters.
  TGraph * bandpassFilter(TGraph *ingr, double low, double high);
  //a brick wall frequency domain filter
  TGraph * brickWallFilter(TGraph *ingr, double low, double high);
    //will take the first chunk of the signal graph (equal to to t_high-t_low)
  vector<double> linspace(double start, double stop, int N);
  // return the value of a window over sample numbers n. types are:
  /*
    0=bartlett (triangle) (default)
    1=welch (parabolic)
    2=hann (gaussian ish)
    3=blackman-nuttall (gaussian ish);
   */
  double window(int i, int n, int type=0);
  //return the value of a bartlett window over sample numbers n
  double bartlettWindow(int i, int n);
  //return the value of a welch window over sample numbers n
  double welchWindow(int i, int n);
  //return the value of a hann window over sample numbers n
  double hannWindow(int i, int n);
  //return the value of a blackman-nuttall window
  double blackmanNuttallWindow(int i, int n);
  //apply a window of the selected type to the graph inGr in time window.
  //and add it to the indicated region of the background graph.
  TGraph * applyWindow(TGraph* inGr, double startt, double endt, int type=0);
  //plot a window with the given parameters.
  TGraph * plotWindow(double peakAmplitude, double len, double GSs, double startt, double endt, int type);
  TGraph * makeNullData(TGraph *sig, TGraph * back, double t_min, double t_max, double scale=1.);
  TGraph * makeNullDataFixedLength(TGraph *sig, TGraph *back, double t_min, int nSamps);
  double sidebandSubtraction2DWithErrors(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, double & err, int draw=0, Color_t color=kRed, double alpha=1.);
  double sidebandSubtractionXAxisWithErrors(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, double &err, int draw=0, Color_t color=kRed);
  double sidebandSubtractionYAxisWithErrors(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, double &err, int draw=0, Color_t color=kRed);
  double sidebandSubtraction2DDev(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, double & err, int draw=0, Color_t color=kRed);
  double sidebandSubtraction2D(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, int draw=0, Color_t color=kRed);

  //degrees to radians
  double deg2Rad(double deg);
  //radians to degrees
  double rad2Deg(double rad);
  //add some noise to a graph
  TGraph * addNoise(TGraph * inGr, double level);
  //find the x values of zero crossings. can also plot the relative time between subsequent zero crossings if relative = 1
  TGraph * getZeroCrossGraph(TGraph * inGr, int relative=0);
  //find the x values of zero crossings and put them in a histogram. returns the number of zero crossings in that graph.
  int fillZeroCrossHist(TGraph * inGr, TH1D* hist, double weight=1., double threshold=0.);
  //fill a histogram with an array
  TH1F * histogram(TGraph *gr, int nbins=20, TString title="hist", TString name="hist");

  //get the time of the first threshold crossing after after
  double getFirstThresholdCrossing(TGraph *inGr, double thresh, double after=0.);
  double getLastThresholdCrossing(TGraph *inGr, double thresh, double after=0.);


  //drawing things
  void style(TGraph *inGr, Color_t color, double lineWidth=1., int lineStyle=1);
  //utilities for common stuff such as titles and axes
  void titles(TGraph *inGr, TString title, TString xtitle, TString ytitle);
  void ranges(TGraph *inGr, double x1, double x2, double y1, double y2);
  void xrange(TGraph *inGr, double x1, double x2);
  void yrange(TGraph *inGr, double y1, double y2);

  void titles(TGraphErrors *inGr, TString title, TString xtitle, TString ytitle);
  void ranges(TGraphErrors *inGr, double x1, double x2, double y1, double y2);
  void xrange(TGraphErrors *inGr, double x1, double x2);
  void yrange(TGraphErrors *inGr, double y1, double y2);
  
  void titles(TProfile *inGr, TString title, TString xtitle, TString ytitle);
  void ranges(TProfile *inGr, double x1, double x2, double y1, double y2);
  void xrange(TProfile *inGr, double x1, double x2);
  void yrange(TProfile *inGr, double y1, double y2);
  void titles(TH1F *inGr, TString title, TString xtitle, TString ytitle);
  void ranges(TH1F *inGr, double x1, double x2, double y1, double y2);
  void xrange(TH1F *inGr, double x1, double x2);
  void yrange(TH1F *inGr, double y1, double y2);

  void titles(TH2D *inGr, TString title, TString xtitle, TString ytitle, TString ztitle);
  void ranges(TH2D *inGr, double x1, double x2, double y1, double y2, double z1=0, double z2=0);
  void xrange(TH2D *inGr, double x1, double x2);
  void yrange(TH2D *inGr, double y1, double y2);
  void zrange(TH2D *inGr, double z1, double z2);

  void titles(TH2F *inGr, TString title, TString xtitle, TString ytitle, TString ztitle);
  void ranges(TH2F *inGr, double x1, double x2, double y1, double y2);

  void xrange(TH2F *inGr, double x1, double x2);
  void yrange(TH2F *inGr, double y1, double y2);
  //draw a bunch of graphs
  void draw(vector<TGraph*> inGr, TString option="");
    //draw a bunch of graphs
  void draw(int nGraphs, TGraph ** inGr, TString option="");
  //draw a bunch of hists
  void draw(vector<TH1D*> inGr, TString drawOption1);
  //cast a th1d to a tgraph
  TGraph * toGraph(TH1D * hist);
  TGraph * toGraph(TH1F * hist);
  //draw cursors centered on the peak of a spectrogram.
  //before use, must draw the canvas or update with can->Update()
  //so that the gPad calls work correctly.
  vector<TLine*> drawPeakCursorXY(TH2D* inHist, Color_t color);
  //get a slice of a spectrogram from ylow to yhigh
  TGraph * getSliceY(TH2D* hist, double ylow, double yhigh);
  //get many slices
  vector<TGraph*> getSlicesY(TH2D * hist, int nSlices, double ylow, double yhigh, TString name="name", TString title="title");
  //a pretty warm palette
  void setWarmPalette();
  //a pretty cool palette
  void setCoolPalette();
  //cold palette
  void setColdPalette();
  void setHotPalette();
  void set2DPalette();
  TCanvas * canvas(TString title="can", TString name="can", int xdim=700, int ydim=600);
  //miscellaneous TGraph things


  //return a TGraph, evenly sampled along dt
  TGraph * evenSample(TGraph *inGr, double dt);
  //zero pad a tgraph. requires an evenly sampled tgraph.
  TGraph * zeroPad(TGraph *inGr, int num, int whichEnd=1);
  //get the distance between tvectors
  double distance3(TVector3 one, TVector3 two);
  //get the distance from one vector to several others
  double * distance3(int N, TVector3 one, TVector3 * two);
  //time of flight between two tvectors
  double timeOfFlight(TVector3 one, TVector3 two, double n=1.);
  //time of flight from one vector to several others
  double * timeOfFlight(int N, TVector3 one, TVector3 * two, double n=1.);
  //delta t between 2 antennas for a single source. it is the time a signal hits number 1 minus the time it hits number 2;
  double  dTimeOfFlight(TVector3 source,TVector3 one, TVector3 two, double n=1.);
  //delta t's between different antennas for a single source
  double ** dTimeOfFlight(int N, TVector3 one, TVector3 * two, double n=1.);
  //offset a bunch of antennas for the correct delta t's for single source
  //  TGraph ** delayGraphs

  //return the neutrino interaction length for a neutrino with log10(E) in GeV logE=6 is 1PeV, for example
  double lInt(double logE);

  /*******************************************************

The FFT namespace, for everything to do with FFTs.

  *******************************************************/

  
  namespace FFT{

    //returns a tgraph2d, x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.
    TGraph2D * fft(TGraph *inGr);
    //this version returns a complex vector from real input
    int fft(int n, double * in, complex<double> * out);
    //returns the inverse fft of the data. must be structured as above (x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.)
    TGraph * ifft(TGraph2D *inGr);
    //complex array to double
    int ifft(int n, complex<double> * in, double * out);
    //complex valued forward FFT.
    int cfft(int n, complex<double> * inVec, complex<double> * outVec);
    //complex valued backward FFT.
    int cifft(int n, complex<double> * inVec, complex<double> * outVec);
    //sine transform
    TGraph * sineTransform(TGraph * inGr);
    double * sineTransform(int n, double * in);

    void fftshift(int N, complex<double>* in);
    //return the Hilbert transform
    TGraph * hilbertTransform(TGraph *inGr);
    //plot the phase of the full graph 
    TGraph * plotPhase(TGraph *inGr);
    //get the phase at one frequency
    double getPhaseAt(TGraph *inGr, double freq);
    //zero the phase at the given frequency
    TGraph * zeroPhaseAt(TGraph *inGr, double freq, int debug=0);
    //set the phase angle at the given frequency
    TGraph * setPhaseAt(TGraph *inGr, double freq, double phaseAng, int debug=0);
    //get the power magnitude at some frequency
    double getMagAt(TGraph *inGr, double freq);
    //return the Hilbert envelope
    TGraph * hilbertEnvelope(TGraph *inGr);
    //return the power spectral density in dBm/Hz, rBW is the resolution bandwith of the system, used to calculate the density. defaults to Nyquist.

    TGraph * psd(TGraph *inGr, int dbFlag=1, double rBW=0.);


/*
      return the spectrogram, with various options:

binsize is what python calls nfft. it's the number of samples in the chunk over which the FFT is calculated. this is one spectrogram 'bin'

overlap is how many samples each bin overlaps with the next. helps somewhat with smoothing.

zero pad length is the length to which the chunk is symmetrically zero-padded.

win_type is an enumeration of window types to be applied to each bin. this helps avoid discontinuities and noise in the spectrogram. see the window function for the window types.
     */
    TH2D* spectrogram(TGraph *gr, Int_t binsize = 128, Int_t overlap=32, Int_t zero_pad_length=128, int win_type=0, int dbFlag=1, double ymin=0, double ymax=3.);
    TGraph2D* spectrogramGraph(TGraph *gr, Int_t binsize=64 , Int_t overlap=32, Int_t zero_pad_length=128, int win_type=0, int dbFlag=1, double ymin=0, double ymax=3.);
    //averages a vector of spectrograms. must be the same size.
    TH2D* avgSpectrograms(vector<TH2D*>  inh);

//plot the peak frequency for each bin, with the options as above

TGraph* peakFreqGraph(TGraph *gr, Int_t binsize , Int_t overlap, Int_t zero_pad_length, int win_type, double thresh=0.);
  }






  /*****************************************************

the SVD namespace, which has useful utilities for SVD filtration methods

  ******************************************************/

  namespace SVD{

    //normmalize a tvector
    TVectorD normalize(TVectorD vec);
    //return the norm of a vector
    double norm(TVectorD vec);
    //build a matrix out of events, with each event a row in the matrix.
    //must all be the same length.
    TMatrixD eventMatrix(vector<TGraph*> vecs);
    //make a density matrix partitioned along D
    TMatrixD densityMatrix(TGraph *vec, int D=1, int xlim=0);
    //returns a matrix after removing all coffs higher than below.
    //if above is provided, will return the matrix as constructed from
    //the modes between above and below.
    TMatrixD truncateSVD(TDecompSVD svd, int below, int above=0);
    //returns a matrix constructed from a single singular value
    TMatrixD reconstructSingle(TDecompSVD svd, int val);
    //cast a tgraph to a tvectord
    TVectorD toVector(TGraph *vec);
    //makes a row average of matrix m
    TVectorD avgVector(TMatrixD m);
    //builds a basis of patterns from the matrix m, up to the number of
    //patterns num
    TMatrixD buildBasis(TMatrixD m, int num=10);
    //returns the expansion coeffs for a vector in a basis B
    TVectorD getCoefficients(TVectorD V, TMatrixD B);
    //expands a vector in a basis B
    TVectorD expandInBasis(TVectorD V, TMatrixD B, int num=10);
    TGraph * expandInBasis(TGraph * G, TMatrixD B, int num=10);

    //filter a vector using a basis. the expansion of the vector in the basis B (to
    //order num) will be removed from the vector.
    TVectorD filter(TVectorD V, TMatrixD B, int num);
    TGraph * filter(TGraph *G, TMatrixD B, int num);
    //will align a signal matrix and a background matrix to their respective
    //reference matrices, build a basis out of backs ,
    //and filter sigs using this basis to order num
    TGraph * alignAndFilter(vector<TGraph*> sigs, vector<TGraph*> sigsref, vector<TGraph*> backs, vector<TGraph*> backsref, int num);

    //must be for square matrix, a Ralston-style filter matrix
    TMatrixD makeFilter(TDecompSVD svd, int below, int above=0);
    //flatten the indices of a matrix in row major.
    TVectorD flatten(TMatrixD m);
    //cast a vector to graph for plotting.
    TGraph * toGraph(TVectorD v, double samplerate=1., double delay=0., TString name="");
    //makes a 2d histo from a matrix, for visualization.
    TH2F * matrixMap(TMatrixD M, TString name="");


  }

  namespace SIM{
    //the shower age expression for the NKG approximation  from arXiv:1503.02808
    double ss(double x, double E, double x_0, double e_0);
    //number of leptons as a function of depth in radiation lengths in the NKG approximation from  arXiv:1503.02808
    double n(double x, double E, double x_0, double e_0);
  }


  namespace DFT{
    TGraph2D * udft(TGraph * inGr, double fSampMean=1.);
    TGraph * upsd(TGraph * inGr, double fSampMean=1., int log=0);
    TGraph * iudft(TGraph2D * inGr, double GSs);
  }

  namespace demod{
    void normalizedDeChirp(Int_t * one, Int_t * two, int offset, int insize, Float_t * out);
    
    TGraph * normalizedDeChirp(TGraph * one, int offset);
    TGraph * deChirp(TGraph * one, int offset);
  }
}


#endif
