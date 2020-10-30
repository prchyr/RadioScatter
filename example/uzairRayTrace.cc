#include <RadioScatter/IceRayTracing.hh>

/*To compile this code just do:*/
/* g++ -std=c++0x uzairRayTrace.cc -o uzairRayTrace `root-config --cflags --glibs --libs`  -lRadioScatter */

/* This is an example script to call and play with Uzair's Ray Tracer */

int main(int argc, char**argv){

  if(argc<7){
    cout<<"More Arguments needed!"<<endl;
    cout<<"Example run command: ./IceRayTracing 0 0 -100 0 100 -5"<<endl;
    cout<<"Here the first three arguments are Tx coordinates i.e. x_Tx=0 m, y_Tx=0 m, z_Tx=-100 m  and the second three arguments are Rx coordinates i.e. x_Rx=0 m, y_Rx=100 m, z_Rx=-5 m "<<endl;
    return 0;
  }
  if(argc==7){
    cout<<"Tx set at X="<<atof(argv[1])<<" m, Y="<<atof(argv[2])<<" m, Z="<<atof(argv[3])<<" m"<<endl;
    cout<<"Rx set at X="<<atof(argv[4])<<" m, Y="<<atof(argv[5])<<" m, Z="<<atof(argv[6])<<" m"<<endl;
  } 
  if(argc>7){
    cout<<"More Arguments than needed!"<<endl;
    cout<<"Example run command: ./IceRayTracing 0 0 -100 0 100 -5"<<endl;
    cout<<"Here the first three arguments are Tx coordinates i.e. x_Tx=0 m, y_Tx=0 m, z_Tx=-100 m  and the second three arguments are Rx coordinates i.e. x_Rx=0 m, y_Rx=100 m, z_Rx=-5 m "<<endl;
    return 0;
  }

  double TxCor[3]={atof(argv[1]),atof(argv[2]),atof(argv[3])};
  double RxCor[3]={atof(argv[4]),atof(argv[5]),atof(argv[6])};
  
  /* An example for Direct and Reflected Rays */
  // double TxCor[3]={0,0,-100};
  // double RxCor[3]={0,100,-5};

  /* An example for Direct and Refracted Rays */
  // double TxCor[3]={0,0,-67.5489};
  // double RxCor[3]={0,1338.3,-800};
  
  double x0=0;
  double z0=TxCor[2];
  double x1=sqrt(pow(TxCor[0]-RxCor[0],2)+pow(TxCor[1]-RxCor[1],2));
  double z1=RxCor[2];

  double * getresults=IceRayTracing::IceRayTracing(x0,z0,x1,z1);

  cout<<" "<<endl;
  cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
  if(getresults[6]!=0 && getresults[7]!=0 && getresults[8]!=0){
    cout<<"No Possible Ray Paths between Tx and Rx!!! :("<<endl;
  }

  double A0, frequency, Lvalue;
  if(getresults[6]!=0 && getresults[8]!=0){ 
    cout<<"Direct and Refracted: dt(D,R)="<<(getresults[5]-getresults[3])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[6]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[3]*pow(10,9)<<" ns"<<endl;
    A0=1;
    frequency=0.1;//in GHz
    Lvalue=getresults[12];
    cout<<"Total Attenuation for Direct is (at "<<frequency<<" GHz): "<<IceRayTracing::GetTotalAttenuationDirect (A0,frequency,z0,z1,Lvalue)<<endl;

    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;
    cout<<"zmax: "<<getresults[15]<<" m"<<endl;
    A0=1;
    frequency=0.1;//in GHz
    Lvalue=getresults[14];
    cout<<"Total Attenuation for Refracted is (at "<<frequency<<" GHz): "<<IceRayTracing::GetTotalAttenuationRefracted (A0,frequency,z0,z1,getresults[15],Lvalue)<<endl;
  }
      
  if(getresults[6]!=0 && getresults[7]!=0){
    cout<<"Direct and Reflected: dt(D,R)="<<(getresults[4]-getresults[3])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[6]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[3]*pow(10,9)<<" ns"<<endl;
    A0=1;
    frequency=0.1;//in GHz
    Lvalue=getresults[12];
    cout<<"Total Attenuation for Direct is (at "<<frequency<<" GHz): "<<IceRayTracing::GetTotalAttenuationDirect (A0,frequency,z0,z1,Lvalue)<<endl;
    
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[7]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;   
    cout<<"Incident Angle in Ice on the Surface: "<<getresults[11]<<" deg"<<endl;
    A0=1;
    frequency=0.1;//in GHz
    Lvalue=getresults[13];
    cout<<"Total Attenuation for Reflected is (at "<<frequency<<" GHz): "<<IceRayTracing::GetTotalAttenuationReflected (A0,frequency,z0,z1,Lvalue)<<endl;
  }
      
  if(getresults[8]!=0 && getresults[7]!=0 && getresults[6]==0){
    cout<<"Refracted and Reflected: dt(D,R)="<<(getresults[4]-getresults[5])*pow(10,9)<<" ns "<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[7]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    cout<<"Incident Angle in Ice on the Surface: "<<getresults[11]<<" deg"<<endl;
    A0=1;
    frequency=0.1;//in GHz
    Lvalue=getresults[13];
    
    cout<<"Total Attenuation for Reflected is (at "<<frequency<<" GHz): "<<IceRayTracing::GetTotalAttenuationReflected (A0,frequency,z0,z1,Lvalue)<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;
    cout<<"zmax: "<<getresults[15]<<" m"<<endl;
    A0=1;
    frequency=0.1;//in GHz
    Lvalue=getresults[14];
    cout<<"Total Attenuation for Refracted is (at "<<frequency<<" GHz): "<<IceRayTracing::GetTotalAttenuationRefracted (A0,frequency,z0,z1,getresults[15],Lvalue)<<endl;
  }
  
  delete []getresults;  

  return 0;
  
}

