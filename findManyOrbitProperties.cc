/*******************************************************************************
*                                                                              *
*  findManyOrbitProperties.cc                                                  *
*                                                                              *
*  C++ code written by Paul McMillan, 2007-                                    *
*  Lund Observatory, Lund University.                                          *
*  address: Box 43, SE-221 00 Lund, Sweden                                     *
*  e-mail:  paul@astro.lu.se                                                   *
*                                                                              *
*******************************************************************************/


#include "GalPot.h"
#include "OrbitIntegrator.h"
#include <ctime>
#include <omp.h>

int main(int argc,char *argv[])  
{

    double intTime=8000;
    int nData;
    
    
  ifstream file,data;
  ofstream output;
  string potfile = "pot/test.Tpot",line;
  GalaxyPotential *Phi;


 // Read potential from file 
  file.open(potfile.c_str());
  if(!file){
    cerr << "Input file does not exist. Filename: " << potfile << "\n";
    return 1;
  }
    
   Phi = new GalaxyPotential(file);

  //file.clear( );
  //file.seekg( 0, std::ios::beg );
    
    
  //GalaxyPotential Phi2(file);
 file.close();
  
    int nDisk=Phi->NumberofDisks() ;
    int nSphr=Phi->NumberofSpheroids();
    DiskPar * Dpar = new  DiskPar[nDisk];
    SphrPar * Spar = new SphrPar[nSphr];
    
    for(int i=0; i!=nDisk; i++)
        Dpar[i]=Phi->DiskParameter(i);
    for(int i=0; i!=nSphr; i++)
        Spar[i]=Phi->SpheroidParameter(i);
    
    GalaxyPotential Phi2(nDisk, Dpar, nSphr, Spar);

    
  if(argc<3) {
    cerr << "Input: input_file output_file\n";
    cerr << "input_file must contain columns: R z v_R v_z v_phi\n";
    cerr << "   (distance in kpc, velocity in km/s)\n";
    cerr << " output is in kpc (distances), km^2/s^2 (energy), kpc km/s (angular momentum)\n";
    return 1;
  }

  // Open file for input
  data.open(argv[1]);
  if(!data){
    cerr << "Input file does not exist. Filename: " << string(argv[1]) << "\n";
    return 0;
  }



  
  
  
  int nline;
    int count=0;
    std::cout<<"Counting line"<<std::endl << std::flush;
    while (getline(data,line)) {
        count+=1;
    }
    nData=count;
    data.clear( );
    data.seekg( 0, std::ios::beg );
    std::cout<<"Ndata"<<"  "<<nData<<std::endl<<std::flush;
    
    
    // Setup class
    Vector <double,6> XV=1.;
    Vector <double,6> XVSun=1.;
    Vector<double,6> *XVL = new Vector<double,6> [nData];
    double Rsun=8.2;
    
    XVSun[0] = Rsun*Units::kpc;
    XVSun[1] = 0;
    XVSun[2] = 0;
    XVSun[3] = 0;
    XVSun[4] = 0;
    XVSun[5] = sqrt(Phi2.vcsquare(Rsun));
    
    
  // read file
    std::cout<<"Reading file"<<std::endl << std::flush;
    for (int i=0; i<nData; i++) {
        getline(data,line);
    if(line[0] != '#') {
      std::stringstream ss(line);
    
        for(int j=0;j!=6;j++)
            ss >> XVL[i][j];
        
    double Vc=sqrt(Phi2.vcsquare(XVL[i][0]));
    // Convert input to code coordinates
      XVL[i][0] *= Units::kpc;
      XVL[i][1] *= Units::kpc;
      XVL[i][2] *= Units::degree;
      XVL[i][3] *= Units::kms;
      XVL[i][4] *= Units::kms;
      XVL[i][5] =  XVL[i][5]*Units::kms+Vc;
      // Convert input to code coordinates
      
        
        
    }
    }
    data.close();
    
    
    
    int nOut=20;
    int IntegrationFail;
    Vector<double,7> *OVL = new Vector<double,7> [nData];
    OrbitIntegratorWithStats* OIVEC= new  OrbitIntegratorWithStats[nData];
    
    
    
   
    // read file
    std::cout<<"Kicking stars...../"<<std::endl;
	OrbitIntegratorWithStats OISun(XVSun, Phi, intTime);
    double *tOutS = new double[nOut];
    Vector<double,6> *OrbOutS = new Vector<double,6>[nOut];
    OISun.runWithOutputIncludingTime(OrbOutS,tOutS,nOut);
	double Phi_Sun=OrbOutS[nOut-1][2];
    OrbitIntegratorWithStats OI;
    

    
    #pragma omp parallel
    {
        
    OrbitIntegratorWithStats OI;
    double *tOut = new double[nOut];
    Vector<double,6> *OrbOut = new Vector<double,6>[nOut];

    Potential *Phic = new GalaxyPotential(nDisk, Dpar, nSphr, Spar);

    
    #pragma omp for
    for (int i=0; i<nData; i++) {


        OI.setup(XVL[i], Phic, intTime);
        OI.runWithOutputIncludingTime(OrbOut,tOut,nOut);
 

        
            for(int j=0;j!=6;j++)
            {
				if(j==2)
				{
                OVL[i][j]=OrbOut[nOut-1][j]-Phi_Sun;	
				}
				else
				{
                OVL[i][j]=OrbOut[nOut-1][j];
				}
            }
        
            OVL[i][6]=tOut[nOut-1];

    

        
    }
}
    
    
    
    //Writing output
    std::cout<<"Writing results..../"<<std::endl;
    
    // Open for output
    output.open(argv[2]);
    
    // Write header
    output << "#R z Phi VR Vz Vphi T  Rini   zini Phiini VRini Vzini Vphiini \n" << std::flush;
    
    for (int i=0; i<nData; i++) {
        
        output <<  OVL[i][0] / Units::kpc  << ' '
        << OVL[i][1] / Units::kpc  << ' '
        << OVL[i][2] / Units::degree  << ' '
        << OVL[i][3] / Units::kms  << ' '
        << OVL[i][4] / Units::kms  << ' '
        << OVL[i][5] / Units::kms  << ' '
        <<  OVL[i][6] << ' '
        << XVL[i][0] / Units::kpc  << ' '
        << XVL[i][1] / Units::kpc  << ' '
        << XVL[i][2] / Units::degree  << ' '
        << XVL[i][3] / Units::kms  << ' '
        << XVL[i][4] / Units::kms  << ' '
        << XVL[i][5] / Units::kms << ' '
        << std::endl << std::flush;
        
    }
    
    output.close();
  

  


  return 0;

}
