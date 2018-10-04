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
#include <iomanip>
#include <sstream>

int main(int argc,char *argv[])  
{

    double intTime=8000;
    int nData;
    
  const double PI = 3.141592653589793;
  const double torad = 180/PI;
  ifstream file,data;
  ofstream output, outputor;
  string potfile = "pot/PJM17_best.Tpot",line;
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
    Vector<double,6> *XVL = new Vector<double,6> [nData];
    double Rsun=8.2;
    
    
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
      XVL[i][5] *= Units::kms;
      // Convert input to code coordinates
      
        
        
    }
    }
    data.close();
    
    std::istringstream ss(argv[4]);
    int nOut;
    ss>>nOut;
    

    int IntegrationFail;
    Vector<double,7> *OVL = new Vector<double,7> [nData];
    Vector<double,12> *EVL = new Vector<double,12> [nData];
    Vector<double, 8> *PSL = new Vector<double,8> [nData*nOut];
    OrbitIntegratorWithStats* OIVEC= new  OrbitIntegratorWithStats[nData];
    
    
    
   
    // read file
    std::cout<<"Norbits  "<<nOut <<std::endl;
    std::cout<<"Wait..... stars are orbiting/"<<std::endl;
    OrbitIntegratorWithStats OI;
    

    
    #pragma omp parallel
    {
        
    OrbitIntegratorWithStats OI;
    double *tOut = new double[nOut];
    Vector<double,6> *OrbOut = new Vector<double,6>[nOut];

    Potential *Phic = new GalaxyPotential(nDisk, Dpar, nSphr, Spar);
	double L,Lz,LR,LT,inc,incmax,incmin,incmean,Lzmin,Lzmax,Lmin,Lmax,thetamax;
	int iini, ifin;
		
	
    #pragma omp for
    for (int i=0; i<nData; i++) {
		

        OI.setup(XVL[i], Phic, intTime);
        OI.runWithOutputIncludingTime(OrbOut,tOut,nOut);
		
		iini=i*nOut;
		ifin=iini+(nOut-1);
		
		
        if(OI.run() == 0) {
			
			incmin=2*PI;
			incmax=0;
			incmean=0;
            for(int j=0;j!=nOut;j++)
            {	
				//MOM ang stuff
				Lz=OrbOut[j][0]*OrbOut[j][5]; //R*Vt
				LR=-OrbOut[j][1]*OrbOut[j][5]; //-z*Vt
				LT=OrbOut[j][1]*OrbOut[j][3]-OrbOut[j][0]*OrbOut[j][4]; //z*VR-R*Vz
				L=sqrt(LR*LR+Lz*Lz+LT*LT);
				inc=acos(Lz/L);
				incmean+=inc;
				if (inc<incmin){
					incmin=inc;
				}
				if (inc>incmax){
					incmax=inc;
				}	
				
				//Save orbits stuff
				PSL[iini+j][0]=i;
				PSL[iini+j][7]=tOut[j];
				for(int k=1; k!=7; k++){
					PSL[iini+j][k]=OrbOut[j][k-1];
				}
					
				
				
            }

			incmean=(incmean/nOut)*torad;
			incmin=incmin*torad;
			incmax=incmax*torad;
			thetamax=atan(OI.Maxz/OI.MaxR)*torad;
			
            OVL[i][6]=tOut[nOut-1];
		EVL[i][0]=OI.Minr/Units::kpc; //rperi
	        EVL[i][1]=OI.Maxr/Units::kpc; //rapo
	        EVL[i][2]=OI.MinR/Units::kpc; //MaxR
	        EVL[i][3]=OI.MaxR/Units::kpc; //MinR
		EVL[i][4]=OI.Maxz/Units::kpc; //zmax
	        EVL[i][5]=(OI.Maxr-OI.Minr)/(OI.Maxr+OI.Minr); //ecc
	        EVL[i][6]=OI.Energy/(Units::kms*Units::kms);  //E
	        EVL[i][7]=OI.Lz/(Units::kpc*Units::kms); //Lz
	        EVL[i][8]=incmin; //imin
	        EVL[i][9]=incmax; //imax
	        EVL[i][10]=incmean; //imean
	        EVL[i][11]=thetamax; //thetamax
		} else {
				
				for(int j=0;j!=nOut;j++){
				//Save orbits stuff
				PSL[iini+j][0]=i;
				for(int k=1; k!=7; k++){
					PSL[iini+j][k]=OrbOut[j][k-1];
					}
				}
			
		        EVL[i][0]=1.e10;
		        EVL[i][1]=1.e10;
		        EVL[i][2]=1.e10;
		        EVL[i][3]=1.e10;
		        EVL[i][4]=1.e10;
		        EVL[i][5]=1.e10;
		        EVL[i][6]=1.e10;
		        EVL[i][7]=1.e10;
		        EVL[i][8]=1.e10;
		        EVL[i][9]=1.e10;
		        EVL[i][10]=1.e10;
		        EVL[i][11]=1.e10;
			
		}
			

        
    }
}
    

    
    //Writing output
    std::cout<<"Writing results..../"<<std::endl;
    
    // Open for output
    output.open(argv[2]);

    
    // Write header
    //output << "#R z Phi VR Vz Vphi T  Rini   zini Phiini VRini Vzini Vphiini \n" << std::flush;
    output << "#0-id 1-rmin 2-rmax 3-Rmin 4-Rmax 5-zmax 6-ecc 7-E 8-Lz 9-incmin 10-incmax 11-incmean 12-thetamax \n" << std::flush;
    
	//Output properties
    for (int i=0; i<nData; i++) {
        
        output
        << i << ' '  //id
		<< EVL[i][0] << ' '  //rperi
        << EVL[i][1] << ' '  //rapo
        << EVL[i][2] << ' '  //MinR
        << EVL[i][3]  << ' ' //MaxR
        << EVL[i][4]  << ' '  //zmax
        << EVL[i][5] << ' '  //ecc
        << EVL[i][6]  << ' ' //E
		<< EVL[i][7]  << ' '  //Lz
		<< EVL[i][8]  << ' '  //incin
		<< EVL[i][9]  << ' '  //incmax
		<< EVL[i][10]  << ' '  //incmean
		<< EVL[i][11]  << ' '  //thetamax
        << std::endl << std::flush;
         
        //output<< EVL[i][2] << std::endl << std::flush;
        
    }
	
    output.close();
	
	
	outputor.open(argv[3]);
	outputor << "#0-id 1-R 2-z 3-Phi 4-VR 5-Vz 6-Vphi 7-t \n" << std::flush;
	
	int j;	
	for (int i=0; i<nData*nOut; i++) {
	    j=floor(i/nOut);
            outputor
        	<<PSL[i][0]  << ' ' //id
        	<< PSL[i][1] / Units::kpc  << ' '   //R
        	<< PSL[i][2] / Units::kpc  << ' '   //z
        	<< PSL[i][3] / Units::degree  << ' '   //Phi
        	<< PSL[i][4] / Units::kms  << ' '   //VR
        	<< PSL[i][5] / Units::kms << ' '   //Vz
		<< PSL[i][6] / Units::kms  << ' '   //VPhi
		<< PSL[i][7]  / Units::Gyr << ' '   //time
		<< EVL[j][0] << ' '  //rperi
        	<< EVL[j][1] << ' '  //rapo
        	<< EVL[j][2] << ' '  //MinR
        	<< EVL[j][3]  << ' ' //MaxR
        	<< EVL[j][4]  << ' '  //zmax
        	<< EVL[j][5] << ' '  //ecc
        	<< EVL[j][6]  << ' ' //E
		<< EVL[j][7]  << ' '  //Lz
		<< EVL[j][8]  << ' '  //incin
		<< EVL[j][9]  << ' '  //incmax
		<< EVL[j][10]  << ' '  //incmean
		<< EVL[j][11]  << ' '  //thetamax
        	<< std::endl << std::flush;
				
	}
	
    outputor.close();
  

  


  return 0;

}
