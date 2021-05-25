#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

struct sphere
{
        double cx;
        double cy;
        double r;
        double charge;
};

void ComputePotential(sphere , sphere , double*, double*, double*, double*, double*, double , double , double , double , int , int , double , double  );

void ComputeElectricField(int , int ,double* ,double , double , double* , double* );

void spheres()
{
    sphere s1,s2;
    //setting values
    s1.cx=3;
    s1.cy=3;
    s1.r=1;
    s1.charge=60;
    s2.cx=7;
    s2.cy=7;
    s2.r=1;
    s2.charge=60;
    double Lx=10, Ly=10;
    int Nx=100, Ny=100;
    
    double thrs=0.1; //the value that needed to be reached in order to find the satisfying potential result. Clearly, the lower the threshold, the more precise will be the result (and the longer it will take). 
    double border=0.;
    //the followeing will help while computing potential 
    double* V0= new double[Nx*Ny];
    double* V1= new double[Nx*Ny];
    double* V2= new double[Nx*Ny];
    double* func=new double[Nx*Ny];//to charge and border distributions
    double*  mV= new double[Nx*Ny];
    
    double hx=Lx/(Nx-1.);
    double hy=Ly/(Ny-1.);
    std::ofstream fphi,fcharge,felectric;
    
    
    ComputePotential( s1,  s2, V0, V1, V2, func, mV,  hx,  hy,  Lx,  Ly,  Nx,  Ny,  thrs,  border);
    
  
    
    
    // Now, I need to compute Electric field as E=-grad(V):
     double* Ex= new double[Nx*Ny];
     double* Ey= new double[Ny*Nx];
    
      ComputeElectricField(Nx,Ny,V2,hx,hy, Ex, Ey);
  
    
    
    //scrive grafici
    fcharge.open("carica.txt");
    fphi.open("V.txt");
    felectric.open("electricField.txt");
    fcharge.precision(5);
    fphi.precision(5);
    felectric.precision(5);
    int iCount=0;
    for(int j=0;j<Ny;j++){
        double y;
        y=Ly*j/(Ny-1.);
        for(int i=0;i<Nx;i++){
            double x;
            x=Lx*i/(Nx-1.);
            
            fcharge << x << "  " << y << "  " << func[iCount]*(-2.*(1./pow(hx,2)+1./pow(hy,2))) << '\n';
            fphi << x << "  " << y << "  " << -V1[iCount] << '\n';
            felectric << x << "  " << y << "  " << Ex[iCount] << "  " <<Ey[iCount]<< '\n';
            iCount++;
        }
        fphi << '\n';
        fcharge << '\n';
    }
    
    fcharge.close();
    fphi.close();
    felectric.close();


    delete [] V0;
    delete [] V1;
    delete [] V2;
    delete [] func;
    delete [] mV;
    delete [] Ex;
    delete [] Ey;

    }

void ComputePotential(sphere s1, sphere s2, double*V0, double*V1, double*V2, double*func, double*mV, double hx, double hy, double Lx, double Ly, int Nx, int Ny, double thrs, double border )
{
	
	int iCount=0;
    for(int j=0;j<Ny;j++){
        double y;
        y=Ly*j/(Ny-1.);
        for(int i=0;i<Nx;i++){
            double x;
            x=Lx*i/(Nx-1.);
            //set charge distribution, according to position in grid we set charge 1, charge2 or 0
            if(sqrt(pow(x-s1.cx,2)+pow(y-s1.cy,2))<=s1.r) 
                     func[iCount]=s1.charge;
                    
            else if(sqrt(pow(x-s2.cx,2)+pow(y-s2.cy,2))<=s2.r)
                    func[iCount]=s2.charge;
                    
            else
                     func[iCount]=0.;


            //set border condition 
            if(i==0 || i==(Nx-1) || j==1 || j==(Ny-1))
                func[iCount]+=border;
              
            
            //now that we set charge and border, we can compute f inside border 
            if(i!=0 && i!=(Nx-1) && j!=0 && j !=(Ny-1)){
                func[iCount]/=(-2.*(1./pow(hx,2)+1./pow(hy,2)));
            }
            
            //update counter used to set func
            iCount++;
        }
    }
    
    
    //set border conditions
    for(int i=0;i<Nx*Ny;i++){
        V1[i]=func[i];
        V0[i]=func[i];
        V2[i]=0.;
      
    }
 
    double distance=1e10;//setting a very high value for distance so that we are sure that it will be for sure higher than threshold
    while( distance>thrs){
         for(int i=0;i<Nx*Ny;i++){
             mV[i]=0;
             
         }
        iCount=0;
        //we now have to set mV
        for(int j=0;j<Ny;j++){
            for(int i=0;i<Nx;i++){
                int jCount;
                double div;
 
                 if(i!=0 && i!=(Nx-1) && j!=0 && j !=(Ny-1)){
                 
                     div=(-2.*(1./pow(hx,2)+1./pow(hy,2))); //this will make the following divisions much more compact
                  
                     if(i>0){
                         jCount=j*Nx+(i-1);
                         mV[iCount]+=-V0[jCount]/pow(hx,2)/div;
                        
                     }
                     if(i<Nx-1){
                         jCount=j*Nx+(i+1);
                         mV[iCount]+=-V0[jCount]/pow(hx,2)/div;
  
                     }
                     if(j>0){
                         jCount=(j-1)*Nx+i;
                         mV[iCount]+=-V0[jCount]/pow(hy,2)/div;
                         
                    }
                     if(j<Ny-1){
                         jCount=(j+1)*Nx+i;
                          mV[iCount]+=-V0[jCount]/pow(hy,2)/div;
                         
                     }
                }
                 iCount++;

            }
        }
     //add results for mV to V1 and store mV in V0
        for(int i=0;i<Nx*Ny;i++){
            
            V1[i]+=mV[i];
            V0[i]=mV[i];
            
        }

        //set distance to zero and compute it to check if it's needed another iteration to reach threshold
        distance=0;
        for(int i=0;i<Nx*Ny;i++)
        	distance+=pow(V1[i]-V2[i],2);
        
        //we hade a ^2 result, the real distance si the squared root
        distance=sqrt(distance);
        //lets print it just to let the user know how long will it take to finish
        std::cout<<distance<<std::endl;
        //store V1 in V2
        for(int i=0;i<Nx*Ny;i++){
            
            V2[i]=V1[i];
        }
    }
    //at the end of the while cycle, V2 is the final potential that satisfies our needs. We can use it to compute electic field!! 
}

void ComputeElectricField(int Nx, int Ny,double* V1,double hx, double hy, double* Ex, double* Ey)
{
// just move step by step in x+ direction, computing for "vertical lines" Ey
	for (int i=0;i<(Nx);i++)
    {
    	
    	for  (int j=0;j<Ny; j++)
    	{
    	
    		if (i==0 || i==(Nx-1) || j==0 || j==(Ny-1)) 
                 Ey[Nx*i+j]=0;
                else
                 Ey[Nx*i+j]=((V1[(i+1)*Nx+j]-V1[(i-1)*Nx+j])/(2*hy));
	
        }
    }
    
    // lets compute symmetrically Ex, moving step by step in y+ direction and computing Ex for "horizontal lines"
     for (int i=0;i<(Ny);i++)
    {
    	
    	for  (int j=0;j<Nx; j++)
    	{
    	
    		if (i==0 || i==(Ny-1) || j==0 || j==(Nx-1)) 
                 Ex[Ny*i+j]=0;
                else
                 Ex[Ny*i+j]=((V1[(i)*Ny+(j+1)]-V1[(i)*Ny+(j-1)])/(2*hx));
	
        }
    }

}

