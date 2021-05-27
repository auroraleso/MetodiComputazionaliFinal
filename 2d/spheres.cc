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
    s1.charge=10;
    s2.cx=7;
    s2.cy=7;
    s2.r=1;
    s2.charge=60;
    double Lx=10, Ly=10;
    int Nx=100, Ny=100;
    
    double threshold=0.1; //the value that needed to be reached in order to find the satisfying potential result. Clearly, the lower the threshold, the more precise will be the result (and the longer it will take). 
    double border=0.;
    //the followeing will help while computing potential 
    double* V0= new double[Nx*Ny];
    double* V1= new double[Nx*Ny];
    double* V2= new double[Nx*Ny];
    double* func=new double[Nx*Ny];//to charge and border distributions, its the "ftilde" according to explainations
    double*  MV= new double[Nx*Ny]; //is the product matrix * potential, (M*V)
    
    double hx=Lx/(Nx-1.);
    double hy=Ly/(Ny-1.);
    std::ofstream fV,fcharge,felectric;
    
    
    ComputePotential( s1,  s2, V0, V1, V2, func, MV,  hx,  hy,  Lx,  Ly,  Nx,  Ny,  threshold,  border);
    
  
    
    
    // Now, I need to compute Electric field as E=-grad(V):
     double* Ex= new double[Nx*Ny];
     double* Ey= new double[Ny*Nx];
    
      ComputeElectricField(Nx,Ny,V2,hx,hy, Ex, Ey);
  
    
    
    //let's save all data computed in txt so thate we'll be able to plot them
    fcharge.open("carica.txt");
    fV.open("V.txt");
    felectric.open("electricField.txt");
    //set precision for saving
    fcharge.precision(5);
    fV.precision(5);
    felectric.precision(5);
    int iCount=0;
    //write to txt in columns
    for(int j=0;j<Ny;j++){
        double y;
        y=Ly*j/(Ny-1.);
        for(int i=0;i<Nx;i++){
            double x;
            x=Lx*i/(Nx-1.);
            
            fcharge << x << "  " << y << "  " << func[iCount]*(-2.*(1./pow(hx,2)+1./pow(hy,2))) << '\n';
            fV << x << "  " << y << "  " << -V1[iCount] << '\n';
            felectric << x << "  " << y << "  " << Ex[iCount] << "  " <<Ey[iCount]<< '\n';
            iCount++;
        }
        fV << '\n';
        fcharge << '\n';
    }
    
    
    //close all opened files 
    fcharge.close();
    fV.close();
    felectric.close();

   //deleted the dynamical memory we allocated for arrays
    delete [] V0;
    delete [] V1;
    delete [] V2;
    delete [] func;
    delete [] MV;
    delete [] Ex;
    delete [] Ey;

    }
    
    
    
    
//function called to compute potential
void ComputePotential(sphere s1, sphere s2, double*V0, double*V1, double*V2, double*func, double*MV, double hx, double hy, double Lx, double Ly, int Nx, int Ny, double threshold, double border )
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
              
            
            //now that we set charge and border, we can compute f inside border by multiplying the setted value for -2(1/hx^2+1/hy^2)
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
    while( distance>threshold){
         for(int i=0;i<Nx*Ny;i++){
             MV[i]=0;
             
         }
        iCount=0;
        //we now have to set MV: according to the conditions we set for each case, we obtain different values
        for(int j=0;j<Ny;j++){
            for(int i=0;i<Nx;i++){
                int jCount;
                double div=(-2.*(1./pow(hx,2)+1./pow(hy,2)));//this will make the following divisions much more compact
 
                 if(i!=0 && i!=(Nx-1) && j!=0 && j !=(Ny-1))
                 {
		     if(i>0)
		     {
                         jCount=j*Nx+(i-1);
                         MV[iCount]+=-V0[jCount]/pow(hx,2)/div;
                        
                     }
                     if(i<Nx-1)
                     {
                         jCount=j*Nx+(i+1);
                         MV[iCount]+=-V0[jCount]/pow(hx,2)/div;
  
                     }
                     if(j>0)
                     {
                         jCount=(j-1)*Nx+i;
                         MV[iCount]+=-V0[jCount]/pow(hy,2)/div;
                         
                    }
                     if(j<Ny-1)
                     {
                         jCount=(j+1)*Nx+i;
                          MV[iCount]+=-V0[jCount]/pow(hy,2)/div;
                         
                     }
                }
                 iCount++;

            }
        }
     //add results for MV to V1 and store MV in V0
        for(int i=0;i<Nx*Ny;i++){
            
            V1[i]+=MV[i];
            V0[i]=MV[i];
            
        }

        //set distance to zero and compute it to check if it's needed another iteration to reach threshold
        distance=0;
        for(int i=0;i<Nx*Ny;i++)
        	distance+=pow(V1[i]-V2[i],2);
        
        //we hade a ^2 result, the real distance is the squared root
        distance=sqrt(distance);
        //lets print it just to let the user know how long will it take to finish
        std::cout<<distance<<std::endl;
        //store V1 in V2
        for(int i=0;i<Nx*Ny;i++){
            
            V2[i]=V1[i]; 
            
            /*doing this step, at the following iteration (Gauss-Seidel) we check if the difference between V^i(V2) and 			V^(i+1)(the new V1) is underneath the threshold, this means we reached convergence (Richardson) and we're
            	satisfied with the last vector obtained for V. Continuing to iterate is not catastrophic since 
                we'll produce the same vector again and again, but is not convenient in terms of time and memory 			       	used. */	 
        }
    }
    //at the end of the while cycle, V2 is the final potential that satisfies our needs. We can use it to compute electic field!! 
}



//function called to derive potential and obtain electric field
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

