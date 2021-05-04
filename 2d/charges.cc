#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

//1/(4*pi*epsilon0)=k
 long int k=8987551788;


struct charge
{
	double q;
	double x;
	double y;
};
void ComputePotential(charge , charge , double* , int , int , int , int , double , double );

void ComputeElectricField(charge, charge, double*, double* , double* , int , int , double , double, double, double );
void charges()
{
	
	std::cout<<"You have chosen CHARGES!\n";
	
	
	double Lx=3, Ly=3;
    	int Nx=100, Ny=100;
    	double hx=(Lx)/(Nx-1.), hy=Ly/(Nx-1.);
    	
    	
    	
	double* V= new double [Nx*Ny];
	double* Ex= new double[Nx*Ny];
     	double* Ey= new double[Ny*Nx];
     	
     	
     	charge q1, q2;
     	q1.q=10*pow(10,-6), q1.x=1.5, q1.y=2;
     	q2.q=-10*pow(10,-6), q2.x=1.6, q2.y=2;
     
     	
     	
     	ComputePotential(q1,q2, V, Nx, Ny, hx, hy,Lx,Ly);
     	ComputeElectricField(q1,q2,V, Ex, Ey, Nx, Ny, hx, hy, Lx, Ly);
     	
     	
     	
    std::ofstream fV, fElectric;	
    //fcharge.open("carica.txt");
    fV.open("Vcharges.txt");
    fElectric.open("electricFieldcharges.txt");
    //fcharge.precision(5);
    fV.precision(5);
    fElectric.precision(5);
    int iCount=0;
    for(int j=0;j<Ny;j++){
        double y;
        y=Ly*j/(Ny-1.);
        for(int i=0;i<Nx;i++){
            double x;
            x=Lx*i/(Nx-1.);
            
           // fcharge << x << "  " << y << "  " << V[iCount]*(-2.*(1./pow(hx,2)+1./pow(hy,2))) << '\n';
            fV << x << "  " << y << "  " << V[iCount] << '\n';
            fElectric << x << "  " << y << "  " << Ex[iCount] << "  " <<Ey[iCount]<< '\n';
            iCount++;
        }
        fV << '\n';
       // fcharge << '\n';
    }
    
   // fcharge.close();
    fV.close();
    fElectric.close();
    
    
     	
     	delete[] Ex;
     	delete[] Ey;
     	delete[] V;
     	
}
void ComputePotential(charge q1, charge q2, double* V, int Nx, int Ny, int hx, int hy, double Lx, double Ly)
{
	double x,y,V1,V2;
	for (int i=0; i<Ny; i++)
	{
		y=Ly*i/(Ny-1.);
		for (int j=0; j<Nx;j++)
		{
			x=Lx*j/(Nx-1.);
			
			if(sqrt(pow(x,2)+pow(y,2))!=sqrt(pow(q1.x,2)+pow(q1.y,2)) && sqrt(pow(x,2)+pow(y,2))!=sqrt(pow(q2.x,2)+pow(q2.y,2)) )
			
			{
				V1=1*(q1.q/(sqrt(pow((x-q1.x),2)+pow((y-q1.y),2))));
				V2=1*(q2.q/((sqrt(pow((x-q2.x),2)+pow((y-q2.y),2)))));
				V[i*Ny+j]=V1+V2;
			}
			else
			{
				V[i*Ny+j]=0;
			}
			
		}
	}
		
}
void ComputeElectricField(charge q1, charge q2, double*V, double* Ex, double* Ey, int Nx, int Ny, double hx, double hy, double Lx, double Ly)
{
	//E=-gradV
	for (int i=0;i<(Ny);i++)
    {
    	
    	for  (int j=0;j<Nx; j++)
    	{
    	
    		if (i==0 || i==(Nx-1) || j==0 || j==(Ny-1) || V[Ny*i+j]==0) 
                 Ey[Ny*i+j]=0;
                else
                 Ey[Ny*i+j]=-((V[(i+1)*Ny+j]-V[(i-1)*Ny+j])/(2*hy));
	
        }
    }
    
    
     for (int i=0;i<(Nx);i++)
    {
    	
    	for  (int j=0;j<Ny; j++)
    	{
    	
    		if (i==0 || i==(Ny-1) || j==0 || j==(Nx-1)) 
                 Ex[Nx*i+j]=0;
                else
                 Ex[Nx*i+j]=-((V[(i)*Nx+(j+1)]-V[(i)*Nx+(j-1)])/(2*hx));
	
        }
    }
    
    
    //otherwise
    
   /* double x,y,E1,E2;
    double* E=new double [Nx*Ny];
	for (int i=0; i<Nx; i++)
	{
		x=Lx*i/(Nx-1.);
		for (int j=0; j<Ny;j++)
		{
			y=Ly*j/(Ny-1.);
			
			if(sqrt(pow(x,2)+pow(y,2))!=sqrt(pow(q1.x,2)+pow(q1.y,2)) && sqrt(pow(x,2)+pow(y,2))!=sqrt(pow(q2.x,2)+pow(q2.y,2)) )
			
			{
				E1=k*(q1.q/pow((y-q1.y),2)+pow(x-q1.x,2));
				E2=k*(q2.q/pow((y-q2.y),2)+pow(x-q2.x,2));
				E[i*Nx+j]=E1+E2;
				Ex[i*Nx+j]=
			}
			
		}
	}*/
	

    
}







