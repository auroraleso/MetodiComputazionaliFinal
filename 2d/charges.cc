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

void ComputeElectricField(double*, double* , double* , int , int , double , double );
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
     	q1.q=10, q1.x=1.5, q1.y=2;
     	q2.q=-10, q2.x=2, q2.y=2;
     
     	
     	
     	ComputePotential(q1,q2, V, Nx, Ny, hx, hy,Lx,Ly);
     	ComputeElectricField(V, Ex, Ey, Nx, Ny, hx, hy);
     	
     	
     	
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
	for (int i=0; i<Nx; i++)
	{
		x=Lx*i/(Nx-1.);
		for (int j=0; j<Ny;j++)
		{
			y=Ly*j/(Ny-1.);
			
			if(sqrt(pow(x,2)+pow(y,2))!=sqrt(pow(q1.x,2)+pow(q1.y,2)) && sqrt(pow(x,2)+pow(y,2))!=sqrt(pow(q2.x,2)+pow(q2.y,2)) )
			
			{
				V1=k*(q1.q/(sqrt(pow((x-q1.x),2)+pow((y-q1.y),2))));
				V2=k*(q2.q/((sqrt(pow((x-q2.x),2)+pow((y-q2.y),2)))));
				V[i*Nx+j]=V1+V2;
			}
			
		}
	}
		
}
void ComputeElectricField(double*V, double* Ex, double* Ey, int Nx, int Ny, double hx, double hy)
{
	//E=-gradV
	for (int i=0;i<(Nx);i++)
    {
    	
    	for  (int j=0;j<Ny; j++)
    	{
    	
    		if (i==0 || i==(Nx-1) || j==0 || j==(Ny-1)) 
                 Ey[Nx*i+j]=0;
                else
                 Ey[Nx*i+j]=-((V[(i+1)*Nx+j]-V[(i-1)*Nx+j])/(2*hx));
	
        }
    }
    
    
     for (int i=0;i<(Ny);i++)
    {
    	
    	for  (int j=0;j<Nx; j++)
    	{
    	
    		if (i==0 || i==(Ny-1) || j==0 || j==(Nx-1)) 
                 Ex[Ny*i+j]=0;
                else
                 Ex[Ny*i+j]=-((V[(i)*Ny+(j+1)]-V[(i)*Ny+(j-1)])/(2*hy));
	
        }
    }
}







