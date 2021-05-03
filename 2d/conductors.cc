#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

struct conductor
{
	double cx;
	double cy;
	double r;
	double q;
};

void ComputePotentialC(conductor , conductor , double* , int , int , int , int , double , double );
void ComputeElectricFieldC(conductor, conductor, double*, double* , double* , int , int , double , double, double, double );


void conductors ()
{
	double Lx=10, Ly=10;
    	int Nx=100, Ny=100;
    	double hx=(Lx)/(Nx-1.), hy=Ly/(Ny-1.);
    	
	std::cout<<"You have chosen CONDUCTORS!\n";
	conductor c1;
	conductor c2;
	//setting values 
	c1.cx=3, c1.cy=5, c1.r=1, c1.q=80;
	c2.cx=8, c2.cy=2, c2.r=1, c2.q=-80;
	
	
	double* V= new double [Nx*Ny];
	double* Ex= new double[Nx*Ny];
     	double* Ey= new double[Ny*Nx];
     	
     	
     	
     	ComputePotentialC(c1,c2, V, Nx, Ny, hx, hy,Lx,Ly);
     	ComputeElectricFieldC(c1,c2,V, Ex, Ey, Nx, Ny, hx, hy, Lx, Ly);
     	
     	std::ofstream fV, fElectric;	
    //fcharge.open("carica.txt");
    fV.open("Vconductors.txt");
    fElectric.open("electricFieldconductors.txt");
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



void ComputePotentialC(conductor c1, conductor c2, double* V, int Nx, int Ny, int hx, int hy, double Lx, double Ly)
{
	double x,y;
	for (int i=0; i<Nx; i++)
	{
		x=Lx*i/(Nx-1.);
		for (int j=0; j<Ny;j++)
		{
			y=Ly*j/(Ny-1.);
			
			if(sqrt(pow(x-c1.cx,2)+pow(y-c1.cy,2))>=c1.r && sqrt(pow(x-c2.cx,2)+pow(y-c2.cy,2))>=c2.r)
			
			{
				V[Nx*i+j]=8987551788*(c1.q/(sqrt(pow(x-c1.cx,2)+pow(y-c1.cy,2)))+c2.q/(pow(x-c2.cx,2)+pow(y-c2.cy,2)));
			}
			else if (sqrt(pow(x-c1.cx,2)+pow(y-c1.cy,2))<=c1.r)
			{
			V[Nx*i+j]=8987551788*(c1.q/(c1.r));
			}
			else if (sqrt(pow(x-c2.cx,2)+pow(y-c2.cy,2))<=c2.r)
			{
			V[Nx*i+j]=8987551788*(c2.q/(c2.r));
			}
			
		}
	}
		
}
void ComputeElectricFieldC(conductor c1, conductor c2,double*V,  double* Ex, double* Ey, int Nx, int Ny, double hx, double hy, double Lx, double Ly)
{
	//E=-gradV
	double x,y;
	for (int i=0; i<Nx; i++)
	{
		x=Lx*i/(Nx-1.);
		for (int j=0; j<Ny;j++)
		{
			y=Ly*j/(Ny-1.);
    	
    		if (i==0 || i==(Nx-1) || j==0 || j==(Ny-1)) 
                 Ey[Nx*i+j]=0;
                else if (sqrt(pow(x-c1.cx,2)+pow(y-c1.cy,2))<=c1.r || sqrt(pow(x-c2.cx,2)+pow(y-c2.cy,2))<=c2.r)
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
                else if (sqrt(pow(x-c1.cx,2)+pow(y-c1.cy,2))<=c1.r || sqrt(pow(x-c2.cx,2)+pow(y-c2.cy,2))<=c2.r)
                	Ex[Ny*i+j]=0;
                else
                 Ex[Ny*i+j]=-((V[(i)*Ny+(j+1)]-V[(i)*Ny+(j-1)])/(2*hy));
	
        }
    }
}




