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
    s1.cx=4;
    s1.cy=2;
    s1.r=1;
    s1.charge=60;
    s2.cx=7;
    s2.cy=7;
    s2.r=1;
    s2.charge=-60;
    double Lx=10, Ly=10;
    int Nx=100, Ny=100;
    
    double thrs=0.1;
    double border=0.;
    double* phi0= new double[Nx*Ny];//funzione da trovare
    double* phi1= new double[Nx*Ny];//funzione da trovare
    double* phi2= new double[Nx*Ny];//funzione da trovare
    double* ftilde=new double[Nx*Ny];//distribuzione di carica e fattori bordo
    double*  mphi= new double[Nx*Ny];
    
    double hx=Lx/(Nx-1.);
    double hy=Ly/(Ny-1.);
    std::ofstream fphi,fcharge,felectric;
    
    
    ComputePotential( s1,  s2, phi0, phi1, phi2, ftilde, mphi,  hx,  hy,  Lx,  Ly,  Nx,  Ny,  thrs,  border);
    
  
    
    
    // Now, I need to compute Electric field as E=-grad(V):
     double* Ex= new double[Nx*Ny];
     double* Ey= new double[Ny*Nx];
    
      ComputeElectricField(Nx,Ny,phi1,hx,hy, Ex, Ey);
  
    
    
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
            
            fcharge << x << "  " << y << "  " << ftilde[iCount]*(-2.*(1./pow(hx,2)+1./pow(hy,2))) << '\n';
            fphi << x << "  " << y << "  " << -phi1[iCount] << '\n';
            felectric << x << "  " << y << "  " << Ex[iCount] << "  " <<Ey[iCount]<< '\n';
            iCount++;
        }
        fphi << '\n';
        fcharge << '\n';
    }
    
    fcharge.close();
    fphi.close();
    felectric.close();


    delete [] phi0;
    delete [] phi1;
    delete [] phi2;
    delete [] ftilde;
    delete [] mphi;
    delete [] Ex;
    delete [] Ey;

    }

void ComputePotential(sphere s1, sphere s2, double*phi0, double*phi1, double*phi2, double*ftilde, double*mphi, double hx, double hy, double Lx, double Ly, int Nx, int Ny, double thrs, double border )
{
	int dist_type=2;
	int iCount=0;
    for(int j=0;j<Ny;j++){
        double y;
        y=Ly*j/(Ny-1.);
        for(int i=0;i<Nx;i++){
            //std::cout<<i<<" "<<iCount<<std::endl;
            double x;
            x=Lx*i/(Nx-1.);
            //add charge distribution
        
                     switch(dist_type)
            {
                case 1:
                    if(sqrt(pow(x-s2.cx,2)+pow(y-s2.cy,2))<=s2.r) {
                        ftilde[iCount]=s2.charge;
                       
                    }else{
                        ftilde[iCount]=0.;
                    
                    }
                    break;

                case 2:
                    if(sqrt(pow(x-s1.cx,2)+pow(y-s1.cy,2))<=s1.r) {
                                    ftilde[iCount]=s1.charge;
                    }else if(sqrt(pow(x-s2.cx,2)+pow(y-s2.cy,2))<=s2.r){
                                    ftilde[iCount]=s2.charge;;
                    }else{
                        ftilde[iCount]=0.;
                                }
                            break;
	        
		
            }
                           
	        
		
            
            //metti il bordo
            if(i==0 || i==(Nx-1) || j==1 || j==(Ny-1)){
                ftilde[iCount]+=border;
                //std::cout<<"border\n";
            }
            //trasforma da f a f tilde
            if(i!=0 && i!=(Nx-1) && j!=0 && j !=(Ny-1)){
                ftilde[iCount]/=(-2.*(1./pow(hx,2)+1./pow(hy,2)));
                    //std::cout<<"inside\n";
            }
            
            iCount++;
                        
            
        }
    }
    
    
    //metti in phi0 condioni iniziali ossia pari a 0
    for(int i=0;i<Nx*Ny;i++){
        phi1[i]=ftilde[i];
        phi0[i]=ftilde[i];
        phi2[i]=0.;
      
    }

    double diff=1e10;
    while( diff>thrs){
         for(int i=0;i<Nx*Ny;i++){
             mphi[i]=0;
             
         }
        iCount=0;
        //aggiungi Mtilde/phi0
        for(int j=0;j<Ny;j++){
            for(int i=0;i<Nx;i++){
                int jCount;
                double div;
 
                 if(i!=0 && i!=(Nx-1) && j!=0 && j !=(Ny-1)){
                     div=(-2.*(1./pow(hx,2)+1./pow(hy,2)));
                   // std::cout<<"div "<<div<<std::endl;
                     if(i>0){
                         jCount=j*Nx+(i-1);
                         mphi[iCount]+=-phi0[jCount]/pow(hx,2)/div;
                          
                         
                     }
                     if(i<Nx-1){
                         jCount=j*Nx+(i+1);
                         mphi[iCount]+=-phi0[jCount]/pow(hx,2)/div;
                         
                        
                        
                     }
                     if(j>0){
                         jCount=(j-1)*Nx+i;
                         mphi[iCount]+=-phi0[jCount]/pow(hy,2)/div;
                         
                    }
                     if(j<Ny-1){
                         jCount=(j+1)*Nx+i;
                          mphi[iCount]+=-phi0[jCount]/pow(hy,2)/div;
                         
                     }
                }
                 iCount++;

            }
        }
     //metti in phi1
        for(int i=0;i<Nx*Ny;i++){
            //std::cout<<i<<") "<<mphi[i]<<std::endl;
            phi1[i]+=mphi[i];
            phi0[i]=mphi[i];
            
        }

        //calcola diffrenza
        diff=0;
        for(int i=0;i<Nx*Ny;i++){
           // std::cout<<"i:" << i<<std::endl;
            diff+=pow(phi1[i]-phi2[i],2);
            //std::cout << "Differenza: " << phi1[i] <<","<< phi2[i]<<"\n";
        }
        diff=sqrt(diff);
        std::cout<<diff<<std::endl;
        //copia phi0 in phi1
        for(int i=0;i<Nx*Ny;i++){
            
            phi2[i]=phi1[i];
        }
    }
 
    
    
}

void ComputeElectricField(int Nx, int Ny,double* phi1,double hx, double hy, double* Ex, double* Ey)
{
	for (int i=0;i<(Nx);i++)
    {
    	
    	for  (int j=0;j<Ny; j++)
    	{
    	
    		if (i==0 || i==(Nx-1) || j==0 || j==(Ny-1)) 
                 Ey[Nx*i+j]=0;
                else
                 Ey[Nx*i+j]=((phi1[(i+1)*Nx+j]-phi1[(i-1)*Nx+j])/(2*hy));
	
        }
    }
    
    
     for (int i=0;i<(Ny);i++)
    {
    	
    	for  (int j=0;j<Nx; j++)
    	{
    	
    		if (i==0 || i==(Ny-1) || j==0 || j==(Nx-1)) 
                 Ex[Ny*i+j]=0;
                else
                 Ex[Ny*i+j]=((phi1[(i)*Ny+(j+1)]-phi1[(i)*Ny+(j-1)])/(2*hx));
	
        }
    }

}

