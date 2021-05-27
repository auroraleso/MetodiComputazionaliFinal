#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
//#include "charges.cc"
//#include "conductors.cc"
//#include "spheres.cc"
void charges();
void conductors();
void spheres();
using namespace std;
int main()
{
	int choice;
	bool accepted=false;
	while (!accepted)
	{
		cout<<"Select an option:\n"<<
		"1) Set two charges;\n"<<
		"2) Set two conductors;\n"<<
		"3) Set two isolated spheres;\n";
		cout<<"Your choice (1,2,3): ";
		cin>>choice;
		if (choice==1)
		{
			charges();
			accepted=true;
		}
		else if (choice==2)
		{	
			conductors();
			accepted=true;
		}
		else if (choice==3)
		{	
			spheres();
			accepted=true;
		}
		else
		{
			//go back in the while, not allowed choice!
			cout<<"Not allowed request! Let me remind you...\n";
			
		}
	}
	
	

	return 0;
}

