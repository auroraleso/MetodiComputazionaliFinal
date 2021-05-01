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
	bool accepted=0;
	while (accepted==0)
	{
		cout<<"Select an option, 3D plots:\n"<<
		"1) Set two charges;\n"<<
		"2) Set two conductors;\n"<<
		"3) Set two isolated spheres;\n";
		cout<<"Your choice (1,2,3): ";
		cin>>choice;
		if (choice==1)
		{
			charges();
			accepted=1;
		}
		else if (choice==2)
		{	
			conductors();
			accepted=1;
		}
		else if (choice==3)
		{	
			spheres();
			accepted=1;
		}
		else
		{
			//go back in the while, not allowed choice!
			cout<<"Not allowed request! Let me remind you...\n";
			accepted=0;
		}
	}
	
	

	return 0;
}

