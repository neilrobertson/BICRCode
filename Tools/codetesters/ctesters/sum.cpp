// file sum.cpp
#include <iostream>
#include <fstream>

int main()
{
   fstream ffile("f.dat");
   int nv;
   float sum=0;
   ffile >> nv; 
   for(int i=0; i<nv; i++)
   {  ffile >> fi;
      sum += fi;
  }
  cout << "sum=" << sum << endl;
}
