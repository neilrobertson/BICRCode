#include <iostream>
#include <fstream>

using namespace std;

int main()
{
	string temp;
     ifstream fin( "/home/pzs/histone/HISTONE_DATA/raw_affy_data/K562_H3K36me1_Normtogether_1_signal.txt-processed.txt" );
    ofstream fout( "/tmp/results.txt", ios::out );
    while (! fin.eof() )
	{
      getline(fin, temp);
	  fout << temp  << endl;
	}

}
