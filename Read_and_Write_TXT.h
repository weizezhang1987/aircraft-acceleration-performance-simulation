#include<iostream>
#include<fstream>
using namespace std;

class read_and_write_array  //sample read and write.
{
public:
	float temp_read[1000];  //temorary array to store data read from txt file
	int number_to_write;
	int number_to_read;
	string address;
public:
	float write_array(float array_to_write[])
	{
	ofstream in;
	in.open(address); 
	for (int m = 0; m < number_to_write; m++) 
	{
		in << array_to_write[m] << endl;
	}
	in.close();
	return true;
	}


	float read_array()
	{
	ifstream out;
	out.open(address);
	for (int m = 0; m < number_to_read; m++)
	{
		out >> temp_read[m];
	}
	out.close();
	return true;
	}

	/************
	for (i = 0; i < 1000; i++)
	{
		cout << input[i] <<" "<<out[i]<< endl;
	}
	*///////////
	//return 0;
};
