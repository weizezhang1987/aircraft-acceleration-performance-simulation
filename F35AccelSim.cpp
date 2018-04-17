// zwztest1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "math.h"
#include "RunKu.h"
#include "LinearRegression.h"
#include "Read_and_Write_TXT.h"
#include <vector>
#include <iostream>
#include <stdlib.h>     
#include <time.h>       
//using namespace std;

int main()
{	
	srand (time(NULL));
	// first R-K initial
	Rungue_Kutta R_K1;
	R_K1.timestep=0.005;
	float test_value_initial[10]={0};
	int i=1;

	for (i=1;i<10;i++)
	{
	test_value_initial[i]=1.000; //initial value 
	}

	for(R_K1.ii=1;R_K1.ii<10;R_K1.ii++)
	{
	R_K1.state_real[R_K1.ii]=test_value_initial[R_K1.ii];
	}
	// first R-K initial finishes

	// second T-K initial
	Rungue_Kutta_2 R_K2;
	R_K2.timestep=0.002;
	//float test_value_initial[10]={0};
	//int i=1;

	for (i=1;i<10;i++)
	{
	test_value_initial[i]=1.000; //initial value 
	}

	for(R_K2.ii=1;R_K2.ii<10;R_K2.ii++)
	{
	R_K2.state_real[R_K2.ii]=test_value_initial[R_K2.ii];
	}
	// second R-K initial finishes

	// first order sim initial
	Rungue_Kutta_first_order_sim R_K_1st_order;
	R_K_1st_order.timestep=0.002;

	R_K_1st_order.risingtime_max=0.29;
	R_K_1st_order.risingtime_min=0.15;
	R_K_1st_order.risingtime=R_K_1st_order.risingtime_min;

	R_K_1st_order.overshoot_max=0.1;
	R_K_1st_order.overshoot_min=0.05;
	R_K_1st_order.overshoot=R_K_1st_order.overshoot_min;

	for (i=1;i<10;i++)
	{
	test_value_initial[i]=0.000; //initial value 
	}

	for(R_K_1st_order.ii=1;R_K_1st_order.ii<10;R_K_1st_order.ii++)
	{
	R_K_1st_order.state_real[R_K_1st_order.ii]=test_value_initial[R_K_1st_order.ii];
	}

	ofstream file_1;
	file_1.open("D:\\read_and_write_test_2.txt");
	// first order sim initial finishes

	//********RunKu starts here********

	for (i=0;i<int(10.000/R_K_1st_order.timestep);i++) //Change iteration numbers here
	{
		/*
		for(R_K1.ii=1;R_K1.ii<10;R_K1.ii++)
	    {
	    R_K1.state_real[R_K1.ii]=update the input.
	    }
		*/

		// first R-K 
		R_K1.calculate_gradient();

		for (R_K1.ii=1;R_K1.ii<10;R_K1.ii++)
		{
			R_K1.state_predict[R_K1.ii]=R_K1.state_real[R_K1.ii]+R_K1.k_final[R_K1.ii]*R_K1.timestep;
		}

		for (R_K1.jj=1;R_K1.jj<10;R_K1.jj++)
		{
			R_K1.state_real[R_K1.jj]=R_K1.state_predict[R_K1.jj];
		}
		// first R-K finishes

		// second R-K 
		R_K2.input_real1=sin(3.1415926*(float(i)*R_K2.timestep));
		
		R_K2.calculate_gradient();

		for (R_K2.ii=1;R_K2.ii<10;R_K2.ii++)
		{
			R_K2.state_predict[R_K2.ii]=R_K2.state_real[R_K2.ii]+R_K2.k_final[R_K2.ii]*R_K2.timestep;
		}

		for (R_K2.jj=1;R_K2.jj<10;R_K2.jj++)
		{
			R_K2.state_real[R_K2.jj]=R_K2.state_predict[R_K2.jj];
		}
		// second R-K finishes

		// first order sim R-K
		if (i<int(1.00/R_K_1st_order.timestep))
		    R_K_1st_order.input_real1=0.25*1;
		if ((i>=int(1.00/R_K_1st_order.timestep))&&(i<int(2.00/R_K_1st_order.timestep)))
			R_K_1st_order.input_real1=0.25*2;
		if ((i>=int(2.00/R_K_1st_order.timestep))&&(i<int(3.00/R_K_1st_order.timestep)))
			R_K_1st_order.input_real1=0.25*3;
		if ((i>=int(3.00/R_K_1st_order.timestep))&&(i<int(4.00/R_K_1st_order.timestep)))
			R_K_1st_order.input_real1=0.25*4;
		if ((i>=int(4.00/R_K_1st_order.timestep))&&(i<int(5.00/R_K_1st_order.timestep)))
			R_K_1st_order.input_real1=0.25*3;
		if ((i>=int(5.00/R_K_1st_order.timestep))&&(i<int(6.00/R_K_1st_order.timestep)))
			R_K_1st_order.input_real1=0.25*2;
		if ((i>=int(6.00/R_K_1st_order.timestep))&&(i<int(7.00/R_K_1st_order.timestep)))
			R_K_1st_order.input_real1=0.25*1;
		if ((i>=int(7.00/R_K_1st_order.timestep))&&(i<int(8.00/R_K_1st_order.timestep)))
			R_K_1st_order.input_real1=0.25*0;

		if (int(10.0*float(i)*R_K_1st_order.timestep)%2==0)   //Trigger random every 0.2 second
			R_K_1st_order.rand_trigger=1;
		else R_K_1st_order.rand_trigger=0;   //Trigger no random 

		//if ((i>=int(1.00/R_K_1st_order.timestep))&&(i<int(2.00/R_K_1st_order.timestep)))   //Trigger at a specific step 
			//R_K_1st_order.rand_trigger=2; //Trigger override to fixed value
		
		R_K_1st_order.calculate_gradient();

		for (R_K_1st_order.ii=1;R_K_1st_order.ii<10;R_K_1st_order.ii++)
		{
			R_K_1st_order.state_predict[R_K_1st_order.ii]=R_K_1st_order.state_real[R_K_1st_order.ii]+R_K_1st_order.k_final[R_K_1st_order.ii]*R_K_1st_order.timestep;
		}

		for (R_K_1st_order.jj=1;R_K_1st_order.jj<10;R_K_1st_order.jj++)
		{
			R_K_1st_order.state_real[R_K_1st_order.jj]=R_K_1st_order.state_predict[R_K_1st_order.jj];
		}

		file_1<<float(i)*R_K_1st_order.timestep<<" "<<R_K_1st_order.input_real1<<" "<<R_K_1st_order.state_real[1]<<" "<<R_K_1st_order.state_real[1]+R_K_1st_order.noise<<endl;
		// first order sim R-K finishes

	}

	std::cout << R_K1.state_real[1] << " " << R_K1.state_real[2] << " " << R_K1.state_real[3] << " " << R_K1.state_real[4] << " " << R_K1.state_real[5] << " " << R_K1.state_real[6] << " " << R_K1.state_real[7] << " " << R_K1.state_real[8] << " " << R_K1.state_real[9] << endl;
	std::cout << R_K2.state_real[1] << " " << R_K2.state_real[2] << " " << R_K2.state_real[3] << " " << R_K2.state_real[4] << " " << R_K2.state_real[5] << " " << R_K2.state_real[6] << " " << R_K2.state_real[7] << " " << R_K2.state_real[8] << " " << R_K2.state_real[9] << endl;
	//************Runku finishes here********
	getchar();
	//************Linear regression here********
	vector<float> array_x(10);
	vector<float> array_y(10);
    // vector initialization
	int ii_vec = 0;
    for(ii_vec=0;ii_vec<10;ii_vec++) 
	{
	array_x[ii_vec]=ii_vec+1;
	array_y[ii_vec]=ii_vec+1+sin(0.8*ii_vec);
	}
	// vector initialization finishes
	LeastSquare ls(array_x, array_y);
	ls.print();
	//************Linear regression finishes********
	getchar();
	//************read and write here********
	float array_to_write_sample[]={0,1,2,3,4};
	
	read_and_write_array R_W_1;
	
	R_W_1.number_to_read=3;
	R_W_1.number_to_write=3;
	R_W_1.address="D:\\read_and_write_test_1.txt";
	
	R_W_1.write_array(array_to_write_sample);
	
	R_W_1.read_array();
	
	for (i=0;i<R_W_1.number_to_read;i++)
	{
		std::cout << array_to_write_sample[i] <<" "<<R_W_1.temp_read[i]<< endl;
	}

	if (R_W_1.temp_read[0]==0)
	{std::cout << "Matched" << endl;}
	//************read and write finishes********
	getchar();

	

}

