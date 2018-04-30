// zwztest1.cpp : Defines the entry point for the console application.
//

#include "math.h"
#include "RunKu.h"
#include "LinearRegression.h"
#include "Read_and_Write_TXT.h"
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <time.h>
//using namespace std;

struct A_T1_combination
{
    double A = 0.0;
    double T1 = 0.0;
};

int main()
{
	// F-35 accel initial
    double g = 9.8;
	double m1 = 19000;
	double m2 = 15600;
	double density1 = 0.774;
	double density2 = 1.111;
	double cla = 0.07;
	double cd0 = 0.02;
	double S = 42.7;
	double A = 0.15; // not accurate. to be corrected
	double T1 = 13000 * 9.8; // to be corrected
	double T2 = 18000 * 9.8; // to be corrected
    std::vector<A_T1_combination> possible_combination;

	Rungue_Kutta_simple_integration F35_accel;
	F35_accel.timestep = 0.1;
	float F35_initial[10] = {0.0};
	F35_initial[1] = 193.3;

	for(F35_accel.ii = 1; F35_accel.ii < 10; F35_accel.ii++)
	{
		F35_accel.state_real[F35_accel.ii] = F35_initial[F35_accel.ii];
	}
	// F-35 accel initial finishes

	double A_table[20] = {0.0};
	double T1_table[20] = {0.0};
	double v_final_table[20][20] = {{0.0}};
	for (int i = 0; i < 20; i++)
	{
		A_table[i] = 0.05 + 0.01 * i;
		T1_table[i] = 10000 * g + 500 * g * i;
	}

	std::cout << "Official data: at 4527 m, from 193.3m/s and accelerate for 17.9 s, F-35 should achieve 306 m/s" << std::endl;
	std::cout << "Official data verification at 4527 m, after accelerating for 17.9 s:" << std::endl;
	for (int j = 0; j < 20; j++)
	{
		for (int k = 0; k < 20; k++)
		{
			A = A_table[j];
			T1 = T1_table[k];

			for(F35_accel.ii = 1; F35_accel.ii < 10; F35_accel.ii++)
            {
                F35_accel.state_real[F35_accel.ii] = F35_initial[F35_accel.ii];
            }

			for (int i = 0; i < int(17.9 / F35_accel.timestep); i++) // Accelerate until 17.9 s
			{
				double cl = m1 * g / (0.5 * density1 * F35_accel.state_real[1] * F35_accel.state_real[1] * S);
				double cd = cd0 + A * cl * cl;
				double D = cd / cl * m1 * g;
				F35_accel.input_real1 = 1.0 / m1 * (T1 - D); // input_real1 is the only external parameter that is used to calculate gradient.
				F35_accel.calculate_gradient();

				for (F35_accel.ii = 1; F35_accel.ii < 10; F35_accel.ii++)
				{
					F35_accel.state_predict[F35_accel.ii] = F35_accel.state_real[F35_accel.ii] + F35_accel.k_final[F35_accel.ii] * F35_accel.timestep;
				}

				for (F35_accel.jj = 1; F35_accel.jj < 10; F35_accel.jj++)
				{
					F35_accel.state_real[F35_accel.jj] = F35_accel.state_predict[F35_accel.jj];
				}
			}
			v_final_table[j][k] = F35_accel.state_real[1];
			if (v_final_table[j][k] < 308.0 && v_final_table[j][k] > 304.0)
            {
                std::cout << "A: " << A_table[j] << " T1: " << T1_table[k] << " v_final: " << v_final_table[j][k] << std::endl;
                A_T1_combination combo_temp;
                combo_temp.A = A_table[j];
                combo_temp.T1 = T1_table[k];
                possible_combination.push_back(combo_temp);
            }
		}
	}

	// verification finished. Now initialize prediction.
	std::cout << "////// ****** //////" << std::endl;
	std::cout << "F35 sim v3.02 at 1000 m, 193.3m/s to 306m/s acceleration:" << std::endl;
	F35_initial[1] = 193.3;
	for (int i = 0; i < possible_combination.size(); i++)
    {
        A = possible_combination[i].A;
        T1 = possible_combination[i].T1;
        T2 = density2 / density1 * T1;

        for(F35_accel.ii = 1; F35_accel.ii < 10; F35_accel.ii++)
        {
            F35_accel.state_real[F35_accel.ii] = F35_initial[F35_accel.ii];
        }

        for (int j = 0; j < int(20.0 / F35_accel.timestep); j++)
        {
            double cl = m2 * g / (0.5 * density2 * F35_accel.state_real[1] * F35_accel.state_real[1] * S);
            double cd = cd0 + A * cl * cl;
            double D = cd / cl * m2 * g;
            F35_accel.input_real1 = 1.0 / m2 * (T2 - D); // input_real1 is the only external parameter that is used to calculate gradient.
            F35_accel.calculate_gradient();

            for (F35_accel.ii = 1; F35_accel.ii < 10; F35_accel.ii++)
            {
                F35_accel.state_predict[F35_accel.ii] = F35_accel.state_real[F35_accel.ii] + F35_accel.k_final[F35_accel.ii] * F35_accel.timestep;
            }

            for (F35_accel.jj = 1; F35_accel.jj < 10; F35_accel.jj++)
            {
                F35_accel.state_real[F35_accel.jj] = F35_accel.state_predict[F35_accel.jj];
            }

            if (F35_accel.state_real[1] > 306)
            {
                std::cout << "A: " << A << " T2: " << T2 << " v_ini: " << F35_initial[1] << " v_ final: " << F35_accel.state_real[1] << " time: " << j * F35_accel.timestep << std::endl;
                break;
            }
        }
    }

    std::cout << "////// ****** //////" << std::endl;
    std::cout << "F35 sim v3.02 at 1000 m, 166.6m/s (600km/h) to 305.5m/s (1100km/h) acceleration:" << std::endl;
	F35_initial[1] = 166.6;
	for (int i = 0; i < possible_combination.size(); i++)
    {
        A = possible_combination[i].A;
        T1 = possible_combination[i].T1;
        T2 = density2 / density1 * T1;

        for(F35_accel.ii = 1; F35_accel.ii < 10; F35_accel.ii++)
        {
            F35_accel.state_real[F35_accel.ii] = F35_initial[F35_accel.ii];
        }

        for (int j = 0; j < int(20.0 / F35_accel.timestep); j++)
        {
            double cl = m2 * g / (0.5 * density2 * F35_accel.state_real[1] * F35_accel.state_real[1] * S);
            double cd = cd0 + A * cl * cl;
            double D = cd / cl * m2 * g;
            F35_accel.input_real1 = 1.0 / m2 * (T2 - D); // input_real1 is the only external parameter that is used to calculate gradient.
            F35_accel.calculate_gradient();

            for (F35_accel.ii = 1; F35_accel.ii < 10; F35_accel.ii++)
            {
                F35_accel.state_predict[F35_accel.ii] = F35_accel.state_real[F35_accel.ii] + F35_accel.k_final[F35_accel.ii] * F35_accel.timestep;
            }

            for (F35_accel.jj = 1; F35_accel.jj < 10; F35_accel.jj++)
            {
                F35_accel.state_real[F35_accel.jj] = F35_accel.state_predict[F35_accel.jj];
            }

            if (F35_accel.state_real[1] > 305.5)
            {
                std::cout << "A: " << A << " T2: " << T2 << " v_ini: " << F35_initial[1] << " v_ final: " << F35_accel.state_real[1] << " time: " << j * F35_accel.timestep << std::endl;
                break;
            }
        }
    }

}

