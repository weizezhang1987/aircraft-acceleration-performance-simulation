#include <complex>
#include <math.h>
using namespace std;

class Rungue_Kutta  //sample Rungue_Kutta
{
public:
	
	int dim_state;
	int i;
	int ii;
	int iii;
	int j;
	int jj;
	int jjj;

	float timestep;

	float k1[10];
	float k2[10];
	float k3[10];
	float k4[10];
	float k_final[10];
	
	//float input1,input2,input3,input4,input5,input6;
	//float state[5];
	float state_real[10];  //current state
	float state_predict[10];  //predicted future state
	float state_temp[10];  //temporary state used in gradient calculation
	float state_dot[10];   //return value of differential equation
	float input_real1,input_real2,input_real3,input_real4,input_real5,input_real6;  //possible additional parameters, usually unused

	float *state_dot_func(float state[10],float input1,float input2,float input3,float input4,float input5,float input6)
	{
	state_dot[1]=state[1];
	state_dot[2]=state[2];
	state_dot[3]=state[3];
	state_dot[4]=state[4];
	state_dot[5]=state[5];
	state_dot[6]=state[6];
	state_dot[7]=state[7];
	state_dot[8]=state[8];
	state_dot[9]=state[9];
	return state_dot;
	}

	float *calculate_gradient()
	{
	i=1;
    dim_state=10;

	state_dot_func(state_real,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k1[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k1[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k2[i]=state_dot[i];}
	//*****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k2[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k3[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k3[i]*timestep;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k4[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{k_final[i]=(k1[i]+2.00*k2[i]+2.00*k3[i]+k4[i])/6.00;}
	
	return k_final;

	}

};

class Rungue_Kutta_2  //sample Rungue_Kutta, second one, a complex differentiator
{
public:
	
	int dim_state;
	int i;
	int ii;
	int iii;
	int j;
	int jj;
	int jjj;

	float timestep;

	float k1[10];
	float k2[10];
	float k3[10];
	float k4[10];
	float k_final[10];
	
	//float input1,input2,input3,input4,input5,input6;
	//float state[5];
	float state_real[10];  //current state
	float state_predict[10];  //predicted future state
	float state_temp[10];  //temporary state used in gradient calculation
	float state_dot[10];   //return value of differential equation
	float input_real1,input_real2,input_real3,input_real4,input_real5,input_real6;  //possible additional parameters, usually unused

	//******differential equation intial
	float epsilon;

    float k11;
    float k12;
    float k21;
    float k22;
    float k31;
    float k32;

	float z1_dot,z2_dot,z1_est_dot,z2_est_dot,x_est_dot;
	float input,z1,z2,z1_est,z2_est,x_est;
	//******differential equation initialization finishes

	float *state_dot_func(float state[10],float input1,float input2,float input3,float input4,float input5,float input6)
	{
	epsilon=1.00;

    k11=-13.0;
    k12=-13.0;
    k21=-13.0;
    k22=-13.0;
    k31=-13.0;
    k32=-169.0;
	
	z1=state[1];
	z2=state[2];
	z1_est=state[3];
	z2_est=state[4];
	x_est=state[5];

    z1_dot=1.0/epsilon*z2;
    z2_dot=1.0/epsilon*input1;
    z1_est_dot=1.0/epsilon*(z2_est+k11*(z1_est-z1)+k12*(z2_est-z2));
    z2_est_dot=1.0/epsilon*(x_est+k21*(z1_est-z1)+k22*(z2_est-z2));
    x_est_dot=1.0/epsilon*(k31*(z1_est-z1)+k32*(z2_est-z2));

	state_dot[1]=z1_dot;
	state_dot[2]=z2_dot;
	state_dot[3]=z1_est_dot;
	state_dot[4]=z2_est_dot;
	state_dot[5]=x_est_dot;

    state_dot[6]=0;
	state_dot[7]=0;
	state_dot[8]=0;
	state_dot[9]=0;

	return state_dot;
	}

	float *calculate_gradient()
	{
	i=1;
    dim_state=10;

	state_dot_func(state_real,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k1[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k1[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k2[i]=state_dot[i];}
	//*****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k2[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k3[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k3[i]*timestep;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k4[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{k_final[i]=(k1[i]+2.00*k2[i]+2.00*k3[i]+k4[i])/6.00;}
	
	return k_final;

	}

};

class Rungue_Kutta_simple_integration  //sample imtegration of a time-variant value in R-K
{
public:
	
	int dim_state;
	int i;
	int ii;
	int iii;
	int j;
	int jj;
	int jjj;

	float timestep;

	float k1[10];
	float k2[10];
	float k3[10];
	float k4[10];
	float k_final[10];
	
	//float input1,input2,input3,input4,input5,input6;
	//float state[5];
	float state_real[10];  //current state
	float state_predict[10];  //predicted future state
	float state_temp[10];  //temporary state used in gradient calculation
	float state_dot[10];   //return value of differential equation
	float input_real1,input_real2,input_real3,input_real4,input_real5,input_real6,input_real7,input_real8,input_real9;  //possible additional parameters, usually unused

	float *state_dot_func(float state[10],float input1,float input2,float input3,float input4,float input5,float input6,float input7,float input8,float input9)
	{
	state_dot[1]=input1;
	state_dot[2]=input2;
	state_dot[3]=input3;
	state_dot[4]=input4;
	state_dot[5]=input5;
	state_dot[6]=input6;
	state_dot[7]=input7;
	state_dot[8]=input8;
	state_dot[9]=input9;
	return state_dot;
	}

	float *calculate_gradient()
	{
	i=1;
    dim_state=10;

	state_dot_func(state_real,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6,input_real7,input_real8,input_real9);
	for (i=1;i<dim_state;i++)
	{k1[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k1[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6,input_real7,input_real8,input_real9);
	for (i=1;i<dim_state;i++)
	{k2[i]=state_dot[i];}
	//*****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k2[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6,input_real7,input_real8,input_real9);
	for (i=1;i<dim_state;i++)
	{k3[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k3[i]*timestep;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6,input_real7,input_real8,input_real9);
	for (i=1;i<dim_state;i++)
	{k4[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{k_final[i]=(k1[i]+2.00*k2[i]+2.00*k3[i]+k4[i])/6.00;}
	
	return k_final;

	}

};

class Rungue_Kutta_first_order_sim  //first order simulator
{
public:
	//dynamic parameter settings
	float noise;
	int rand_trigger;
	float overshoot,overshoot_max,overshoot_min,risingtime,risingtime_max,risingtime_min,square;// overshoot 0-0.5, rising_time in seconds, square is a temporary value
	float psy,omega_n;
	//dynamic parameter settings fnished
	int dim_state;
	int i;
	int ii;
	int iii;
	int j;
	int jj;
	int jjj;

	float timestep;

	float k1[10];
	float k2[10];
	float k3[10];
	float k4[10];
	float k_final[10];
	
	//float input1,input2,input3,input4,input5,input6;
	//float state[5];
	float state_real[10];  //current state
	float state_predict[10];  //predicted future state
	float state_temp[10];  //temporary state used in gradient calculation
	float state_dot[10];   //return value of differential equation
	float input_real1,input_real2,input_real3,input_real4,input_real5,input_real6;  //possible additional parameters, usually unused

	float *state_dot_func(float state[10],float input1,float input2,float input3,float input4,float input5,float input6)
	{
    //calculate dynamic parameters
		if (rand_trigger==1)
		{
	    overshoot=overshoot_min+(rand()/float(RAND_MAX))*(overshoot_max-overshoot_min);
	    risingtime=risingtime_min+(rand()/float(RAND_MAX))*(risingtime_max-risingtime_min);
		}
		if (rand_trigger==2)  // trigger override
		{overshoot=0.35;}
	psy=(overshoot-0.5625)/(-0.625);
	omega_n=4.0/(psy*risingtime);
    noise=((rand()/float(RAND_MAX))-0.5)*2.0*0.032; //noise level to 0.5, 0.4...
    // dynamic parameters finished
	state_dot[1]=state[2];
	state_dot[2]=-omega_n*omega_n*state[1]-2*psy*omega_n*state[2]+omega_n*omega_n*input1;
	state_dot[3]=0;
	state_dot[4]=0;
	state_dot[5]=0;
	state_dot[6]=0;
	state_dot[7]=0;
	state_dot[8]=0;
	state_dot[9]=0;
	return state_dot;
	}

	float *calculate_gradient()
	{
	i=1;
    dim_state=10;

	state_dot_func(state_real,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k1[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k1[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k2[i]=state_dot[i];}
	//*****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k2[i]*timestep/2.00;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k3[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{state_temp[i]=state_real[i]+k3[i]*timestep;}
	state_dot_func(state_temp,input_real1,input_real2,input_real3,input_real4,input_real5,input_real6);
	for (i=1;i<dim_state;i++)
	{k4[i]=state_dot[i];}
	//****************************
	for (i=1;i<dim_state;i++)
	{k_final[i]=(k1[i]+2.00*k2[i]+2.00*k3[i]+k4[i])/6.00;}
	
	return k_final;

	}

};