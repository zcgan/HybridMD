#include "MDpara.h"
#include "Precisionsetting.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void read_para( )
{
	int i,j;
	double para[44];
	char name[400];
	FILE *input, *prec_header;
	double eps=0.001;
	double prec_setting[4][8]=prec_para;

	if((input=fopen(paraname,"r"))==NULL)
	{
		printf("can not find the input file!\n");
		exit (0);
	}
	i=0;
	while(!feof(input))		
	{ 
        if(i!=0&&i!=28)
            fscanf(input,"%s%lf",name,&para[i]);
        else if(i==0||i==28)
            fscanf(input,"%s",name);
        if(i!=0&&i!=28&&i<44)
            printf("%s is:  %f\n",name,para[i]);
        else if(i==0||i==28)
            printf("%s\n",name);		
        i++;
		

	}
	fclose(input);
	
	N=int(para[1]);
	Ntype=int(para[2]);
    if (Ntype>4)
    {
        printf("Sorry, the current code only allows at most 4 types:"
		" colloid1 col2 and ion1, ion2. Please change "
		"Ntype in the input file and rerun the code.\n");
        exit(0);
    }
	N_col=int(para[3]);
	N_ion=int(para[4]);
    if (N_col+N_ion!=N)
    {
        printf("Sorry, N_col+N_ion have to be equal to the total num of"
	" particles N, please fix the input file and rerun the code.\n");
        exit(0);
    }

    N_col1=int(para[5]);
    N_col2=int(para[6]);
    if (N_col1+N_col2!=N_col)
    {
        printf("Sorry, N_col1+N_col2 have to be equal to the total num of"
	" colloids N_col, please fix the input file and rerun the code.\n");
        exit(0);
    }
    
    N_ion1=int(para[7]);
    N_ion2=int(para[8]);
    if (N_ion1+N_ion2!=N_ion)
    {
        printf("Sorry, N_ion1+N_ion2 have to be equal to the total num of ions" 
		"N_ion, please fix the input file and rerun the code.\n");
        exit(0);
    }


    q_col1=para[9];
    q_col2=para[10];
    q_ion1=para[11];
    q_ion2=para[12];
   if (fabs(q_col1*N_col1+q_col2*N_col2+q_ion1*N_ion1+q_ion2*N_ion2)>eps)
    {
        printf("Sorry, the total charge nuetrality condition is not satisfied"
	", please fix the input file and rerun the code.\n");
        exit(0);
    }

    r_col1=para[13];
    r_col2=para[14];
    r_ion=para[15];
    c_lj=para[16];
    ei1=para[17];
    ei2=para[18];
    epsi_ion=para[19];
    m_col1=para[20];
    m_col2=para[21];
    m_ion=para[22];

    QM=para[23];
    RM=para[24];
    epsi_M=para[25];
    Rshell=para[26];
    temperature=para[27];
    KbT=temperature; /* set Kb=1.0 */

    dt_initial=para[29];
    dt_final=para[30];
    dt_inc_interval=int(para[31]);
    rescale_cycle=int(para[32]);
    FricCoef=para[33];
    Langevin_equi_steps=int(para[34]);
    N_lowtemp_annealing=int(para[35]);
    N_hightemp_annealing=int(para[36]);
    N_annealing_cycle=int(para[37]);
    Temprature_annealing=para[38];
    Production_steps=int(para[39]);
    sample_cycle=int(para[40]);
    bin_size=para[41];
    prec_set=int(para[42]);
    read_config=int(para[43]);
    strcpy(config_name,name);

if(read_config==1)
	printf("Reading the initial configuration from file %s\n",config_name);
	
	/*rescale size*/
    /*
	if(fabs(r_col-r_ion-1.0)>eps)
	{
		if(fabs(r_col-r_ion)>eps)
		{
		size_scaling=r_col-r_ion;
		r_col=r_col/size_scaling;
		r_ion=r_ion/size_scaling;
		RM=RM/size_scaling;
		Rshell=Rshell/size_scaling;
		}
	}*/
	
	time_stepratio=ceil(log(dt_final/dt_initial)/log(2.0));
	step_rescale=(time_stepratio+1)*dt_inc_interval;
	step_burnin=step_rescale+Langevin_equi_steps;
	step_max=step_burnin+Production_steps;
    
    if(r_col1>r_col2)
        r_col=r_col1;
    else
        r_col=r_col2;
    
	c_lj=c_lj*2.0;   //the input was half c_lj;
	Rc_lj1=2.0*r_col1+(factor_lj-1.0)*c_lj;
	Rc_lj2=2.0*r_col2+(factor_lj-1.0)*c_lj;
	
	
	/*read the header file for precision settings*/

	for(i=0;i<4;i++)
	{
		precision[i]=prec_setting[i][0];
		order_p[i]=prec_setting[i][1];
		im_num[i]=prec_setting[i][2];
		imm_num[i]=prec_setting[i][3];
		source_tol[i]=prec_setting[i][4];
		sph_tol[i]=prec_setting[i][5];
		gmres_tol[i]=prec_setting[i][6];
		fmm_tol[i]=	prec_setting[i][7];
		
	}
	
/*	if((prec_header=fopen("Precisionsetting.h","r"))==NULL)
	{
		printf("can not find the precision header file!\n");
		exit (0);
	}
	i=0;
	while(!feof(prec_header))		
	{ 
		if(i<4)
		{
		fscanf(prec_header,"%lf %lf %lf %lf %lf %lf %lf %lf\n",\
&precision[i],&order_p[i],&im_num[i],&imm_num[i],\
&source_tol[i],&sph_tol[i],&gmres_tol[i],&fmm_tol[i]);	
		}
		else
			return;
		i++;
		
	}
	fclose(prec_header);*/
}

