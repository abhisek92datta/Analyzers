#include<iostream>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<string.h>
#include<stdlib.h>

using namespace std;

int main()
{

	ifstream fin;
	ofstream fout, fout_plot;

	int col;
	double ulim,dlim,nbins;
	double bin_size;
	char xlabel[100];
	
	char awk_comm[100];
	char file1[100],file2[100];
	
	int N1,N2;
	int is_SL,is_DL;
	int n_sl_mc, n_dl_mc, n_sl_data, n_dl_data;
	
	double SL_mc_var[10000], DL_mc_var[10000];
	double SL_data_var[10000], DL_data_var[10000];
	double temp;
	double x[200],y_sl_mc[200],y_sl_data[200],y_dl_mc[200],y_dl_data[200];
	double sl_data_err[200], dl_data_err[200];
	
	sprintf(file1,"mc.csv");
	sprintf(file2,"data.csv");
	
	cin>>col;
	cin>>dlim>>ulim>>nbins;	
	cin>>xlabel;

	sprintf(awk_comm, "awk -v lc=\"%d\" '{ if(FNR > 1){print $4+$5, $6+$7+$8, $lc } }' %s > temp_mc.txt", col,file1);
	system(awk_comm);
	sprintf(awk_comm, "awk -v lc=\"%d\" '{ if(FNR > 1){print $4+$5, $6+$7+$8, $lc } }' %s > temp_data.txt", col,file2);
	system(awk_comm);
	
	
	N1 = N2 = 0;
	n_sl_mc = n_dl_mc = n_sl_data = n_dl_data = 0;
	
	//MC
	fin.open("temp_mc.txt");
	while(!fin.eof())
	{
		fin>>is_SL>>is_DL;
		if(is_SL == 1)
		{
			fin>>SL_mc_var[n_sl_mc];
			n_sl_mc++;
		}
		else 
		{
			fin>>DL_mc_var[n_dl_mc];
			n_dl_mc++;
		}			
		N1++;
	}
	N1--;
	if(is_SL==1)
		n_sl_mc--;
	else
		n_dl_mc--;
	fin.close();
	
	//DATA
	fin.open("temp_data.txt");
	while(!fin.eof())
	{
		fin>>is_SL>>is_DL;
		if(is_SL == 1)
		{
			fin>>SL_data_var[n_sl_data];
			n_sl_data++;
		}
		else 
		{
			fin>>DL_data_var[n_dl_data];
			n_dl_data++;
		}			
		N2++;
	}
	N2--;
	if(is_SL==1)
		n_sl_data--;
	else
		n_dl_data--;
	fin.close();

	bin_size = (ulim-dlim)/nbins;
	for(int i=0; i<nbins; i++)
	{			
		x[i] = dlim + (i+0.5)*bin_size;
		y_sl_mc[i] = y_sl_data[i] = y_dl_mc[i] = y_dl_data[i] = 0;
	}
	
	//SL MC
	for(int j=0; j<n_sl_mc; j++)
	{
		for(int i=0; i<nbins; i++)
		{			
			temp = dlim + (i+1)*bin_size;
			if(SL_mc_var[j]<=temp)
			{
				y_sl_mc[i]++;
				break;
			}
		}
	}

	// SL DATA
	for(int j=0; j<n_sl_data; j++)
	{
		for(int i=0; i<nbins; i++)
		{			
			temp = dlim + (i+1)*bin_size;
			if(SL_data_var[j]<=temp)
			{
				y_sl_data[i]++;
				break;
			}
		}
	}

	//DL MC
	for(int j=0; j<n_dl_mc; j++)
	{
		for(int i=0; i<nbins; i++)
		{			
			temp = dlim + (i+1)*bin_size;
			if(DL_mc_var[j]<=temp)
			{
				y_dl_mc[i]++;
				break;
			}
		}
	}

	//DL DATA
	for(int j=0; j<n_dl_data; j++)
	{
		for(int i=0; i<nbins; i++)
		{			
			temp = dlim + (i+1)*bin_size;
			if(DL_data_var[j]<=temp)
			{
				y_dl_data[i]++;
				break;
			}
		}
	}

	fout.open("dist_hist.txt");
	for(int i=0; i<nbins; i++)
	{
		y_sl_mc[i] = y_sl_mc[i]/n_sl_mc;
		y_sl_data[i] = y_sl_data[i]/n_sl_data;
		y_dl_mc[i] = y_dl_mc[i]/n_dl_mc;
		y_dl_data[i] = y_dl_data[i]/n_dl_data;
		sl_data_err[i] = (sqrt(y_sl_data[i])/n_sl_data) + (y_sl_data[i]/(n_sl_data*sqrt(n_sl_data)));
		dl_data_err[i] = (sqrt(y_dl_data[i])/n_dl_data) + (y_dl_data[i]/(n_dl_data*sqrt(n_dl_data)));
		
		fout<<x[i]<<"   "<<y_sl_mc[i]<<"   "<<y_sl_data[i]<<"   "<<sl_data_err[i]<<"   "<<y_dl_mc[i]<<"   "<<y_dl_data[i]<<"   "<<dl_data_err[i]<<"\n";
	}
	fout.close();

	fout_plot.open("plot.gnu");	
	
	fout_plot<<"set terminal png font 'Helvetica,20' lw 2 \n";
	fout_plot<<"se ou \"SL.png\"\n";
	fout_plot<<"set title \"Single Lepton\"\n";
	fout_plot<<"set xrange ["<<dlim<<":"<<ulim<<"] \n";
	fout_plot<<"set yrange [0:0.22] \n";
	fout_plot<<"set xlabel \""<<xlabel<<"\"\n";
	fout_plot<<"set ylabel \"Normalized Counts\"\n";
	fout_plot<<"set boxwidth "<<bin_size<<" \n";
	//fout_plot<<"set style fill solid border -1 \n";
	fout_plot<<"set style fill solid \n";
	fout_plot<<"p \"dist_hist.txt\" u 1:2 w boxes lt 1 ti \"MC\", \"dist_hist.txt\" u 1:3:($3+$4):($3-$4) with errorbars pt 7 ps 1 lc -1 ti \"DATA\" \n\n";
	
	fout_plot<<"se ou \"DL.png\"\n";
	fout_plot<<"set title \"Di Lepton\"\n";
	fout_plot<<"set xrange ["<<dlim<<":"<<ulim<<"] \n";
	fout_plot<<"set yrange [0:0.22] \n";
	fout_plot<<"set xlabel \""<<xlabel<<"\"\n";
	fout_plot<<"set ylabel \"Normalized Counts\"\n";
	fout_plot<<"set boxwidth "<<bin_size<<" \n";
	//fout_plot<<"set style fill solid border -1 \n";
	fout_plot<<"set style fill solid \n";
	fout_plot<<"p \"dist_hist.txt\" u 1:5 w boxes lt 1 ti \"MC\", \"dist_hist.txt\" u 1:6:($6+$7):($6-$7) with errorbars pt 7 ps 1 lc -1 ti \"DATA\" \n\n";
	
	fout_plot.close();
	system("gnuplot plot.gnu");
	
	system("rm -rf plot.gnu");
	system("rm -rf temp*");
	system("rm -rf dist_hist.txt");
	return 0;
}
