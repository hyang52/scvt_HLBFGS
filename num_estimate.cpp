/*
 * num_estimate.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: huanhuan
 */
#include <stdlib.h>
#include <math.h>
#include <boost/mpi.hpp>
#include "GetPot.hpp"
#include "initial_points.h"
#include "densities.h"

using namespace std;
namespace mpi = boost::mpi;

//This routine draws MC pts to estimate required num of CVT cells for reaching certain resolution
int main(int argc, char **argv){
	//Processor information and global communicator
	//master id is always 0
	mpi::environment env(argc, argv);
	mpi::communicator world;
	int id = world.rank();
	int num_procs = world.size();

    string dataFileName;
    GetPot command_line (argc, argv);
    dataFileName = command_line.follow ("parameters", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

	//Read in parameters from dataFile.
	double max_resolution = dataFile("scvt/resolution/max_inner_resolution", 120000.0);
	double min_resolution = dataFile("scvt/resolution/min_inner_resolution", 30000.0);

	vector<pnt> points;
	int num_pts=1e6;
	double area_per_sample = 4.0 * M_PI * pow(6371000,2.0) / num_pts;
	double sum_nhexs, nhexs, dl, hex_area;

	if(id==0)
		cout << "Estimating the required number of cells to reach max_inner_resolution = " << max_resolution << endl;

	nhexs = 0.0;
	cout << "On proc "<<id<<": " ;
	makeMCPoints(num_pts, points);

	// find min rho, which corresponds to max h. then use the cell as reference
	double dens[num_pts];
 	for(int i=0; i<num_pts;++i){
		dens[i] = density(points[i]);
	}
	double rho_min = *std::min_element(dens,dens+num_pts);
	for(int i=0; i<num_pts;++i)
	{
		dl = max_resolution * pow(rho_min/dens[i],0.25);
		hex_area = sqrt(3.0) / 2.0 * pow(dl,2.0);
		nhexs += area_per_sample / hex_area;
	}
	mpi::reduce(world, nhexs, sum_nhexs, std::plus<double>(), 0);
	if(id==0)
	{
		cout << "-------------------------------------------------------------------------------------\n";
		cout << "Estimated # hexs to reach max_inner_resolution " << max_resolution << " is " << int(sum_nhexs/num_procs) << endl;
	}

	double rho_max = *std::max_element(dens,dens+num_pts);
	nhexs = 0.0;
	for(int i=0; i<num_pts;++i)
	{
		dl = min_resolution * pow(rho_max/dens[i],0.25);
		hex_area = sqrt(3.0) / 2.0 * pow(dl,2.0);
		nhexs += area_per_sample / hex_area;
	}
	mpi::reduce(world, nhexs, sum_nhexs, std::plus<double>(), 0);
	if(id==0)
	{
		cout << "Estimated # hexs to reach min_inner_resolution " << min_resolution << " is " << int(sum_nhexs/num_procs) << endl;
		//cout << "rho_max = "<<rho_max<<"; rho_min="<<rho_min<<endl;
		cout << "-------------------------------------------------------------------------------------\n";
		cout << "Given refinement step k, the necessary num of initial pts is (N^final - 2)/4^k + 2 " << endl;
		cout << "Given N^0, the necessary num of refinement steps is log_4( (N^final - 2)/(N^0-2) ) \n";
		cout << "-------------------------------------------------------------------------------------" << endl;
	}
	return 0;
}


