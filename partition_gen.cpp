/*
 * partition_gen.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: Huanhuan Yang	huan2yang@outlook.com
 */

#include <stdlib.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <tr1/unordered_set>
#include <vector>
#include <math.h>
#include <assert.h>

#include "triangulation.h"
#include "setup_routines.h"
#include "loop_routines.h"
#include "finalize_routines.h"
#include "initial_points.h"
#include "densities.h"

// Uses namespaces std and tr1. tr1 is used for unordered_set which gives unique triangulation at the end.
using namespace std;
using namespace tr1;

//This routine is using Llyod method to generate small-size SCVT (on unit sphere) that will be used later for partition 
int main(int argc, char **argv){
	//Processor information and global communicator
	//master id is always 0
	mpi::environment env(argc, argv);
	mpi::communicator world;
	int id, num_procs;
	id = world.rank();
	num_procs = world.size();

	if(id==0 && argc<3)	{
		printf("Command argument error! Run with 3 args as: ./partition_gen.exe -p 320\n");
		exit(0);
	}
	int num_pts = atoi(argv[2]);

	//Each processor has a copy of all points, and has a n_points (new points) vector for which points it should update.
	vector<pnt> points;
	vector<pnt> n_points;
	vector<pnt>::iterator point_itr;
	//Each processor has a list of all regions, as well as it's own regions (only one per processor currently)
	vector<region> regions;
	vector<region> my_regions;
	vector<region>::iterator region_itr;
	vector<tri> all_triangles;

	int it, i;
	int stop;
	mpi::request *ave_comms, *max_comms, *l1_comms;
	double *my_ave, *my_max, *my_l1;
	double glob_ave, glob_max, glob_l1;
	double last_glob_ave, last_glob_max, last_glob_l1;
	optional ave_opti, max_opti, l1_opti;

	// Input flags for the Triangle package
	string flags_str = "QBPIOYYiz";
	char * flags = new char[flags_str.size()+1];
	strcpy(flags,flags_str.c_str());

	my_ave = new double[num_procs];
	my_max = new double[num_procs];
	my_l1 = new double[num_procs];
	ave_comms = new mpi::request[num_procs];
	max_comms = new mpi::request[num_procs];
	l1_comms = new mpi::request[num_procs];
	last_glob_ave = 1.0;
	last_glob_max = 1.0;
	last_glob_l1 = 1.0;

	// Was global constants: value given in config.xml. Here we fixed the values
	int conv = 2; //convergence_check by Maximum movement
	int max_it = 5001;
	double eps = 1.0E-6;
	int div_levs = 1; // level for triangle_division in integration
	int quad_rule = 3; //7 Point Rule;  			
	int use_barycenter = 0; // use circumcenter as delaunay_triangle_center
	int sort_method = sort_dot;

	// Read in regions. Setup initial point set
	if(id == 0){
		printf("Initial points being created with Fibonacci Grid.\n");
		makeFibonacciGridPoints(num_pts, points);

		string regFile = "./region/RegionList."+to_string(max(2,num_procs));
		string connFile = "./region/RegionTriangulation."+to_string(max(2,num_procs));
		int read = buildRegions(id, regions, regFile, connFile);
		if(read==1){
			printf("Error: RegionList.n and/or RegionTriangulation.n DON'T EXIST!\n");
			printf("Region size must equal num_procs, if num_procs > 1. You can run partition_gen.exe to generate new partition files\n");
			exit(1);
		}
	}
	mpi::broadcast(world,regions,0);
	mpi::broadcast(world,points,0);

	// Each processor clears my_regions (to make sure it's empty) and add it's own region into it's list.
	// If there is only 1 processor, that processor takes all regions, which is a serial computation.
	my_regions.clear();
	if(num_procs > 1){
		my_regions.push_back(regions.at(id));
		if(num_procs != regions.size()){
			cout << "Region error ---- Region size must equal number of processors." << endl;
			cout << " Or you must only use 1 processor" << endl;
			assert(num_procs == regions.size());
		}
	} else {
		for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
			my_regions.push_back((*region_itr));
		}
	}

	// Each processor needs to setup the quadrature rules.
	Quadrature quadr(quad_rule);

	// Lloyd iterations
	stop = 0;
	if(id==0)
		cout << "\nit ||  glob_ave  glob_l1  glob_max  \n"<< endl;
	for(it = 0; it < max_it && !stop; it++){
		glob_ave = 0.0;
		glob_max = 0.0;
		glob_l1 = 0.0;
		for(i = 0; i < num_procs; i++){
			my_ave[i] = 0.0;
			my_max[i] = 0.0;
			my_l1[i] = 0.0;
		}

		clearRegions(id, my_regions);
		sortPoints(id, regions, points, sort_method, my_regions);
		triangulateRegions(id, flags, my_regions);
		integrateRegions(id, div_levs, quadr, use_barycenter, regions, my_regions, points, n_points);
		computeMetrics(id, points, n_points, my_ave[id],my_max[id], my_l1[id]);

		// Start non-blocking sends and receives of metrics
		if(id == 0){
			for(i = 1; i < num_procs; i++){
				ave_comms[i] = world.irecv(i,msg_ave,my_ave[i]);
				max_comms[i] = world.irecv(i,msg_max,my_max[i]);
				l1_comms[i] = world.irecv(i,msg_l1,my_l1[i]);
			}
		} else {
			ave_comms[id] = world.isend(0,msg_ave,my_ave[id]);
			max_comms[id] = world.isend(0,msg_max,my_max[id]);
			l1_comms[id] = world.isend(0,msg_l1,my_l1[id]);
		}

		transferUpdatedPoints(world, my_regions, n_points, points);

		// Finish metrics sends and receives, check for convergence and broadcast a stop request to all processors.
		if(id == 0)
		{
			glob_ave = my_ave[id];
			glob_max = my_max[id];
			glob_l1 = my_l1[id];

			for(i = 1; i < num_procs; i++){
				ave_opti = ave_comms[i].test();
				max_opti = max_comms[i].test();
				l1_opti = l1_comms[i].test();

				if(!ave_opti) ave_comms[i].wait();
				if(!max_opti) max_comms[i].wait();
				if(!l1_opti) l1_comms[i].wait();

				glob_ave += my_ave[i];
				glob_max = std::max(glob_max, my_max[i]);
				glob_l1 += my_l1[i];
			}
			glob_ave = sqrt(glob_ave);
			glob_l1 = glob_l1/sqrt(points.size()); 
			glob_max = glob_max*sqrt(points.size());

			if(id==0 && it%100==0)
				cout << it << " " << glob_ave << " " << glob_l1 << " " << glob_max << endl;

			if(conv == 1 && glob_ave < eps){
				cout << "Converged on average movement." << endl;
				stop = 1;
			} else if(conv == 2 && glob_max < eps){
				cout << "Converged on maximum movement." << endl;
				stop = 1;
			}

			last_glob_ave = glob_ave;
			last_glob_max = glob_max;
			last_glob_l1 = glob_l1;
		}

		mpi::broadcast(world,stop,0);
	}

	//Gather all updated points onto master processor, for printing to RegionList.n file
	gatherAllUpdatedPoints(world, n_points, points);

	// Compute final triangulation by merging all triangulations from each processor into an
	// unordered_set, and then ordering them ccw before printing them out.
	// write triangles to RegionTriangulation.n file
	clearRegions(id, my_regions);
	sortPoints(id, regions, points, sort_vor, my_regions);
	makeFinalTriangulations(id, flags, regions, my_regions);

	printMyFinalTriangulation(world, my_regions, all_triangles, "./region/RegionTriangulation."+to_string(num_pts));

	if(id == 0){
		ofstream end_pts("./region/RegionList."+to_string(num_pts));
		for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
			end_pts << (*point_itr) << endl;
		}
		end_pts.close();
	}

	return 0;
}





