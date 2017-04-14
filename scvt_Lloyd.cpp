/*
 * scvt_Lloyd.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: Huanhuan Yang	huan2yang@outlook.com
 *
 *  This programs reshapes the data structures from 
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff
 *
 *  Modification is made on the convergence control
 *  Add option for energy function evaluation
 */

// Enable _DEBUG for output of routines. Useful for debugging any issues in mpi.
//#define _DEBUG

#include <stdlib.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <tr1/unordered_set>
#include <vector>
#include <math.h>
#include <assert.h>

#include "GetPot.hpp"
#include "triangulation.h"
#include "setup_routines.h"
#include "loop_routines.h"
#include "finalize_routines.h"
#include "initial_points.h"
#include "densities.h"

// Uses namespaces std and tr1. tr1 is used for unordered_set which gives unique triangulation at the end.
using namespace std;
using namespace tr1;

// Global variables
int id, num_procs;
mpi::communicator world;
//Each processor has a list of all regions, as well as it's own regions (only one per processor currently)
vector<region> regions;
vector<region> my_regions;
vector<pnt> points;
vector<pnt> n_points;
Quadrature quadr;
char * flags;
// Global constant parameters: value given in data file: parameters
int sort_method;
int div_levs;
int use_barycenter;

//This routine is using quasi-Newton optimization for SCVT computation. In particular, LBFGS preconditioned by Lloyd step is used.
int main(int argc, char **argv){
	//Processor information and global communicator
	//master id is always 0
	mpi::environment env(argc, argv);
	id = world.rank();
	num_procs = world.size();

	vector<pnt>::iterator point_itr;
	pnt p;
	vector<region>::iterator region_itr;
	vector<tri> all_triangles;
	//vector<pnt> boundary_points;
	//vector<pnt>::iterator boundary_itr;
	int it, i, it_bisect, num_avePoints, num_myPoints;
	double my_l2, my_max, my_l1, glob_l2, glob_max, glob_l1;	//for Lloyd step computation
	//double proj_alpha;

    string dataFileName;
    GetPot command_line (argc, argv);
    dataFileName = command_line.follow ("parameters", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

	// constant parameters reading from dataFile
	int points_begin = dataFile("scvt/initial/initial_point_set", 3); 
	int num_pts = dataFile("scvt/initial/number_of_generated_points", 12);
	sort_method = dataFile("partition/sort_method", 0);
	int num_bisections = dataFile("scvt/num_bisections", 0);
	int quad_rule = dataFile("scvt/integration/quadrature_rule", 0);	
	use_barycenter = dataFile("scvt/integration/delaunay_triangle_center", 0);
	div_levs = dataFile("scvt/integration/division_levs", 1);
	bool save_before_bisect = dataFile("fileIO/save_before_bisect", false);
	double dx_tol = dataFile("Lloyd/dx_tol", 1.0);
	int max_itr = dataFile("Lloyd/max_itr", 10000);
	bool eval_f = dataFile("Lloyd/eval_f", false);
	string fileName = dataFile("Lloyd/itr_polarFile", "");
	ofstream itrPolarFile;
	//double max_bdryResol = dataFile("scvt/resolution/max_boundary_resolution", 40000.0);

	// Input flags for the Triangle package
	string flags_str = "QBPIOYYiz";
	flags = new char[flags_str.size()+1];
	strcpy(flags,flags_str.c_str());

	//Define timers for performance studies
	const int num_timers = 4;
	mpi_timer timers[num_timers];
	string t_names[num_timers] = {"Cumulative Time", "Global Time", "Final Triangulation", "Final Bisection"};
	for(i = 0; i < num_timers; i++){
		timers[i] = mpi_timer(t_names[i]);
	}
	timers[1].start(); // Global Time Timer

	// Setup initial point set and build regions
	if(id == 0){
		if (points_begin == 0){
			readPoints(num_pts, points);
			cout << "\n" << num_pts <<" points being read in from SaveVertices." << endl;
		} else if(points_begin == 1){
			makeMCPoints(num_pts, points);
			cout << "\n" << num_pts <<" points being created with Monte Carlo." << endl;
		} else if(points_begin == 2){
			makeGeneralizedSpiralPoints(num_pts, points);
			cout << "\n" << num_pts << " points being created with Generalized Spiral." << endl;
		} else if(points_begin == 3){
			makeFibonacciGridPoints(num_pts, points);
			cout << "\n" << num_pts << " points being created with Fibonacci Grid." << endl;
		}

		//readBoundaries(max_bdryResol, boundary_points);

		string regFile = "./region/RegionList."+to_string(max(2,num_procs));
		string connFile = "./region/RegionTriangulation."+to_string(max(2,num_procs));
		int read = buildRegions(id, regions, regFile, connFile);
		if(read==1){
			printf("Error: RegionList.n and/or RegionTriangulation.n DON'T EXIST!\n");
			printf("Region size must equal num_procs, if num_procs > 1. You can run partition_gen.exe to generate new partition files\n");
			exit(1);
		}
		
		ofstream pts_out("point_initial.dat");
		for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
			pts_out << (*point_itr) << endl;
		}
		pts_out.close();
	}

	// Broadcast regions, and initial point set to each processor.
	mpi::broadcast(world,regions,0);
	mpi::broadcast(world,num_pts,0);
	mpi::broadcast(world,points,0);
	//mpi::broadcast(world,boundary_points,0);

	// Each processor needs to setup the quadrature rules.
	quadr.setQuadRule(quad_rule);

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
		vector<region>::iterator region_itr;
		for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
			my_regions.push_back((*region_itr));
		}
	}

	// Loop over it_bisect
	for(it_bisect = 0; it_bisect <= num_bisections; it_bisect++)
	{
		if(it_bisect==num_bisections && id==0)
		{
			itrPolarFile.open(fileName.c_str());
			//itrPolarFile.open(fileName.c_str(),ios::app);
		}

		if(id==0)
			cout << "\nLloyd itr ||  dx_l2  |  f  |  df_normPolar\n"<< endl;

		bool stop = false;
		
		for(it = 0; it < max_itr && !stop; it++)
		{
			//timers[0].init();

			if(it_bisect==num_bisections)
				timers[0].start();
			clearRegions(id, my_regions);
			sortPoints(id, regions, points, sort_method, my_regions);
			triangulateRegions(id, flags, my_regions);
			integrateRegions(id, div_levs, quadr, use_barycenter, regions, my_regions, points, n_points);
			
			if(it_bisect==num_bisections)
				timers[0].stop();

/*			if(it > max_it_no_proj){
				proj_alpha = max((double)(it-max_it_no_proj), 0.0)/max((double)max_it_scale_alpha, 1.0);
				projectToBoundary(proj_alpha, points, boundary_points, n_points, my_regions);
			}
*/

			double fval=0.0, grad_norm=0.0, grad_normPolar=0.0;
			if(eval_f)
			{
				double my_energy;
				inteEnergy(id, div_levs, quadr, use_barycenter, regions, my_regions, points, my_energy);
				mpi::reduce(world, my_energy, fval, std::plus<double>(), 0);

				vector<pnt>	distr_grad, gradients;
				inteGradient(id, div_levs, quadr, use_barycenter, regions, my_regions, points, distr_grad);
				gradients.resize(points.size());
				gatherAllUpdatedPoints(world, distr_grad, gradients);
				for(int k=0; k<gradients.size(); ++k)
					grad_norm += gradients[k].magnitude2();
				grad_norm = sqrt(grad_norm);

				double p_lat, p_lon, g_lat, g_lon;
				for(int k=0; k<points.size(); ++k){
					p_lat = points[k].getLat();
					p_lon = points[k].getLon();
					g_lat = -1.0*gradients[k].x*sin(p_lat)*cos(p_lon) - gradients[k].y*sin(p_lat)*sin(p_lon)
							+ gradients[k].z*cos(p_lat);
					g_lon = -1.0*gradients[k].x*points[k].y + gradients[k].y*points[k].x;

					grad_normPolar += g_lat*g_lat + g_lon*g_lon;
				}
				grad_normPolar = sqrt(grad_normPolar);		
			}

			if(it_bisect==num_bisections)
				timers[0].start();
			computeMetrics(id, points, n_points, my_l2, my_max, my_l1);
			mpi::reduce(world, my_l2, glob_l2, std::plus<double>(), 0);
			//mpi::reduce(world, my_max, glob_max, mpi::maximum<double>(), 0);
			//mpi::reduce(world, my_l1, glob_l1, std::plus<double>(), 0);
			glob_l2 = sqrt(glob_l2);
			//glob_max = glob_max*sqrt(points.size());
			//glob_l1 = glob_l1/sqrt(points.size()); 
			mpi::broadcast(world, glob_l2, 0);
			//mpi::broadcast(world, glob_max, 0);
			//mpi::broadcast(world, glob_l1, 0);

			if(glob_l2 > dx_tol){ 
				// the two functions have no difference in serial
				if(eval_f)	gatherAllUpdatedPoints(world, n_points, points);
				else transferUpdatedPoints(world, my_regions, n_points, points);
			}else{
				gatherAllUpdatedPoints(world, n_points, points);
				stop = true;
			}
			if(it_bisect==num_bisections)
				timers[0].stop();

			if(id==0) 
			{	
				cout << it << " " << glob_l2 << " " << fval << " " << grad_normPolar << endl;

				if(it_bisect==num_bisections)
					itrPolarFile << it << " " << glob_l2 << " " << fval << " " << grad_normPolar << " " << timers[0].total_time <<"\n" ;
				if(stop)
					itrPolarFile.close();
			}
		}

		if(save_before_bisect && it_bisect < num_bisections)
		{
			string postNm = to_string(points.size());

			clearRegions(id, my_regions);
			sortPoints(id, regions, points, sort_vor, my_regions);
			makeFinalTriangulations(id, flags, regions, my_regions);
			printMyFinalTriangulation(world, my_regions, all_triangles, "triangles.dat."+postNm);

			if(id == 0){
				ofstream end_pts("end_points.dat."+postNm);
				ofstream pt_dens("point_density.dat."+postNm);
				//ofstream bdry_pts("boundary_points.dat");
				for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
					end_pts << (*point_itr) << endl;
					pt_dens << density((*point_itr)) << endl;
				}
				//for(boundary_itr = boundary_points.begin(); boundary_itr != boundary_points.end(); boundary_itr++){
				//	bdry_pts << (*boundary_itr) << endl;
				//}
				end_pts.close();
				pt_dens.close();
				//bdry_pts.close();
			}
		}

		// Bisect if needed
		if(it_bisect < num_bisections){
			bisectTriangulation(flags, world, my_regions, all_triangles, regions, points, 0);
		} else {
			if(id == 0 && num_bisections>0){
				cout << "\nNo more bisections for convergence" << endl;
			}
		}

	}
	timers[1].stop();

	// Compute final triangulation by merging all triangulations from each processor into an
	// unordered_set, and then ordering them ccw before printing them out.
	// write triangles to triangles.dat
	timers[2].start(); // Final Triangulation Timer
	clearRegions(id, my_regions);
	sortPoints(id, regions, points, sort_vor, my_regions);
	makeFinalTriangulations(id, flags, regions, my_regions);
	timers[2].stop();

	//Compute average points per region for diagnostics
	num_avePoints = 0;
	num_myPoints = 0;
	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		num_myPoints += (*region_itr).points.size();
	}
	mpi::reduce(world, num_myPoints, num_avePoints, std::plus<int>(), 0);
	num_avePoints = num_avePoints / regions.size();
	if(id == 0){
		cout << "\nAverage points per region: " << num_avePoints << endl;
	}


	printMyFinalTriangulation(world, my_regions, all_triangles);

	if(id == 0){
		ofstream end_pts("end_points.dat");
		ofstream pt_dens("point_density.dat");
		//ofstream bdry_pts("boundary_points.dat");
		for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
			end_pts << (*point_itr) << endl;
			pt_dens << density((*point_itr)) << endl;
		}
		//for(boundary_itr = boundary_points.begin(); boundary_itr != boundary_points.end(); boundary_itr++){
		//	bdry_pts << (*boundary_itr) << endl;
		//}
		end_pts.close();
		pt_dens.close();
		//bdry_pts.close();
	}

	//Bisect all edges of all triangles to give an extra point set at the end, SaveVertices
	timers[3].start(); // Final Bisection Timer
	if(num_bisections>0) 
		bisectTriangulation(flags, world, my_regions, all_triangles, regions, points, 1, "bisected_points.dat");
	timers[3].stop();

	//Print out final timers, for global times.
	if(id == 0){
		cout << endl << " ---- Final Timers ---- " << endl;
		for(i = 1; i < num_timers; i++){
			cout << timers[i];
		}

	}

	return 0;
}










