/*
 * scvt_HLBFGS.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: Huanhuan Yang	huan2yang@outlook.com
 *
 *  The Overlapping domain decomposition method for parallelization is following
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff
 *
 *	This program is using non-overlapping local vector for storing points, 
 *  as compared with another version which use disjointly distributed vector for updating points during optimization
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
#include "HLBFGS/HLBFGS.h"
#include "HLBFGS/HLBFGS_BLAS.h"
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
vector<double> my_lloyds;
vector<double> my_bots;
vector<int> disjDistrIdx;
Quadrature quadr;
char * flags;
// Global constant parameters: value given in data file: parameters
int sort_method;
int div_levs;
int use_barycenter;




//////////////////////////////////////////////////////////////////////////
void evalfunc(int my_N, double* my_x, double *my_prev_x, double* f, double* my_g)
{
	*f = 0;
	double my_energy;
	vector<pnt> my_points;
	vector<pnt>	disj_grad;
	vector<pnt>	disj_lloyd;
	vector<pnt> my_gradPts;
	vector<double> disj_bots;

	pnt p;
	for(int i=0; i<my_N/2; ++i)
	{
		my_x[2*i] = fmod(my_x[2*i],M_PI);
		if(my_x[2*i]>M_PI/2.0 )	my_x[2*i] -= M_PI;
		if(my_x[2*i]<-M_PI/2.0 )	my_x[2*i] += M_PI;
		p = pntFromLatLon(my_x[2*i],my_x[2*i+1]);
		p.idx = disjDistrIdx[i];
		my_points.push_back(p);
	}
	transferUpdatedPoints(world, my_regions, my_points, points);
	clearRegions(id, my_regions);
	sortPoints(id, regions, points, sort_method, my_regions);
	triangulateRegions(id, flags, my_regions);
	inteEnergGrad(id, div_levs, quadr, use_barycenter, regions, my_regions, points, 
					my_energy, disj_grad, disj_lloyd, disj_bots);
	mpi::reduce(world, my_energy, *f, std::plus<double>(), 0);
	mpi::broadcast(world, *f, 0);
	if(id==0)	cout << "\n f ="<< *f <<endl;

	transferByDisjDistrIdx(world, my_regions, points.size(), disjDistrIdx, 
							disj_grad, my_gradPts, disj_lloyd, my_lloyds, disj_bots, my_bots);

	double lat, lon;
	for(int i=0; i<my_N/2; ++i){
		p = points.at(disjDistrIdx[i]);
		//lat = -1.0*my_gradPts[i].x*p.z*p.x/sqrt(1-p.z*p.z+1e-100)
		//		-my_gradPts[i].y*p.z*p.y/sqrt(1-p.z*p.z+1e-100)
		//		+my_gradPts[i].z*sqrt(1-p.z*p.z);
		lat = -1.0*my_gradPts[i].x*sin(my_x[2*i])*cos(my_x[2*i+1])
				-my_gradPts[i].y*sin(my_x[2*i])*sin(my_x[2*i+1])
				+my_gradPts[i].z*cos(my_x[2*i]);
		lon = -1.0*my_gradPts[i].x*p.y + my_gradPts[i].y*p.x;

		my_g[2*i] = lat;
		my_g[2*i+1] = lon;
	}

}


void checkGrad(int my_N, double* my_x, double *my_prev_x, double* f, double* my_g) 
{
	double h=1e-8, f1, f2;
	double tmp[my_N];
	
	for(int i=0; i<my_N/2; ++i)
	{
		evalfunc(my_N, my_x, 0, &f1, tmp);
		my_x[2*i] = my_x[2*i]+h;		
		evalfunc(my_N, my_x, 0, &f2, tmp);
		my_x[2*i] = my_x[2*i]-h;		
		my_g[2*i] = (f2 - f1)/h;
	
		my_x[2*i+1] = my_x[2*i+1]+h;
		evalfunc(my_N, my_x, 0, &f2, tmp);
		my_x[2*i+1] = my_x[2*i+1]-h;
		my_g[2*i+1] = (f2 - f1)/h;
	}

}

//////////////////////////////////////////////////////////////////////////
void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
	if(id==0){
		std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
	}
}


void HLBFGS_DSCALDV(const int n, const double *a, double *y)
{
	int i = 0;
#ifdef USE_OPENMP
#pragma omp parallel for private(i) reduction(+:result)
#endif
	for (i = 0; i < n; i++)
	{
		y[i] = y[i] / a[i];
	}
}


//////////////////////////////////////////////////////////////////////////
void defined_update_Hessian(int N, int M, double *q, double *s, double *y,
						   int cur_pos, double *diag, int INFO[], double *x, mpi::communicator* comm=0)
{
	if (M <= 0)
	{
		return;
	}

	if (INFO[2] == 0)
	{
		for(int j=0; j<N; ++j)
		{
			q[j] = my_lloyds[j] - x[j];
		}
	}
	else if (INFO[2] > 0)
	{
		if (INFO[3] == 0 || INFO[3] == 1)
			HLBFGS_UPDATE_Hessian(N, M, q, s, y, cur_pos, diag, INFO, x, &world);
		if (INFO[3] == 2)
		{
			HLBFGS_DSCALDV(N, &my_bots[0], q);
			HLBFGS_DSCAL(N, 0.5, q);
		}
	}
}






//This routine is using quasi-Newton optimization for SCVT computation. In particular, LBFGS preconditioned by Lloyd step is used.
int main(int argc, char **argv){
	//Processor information and global communicator
	//master id is always 0
	mpi::environment env(argc, argv);
	id = world.rank();
	num_procs = world.size();

	vector<pnt> n_points;
	vector<pnt>::iterator point_itr;
	vector<region>::iterator region_itr;
	vector<tri> all_triangles;
	pnt p;
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
	//double max_bdryResol = dataFile("scvt/resolution/max_boundary_resolution", 40000.0);

	// Input flags for the Triangle package
	string flags_str = "QBPIOYYiz";
	flags = new char[flags_str.size()+1];
	strcpy(flags,flags_str.c_str());

	// constant parameters for HLBFGS, reading from dataFile
	double Lloyd_tol = dataFile("HLBFGS/real_tol/Lloyd_tol", 1.0);
	int mem_num = dataFile("HLBFGS/memory_number", 7);
	double parameter[10];
	int info[15];
	//initialize
	INIT_HLBFGS(parameter, info);
	parameter[0] = dataFile("HLBFGS/real_tol/c_1", 1e-4);
	parameter[2] = dataFile("HLBFGS/real_tol/c_2", 0.9); 
	parameter[7] = dataFile("HLBFGS/real_tol/g_tol", 1e-10); 
	parameter[8] = dataFile("HLBFGS/real_tol/f-f_tol", 1e-16); 
	parameter[9] = dataFile("HLBFGS/real_tol/dx_tol", 1e-16); 
	info[0] = dataFile("HLBFGS/int_tol/max_LS", 20);
	info[3] = dataFile("HLBFGS/int_tol/invH0_strategy", 0);
	info[4] = dataFile("HLBFGS/int_tol/max_itr", 100000);								
	info[14] = dataFile("HLBFGS/int_tol/max_reset", 0);
	if(info[3]==2)	info[12] = 1;

	//Define timers for performance studies
	const int num_global_timers = 4;
	mpi_timer global_timers[num_global_timers];
	string global_names[num_global_timers] = {"Global Time", "Final Gather", "Final Triangulation", "Final Bisection"};
	for(i = 0; i < num_global_timers; i++){
		global_timers[i] = mpi_timer(global_names[i]);
	}
	global_timers[0].start(); // Global Time Timer

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
		/////////////////////////////////////////////////////
			// Lloyd iteration as starting up
		it = 0;
		bool stop = false;
		do{
			clearRegions(id, my_regions);
			sortPoints(id, regions, points, sort_method, my_regions);
			triangulateRegions(id, flags, my_regions);
			integrateRegions(id, div_levs, quadr, use_barycenter, regions, my_regions, points, n_points);

/*			if(it > max_it_no_proj){
				proj_alpha = max((double)(it-max_it_no_proj), 0.0)/max((double)max_it_scale_alpha, 1.0);
				projectToBoundary(proj_alpha, points, boundary_points, n_points, my_regions);
			}
*/
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
			if(id==0) cout << "Initial Lloyd itr "<< ++it <<": dx_l2 = "<< glob_l2 << endl;

			if(glob_l2 > Lloyd_tol){ 
				transferUpdatedPoints(world, my_regions, n_points, points);
			}else{
				gatherAllUpdatedPoints(world, n_points, points);
				stop = true;
			}
		}while(!stop);

		/////////////////////////////////////////////////////
		// Quasi-Newton iteration when dx <= Lloyd_tol

		clearRegions(id, my_regions);
		sortPoints(id, regions, points, sort_method, my_regions);
		getDisjointIndex(regions, my_regions, disjDistrIdx);

		std::vector<double> x(disjDistrIdx.size()*2);
		for(i=0; i<disjDistrIdx.size(); ++i)
		{
			x[2*i] = points.at(disjDistrIdx[i]).getLat();
			x[2*i+1] = points.at(disjDistrIdx[i]).getLon();
		}

/*		//check gradient
		if(num_procs==1)
		{
			double f;
			double g1[x.size()];
			double g2[x.size()];
			evalfunc(x.size(), &x[0], 0, &f, g1);
			checkGrad(x.size(), &x[0], 0, &f, g2);
			if(id==0)	cout<<"\n grad =\n ";
			for(i=0; i<x.size(); ++i)
				if(id==0)	cout<<g1[i]<<" ";

			if(id==0)	cout<<"\n\n check FD grad =\n ";
			for(i=0; i<x.size(); ++i)
				if(id==0)	cout<<g2[i]<<" ";
			if(id==0)	cout<<"\n\n check error grad =\n ";
			for(i=0; i<x.size(); ++i)
				if(id==0)	cout<<g1[i]-g2[i]<<" ";
		}else
			if(id==0)	cout<<"\n Warning: check gradient only when using 1 processor ! ";
*/
		if(id==0)
			cout << "\n itr(Quasi-Newton): num_feval |  f_val  |  g_norm  |\n"<< endl;

		HLBFGS(disjDistrIdx.size()*2, mem_num, &x[0], evalfunc, 0, defined_update_Hessian, newiteration, parameter, info, &world);

		if(id==0)
			cout << "-----------Summary: totoal number of f (and g) evals: "<< it+info[1]<< endl;
		vector<pnt> my_points;
		pnt p;
		for(i=0; i<x.size()/2; ++i)
		{
			p = pntFromLatLon(x[2*i],x[2*i+1]);
			p.idx = disjDistrIdx[i];
			my_points.push_back(p);
		}
		gatherAllUpdatedPoints(world, my_points, points);

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
	global_timers[0].stop();

	//Compute average points per region for diagnostics
	num_avePoints = 0;
	num_myPoints = 0;
	clearRegions(id, my_regions);
	sortPoints(id, regions, points, sort_method, my_regions);
	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		num_myPoints += (*region_itr).points.size();
	}

	mpi::reduce(world, num_myPoints, num_avePoints, std::plus<int>(), 0);
	num_avePoints = num_avePoints / regions.size();
	if(id == 0){
		cout << "\nAverage points per region: " << num_avePoints << endl;
	}


	// Compute final triangulation by merging all triangulations from each processor into an
	// unordered_set, and then ordering them ccw before printing them out.
	// write triangles to triangles.dat
	global_timers[2].start(); // Final Triangulation Timer
	clearRegions(id, my_regions);
	sortPoints(id, regions, points, sort_vor, my_regions);
	makeFinalTriangulations(id, flags, regions, my_regions);
	global_timers[2].stop();

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
	global_timers[3].start(); // Final Bisection Timer
	if(num_bisections>0) 
		bisectTriangulation(flags, world, my_regions, all_triangles, regions, points, 1, "bisected_points.dat");
	global_timers[3].stop();

	//Print out final timers, for global times.
	if(id == 0){
		cout << endl << " ---- Final Timers ---- " << endl;
		for(i = 0; i < num_global_timers; i++){
			cout << global_timers[i];
		}

	}

	return 0;
}










