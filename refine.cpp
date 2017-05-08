/*
 * refine.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: Huanhuan Yang	huan2yang@outlook.com
 *
 */

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
namespace mpi = boost::mpi;

// Global variables
int id, num_procs;
mpi::communicator world;
//Each processor has a list of all regions, as well as it's own regions (only one per processor currently)
vector<region> regions;
vector<region> my_regions;
vector<pnt> points;
char * flags;
int it_bisect;
int sort_method = sort_vor;
int num_bisections;
int num_pts;

// this program read points from SaveVertices file, do continuous refinement, 
// and save the refined mesh to "end_points.dat."+postNm
int main(int argc, char **argv){
	//Processor information and global communicator
	//master id is always 0
	mpi::environment env(argc, argv);
	id = world.rank();
	num_procs = world.size();

	if(id==0 && argc<3)	{
		printf("Command argument error! Run with 3 args as: ./refine.exe -r 2\n");
		exit(0);
	}
	num_bisections = atoi(argv[2]);

	vector<pnt> n_points;
	vector<pnt>::iterator point_itr;
	vector<region>::iterator region_itr;
	vector<tri> all_triangles;
	pnt p;

	// Input flags for the Triangle package
	string flags_str = "QBPIOYYiz";
	flags = new char[flags_str.size()+1];
	strcpy(flags,flags_str.c_str());

	// Setup initial point set and build regions
	if(id == 0){
		readPoints(num_pts, points);
		cout << "\n" << num_pts <<" points being read in from SaveVertices." << endl;

		string regFile = "./region/RegionList."+to_string(max(2,num_procs));
		string connFile = "./region/RegionTriangulation."+to_string(max(2,num_procs));
		int read = buildRegions(id, regions, regFile, connFile);
		if(read==1){
			printf("Error: RegionList.n and/or RegionTriangulation.n DON'T EXIST!\n");
			printf("Region size must equal num_procs, if num_procs > 1. You can run partition_gen.exe to generate new partition files\n");
			exit(1);
		}
	}
	// Broadcast regions, and initial point set to each processor.
	mpi::broadcast(world,regions,0);
	mpi::broadcast(world,num_pts,0);
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
		vector<region>::iterator region_itr;
		for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
			my_regions.push_back((*region_itr));
		}
	}

	clearRegions(id, my_regions);
	sortPoints(id, regions, points, sort_method, my_regions);

	////////////////  If print out end_points only /////////////////////////////
	for(it_bisect = 0; it_bisect < num_bisections; it_bisect++)
	{
		bisectTriangulation(flags, world, my_regions, all_triangles, regions, points, 0);

		string postNm = to_string(points.size());
		if(id == 0)
		{
			ofstream end_pts("end_points.dat."+postNm);
			for(point_itr = points.begin(); point_itr != points.end(); ++point_itr)		
				end_pts << (*point_itr) << endl;			
			end_pts.close();
		}

		clearRegions(id, my_regions);
		sortPoints(id, regions, points, sort_method, my_regions);
	}


/*	////////////////  If print out end_points triangles point_density /////////////////////////////
	for(it_bisect = 0; it_bisect < num_bisections; it_bisect++)
	{
		bisectTriangulation(flags, world, my_regions, all_triangles, regions, points, 0);

		string postNm = to_string(points.size());

		clearRegions(id, my_regions);
		sortPoints(id, regions, points, sort_method, my_regions);
		makeFinalTriangulations(id, flags, regions, my_regions);
		printMyFinalTriangulation(world, my_regions, all_triangles, "triangles.dat."+postNm);

		if(id == 0)
		{
			ofstream end_pts("end_points.dat."+postNm);
			ofstream pt_dens("point_density.dat."+postNm);

			for(point_itr = points.begin(); point_itr != points.end(); ++point_itr)
			{
				end_pts << (*point_itr) << endl;
				pt_dens << density((*point_itr)) << endl;
			}

			end_pts.close();
			pt_dens.close();

		}
	}
*/



	return 0;
}










