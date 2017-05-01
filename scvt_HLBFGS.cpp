/*
 * scvt_HLBFGS.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: Huanhuan Yang huan2yang@outlook.com
 *
 *  The Overlapping domain decomposition method for parallelization is following
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff
 *
 *  This program is using non-overlapping local vector for storing points,
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
ofstream itrFile;
int it_bisect;
boost::shared_ptr<mpi_timer> itr_timer;
double initial_time;
//int  my_numPts=0;
//int  num_sorts=0;
// Global constant parameters: value given in data file: parameters
int sort_method;
int div_levs;
int use_barycenter;
int num_bisections;
string itrFileName;
bool save_bisectItr;


//////////////////////////////////////////////////////////////////////////
int evalfunc(int my_N, double* my_x, double *my_prev_x, double* f, double* my_g)
{
  *f = 0;
  double my_energy;
  vector<pnt> my_points;
  vector<pnt> disj_grad;
  vector<pnt> disj_lloyd;
  vector<pnt> my_gradPts;
  vector<double> disj_bots;
  double pNm;
  pnt p;
  for(int i=0; i<my_N/3; ++i)
  {
    pNm = sqrt(my_x[3*i]*my_x[3*i] + my_x[3*i+1]*my_x[3*i+1] + my_x[3*i+2]*my_x[3*i+2]);
    my_x[3*i] /= pNm;
    my_x[3*i+1] /= pNm;
    my_x[3*i+2] /= pNm;
    p.x = my_x[3*i]; p.y = my_x[3*i+1]; p.z = my_x[3*i+2];
    p.idx = disjDistrIdx[i];
    my_points.push_back(p);
  }
  transferUpdatedPoints(world, my_regions, my_points, points);
  clearRegions(id, my_regions);
  sortPoints(id, regions, points, sort_method, my_regions);
    //for(vector<region>::iterator region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
    //    my_numPts += (*region_itr).points.size();
    //}
  //num_sorts++;
  triangulateRegions(id, flags, my_regions);
  inteEnergGrad(id, div_levs, quadr, use_barycenter, regions, my_regions, points,
          my_energy, disj_grad, disj_lloyd, disj_bots);
  mpi::reduce(world, my_energy, *f, std::plus<double>(), 0);
  mpi::broadcast(world, *f, 0);
  //if(id==0) cout << "\n f ="<< *f <<endl;

  int my_ret = transferByDisjDistrIdx(world, my_regions, points, &my_x[0], disjDistrIdx,
                   disj_grad, &my_g[0], disj_lloyd, my_lloyds, disj_bots, my_bots);
  int ret;
  mpi::reduce(world, my_ret, ret, mpi::maximum<double>(), 0);
  mpi::broadcast(world, ret, 0);

  return ret;
}



//////////////////////////////////////////////////////////////////////////
void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
  if(id==0){
    if(it_bisect==num_bisections || save_bisectItr)
    {
      itr_timer->stop();
      itrFile << iter << " " << call_iter << " " << *f << " " << *gnorm << " "
          << initial_time + itr_timer->total_time <<"\n" ;
      itr_timer->start();
    }
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
    //return;
    /*for(int j=0; j<N; ++j)
    {
      q[j] = my_lloyds[j] - x[j];
    }*/
    HLBFGS_DSCALDV(N, &my_bots[0], q);
    HLBFGS_DSCAL(N, 0.5, q);
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
  int it, i, it_restr, count, num_avePoints, num_myPoints;
  double my_l2, my_max, my_l1, glob_l2, glob_max, glob_l1;  //for Lloyd step computation
  //double proj_alpha;

    string dataFileName;
    GetPot command_line (argc, argv);
    dataFileName = command_line.follow ("parameters", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

  // constant parameters reading from dataFile
  int points_begin = dataFile("scvt/initial/initial_point_set", 3);
  int num_pts = dataFile("scvt/initial/number_of_generated_points", 12);
  sort_method = dataFile("partition/sort_method", 0);
  num_bisections = dataFile("scvt/num_bisections", 0);
  bool bisect_final = dataFile("scvt/bisect_final", false);
  int quad_rule = dataFile("scvt/integration/quadrature_rule", 0);
  use_barycenter = dataFile("scvt/integration/delaunay_triangle_center", 0);
  div_levs = dataFile("scvt/integration/division_levs", 1);
  bool save_before_bisect = dataFile("fileIO/save_before_bisect", false);
  itrFileName = dataFile("HLBFGS/itr_File", "");
  save_bisectItr = dataFile("fileIO/save_bisect_itrFile", false);
  //double max_bdryResol = dataFile("scvt/resolution/max_boundary_resolution", 40000.0);

  // Input flags for the Triangle package
  string flags_str = "QBPIOYYiz";
  flags = new char[flags_str.size()+1];
  strcpy(flags,flags_str.c_str());

  // constant parameters for HLBFGS, reading from dataFile
  double Lloyd_tol = dataFile("HLBFGS/real_tol/Lloyd_tol", 1.0);
  double Lloyd_tolPerc = dataFile("HLBFGS/real_tol/Lloyd_tolPerc", 1.0);
  int max_restart = dataFile("HLBFGS/int_tol/max_restart", 0);
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
  if(info[3]==2)  info[12] = 1;

  //Define timers for performance studies
  const int num_timers = 5;
  mpi_timer timers[num_timers];
  string global_names[num_timers] = {"Global Time", "Final Gather", "Final Triangulation", "Initialization", "Final Bisection"};
  for(i = 0; i < num_timers; i++){
    timers[i] = mpi_timer(global_names[i]);
  }
  timers[0].start(); // Global Time Timer
  timers[3].start(); // Initialization timer
  itr_timer.reset(new mpi_timer("Cumulative Time"));

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
    } else if(points_begin == 4){
        makeMCPoints_rejection(num_pts, points);
        cout << "\n" << num_pts <<" points being created with Monte Carlo with rejection." << endl;
        num_pts << points.size();
    } else if(points_begin == 5){
        makeGeneralizedSpiralPoints_rejection(num_pts, points);
        cout << "\n" << num_pts << " points being created with Generalized Spiral with rejection." << endl;
        num_pts << points.size();
    } else if(points_begin == 6){
        makeFibonacciGridPoints_rejection(num_pts, points);
        cout << "\n" << num_pts << " points being created with Fibonacci Grid with rejection." << endl;
        num_pts << points.size();
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


  if(save_bisectItr && id==0)
  {
    timers[3].stop();
    initial_time = timers[3].total_time;
    itrFile.open(itrFileName.c_str());
    itr_timer->start();
  }

  // Loop over it_bisect
  for(it_bisect = 0; it_bisect <= num_bisections; it_bisect++)
  {
    if((!save_bisectItr) && it_bisect==num_bisections && id==0)
    {
      timers[3].stop();
      initial_time = timers[3].total_time;
      itrFile.open(itrFileName.c_str());
      itr_timer->start();
    }

    // Loop to restart Lloyd-HLBFGS in case bad communication happends in transferByDisjDistrIdx
    count = 0;
    for(it_restr = 0; it_restr <= max_restart; it_restr++)
    {
      it = 0;
      /////////////////////////////////////////////////////
/*        // Lloyd iteration as starting up
      bool stop = false;
      do{
        clearRegions(id, my_regions);
        sortPoints(id, regions, points, sort_method, my_regions);
        triangulateRegions(id, flags, my_regions);
        integrateRegions(id, div_levs, quadr, use_barycenter, regions, my_regions, points, n_points);

  //      if(it > max_it_no_proj){
  //        proj_alpha = max((double)(it-max_it_no_proj), 0.0)/max((double)max_it_scale_alpha, 1.0);
  //        projectToBoundary(proj_alpha, points, boundary_points, n_points, my_regions);
  //      }

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
        if(it == 0 && Lloyd_tolPerc < 1.0)
          Lloyd_tol = glob_l2 * Lloyd_tolPerc;
        it++;
        if(id == 0)
          cout << "Lloyd itr "<< it <<": dx_l2 = "<< glob_l2 << endl;

        if(glob_l2 > Lloyd_tol){
          transferUpdatedPoints(world, my_regions, n_points, points);
        }else{
          gatherAllUpdatedPoints(world, n_points, points);
          stop = true;
        }
      }while(!stop);
*/
      /////////////////////////////////////////////////////
      // Quasi-Newton iteration when dx <= Lloyd_tol
      clearRegions(id, my_regions);
      sortPoints(id, regions, points, sort_method, my_regions);
      getDisjointIndex(regions, my_regions, disjDistrIdx);
      std::vector<double> x(disjDistrIdx.size()*3);
      for(i=0; i<disjDistrIdx.size(); ++i)
      {
        x[3*i] = points.at(disjDistrIdx[i]).x;
        x[3*i+1] = points.at(disjDistrIdx[i]).y;
        x[3*i+2] = points.at(disjDistrIdx[i]).z;
      }

      if(id==0)
        cout << "LBFGS itr: num_feval |  f_val  |  g_norm  |"<< endl;


      int ret = HLBFGS(disjDistrIdx.size()*3, mem_num, &x[0], evalfunc, 0,
               defined_update_Hessian, newiteration, parameter, info, &world);


      count += it+info[1];
      if(id==0){
        cout << "-----------Summary: totoal number of f (and g) evals: "<< it+info[1]<< endl;
        cout << "---------Cumulative totoal number of f (and g) evals: "<< count << endl;
      }
      vector<pnt> my_points;
      pnt p;
      for(i=0; i<x.size()/3; ++i)
      {
        p.x=x[3*i]; p.y=x[3*i+1]; p.z=x[3*i+2];
        p.normalize();
        p.idx = disjDistrIdx[i];
        my_points.push_back(p);
      }
      gatherAllUpdatedPoints(world, my_points, points);

      if(ret==1 && it_restr<max_restart)
      {
        if(id==0)
          cout << "\n******Warning : bad communication in transferByDisjDistrIdx. Restarting Lloyd-LBFGS......" << endl;
        continue;
      }else if(ret==1 && it_restr==max_restart){
        if(id==0){
          cout << "\n***********************************************************************************";
          cout << "\nHave restarted Lloyd-LBFGS "<<max_restart<<" times, but still get bad communication.";
          cout << "\nYou may need to use refinement from coarse mesh, or better initialization ! " << endl;
        }
        exit(1);
      }else
      {
        timers[0].stop(); // Global Time Timer
        if(!save_bisectItr)
          timers[3].stop();
        if(save_before_bisect && it_bisect < num_bisections)
        {
          string postNm = to_string(points.size());

          clearRegions(id, my_regions);
          sortPoints(id, regions, points, sort_method, my_regions);
          makeFinalTriangulations(id, flags, regions, my_regions);
          printMyFinalTriangulation(world, my_regions, all_triangles, "triangles.dat."+postNm);

          if(id == 0)
          {
            ofstream end_pts("end_points.dat."+postNm);
            ofstream pt_dens("point_density.dat."+postNm);
            //ofstream bdry_pts("boundary_points.dat");
            for(point_itr = points.begin(); point_itr != points.end(); ++point_itr)
            {
              end_pts << (*point_itr) << endl;
              pt_dens << density((*point_itr)) << endl;
            }
            //for(boundary_itr = boundary_points.begin(); boundary_itr != boundary_points.end(); boundary_itr++){
            //  bdry_pts << (*boundary_itr) << endl;
            //}
            end_pts.close();
            pt_dens.close();
            //bdry_pts.close();
          }
        }
        timers[0].start(); // Global Time Timer
        if(!save_bisectItr)
          timers[3].start();
        break;
      }
    }
    if(it_bisect==num_bisections)
      itrFile.close();

    // Bisect if needed
    if(it_bisect < num_bisections)
    {
      if(!save_before_bisect)
      {
            clearRegions(id, my_regions);
            sortPoints(id, regions, points, sort_method, my_regions);
      }
      bisectTriangulation(flags, world, my_regions, all_triangles, regions, points, 0);
    } else {
      if(id == 0 && num_bisections>0){
        cout << "\nNo more bisections for convergence" << endl;
      }
    }
  }
  timers[0].stop();

  // Compute final triangulation by merging all triangulations from each processor into an
  // unordered_set, and then ordering them ccw before printing them out.
  // write triangles to triangles.dat
  timers[2].start(); // Final Triangulation Timer
  clearRegions(id, my_regions);
  sortPoints(id, regions, points, sort_vor, my_regions);
  makeFinalTriangulations(id, flags, regions, my_regions);
  timers[2].stop();

  printMyFinalTriangulation(world, my_regions, all_triangles, "triangles.dat."+itrFileName);

  if(id == 0){
    ofstream end_pts("end_points.dat."+itrFileName);
    ofstream pt_dens("point_density.dat."+itrFileName);
    //ofstream bdry_pts("boundary_points.dat");
    for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
      end_pts << (*point_itr) << endl;
      pt_dens << density((*point_itr)) << endl;
    }
    //for(boundary_itr = boundary_points.begin(); boundary_itr != boundary_points.end(); boundary_itr++){
    //  bdry_pts << (*boundary_itr) << endl;
    //}
    end_pts.close();
    pt_dens.close();
    //bdry_pts.close();
  }

/*  //Compute average points per region for diagnostics
    cout << "\nAverage points per iteration: " << my_numPts/num_sorts << endl;
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
*/
  //Bisect all edges of all triangles to give an extra point set at the end, SaveVertices
  timers[4].start(); // Final Bisection Timer
  if(bisect_final)
    bisectTriangulation(flags, world, my_regions, all_triangles, regions, points, 1, "bisected_points.dat."+itrFileName);
  timers[4].stop();

  //Print out final timers, for global times.
  if(id == 0){
    cout << endl << " ---- Final Timers ---- " << endl;
    for(i = 0; i < num_timers; i++){
      cout << timers[i];
    }

  }

  return 0;
}










