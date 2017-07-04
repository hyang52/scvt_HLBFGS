/*
 * finalize_routines.h
 *
 *  Created on: Jan 13, 2017
 *  Author: Huanhuan Yang	huan2yang@outlook.com
 *
 *  This file inherits and reshapes some data structures from 
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff
 */

#ifndef FINALIZE_ROUTINES_H_
#define FINALIZE_ROUTINES_H_

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#include "loop_routines.h"

class bpt {/*{{{*/
	public:
struct bisect_hasher {/*{{{*/
	size_t operator()(const pair<int,int> &p) const {
		uint32_t hash;
		size_t i, key[2] = { (size_t)p.first, (size_t)p.second};
		for(hash = i = 0; i < sizeof(key); ++i) {
			hash += ((uint8_t *)key)[i];
			hash += (hash << 10);
			hash ^= (hash >> 6);
		}
		hash += (hash << 3);
		hash ^= (hash >> 11);
		hash += (hash << 15);
		return hash;
	}
};/*}}}*/
};/*}}}*/


void clearRegions(const int id, vector<region> &region_vec){/*{{{*/
	//Clear all of my regions
#ifdef _DEBUG
	cerr << "Clearing local regions " << id << endl;
#endif
	vector<region>::iterator region_itr;
	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		(*region_itr).points.clear();
		(*region_itr).triangles.clear();

		assert((*region_itr).points.empty());
		assert((*region_itr).triangles.empty());
	}
#ifdef _DEBUG
	cerr << "Done Clearing local regions " << id << endl;
#endif
}/*}}}*/


void printRegions(const int id, vector<pnt>& points, int sort_method, vector<region> &regions){/*{{{*/
	// This function is only for debugging purposes.
	// It's used to print out the information associated with each region
	// to verify the region information gets read correctly.
	ofstream r_out("RegionCenters.dat");
	ofstream rr_out("RegionRadii.dat");
	ofstream rp_out("RegionProperties.dat");
	double radius;
	vector<region>::iterator region_itr;

	clearRegions(id, regions);
	sortPoints(id, regions, points, sort_method, regions);

	for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
		radius = (*region_itr).radius;

		radius = acos(-radius + 1.0);
		r_out << (*region_itr).center << endl;
		rr_out << radius << endl;
		rp_out << (*region_itr).points.size() << endl;
	}

	r_out.close();
	rr_out.close();

	rp_out.close();
	clearRegions(id, regions);
}/*}}}*/


void makeFinalTriangulations(const int id, char* flags, vector<region>& regions, vector<region> &region_vec){/*{{{*/
	//Make final triangulations, triangulates regions as in the function triangulateRegions, however
	//here the triangles are push_back to regions' triangle vectors non-repeatedly
	int i;
	int vi1, vi2, vi3;
	int swp_vi;
	struct triangulateio in, out, vorout;
	tri t;
	double cradius;
	pnt ccenter;
	pnt a, b, c;
	pnt ac, bc;
	pnt x_hat, y_hat, axis;
	pnt Q;
	double s, min_dir;
	double a_dist, b_dist, c_dist;
	double a_dist_min, b_dist_min, c_dist_min;
	double dist_temp;
	double x, y, z;
	double criteria;
	vector<region>::iterator region_itr;
	vector<pnt>::iterator point_itr;
	vector<int>::iterator neighbor_itr;


#ifdef _DEBUG
	cerr << "Triangulating points (local) " << id << endl;
#endif

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		(*region_itr).triangles.clear();

		in.numberofpoints = (*region_itr).points.size();
		in.numberofpointattributes = 0;
		in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
		in.numberofsegments = 0;
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.regionlist = (double *) NULL;
		in.pointmarkerlist = (int *) NULL;

		axis = (*region_itr).center;
		min_dir = min(fabs(axis.x),min(fabs(axis.y),fabs(axis.z)));
		if(min_dir == fabs(axis.x)){
			axis.x = 1.0;
			axis.y = 0.0;
			axis.z = 0.0;
		} else if (min_dir == fabs(axis.y)){
			axis.x = 0.0;
			axis.y = 1.0;
			axis.z = 0.0;
		} else if (min_dir == fabs(axis.z)){
			axis.x = 0.0;
			axis.y = 0.0;
			axis.z = 1.0;
		}

		x_hat = (*region_itr).center.cross(axis);
		x_hat.normalize();
		y_hat = (*region_itr).center.cross(x_hat);
		y_hat.normalize();

		i = 0;
		for(point_itr = (*region_itr).points.begin(); point_itr != (*region_itr).points.end(); ++point_itr){
			s = 2.0/(*region_itr).center.dot((*point_itr) + (*region_itr).center);
			Q = s*(*point_itr) + (s - 1.0) * (*region_itr).center;
			in.pointlist[2*i] = x_hat.dot(Q);
			in.pointlist[2*i+1] = y_hat.dot(Q);
			i++;
		}

		out.pointlist = (double *)NULL;
		out.trianglelist = (int *)NULL;

		triangulate(flags,&in,&out,&vorout);

		for(i = 0; i < out.numberoftriangles; i++){
			vi1 = out.trianglelist[3*i];
			vi2 = out.trianglelist[3*i+1];
			vi3 = out.trianglelist[3*i+2];

			a = (*region_itr).points.at(vi1);
			b = (*region_itr).points.at(vi2);
			c = (*region_itr).points.at(vi3);

			if(isCcw(a,b,c)==-1){
				b = (*region_itr).points.at(vi3);
				c = (*region_itr).points.at(vi2);
			}

			vi1 = a.idx;
			vi2 = b.idx;
			vi3 = c.idx;

			try {
				circumcenter(a,b,c,ccenter);
				ccenter.normalize();
			} catch ( int n ) {
				ccenter = (a + b + c) / 3.0;
			}

			c_dist = ccenter.dotForAngle((*region_itr).center);
			c_dist_min = 10.0;

			for(neighbor_itr = (*region_itr).neighbors.begin(); neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
				if((*neighbor_itr) != (*region_itr).center.idx){
					dist_temp = ccenter.dotForAngle(regions.at((*neighbor_itr)).center);
					c_dist_min = min(c_dist_min, dist_temp);
				}
			}

			if(c_dist < c_dist_min){
				t = tri(vi1, vi2, vi3);
				(*region_itr).triangles.push_back(t);
			}
		}
		free(in.pointlist);
		free(in.regionlist);
		free(in.pointmarkerlist);
		free(out.pointlist);
		free(out.trianglelist);
	}
#ifdef _DEBUG
	cerr << "Done triangulating points (local) " << id << endl;
#endif
	return;
}/*}}}*/




void gatherFinalTriangles(const mpi::communicator& world, vector<region>& my_regions, vector<tri>& all_triangles){/*{{{*/
	//Merge all finalTriangulations made with makeMyFinalTriangulations onto master processor.
	//Master processor inserts all triangles into an unordered_set to create a list of unique triangles
	mpi::request mycomm;
	vector<tri> temp_tris_out;
	vector<tri> temp_tris_in;
	tri t;
	vector<region>::iterator region_itr;
	vector<tri>::iterator tri_itr;
	int id=world.rank();

	#ifdef _DEBUG
		cerr << " Store my final triangulation " << id << endl;
	#endif

    all_triangles.clear();

    if(id == 0){
        temp_tris_in.clear();
        for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
            for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
                t = (*tri_itr);

                all_triangles.push_back(t);
            }
        }
    } else {
        for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
            for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
                t = (*tri_itr);

                temp_tris_out.push_back(t);
            }
        }
    }

	for(int i = 1; i < world.size(); i++){
		if(id == i){
			mycomm = world.isend(0,i,temp_tris_out);
		} else if (id == 0){
			world.recv(i, i, temp_tris_in);

            for(tri_itr = temp_tris_in.begin(); tri_itr != temp_tris_in.end(); ++tri_itr){
                t = (*tri_itr);

                all_triangles.push_back(t);
            }

            temp_tris_in.clear();
		}
	}

	if(id != 0){
		mycomm.wait();
		temp_tris_out.clear();
	}

	#ifdef _DEBUG
		cerr << " Store my final triangulation done " << id << endl;
	#endif
}/*}}}*/
/*}}}*/


void printMyFinalTriangulation(const mpi::communicator& world, vector<region>& my_regions, vector<tri>& all_triangles, 
							   const string triFile="triangles.dat"){/*{{{*/
	//Merge all finalTriangulations made with makeFinalTriangulations onto master processor.
	//Master processor inserts all triangles into an unordered_set to create a list of unique triangles
	//Triangles are then written into a file triangles.dat in ccw order
	mpi::request mycomm;
	vector<tri> temp_tris_out;
	vector<tri> temp_tris_in;
	unordered_set<tri, tri::hasher> unique_tris;
	unordered_set<tri, tri::hasher>::iterator utri_itr;
	tri t;
	vector<tri>::iterator tri_itr;

	#ifdef _DEBUG
		cerr << " Print my final triangulation " << world.rank() << endl;
	#endif

	gatherFinalTriangles(world, my_regions, all_triangles);
	if(world.rank() == 0){
		ofstream tris_out(triFile);

		for(tri_itr = all_triangles.begin(); tri_itr != all_triangles.end(); ++tri_itr){
			t = (*tri_itr);
			tris_out << t << endl;
		}

		tris_out.close();
	}

	all_triangles.clear();

	#ifdef _DEBUG
		cerr << " Print my final triangulation done " << world.rank() << endl;
	#endif
}/*}}}*/


void printAllFinalTriangulation( vector<region>& regions, vector<pnt>& points, const string triFile="triangles.dat"){/*{{{*/
	//Single processor version of printMyFinalTriangulation
	vector<tri> temp_tris_out;
	vector<tri> temp_tris_in;
	unordered_set<tri, tri::hasher> unique_tris;
	unordered_set<tri, tri::hasher>::iterator utri_itr;
	tri t;
	vector<region>::iterator region_itr;
	vector<tri>::iterator tri_itr;

	for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			temp_tris_out.push_back((*tri_itr));
			unique_tris.insert((*tri_itr).sortedTri());
		}
	}

	ofstream utris_out(triFile);

	for(utri_itr = unique_tris.begin(); utri_itr != unique_tris.end(); ++utri_itr){
		t = (*utri_itr);

		if(isCcw(points.at(t.vi1),points.at(t.vi2),points.at(t.vi3))==-1){
			int swp_v;

			swp_v = t.vi2;
			t.vi2 = t.vi3;
			t.vi3 = swp_v;
		}
		utris_out << t << endl;
	}
	utris_out.close();
}/*}}}*/




/* ***** Routines for Points ***** {{{*/

#ifdef USE_NETCDF
int writeRestartFileOverwriteNC( const int it, vector<pnt> &points ) {/*{{{*/

	// set up the file and create
	static const int NC_ERR = 2;
	NcError err(NcError::verbose_nonfatal);
	NcFile grid("point_restart.nc", NcFile::Replace);
	if(!grid.is_valid()) return NC_ERR;

	// define dimensions
	NcDim *nPointsDim;
	NcDim *itDim;

	// store dimensions
	if (!(nPointsDim = grid.add_dim( "nPoints",	points.size() ))) return NC_ERR;
	if (!(itDim = grid.add_dim( "it", it ))) return NC_ERR;

	// create/populate point arrays
	double *xPoint, *yPoint, *zPoint;
	xPoint = new double[nPointsDim->size()];
	yPoint = new double[nPointsDim->size()];
	zPoint = new double[nPointsDim->size()];

	for ( int ii = 0; ii < nPointsDim->size(); ii++ ) {
		xPoint[ii] = points.at(ii).x;
		yPoint[ii] = points.at(ii).y;
		zPoint[ii] = points.at(ii).z;
	}

	// define coordinate ncvars
	NcVar *xPointVar, *yPointVar, *zPointVar;
	if (!(xPointVar = grid.add_var("xPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(yPointVar = grid.add_var("yPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(zPointVar = grid.add_var("zPoint", ncDouble, nPointsDim))) return NC_ERR;

	// write coordinate data
	if (!xPointVar->put(xPoint,nPointsDim->size())) return NC_ERR;
	if (!yPointVar->put(yPoint,nPointsDim->size())) return NC_ERR;
	if (!zPointVar->put(zPoint,nPointsDim->size())) return NC_ERR;

	// clear memory
	delete[] xPoint;
	delete[] yPoint;
	delete[] zPoint;

	// scope closure destroys NC objs
	return 0;

}/*}}}*/

int writeRestartFileRetainNC( const int it, const vector<pnt> &points ) {/*{{{*/

	// set up the file and create
	static const int NC_ERR = 2;
	NcError err(NcError::verbose_nonfatal);
	std::ostringstream int_hole;
	int_hole << it;
	std::string restart_name = "point_restart_" + int_hole.str() + ".nc";
	NcFile grid(restart_name.c_str(), NcFile::Replace);
	if(!grid.is_valid()) return NC_ERR;

	// define dimensions
	NcDim *nPointsDim;
	NcDim *itDim;

	// store dimensions
	if (!(nPointsDim = grid.add_dim( "nPoints",	points.size() ))) return NC_ERR;
	if (!(itDim = grid.add_dim( "it", it ))) return NC_ERR;

	// create/populate point arrays
	double *xPoint, *yPoint, *zPoint;
	xPoint = new double[nPointsDim->size()];
	yPoint = new double[nPointsDim->size()];
	zPoint = new double[nPointsDim->size()];

	for ( int ii = 0; ii < nPointsDim->size(); ii++ ) {
		xPoint[ii] = points.at(ii).x;
		yPoint[ii] = points.at(ii).y;
		zPoint[ii] = points.at(ii).z;
	}

	// define coordinate ncvars
	NcVar *xPointVar, *yPointVar, *zPointVar;
	if (!(xPointVar = grid.add_var("xPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(yPointVar = grid.add_var("yPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(zPointVar = grid.add_var("zPoint", ncDouble, nPointsDim))) return NC_ERR;

	// write coordinate data
	if (!xPointVar->put(xPoint,nPointsDim->size())) return NC_ERR;
	if (!yPointVar->put(yPoint,nPointsDim->size())) return NC_ERR;
	if (!zPointVar->put(zPoint,nPointsDim->size())) return NC_ERR;

	// clear memory
	delete[] xPoint;
	delete[] yPoint;
	delete[] zPoint;

	// scope closure destroys NC objs
	return 0;

}/*}}}*/
#endif

int writeRestartFileOverwriteTXT( const int it, vector<pnt>& points) {/*{{{*/

	char temp[32];
	sprintf(temp,"point_restart.dat");

	ofstream pts_out(temp);
	vector<pnt>::iterator point_itr;
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		pts_out << (*point_itr) << endl;
	}

	pts_out.close();

}/*}}}*/

int writeRestartFileRetainTXT( const int it, vector<pnt>& points) {/*{{{*/

	char temp[32];
	sprintf(temp,"point_restart_%d.dat",it);

	ofstream pts_out(temp);
	vector<pnt>::iterator point_itr;
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		pts_out << (*point_itr) << endl;
	}

	pts_out.close();

}/*}}}*/

void writePointsAsRestart(const mpi::communicator& world, vector<pnt>& n_points, vector<pnt>& points, const int it, const restart_mode_type restart_mode, const fileio_mode_type fileio_mode){/*{{{*/

	#ifdef _DEBUG
		cerr << "Writing restart file " << world.rank() << endl;
	#endif

	if(it > 0){
		gatherAllUpdatedPoints(world, n_points, points);

		if(world.rank() == 0){

		#ifdef USE_NETCDF
			cout << "restart mode " << restart_mode << " fileio_mode " << fileio_mode << endl;
			if ( restart_mode == RESTART_OVERWRITE && fileio_mode == FILEIO_NETCDF ) {
				writeRestartFileOverwriteNC( it, points );
			} else if ( restart_mode == RESTART_OVERWRITE && fileio_mode == FILEIO_TXT ) {
				writeRestartFileOverwriteTXT( it, points );
		    } else if ( restart_mode == RESTART_OVERWRITE && fileio_mode == FILEIO_BOTH ) {
		    	writeRestartFileOverwriteNC( it, points );
		    	writeRestartFileOverwriteTXT( it, points );
		    } else if ( restart_mode == RESTART_RETAIN  && fileio_mode == FILEIO_NETCDF ) {
		    	writeRestartFileRetainNC( it, points );
			} else if ( restart_mode == RESTART_RETAIN && fileio_mode == FILEIO_TXT ) {
				writeRestartFileRetainTXT( it, points );
			} else if ( restart_mode == RESTART_RETAIN && fileio_mode == FILEIO_BOTH ) {
				writeRestartFileRetainNC( it, points );
				writeRestartFileRetainTXT( it, points );
			}
		#else
			if ( restart_mode == RESTART_OVERWRITE ) {
				writeRestartFileOverwriteTXT( it, points );
			} else if ( restart_mode == RESTART_RETAIN ) {
				writeRestartFileRetainTXT( it, points );
		    }
		#endif
		}

		/*
		if(world.rank() == 0){
			ofstream pts_out(temp);
			for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
				pts_out << (*point_itr) << endl;
			}
			pts_out.close();

			if(world.size() > 1){
				world.isend(world.rank()+1, msg_restart, junk);
			}
		} else if(world.rank() == world.size()-1){
			world.recv(world.rank()-1, msg_restart, junk);
			ofstream pts_out(temp, fstream::app);

			for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
				pts_out << (*point_itr) << endl;
			}

			pts_out.close();
		} else {
			world.recv(world.rank()-1, msg_restart, junk);
			ofstream pts_out(temp, fstream::app);

			for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
				pts_out << (*point_itr) << endl;
			}

			pts_out.close();
			world.isend(world.rank()+1, msg_restart, junk);
		}
		*/
	} else {
		return;
	}

	#ifdef _DEBUG
		cerr << "Done Writing restart file " << world.rank() << endl;
	#endif
}/*}}}*/



/* ***** Bisect Edges Routines ***** {{{ */
void bisectEdges(const mpi::communicator& world, vector<tri> & all_triangles, vector<region>& regions, vector<pnt>& points, int end){/*{{{*/
	//void bisectEdges(int end)
	// This function actually performs the bisection of the point set.
	// All bisections are performed on the master process
	// INPUT:
	// 		int end - detemines if the routine prints or not.
	//
	// 	No output.
	// 	All new points are added directly into the points vector
	unordered_set<pair<int,int>, bpt::bisect_hasher> bisection_pts;
	unordered_set<pair<int,int>, bpt::bisect_hasher>::iterator bi_pts_itr;
	pair<int, int> in_pair;
	tri t;
	int vi1, vi2, vi3;
	int i;
	vector<tri>::iterator tri_itr;

	world.barrier();
	if(world.rank() == 0){
		if(!end)
			cout << "\n\n********************************************************************************************\n"
				 << "Bisecting from " << points.size() << " points to ";

		i = points.size();
		for(tri_itr = all_triangles.begin(); tri_itr != all_triangles.end(); ++tri_itr){
			vi1 = (*tri_itr).vi1;
			vi2 = (*tri_itr).vi2;
			vi3 = (*tri_itr).vi3;
			pnt new_p;

			if(vi1 < vi2){
				in_pair.first = vi1;
				in_pair.second = vi2;
				bisection_pts.insert(in_pair);
				i++;
			}

			if(vi2 < vi3){
				in_pair.first = vi2;
				in_pair.second = vi3;
				bisection_pts.insert(in_pair);
				i++;
			}

			if(vi3 < vi1){
				in_pair.first = vi3;
				in_pair.second = vi1;
				bisection_pts.insert(in_pair);
				i++;
			}
		}

		i = points.size();
		for(bi_pts_itr = bisection_pts.begin(); bi_pts_itr != bisection_pts.end(); ++bi_pts_itr){
			pnt p = 0.5*(points.at((*bi_pts_itr).first) + points.at((*bi_pts_itr).second));
			p.idx = i;
			p.isBdry = 0;
			p.normalize();
			points.push_back(p);
			i++;
		}
		if(!end)
			cout << points.size() << " points." << endl;
	} else {
		points.clear();
	}

	clearRegions(world.rank(), regions);
	world.barrier();
	mpi::broadcast(world, points, 0);
	return;
}/*}}}*/


void bisectTriangulation(char* flags, const mpi::communicator& world, vector<region>& my_regions, vector<tri>& all_triangles, vector<region>& regions, vector<pnt>& points, int output, const string fileName="bisected_points.dat"){/*{{{*/
	// void bisectTriangulation(int end)
	//
	// Sets up the required datastructures to perform the bisection with
	// bisectEdges
	//
	// Input
	if(!output)
		makeFinalTriangulations(world.rank(), flags, regions, my_regions);

	gatherFinalTriangles(world, my_regions, all_triangles);
	bisectEdges(world, all_triangles, regions, points, output);

	vector<pnt>::iterator point_itr;

	if(world.rank() == 0 && output){
		ofstream pts_out(fileName);
		for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
			pts_out << (*point_itr) << endl;
		}
		pts_out.close();
	}

	return;
}/*}}}*/
/* }}} */



#endif /* FINALIZE_ROUTINES_H_ */
