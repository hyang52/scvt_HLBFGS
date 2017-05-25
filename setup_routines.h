/*
 * setup_routines.h
 *
 *  Created on: Jan 13, 2017
 *  Author: Huanhuan Yang   huan2yang@outlook.com
 *
 *  This file inherits and reshapes some data structures from
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff
 */

#ifndef SETUP_ROUTINES_H_
#define SETUP_ROUTINES_H_

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include "triangulation.h"
#include "densities.h"

using namespace std;
using namespace tr1;
namespace mpi = boost::mpi;

class region{/*{{{*/
    private:
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive & ar, const unsigned int version)
            {
                ar & center;
                ar & radius;
                ar & input_radius;
                ar & triangles;
                ar & neighbors;
                ar & neighbors1;
                ar & neighbors2;
                ar & boundary_points;
                ar & loop_start;
                ar & loop_stop;
            }

    public:
        pnt center;
        double radius;
        double input_radius;
        vector<pnt> points;
		vector<double> overlap_dists;
        vector<tri> triangles;
        vector<int> neighbors; // First Level of Neighbors
        vector<int> neighbors1; // First Level of Neighbors + Self
        vector<int> neighbors2; // Second Level of Neighbors + First Level of Neighbors
        vector<pnt> boundary_points;
        vector<int> loop_start; // beginning point in loop
        vector<int> loop_stop; // ending point in loop
};/*}}}*/

struct int_hasher {/*{{{*/
      size_t operator()(const int v) const { return v; }
};/*}}}*/

class mpi_timer{/*{{{*/
    public:
        mpi::timer my_timer;
        double total_time;
        int num_calls;
        string name;

        mpi_timer() : total_time(0), num_calls(0), name("Default") { };
        mpi_timer(string in_name) : total_time(0), num_calls(0), name(in_name) { };

        mpi_timer operator+(const mpi_timer &t) const {/*{{{*/
            mpi_timer nt;
            nt.init("Sum");
            nt.num_calls = t.num_calls;
            nt.total_time = t.total_time + total_time;

            return nt;
        }/*}}}*/
        void init(string in_name){/*{{{*/
            total_time = 0.0;
            num_calls = 0;
            name = in_name;
        }/*}}}*/
        void init(){/*{{{*/
            total_time = 0;
            num_calls = 0;
        }/*}}}*/
        void start(){/*{{{*/
            my_timer.restart();
        }/*}}}*/
        void stop(){/*{{{*/
            total_time += my_timer.elapsed();
            num_calls++;
        }/*}}}*/
};/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const mpi_timer &t){/*{{{*/
    if(t.num_calls > 0){
    os << t.name << ": " << t.total_time*1e3 << " (ms), " << (t.total_time/t.num_calls)*1e3 << " (ms). Called " << t.num_calls << " times." << endl;
    } else {
    os << t.name << ": Never called. Time = 0.0 (ms)" << endl;
    }
    return os;
}/*}}}*/

// Message tags for sending and receiving non-blocking messages
enum {msg_points, msg_tri_print , msg_restart, msg_ave, msg_max, msg_l1};
// Sort types
enum {sort_dot, sort_vor, sort_vorMD};
enum restart_mode_type { RESTART_OVERWRITE, RESTART_RETAIN };
enum fileio_mode_type { FILEIO_TXT, FILEIO_NETCDF, FILEIO_BOTH };


/* ***** Setup Routines ***** {{{*/

void readBoundaries(const double max_resolution, vector<pnt>& boundary_points){/*{{{*/
    int count_count, bdry_count, fill_count, bdry_total;
    int add_count;
    int count_start, count_stop;
    int cur_loop_start, cur_loop_length;
    vector<int> loop_start;
    vector<int> loop_stop;
    int p0_idx, p1_idx;
    //pnt p;
    //pnt p_b, p_e;
    double bdry_lon, bdry_lat;
    double dtr = M_PI/180.0;
    double point_delta;
    double add_spacing;
    double denom;
    double c0, c1;
    double t, omega;
    double r_earth = 6371.0; //avg radius in km

    // gw: read boundary points file
    bdry_count = 0;
    ifstream bdry_in("SaveBoundaries");
    if(!bdry_in)
        return;

    while(!bdry_in.eof()){
        bdry_in >> bdry_lon >> bdry_lat;
        bdry_in.ignore(10000,'\n');

        if(bdry_in.good()){

            bdry_lon *= dtr;
            bdry_lat *= dtr;
            pnt p = pntFromLatLon(bdry_lat, bdry_lon);

            p.normalize();
            p.idx = bdry_count;
            p.isBdry = 0;
            bdry_count++;
            boundary_points.push_back(p);
        } // end if input good

    } // end while not eof
    bdry_in.close();

    // gw: read loop counts file
    count_count = 0;
    ifstream count_in("SaveLoopCounts");
    if(!count_in)
        return;

    while(!count_in.eof()){
        count_in >> count_start >> count_stop;
        count_in.ignore(10000,'\n');

        if(count_in.good()){
            loop_start.push_back(count_start);
            loop_stop.push_back(count_stop);
            count_count++;

        } // end if input good
    } // end while not eof
    count_in.close();
    // gw: loop over loops
    fill_count = 0;
    for(int cur_loop = 0; cur_loop < count_count; cur_loop++){
        cur_loop_start = loop_start.at(cur_loop) - 1;
        cur_loop_length = loop_stop.at(cur_loop) - cur_loop_start;

        // gw: loop over point pairs in current loop
        for(int cur_pair = 0; cur_pair < cur_loop_length; cur_pair++){

            p0_idx = cur_loop_start + cur_pair;
            p1_idx = cur_loop_start + (cur_pair+1)%cur_loop_length;
            pnt p0 = boundary_points.at(p0_idx);
            pnt p1 = boundary_points.at(p1_idx);
            point_delta = p1.dotForAngle(p0);

            // gw: if distance between pair is greater than allowed amount
            //     then add some additional points
            if ( (point_delta * r_earth) > max_resolution ) {

                // gw: figure out how many points to added
                add_count = (int)ceil( (point_delta * r_earth) ) / max_resolution;
                add_spacing = 1.0 / ((double)add_count + 1);
                denom = sin( point_delta );
                bdry_lat = p1.getLat() - p0.getLat();
                bdry_lon = p1.getLon() - p0.getLon();

                // gw: loop for adding point(s)
                for(int cur_add = 1; cur_add <= add_count; cur_add++){

                    /* Slerp - Doesn't work for constant lat, lon curves
                    // http://en.wikipedia.org/wiki/Slerp
                    t = add_spacing * ( (double)cur_add + 1 );
                    c0 = sin( (1.0 - t) * point_delta ) / denom;
                    c1 = sin( t * point_delta ) / denom;
                    pnt temp_point = c0 * p0 + c1 * p1; // */

                    // Linear interpolation, using Lat Lon.
                    pnt temp_point = pntFromLatLon( p0.getLat() + cur_add * add_spacing * bdry_lat,
                                                    p0.getLon() + cur_add * add_spacing * bdry_lon);

                    temp_point.normalize();
                    temp_point.idx = bdry_count + fill_count;
                    boundary_points.push_back(temp_point);

                    fill_count++;
                }

            } // end if need to add fill points

        } // end loop over point pairs in current loop

    } // end loop over loops
    cout << "Read in " << bdry_count << " boundary points." << endl;
    cout << "Made " << fill_count << " fill points." << endl;
    cout << "There are " << boundary_points.size() << " boundary points total." << endl;

    //num_bdry = boundary_points.size();

}/*}}}*/


int buildRegions(const int id, vector<region>& regions, const string regList="RegionList", const string regConn="RegionTriangulation"){/*{{{*/
    //Read in region centers, and connectivity (triangulation) from files
    //RegionList, and RegionTriangulation
    //From these, regions are setup to have a center and a radius.
    ifstream region_list(regList);
    ifstream region_connectivity(regConn);
    vector<region>::iterator region_itr;
    vector<int>::iterator neighbor_itr;
    vector<unordered_set<int, int_hasher> > region_neighbors;
    unordered_set<int, int_hasher>::iterator region_neigh_itr;
    unordered_set<int, int_hasher> neighbors1, neighbors2;
    unordered_set<int, int_hasher>::iterator neigh_itr;
    vector<int>::iterator cur_neigh_itr;
    region r;
    pnt p;
    double loc_radius, max_radius;
    double alpha, beta;
    int min_connectivity;
    int t1, t2, t3;
    int i;

    alpha = 2.0/3.0;
    beta = 4;

#ifdef _DEBUG
    cerr << "Building regions " << id << endl;
#endif


    if(!region_list){
        cout << "\nFailed to open file " << regList << endl;
        return 1;
    }

    // Read in region centers
    i = 0;
    while(!region_list.eof()){
        region_list >> p;
        p.idx = i;
        i++;
        p.normalize();
        r.center = p;
        r.radius = 0.0;

        if(region_list.good()){
            regions.push_back(r);
        }
    }
    region_list.close();

    region_neighbors.resize(regions.size());

    if(!region_connectivity){
        cout << "\nFailed to open file " << regConn << endl;
        return 1;
    }

    //Read in region triangulation, and insert into hash table to get unique connectivity.
    min_connectivity = regions.size();
    while(!region_connectivity.eof()){
        region_connectivity >> t1 >> t2 >> t3;

        if(t1 < min_connectivity)
            min_connectivity = t1;
        if(t2 < min_connectivity)
            min_connectivity = t2;
        if(t3 < min_connectivity)
            min_connectivity = t3;
    }
    region_connectivity.close();

    region_connectivity.open(regConn);
    while(!region_connectivity.eof()){
        region_connectivity >> t1 >> t2 >> t3;

        t1 = t1 - min_connectivity;
        t2 = t2 - min_connectivity;
        t3 = t3 - min_connectivity;

        if(region_connectivity.good()){
            region_neighbors[t1].insert(t2);
            region_neighbors[t1].insert(t3);

            region_neighbors[t2].insert(t1);
            region_neighbors[t2].insert(t3);

            region_neighbors[t3].insert(t1);
            region_neighbors[t3].insert(t2);
        }
    }
    region_connectivity.close();

    //Compute region radii by dotting with each neighbor, and taking the max distance.
    for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
        max_radius = 0.0;
        loc_radius = 0.0;

        neighbors1.insert((*region_itr).center.idx);
        for(region_neigh_itr = region_neighbors[(*region_itr).center.idx].begin();
                    region_neigh_itr != region_neighbors[(*region_itr).center.idx].end();
                    ++region_neigh_itr){

            loc_radius = (*region_itr).center.dotForAngle(regions[(*region_neigh_itr)].center);
            //cerr << "  loc_radius = " << loc_radius << endl;
            if(loc_radius > max_radius){
                max_radius = loc_radius;
            }
            if((*region_neigh_itr) != (*region_itr).center.idx)
            	(*region_itr).neighbors.push_back((*region_neigh_itr));
        }

        (*region_itr).radius = std::min(max_radius,M_PI);
        (*region_itr).input_radius = std::min(max_radius,M_PI);

        (*region_itr).points.clear();
        (*region_itr).triangles.clear();
        (*region_itr).boundary_points.clear();
    }

    // Build first and second levels of neighbors, for use in more complicated sort method.
    for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
        neighbors1.clear();
        neighbors2.clear();

        neighbors1.insert((*region_itr).center.idx);
        //neighbors2.insert((*region_itr).center.idx);
        for(cur_neigh_itr = (*region_itr).neighbors.begin(); cur_neigh_itr != (*region_itr).neighbors.end(); ++cur_neigh_itr){
            neighbors1.insert((*cur_neigh_itr));
            neighbors2.insert((*cur_neigh_itr));

            for(neighbor_itr = regions.at((*cur_neigh_itr)).neighbors.begin();
                neighbor_itr != regions.at((*cur_neigh_itr)).neighbors.end(); ++neighbor_itr){
                if((*neighbor_itr)!=(*region_itr).center.idx)
                    neighbors2.insert((*neighbor_itr));
            }
        }

        for(neigh_itr = neighbors1.begin(); neigh_itr != neighbors1.end(); ++neigh_itr){
            (*region_itr).neighbors1.push_back((*neigh_itr));
        }

        for(neigh_itr = neighbors2.begin(); neigh_itr != neighbors2.end(); ++neigh_itr){
            (*region_itr).neighbors2.push_back((*neigh_itr));
        }
    }

    region_neighbors.clear();

    return 0;

}/*}}}*/



int buildOverlap(const int id, const vector<pnt>& points, vector<region>& regions, vector<region> &region_vec){
    vector<region>::iterator region_itr;
    vector<pnt>::const_iterator point_itr;
    vector<int>::iterator neighbor_itr;

    double dist_temp;
    double my_dist;
    double min_dist;
    int min_region;
	int my_DisjNumPts;
	pnt neibPt, midPt, tmp, overlap_nbPt; 
	double shifted_dist;
	double coeff;

    for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		my_DisjNumPts = 0;
        for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
            min_dist = M_PI;
            my_dist = (*point_itr).dotForAngle((*region_itr).center);

            for(neighbor_itr = (*region_itr).neighbors.begin();
                    neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
                dist_temp = (*point_itr).dotForAngle(regions.at((*neighbor_itr)).center);

                if(dist_temp < min_dist){
                    min_region = (*neighbor_itr);
                    min_dist = dist_temp;
                }
            }
            if(my_dist < min_dist || (my_dist == min_dist && (*region_itr).center.idx < min_region)){
				my_DisjNumPts++;
            }
        }
		(*region_itr).overlap_dists.clear();
		//(*region_itr).overlap_nbPts.clear();
        for(neighbor_itr = (*region_itr).neighbors.begin();
                neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
			/*neibPt = regions.at((*neighbor_itr)).center;
			midPt = ( (*region_itr).center + neibPt ) / 2.0;			
			try{
				midPt.normalize();
			}catch(int j){
				cout << "throwing error from zero norm, replacing by polar coor for middle point" << endl;
				double lat = ( (*region_itr).center.getLat() + neibPt.getLat() )/2.0;
				double lon = ( (*region_itr).center.getLon() + neibPt.getLon() )/2.0;
				midPt = pntFromLatLon(lat,lon);
			}
			tmp = midPt - (*region_itr).center;
			tmp = tmp - midPt * tmp.dot(midPt);
			overlap_nbPt = midPt + tmp*shifted_dist;
			overlap_nbPt.normalize();
			(*region_itr).overlap_nbPts.push_back(overlap_nbPt);
			*/
			if(my_DisjNumPts <= 1600)			
				coeff = 1.0 - sqrt(my_DisjNumPts)/80.0;
			else
				coeff = pow(min(density((*region_itr).center)/density(midPt),16.0), 1.0/4.0) * 20.0/ sqrt(my_DisjNumPts);
			shifted_dist = neibPt.dotForAngle((*region_itr).center)
						   * min(coeff, 1.0);
			//cout << "shifted_dist = " << shifted_dist << ", radius = " << (*region_itr).input_radius << endl;
			(*region_itr).overlap_dists.push_back(shifted_dist);
        }

    }

}/*}}}*/




/* ***** Generic Region Routines ***** {{{ */
void sortPoints(const int id, vector<region>& regions, const vector<pnt>& points, int sort_type, vector<region> &region_vec){/*{{{*/
    //Sort points into my region(s).
    //This is done using a dot product and checking if the dot product is inside of the region radius
    double val;
    vector<region>::iterator region_itr;
    vector<pnt>::const_iterator point_itr;
    vector<int>::iterator neighbor_itr;

#ifdef _DEBUG
    cerr << "Sorting Points " << id << endl;
#endif
    if(sort_type == sort_dot){
        //Simple Dot Product Sort using the radius based decomposition.
        for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
            (*region_itr).radius = (*region_itr).input_radius;
            for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
                val = (*point_itr).dotForAngle((*region_itr).center);

                if(val < (*region_itr).radius){
                    (*region_itr).points.push_back((*point_itr));
                }
            }
    //      cerr << ": rad = " << (*region_itr).radius << ", val =" << val << endl;
        }
    } else if (sort_type == sort_vor){
        //More complicated sort, that sorts by Voronoi cells keeping current regions points, as well as neighboring regions points.
        //Should handle variable resolution meshes better than the more simple dot product sorting.
        double my_dist;
        int added;
        double min_dist;
        double max_dist;
        int min_region;
        vector<int>::iterator cur_neigh_itr;

        for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
            max_dist = 0.0;
            for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
                min_dist = M_PI;
                my_dist = (*point_itr).dotForAngle((*region_itr).center);

                //if(my_dist < (*region_itr).input_radius)
                {
                    for(neighbor_itr = (*region_itr).neighbors2.begin();
                            neighbor_itr != (*region_itr).neighbors2.end(); ++neighbor_itr){

                        val = (*point_itr).dotForAngle(regions.at((*neighbor_itr)).center);

                        if(val < min_dist){
                            min_region = (*neighbor_itr);
                            min_dist = val;
                        }
                    }

                    added = 0;

                    if(my_dist < min_dist){
                        (*region_itr).points.push_back((*point_itr));
                        added = 1;
                    }
                    for(neighbor_itr = (*region_itr).neighbors1.begin();
                            neighbor_itr != (*region_itr).neighbors1.end() && added == 0;
                            ++neighbor_itr){
                        if(min_region == (*neighbor_itr)){
                            //val = (*region_itr).center.dotForAngle(regions[min_region].center);
                            //if(my_dist < val)
                            {
                                (*region_itr).points.push_back((*point_itr));
                                added = 1;
                            }
                        }
                    }
                }
                if(my_dist > max_dist && added==1){
                    max_dist = my_dist;
                }

            }
            (*region_itr).radius = max_dist;
        }
    } else if (sort_type == sort_vorMD){
        double my_dist;
        double min_dist;

        for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
            for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
                min_dist = M_PI;
                my_dist = (*point_itr).dotForAngle((*region_itr).center);
            
				for(int i=0; i<(*region_itr).neighbors.size(); ++i){
                    /*val = (*point_itr).dotForAngle((*region_itr).overlap_nbPts[i]);
                    if(val < min_dist){
                        min_dist = val;
                    }*/
                    val = (*point_itr).dotForAngle(regions.at((*region_itr).neighbors[i]).center) + (*region_itr).overlap_dists[i];
                    if(val < min_dist){
                        min_dist = val;
                    }
				}
				min_dist = min(min_dist,M_PI*0.75);
                if(my_dist <= min_dist){
                    (*region_itr).points.push_back((*point_itr));
                }

            }
			if((*region_itr).points.size()==0){
				for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		            min_dist = M_PI;
		            my_dist = (*point_itr).dotForAngle((*region_itr).center);
		        
					for(int i=0; i<(*region_itr).neighbors.size(); ++i){
		                /*val = (*point_itr).dotForAngle((*region_itr).overlap_nbPts[i]);
		                if(val < min_dist){
		                    min_dist = val;
		                }*/
		                val = (*point_itr).dotForAngle(regions.at((*region_itr).neighbors[i]).center) + (*region_itr).overlap_dists[i];
		                if(val < min_dist){
		                    min_dist = val;
		                }
					}
					min_dist = min(min_dist,M_PI*0.75);
		            if(my_dist <= min_dist){
		                (*region_itr).points.push_back((*point_itr));
		            }
            	}
			}

        }
    }
#ifdef _DEBUG
    cerr << "Done Sorting Points (Local) " << id << endl;
#endif
    return;
}/*}}}*/


void sortBoundaryPoints(const int id, vector<region>& regions, vector<pnt>& boundary_points, vector<region> &region_vec){/*{{{*/
    //Sort points into my region(s).
    //This is done using a dot product and checking if the dot product is inside of the region radius
    double val, min_val;
    int index;
    vector<pnt>::iterator point_itr;
    vector<region>::iterator region_itr;
#ifdef _DEBUG
    cerr << "Sorting Boundary Points " << id << endl;
#endif
    //simple dot product sort using the radius based decomposition.
    for(point_itr = boundary_points.begin(); point_itr != boundary_points.end(); point_itr++){
        min_val = M_PI;
        for(region_itr = regions.begin(); region_itr != regions.end(); region_itr++){
            val = (*point_itr).dotForAngle((*region_itr).center);

            if(val < min_val){
                min_val = val;
                index = (*region_itr).center.idx;
            }
        }

        for(region_itr = region_vec.begin(); region_itr != region_vec.end(); region_itr++){
            if((*region_itr).center.idx == index){
                (*region_itr).boundary_points.push_back((*point_itr));
            }
        }
    }
#ifdef _DEBUG
    cerr << "Done Sorting Boundary Points (Local) " << id << endl;
#endif
    return;
}/*}}}*/


#endif /* SETUP_ROUTINES_H_ */
