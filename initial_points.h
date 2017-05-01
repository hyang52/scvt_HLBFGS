/*
 * initial_routines.h
 *
 *  Created on: Jan 13, 2017
 *      Author: Huanhuan Yang	huan2yang@outlook.com
 *      Update by: Lindley Graham
 *
 *  This file inherits and reshapes some data structures from 
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff
 */

#ifndef INITIALIZE_ROUTINES_H_
#define INITIALIZE_ROUTINES_H_

#include <iostream>
#include <fstream>
#include "triangulation.h"
#include "densities.h"

#define SEED	3729

using namespace std;

/* ***** Point Init Routines ***** {{{*/
void readPoints(int& n, vector<pnt>& points){/*{{{*/
	//Read in initial point set from SaveVertices
	ifstream points_in("SaveVertices");
	pnt p;
	int i;

	i = 0;
	while(!points_in.eof()){
		points_in >> p;
		points_in.ignore(10000,'\n');
		p.idx = i;
		p.isBdry = 0;
		p.normalize();

		if(points_in.good()){
			points.push_back(p);
		}
		i++;
	}

	n = points.size();

	cout << "Read in " << n << " points from SaveVertices" << endl;

}/*}}}*/


void makeMCPoints(int n, vector<pnt>& points){/*{{{*/
	//Create Monte Carlo random point set
	int i, j;
	srand48(time(NULL));
	double dlon, dlat;
	double dtr;
	double lat, lon;
	double x, y, z;
	double dens_check, dens_comp;
	pnt p;

	dtr = M_PI/180.0;

	// Uniform Spherical Points
	for(i = 0; i < n; i++){
		x = drand48()*2.0-1.0;
		y = drand48()*2.0-1.0;
		z = drand48()*2.0-1.0;

		p = pnt(x,y,z);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();

		points.push_back(p);
	}

	// */
	/* Uniform points on equator between -20 and +20 latitude
	lat = 20.0 * dtr;
	dlon = 2.0*M_PI/(n/10);
	j = 0;
	for(i = 0; i < (n / 10); i++){
		lon = i*dlon;
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
	}
	lat = -1.0*lat;
	for(i = n/10; i < (2 * n / 10); i++){
		lon = (i-(n/10))*dlon;
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
	}
	for(i = (2*n/10); i < n; i++){
		lon = drand48()*2*M_PI;
//		do{
			lat = (drand48()*2.0 - 1.0)*(15.0*dtr);
//		}while(!(fabs(lat) < 15.0*dtr));
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 0;
		points.push_back(p);
		j++;
	}
	// */
	/* Uniform points in rectangle, from 0 to 45 north lat, and 0 to 90 lon
	lat = 0.0;
	dlon = 90.0*dtr/(n/10);
	lon = 0.0;
	j = 0;
	for(i = 0; i < (n/10)+1; i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lon += dlon;
	}
	lat = 45.0*dtr;
//	dlon = 90.0*dtr/(2*n/30);
	lon = 0.0;
	for(i = 0; i < (n/10)+1; i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lon += dlon;
	}
	lon = 0.0;
	dlat = 45.0*dtr/(n/20);
	lat = dlat;
	for(i = 1; i < (n/20); i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lat += dlat;
	}
	lon = M_PI/2.0;
	lat = dlat;
	for(i = 1; i < (n/20); i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lat += dlat;
	}
	for(i = j; i < n; i++){
		lon = (drand48()*0.85 + 0.05)*M_PI/2.0;
		lat = (drand48()*0.85 + 0.05)*M_PI/4.0;
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = i;
		p.isBdry = 0;
		points.push_back(p);
	}
	// */
	/* Uniform points in a triangle from 1,0,0, 0,1,0, and 0,0,1
	lat = 0.0;
	dlon = (M_PI/2.0)/(n/10);
	j = 0;
	for(i = 0; i < (n/10)+1; i++){
		lon = i * dlon;
		p = pntFromLatLon(lat,lon);
		p.idx = j;
		p.isBdry = 2;
		p.normalize();
		points.push_back(p);
		j++;
	}
	lon = 0.0;
	dlat = (M_PI/2.0)/(n/10);
	lat = dlat;
	for(i = 1; i < (n/10)+1; i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 2;
		points.push_back(p);
		lat += dlat;
		j++;
	}
	lon = M_PI/2.0;
	lat = dlat;
	for(i = 1; i < (n/10); i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 2;
		points.push_back(p);
		lat += dlat;
		j++;
	}
	for(i = j; i < n; i++){
		lat = (drand48()*0.85 + 0.05)*M_PI/2.0;
		lon = (drand48()*0.85 + 0.05)*M_PI/2.0;
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = i;
		p.isBdry = 0;
		points.push_back(p);
	}
	// */
	cout << "Created " << points.size() << " points using monte carlo." << endl;
}/*}}}*/

void makeMCPoints_rejection(int n, vector<pnt>& points){
	//Create Monte Carlo random point set
	int i, j;
	srand48(time(NULL));
	double dlon, dlat;
	double dtr;
	double lat, lon;
	double x, y, z;
	double dens_check, dens_comp;
	pnt p;

	dtr = M_PI/180.0;

	// Non-Uniform Spherical Points w/ rejection sampling
	i = 0;
	do{
		x = drand48()*2.0-1.0;
		y = drand48()*2.0-1.0;
		z = drand48()*2.0-1.0;

		p = pnt(x,y,z);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();

		dens_check = density(p);
		dens_comp = drand48() * density_max;

		if (dens_comp <= dens_check)
		{
			points.push_back(p);
			++i;
		}
		}
	while(i < n);

	cout << "Created " << points.size() << " points using monte carlo." << endl;
	}



void makeGeneralizedSpiralPoints(int n, vector<pnt>& points){/*{{{*/
	//Create Generalize Spiral point set
	int i;
	int idx;
	pnt p;
	double phi_curr, h, theta, aa, bb, cc;
	double gsC = 3.809;
	double twopi_dp;

	twopi_dp = 2.0*M_PI;

	p = pnt(0.0,0.0,0.0,0,0);

	//	first pt, loop primer
	i = 0;
	h = -1.0;
	theta = acos(h);
	phi_curr = 0;

	p.x = cos( phi_curr ) * sin( theta );
	p.y = sin( phi_curr ) * sin( theta );
	p.z = cos( theta );
	p.idx = i;
	p.isBdry = 0;
	points.push_back(p);

	for(i = 1; i < n-1; i++){
		h = -1.0 + (2.0*(double)(i))/(double)(n-1);
		theta = acos(h);
		aa = phi_curr;
		bb = gsC / sqrt((double)n);
		cc = 1.0 / sqrt(1.0-h*h);
		phi_curr =  fmod(aa + bb * cc,twopi_dp);

		p.x = cos( phi_curr ) * sin( theta );
		p.y = sin( phi_curr ) * sin( theta );
		p.z = cos( theta );
		p.idx = i;
		p.isBdry = 0;
		points.push_back(p);
	}

	p.x = 0.0;
	p.y = 0.0;
	p.z = 1.0;
	p.idx = n-1;
	p.isBdry = 0;
	points.push_back(p);

	cout << "Created " << points.size() << " points using generalized spiral." << endl;
}/*}}}*/


void makeGeneralizedSpiralPoints_rejection(int n, vector<pnt>& points){/*{{{*/
	//Create Generalize Spiral point set
	int i;
	int idx;
	int n_accept = 0;
	pnt p;
	double phi_curr, h, theta, aa, bb, cc, ratio;
	double gsC = 3.809;
	double twopi_dp;
	double dens_check, dens_comp;

	twopi_dp = 2.0*M_PI;

	// Do an inital sampling to estiamte the accept/reject ratio
	//	first pt, loop primer
	p = pnt(0.0,0.0,0.0,0,0);

	i = 0;
	h = -1.0;
	theta = acos(h);
	phi_curr = 0;

	p.x = cos( phi_curr ) * sin( theta );
	p.y = sin( phi_curr ) * sin( theta );
	p.z = cos( theta );
	p.idx = i;
	p.isBdry = 0;
	points.push_back(p);
	dens_check = density(p);
	dens_comp = drand48() * density_max;

	if (dens_comp <= dens_check)
		++n_accept;

	for(i = 1; i < n-1; i++){
		h = -1.0 + (2.0*(double)(i))/(double)(n-1);
		theta = acos(h);
		aa = phi_curr;
		bb = gsC / sqrt((double)n);
		cc = 1.0 / sqrt(1.0-h*h);
		phi_curr =  fmod(aa + bb * cc,twopi_dp);

		p.x = cos( phi_curr ) * sin( theta );
		p.y = sin( phi_curr ) * sin( theta );
		p.z = cos( theta );
		p.idx = i;
		p.isBdry = 0;

		dens_check = density(p);
		dens_comp = drand48() * density_max;

		if (dens_comp <= dens_check)
			++n_accept;
	}


	p.x = 0.0;
	p.y = 0.0;
	p.z = 1.0;
	p.idx = n-1;
	p.isBdry = 0;
	dens_check = density(p);
	dens_comp = drand48() * density_max;

	if (dens_comp <= dens_check)
		++n_accept;

	// Increase the number of points based on an approximate accept/reject ratio
	// estiamte the accept/reject ratio
	ratio = (double) n / (double) n_accept;
	cout << "Acceptance ratio estimate " << ratio << endl;
	cout << "Goal " << n << endl;
	cout << "Accepted " << n_accept << endl;
	n *= (int) ratio * 1.5;
	// Sample generalized sphere points accordinly
	//	first pt, loop primer
	p = pnt(0.0,0.0,0.0,0,0);

	i = 0;
	h = -1.0;
	theta = acos(h);
	phi_curr = 0;

	p.x = cos( phi_curr ) * sin( theta );
	p.y = sin( phi_curr ) * sin( theta );
	p.z = cos( theta );
	p.idx = i;
	p.isBdry = 0;
	dens_comp = drand48() * density_max;

	if (dens_comp <= dens_check)
		points.push_back(p);



	for(i = 1; i < n-1; i++){
		h = -1.0 + (2.0*(double)(i))/(double)(n-1);
		theta = acos(h);
		aa = phi_curr;
		bb = gsC / sqrt((double)n);
		cc = 1.0 / sqrt(1.0-h*h);
		phi_curr =  fmod(aa + bb * cc,twopi_dp);

		p.x = cos( phi_curr ) * sin( theta );
		p.y = sin( phi_curr ) * sin( theta );
		p.z = cos( theta );
		p.idx = i;
		p.isBdry = 0;

		dens_check = density(p);
		dens_comp = drand48() * density_max;

		if (dens_comp <= dens_check)
			points.push_back(p);
	}


	p.x = 0.0;
	p.y = 0.0;
	p.z = 1.0;
	p.idx = n-1;
	p.isBdry = 0;
	dens_comp = drand48() * density_max;

	if (dens_comp <= dens_check)
		points.push_back(p);


	cout << "Created " << points.size() << " points using generalized spiral with rejection." << endl;
}/*}}}*/

void makeFibonacciGridPoints_rejection(int n, vector<pnt>& points){/*{{{*/
	const double g_ratio = (1.0+sqrt(5)) / 2.0;
	double lambda, phi, sinphi, x, y, z, ratio;
	int i, j, m;
	int n_accept = 0;
	pnt p;
	double dens_check, dens_comp;

	m = n/2;

	// Do an inital sampling to estiamte the accept/reject ratio
	for (i = 0; i < n; i++){
		j = i - m;

		lambda = i / g_ratio;
		sinphi = (2.0*j) / (n + 1);
		phi = asin(sinphi);
		x = cos(phi) * sin(lambda);
		y = cos(phi) * cos(lambda);
		z = sin(phi);
		p = pntFromLatLon(phi, lambda);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();
		dens_check = density(p);
		dens_comp = drand48() * density_max;
		
		if (dens_comp <= dens_check)
		{n_accept++;

		}
	}

	// Increase the number of points based on an approximate accept/reject ratio
	// estiamte the accept/reject ratio
	ratio = (double) n / (double) n_accept;
	cout << "Acceptance ratio estimate " << ratio << endl;
	n *= (int) ratio * 1.2;
	m = n/2;

	// Sample Fibonacci Grid points accordingly
	for (i = 0; i < n; i++){
		j = i - m;

		lambda = i / g_ratio;
		sinphi = (2.0*j) / (n + 1);
		phi = asin(sinphi);
		p.x = cos(phi) * sin(lambda);
		p.y = cos(phi) * cos(lambda);
		p.z = sin(phi);
		//p = pntFromLatLon(phi, lambda);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();
		dens_check = density(p);
		dens_comp = drand48() * density_max;

		if (dens_comp <= dens_check)
			points.push_back(p);
		
	}

	cout << "Created " << points.size() << " points using Fibonacci Grid with rejection." << endl;
}/*}}}*/

void makeFibonacciGridPoints(int n, vector<pnt>& points){/*{{{*/
	const double g_ratio = (1.0+sqrt(5)) / 2.0;
	double lambda, phi, sinphi, x, y, z;
	int i, j, m;
	pnt p;

	m = n/2;

	for (i = 0; i < n; i++){
		j = i - m;

		lambda = i / g_ratio;
		sinphi = (2.0*j) / (n + 1);
		phi = asin(sinphi);
		x = cos(phi) * sin(lambda);
		y = cos(phi) * cos(lambda);
		z = sin(phi);
		p = pntFromLatLon(phi, lambda);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();
		points.push_back(p);
	}
	cout << "Created " << points.size() << " points using Fibonacci Grid." << endl;
}/*}}}*/
/*}}}*/


#endif /* INITIALIZE_ROUTINES_H_ */
