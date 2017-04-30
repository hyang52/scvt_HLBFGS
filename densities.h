/*
 * density_routines.h
 *
 *  Created on: Jan 13, 2017
 *      Author: Huanhuan Yang	huan2yang@outlook.com
 *
 *  This file inherits and reshapes some data structures from 
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff
 */

#ifndef DENSITY_ROUTINES_H_
#define DENSITY_ROUTINES_H_

#include <stdlib.h>     
#include "triangulation.h"

double ellipse_density(const pnt &p, double lat_c, double lon_c, double lat_width, double lon_width){/*{{{*/
	//density returns the value of the density function at point p
	//	return 1.0; // Uniform density
	//	/* Ellipse Density function
	pnt work;
	double r1, r2, r;
	double dtr;
	double width, trans_cent, min_val, norm, density;

	dtr = M_PI/180.0;

	work = pntFromLatLon(p.getLat(), lon_c*dtr);
	r1 = work.dotForAngle(p);

	work = pntFromLatLon(lat_c*dtr, p.getLon());
	r2 = work.dotForAngle(p);

	r1 = r1/lon_width;
	r2 = r2/lat_width;
	r = sqrt (r1*r1 + r2*r2);

	width = 0.15;
	trans_cent = 30.0 * dtr;
	min_val = 1.0/12.0;
	min_val = pow(min_val,4);
	norm = 1.0/(1.0-min_val);
	density = ((tanh((trans_cent-r)*(1.0/width))+1.0)/2)/norm + min_val;

	return density;
	// */

}/*}}}*/
double pop_lowres_density(const pnt &p){/*{{{*/
    double dtr;
    double density, gamma, lat, lat_cent, width;

	dtr = M_PI/180.0;

    lat_cent = 25.0 * dtr;
    width = 6.0 * dtr;
    gamma = powf(1.0 / 1.95, 4.0);
    lat = p.getLat();

    density = ((1.0-gamma) * (1.0 + tanh( (lat_cent - fabs(lat)) / width)) / 2.0) + gamma;

    return density;
}/*}}}*/
double pop_highres_density(const pnt &p){/*{{{*/
    double dtr;
    double density, lat, gamma;

    dtr = M_PI/180.0;

    gamma = pow(1.0 / 3.0, 4);
    lat = p.getLat();

    density = (1-gamma) * pow( sin(lat), 4) + gamma;

    return density;
}/*}}}*/
/*}}}*/

double density(const pnt &p){/*{{{*/
	//density returns the value of the density function at point p
	//return 1.0; // Uniform density

	 //Density function for Shallow Water Test Case 5
	pnt cent;
	double r;
	double norm;
	double density;
	double min_val;
	double width, trans_cent;

	cent = pnt(0.0, -0.866025403784439, 0.5);
	cent.normalize();

	width = 0.15;
	min_val = 1.0/8.0;
	min_val = pow(min_val,4);
	trans_cent = 30.0*M_PI/180.0;
	norm = 1.0/(1.0-min_val);

	r = cent.dotForAngle(p);

	density = ((tanh((trans_cent-r)*(1.0/width))+1.0)/2.0)/norm + min_val;

	return density;
	// */

	/* Ellipse density function.

	return ellipse_density(p, 40.0, 0.0, 1.0, 0.5);
	// */

    /* Pop low resolution density function.
    return pop_lowres_density(p);
    // */

    // /* Pop high resolution density function.
    //return pop_highres_density(p); 
    // */

}/*}}}*/

#endif /* DENSITY_ROUTINES_H_ */
