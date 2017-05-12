/*
 * loop_routines.h
 *
 *  Created on: Jan 13, 2017
 *      Author: Huanhuan Yang	huan2yang@outlook.com
 *
 *  This file play the key role in scvt computation: Lloyd step, energy function, gradient.  
 *  It contains some data structures from 
 *  the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff.
 *  In particular, the routine for integration is larged changed.
 */

#ifndef LOOP_ROUTINES_H_
#define LOOP_ROUTINES_H_

#include <tr1/unordered_set>
#include "Triangle/triangle.h"
#include "setup_routines.h"
#include "densities.h"

// Uses boost mpi to make use of serialization routines for packing and unpacking of classes
namespace mpi = boost::mpi;
typedef boost::optional<mpi::status> optional;

class Quadrature
{
public:
    //! Constructor
	Quadrature(const int quadRule=2);
    //! Virtual Destructor
    virtual ~Quadrature () {}

	void integTopBot(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const{/*{{{*/
		//Call the appropriate quadrature rule that has been requested
		switch(quad_rule){
			case 0:
				integTopBot_CR(A,B,C,top,bot);
				break;
			case 1:
				integTopBot_VR(A,B,C,top,bot);
				break;
			case 2:
				integTopBot_MP(A,B,C,top,bot);
				break;
			case 3:
				integTopBot_7P(A,B,C,top,bot);
				break;
			case 4:
				integTopBot_13P(A,B,C,top,bot);
				break;
			case 5:
				integTopBot_19P(A,B,C,top,bot);
				break;
			default:
				integTopBot_MP(A,B,C,top,bot);
				break;
		}
	}/*}}}*/

	void integEnerg(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const{/*{{{*/
		//Call the appropriate quadrature rule that has been requested
		switch(quad_rule){
			case 0:
				integEnerg_CR(c,A,B,C,energ);
				break;
			case 1:
				integEnerg_VR(c,A,B,C,energ);
				break;
			case 2:
				integEnerg_MP(c,A,B,C,energ);
				break;
			case 3:
				integEnerg_7P(c,A,B,C,energ);
				break;
			case 4:
				integEnerg_13P(c,A,B,C,energ);
				break;
			case 5:
				integEnerg_19P(c,A,B,C,energ);
				break;
			default:
				integEnerg_MP(c,A,B,C,energ);
				break;
		}
	}/*}}}*/

	void integrate(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const{/*{{{*/
		//Call the appropriate quadrature rule that has been requested
		switch(quad_rule){
			case 0:
				integrate_CR(c,A,B,C,energ,top,bot);
				break;
			case 1:
				integrate_VR(c,A,B,C,energ,top,bot);
				break;
			case 2:
				integrate_MP(c,A,B,C,energ,top,bot);
				break;
			case 3:
				integrate_7P(c,A,B,C,energ,top,bot);
				break;
			case 4:
				integrate_13P(c,A,B,C,energ,top,bot);
				break;
			case 5:
				integrate_19P(c,A,B,C,energ,top,bot);
				break;
			default:
				integrate_MP(c,A,B,C,energ,top,bot);
				break;
		}
	}/*}}}*/

	void setQuadRule(const int quadRule);
private:
	void integTopBot_CR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const;
	void integTopBot_VR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const;
	void integTopBot_MP(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const;
	void integTopBot_7P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const;
	void integTopBot_13P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const;
	void integTopBot_19P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const;

	void integEnerg_CR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const;
	void integEnerg_VR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const;
	void integEnerg_MP(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const;
	void integEnerg_7P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const;
	void integEnerg_13P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const;
	void integEnerg_19P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const;

	void integrate_CR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const;
	void integrate_VR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const;
	void integrate_MP(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const;
	void integrate_7P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const;
	void integrate_13P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const;
	void integrate_19P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const;

	int quad_rule;
	double wq7[7], q7[4];
	double wq13[13], q13[8];
	double wq19[19], q19[12];

};


Quadrature::Quadrature(const int quadRule) {
	// This constructor initializes the quadrature weights (wq) and combination factors (q) for different quadrature rules on triangles.
	setQuadRule(quadRule);

}



void Quadrature::setQuadRule(const int quadRule) {
	quad_rule = quadRule;
	switch(quad_rule){
		case 0:
		case 1:
		case 2:
			break;
		case 3:
			q7[0] = (6.0 - sqrt(15.0))/21.0;
			q7[1] = (9.0 + 2.0 * sqrt(15.0))/21.0;
			q7[2] = (6.0 + sqrt(15.0))/21.0;
			q7[3] = (9.0 - 2.0 * sqrt(15.0))/21.0;

			wq7[0] = (155.0 - sqrt(15.0))/1200.0;
			wq7[1] = (155.0 - sqrt(15.0))/1200.0;
			wq7[2] = (155.0 - sqrt(15.0))/1200.0;
			wq7[3] = (155.0 + sqrt(15.0))/1200.0;
			wq7[4] = (155.0 + sqrt(15.0))/1200.0;
			wq7[5] = (155.0 + sqrt(15.0))/1200.0;
			wq7[6] = 9.0/40.0;

			break;
		case 4:
			q13[0] = 1.0/3.0;
			q13[1] = 0.479308067841923;
			q13[2] = 0.260345966079038;
			q13[3] = 0.869739794195568;
			q13[4] = 0.065130102902216;
			q13[5] = 0.638444188569809;
			q13[6] = 0.312865496004875;
			q13[7] = 0.048690315425316;

			wq13[0] = -0.149570044467670;
			wq13[1] = 0.175615257433204;
			wq13[2] = 0.175615257433204;
			wq13[3] = 0.175615257433204;
			wq13[4] = 0.053347235608839;
			wq13[5] = 0.053347235608839;
			wq13[6] = 0.053347235608839;
			wq13[7] = 0.077113760890257;
			wq13[8] = 0.077113760890257;
			wq13[9] = 0.077113760890257;
			wq13[10] = 0.077113760890257;
			wq13[11] = 0.077113760890257;
			wq13[12] = 0.077113760890257;

			break;
		case 5:
			q19[0] = 1.0 / 3.0;
			q19[1] = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
			q19[2] = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
			q19[3] = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
			q19[4] = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
			q19[5] = ( 40.0 - 10.0 * sqrt ( 15.0 )+ 10.0 * sqrt ( 7.0 ) + 2.0 * sqrt ( 105.0 ) ) / 90.0;
			q19[6] = ( 25.0 +  5.0 * sqrt ( 15.0 )-  5.0 * sqrt ( 7.0 ) - sqrt ( 105.0 ) ) / 90.0;
			q19[7] = ( 40.0 + 10.0 * sqrt ( 15.0 )+ 10.0 * sqrt ( 7.0 ) - 2.0 * sqrt ( 105.0 ) ) / 90.0;
			q19[8] = ( 25.0 -  5.0 * sqrt ( 15.0 )-  5.0 * sqrt ( 7.0 ) + sqrt ( 105.0 ) ) / 90.0;
			q19[9] = ( 40.0 + 10.0 * sqrt ( 7.0 ) ) / 90.0;
			q19[10] = ( 25.0 +  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 )- sqrt ( 105.0 ) ) / 90.0;
			q19[11] = ( 25.0 -  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 )+ sqrt ( 105.0 ) ) / 90.0;

			wq19[0] = ( 7137.0 - 1800.0 * sqrt ( 7.0 ) ) / 62720.0;
			wq19[1] = (-9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 )
				/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
				+ 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
			wq19[2] = (-9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 )
				/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
				+ 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
			wq19[3] = (-9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 )
				/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
				+ 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
			wq19[4] = (-9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 )
				/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
				- 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
			wq19[5] = (-9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 )
				/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
				- 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
			wq19[6] = (-9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 )
				/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
				- 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
			wq19[7] = (( 102791225.0 - 23876225.0 * sqrt ( 15.0 )
					- 34500875.0 * sqrt ( 7.0 )	+ 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0)
					/3.0;
			wq19[8] = (( 102791225.0 - 23876225.0 * sqrt ( 15.0 )
					- 34500875.0 * sqrt ( 7.0 )	+ 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0)
					/3.0;
			wq19[9] = (( 102791225.0 - 23876225.0 * sqrt ( 15.0 )
					- 34500875.0 * sqrt ( 7.0 )	+ 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0)
					/3.0;
			wq19[10] = (( 102791225.0 + 23876225.0 * sqrt ( 15.0 )
					- 34500875.0 * sqrt ( 7.0 )	- 9914825 * sqrt ( 105.0 ) ) / 59157504.0)/3.0;
			wq19[11] = (( 102791225.0 + 23876225.0 * sqrt ( 15.0 )
					- 34500875.0 * sqrt ( 7.0 )	- 9914825 * sqrt ( 105.0 ) ) / 59157504.0)/3.0;
			wq19[12] = (( 102791225.0 + 23876225.0 * sqrt ( 15.0 )
					- 34500875.0 * sqrt ( 7.0 )	- 9914825 * sqrt ( 105.0 ) ) / 59157504.0)/3.0;
			wq19[13] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
			wq19[14] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
			wq19[15] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
			wq19[16] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
			wq19[17] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
			wq19[18] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;

			break;
		default:
			break;
	}
}



/* ***** Integration Routines *****  {{{*/
void Quadrature::integTopBot_CR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const{/*{{{*/
	//Integrate a triangle using the Circumcenter rule
	pnt cent;
	double d;

	circumcenter(A,B,C,cent);
	cent.normalize();
	d = density(cent);

	bot = d;
	top = d*cent;
}/*}}}*/


void Quadrature::integTopBot_VR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const{/*{{{*/
	//Integrate a triangle using the vertex rule
	double d1, d2, d3;

	d1 = density(A)/3.0;
	d2 = density(B)/3.0;
	d3 = density(C)/3.0;

	bot = d1+d2+d3;
	top = d1*A + d2*B + d3*C;
}/*}}}*/


void Quadrature::integTopBot_MP(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points p1, p2, and p3.
	 *
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: Cent and dens
	 *		  Cent: Cent contains the three components of the first integral
	 *        dens: Dens contains the value of the bottom integral.
	 *****************************************************************************/
	pnt p1, p2, p3;
	double da, db, dc;
	double d1, d2, d3;
	int i, j; // Loop Index

	top = pnt(0.0,0.0,0.0,0,0);
	bot = 0.0;

	p1 = (A+B)/2.0;
	p2 = (B+C)/2.0;
	p3 = (C+A)/2.0;

	p1.normalize();
	p2.normalize();
	p3.normalize();

	da = density(A);
	db = density(B);
	dc = density(C);

	d1 = (da+db)/2.0;
	d2 = (db+dc)/2.0;
	d3 = (dc+da)/2.0;

	bot = (d1+d2+d3)/3.0;;
	top = (p1*d1 + p2*d2 + p3*d3)/3.0;
}/*}}}*/


void Quadrature::integTopBot_7P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points p1, p2, and p3.
	 *
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: Cent and dens
	 *		  Cent: Cent contains the three components of the first integral
	 *        dens: Dens contains the value of the bottom integral.
	 *****************************************************************************/
	double d;
	pnt z[7];
	int i; // Loop Index
	const double *wq, *q;
	wq = wq7;
	q = q7;

	top = pnt(0.0,0.0,0.0,0,0);
	bot = 0.0;

	z[0] = A*q[0] + B*q[0] + C*q[1];
	z[1] = A*q[0] + B*q[1] + C*q[0];
	z[2] = A*q[1] + B*q[0] + C*q[0];
	z[3] = A*q[2] + B*q[2] + C*q[3];
	z[4] = A*q[2] + B*q[3] + C*q[2];
	z[5] = A*q[3] + B*q[2] + C*q[2];
	z[6] = (A + B + C)/3.0;

	for(i = 0; i < 7; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		top += z[i]*d;
		bot += d;
	}

}/*}}}*/


void Quadrature::integTopBot_13P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points A, B, and C.
	 *
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: top and bot
	 *		  top: top contains the three components of the first integral
	 *        bot: bot contains the value of the bottom integral.
	 *****************************************************************************/
	double d;
	pnt z[13];
	int i;
	const double *wq, *q;
	wq = wq13;
	q = q13;

	top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[1]*C;
	z[2] = q[2]*A+q[1]*B+q[1]*C;
	z[3] = q[1]*A+q[1]*B+q[2]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[7]*C;
	z[8] = q[5]*A+q[7]*B+q[6]*C;
	z[9] = q[6]*A+q[5]*B+q[7]*C;
	z[10] = q[6]*A+q[7]*B+q[5]*C;
	z[11] = q[7]*A+q[5]*B+q[6]*C;
	z[12] = q[7]*A+q[6]*B+q[5]*C;

	for(i = 0; i < 13; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		top += z[i]*d;
		bot += d;
	}

}/*}}}*/


void Quadrature::integTopBot_19P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points A, B, and C.
	 *
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: top and bot
	 *		  top: top contains the three components of the first integral
	 *        bot: bot contains the value of the bottom integral.
	 *****************************************************************************/
	double d;
	pnt z[19];
	int i;
	const double *wq, *q;
	wq = wq19;
	q = q19;

	top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[2]*C;
	z[2] = q[2]*A+q[1]*B+q[2]*C;
	z[3] = q[2]*A+q[2]*B+q[1]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[6]*C;
	z[8] = q[6]*A+q[5]*B+q[6]*C;
	z[9] = q[6]*A+q[6]*B+q[5]*C;
	z[10] = q[7]*A+q[8]*B+q[8]*C;
	z[11] = q[8]*A+q[7]*B+q[8]*C;
	z[12] = q[8]*A+q[8]*B+q[7]*C;
	z[13] = q[9]*A+q[10]*B+q[11]*C;
	z[14] = q[9]*A+q[11]*B+q[10]*C;
	z[15] = q[10]*A+q[9]*B+q[11]*C;
	z[16] = q[10]*A+q[11]*B+q[9]*C;
	z[17] = q[11]*A+q[9]*B+q[10]*C;
	z[18] = q[11]*A+q[10]*B+q[9]*C;

	for(i = 0; i < 19; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		top += z[i]*d;
		bot += d;
	}

}/*}}}*/






//Revised version by Huanhuan, 	using area = 4.0*atan(tanqe) in triArea()
void divideIntegrate(const int levs, const Quadrature& quadr, const pnt &A, const pnt &B, const pnt &C, pnt &Top, double &bot){/*{{{*/
	//divideIntegrate integrates a triangle, with an arbitrary number of subdivisions.
	//Each subdivision produces 4 triangles from 1 triangle.
	//Master triangle vertices are A, B, C.
	//Number of divisions is levs
	//Top and bot are returned, being the top and bottom of the centroid calculation.

	pnt temp;
	double bot_temp;
	double area;

	Top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	if(levs==0)
	{
		quadr.integTopBot(A,B,C,temp,bot_temp);
		area = triArea(A,B,C);
		Top += temp*area;
		bot += bot_temp*area;
	}else{
		pnt z1, z2, z3;

		z1 = A;
		z2 = (A+B)/2.0;
		z3 = (A+C)/2.0;
		z2.normalize();
		z3.normalize();
		divideIntegrate(levs-1, quadr, z1, z2, z3, temp, bot_temp);
		Top += temp;
		bot += bot_temp;

		z1 = B;
		z3 = (B+C)/2.0;
		z3.normalize();
		divideIntegrate(levs-1, quadr, z1, z3, z2, temp, bot_temp);
		Top += temp;
		bot += bot_temp;

		z1 = C;
		z2 = (A+C)/2.0;
		z2.normalize();
		divideIntegrate(levs-1, quadr, z1, z2, z3, temp, bot_temp);
		Top += temp;
		bot += bot_temp;

		z1 = (A+B)/2.0;
		z1.normalize();
		divideIntegrate(levs-1, quadr, z1, z3, z2, temp, bot_temp);
		Top += temp;
		bot += bot_temp;
	}
}




void Quadrature::integEnerg_CR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const{/*{{{*/
	//Integrate a triangle using the Circumcenter rule
	//The first pnt c is the center of corresponding Voronoi cell
	pnt cent;
	double d;

	circumcenter(A,B,C,cent);
	cent.normalize();
	d = density(cent);

	energ = (cent-c).magnitude2()*d;

}/*}}}*/


void Quadrature::integEnerg_VR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const{/*{{{*/
	//Integrate a triangle using the vertex rule
	double d1, d2, d3;

	d1 = density(A)/3.0;
	d2 = density(B)/3.0;
	d3 = density(C)/3.0;

	energ = (A-c).magnitude2()*d1 + (B-c).magnitude2()*d2 + (C-c).magnitude2()*d3;

}/*}}}*/


void Quadrature::integEnerg_MP(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points p1, p2, and p3.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	pnt p1, p2, p3;
	double da, db, dc;
	double d1, d2, d3;
	int i, j; // Loop Index

	p1 = (A+B)/2.0;
	p2 = (B+C)/2.0;
	p3 = (C+A)/2.0;

	p1.normalize();
	p2.normalize();
	p3.normalize();

	da = density(A);
	db = density(B);
	dc = density(C);

	d1 = (da+db)/2.0;
	d2 = (db+dc)/2.0;
	d3 = (dc+da)/2.0;

	energ = ((p1-c).magnitude2()*d1 + (p2-c).magnitude2()*d2 + (p3-c).magnitude2()*d3)/3.0;

}/*}}}*/


void Quadrature::integEnerg_7P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points p1, p2, and p3.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	double d;
	pnt z[7];
	int i; // Loop Index
	const double *wq, *q;
	wq = wq7;
	q = q7;

	energ = 0.0;

	z[0] = A*q[0] + B*q[0] + C*q[1];
	z[1] = A*q[0] + B*q[1] + C*q[0];
	z[2] = A*q[1] + B*q[0] + C*q[0];
	z[3] = A*q[2] + B*q[2] + C*q[3];
	z[4] = A*q[2] + B*q[3] + C*q[2];
	z[5] = A*q[3] + B*q[2] + C*q[2];
	z[6] = (A + B + C)/3.0;

	for(i = 0; i < 7; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		energ += (z[i]-c).magnitude2()*d;
	}

}/*}}}*/


void Quadrature::integEnerg_13P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points A, B, and C.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	double d;
	pnt z[13];
	int i;
	const double *wq, *q;
	wq = wq13;
	q = q13;

	energ = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[1]*C;
	z[2] = q[2]*A+q[1]*B+q[1]*C;
	z[3] = q[1]*A+q[1]*B+q[2]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[7]*C;
	z[8] = q[5]*A+q[7]*B+q[6]*C;
	z[9] = q[6]*A+q[5]*B+q[7]*C;
	z[10] = q[6]*A+q[7]*B+q[5]*C;
	z[11] = q[7]*A+q[5]*B+q[6]*C;
	z[12] = q[7]*A+q[6]*B+q[5]*C;

	for(i = 0; i < 13; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		energ += (z[i]-c).magnitude2()*d;
	}

}/*}}}*/


void Quadrature::integEnerg_19P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points A, B, and C.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	double d;
	pnt z[19];
	int i;
	const double *wq, *q;
	wq = wq19;
	q = q19;

	energ = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[2]*C;
	z[2] = q[2]*A+q[1]*B+q[2]*C;
	z[3] = q[2]*A+q[2]*B+q[1]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[6]*C;
	z[8] = q[6]*A+q[5]*B+q[6]*C;
	z[9] = q[6]*A+q[6]*B+q[5]*C;
	z[10] = q[7]*A+q[8]*B+q[8]*C;
	z[11] = q[8]*A+q[7]*B+q[8]*C;
	z[12] = q[8]*A+q[8]*B+q[7]*C;
	z[13] = q[9]*A+q[10]*B+q[11]*C;
	z[14] = q[9]*A+q[11]*B+q[10]*C;
	z[15] = q[10]*A+q[9]*B+q[11]*C;
	z[16] = q[10]*A+q[11]*B+q[9]*C;
	z[17] = q[11]*A+q[9]*B+q[10]*C;
	z[18] = q[11]*A+q[10]*B+q[9]*C;

	for(i = 0; i < 19; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		energ += (z[i]-c).magnitude2()*d;
	}

}/*}}}*/






//Revised version by Huanhuan, 	using area = 4.0*atan(tanqe) in triArea()
void divideIntegEnerg(const int levs, const Quadrature& quadr, const pnt& VorC, const pnt &A, const pnt &B, const pnt &C, double &energ){
	//divideIntegEnerg integrates a triangle, with an arbitrary number of subdivisions.
	//Master triangle vertices are A, B, C. Not clear about its orientation
	//Each subdivision produces 4 triangles from 1 triangle. Keep the orientation
	//The first pnt VorC is the center of corresponding Voronoi cell
	//Number of divisions is levs

	pnt temp;
	double energ_temp;
	double area;

	energ = 0.0;

	if(levs==0)
	{
		quadr.integEnerg(VorC,A,B,C,energ_temp);
		area = triArea(A,B,C);
		energ += energ_temp*area;
	}else{
		pnt z1, z2, z3;

		z1 = A;
		z2 = (A+B)/2.0;
		z3 = (A+C)/2.0;
		z2.normalize();
		z3.normalize();
		divideIntegEnerg(levs-1, quadr, VorC, z1, z2, z3, energ_temp);
		energ += energ_temp;

		z1 = B;
		z3 = (B+C)/2.0;
		z3.normalize();
		divideIntegEnerg(levs-1, quadr, VorC, z1, z3, z2, energ_temp);
		energ += energ_temp;

		z1 = C;
		z2 = (A+C)/2.0;
		z2.normalize();
		divideIntegEnerg(levs-1, quadr, VorC, z1, z2, z3, energ_temp);
		energ += energ_temp;

		z1 = (A+B)/2.0;
		z1.normalize();
		divideIntegEnerg(levs-1, quadr, VorC, z1, z3, z2, energ_temp);
		energ += energ_temp;
	}
}








void Quadrature::integrate_CR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const{/*{{{*/
	//Integrate a triangle using the Circumcenter rule
	//The first pnt c is the center of corresponding Voronoi cell
	pnt cent;
	double d;

	circumcenter(A,B,C,cent);
	cent.normalize();
	d = density(cent);

	energ = (cent-c).magnitude2()*d;
	bot = d;
	top = d*cent;
}/*}}}*/


void Quadrature::integrate_VR(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const{/*{{{*/
	//Integrate a triangle using the vertex rule
	double d1, d2, d3;

	d1 = density(A)/3.0;
	d2 = density(B)/3.0;
	d3 = density(C)/3.0;

	energ = (A-c).magnitude2()*d1 + (B-c).magnitude2()*d2 + (C-c).magnitude2()*d3;
	bot = d1+d2+d3;
	top = d1*A + d2*B + d3*C;
}/*}}}*/


void Quadrature::integrate_MP(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points p1, p2, and p3.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	pnt p1, p2, p3;
	double da, db, dc;
	double d1, d2, d3;
	int i, j; // Loop Index

	p1 = (A+B)/2.0;
	p2 = (B+C)/2.0;
	p3 = (C+A)/2.0;

	p1.normalize();
	p2.normalize();
	p3.normalize();

	da = density(A);
	db = density(B);
	dc = density(C);

	d1 = (da+db)/2.0;
	d2 = (db+dc)/2.0;
	d3 = (dc+da)/2.0;

	energ = ((p1-c).magnitude2()*d1 + (p2-c).magnitude2()*d2 + (p3-c).magnitude2()*d3)/3.0;
	bot = (d1+d2+d3)/3.0;;
	top = (p1*d1 + p2*d2 + p3*d3)/3.0;
}/*}}}*/


void Quadrature::integrate_7P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points p1, p2, and p3.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	double d;
	pnt z[7];
	int i; // Loop Index
	const double *wq, *q;
	wq = wq7;
	q = q7;

	energ = 0.0;
	top = pnt(0.0,0.0,0.0,0,0);
	bot = 0.0;

	z[0] = A*q[0] + B*q[0] + C*q[1];
	z[1] = A*q[0] + B*q[1] + C*q[0];
	z[2] = A*q[1] + B*q[0] + C*q[0];
	z[3] = A*q[2] + B*q[2] + C*q[3];
	z[4] = A*q[2] + B*q[3] + C*q[2];
	z[5] = A*q[3] + B*q[2] + C*q[2];
	z[6] = (A + B + C)/3.0;

	for(i = 0; i < 7; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		energ += (z[i]-c).magnitude2()*d;
		top += z[i]*d;
		bot += d;
	}

}/*}}}*/


void Quadrature::integrate_13P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points A, B, and C.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	double d;
	pnt z[13];
	int i;
	const double *wq, *q;
	wq = wq13;
	q = q13;

	energ = 0.0;
	top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[1]*C;
	z[2] = q[2]*A+q[1]*B+q[1]*C;
	z[3] = q[1]*A+q[1]*B+q[2]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[7]*C;
	z[8] = q[5]*A+q[7]*B+q[6]*C;
	z[9] = q[6]*A+q[5]*B+q[7]*C;
	z[10] = q[6]*A+q[7]*B+q[5]*C;
	z[11] = q[7]*A+q[5]*B+q[6]*C;
	z[12] = q[7]*A+q[6]*B+q[5]*C;

	for(i = 0; i < 13; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		energ += (z[i]-c).magnitude2()*d;
		top += z[i]*d;
		bot += d;
	}

}/*}}}*/


void Quadrature::integrate_19P(const pnt &c, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot) const{/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle
	 *      v defined by the points A, B, and C.
	 *	  - The first pnt c is the center of corresponding Voronoi cell
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: energ
	 *****************************************************************************/
	double d;
	pnt z[19];
	int i;
	const double *wq, *q;
	wq = wq19;
	q = q19;

	energ = 0.0;
	top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[2]*C;
	z[2] = q[2]*A+q[1]*B+q[2]*C;
	z[3] = q[2]*A+q[2]*B+q[1]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[6]*C;
	z[8] = q[6]*A+q[5]*B+q[6]*C;
	z[9] = q[6]*A+q[6]*B+q[5]*C;
	z[10] = q[7]*A+q[8]*B+q[8]*C;
	z[11] = q[8]*A+q[7]*B+q[8]*C;
	z[12] = q[8]*A+q[8]*B+q[7]*C;
	z[13] = q[9]*A+q[10]*B+q[11]*C;
	z[14] = q[9]*A+q[11]*B+q[10]*C;
	z[15] = q[10]*A+q[9]*B+q[11]*C;
	z[16] = q[10]*A+q[11]*B+q[9]*C;
	z[17] = q[11]*A+q[9]*B+q[10]*C;
	z[18] = q[11]*A+q[10]*B+q[9]*C;

	for(i = 0; i < 19; i ++){
		z[i].normalize();
		d = density(z[i])*wq[i];

		energ += (z[i]-c).magnitude2()*d;
		top += z[i]*d;
		bot += d;
	}

}/*}}}*/







//Revised version by Huanhuan, 	using area = 4.0*atan(tanqe) in triArea()
void divideIntegrate(const int levs, const Quadrature& quadr, const pnt& VorC, const pnt &A, const pnt &B, const pnt &C, double &energ, pnt &top, double &bot){
	//divideIntegEnerg integrates a triangle, with an arbitrary number of subdivisions.
	//Master triangle vertices are A, B, C. Not clear about its orientation
	//Each subdivision produces 4 triangles from 1 triangle. Keep the orientation
	//The first pnt VorC is the center of corresponding Voronoi cell
	//Number of divisions is levs

	pnt temp;
	double energ_temp, bot_temp;
	double area;

	energ = 0.0;
	top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	if(levs==0)
	{
		quadr.integrate(VorC,A,B,C,energ_temp,temp,bot_temp);
		area = triArea(A,B,C);
		energ += energ_temp*area;
		top += temp*area;
		bot += bot_temp*area;
	}else{
		pnt z1, z2, z3;

		z1 = A;
		z2 = (A+B)/2.0;
		z3 = (A+C)/2.0;
		z2.normalize();
		z3.normalize();
		divideIntegrate(levs-1, quadr, VorC, z1, z2, z3, energ_temp, temp, bot_temp);
		energ += energ_temp;
		top += temp;
		bot += bot_temp;

		z1 = B;
		z3 = (B+C)/2.0;
		z3.normalize();
		divideIntegrate(levs-1, quadr, VorC, z1, z3, z2, energ_temp, temp, bot_temp);
		energ += energ_temp;
		top += temp;
		bot += bot_temp;

		z1 = C;
		z2 = (A+C)/2.0;
		z2.normalize();
		divideIntegrate(levs-1, quadr, VorC, z1, z2, z3, energ_temp, temp, bot_temp);
		energ += energ_temp;
		top += temp;
		bot += bot_temp;

		z1 = (A+B)/2.0;
		z1.normalize();
		divideIntegrate(levs-1, quadr, VorC, z1, z3, z2, energ_temp, temp, bot_temp);
		energ += energ_temp;
		top += temp;
		bot += bot_temp;
	}
}








/* ***** Generic Region Routines ***** {{{ */


void triangulateRegions(const int id, char* flags, vector<region> &region_vec){/*{{{*/
	//Triangulate my region(s) points
	//Points are first stereographically projected into a plane tanget to region center
	//Projected point set is then triangulated using Triangle
	//Indexing is mapped from local (projected) point set to global (spherical) point set
	//Triangles whose circumcircles extend outside of the region radius are considered "bad" and are not kept
	struct triangulateio in, out, vorout;
	pnt ccenter;
	pnt a, b, c;
	pnt ac, bc;
	pnt x_hat, y_hat, axis;
	pnt Q;
	double cradius;
	double scale;
	double x, y, z;
	double criteria;
	double s;
	double min_dir;
	tri t;
	int i;
	int vi1, vi2, vi3;
	int swp_vi;
	vector<region>::iterator region_itr;
	vector<pnt>::iterator point_itr;

#ifdef _DEBUG
	cerr << "Triangulating points (local) " << id << endl;
#endif

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		in.numberofpoints = (*region_itr).points.size();
		in.numberofpointattributes = 0;
		in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
		in.numberofsegments = 0;
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.regionlist = (double *) NULL;
		in.pointmarkerlist = (int *) NULL;

		//cerr << "In Triangulation (local)" << id << ", num pts ="<< in.numberofpoints << endl;

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

			//Counterclockwise tri on plan maybe clockwise on sphere. 
			//If tri is CW, no matter which semisphere point a belongs to, isCCw is false
			if(isCcw(a,b,c)==-1){								
				b = (*region_itr).points.at(vi3);
				c = (*region_itr).points.at(vi2);
			}

			vi1 = a.idx;
			vi2 = b.idx;
			vi3 = c.idx;

			t = tri(vi1, vi2, vi3);
			(*region_itr).triangles.push_back(t);
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


void integrateRegions(const int id, const int div_levs, const Quadrature& quadr, const int use_barycenter, vector<region>& regions, vector<region> &region_vec, vector<pnt>& points, vector<pnt>& n_points){/*{{{*/
	// Integrate Voronoi cells inside of my region
	// Every region updates all points that are closer to their region center than any other region center.
	// This ensures that each point is only updated once.
	pnt *tops;
	double *bots;
	int *sides;
	int vi1, vi2, vi3;
	pnt a, b, c;
	pnt ab, bc, ca;
	pnt ccenter;
	pnt np;
	double dist_temp;
	double a_dist_to_region, a_min_dist;
	double b_dist_to_region, b_min_dist;
	double c_dist_to_region, c_min_dist;
	int a_min_region, b_min_region, c_min_region;
	pnt top_val;
	double bot_val;
	int i;
	int sign;
	vector<region>::iterator region_itr;
	vector<tri>::iterator tri_itr;
	vector<pnt>::iterator point_itr;
	vector<int>::iterator neighbor_itr;


	#ifdef _DEBUG
		cerr << "Integrating regions " << id << endl;
	#endif

	tops = new pnt[points.size()];
	bots = new double[points.size()];
	sides = new int[points.size()];

	for(i = 0; i < points.size(); i++){
		tops[i] = pnt(0.0,0.0,0.0,0,i);
		bots[i] = 0.0;
		sides[i] = 0;
	}

	i = 0;

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			vi1 = (*tri_itr).vi1;
			vi2 = (*tri_itr).vi2;
			vi3 = (*tri_itr).vi3;

			a = points[vi1];
			b = points[vi2];
			c = points[vi3];

			sides[vi1]++;
			sides[vi2]++;
			sides[vi3]++;

			if ( !a.isBdry || !b.isBdry || !c.isBdry ) {
				ab = (a+b)/2.0;
				bc = (b+c)/2.0;
				ca = (c+a)/2.0;

				ab.normalize();
				bc.normalize();
				ca.normalize();

				if(use_barycenter){
					ccenter = (a+b+c)/3.0;
				} else {
					try {
						circumcenter(a,b,c,ccenter);
					} catch ( int n ) {
						ccenter = (a + b + c) / 3.0;
					}
				}
				ccenter.normalize();

				a_dist_to_region = a.dotForAngle((*region_itr).center);
				b_dist_to_region = b.dotForAngle((*region_itr).center);
				c_dist_to_region = c.dotForAngle((*region_itr).center);

				a_min_dist = 10.0;
				b_min_dist = 10.0;
				c_min_dist = 10.0;

				for(neighbor_itr = (*region_itr).neighbors.begin(); neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
					if(regions[(*neighbor_itr)].center.idx != (*region_itr).center.idx){
						dist_temp = a.dotForAngle(regions[(*neighbor_itr)].center);
						if(a_min_dist > dist_temp){
							a_min_dist = dist_temp;
							a_min_region = (*neighbor_itr);
						}
						dist_temp = b.dotForAngle(regions[(*neighbor_itr)].center);
						if(b_min_dist > dist_temp){
							b_min_dist = dist_temp;
							b_min_region = (*neighbor_itr);
						}
						dist_temp = c.dotForAngle(regions[(*neighbor_itr)].center);
						if(c_min_dist > dist_temp){
							c_min_dist = dist_temp;
							c_min_region = (*neighbor_itr);
						}
					}
				}

				if(a_dist_to_region < a_min_dist || (a_dist_to_region == a_min_dist && (*region_itr).center.idx < a_min_region)){
					//Triangle 1 - a ab ccenter
					divideIntegrate(div_levs,quadr,a,ab,ccenter,top_val,bot_val);
					sign = isCcw(a,ab,ccenter);
					tops[a.idx] += top_val*sign;
					bots[a.idx] += bot_val*sign;

					//Triangle 2 - a ccenter ca
					divideIntegrate(div_levs,quadr,a,ccenter,ca,top_val,bot_val);
					sign = isCcw(a,ccenter,ca);
					tops[a.idx] += top_val*sign;
					bots[a.idx] += bot_val*sign;				
				}

				if(b_dist_to_region < b_min_dist || (b_dist_to_region == b_min_dist && (*region_itr).center.idx < b_min_region)){
					//Triangle 1 - b bc ccenter
					divideIntegrate(div_levs,quadr,b,bc,ccenter,top_val,bot_val);
					sign = isCcw(b,bc,ccenter);
					tops[b.idx] += top_val*sign;
					bots[b.idx] += bot_val*sign;

					//Triangle 2 - b ccenter ab
					divideIntegrate(div_levs,quadr,b,ccenter,ab,top_val,bot_val);
					sign = isCcw(b,ccenter,ab);
					tops[b.idx] += top_val*sign;
					bots[b.idx] += bot_val*sign;
				} 

				if(c_dist_to_region < c_min_dist || (c_dist_to_region == c_min_dist && (*region_itr).center.idx < c_min_region)){
					//Triangle 1 - c ca ccenter
					divideIntegrate(div_levs,quadr,c,ca,ccenter,top_val,bot_val);
					sign = isCcw(c,ca,ccenter);
					tops[c.idx] += top_val*sign;
					bots[c.idx] += bot_val*sign;

					//Triangle 2 - c ccenter bc
					divideIntegrate(div_levs,quadr,c,ccenter,bc,top_val,bot_val);
					sign = isCcw(c,ccenter,bc);
					tops[c.idx] += top_val*sign;
					bots[c.idx] += bot_val*sign;
				} 
			}
		}
	}

	n_points.clear();
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		if(!(*point_itr).isBdry){
			if(bots[(*point_itr).idx] != 0.0){
				np = tops[(*point_itr).idx]/bots[(*point_itr).idx];
				np.idx = (*point_itr).idx;
				np.sides = sides[(*point_itr).idx];
				np.isBdry = (*point_itr).isBdry;
				np.normalize();
				n_points.push_back(np);
			}
		} else if((*point_itr).isBdry){
			if((*point_itr).isBdry == 2){
				(*point_itr).isBdry = 0;
			}

			n_points.push_back((*point_itr));
		}
	}

	delete(tops);
	delete(bots);

	#ifdef _DEBUG
		cerr << "Done Integrating regions " << id << endl;
	#endif
	return;
}/*}}}*/






void inteGradient(const int id, const int div_levs, const Quadrature& quadr, const int use_barycenter, vector<region>& regions, vector<region> &my_regions, const vector<pnt>& points, vector<pnt>& distr_grad){/*{{{*/
	// Integrate Voronoi cells inside of my region
	// Every region updates all points that are closer to their region center than any other region center.
	// This ensures that each point is only updated once.
	pnt *tops;
	double *bots;
	int vi1, vi2, vi3;
	pnt a, b, c;
	pnt ab, bc, ca;
	pnt ccenter;
	pnt grad_i;
	double dist_temp;
	double a_dist_to_region, a_min_dist;
	double b_dist_to_region, b_min_dist;
	double c_dist_to_region, c_min_dist;
	int a_min_region, b_min_region, c_min_region;
	pnt top_val;
	double bot_val;
	int i;
	int sign;
	vector<region>::iterator region_itr;
	vector<tri>::iterator tri_itr;
	vector<pnt>::const_iterator point_itr;
	vector<int>::iterator neighbor_itr;


	#ifdef _DEBUG
		cerr << "Integrating regions " << id << endl;
	#endif

	tops = new pnt[points.size()];
	bots = new double[points.size()];

	for(i = 0; i < points.size(); i++){
		tops[i] = pnt(0.0,0.0,0.0,0,i);
		bots[i] = 0.0;
	}

	i = 0;

	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			vi1 = (*tri_itr).vi1;
			vi2 = (*tri_itr).vi2;
			vi3 = (*tri_itr).vi3;

			a = points[vi1];
			b = points[vi2];
			c = points[vi3];

			if ( !a.isBdry || !b.isBdry || !c.isBdry ) {
				ab = (a+b)/2.0;
				bc = (b+c)/2.0;
				ca = (c+a)/2.0;

				ab.normalize();
				bc.normalize();
				ca.normalize();

				if(use_barycenter){
					ccenter = (a+b+c)/3.0;
				} else {
					try {
						circumcenter(a,b,c,ccenter);
					} catch ( int n ) {
						ccenter = (a + b + c) / 3.0;
					}
				}
				ccenter.normalize();

				a_dist_to_region = a.dotForAngle((*region_itr).center);
				b_dist_to_region = b.dotForAngle((*region_itr).center);
				c_dist_to_region = c.dotForAngle((*region_itr).center);

				a_min_dist = 10.0;
				b_min_dist = 10.0;
				c_min_dist = 10.0;

				for(neighbor_itr = (*region_itr).neighbors.begin(); neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
					if(regions[(*neighbor_itr)].center.idx != (*region_itr).center.idx){
						dist_temp = a.dotForAngle(regions[(*neighbor_itr)].center);
						if(a_min_dist > dist_temp){
							a_min_dist = dist_temp;
							a_min_region = (*neighbor_itr);
						}
						dist_temp = b.dotForAngle(regions[(*neighbor_itr)].center);
						if(b_min_dist > dist_temp){
							b_min_dist = dist_temp;
							b_min_region = (*neighbor_itr);
						}
						dist_temp = c.dotForAngle(regions[(*neighbor_itr)].center);
						if(c_min_dist > dist_temp){
							c_min_dist = dist_temp;
							c_min_region = (*neighbor_itr);
						}
					}
				}

				if(a_dist_to_region < a_min_dist || (a_dist_to_region == a_min_dist && (*region_itr).center.idx < a_min_region)){
					//Triangle 1 - a ab ccenter
					divideIntegrate(div_levs,quadr,a,ab,ccenter,top_val,bot_val);
					sign = isCcw(a,ab,ccenter);
					tops[a.idx] += top_val*sign;
					bots[a.idx] += bot_val*sign;
					//Triangle 2 - a ccenter ca
					divideIntegrate(div_levs,quadr,a,ccenter,ca,top_val,bot_val);
					sign = isCcw(a,ccenter,ca);
					tops[a.idx] += top_val*sign;
					bots[a.idx] += bot_val*sign;			
				}

				if(b_dist_to_region < b_min_dist || (b_dist_to_region == b_min_dist && (*region_itr).center.idx < b_min_region)){
					//Triangle 1 - b bc ccenter
					divideIntegrate(div_levs,quadr,b,bc,ccenter,top_val,bot_val);
					sign = isCcw(b,bc,ccenter);
					tops[b.idx] += top_val*sign;
					bots[b.idx] += bot_val*sign;
					//Triangle 2 - b ccenter ab
					divideIntegrate(div_levs,quadr,b,ccenter,ab,top_val,bot_val);
					sign = isCcw(b,ccenter,ab);
					tops[b.idx] += top_val*sign;
					bots[b.idx] += bot_val*sign;
				} 

				if(c_dist_to_region < c_min_dist || (c_dist_to_region == c_min_dist && (*region_itr).center.idx < c_min_region)){
					//Triangle 1 - c ca ccenter
					divideIntegrate(div_levs,quadr,c,ca,ccenter,top_val,bot_val);
					sign = isCcw(c,ca,ccenter);
					tops[c.idx] += top_val*sign;
					bots[c.idx] += bot_val*sign;
					//Triangle 2 - c ccenter bc
					divideIntegrate(div_levs,quadr,c,ccenter,bc,top_val,bot_val);
					sign = isCcw(c,ccenter,bc);
					tops[c.idx] += top_val*sign;
					bots[c.idx] += bot_val*sign;
				} 
			}
		}
	}

	distr_grad.clear();
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		if(!(*point_itr).isBdry){
			if(bots[(*point_itr).idx] != 0.0){
				grad_i = 2.0*(*point_itr)*bots[(*point_itr).idx] - 2.0*tops[(*point_itr).idx];
				grad_i.idx = (*point_itr).idx;
				distr_grad.push_back(grad_i);
			}
		} else if((*point_itr).isBdry){
			/*if((*point_itr).isBdry == 2){
				(*point_itr).isBdry = 0;
			}
			pnt zero(0.0,0.0,0.0);
			distr_grad.push_back(zero); */
		}
	}

	delete(tops);
	delete(bots);

	#ifdef _DEBUG
		cerr << "Done Integrating regions " << id << endl;
	#endif
	return;
}/*}}}*/





void inteEnergy(const int id, const int div_levs, const Quadrature& quadr, const int use_barycenter, vector<region>& regions, vector<region> &my_regions, const vector<pnt>& points, double& my_energy){/*{{{*/
	// Integrate Voronoi cells inside of my region
	// Every region updates all points that are closer to their region center than any other region center.
	// This ensures that each point is only updated once.
	double *energs;
	int vi1, vi2, vi3;
	pnt a, b, c;
	pnt ab, bc, ca;
	pnt ccenter;
	pnt grad_i;
	double dist_temp;
	double a_dist_to_region, a_min_dist;
	double b_dist_to_region, b_min_dist;
	double c_dist_to_region, c_min_dist;
	int a_min_region, b_min_region, c_min_region;
	double energ_val;
	int i;
	int sign;
	vector<region>::iterator region_itr;
	vector<tri>::iterator tri_itr;
	vector<pnt>::const_iterator point_itr;
	vector<int>::iterator neighbor_itr;


	#ifdef _DEBUG
		cerr << "Integrating regions " << id << endl;
	#endif

	my_energy = 0.0;
	energs = new double[points.size()];

	for(i = 0; i < points.size(); i++){
		energs[i] = 0.0;
	}

	i = 0;

	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			vi1 = (*tri_itr).vi1;
			vi2 = (*tri_itr).vi2;
			vi3 = (*tri_itr).vi3;

			a = points[vi1];
			b = points[vi2];
			c = points[vi3];

			if ( !a.isBdry || !b.isBdry || !c.isBdry ) {
				ab = (a+b)/2.0;
				bc = (b+c)/2.0;
				ca = (c+a)/2.0;

				ab.normalize();
				bc.normalize();
				ca.normalize();

				if(use_barycenter){
					ccenter = (a+b+c)/3.0;
				} else {
					try {
						circumcenter(a,b,c,ccenter);
					} catch ( int n ) {
						ccenter = (a + b + c) / 3.0;
					}
				}
				ccenter.normalize();

				a_dist_to_region = a.dotForAngle((*region_itr).center);
				b_dist_to_region = b.dotForAngle((*region_itr).center);
				c_dist_to_region = c.dotForAngle((*region_itr).center);

				a_min_dist = 10.0;
				b_min_dist = 10.0;
				c_min_dist = 10.0;

				for(neighbor_itr = (*region_itr).neighbors.begin(); neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
					if(regions[(*neighbor_itr)].center.idx != (*region_itr).center.idx){
						dist_temp = a.dotForAngle(regions[(*neighbor_itr)].center);
						if(a_min_dist > dist_temp){
							a_min_dist = dist_temp;
							a_min_region = (*neighbor_itr);
						}
						dist_temp = b.dotForAngle(regions[(*neighbor_itr)].center);
						if(b_min_dist > dist_temp){
							b_min_dist = dist_temp;
							b_min_region = (*neighbor_itr);
						}
						dist_temp = c.dotForAngle(regions[(*neighbor_itr)].center);
						if(c_min_dist > dist_temp){
							c_min_dist = dist_temp;
							c_min_region = (*neighbor_itr);
						}
					}
				}

				if(a_dist_to_region < a_min_dist || (a_dist_to_region == a_min_dist && (*region_itr).center.idx < a_min_region)){
					//Triangle 1 - a ab ccenter
					divideIntegEnerg(div_levs,quadr,a,a,ab,ccenter,energ_val);
					sign = isCcw(a,ab,ccenter);
					energs[a.idx] += energ_val*sign;
					//Triangle 2 - a ccenter ca
					divideIntegEnerg(div_levs,quadr,a,a,ccenter,ca,energ_val);
					sign = isCcw(a,ccenter,ca);
					energs[a.idx] += energ_val*sign;
				} 

				if(b_dist_to_region < b_min_dist || (b_dist_to_region == b_min_dist && (*region_itr).center.idx < b_min_region)){
					//Triangle 1 - b bc ccenter
					divideIntegEnerg(div_levs,quadr,b,b,bc,ccenter,energ_val);
					sign = isCcw(b,bc,ccenter);
					energs[b.idx] += energ_val*sign;
					//Triangle 2 - b ccenter ab
					divideIntegEnerg(div_levs,quadr,b,b,ccenter,ab,energ_val);
					sign = isCcw(b,ccenter,ab);
					energs[b.idx] += energ_val*sign;
				}

				if(c_dist_to_region < c_min_dist || (c_dist_to_region == c_min_dist && (*region_itr).center.idx < c_min_region)){
					//Triangle 1 - c ca ccenter
					divideIntegEnerg(div_levs,quadr,c,c,ca,ccenter,energ_val);
					sign = isCcw(c,ca,ccenter);
					energs[c.idx] += energ_val*sign;
					//Triangle 2 - c ccenter bc
					divideIntegEnerg(div_levs,quadr,c,c,ccenter,bc,energ_val);
					sign = isCcw(c,ccenter,bc);
					energs[c.idx] += energ_val*sign;
				} 
			}
		}
	}

	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		// Ignoring boundary conditions at this moment
		my_energy += energs[(*point_itr).idx];
	}

	delete(energs);

	#ifdef _DEBUG
		cerr << "Done Integrating regions " << id << endl;
	#endif
	return;
}/*}}}*/






void inteEnergGrad(const int id, const int div_levs, const Quadrature& quadr, const int use_barycenter, vector<region>& regions, vector<region> &my_regions, const vector<pnt>& points, double& my_energy, vector<pnt>& distr_grad, vector<pnt>& distr_nPoints, double* my_bots){/*{{{*/
	// Integrate Voronoi cells inside of my region
	// Every region updates all points that are closer to their region center than any other region center.
	// This ensures that each point is only updated once.
	double *energs;
	pnt *tops;
	double *bots;
	int vi1, vi2, vi3;
	pnt a, b, c;
	pnt ab, bc, ca;
	pnt ccenter;
	pnt grad_i, pnt_i;
	double dist_temp;
	double a_dist_to_region, a_min_dist;
	double b_dist_to_region, b_min_dist;
	double c_dist_to_region, c_min_dist;
	int a_min_region, b_min_region, c_min_region;
	double energ_val;
	pnt top_val;
	double bot_val;
	int sign;
	vector<region>::iterator region_itr;
	vector<tri>::iterator tri_itr;
	vector<pnt>::const_iterator point_itr;
	vector<int>::iterator neighbor_itr;

	#ifdef _DEBUG
		cerr << "Integrating regions " << id << endl;
	#endif

	my_energy = 0.0;
	energs = new double[points.size()];
	tops = new pnt[points.size()];
	bots = new double[points.size()];
	for(int i = 0; i < points.size(); i++){
		energs[i] = 0.0;
		tops[i] = pnt(0.0,0.0,0.0,0,i);
		bots[i] = 0.0;
	}

	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			vi1 = (*tri_itr).vi1;
			vi2 = (*tri_itr).vi2;
			vi3 = (*tri_itr).vi3;

			a = points[vi1];
			b = points[vi2];
			c = points[vi3];

			if ( !a.isBdry || !b.isBdry || !c.isBdry ) {
				ab = (a+b)/2.0;
				bc = (b+c)/2.0;
				ca = (c+a)/2.0;

				ab.normalize();
				bc.normalize();
				ca.normalize();

				if(use_barycenter){
					ccenter = (a+b+c)/3.0;
				} else {
					try {
						circumcenter(a,b,c,ccenter);
					} catch ( int n ) {
						ccenter = (a + b + c) / 3.0;
					}
				}
				ccenter.normalize();

				a_dist_to_region = a.dotForAngle((*region_itr).center);
				b_dist_to_region = b.dotForAngle((*region_itr).center);
				c_dist_to_region = c.dotForAngle((*region_itr).center);

				a_min_dist = 10.0;
				b_min_dist = 10.0;
				c_min_dist = 10.0;

				for(neighbor_itr = (*region_itr).neighbors.begin(); neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
					if(regions[(*neighbor_itr)].center.idx != (*region_itr).center.idx){
						dist_temp = a.dotForAngle(regions[(*neighbor_itr)].center);
						if(a_min_dist > dist_temp){
							a_min_dist = dist_temp;
							a_min_region = (*neighbor_itr);
						}
						dist_temp = b.dotForAngle(regions[(*neighbor_itr)].center);
						if(b_min_dist > dist_temp){
							b_min_dist = dist_temp;
							b_min_region = (*neighbor_itr);
						}
						dist_temp = c.dotForAngle(regions[(*neighbor_itr)].center);
						if(c_min_dist > dist_temp){
							c_min_dist = dist_temp;
							c_min_region = (*neighbor_itr);
						}
					}
				}

				if(a_dist_to_region < a_min_dist || (a_dist_to_region == a_min_dist && (*region_itr).center.idx < a_min_region)){
					//Triangle 1 - a ab ccenter
					divideIntegrate(div_levs,quadr,a,a,ab,ccenter,energ_val,top_val,bot_val);
					sign = isCcw(a,ab,ccenter);
					energs[a.idx] += energ_val*sign;
					tops[a.idx] += top_val*sign;
					bots[a.idx] += bot_val*sign;
					//Triangle 2 - a ccenter ca
					divideIntegrate(div_levs,quadr,a,a,ccenter,ca,energ_val,top_val,bot_val);
					sign = isCcw(a,ccenter,ca);
					energs[a.idx] += energ_val*sign;
					tops[a.idx] += top_val*sign;
					bots[a.idx] += bot_val*sign;
				} 

				if(b_dist_to_region < b_min_dist || (b_dist_to_region == b_min_dist && (*region_itr).center.idx < b_min_region)){
					//Triangle 1 - b bc ccenter
					divideIntegrate(div_levs,quadr,b,b,bc,ccenter,energ_val,top_val,bot_val);
					sign = isCcw(b,bc,ccenter);
					energs[b.idx] += energ_val*sign;
					tops[b.idx] += top_val*sign;
					bots[b.idx] += bot_val*sign;
					//Triangle 2 - b ccenter ab
					divideIntegrate(div_levs,quadr,b,b,ccenter,ab,energ_val,top_val,bot_val);
					sign = isCcw(b,ccenter,ab);
					energs[b.idx] += energ_val*sign;
					tops[b.idx] += top_val*sign;
					bots[b.idx] += bot_val*sign;
				}

				if(c_dist_to_region < c_min_dist || (c_dist_to_region == c_min_dist && (*region_itr).center.idx < c_min_region)){
					//Triangle 1 - c ca ccenter
					divideIntegrate(div_levs,quadr,c,c,ca,ccenter,energ_val,top_val,bot_val);
					sign = isCcw(c,ca,ccenter);
					energs[c.idx] += energ_val*sign;
					tops[c.idx] += top_val*sign;
					bots[c.idx] += bot_val*sign;
					//Triangle 2 - c ccenter bc
					divideIntegrate(div_levs,quadr,c,c,ccenter,bc,energ_val,top_val,bot_val);
					sign = isCcw(c,ccenter,bc);
					energs[c.idx] += energ_val*sign;
					tops[c.idx] += top_val*sign;
					bots[c.idx] += bot_val*sign;
				} 
			}
		}
	}

	distr_grad.clear();
	distr_nPoints.clear();
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		//my_bots[(*point_itr).idx] = bots[(*point_itr).idx];
		my_bots[(*point_itr).idx] = tops[(*point_itr).idx].dot(*point_itr);
		// Ignoring boundary conditions at this moment
		my_energy += energs[(*point_itr).idx];
		if(!(*point_itr).isBdry){
			if(bots[(*point_itr).idx] != 0.0){
				grad_i = 2.0*(*point_itr)*bots[(*point_itr).idx] - 2.0*tops[(*point_itr).idx];
				grad_i.idx = (*point_itr).idx;
				distr_grad.push_back(grad_i);

				pnt_i = tops[(*point_itr).idx] / bots[(*point_itr).idx];
				pnt_i.normalize();
				pnt_i.idx = (*point_itr).idx;
				distr_nPoints.push_back(pnt_i);
			}
		} else if((*point_itr).isBdry){
			/*if((*point_itr).isBdry == 2){
				(*point_itr).isBdry = 0;
			}
			pnt zero(0.0,0.0,0.0);
			distr_grad.push_back(zero); */
		}
	}

	delete(tops);
	delete(bots);
	delete(energs);

	#ifdef _DEBUG
		cerr << "Done Integrating regions " << id << endl;
	#endif
	return;
}/*}}}*/











void computeMetrics(const int id, vector<pnt>& points, vector<pnt>& n_points, double &ave, double &max, double &l1){/*{{{*/
	//Metrics are computed for my updated points
	pnt norm_pt;
	double val;
	vector<pnt>::iterator point_itr;

	max = 0.0;
	ave = 0.0;
	l1 = 0.0;

#ifdef _DEBUG
	cerr << "Computing local metrics " << id << endl;
#endif

	for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
		norm_pt = points.at((*point_itr).idx) - (*point_itr);
		val = norm_pt.magnitude();

		ave += val*val;
		max = std::max(val,max);
		l1 += val;
	}
#ifdef _DEBUG
	cerr << "Done Computing local metrics " << id << endl;
#endif
	return;
}/*}}}*/




void projectToBoundary(const double proj_alpha, vector<pnt>& points, vector<pnt>& boundary_points, vector<pnt>& n_points, vector<region> &region_vec){/*{{{*/
	double new_min_dist, min_dist, dist, alpha, beta;
	pnt p, p_n;
	pnt a, b, c;
	pnt v1, v2, v3;
	vector<int> closest_cell;
	vector<int> proj_1_point;
	vector<int> proj_2_point;
	vector<int> search_flags;
	int i, index1, index2;
	vector<region>::iterator region_itr;
	vector<pnt>::iterator boundary_itr;
	vector<pnt>::iterator point_itr;


	for(i = 0; i < boundary_points.size(); i++){
		closest_cell.push_back(-1);
	}
	for(i = 0; i < points.size(); i++){
		proj_1_point.push_back(-1);
		proj_2_point.push_back(-1);
		search_flags.push_back(1);
	}

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); region_itr++){
		for(boundary_itr = (*region_itr).boundary_points.begin(); boundary_itr != (*region_itr).boundary_points.end(); boundary_itr++){

			min_dist = M_PI;
			for(point_itr = n_points.begin(); point_itr != n_points.end(); point_itr++){
				dist = (*boundary_itr).dotForAngle((*point_itr));

				if(dist < min_dist){
					min_dist = dist;
					closest_cell.at((*boundary_itr).idx) = (*point_itr).idx;
				}
			}

		}
	}

	for(point_itr = n_points.begin(); point_itr != n_points.end(); point_itr++){
		min_dist = M_PI;

		for(i = 0; i < closest_cell.size(); i++){
			if(closest_cell.at(i) == (*point_itr).idx){
				dist = (*point_itr).dotForAngle(boundary_points.at(i));

				if(dist < min_dist){
					min_dist = dist;
					proj_1_point.at((*point_itr).idx) = i;
				}
			}
		}

		if(proj_1_point.at((*point_itr).idx) > -1){
			closest_cell.at(proj_1_point.at((*point_itr).idx)) = -1;
		}

		min_dist = M_PI;

		for(i = 0; i < closest_cell.size(); i++){
			if(closest_cell.at(i) == (*point_itr).idx){
				dist = (*point_itr).dotForAngle(boundary_points.at(i));

				if(dist < min_dist){
					min_dist = dist;
					proj_2_point.at((*point_itr).idx) = i;
				}
			}
		}

	}

	for(point_itr = n_points.begin(); point_itr != n_points.end(); point_itr++){
		if(proj_1_point.at((*point_itr).idx) > -1){
			index1 = proj_1_point.at((*point_itr).idx);
			index2 = proj_2_point.at((*point_itr).idx);
			if(index2 > -1){
				a = boundary_points.at(index1);
				b = boundary_points.at(index2);

				v1 = a - (*point_itr);
				v2 = b - a;

				alpha = v2.magnitude();
				alpha = alpha*alpha;
				alpha = -v1.dot(v2)/alpha;

				if(alpha < 0.0){
					p = boundary_points.at(index1);
					p.idx = (*point_itr).idx;
					p.isBdry = 0;
					p.normalize();
				} else {
					p = a + alpha * (b - a);

					p.idx = (*point_itr).idx;
					p.isBdry = 0;
					p.normalize();
				}
			} else {
				p = boundary_points.at(index1);
				p.idx = (*point_itr).idx;
				p.isBdry = 0;
				p.normalize();
			}

			alpha = std::min(1.0, proj_alpha);
			p_n = alpha*p + (1.0 - alpha) * (*point_itr);
			p_n.normalize();
			p_n.idx = (*point_itr).idx;
			p_n.isBdry = 0;


			(*point_itr).x = p_n.x;
			(*point_itr).y = p_n.y;
			(*point_itr).z = p_n.z;
			(*point_itr).idx = p_n.idx;
			(*point_itr).isBdry = p_n.isBdry;
		}
	}

}/*}}}*/


void annealPoints(const int anneal_shape_control, const double anneal_percent, vector<region> &region_vec){/*{{{*/
	double rand_x, rand_y, rand_z;
	vector<region>::iterator region_itr;
	vector<pnt>::iterator point_itr;

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); region_itr++){
		for(point_itr = (*region_itr).points.begin(); point_itr != (*region_itr).points.end(); point_itr++){
			if(anneal_shape_control == 1){
				if((*point_itr).sides != 6){
					rand_x = drand48() * 2.0 * anneal_percent - anneal_percent;
					rand_y = drand48() * 2.0 * anneal_percent - anneal_percent;
					rand_z = drand48() * 2.0 * anneal_percent - anneal_percent;

					(*point_itr).x = (*point_itr).x + rand_x;
					(*point_itr).y = (*point_itr).y + rand_y;
					(*point_itr).z = (*point_itr).z + rand_z;

					(*point_itr).normalize();
				}
			} else {
				rand_x = drand48() * 2.0 * anneal_percent - anneal_percent;
				rand_y = drand48() * 2.0 * anneal_percent - anneal_percent;
				rand_z = drand48() * 2.0 * anneal_percent - anneal_percent;

				(*point_itr).x = (*point_itr).x + rand_x;
				(*point_itr).y = (*point_itr).y + rand_y;
				(*point_itr).z = (*point_itr).z + rand_z;

				(*point_itr).normalize();
			}
		}
	}

}/*}}}*/
/* }}} */


/* ***** Communication Routines ***** {{{*/
void transferUpdatedPoints(const mpi::communicator& world, vector<region>& my_regions, vector<pnt>& n_points, vector<pnt>& points){/*{{{*/
	//Each processor transfers it's updated point set (stored in n_points) to it's region neighbors
	//as defined in RegionTriangulation
	//
	//This keeps communications minimal and only updates the points that need to be updated for each region.
	vector<pnt> temp_points_in;
	vector<pnt> temp_points_out;
	vector<mpi::request> comms;
	optional options;
	vector<pnt>::iterator point_itr;
	vector<region>::iterator region_itr;

#ifdef _DEBUG
	cerr << "Transfering updated points " << world.rank() << endl;
#endif

	if(world.size() > 1){
		temp_points_out.clear();
		for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
			points.at((*point_itr).idx) = (*point_itr);
			temp_points_out.push_back((*point_itr));
		}

		for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
			comms.resize((*region_itr).neighbors.size());

			for(int i = 0; i < (*region_itr).neighbors.size(); i++){
				comms[i] = world.isend((*region_itr).neighbors.at(i), msg_points, temp_points_out);
			}

			for(int i = 0; i < (*region_itr).neighbors.size(); i++){
				temp_points_in.clear();
				world.recv((*region_itr).neighbors.at(i), msg_points, temp_points_in);

				for(point_itr = temp_points_in.begin(); point_itr != temp_points_in.end(); ++point_itr){
					points.at((*point_itr).idx) = (*point_itr);
				}
			}

			mpi::wait_all(&comms[0],&(comms[(*region_itr).neighbors.size()]));

			comms.clear();
		}
	} else {
		for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
			points.at((*point_itr).idx) = (*point_itr);
		}
	}

	temp_points_in.clear();
	temp_points_out.clear();
#ifdef _DEBUG
	cerr << "Done Transfering updated points " << world.rank() << endl;
#endif
}/*}}}*/


void gatherAllUpdatedPoints(const mpi::communicator& world, vector<pnt>& n_points, vector<pnt>& points){/*{{{*/
	//All updated points (stored in n_points) are gathered onto the master processor and broadcasted to every processor
	vector<pnt> temp_points;
	mpi::request mycomm;
	optional options;
	vector<pnt>::iterator point_itr;
	int id=world.rank(), num_procs=world.size();

#ifdef _DEBUG
	cerr << "Gatering all updated points " << id << endl;
#endif

	if(num_procs > 1){
		temp_points.clear();
		for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
			temp_points.push_back((*point_itr));
			if(id == 0){
				points.at((*point_itr).idx) = (*point_itr);
			}
		}
		for(int i = 1; i < num_procs; i++){
			if(id == i){
				mycomm = world.isend(0, msg_points, temp_points);
				points.clear();
			}else if(id == 0){
				world.recv(i, msg_points, temp_points);
				for(point_itr = temp_points.begin(); point_itr != temp_points.end(); ++point_itr){
					points.at((*point_itr).idx) = (*point_itr);
				}
				temp_points.clear();
			}
		}
		if(id != 0){
			points.clear();
			options = mycomm.test();
			if(!options){
				mycomm.wait();
			}
			temp_points.clear();
		}
		mpi::broadcast(world,points,0);
	} else {
		for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
			points.at((*point_itr).idx) = (*point_itr);
		}
	}
#ifdef _DEBUG
	cerr << "Done Gatering all updated points " << id << endl;
#endif
}/*}}}*/
/*}}}*/




#endif /* LOOP_ROUTINES_H_ */
