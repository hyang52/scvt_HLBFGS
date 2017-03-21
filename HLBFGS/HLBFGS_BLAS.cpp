///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLBFGS                                                                    //
// http://www.loria.fr/~liuyang/software/HLBFGS/							 //
//                                                                           //
// HLBFGS is a hybrid L-BFGS optimization framework which unifies L-BFGS     //
// method, Preconditioned L-BFGS method and                                  //
// Preconditioned Conjugate Gradient method.                                 //
//                                                                           //
// Version 1.2                                                               //
// March 09, 2010                                                            //
//                                                                           //
// Copyright (C) 2009--2010                                                  //
// Yang Liu                                                                  //
//																			 //
// xueyuhanlang@gmail.com                                                    //
//                                                                           //
// HLBFGS is HLBFGS is freely available for non-commercial purposes.		 //
//                                                                           //
// Maintainer: Huanhuan Yang. huan2yang@outlook.com							 //
//			   Customize for mpi handeling									 //
///////////////////////////////////////////////////////////////////////////////
#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "HLBFGS_BLAS.h"
#include <cmath>

double HLBFGS_DDOT(const int n, const double *x, const double *y, mpi::communicator* world)
{
	double result = 0;
	int i = 0;
#ifdef USE_OPENMP
#pragma omp parallel for private(i) reduction(+:result)
#endif
	for (i = 0; i < n; i++)
	{
		result += x[i] * y[i];
	}

	if(world==0)
	{
		return result;
	}else{
		double g_result;
		mpi::reduce(*world, result, g_result, std::plus<double>(), 0);
		mpi::broadcast(*world, g_result, 0);	
		return g_result;
	}
}

void HLBFGS_DAXPY(const int n, const double alpha, const double *x, double *y)
{
	int i = 0;
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
	for (i = 0; i < n; i++)
	{
		y[i] += alpha * x[i];
	}
}

double HLBFGS_DNRM2(const int n, const double *x, mpi::communicator* world)
{
	double result = 0;
	int i = 0;
#ifdef USE_OPENMP
#pragma omp parallel for private(i) reduction(+:result)
#endif
	for (i = 0; i < n; i++)
	{
		result += x[i] * x[i];
	}

	if(world==0)
	{
		return std::sqrt(result);
	}else{
		double g_result;
		mpi::reduce(*world, result, g_result, std::plus<double>(), 0);
		mpi::broadcast(*world, g_result, 0);	
		return std::sqrt(g_result);
	}
}

void HLBFGS_DSCAL(const int n, const double a, double *x)
{
	int i = 0;
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
	for (i = 0; i < n; i++)
	{
		x[i] *= a;
	}
}

