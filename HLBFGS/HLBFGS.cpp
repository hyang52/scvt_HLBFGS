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
//			   Customize the stopping criteria, mpi handeling, and resetting //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "HLBFGS.h"
#include "LineSearch.h"
#include "ICFS.h"

//////////////////////////////////////////////////////////////////////////
void INIT_HLBFGS(double PARAMETERS[], int INFO[])
{
	PARAMETERS[0] = 1.0e-4; //ftol
	PARAMETERS[1] = 1.0e-16; //xtol
	PARAMETERS[2] = 0.9; //gtol
	PARAMETERS[3] = 1.0e-20; //stpmin
	PARAMETERS[4] = 1.0e+20; //stpmax
	PARAMETERS[5] = 1.0e-16; // ||g||/max(1,||x||)
	PARAMETERS[6] = 1.0e-10; // ||g||
	PARAMETERS[7] = 1.0e-10; // ||g^k||/||f^k||
	PARAMETERS[8] = 1.0e-16; // ||f^k+1 - f^k||/||f^k||
	PARAMETERS[9] = 1.0e-16; // ||x^k+1 - x^k||_l2

	INFO[0] = 20; //max_fev_in_linesearch
	INFO[1] = 0; //total_num_fev
	INFO[2] = 0; //iter
	INFO[3] = 0; //update strategy. 0: standard lbfgs, 1: m1qn3, 2: Lloyd precondition
	INFO[4] = 100000; // max iterations
	INFO[5] = 1; //1: print message, 0: do nothing
	INFO[6] = 10; // T: update interval of Hessian
	INFO[7] = 0; // 0: without hessian, 1: with accurate hessian
	INFO[8] = 15; // icfs parameter
	INFO[9] = 0; // 0: linesearch 1: modified linesearch (do not set 1 in pratice !)
	INFO[10] = 0; // 0: Disable preconditioned CG 1: preconditioned CG
	INFO[11] = 1; // different methods for choosing beta in CG.
	INFO[12] = 1; //internal usage. 0: only update diag in USER_DEFINED_HLBFGS_UPDATE_H
	INFO[13] = 0; // 0: standard lbfgs update, 1: Biggs's update, 2: Yuan's update; 3: Zhang and Xu's update
	INFO[14] = 0; // max num of resetting LBFGS, i.e. setting iter=0 and total_num_fev=0
}

//////////////////////////////////////////////////////////////////////////
void HLBFGS_MESSAGE(bool print, int id, const double PARAMETERS[], mpi::communicator* world)
{
	if (!print)
	{
		return;
	}

	if(world==0 || (world!=0 && world->rank()==0))
	{
		switch (id)
		{
		case 0:
			std::cout << "Please check your input parameters !\n";
			break;
		case 1:
			std::cout << "Convergence : ||g||/max(1,||x||) <= " << PARAMETERS[5]
			<< std::endl;
			break;
		case 2:
			std::cout << "Convergence : ||g|| <=  " << PARAMETERS[6] << std::endl;
			break;
		case 3:
			std::cout << "Convergence : ||g^k||/||f^k|| <=  " << PARAMETERS[7] << std::endl;
			break;
		case 4:
			std::cout << "Warning: linesearch cannot improve anymore!\n";
			break;
		case 5:
			std::cout << "************Exceeds max resetting." << std::flush;
			std::cout <<" If still not convergent, you may need to relax Lloyd_tol or set larger max_reset. \n";
			break;
		case 6:
			std::cout << "Exceeds max iteration \n"<< std::endl;
			break;
		case 7:
			std::cout << "Convergence : ||x^k+1 - x^k||_l2 <=  " << PARAMETERS[9] << std::endl;
			break;
		case 8:
			std::cout << "Convergence : ||f^k+1 - f^k||/||f_^k|| <=  " << PARAMETERS[8] << std::endl;
			break;
		default:
			break;
		}
	}
}

//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_Hessian(int N, int M, double *q, double *s, double *y,
						   int cur_pos, double *diag, int INFO[], double *x, mpi::communicator* world)
{
	if (M <= 0 || INFO[2] == 0)
	{
		return;
	}

	int start = cur_pos * N;

	double *y_start = &y[start];
	double *s_start = &s[start];

	double ys = HLBFGS_DDOT(N, y_start, s_start, &(*world));

	if (INFO[3] == 0)
	{
		double yy = HLBFGS_DDOT(N, y_start, y_start, &(*world));
		double factor = ys / yy;
		if (INFO[12] == 1)
		{
			HLBFGS_DSCAL(N, factor, q);
		}
		else
		{
			diag[0] = factor;
		}

	}
	else if (INFO[3] == 1)
	{

		//m1qn3 update
		double dyy, my_dyy= 0;
		double dinvss, my_dinvss = 0;
		int i = 0;
#ifdef USE_OPENMP
#pragma omp parallel for private(i) reduction(+:dinvss) reduction(+:dyy)
#endif
		for (i = 0; i < N; i++)
		{
			my_dinvss += s_start[i] * s_start[i] / diag[i];
			my_dyy += diag[i] * y_start[i] * y_start[i];
		}
		mpi::reduce(*world, my_dyy, dyy, std::plus<double>(), 0);
		mpi::reduce(*world, my_dinvss, dinvss, std::plus<double>(), 0);
		mpi::broadcast(*world, dyy, 0);
		mpi::broadcast(*world, dinvss, 0);
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
		for (i = 0; i < N; i++)
		{
			diag[i] = 1.0 / (dyy / (ys * diag[i]) + y_start[i] * y_start[i]
			/ ys - dyy * s_start[i] * s_start[i] / (ys * dinvss
				* diag[i] * diag[i]));
		}
		if (INFO[12] == 1)
		{
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
			for (i = 0; i < N; i++)
			{
				q[i] *= diag[i];
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_First_Step(int N, int M, double *q, double *s, double *y,
							  double *rho, double *alpha, int bound, int cur_pos, int iter, mpi::communicator* world)
{
	if (M <= 0)
	{
		return;
	}

	int start;
	double tmp;

	for (int i = bound; i >= 0; i--)
	{
		start = iter <= M ? cur_pos - bound + i : (cur_pos - (bound - i) + M)
			% M;
		alpha[i] = rho[start] * HLBFGS_DDOT(N, q, &s[start * N], &(*world));
		tmp = -alpha[i];
		HLBFGS_DAXPY(N, tmp, &y[start * N], q);
	}
}

//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_Second_Step(int N, int M, double *q, double *s, double *y,
							   double *rho, double *alpha, int bound, int cur_pos, int iter, mpi::communicator* world)
{
	if (M <= 0)
	{
		return;
	}

	int start;
	double tmp;

	for (int i = 0; i <= bound; i++)
	{
		start = iter <= M ? i : (cur_pos + 1 + i) % M;
		tmp = alpha[i] - rho[start] * HLBFGS_DDOT(N, &y[start * N], q, &(*world));
		HLBFGS_DAXPY(N, tmp, &s[start * N], q);
	}
}

//////////////////////////////////////////////////////////////////////////
void HLBFGS_BUILD_HESSIAN_INFO(HESSIAN_MATRIX& m_hessian, int INFO[])
{
	ICFS_INFO& l_info = m_hessian.get_icfs_info();
	l_info.get_p() = INFO[8];
	l_info.set_lrow_ind_size(m_hessian.get_nonzeros()
		+ m_hessian.get_dimension() * l_info.get_p());
	l_info.set_l_size(m_hessian.get_nonzeros() + m_hessian.get_dimension()
		* l_info.get_p());
	l_info.get_icfs_alpha() = 0;
	int n = m_hessian.get_dimension();
	int nnz = m_hessian.get_nonzeros();
	dicfs_(&n, &nnz, m_hessian.get_values(), m_hessian.get_diag(),
		m_hessian.get_colptr(), m_hessian.get_rowind(), l_info.get_l(),
		l_info.get_ldiag(), l_info.get_lcol_ptr(), l_info.get_lrow_ind(),
		&l_info.get_p(), &l_info.get_icfs_alpha(), l_info.get_iwa(),
		l_info.get_wa1(), l_info.get_wa2());
}

//////////////////////////////////////////////////////////////////////////
void CONJUGATE_GRADIENT_UPDATE(int N, double *q, double *prev_q_update,
							   double *prev_q_first_stage, int INFO[], mpi::communicator* world)
{
	//determine beta
	double cg_beta = 1.0;
	if (INFO[11] == 1)
	{
		if (INFO[2] == 0)
		{
			std::copy(q, &q[N], prev_q_first_stage);
			std::copy(q, &q[N], prev_q_update);
			return;
		}
		else
		{
			cg_beta = HLBFGS_DDOT(N, q, q, &(*world));
			cg_beta /= std::fabs(cg_beta
				- HLBFGS_DDOT(N, q, prev_q_first_stage, &(*world)));
			std::copy(q, &q[N], prev_q_first_stage);
		}
	}
	else
	{
		if (INFO[2] == 0)
		{
			std::copy(q, &q[N], prev_q_update);
			return;
		}
	}
	//determine new q

	const double minus_one = -1.0;
	if (cg_beta != 1.0)
		HLBFGS_DSCAL(N, cg_beta, prev_q_update);

	int i = 0;
//#ifdef USE_OPENMP
//#pragma omp parallel for private(i)
//#endif
	for (i = 0; i < N; i++)
	{
		q[i] -= prev_q_update[i];
	}

	double quad_a = HLBFGS_DDOT(N, q, q, &(*world));
	double quad_b = HLBFGS_DDOT(N, q, prev_q_update, &(*world));
	double cg_lambda = -quad_b / quad_a;
	if (cg_lambda > 1)
		cg_lambda = 1;
	else if (cg_lambda < 0)
		cg_lambda = 0;

#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
	for (i = 0; i < N; i++)
	{
		q[i] = cg_lambda * q[i] + prev_q_update[i];
	}
	std::copy(q, &q[N], prev_q_update);
}

//////////////////////////////////////////////////////////////////////////
int HLBFGS(int N, int M, double *x, int EVALFUNC(int, double*, double*,
			double*, double*), void EVALFUNC_H(int, double*, double*, double*,
			double*, HESSIAN_MATRIX&), void USER_DEFINED_HLBFGS_UPDATE_H(int, int,
			double*, double*, double*, int, double*, int[], double*, mpi::communicator*), void NEWITERATION(int,
			int, double*, double*, double*, double*), double PARAMETERS[],
			int INFO[], mpi::communicator* WORLD)
{
	int T = INFO[6];
	if (N < 1 || M < 0 || T < 0 || INFO[4] < 1)
	{
		HLBFGS_MESSAGE(INFO[5] != 0, 0, PARAMETERS, &(*WORLD));
		exit(1);
	}
	//allocate mem
	std::vector<double> q_vec(N), g_vec(N), alpha_vec(M <= 0 ? 0 : M), rho_vec(
		M <= 0 ? 0 : M), s_vec(M <= 0 ? 0 : M * N), y_vec(M <= 0 ? 0 : M
		* N), prev_x_vec(N), prev_g_vec(N), diag_vec, wa_vec(N);
	double *q = &q_vec[0];
	double *g = &g_vec[0];
	double *alpha = M <= 0 ? 0 : &alpha_vec[0];
	double *rho = M <= 0 ? 0 : &rho_vec[0];
	double *s = M <= 0 ? 0 : &s_vec[0];
	double *y = M <= 0 ? 0 : &y_vec[0];
	double *prev_x = &prev_x_vec[0];
	double *prev_g = &prev_g_vec[0];
	double *diag = 0;
	double *wa = &wa_vec[0];
	double update_alpha = 1;
	HESSIAN_MATRIX m_hessian(N);
	if (INFO[3] == 1)
	{
		diag_vec.resize(N);
		diag = &diag_vec[0];
		std::fill(diag_vec.begin(), diag_vec.end(), 1.0);
	}
	std::vector<double> prev_q_first_stage_vec, prev_q_update_vec;
	double *prev_q_first_stage = 0;
	double *prev_q_update = 0;
	double scale = 0.0;
	double cg_dginit = 0;
	if (INFO[10] == 1)
	{
		if (INFO[11] == 1)
		{
			prev_q_first_stage_vec.resize(N);
			prev_q_first_stage = &prev_q_first_stage_vec[0];
		}
		prev_q_update_vec.resize(N);
		prev_q_update = &prev_q_update_vec[0];
	}

	//initialize
	INFO[1] = 0;
	INFO[2] = 0;
	int num_reset = 0;
	int ret;
	double f = 0;
	int maxfev = INFO[0], bound = 0, nfev = 0, cur_pos = 0, start = 0;
	//line search parameters
	double stp, ftol = PARAMETERS[0], xtol = PARAMETERS[1], gtol =
		PARAMETERS[2], stpmin = PARAMETERS[3], stpmax = PARAMETERS[4];
	int info, keep[20];
	double gnorm, rkeep[40];
	std::fill(&rkeep[0], &rkeep[40], 0.0);
	std::fill(&keep[0], &keep[20], 0);

	m_hessian.get_icfs_info().allocate_mem(N);
	char task1 = 'N';
	char task2 = 'T';
	double prev_f;
	int i;

	//////////////////////////////////////////////////////////////////////////
	do
	{
		if (INFO[7] == 1 && ((T == 0) || (INFO[2] % T == 0)))
		{
			EVALFUNC_H(N, x, INFO[2] == 0 ? 0 : prev_x, &f, g, m_hessian);
			HLBFGS_BUILD_HESSIAN_INFO(m_hessian, INFO);
		}
		else if (INFO[2] == 0)
		{
			EVALFUNC(N, x, 0, &f, g);
			INFO[1]++;
		}

		if (INFO[2] > 0 && M > 0)
		{
			//compute s and y
			start = cur_pos * N;
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
			for (i = 0; i < N; i++)
			{
				s[start + i] = x[i] - prev_x[i];
				y[start + i] = g[i] - prev_g[i];
			}
			rho[cur_pos] = 1.0 / HLBFGS_DDOT(N, &y[start], &s[start], &(*WORLD));
			if (INFO[13] == 1)
			{
				update_alpha = 1.0 / (rho[cur_pos] * 6 * (prev_f - f
					+ HLBFGS_DDOT(N, g, &s[start], &(*WORLD))) - 2.0);
			}
			else if (INFO[13] == 2)
			{
				update_alpha = 1.0 / (rho[cur_pos] * 2 * (prev_f - f
					+ HLBFGS_DDOT(N, g, &s[start], &(*WORLD))));
			}
			else if (INFO[13] == 3)
			{
				update_alpha = 1.0 / (1 + rho[cur_pos] * (6 * (prev_f - f) + 3
					* (HLBFGS_DDOT(N, g, &s[start], &(*WORLD)) + HLBFGS_DDOT(N,
					prev_g, &s[start], &(*WORLD)))));
			}
			if (INFO[13] != 0)
			{
				if (update_alpha < 0.01)
				{
					update_alpha = 0.01;
				}
				else if (update_alpha > 100)
				{
					update_alpha = 100;
				}
				rho[cur_pos] *= update_alpha;
			}
		}
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
		for (i = 0; i < N; i++)
		{
			q[i] = -g[i];
		}

		if (INFO[2] > 0 && M > 0)
		{
			bound = INFO[2] > M ? M - 1 : INFO[2] - 1;
			HLBFGS_UPDATE_First_Step(N, M, q, s, y, rho, alpha, bound, cur_pos,
				INFO[2], &(*WORLD));
		}

		if (INFO[10] == 0)
		{
			if (INFO[7] == 1)
			{
				ICFS_INFO& l_info = m_hessian.get_icfs_info();
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task1);
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task2);
			}
			else
			{
				USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, diag, INFO, x, &(*WORLD));
			}
		}
		else
		{
			if (INFO[7] == 1)
			{
				ICFS_INFO& l_info = m_hessian.get_icfs_info();
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task1);
				CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update,
					prev_q_first_stage, INFO);
				cg_dginit = -HLBFGS_DDOT(N, q, q, &(*WORLD));
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task2);
			}
			else
			{
				INFO[12] = 0;
				USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, INFO[3]
				== 0 ? (&scale) : diag, INFO, x, &(*WORLD));
				if (INFO[3] == 0)
				{
					if (M > 0 && INFO[2] > 0 && scale != 1.0)
					{
						scale = std::sqrt(scale);
						HLBFGS_DSCAL(N, scale, q);
					}
					CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update,
						prev_q_first_stage, INFO);
					cg_dginit = -HLBFGS_DDOT(N, q, q, &(*WORLD));
					if (M > 0 && INFO[2] > 0 && scale != 1.0)
						HLBFGS_DSCAL(N, scale, q);
				}
				else
				{
					if (M > 0 && INFO[2] > 0)
					{
						//use prev_g as temporary array
						for (int i = 0; i < N; i++)
						{
							prev_g[i] = std::sqrt(diag[i]);
							q[i] *= prev_g[i];
						}
					}
					CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update,
						prev_q_first_stage, INFO);
					cg_dginit = -HLBFGS_DDOT(N, q, q, &(*WORLD));
					if (M > 0 && INFO[2] > 0)
					{
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
						for (i = 0; i < N; i++)
						{
							q[i] *= prev_g[i];
						}
					}

				}
				INFO[12] = 1;
			}
		}

		if (INFO[2] == 0 && num_reset > 0 && f >= prev_f)
		{
			if(WORLD->rank()==0)	std::cout << "Convergence: cannot improve anymore!\n"<< std::endl;
			return 0;
		}
		if (INFO[2] == 0 && HLBFGS_DDOT(N, g, q, &(*WORLD)) >= 0.)
		{
		 	HLBFGS_DAXPY(N, 1.0, q, x);
		 	continue;
		}

		if (INFO[2] > 0 && M > 0)
		{
			HLBFGS_UPDATE_Second_Step(N, M, q, s, y, rho, alpha, bound,
				cur_pos, INFO[2], &(*WORLD));
			cur_pos = (cur_pos + 1) % M;
		}

		//store g and x
		std::copy(&x[0], &x[N], prev_x_vec.begin());
		std::copy(&g[0], &g[N], prev_g_vec.begin());
		prev_f = f;
		//linesearch, find new x
		bool blinesearch = true;
		if (INFO[2] == 0)
		{
			gnorm = HLBFGS_DNRM2(N, g, &(*WORLD));
			stp = 1.0;// / gnorm;
		}
		else
		{
			stp = 1;
		}

		info = 0;

		do
		{
			MCSRCH(&N, x, &f, g, q, &stp, &ftol, &gtol, &xtol, &stpmin,
				&stpmax, &maxfev, &info, &nfev, wa, keep, rkeep, INFO[10]
			== 0 ? 0 : (&cg_dginit), &(*WORLD));
			blinesearch = (info == -1);
			if (blinesearch)
			{
				ret = EVALFUNC(N, x, prev_x, &f, g);
				if(ret == 1)
					return 1;
				INFO[1]++;
			}

			if (INFO[9] == 1 && prev_f > f) //modify line search to avoid too many function calls
			{
				info = 1;
				break;
			}

		} while (blinesearch);

		gnorm = HLBFGS_DNRM2(N, g, &(*WORLD));
		INFO[2]++;
		NEWITERATION(INFO[2], INFO[1], x, &f, g, &gnorm);
		double xnorm = HLBFGS_DNRM2(N, x, &(*WORLD));
		xnorm = 1 > xnorm ? 1 : xnorm;
		double dxnorm = HLBFGS_DNRM2(N, q, &(*WORLD))*stp;
		//if(WORLD->rank()==0)	std::cout << "polar dxnorm = " << dxnorm << std::endl;
		rkeep[2] = gnorm;
		rkeep[8] = xnorm;

		// The following stopping criteria should be checked in order.
		if (gnorm / xnorm <= PARAMETERS[5])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 1, PARAMETERS, &(*WORLD));
			break;
		}
		if (gnorm < PARAMETERS[6])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 2, PARAMETERS, &(*WORLD));
			break;
		}
		if (gnorm/f <= PARAMETERS[7])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 3, PARAMETERS, &(*WORLD));
			break;
		}
		if (info != 1 || stp < stpmin || stp > stpmax)
		{
			//if(WORLD->rank()==0)	std::cout << "info from LS:" << info << std::endl;
			HLBFGS_MESSAGE(INFO[5] != 0, 4, PARAMETERS, &(*WORLD));

			if(num_reset < INFO[14]){
				// if fail in line search, reset LBFGS
				if(WORLD->rank()==0)	std::cout << "--------------Reseting HLBFGS: let itr=0 -----------------\n";
				cur_pos = 0;
				INFO[2] = 0;
				prev_f = f;
				USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, diag, INFO, x, &(*WORLD));
				HLBFGS_DAXPY(N, 1.0, q, x);
				num_reset++;
				continue;
			}else{
				HLBFGS_MESSAGE(INFO[5] != 0, 5, PARAMETERS, &(*WORLD));
				break;
			}
		}
		if (INFO[2] > INFO[4])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 6, PARAMETERS, &(*WORLD));

			if(num_reset < INFO[14]){
				// if reach max_itr, reset LBFGS
				if(WORLD->rank()==0)	std::cout << "--------------Reseting HLBFGS: let itr=0 -----------------\n";
				cur_pos = 0;
				INFO[2] = 0;
				prev_f = f;
				USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, diag, INFO, x, &(*WORLD));
				HLBFGS_DAXPY(N, 1.0, q, x);
				num_reset++;
				continue;
			}else{
				HLBFGS_MESSAGE(INFO[5] != 0, 5, PARAMETERS, &(*WORLD));
				break;
			}
		}
		if (dxnorm <= PARAMETERS[9])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 7, PARAMETERS, &(*WORLD));
			break;
		}
		if (std::abs( (f-prev_f)/prev_f ) <= PARAMETERS[8])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 8, PARAMETERS, &(*WORLD));
			break;
		}
	} while (true);

	return 0;
}
