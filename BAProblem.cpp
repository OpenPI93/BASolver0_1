/**
* This file is part of MT-BA.
*
* Copyright 2018 Cao Lin.
* Developed by Cao Lin ,
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* MT-BA is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* MT-BA is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with MT-BA. If not, see <http://www.gnu.org/licenses/>.
*/

#include"BAProblem.h"
#include <opencv.hpp>
#include <sophus\se3.hpp>
#include <eigen3\Eigen\Core>
#include <omp.h>
#include <thread>

BASolver::BASolver(double* _K)
{ 
	K[0] = _K[0];	//fx
	K[1] = _K[1];	//fy
	K[2] = _K[2];	//cx
	K[3] = _K[3];	//cy
	//init mu as anything < 0, because the mu in LMA is always > 0 so that we can compute mu = tao * max{H[i][i]} before the first iteration where tao = 0.000006.
	mu = -1;
}

void BASolver::sortObservers()
{
	std::sort(observers.begin(), observers.end(), [&](observer a, observer b) {return a.cam < b.cam; });
}

clMat* BASolver::getSolveMatrix(double& cost, double& Q_down)
{
	int cameracount = cameras.size();
	int pointcount = points.size();
	int observercount = observers.size();
	// int the equivalent H * x = b ,
	//the H Matrix has special structure, we define 1 as the pose of camera and 2 means the location of the point. so that
	//
	//		H11   H12         x1         b1
	//					  *          =
	//	   H12^T  H22         x2         b2
	//
	
	clMat* H11 = new clMat(cameracount * 6, cameracount * 6);
	clMat* H12 = new clMat(cameracount * 6, pointcount * 3);
	clMat* H22 = new clMat(pointcount * 3, pointcount * 3);

	clMat* b1 = new clMat(cameracount * 6, 1);
	clMat* b2 = new clMat(pointcount * 3, 1);

	/////////////////////////////get Matrix////////////////////////////////////////
	
	//compute jacobian from the observer info
	int index = 0;
//#pragma omp parallel for
	for (auto edge : observers)
	{
		//we will use Sophus library to solve the SE3 and Eigen to store data;
		//step 1: prepare data
		Sophus::SE3d T = Sophus::SE3d::exp(Eigen::Matrix<double, 6, 1>(cameras[edge.cam].data()));
		Eigen::Vector3d ori_pt(points[edge.pt].data());
		Eigen::Vector3d cur_pt = T * ori_pt;

		double fx = K[0];
		double fy = K[1];
		double x = cur_pt[0];
		double y = cur_pt[1];
		double z = cur_pt[2];
		double inv_z = 1. / cur_pt[2];
		double inv_z2 = inv_z * inv_z;

		//step 2: compute Jacobian for each observe

		//∂e/∂ξ = ∂e/∂p' * ∂p'/∂ξ
		//∂e/∂p = ∂e/∂p' * R
		//∂e/∂p'
		Eigen::Matrix<double, 2, 3> tmp1;
		tmp1(0, 0) = fx * inv_z;
		tmp1(0, 1) = 0;
		tmp1(0, 2) = -fx * x * inv_z2;

		tmp1(1, 0) = 0;
		tmp1(1, 1) = fy * inv_z;
		tmp1(1, 2) = -fy * y * inv_z2;

		//∂p'/∂ξ
		Eigen::Matrix<double, 3, 6> tmp2;
		tmp2(0, 0) = 1;
		tmp2(0, 1) = 0;
		tmp2(0, 2) = 0;
		tmp2(0, 3) = 0;
		tmp2(0, 4) = z;
		tmp2(0, 5) = -y;

		tmp2(1, 0) = 0;
		tmp2(1, 1) = 1;
		tmp2(1, 2) = 0;
		tmp2(1, 3) = -z;
		tmp2(1, 4) = 0;
		tmp2(1, 5) = x;

		tmp2(2, 0) = 0;
		tmp2(2, 1) = 0;
		tmp2(2, 2) = 1;
		tmp2(2, 3) = y;
		tmp2(2, 4) = -x;
		tmp2(2, 5) = 0;

		//step 3: compute the Jacobian Matrix for this observe
		Eigen::Matrix<double, 2, 6> Jc = -tmp1 * tmp2;
		Eigen::Matrix<double, 2, 3> Jp = -tmp1 * T.rotationMatrix();

		//step4: write to the H Matrix
		int cam(edge.cam);
		int pt(edge.pt);

		//H11
#pragma UNROLL(6)
		for(int r = 0 ; r < 6 ; ++r)
		{			
			H11->data[cam * 6 + r][cam * 6 + 0] += Jc(0, r) * Jc(0, 0) + Jc(1, r) * Jc(1, 0);
			H11->data[cam * 6 + r][cam * 6 + 1] += Jc(0, r) * Jc(0, 1) + Jc(1, r) * Jc(1, 1);
			H11->data[cam * 6 + r][cam * 6 + 2] += Jc(0, r) * Jc(0, 2) + Jc(1, r) * Jc(1, 2);
			H11->data[cam * 6 + r][cam * 6 + 3] += Jc(0, r) * Jc(0, 3) + Jc(1, r) * Jc(1, 3);
			H11->data[cam * 6 + r][cam * 6 + 4] += Jc(0, r) * Jc(0, 4) + Jc(1, r) * Jc(1, 4);
			H11->data[cam * 6 + r][cam * 6 + 5] += Jc(0, r) * Jc(0, 5) + Jc(1, r) * Jc(1, 5);
		}
	
		//H22
#pragma UNROLL(3)
		for (int r = 0; r < 3; ++r)
		{
			H22->data[pt * 3 + r][pt * 3 + 0] += Jp(0, r) * Jp(0, 0) + Jp(1, r) * Jp(1, 0);
			H22->data[pt * 3 + r][pt * 3 + 1] += Jp(0, r) * Jp(0, 1) + Jp(1, r) * Jp(1, 1);
			H22->data[pt * 3 + r][pt * 3 + 2] += Jp(0, r) * Jp(0, 2) + Jp(1, r) * Jp(1, 2);
		}

		//H12
#pragma UNROLL(6)
		for (int r = 0; r < 6; ++r)
		{
			H12->data[cam * 6 + r][pt * 3 + 0] = Jc(0, r) * Jp(0, 0) + Jc(1, r) * Jp(1, 0);
			H12->data[cam * 6 + r][pt * 3 + 1] = Jc(0, r) * Jp(0, 1) + Jc(1, r) * Jp(1, 1);
			H12->data[cam * 6 + r][pt * 3 + 2] = Jc(0, r) * Jp(0, 2) + Jc(1, r) * Jp(1, 2);
		}
		//step 5: compute error
		cur_pt = cur_pt * inv_z;
		Eigen::Vector2d uv;
		//we need to compute b = -e * J, this part is -e
		uv[0] =	cur_pt[0] * fx + K[2] - edge.u;
		uv[1] = cur_pt[1] * fy + K[3] - edge.v;

		// step 6: write to the b Matrix
		
		b1->data[cam * 6 + 0][0] += Jc(0, 0) * uv[0] + Jc(1, 0) * uv[1];
		b1->data[cam * 6 + 1][0] += Jc(0, 1) * uv[0] + Jc(1, 1) * uv[1];
		b1->data[cam * 6 + 2][0] += Jc(0, 2) * uv[0] + Jc(1, 2) * uv[1];
		b1->data[cam * 6 + 3][0] += Jc(0, 3) * uv[0] + Jc(1, 3) * uv[1];
		b1->data[cam * 6 + 4][0] += Jc(0, 4) * uv[0] + Jc(1, 4) * uv[1];
		b1->data[cam * 6 + 5][0] += Jc(0, 5) * uv[0] + Jc(1, 5) * uv[1];

		b2->data[pt * 3 + 0][0] += Jp(0, 0) * uv[0] + Jp(1, 0) * uv[1];
		b2->data[pt * 3 + 1][0] += Jp(0, 1) * uv[0] + Jp(1, 1) * uv[1];
		b2->data[pt * 3 + 2][0] += Jp(0, 2) * uv[0] + Jp(1, 2) * uv[1];

		cost += uv[0] * uv[0] + uv[1] * uv[1];
	}
	/////////////////end of get Matrix/////////////////////////
	//init mu
	if (mu < 0)
	{
		double max = -999999999999999;
		for (int i = 0; i < H11->getRow(); ++i)
			if (max < H11->data[i][i])
				max = H11->data[i][i];
		for (int i = 0; i < H22->getRow(); ++i)
			if (max < H22->data[i][i])
				max = H22->data[i][i];
		mu = 0.000001 * max;
	}
	for (int i = 0; i < H11->getRow(); ++i)
		H11->data[i][i] += mu;
	for (int i = 0; i < H22->getRow(); ++i)
		H22->data[i][i] += mu;

	clMat* H22inverse = getH22Inverse(H22);
	clMat* T = getMatrixT(H12, H22inverse);
	getNewb1(T, b2, b1);
	clMat* A = getMatrixA(T, H12);
	clMat* x1 = getMatrix_x1(A, H11, b1);
	getNewb2(H12, x1, b2);
	clMat* x2 = getMatrix_x2(H22inverse, b2);

	clMat* aimMatrix = new clMat(x1->getRow() + x2->getRow(), 1);
	for (int i = 0; i < x1->getRow(); ++i) 
	{
		aimMatrix->data[i][0] = x1->data[i][0];
		Q_down += x1->data[i][0] * (mu * x1->data[i][0] - b1->data[i][0]);
	}
	int x1row = x1->getRow();
	for (int i = 0; i < x2->getRow(); ++i)
	{
		aimMatrix->data[i + x1row][0] = x2->data[i][0];
		Q_down += x2->data[i][0] * (mu * x2->data[i][0] + b2->data[i][0]);
	}

//#ifdef _DEBUG
//	showMatrixStructure(*H22);
//#endif // _DEBUG
	delete H22;
	delete H11;
	delete H12;
	delete H22inverse;
	delete b1;
	delete b2;
	delete A;
	delete x1;
	delete x2;
	delete T;

	return aimMatrix;
}

void BASolver::solveProblem(const int iteration)
{
	cout << "there are " << cameras.size() << " camera";
	if (cameras.size() > 1)cout << "s, ";
	else cout << ", ";
	cout << points.size() << " point";
	if (points.size() > 1)cout << "s ";
	else cout << " ";
	cout << "and " << observers.size() << " observer";
	if (observers.size() > 1)cout << "s.\n";
	else cout << ".\n";

	double cost = 0.0, lastcost = 0.0;
	int v = 2;
	//Q = (F(x) - F(x_new)) / (0.5 * x ^ T * (mu * x - b)), Q_down means (x ^ T * (mu * x - b)) which we can compute in getSolveMatrix, and F(x) = ||error||2
	double Q_down;

	for (int iter = 0; iter < 100; iter++)
	{
		Q_down = 0.0;
		clMat* update = getSolveMatrix(cost, Q_down);
		
		if (iter > 0 && cost > lastcost) {
			cout << "lastcost > cost\n";
			break;
		}
		computeUpdate(update, cost, Q_down, v);
		lastcost = cost;	
		cout << "iteration " << iter << " cost=" << cout.precision(12) << cost << endl;	
		cost = 0.0;
	}
	Eigen::Matrix<double, 6, 1> data;
	for (int i = 0; i < 6; ++i)
		data[i] = cameras[0][i];
	cout << Sophus::SE3d::exp(data).matrix() << endl;;
}

clMat* BASolver::getH22Inverse(const clMat* H22)
{
	clMat* H22inverse = new clMat(H22->getRow(), H22->getCol());
	int cols = points.size() * 3;
	//calculating the inverse of the three order matrix by multi-threading
	auto f = [&](int pt_start, int pt_end) 
	{
		// a0  a1  a3
		// a1  a2  a4
		// a3  a4  a5
		for (int pt = pt_start; pt < pt_end; ++pt) 
		{
			double a0 = H22->data[pt * 3][pt * 3];
			double a1 = H22->data[pt * 3][pt * 3 + 1];
			double a2 = H22->data[pt * 3 + 1][pt * 3 + 1];
			double a3 = H22->data[pt * 3][pt * 3 + 2];
			double a4 = H22->data[pt * 3 + 1][pt * 3 + 2];
			double a5 = H22->data[pt * 3 + 2][pt * 3 + 2];
			double det_a = a0 * a2 * a5 + 2 * a1 * a4 * a3 - a3 * a3* a2 - a5 * a1 * a1 - a0 * a4 * a4;
			double inv_det_a = 1. / det_a;
			
			H22inverse->data[pt * 3][pt * 3] = (a2 * a5 - a4 * a4) * inv_det_a;
			H22inverse->data[pt * 3][pt * 3 + 1] =
				H22inverse->data[pt * 3 + 1][pt * 3] =
				(a4 * a3 - a1 * a5) * inv_det_a;
			H22inverse->data[pt * 3 + 1][pt * 3 + 1] = (a0 * a5 - a3 * a3) * inv_det_a;
			H22inverse->data[pt * 3][pt * 3 + 2] =
				H22inverse->data[pt * 3 + 2][pt * 3] =
				(a1 * a4 - a2 * a3) * inv_det_a;
			H22inverse->data[pt * 3 + 1][pt * 3 + 2] =
				H22inverse->data[pt * 3 + 2][pt * 3 + 1] =
				(a1 * a3 - a0 * a4) * inv_det_a;
			H22inverse->data[pt * 3 + 2][pt * 3 + 2] = (a0 * a2 - a1 * a1) * inv_det_a;
		}
	};
#ifdef MULTI_THREADING
	int temp = points.size() / 4;
	std::array<std::thread, 4> threads;
	threads[0] = (std::thread(f, temp * 0, temp * 1));
	threads[1] = (std::thread(f, temp * 1, temp * 2));
	threads[2] = (std::thread(f, temp * 2, temp * 3));
	threads[3] = (std::thread(f, temp * 3, points.size()));

	threads[0].join();
	threads[1].join();
	threads[2].join();
	threads[3].join();
#else
	f(0, points.size());
#endif
	
	return H22inverse;
}

clMat* BASolver::getMatrixT(const clMat* H12, const clMat* H22inverse)
{
	clMat* T = new clMat(H12->getRow(), H12->getCol());

	//compute Matrix T by multi-threading
	auto f = [&](int index_start, int index_end)
	{
		//for each observer
		for (int i = index_start; i < index_end; ++i)
		{
			int cam = observers[i].cam;
			int pt = observers[i].pt;

			for (int r = cam * 6; r < (cam + 1) * 6; ++r)
			{
				for (int c = pt * 3; c < (pt + 1) * 3; ++c)
				{
					T->data[r][c] = 
						H12->data[r][pt * 3 + 0] * H22inverse->data[pt * 3 + 0][c] +
						H12->data[r][pt * 3 + 1] * H22inverse->data[pt * 3 + 1][c] +
						H12->data[r][pt * 3 + 2] * H22inverse->data[pt * 3 + 2][c];
				}
			}
		}
	};
#ifdef MULTI_THREADING
	int temp = observers.size() / 4;
	std::array<std::thread, 4> threads;
	threads[0] = (std::thread(f, temp * 0, temp * 1));
	threads[1] = (std::thread(f, temp * 1, temp * 2));
	threads[2] = (std::thread(f, temp * 2, temp * 3));
	threads[3] = (std::thread(f, temp * 3, observers.size()));

	threads[0].join();
	threads[1].join();
	threads[2].join();
	threads[3].join();
#else
	f(0, observers.size());
#endif
	return T;
}

void BASolver::getNewb1(const clMat* T, const clMat* b2, clMat* b1)
{
	//T->getRow()
	int Trow = T->getRow();
	//T->getCol()
	int Tcol = T->getCol();

	auto f = [&](int start, int end)
	{
		for (int i = start; i < end; ++i)
		{
			for (int j = 0; j < Tcol; ++j)
			{
				b1->data[i][0] -= T->data[i][j] * b2->data[j][0];
			}
		}
	};

#ifdef MULTI_THREADING
	int temp = Trow / 4;
	std::array<std::thread, 4> threads;
	threads[0] = (std::thread(f, temp * 0, temp * 1));
	threads[1] = (std::thread(f, temp * 1, temp * 2));
	threads[2] = (std::thread(f, temp * 2, temp * 3));
	threads[3] = (std::thread(f, temp * 3, Trow));

	threads[0].join();
	threads[1].join();
	threads[2].join();
	threads[3].join();
#else
	f(0, Trow);
#endif
}

clMat* BASolver::getMatrixA(const clMat* T, const clMat* H12)
{
	//T->getRow()
	int Trow = T->getRow();
	//T->getCol()
	int Tcol = T->getCol();
	clMat *A = new clMat(Trow, Trow);

	//  Matrix T and H12
	//					c
	//	start	*****************
	//	end		*****************
	//			*****************		i
	//			*****************
	//			*****************
	//			
	// the resault for each f
	//					i
	//	start	*****************
	//	end		*****************
	auto f = [&](int start, int end)
	{
		for (int r = start; r < end; ++r)
		{
			for (int i = 0; i < Trow; ++i)
			{
				for (int c = 0; c < Tcol; ++c)
				{
					A->data[r][i] += T->data[r][c] * H12->data[i][c];
				}
			}
		}
	};
#ifdef MULTI_THREADING
	int temp = Trow / 4;
	std::array<std::thread, 4> threads;
	threads[0] = (std::thread(f, temp * 0, temp * 1));
	threads[1] = (std::thread(f, temp * 1, temp * 2));
	threads[2] = (std::thread(f, temp * 2, temp * 3));
	threads[3] = (std::thread(f, temp * 3, Trow));

	threads[0].join();
	threads[1].join();
	threads[2].join();
	threads[3].join();
#else
	f(0, Trow);
#endif

	return A;
}

clMat* BASolver::getMatrix_x1(const clMat* A, const clMat* H11, const clMat* b1)
{
	clMat solver(A->getRow(), A->getCol() + 1);
	int Arow = A->getRow();
	int Acol = A->getCol();
	//we use elementary transformation to solve the matrix equation, f1 to merge H11 - A | b1.
	auto f = [&](int start, int end)
	{
		for (int r = start; r < end; ++r)
		{
			for (int c = 0; c < Acol; ++c)
				solver.data[r][c] = H11->data[r][c] - A->data[r][c];
			solver.data[r][Acol] = b1->data[r][0];
		}
	};
#ifdef MULTI_THREADING
	int temp = Arow / 4;
	std::array<std::thread, 4> threads;
	threads[0] = (std::thread(f, temp * 0, temp * 1));
	threads[1] = (std::thread(f, temp * 1, temp * 2));
	threads[2] = (std::thread(f, temp * 2, temp * 3));
	threads[3] = (std::thread(f, temp * 3, Arow));

	threads[0].join();
	threads[1].join();
	threads[2].join();
	threads[3].join();
#else
	f(0, Arow);
#endif
	
	//solve matrix equation by elementary tansformation. 
	int N = Arow;
	int M = Arow + 1;
	for (int i = 0; i < N; ++i)
	{
		for (int c = i + 1; c < M; ++c)
			solver.data[i][c] /= solver.data[i][i];
		for (int r = i + 1; r < N; ++r)
		{
#ifdef MULTI_THREADING
#ifdef _OPENMP
#pragma omp parallel for num_threads(4)
#endif
#endif
			for (int c = i + 1; c < M; ++c)
			{
				if (solver.data[r][i]) {
					solver.data[r][c] /= solver.data[r][i];
					solver.data[r][c] -= solver.data[i][c];
				}
			}
		}
	}
	clMat* x1 = new clMat(N, 1);
	x1->data[N - 1][0] = solver.data[N - 1][N];
	for (int r = N - 2; r > -1; --r)
	{		
#ifdef MULTI_THREADING
#ifdef _OPENMP
#pragma omp parallel for num_threads(4)
#endif
#endif
		for (int i = N - 1; i > r; --i)
		{
			solver.data[r][N] -= solver.data[i][N] * solver.data[r][i];
		}
		x1->data[r][0] = solver.data[r][N];
	}

	return x1;
}

void BASolver::getNewb2(const clMat* H12, const clMat* x1, clMat* b2)
{
	//b2 -= H12.transpose * x1;
	int row = H12->getRow();
	int col = H12->getCol();

	auto f = [&](int start, int end)
	{
		for (int r = start; r < end; ++r)
		{
			for (int c = 0; c < row; ++c)
				b2->data[r][0] -= H12->data[c][r] * x1->data[c][0];
		}
	};
#ifdef MULTI_THREADING
	int temp = col / 4;
	std::array<std::thread, 4> threads;
	threads[0] = (std::thread(f, temp * 0, temp * 1));
	threads[1] = (std::thread(f, temp * 1, temp * 2));
	threads[2] = (std::thread(f, temp * 2, temp * 3));
	threads[3] = (std::thread(f, temp * 3, col));

	threads[0].join();
	threads[1].join();
	threads[2].join();
	threads[3].join();
#else
	f(0, col);
#endif
}

clMat* BASolver::getMatrix_x2(const clMat* H22inverse, const clMat* b2)
{
	int pt = points.size();
	clMat* x2 = new clMat(pt * 3, 1);

	auto f = [&](int start, int end)
	{
		for (int i = start; i < end; ++i)
		{
#pragma UNROLL(3)
			for (int r = 0; r < 3; ++r) 
			{
#pragma UNROLL(3)
				for (int c = 0; c < 3; ++c)
				{
					x2->data[i * 3 + r][0] += 
						H22inverse->data[i * 3 + r][i * 3 + c] * b2->data[i * 3 + c][0];
				}
			}				
		}
	};
#ifdef MULTI_THREADING
	int temp = points.size() / 4;
	std::array<std::thread, 4> threads;
	threads[0] = (std::thread(f, temp * 0, temp * 1));
	threads[1] = (std::thread(f, temp * 1, temp * 2));
	threads[2] = (std::thread(f, temp * 2, temp * 3));
	threads[3] = (std::thread(f, temp * 3, points.size()));

	threads[0].join();
	threads[1].join();
	threads[2].join();
	threads[3].join();
#else
	f(0, points.size());
#endif

	return x2;
}

void BASolver::showMatrixStructure(clMat& mat)
{
	cv::Mat outputH(mat.getRow(), mat.getRow(), CV_8UC1);
	for (int i = 0; i < outputH.rows; ++i)
	{
		for (int j = 0; j < outputH.cols; ++j)
		{
			if (mat.data[i][j])
				outputH.at<uchar>(i, j) = 0;
			else
				outputH.at<uchar>(i, j) = 255;
		}
	}
	cv::imshow("Matrix", outputH);
	cv::imwrite(".\Matrix.png", outputH);
	cv::waitKey(1);
}

void BASolver::computeUpdate(clMat* update, double cost, double Q_down, int v)
{
	double cost_new = 0.0;
	std::vector< Eigen::Matrix<double, 6, 1> > cam_new;
	std::vector<point> point_new;
	for (int i = 0; i < cameras.size(); ++i)
	{
		Sophus::SE3d T = Sophus::SE3d::exp(Eigen::Matrix<double, 6, 1>(cameras[i].data()));
		double cur_cam[6];
		for (int j = 0; j < 6; ++j)
			cur_cam[j] = update->data[6 * i + j][0];
		Sophus::SE3d T_update = Sophus::SE3d::exp(Eigen::Matrix<double, 6, 1>(cur_cam));
		T = T_update * T;		
		Eigen::Matrix<double, 6, 1> value = T.log();
		cam_new.push_back(value);
	}
	int temp = cameras.size() * 6;
	for (int i = 0; i < points.size(); ++i)
	{
		point temp_point;
		temp_point[0] = points[i][0] + update->data[temp + i * 3 + 0][0];
		temp_point[1] = points[i][1] + update->data[temp + i * 3 + 1][0];
		temp_point[2] = points[i][2] + update->data[temp + i * 3 + 2][0];
		point_new.push_back(temp_point);
	}

	//compute the new cost
	for (auto edge : observers)
	{		
		Sophus::SE3d T = Sophus::SE3d::exp(cam_new[edge.cam]);
		Eigen::Vector3d ori_pt(point_new[edge.pt].data());
		Eigen::Vector3d cur_pt = T * ori_pt;

		double fx = K[0];
		double fy = K[1];
		double z = cur_pt[2];
		double inv_z = 1. / cur_pt[2];

		cur_pt = cur_pt * inv_z;
		Eigen::Vector2d uv;
		uv[0] = cur_pt[0] * fx + K[2] - edge.u;
		uv[1] = cur_pt[1] * fy + K[3] - edge.v;

		cost_new += uv[0] * uv[0] + uv[1] * uv[1];
	}

	double Q = (cost - cost_new) / Q_down;

	if (Q > 0)
	{
		for (int i = 0; i < cameras.size(); ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				cameras[i][j] = cam_new[i][j];
			}
		}
		for (int i = 0; i < points.size(); ++i)
		{
			points[i][0] = point_new[i][0];
			points[i][1] = point_new[i][1];
			points[i][2] = point_new[i][2];
		}
		double temp = 1 - pow((2 * Q - 1), 3);
		mu *= (0.3333333 > temp ? 0.3333333 : temp);
		v = 2;
	}
	else
	{
		mu *= v;
		v *= 2;
	}
	cout << "mu = " << mu << "  ";

	delete update;
}