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

#pragma once
#include "matrix_equation.h"

#include <vector>
#include <array>
using namespace std;

typedef array<double, 6> baseCamera;
typedef array<double, 3> point;
typedef struct _observer{ double u; double v; int cam; int pt; } observer;
typedef array<double, 9> balCamera;

class BASolver
{
public:
	BASolver(double* _K);
	~BASolver() {}
	void addCamera(baseCamera& _mcamera) { cameras.push_back(_mcamera); }//SE3d
	void addPoint(point& _mpoint) { points.push_back(_mpoint); }
	void addObserver(observer& _mobserver) { observers.push_back(_mobserver); }
	//Levenberg¨CMarquardt Method
	void solveProblem(const int iteration);
protected:
	void sortObservers();
	//return value : updata matrix
	clMat* getSolveMatrix(double& cost, double& Q_down);
	//we will use seven steps to get the resule. 
	//step 1: calculating the inverse of H22
	clMat* getH22Inverse(const clMat* H22);
	//step 2: H12 * H22.inverse, we named the result T
	clMat* getMatrixT(const clMat* H12, const clMat* H22inverse);
	//step 3: add T * b2 to b1
	void getNewb1(const clMat* T, const clMat* b2, clMat* b1);
	//step 4: get Matrix A = T * H12.transpose
	clMat* getMatrixA(const clMat* T, const clMat* H12);
	//step 5: solve the equivalent (H11 - A) * x1 = b1
	clMat* getMatrix_x1(const clMat* A, const clMat* H11, const clMat* b1);
	//step 6: b2 -= H12.transpose * x1  
	void getNewb2(const clMat* H12, const clMat* x1, clMat* b2);
	//step 7: x2 = H22.inverse * b2
	clMat* getMatrix_x2(const clMat* H22inverse, const clMat* b2);
	//compute update for every camera and point
	void computeUpdate(clMat* update, double cost, double Q_down, int v);
	//use OpenCV to show the structure of the Matrix mat
	void showMatrixStructure(clMat& mat);
	vector<baseCamera> cameras;
	vector<point> points;
	vector<observer> observers;
	double mu;
	array<double, 4> K;

};