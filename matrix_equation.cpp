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

#include "matrix_equation.h"

void matrix_joint(const clMat& H, const clMat& b, clMat& des)
{
	des.resize(H.getRow(), H.getCol() + b.getCol());
	int row = des.getRow();
	int h_col = H.getCol();
	int b_col = b.getCol();
	for (int i = 0; i < row; ++i)
	{
#ifdef _MSC_VER
	
		for (int j = 0; j < h_col; ++j)
			des.data[i][j] = H.data[i][j];
		for (int j = 0; j < b_col; ++j)
			des.data[i][j + h_col] = b.data[i][j];

#else
		memcpy(des.data[0] + i * (h_col + b_col), H.data[0] + i * h_col, sizeof(double) * row * h_col);
		memcpy(des.data[0] + i * (h_col + b_col) + h_col,
			b.data[0] + i * b_col,
			sizeof(double) * row * b_col);
#endif
	}
}

void first_step(clMat& des)
{
	int N = des.getRow();
	int M = des.getCol();
	for (int i = 0; i < N; ++i)
	{
		//第一行
		for (int c = i + 1; c < M; ++c)
			des.data[i][c] /= des.data[i][i];
		//后面每一行
		for (int r = i + 1; r < N; ++r)
		{
			for (int c = i + 1; c < M; ++c)
			{
				if (des.data[r][i]) {
					des.data[r][c] /= des.data[r][i];
					des.data[r][c] -= des.data[i][c];
				}
			}
		}
	}
	for (int r = N - 2; r > -1; --r)
	{//被求解的矩阵
		for (int c = N; c < M; ++c)
		{
			for (int i = N - 1; i > r; --i)
			{
				des.data[r][c] -= des.data[i][c] * des.data[r][i];
			}			
		}
	}
}
