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

#include "matrix.h"
#include <iomanip>
#include <exception>

clMat::clMat(int _row, int _col) : row(_row), col(_col)
{
	if (row < 1) {
		std::cerr << "row is lower than 1 in " << __FILE__ << " line : " << __LINE__ << endl;
		throw std::runtime_error("row is lower than 1");
	}
	if (col < 1) {
		std::cerr << "col is lower than 1 in " << __FILE__ << " line : " << __LINE__ << endl;
		throw std::runtime_error("col is lower than 1");
	}

	data = new double*[row];
	data[0] = new double[row * col]();

	for (int i = 1; i < row; ++i)
	{
		data[i] = data[0] + i * col;
	}
}

clMat::~clMat()
{
	delete[]data[0];
	for (int i = 0; i < row; ++i)
		data[i] = nullptr;
	delete[]data;
	data = nullptr;
}

void clMat::setData(double* _data)
{
#ifdef _MSC_VER
	memcpy_s(data[0], sizeof(double) * row * col, _data, sizeof(double) * row * col);
#else
	memcpy(data[0], _data, sizeof(double) * row * col);
#endif // check the IDE
	
}

void clMat::resize(int _row, int _col)
{
	delete[]data[0];
	delete[]data;
	data = nullptr;
	row = _row;
	col = _col;

	if (row < 1) {
		std::cerr << "row is lower than 1 in " << __FILE__ << " line : " << __LINE__ << endl;
		throw std::runtime_error("row is lower than 1");
	}
	if (col < 1) {
		std::cerr << "col is lower than 1 in " << __FILE__ << " line : " << __LINE__ << endl;
		throw std::runtime_error("col is lower than 1");
	}

	data = new double*[row];
	data[0] = new double[row * col]();

	for (int i = 1; i < row; ++i)
	{
		data[i] = data[0] + i * col;
	}
}

void clMat::output(std::ostream& out)const
{
	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			out << std::setw(7) << data[i][j] << "  ";
		}
		out << endl;
	}
}

std::ostream& operator << (std::ostream& out, const clMat& x)
{
	x.output(out);
	return out;
}

clMat clMat::Identity(int _size)
{
	clMat* I = new clMat(_size, _size);
	for (int i = 0; i < _size; ++i)
		I->data[i][i] = 1;
	return *I;
}