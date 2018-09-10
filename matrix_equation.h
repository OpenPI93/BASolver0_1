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

#include <iostream>
using std::endl;
using std::cout;
#include "matrix.h"


void matrix_equation(clMat& H, clMat& b);
void matrix_joint(const clMat& H, const clMat& b, clMat& des);
void first_step(clMat& des);
void second_step(clMat& des);
