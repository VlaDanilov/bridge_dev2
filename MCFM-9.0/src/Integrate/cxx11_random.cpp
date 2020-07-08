/*
!  Copyright (C) 2019, respective authors of MCFM.
!
!  This program is free software: you can redistribute it and/or modify it under
!  the terms of the GNU General Public License as published by the Free Software
!  Foundation, either version 3 of the License, or (at your option) any later
!  version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!  PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program. If not, see <http://www.gnu.org/licenses/>
*/
 
#include <random>
#include <omp.h>
#include <iostream>

static std::vector<std::mt19937_64> engines;
static std::vector<std::uniform_real_distribution<double>> distribs;

extern "C" {

void cxx11_init_random(int* seeds)
{
  for (int i=0; i<omp_get_max_threads(); ++i) {
    engines.emplace_back(seeds[i]);
    distribs.emplace_back(0.0,1.0);
  }
}

double cxx11_random_number()
{
  int tid = omp_get_thread_num();
  return (distribs[tid])(engines[tid]);
}

}
