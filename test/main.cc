/********************************************************************************************************************************************
 * Poisson
 * 
 * Copyright (C) 2020, Mahendra K. Verma
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */

/*! \file main.cc
 *
 *  \brief test file for poisoon library
 *
 *  \author Ali Asad and Roshan Samuel
 *  \date Feb 2020
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <poisson.h>

int main()
{
    MPI_Init(NULL, NULL);

    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD, &rnk);

	blitz::Array<int, 1> poissonInitInt(21);
    blitz::Array<real, 1> poissonInitReal(7);
    
    int xIndex = 6;
    int yIndex = 6;
    int zIndex = 6;
    int xMesh = 0;
    int yMesh = 0;
    int zMesh = 0;
    int nopX = 2;
    int nopY = 2;
    int nThreads = 1;
    int xPerInt = 0;
    int yPerInt = 0;
    int zPerInt = 0;
    int vcDpth = 3;
    int vcCnt = 5;
    int preSmth = 10;
    int postSmth = 10;
    std::vector<int> intrSmth;
    intrSmth.resize(5,0);
    intrSmth[0] = 10;
    intrSmth[1] = 10;
    intrSmth[2] = 10;
    intrSmth[3] = 10;
    intrSmth[4] = 0;

    real xL = 1.0;
    real yL = 1.0;
    real zL = 1.0;
    real btX = 1.0;
    real btY = 1.0;
    real btZ = 1.0;
    real tol = 0.00001;

    poissonInitInt(0) = xIndex;
    poissonInitInt(1) = yIndex;
    poissonInitInt(2) = zIndex;
    poissonInitInt(3) = xMesh;
    poissonInitInt(4) = yMesh;
    poissonInitInt(5) = zMesh;
    poissonInitInt(6) = nopX;
    poissonInitInt(7) = nopY;
    poissonInitInt(8) = nThreads;
    poissonInitInt(9) = xPerInt;
    poissonInitInt(10) = yPerInt;
    poissonInitInt(11) = zPerInt;
    poissonInitInt(12) = vcDpth;
    poissonInitInt(13) = vcCnt;
    poissonInitInt(14) = preSmth;
    poissonInitInt(15) = postSmth;
    poissonInitInt(16) = intrSmth[0];
    poissonInitInt(17) = intrSmth[1];
    poissonInitInt(18) = intrSmth[2];
    poissonInitInt(19) = intrSmth[3];
    poissonInitInt(20) = intrSmth[4];
    
    poissonInitReal(0) = xL;
    poissonInitReal(1) = yL;
    poissonInitReal(2) = zL;
    poissonInitReal(3) = btX;
    poissonInitReal(4) = btY;
    poissonInitReal(5) = btZ;
    poissonInitReal(6) = tol;

    
   /*****************************************************************************************************************************************
    * CALLING POISSON LIBRARY
    *****************************************************************************************************************************************
   */ 
 
    if(rnk == 0) {
        std::cout << "Hey let's call poisson library!" << std::endl; 
    }
    
    poisson poisson_lib;
    poisson_lib.Init(poissonInitInt, poissonInitReal);

    if(rnk == 0) {
        std::cout << poisson_lib.pressureData.size() << std::endl;
        std::cout << "Hey poisson worked!" << std::endl;
    }
    

//******************************************************************************************************
 
    if(rnk == 0) {
        std::cout << "Hey let's call multigrid_d2 from poisson library!" << std::endl; 
    } 

   multigrid_d2 multigrid_d2_lib;
   multigrid_d2_lib.Init_d2(poissonInitInt, poissonInitReal); 
   
   if(rnk == 0) {
        std::cout << multigrid_d2_lib.pressureData.shape() << std::endl;
        std::cout << "Hey multigrid_d2 class from poisson library worked!" << std::endl; 
    }


//*******************************************************************************************************

    if(rnk == 0) {
        std::cout << "Hey let's call multigrid_d3 from poisson library!" << std::endl; 
    } 

   multigrid_d3 multigrid_d3_lib;
   multigrid_d3_lib.Init_d3(poissonInitInt, poissonInitReal); 
   
   if(rnk == 0) {
        std::cout << multigrid_d3_lib.pressureData.shape() << std::endl;
        std::cout << "Hey multigrid_d3 class from poisson library worked!" << std::endl; 
    }
    
    MPI_Finalize();
  	return 0;
}
