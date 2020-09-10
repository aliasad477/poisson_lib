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
/*! \file poisson.h
 *
 *  \brief Class declaration of poisson
 *
 *  \author Ali Asad and Roshan Samuel
 *  \date Feb 2020
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef POISSON_H
#define POISSON_H

#include <blitz/array.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include "grid.h"

#ifdef REAL_DOUBLE
#define H5T_NATIVE_REAL H5T_NATIVE_DOUBLE
#define MPI_FP_REAL MPI_DOUBLE
#define real double
#else
#define H5T_NATIVE_REAL H5T_NATIVE_FLOAT
#define MPI_FP_REAL MPI_FLOAT
#define real float
#endif

class poisson {
    protected:
        int vLevel, maxCount;
        int xStr, yStr, zStr;
        int xEnd, yEnd, zEnd;

#ifdef TIME_RUN
        real solveTimeComp;
        real solveTimeTran;
        real smothTimeComp;
        real smothTimeTran;
#endif

        //parser variables..........................................................................
        int xInd, yInd, zInd;
        int vcDepth, vcCount;
        int preSmooth, postSmooth;
        int xGrid, yGrid, zGrid;
        int xPerInt, yPerInt, zPerInt;
        bool xPer, yPer, zPer;
        int nThreads;
        int npY, npX;
        std::vector<int> interSmooth;

        real tolerance;
        real xLen, yLen, zLen;        
        real betaX, betaY, betaZ;

        static inline int pmod(int a, int b) {return (a % b + b) % b;};

        inline int findRank(int xR, int yR) {return pmod(yR, npY)*npX + pmod(xR, npX);};

        //..........................................................................................

        blitz::Array<int, 1> nearRanks;

        blitz::Array <int, 1> gridInitInt;
        blitz::Array <real, 1> gridInitReal;
        
        blitz::Array<real, 3> residualData;
        blitz::Array<real, 3> iteratorTemp;
        blitz::Array<real, 3> smoothedPres;

        blitz::Array<int, 1> mgSizeArray;
        blitz::Array<int, 1> strideValues;

        blitz::TinyVector<int, 3> localSizeIndex;

        blitz::Array<MPI_Request, 1> recvRequest;
        blitz::Array<MPI_Status, 1> recvStatus;

        blitz::Array<real, 1> hx, hy, hz;

        blitz::Array<real, 1> xixx, xix2;
        blitz::Array<real, 1> etyy, ety2;
        blitz::Array<real, 1> ztzz, ztz2;

        blitz::Array<blitz::Range, 1> xMeshRange, yMeshRange, zMeshRange;

        blitz::Array<MPI_Datatype, 1> xMGArray;
        blitz::Array<MPI_Datatype, 1> yMGArray;

        blitz::Array<blitz::TinyVector<int, 3>, 1> mgSendLft, mgSendRgt;
        blitz::Array<blitz::TinyVector<int, 3>, 1> mgRecvLft, mgRecvRgt;

        blitz::Array<blitz::TinyVector<int, 3>, 1> mgSendFrn, mgSendBak;
        blitz::Array<blitz::TinyVector<int, 3>, 1> mgRecvFrn, mgRecvBak;

        virtual void solve();
        virtual void prolong();
        virtual void smooth(const int smoothCount);

        virtual void initMeshRanges();

        virtual void setStagBounds();
        virtual void setLocalSizeIndex();

        virtual void setCoefficients();
        virtual void copyStaggrDerivs();

        virtual void imposeBC();
        virtual void updatePads();
        virtual void createMGSubArrays();

        virtual void vCycle();

        grid grid_lib;

        void initializeArrays();
        void nearRanksDefinition();

    public:
        blitz::Array<real, 3> pressureData;
        blitz::Array<real, 3> inputRHSData;

        blitz::RectDomain<3> stagFull;
        blitz::RectDomain<3> stagCore;

        void Init(blitz::Array <int, 1> parserDataInt, blitz::Array <real, 1> parserDataReal);

        virtual void mgSolve(blitz::Array<real, 3> inFn, const blitz::Array<real, 3> rhs);

        virtual real testTransfer();
        virtual real testProlong();
        virtual real testPeriodic();
        virtual real testSolve();

        virtual ~poisson();
};

/**
 ********************************************************************************************************************************************
 *  \class poisson poisson.h "lib/poisson.h"
 *  \brief The base class poisson and its derived classes multigrid_d2 and multigrid_d3
 *
 *  The class implements the geometric multi-grid method for solving the Poisson equation on a non-uniform grid across MPI decomposed
 *  domains for parallel computations.
 *  The data structure used by the class for computing multi-grid V-cycles across sub-domains is a blitz array with very wide overlap.
 *  When calculating the finite differences at the sub-domain boundaries, at the coarsest level of the V-cycle, data points from very
 *  deep within the neighbouring sub-domains are necessary.
 *  This is the reason for using a wide pad, that spans up to the nearest node of the adjacent sub-domain at the coarsest mesh.
 *  This increases the memory footprint, but doesn't increase the computational time as only a single finite-difference calculation
 *  is being done using the pads at all levels of the V-cycle.
 *
 *  All the necessary functions to perform the V-cycle - prolongation, solving at coarsest mesh, smoothening, etc. are implemented
 *  within the \ref poisson class.
 ********************************************************************************************************************************************
 */

class multigrid_d2: public poisson {
    private:
        blitz::Array<real, 1> hx2, hz2, hzhx;

        void solve();
        void prolong();
        void smooth(const int smoothCount);

        void initMeshRanges();

        void setStagBounds();
        void setLocalSizeIndex();

        void setCoefficients();
        void copyStaggrDerivs();

        void imposeBC();
        void updatePads();
        void createMGSubArrays();

        void vCycle();

    public:
        void Init_d2(blitz::Array <int, 1> parserDataInt, blitz::Array <real, 1> parserDataReal);
        
        void mgSolve(blitz::Array<real, 3> inFn, const blitz::Array<real, 3> rhs);

        real testTransfer();
        real testProlong();
        real testPeriodic();
        real testSolve();

        ~multigrid_d2() {};
};

/**
 ********************************************************************************************************************************************
 *  \class multigrid_d2 poisson.h "lib/poisson.h"
 *  \brief The derived class from poisson to perform multi-grid operations on a 2D grid
 *
 *  The 2D implementation ignores the y-direction component of the computational domain.
 ********************************************************************************************************************************************
 */

class multigrid_d3: public poisson {
    private:
        blitz::Array<real, 1> hxhy, hyhz, hzhx, hxhyhz;

        void solve();
        void prolong();
        void smooth(const int smoothCount);

        void initMeshRanges();

        void setStagBounds();
        void setLocalSizeIndex();

        void setCoefficients();
        void copyStaggrDerivs();

        void imposeBC();
        void updatePads();
        void createMGSubArrays();

        void vCycle();

    public:
        void Init_d3(blitz::Array <int, 1> parserDataInt, blitz::Array <real, 1> parserDataReal);

        void mgSolve(blitz::Array<real, 3> inFn, const blitz::Array<real, 3> rhs);

        real testTransfer();
        real testProlong();
        real testPeriodic();
        real testSolve();

        ~multigrid_d3() {};
};

/**
 ********************************************************************************************************************************************
 *  \class multigrid_d3 poisson.h "lib/poisson.h"
 *  \brief The derived class from poisson to perform multi-grid operations on a 3D grid
 *
 *  The 3D implementation of the multi-grid method differs from the 2D version in that the \ref solve \ref smooth functions use a different
 *  equation with extra terms, and the \ref prolong operation needs to perform extra interpolation steps in the y-direction.
 ********************************************************************************************************************************************
 */

#endif
