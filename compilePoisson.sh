PROC=4
REAL_TYPE="DOUBLE"
#REAL_TYPE="SINGLE"
#PLANAR="PLANAR"
#TIME_RUN="TIME_RUN"
#TEST_RUN="TEST_RUN"
EXECUTE_AFTER_COMPILE="EXECUTE"

# NO USER MODIFICATIONS NECESSARY BELOW THIS LINE

# REMOVE PRE-EXISTING EXECUTATBLES
rm -f test_poisson

# IF build DIRECTORY DOESN'T EXIST, CREATE IT
if [ ! -d build ]; then
    mkdir build
fi

# SWITCH TO build DIRECTORY
cd build

# RUN Cmake WITH NECESSARY FLAGS AS SET BY USER
if [ -z $PLANAR ]; then
    if [ -z $TEST_RUN ]; then
        if [ -z $TIME_RUN ]; then
            if [ "$REAL_TYPE" == "DOUBLE" ]; then
                CC=mpicc CXX=mpicxx cmake ../ -DREAL=FLOAT
            else
                CC=mpicc CXX=mpicxx cmake ../ -DREAL=DOUBLE
            fi
        else
            CC=mpicc CXX=mpicxx cmake ../ -DTIME_RUN=ON
        fi
    else
        CC=mpicc CXX=mpicxx cmake ../ -DTEST_RUN=ON
    fi
else
    if [ -z $TEST_RUN ]; then
        if [ -z $TIME_RUN ]; then
            CC=mpicc CXX=mpicxx cmake ../ -DPLANAR=ON
        else
            CC=mpicc CXX=mpicxx cmake ../ -DPLANAR=ON -DTIME_RUN=ON
        fi
    else
        CC=mpicc CXX=mpicxx cmake ../ -DTEST_RUN=ON -DPLANAR=ON
    fi
fi

# COMPILE
make -j4

#move executable outside build
mv ./test_poisson ../

# SWITCH TO PARENT DIRECTORY
cd ../

# RUN CODE IF REQUESTED BY USER
if ! [ -z $EXECUTE_AFTER_COMPILE ]; then
    if [ -z $TEST_RUN ]; then
        mpirun -np $PROC ./test_poisson
    else
        mpirun -np $PROC ./test_poisson
    fi
fi
