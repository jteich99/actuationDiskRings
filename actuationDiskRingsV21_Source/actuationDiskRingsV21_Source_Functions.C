// #include "actuationDiskRingsV21_Source.H"
#include "volFields.H"
#include <functional>
#include <math.h>
#include "fvc.H"
#include "fvCFD.H"
#include <map>

float posInTableUref2(
    int posI,
    List<scalar> Uref2List_,
    List<scalar> UrefList_,
    List<scalar> rList_,
    scalar rtable,
    vector U_dPointCells_ntr,
    List<scalar> UdiList_
) {
    int difference = 1000;
    int pos = 0;
    int lastPos = 0;

    for (int i = 0; i < (Uref2List_.size() - 1); i = i + 1)
    {
        if (
            (Uref2List_[i] == UrefList_[posI]) and
            (rList_[i] == rtable) and
            (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) &&
            ((mag(U_dPointCells_ntr) - UdiList_[i]) >= 0)) // look into calibration table 2 (Udi_table) for the Uref value, and search for the Udi closest to the Udi of the node
        {
            difference = fabs(mag(U_dPointCells_ntr) - UdiList_[i]);
            pos = i; // the position of the lower value
        }

        if (
            (Uref2List_[i] == UrefList_[posI]) and 
            (rList_[i] == rtable) and 
            (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && 
            ((mag(U_dPointCells_ntr) - UdiList_[i]) < 0)
        )
        {
            lastPos = i;
        }
    }

    if (difference == 1000)
    {
        // just in case the value is outside the table
        pos = lastPos;
    }

    return pos;
}

float gFunction(
    float x,
    float rootDistance
) {
    float a = 2.335;
    int b = 4;
    // Info << "x = " << x << endl;
    // Info << "rootDistance = " << rootDistance << endl;
    // Info << "x/rootDistance = " << x/rootDistance << endl;
    // Info << "pow((x/rootDistance), b) = " << pow((x/rootDistance), b) << endl;
    float g = 1 - std::exp( - a * pow((x/rootDistance), b) );
    return g;
}

float FFunction(
    float x,
    float lambda
){
    float Nb = 3;
    // Info << "1 + pow(lambda,2) = " << 1 + pow(lambda,2) << endl;
    // Info << "std::exp( - (Nb/2) * std::sqrt(1 + pow(lambda,2)) * (1 - x) ) = " << std::exp( - (Nb/2) * std::sqrt(1 + pow(lambda,2)) * (1 - x) ) << endl;
    float F = (2/M_PI) * std::acos( std::exp( - (Nb/2) * std::sqrt(1 + pow(lambda,2)) * (1 - x) ) );
    return F;
}

float a1Function(
    float rootDistance,
    float lambda
){
    float a1 = 0;
    float resolution = 100;
    float gPrev;
    float FPrev;
    float funcPrev;
    float g;
    float F;
    float func;
    for (float x = 2/resolution; x < 1; x += 1/resolution) {
        // Info << "x - 1/resolution = " << x - 1/resolution << endl;
        gPrev = gFunction(x - 1/resolution, rootDistance);
        FPrev = FFunction(x - 1/resolution, lambda);
        funcPrev = pow(gPrev,2) * pow(FPrev,2) / (x - 1/resolution);
        g = gFunction(x, rootDistance);
        F = FFunction(x, lambda);
        // Info << "x = " << x << endl;
        func = pow(g,2) * pow(F,2) / x;
        a1 += (1/resolution) * ( func + funcPrev )/2;          
    }
    return a1;
}

float a2Function(
    float rootDistance,
    float lambda
){
    float a2 = 0;
    float resolution = 100;
    float gPrev;
    float FPrev;
    float funcPrev;
    float g;
    float F;
    float func;
    for (float x = 2/resolution; x < 1; x += 1/resolution) {
        gPrev = gFunction(x - 1/resolution, rootDistance);
        FPrev = FFunction(x - 1/resolution, lambda);
        funcPrev = gPrev * FPrev * (x - 1/resolution);
        g = gFunction(x, rootDistance);
        F = FFunction(x, lambda);
        func = g * F * x;
        a2 += (1/resolution) * ( func + funcPrev )/2;          
    }
    return a2;
}
