// #include "actuationDiskRingsV21_Source.H"
#include "volFields.H"
#include <functional>
#include <math.h>
#include "fvc.H"
#include "fvCFD.H"
#include <map>

vector getNodePosition(
    scalar tita_n_Rad,
    scalar rMed_r,
    scalar yawRad,
    vector diskPoint_,
    int ring,
    int numberRings_
) {
    scalar x_node = 0;
    scalar y_node = 0;
    scalar z_node = 0;

    if (ring != numberRings_)
    {
        x_node = -1 * rMed_r * std::sin(tita_n_Rad) * std::sin(yawRad);
        y_node = rMed_r * std::sin(tita_n_Rad) * std::cos(yawRad);
        z_node = rMed_r * std::cos(tita_n_Rad);
    }

    x_node += diskPoint_[0];
    y_node += diskPoint_[1];
    z_node += diskPoint_[2];

    vector Bi = vector(x_node, y_node, z_node);
    return Bi;
}

vector getNodeVelocity(
    List<label> nodeCellID_,
    int total_nodes_counter,
    int ring,
    int numberRings_,
    vectorField U,
    volTensorField gradU,
    int nodesNumber_,
    int gradInterpolation_,
    vectorField cellCentres,
    vector Bi
) {
    vector U_dPointCells = vector(1000, 1000, 1000);
    if (ring == numberRings_)
    {
        if (nodeCellID_[nodesNumber_-1] != -1) //if the closer cell is in this procesor
        {
            U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
        }
    }
    else
    {
        if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
        {
            U_dPointCells = U[nodeCellID_[total_nodes_counter]];
        }
    }
    if ( 
        (gradInterpolation_ == 1) and
        (ring != numberRings_)
    ){
        // Info << "cell center = " << mesh().cellCentres()[nodeCellID_[total_nodes_counter]] << endl;
        // Info << "cell center = " << cells[nodeCellID_[total_nodes_counter]] << endl;
        // vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
        vector dx = Bi - cellCentres[nodeCellID_[total_nodes_counter]];
        vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
        U_dPointCells += dU;
    }
    reduce(U_dPointCells, minOp<vector>()); // take only normal values of U
    if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
    {
        U_dPointCells = vector(10, 0, 0);
        Info << "OpenFOAM cell Not found" << endl;
        // Info << "ring: " << ring << endl;
        // Info << "node: " << total_nodes_counter << endl;
        // Info << "radius: " << radius << endl;
    }
    return U_dPointCells;
}

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
