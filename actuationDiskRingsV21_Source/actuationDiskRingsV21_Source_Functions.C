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

tensor getNodeTransformTensor(
    vector Bi,
    vector diskPoint,
    vector uniDiskDir
){
    // calculate the blade vector
    vector bladeDir = Bi - diskPoint;
    vector bladeUniDir;
    if (mag(bladeDir) == 0){
        bladeUniDir = vector(0,0,1);
    } else {
        bladeUniDir = bladeDir / mag(bladeDir);
    }

    // calculate the tangential vector
    vector F_tita_dir = vector(
        uniDiskDir[1] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[1],
        -1 * (uniDiskDir[0] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[0]),
        uniDiskDir[0] * bladeUniDir[1] - uniDiskDir[1] * bladeUniDir[0]
    ); 
    F_tita_dir = F_tita_dir / mag(F_tita_dir);
    
    // calculate the tensor transformation of coordinates
    vector vector_n = -1 * uniDiskDir;
    vector vector_t = F_tita_dir;
    vector vector_r = bladeUniDir;
    tensor transform(
        vector_n[0], vector_t[0], vector_r[0],
        vector_n[1], vector_t[1], vector_r[1],
        vector_n[2], vector_t[2], vector_r[2]); 

    return transform;
}

float getNodePhiAngle(
    tensor transform,
    vector U_dPointCells,
    float radius,
    float omega
){
    vector U_dPointCells_ntr = inv(transform) & U_dPointCells; 
    float U_n = -1 * U_dPointCells_ntr[0];
    float U_t = -1 * U_dPointCells_ntr[1];

    float phi = M_PI; // so that center node has sen(phi)=1
    if (radius != 0)
    {
        if (omega > 0)
        {
            phi = Foam::atan(U_n / (U_t - radius * omega));
        }
    }
    
    return phi;
}

float S0Function(
    float Ct,
    float Ct_rated
){
    float S0;
    if (Ct < Ct_rated) {
        S0 = 0.08 * pow((Ct_rated - Ct)/Ct_rated, 3);
    } else if (Ct >= Ct_rated) {
        S0 = 0.05 * (Ct_rated - Ct)/Ct_rated;
    }
    return S0;
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

float tipFactorFunction(
    int tipFactorType,
    float x,
    float lambda,
    float phi
){
    float Nb = 3;
    float F;
    if (tipFactorType == 0){
        // tip factor off
        F = 1;
    }
    else if (tipFactorType == 1){
        // tip factor by Shen 2005
        scalar c1 = 0.125;
        scalar c2 = 27;
        scalar c3 = 0.1;

        if (phi == 0) {
            F = 1; // f = 1/0 = infty, then exp(-infty)=0, then acos(0) = pi/2, then F=1
        } else if (x > 0) {
            scalar g = std::exp(-c1 * (Nb * lambda - c2)) + c3;
            scalar f = (Nb / 2) * (1 - x) / (x * std::sin(phi));
            if (f > 0)
            {
                if (((std::exp(-g * f)) > -1) and ((std::exp(-g * f)) < 1))
                {
                    F = (2 / (M_PI)) * std::acos(std::exp(-g * f));
                }
            }
        } else {
            F = 1;
        }
    }
    else if (tipFactorType == 2){
        // tip factor by Prandtl
        F = (2/M_PI) * std::acos( std::exp( - (Nb/2) * std::sqrt(1 + pow(lambda,2)) * (1 - x) ) );
    }
    else {
        Info << "tipFactor type not valid" << endl;
    }
    
    if ( F < 0 ) {
        Info << "F dio menor a 0!!" << endl;
        F = 0;
    }
    return F;
}

float rootFactorFunction(
    int rootFactorType,
    float x,
    float rootDistance,
    float phi
){
    float g;
    float Nb = 3;
    if (rootFactorType == 0){
        g = 1;
    }
    else if (rootFactorType == 1){
        // root factor by Glauert
        scalar f_tip;
        if (
            (phi == 0) || 
            (x == 0)
         ){
            f_tip = 100; // high value, because if phi=0 f_tip-> infty
        } else {
            f_tip = (Nb / 2) * (1 - x) / (x * std::sin(phi));
        }

        if (x <= rootDistance)
        {
            g = 0;
        }
        else if (phi == 0)
        {
            g = 1;
        }
        else if ((f_tip > 0) and (x < 0.5))
        {
            scalar f = (Nb / 2) * (x - rootDistance) / (x * std::sin(phi));
            if ((std::exp(-f) > -1) and (std::exp(-f) < 1))
            {
                g = (2 / (M_PI)) * std::acos(std::exp(-f));
            }
        }
        else
        {
            g = 1;
        }
    }
    else if (rootFactorType == 2){
        // root factor by Sorensen 2020
        float a = 2.335;
        int b = 4;
        if (x <= rootDistance)
        {
            g = 0;
        }
        else if (x < 0.5)
        {
            g = 1 - std::exp( - a * pow((x/rootDistance), b) );
        }
        else
        {
            g = 1;
        }
    }
    else {
        Info << "rootFactor type not valid" << endl;
    }

    return g;
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
