/*---------------------------------------------------------------------------*\
Author: Dimas Alejandro Barile
CSC - CONICET, Buenos Aires, 2022

-----------------------------------------------------------------------------
OpenFOAM 2.4

Actuation Disk

COMPLETAR

-----------------------------------------------------------------------------

Class
    Foam::fv::actuationDiskRingsV21_Source

Description
    Actuation disk source

    Constant values for momentum source for actuation disk
    \f[
    T = 0.5 \rho A U_{inf}^2 Ct
    \f]

    where:
    \vartable
        A   	| disk area
        U_inf 	| upstream velocity
        Ct	| thrust coeficient
    \endvartable

    \heading Source usage

    Example usage:
    \verbatim
    actuationDiskRingsV21_SourceCoeffs
    {
        fieldNames      (U);        	// names of fields to apply source
        diskDir         (-1 0 0);   	// disk direction
        Ct              0.5;        	// thrust coefficient
    Cp		0.4;	   	// power coefficient
        diskArea        5.0;        	// disk area
        upstreamPoint   (0 0 0);    	// upstream point
    diskPoint   	(1 0 0);    	// upstream point
    }
    \endverbatim


SourceFiles
    actuationDiskRingsV21_Source.C
    actuationDiskRingsV21_SourceTemplates.C

\*---------------------------------------------------------------------------*/

#include "actuationDiskRingsV21_Source.H"
#include "volFields.H"
#include <functional>
#include <math.h>
#include "fvc.H"
#include "fvCFD.H"
#include <map>
// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

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


template <class RhoFieldType>
// void Foam::fv::actuationDiskRingsV21_Source::addactuationDiskRings_AxialInertialResistance(
scalar Foam::fv::actuationDiskRingsV21_Source::addactuationDiskRings_AxialInertialResistance(
    vectorField &Usource,
    const labelList &cells,
    const scalarField &Vcells,
    const RhoFieldType &rho,
    const vectorField &U,
    // vectorField &force) const
    vectorField &force,
    scalar &UrefPrevious) const
{

    //------------------//Information of the actuator line---------------------------

    // Info<< "---Actuator disk N: " << this->name() << endl;
    // int tamano = CtList_.size();
    //

    //------------------List definitions------------------------------------------

    // Scalar acumulation
    scalar Pcells = 0.0;      // Power in each line cell
    scalar Tcells = 0.0;      // Thrust in each line cell
    scalar Ftitacells = 0.0;  // Tangential force in each line cell
    scalar TorqueSects = 0.0; // Torque in each seaction

    std::map<int, float> weightCells; 
    std::map<int, float> weightCells_U;
    std::map<int, float> weightCellsAD;
    // Vector acumulation
    std::map<int, float> weightCellsADCenter;

    vector U_dCenterCells_orient = vector(0, 0, 0);                   // Ud for the center cells of the sphere for orientation calculation
    vector U_dCenterCells = vector(0, 0, 0);                          // Ud for the center cells
    vector U_dCells = vector(0, 0, 0);                                // Ud in each cell
    vector U_infCells = vector(0, 0, 0);                              // Uinf in each cell
    const Field<vector> zoneCellCentres(mesh().cellCentres(), cells); // line cell centers

    // dinamic list (with append) of the AL blade points //Information
    DynamicList<scalar> aziList;    // azimut
    DynamicList<scalar> posList;    // position
    DynamicList<scalar> chordList;  // chord
    DynamicList<scalar> betaList;   // beta angle (pitch + twist)
    DynamicList<scalar> attackList; // angle of attack velocity
    DynamicList<scalar> UnList;     // local normal velocity
    DynamicList<scalar> UtList;     // local tangential velocity
    DynamicList<scalar> UrelList;   // local relative velocity
    DynamicList<scalar> FnList;     // normal force
    DynamicList<scalar> FtList;     // tangential force

    // calculate the radius
    scalar maxR = sqrt(diskArea_ / M_PI);

    //----- YAW ROTATION OPTION  -------------------------------------------------------------
    vector uniDiskDir = vector(0, 0, 0);
    scalar yawRad;
    
    if (yaw_ == 360)
    {
        // volume of the center cells of the sphere
        scalar Vcenter_orient = 0.0;
        forAll(cells, c)
        {
            if (mag(mesh().cellCentres()[cells[c]] - diskPoint_) < (centerRatio_ * maxR))
            {
                Vcenter_orient += Vcells[cells[c]];
            }
        }
        reduce(Vcenter_orient, sumOp<scalar>());

        // Ud vector for the center cells of the sphere
        // aca podria entrar en un loop que ya tenga el filtro hecho a partir del loop anterior, guardando las celdas que pasan el filtro a otra 'listq'
        forAll(cells, c)
        {
            if (mag(mesh().cellCentres()[cells[c]] - diskPoint_) < (centerRatio_ * maxR))
            {
                U_dCenterCells_orient += U[cells[c]] * (Vcells[cells[c]] / Vcenter_orient);
            }
        }
        reduce(U_dCenterCells_orient, sumOp<vector>());

        // no orientation z component
        U_dCenterCells_orient[2] = 0;

        // for self orientation depending on the velocity on the disc
        uniDiskDir = U_dCenterCells_orient / mag(U_dCenterCells_orient);

        // angle between the self orientation vector and the inlet vector
        vector Vinlet = diskDir_ / mag(diskDir_);
        float cosAlpha = (Vinlet[0] * uniDiskDir[0] + Vinlet[1] * uniDiskDir[1]) / (mag(Vinlet) * mag(uniDiskDir));
        float Alpha = acos(cosAlpha) * 180.0 / M_PI;
    }
    else
    {
        yawRad = yaw_ * 2 * M_PI / 360;
        vector diskYawed = vector(diskDir_[0] * cos(yawRad) - diskDir_[1] * sin(yawRad), diskDir_[0] * sin(yawRad) + diskDir_[1] * cos(yawRad), diskDir_[2]);
        uniDiskDir = diskYawed;
    }

    //------------------separation between disc centers in the AD------------------------

    scalar minSep = 0;
    minSep = cellSize_;
    Info << "min horizontal Cell center separation for calculation (deltaX): " << minSep << endl;

    //----create list of cells that belong to the disc-----------------------------------

    // smearing factor recomended when only know the mesh resolution
    scalar E = 2 * minSep;
    scalar En = 2 * minSep;
    scalar Er = 2 * minSep;
    scalar Et = 2 * minSep;
    // for point velocity smaller sphere
    scalar E_U = 0.5 * minSep;

    DynamicList<label> cellsDisc;

    forAll(cells, c)
    {
        // distance from AD plane
        scalar d = mag (uniDiskDir[0]*(mesh().cellCentres()[cells[c]][0]- diskPoint_[0])  
        + uniDiskDir[1]*(mesh().cellCentres()[cells[c]][1]- diskPoint_[1])
        + uniDiskDir[2]*(mesh().cellCentres()[cells[c]][2]- diskPoint_[2]))
        / sqrt(pow (uniDiskDir[0], 2) + pow (uniDiskDir[1], 2) + pow (uniDiskDir[2], 2));

        // distance of the cell center from sphere
        scalar dSphere = mag(mesh().cellCentres()[cells[c]] - diskPoint_);

        if ((dSphere <= maxR * 1.15) and (d <= 3 * E))
        {
            cellsDisc.append(cells[c]);
            // Info <<"appending to cellsDisc" <<endl;
        }
    }
    // Info <<"total cells in the procesor: "<<cellsDisc.size()<<endl;

    //------------------wight of the AD cells depending distance from plane and sphere-----

    scalar centerRatio = centerRatio_; // in the future it has to be determined by user in fvOptions --> it will unlock to make a study of it's influence in the calculation of Ud
    scalar V_AD = 0.0;
    scalar Vcenter = 0.0;

    // loop over the AD cells for weight and volume calculation
    forAll(cellsDisc, c)
    {
        //---weight calculation
        // distance from AD plane
        scalar d = mag(uniDiskDir[0] * (mesh().cellCentres()[cellsDisc[c]][0] - diskPoint_[0]) + uniDiskDir[1] * (mesh().cellCentres()[cellsDisc[c]][1] - diskPoint_[1]) + uniDiskDir[2] * (mesh().cellCentres()[cellsDisc[c]][2] - diskPoint_[2])) / sqrt(pow(uniDiskDir[0], 2) + pow(uniDiskDir[1], 2) + pow(uniDiskDir[2], 2));

        // distance from sphere
        scalar dSphere = mag(mesh().cellCentres()[cellsDisc[c]] - diskPoint_);

        scalar weightADplane = 0;
        scalar weightSphereCenter = 0;
        scalar weightSphereAD = 0;

        // weight calculation
        if (UdCellsMethod_ == 0)
        {
            // sum directly values in cells in AD
            weightADplane = 1;
            weightSphereCenter = 1;
            weightSphereAD = 1;
        }
        else if (UdCellsMethod_ == 1)
        {
            // weight with Gaussian in distance to AD plane
            weightADplane = (1 / (E * sqrt(M_PI))) * exp(-1 * pow((d / E), 2));
            weightSphereCenter = 1;
            weightSphereAD = 1;
        }
        else if (UdCellsMethod_ == 2)
        {
            // weight with Gaussian in distance to AD plane + distance to center of AD
            weightADplane = (1 / (E * sqrt(M_PI))) * exp(-1 * pow((d / E), 2));
            weightSphereCenter = (1 / ((centerRatio * maxR) * sqrt(M_PI))) * exp(-1 * pow((dSphere / (centerRatio * maxR)), 2));
            weightSphereAD = (1 / (maxR * sqrt(M_PI))) * exp(-1 * pow((dSphere / maxR), 2));
        }

        // evaluate weight
        if (dSphere <= maxR)
        {
            weightCellsAD[cellsDisc[c]] = weightADplane * weightSphereAD;
            if (dSphere <=  (centerRatio * maxR))
            {
                weightCellsADCenter[cellsDisc[c]] = weightADplane * weightSphereCenter;
            }
            else
            {
                weightCellsADCenter[cellsDisc[c]] = 0;
            }
        }
        else
        {
            weightCellsAD[cellsDisc[c]] = 0;
            weightCellsADCenter[cellsDisc[c]] = 0;
        }
        //---

        //---volume of the center cells weighted
        Vcenter += Vcells[cellsDisc[c]] * weightCellsADCenter[cellsDisc[c]];
        //---

        //---volume of the AD cells weighted
        V_AD += Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]];
        //---
    }
    reduce(Vcenter, sumOp<scalar>());
    reduce(V_AD, sumOp<scalar>());

    Info << "UdCenterToggle = " << UdCenterToggle_ << endl;
    forAll(cellsDisc, c)
    {
        U_dCenterCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsADCenter[cellsDisc[c]]) / Vcenter);

        if (UdCenterToggle_ == 0)
        {
            U_dCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]]) / V_AD);
        }
        else if (UdCenterToggle_ == 1)
        {
            U_dCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsADCenter[cellsDisc[c]]) / Vcenter);
        }
    }
    reduce(U_dCells, sumOp<vector>());
    reduce(U_dCenterCells, sumOp<vector>());

    // scalar cosU_dCenterCells= (U_dCenterCells[0]*uniDiskDir[0]+U_dCenterCells[1]*uniDiskDir[2])/(mag(U_dCenterCells)*mag(uniDiskDir));
    // scalar U_dCenterCellsYaw=mag(U_dCenterCells)*cosU_dCenterCells;

    scalar cosU_dCells = (U_dCells[0] * uniDiskDir[0] + U_dCells[1] * uniDiskDir[1]) / (mag(U_dCells) * mag(uniDiskDir));
    scalar U_dCellsYaw = mag(U_dCells) * cosU_dCells;

    Info << "U_dCenterCells not yawed: " << mag(U_dCenterCells) << endl;
    Info << "U_dCells not yawed: " << mag(U_dCells) << endl;
    Info << "U_dCells yawed: " << U_dCellsYaw << endl;
    Info << "U_dCells: " << U_dCells << endl;

    int pos1 = 0;
    int pos2 = 0;
    int pos = 0;
    scalar UrefYaw = 0;
    scalar omega = 0;
    scalar pitch = 0;
    scalar Ct = 0;
    scalar Cp = 0;
    float a1;
    float a2;
    float q0;
    DynamicList<float> UavgRings;

    // get the time in this step
    scalar t = mesh().time().value();
    const volVectorField &U_ = mesh().lookupObject<volVectorField>("U");
    volTensorField gradU = fvc::grad(U_);

    if (UrefCalculationMethod_ == 1) {
        //----- Numeric AD - Navarro Diaz 2019 -------------------------------------------------------------------
        //----- Find positions in table 1 (UdAvgList) for interpolation using UdCells ----------------------------
        float difference = GREAT;
        if (mag(U_dCells) < UdAvgList_[0]) // if the U_d in the disc is out the table
        {
            pos = 0; // lower position
        }
        else
        {
            // search the value in the table
            for (int i = 0; i < (UrefList_.size() - 1); i = i + 1)
            {
                if ((fabs(mag(U_dCells) - UdAvgList_[i]) < difference) && ((mag(U_dCells) - UdAvgList_[i]) >= 0))
                {
                    difference = fabs(mag(U_dCells) - UdAvgList_[i]);
                    pos = i; // the position of the lower value
                }
            }
        }

        pos1 = pos;
        pos2 = pos + 1;
        //----- End of find positions in table 1 (UdAvgList) for interpolation using UdCells ----------------------

        //---- Interpolate UrefYaw, omega, pitch, Ct and Cp-----------------------------------------------
        UrefYaw = ((mag(U_dCells) - UdAvgList_[pos]) * ((UrefList_[pos + 1]) - (UrefList_[pos])) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + (UrefList_[pos]);
        Info << "UrefYaw(m/s) interpolated: " << UrefYaw << endl;

        omega = ((mag(U_dCells) - UdAvgList_[pos]) * (omegaList_[pos + 1] - omegaList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + omegaList_[pos];
        Info << "omega(rad/s) interpolated: " << omega << endl;

        pitch = ((UrefYaw - UrefList_[pos]) * ((pitchList_[pos + 1]) - (pitchList_[pos])) / (UrefList_[pos + 1] - UrefList_[pos])) + (pitchList_[pos]);
        Info << "pitch(deg) interpolated " << pitch << endl;
        pitch = pitch * 2 * M_PI / 360;

        Ct = ((mag(U_dCells) - UdAvgList_[pos]) * (CtList_[pos + 1] - CtList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + CtList_[pos];
        Info << "Ct interpolated " << Ct << endl;

        Cp = ((mag(U_dCells) - UdAvgList_[pos]) * (CpList_[pos + 1] - CpList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + CpList_[pos];
        Info << "Cp interpolated " << Cp << endl;
        //---- End of interpolate UrefYaw, omega, pitch, Ct and Cp-----------------------------------------------
        
        float lambda = omega * maxR_ / UrefYaw;
        Info << "lambda = " << lambda << endl;
        a1 = a1Function( rootDistance_, lambda);
        a2 = a2Function( rootDistance_, lambda);
        q0 = ( sqrt(16 * pow(lambda,2) * pow(a2,2) + 8 * a1 * Ct) - 4 * lambda * a2 ) / ( 4 * a1 );
        Info << "q0 = " << q0 << endl; 
    } 
    else if (UrefCalculationMethod_ == 2) {
        // Analytic from Sorensen 2020
        // Ct = from UdAvg and UrefPrevious
        // Ct = 4 * (mag(U_dCells) / UrefPrevious) * ( 1 - (mag(U_dCells) / UrefPrevious) );
        Ct = Ct_;
        Info << "Ct = " << Ct << endl;
        // omega no se de donde sale
        // Info << "maxR_ = " << maxR_ << endl;
        Info << "lambda_ = " << lambda_ << endl;
        omega = lambda_ * UrefPrevious / maxR_;
        Info << "omega = " << omega << endl;
        pitch = 0;

        // q0 calculation considering uniform inflow
        // factores a1 y a2 para calcular q0 del metodo analitico
        a1 = a1Function( rootDistance_, lambda_);
        a2 = a2Function( rootDistance_, lambda_);
        q0 = ( sqrt(16 * pow(lambda_,2) * pow(a2,2) + 8 * a1 * Ct) - 4 * lambda_ * a2 ) / ( 4 * a1 );

        Info << "q0 = " << q0 << endl; 

        // con Udi saco U0i, y el promedio de U0i es Uref
        // porque se vuelvo a sacr Uref con el Ct que calcule voy a obtener de nuevo UrefPrevious porque son la misma expresion
        // Cp calculation
        scalar tita_r = 0;
        scalar tita_n_Rad = 0;
        scalar rMed_r = 0;
        scalar total_nodes_counter = 0;
        total_nodes_counter = 0;
        for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
        {
            // Info << "inside ring loop. ring " << ring << endl;
            tita_r = ringTitaList_[ring];
            rMed_r = ringrMedList_[ring];

            // inicializar velocidad promedio del anillo
            float UavgRing = 0;
            for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
            {
                tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                scalar x_node = 0;
                scalar y_node = 0;
                scalar z_node = 0;

                if (ring != numberRings_)
                {
                    x_node = -1 * rMed_r * sin(tita_n_Rad) * sin(yawRad);
                    y_node = rMed_r * sin(tita_n_Rad) * cos(yawRad);
                    z_node = rMed_r * cos(tita_n_Rad);
                }

                x_node += diskPoint_[0];
                y_node += diskPoint_[1];
                z_node += diskPoint_[2];

                vector Bi = vector(x_node, y_node, z_node);
                float radius = mag(diskPoint_ - Bi);

                vector U_dPointCells = vector(1000, 1000, 1000);
                // Info << "ring " << ring << endl;
                // Info << "nodeIterator " << nodeIterator << endl;
                // Info << "nodeCellID_[total_nodes_counter] = " << nodeCellID_[total_nodes_counter] << endl; 
                // Info << "total_nodes_counter = " << total_nodes_counter << endl; 

                if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
                {
                    if (ring == numberRings_)
                    {
                        U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
                    }
                    else
                    {
                        U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                    }

                    // if (gradInterpolation_ == 1)
                    if ( 
                        (gradInterpolation_ == 1) and
                        (ring != numberRings_)
                    ){
                        vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
                        vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                        U_dPointCells += dU;
                    }
                }

                // if (ring == numberRings_)
                // {
                //     if (nodeCellID_[nodesNumber_-1] != -1) //if the closer cell is in this procesor
                //     {
                //         U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
                //     }
                // }
                // else
                // {
                //     if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
                //     {
                //         U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                //     }
                // }
                // if ( 
                //     (gradInterpolation_ == 1) and
                //     (ring != numberRings_)
                // ){
                //     vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
                //     vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                //     U_dPointCells += dU;
                // }

                reduce(U_dPointCells, minOp<vector>()); // take only normal values of U

                if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
                {
                    U_dPointCells = vector(10, 0, 0);
                    Info << "OpenFOAM cell Not found" << endl;
                    Info << "ring: " << ring << endl;
                    Info << "node: " << total_nodes_counter << endl;
                    Info << "radius: " << radius << endl;
                }
                // Info << "U_dPointCells = " << U_dPointCells << endl;

                total_nodes_counter += 1;

                // vector U_dPointCells = vector(1000, 1000, 1000);
                // Info << "ring " << ring << endl;
                // Info << "nodeIterator " << nodeIterator << endl;
                // Info << "nodeCellID_[total_nodes_counter] = " << nodeCellID_[total_nodes_counter] << endl; 
                // Info << "total_nodes_counter = " << total_nodes_counter << endl;
                // if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
                // {
                //     // medir la velocidad en el nodo
                //     if (ring == numberRings_)
                //     {
                //         U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
                //     }
                //     else
                //     {
                //         U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                //     }

                //     // if (gradInterpolation_ == 1)
                //     if ( 
                //         (gradInterpolation_ == 1) and
                //         (ring != numberRings_)
                //     ){
                //         vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
                //         vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                //         U_dPointCells += dU;
                //     }
                // }
                // reduce(U_dPointCells, sumOp<vector>());
                // if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
                // {
                //     U_dPointCells = vector(10, 0, 0);
                //     Info << "OpenFOAM cell Not found" << endl;
                //     Info << "ring: " << ring << endl;
                //     Info << "node: " << total_nodes_counter << endl;
                //     Info << "radius: " << radius << endl;
                // }
                // sumar a la velocidad del anillo
                UavgRing += mag(U_dPointCells); 
                // reduce(U_dPointCells, minOp<vector>()); // take only normal values of U
                // total_nodes_counter += 1;
            }
            // reduce(UavgRing, sumOp<float>());
            // Info << "ring = " <<  ring << endl;
            // Info << "ringNodesList_[ring] = " <<  ringNodesList_[ring] << endl;
            UavgRing /= ringNodesList_[ring];
            // Info << "UavgRing = " << UavgRing << endl;
            // dividir la suma de velocidades del anillo por la cantidad de nodos en el anillo
            // apendar a la lista de velocidades promedio de anillo
            UavgRings.append(UavgRing);
        }
        // reduce(UavgRings, sumOp<float>());
        // Info << "ok calculating UavgRings" << endl;
        float resolution = 100;
        float integral = 0; 
        int ringCounter = 0;
        float gPrev;
        float FPrev;
        float funcPrev;
        float g;
        float F;
        float func;
        float Udx;
        float UdxPrev;
        float U0x;
        float U0xPrev;
        float r;
        for (float x = 1/resolution; x < 1; x += 1/resolution) {
            Udx = UavgRings[ringCounter]; 
            // Info << "Udx = " << Udx << endl; 
            // Info << "(1 + sqrt(1 - Ct)) = " << (1 + sqrt(1 - Ct)) << endl;
            U0x = 2 * Udx / (1 + sqrt(1 - Ct));

            if (ringrMedList_[ringCounter] - x < 1/resolution) {
                UdxPrev = UavgRings[ringCounter - 1]; 
            } else {
                UdxPrev = UavgRings[ringCounter]; 
            }
            U0xPrev = 2 * UdxPrev / (1 + sqrt(1 - Ct));

            g = gFunction(x, rootDistance_);
            F = FFunction(x, lambda_);

            gPrev = gFunction(x - 1/resolution, rootDistance_);
            FPrev = FFunction(x - 1/resolution, lambda_);

            func = (Udx/U0x) * g * F * x;
            funcPrev = (UdxPrev/U0xPrev) * gPrev * FPrev * (x - 1/resolution);
             
            integral += (1/resolution) * (func + funcPrev)/2;

            r = x * maxR_;
            // Info << "r = " << r << "; ringCounter = " << ringCounter << endl;

            if (
                (x + 1/resolution > ringrMedList_[ringCounter] / maxR_) and
                ( ringCounter < numberRings_ )
            ) {
                ringCounter += 1;
            }
        }
        // reduce(integral, sumOp<float>());
        // Info << "ok integral for Cp" << endl;
        Cp = 4 * lambda_ * q0 * integral;
        Info << "Cp = " << Cp << endl;
        
        // if (mesh().time().value() == 1)
        // {
        //     float UrefAvg = 0;
        //     for (int ring = 0; ring <= (numberRings_); ring = ring + 1) {
        //         UrefAvg += 2 * UavgRings[ring] / (1 + sqrt(1 - Ct));
        //     }
        //     // Info << "numberRings_ = " << numberRings_ << endl;
        //     // reduce(UrefAvg, sumOp<float>());
        //     UrefAvg /= numberRings_;
        //     UrefYaw = UrefAvg;
        // }
        // else
        // {
        //     UrefYaw = 2 * mag(U_dCells) / (1 + sqrt(1 - Ct)); 
        // }
        // float UrefAvg = 0;
        // for (int ring = 0; ring <= (numberRings_); ring = ring + 1) {
        //     UrefAvg += 2 * UavgRings[ring] / (1 + sqrt(1 - Ct));
        // }
        // Info << "numberRings_ = " << numberRings_ << endl;
        // reduce(UrefAvg, sumOp<float>());
        // UrefAvg /= numberRings_;
        // UrefYaw = UrefAvg;
        UrefYaw = 2 * mag(U_dCellsYaw) / (1 + sqrt(1 - Ct));
    } 
    else if (UrefCalculationMethod_ == 3) {
        // generalization of Sorensen 2020, considering discretization in cells to avoid assuming uniform inflow
        // Ct = 4 * (mag(U_dCells) / UrefPrevious) * ( 1 - (mag(U_dCells) / UrefPrevious) );
        Ct = Ct_;
        Info << "Ct = " << Ct << endl;
        // omega no se de donde sale
        // Info << "maxR_ = " << maxR_ << endl;
        Info << "lambda_ = " << lambda_ << endl;
        omega = lambda_ * UrefPrevious / maxR_;
        Info << "omega = " << omega << endl;
        pitch = 0;
        
        // q0 calculation for generalized inflow
        // loop through nodes to calculate the generalized a1 and a2
        a1 = 0;
        a2 = 0;
        float gNode;
        float FNode;
        float xNode;
        float UrefNode;
        float nodeArea;
        scalar tita_r = 0;
        scalar tita_n_Rad = 0;
        scalar rMed_r = 0;
        scalar total_nodes_counter = 0;
        for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
        {
            tita_r = ringTitaList_[ring];
            rMed_r = ringrMedList_[ring];
            for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
            {
                tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                scalar x_node = 0;
                scalar y_node = 0;
                scalar z_node = 0;

                if (ring != numberRings_)
                {
                    x_node = -1 * rMed_r * sin(tita_n_Rad) * sin(yawRad);
                    y_node = rMed_r * sin(tita_n_Rad) * cos(yawRad);
                    z_node = rMed_r * cos(tita_n_Rad);
                }

                x_node += diskPoint_[0];
                y_node += diskPoint_[1];
                z_node += diskPoint_[2];

                vector Bi = vector(x_node, y_node, z_node);
                float radius = mag(diskPoint_ - Bi);

                vector U_dPointCells = vector(1000, 1000, 1000);
                if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
                {
                    if (ring == numberRings_)
                    {
                        U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
                    }
                    else
                    {
                        U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                    }

                    // if (gradInterpolation_ == 1)
                    if ( 
                        (gradInterpolation_ == 1) and
                        (ring != numberRings_)
                    ){
                        vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
                        vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                        U_dPointCells += dU;
                    }
                }
                reduce(U_dPointCells, minOp<vector>()); // take only normal values of U
                if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
                {
                    U_dPointCells = vector(10, 0, 0);
                    Info << "OpenFOAM cell Not found" << endl;
                    Info << "ring: " << ring << endl;
                    Info << "node: " << total_nodes_counter << endl;
                    Info << "radius: " << radius << endl;
                }
                total_nodes_counter += 1;
                
                xNode = rMed_r / maxR_;
                gNode = gFunction(xNode, rootDistance_);  
                FNode = FFunction(xNode, lambda_);  
                nodeArea = ringAreaList_[ring];
                UrefNode = 2 * mag(U_dPointCells) / (1 + sqrt(1 - Ct));
                a1 += pow(UrefNode, 2) * gNode * FNode * nodeArea;
                if (xNode > 0) {
                    a2 += pow(UrefNode, 2) * pow( gNode * FNode / xNode , 2) * nodeArea /2;
                }
            }
        }
        Info << "a1 = " << a1 << endl; 
        Info << "a2 = " << a2 << endl; 
        float density = 1.225;
        UrefYaw = 2 * mag(U_dCellsYaw) / (1 + sqrt(1 - Ct));

        // float Uref_test_node = 2 * mag(U_dCellsYaw) / (1 + sqrt(1 - Ct));
        float thrust = Ct * 0.5 * M_PI * density * pow(maxR_, 2) * pow(UrefYaw, 2);
        // q0 = ( - lambda_ * a1 + sqrt( pow(lambda_ * a1,2) + 4 * a2 * ( Ct * M_PI  * pow(maxR_, 2) / 2 ) ) )/(2 * a2);
        q0 = ( - lambda_ * a1 + sqrt( pow(lambda_ * a1,2) + 4 * a2 * ( thrust / density ) ) )/(2 * a2);

        Info << "q0 = " << q0 << endl; 

        total_nodes_counter = 0;
        for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
        {
            // Info << "inside ring loop. ring " << ring << endl;
            tita_r = ringTitaList_[ring];
            rMed_r = ringrMedList_[ring];

            // inicializar velocidad promedio del anillo
            float UavgRing = 0;
            for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
            {
                tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                scalar x_node = 0;
                scalar y_node = 0;
                scalar z_node = 0;

                if (ring != numberRings_)
                {
                    x_node = -1 * rMed_r * sin(tita_n_Rad) * sin(yawRad);
                    y_node = rMed_r * sin(tita_n_Rad) * cos(yawRad);
                    z_node = rMed_r * cos(tita_n_Rad);
                }

                x_node += diskPoint_[0];
                y_node += diskPoint_[1];
                z_node += diskPoint_[2];

                vector Bi = vector(x_node, y_node, z_node);
                float radius = mag(diskPoint_ - Bi);

                vector U_dPointCells = vector(1000, 1000, 1000);

                if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
                {
                    if (ring == numberRings_)
                    {
                        U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
                    }
                    else
                    {
                        U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                    }

                    // if (gradInterpolation_ == 1)
                    if ( 
                        (gradInterpolation_ == 1) and
                        (ring != numberRings_)
                    ){
                        vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
                        vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                        U_dPointCells += dU;
                    }
                }

                reduce(U_dPointCells, minOp<vector>()); // take only normal values of U

                if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
                {
                    U_dPointCells = vector(10, 0, 0);
                    Info << "OpenFOAM cell Not found" << endl;
                    Info << "ring: " << ring << endl;
                    Info << "node: " << total_nodes_counter << endl;
                    Info << "radius: " << radius << endl;
                }

                total_nodes_counter += 1;
                UavgRing += mag(U_dPointCells); 
            }
            UavgRing /= ringNodesList_[ring];
            // dividir la suma de velocidades del anillo por la cantidad de nodos en el anillo
            // apendar a la lista de velocidades promedio de anillo
            UavgRings.append(UavgRing);
        }
        float resolution = 100;
        float integral = 0; 
        int ringCounter = 0;
        float gPrev;
        float FPrev;
        float funcPrev;
        float g;
        float F;
        float func;
        float Udx;
        float UdxPrev;
        float U0x;
        float U0xPrev;
        float r;
        for (float x = 1/resolution; x < 1; x += 1/resolution) {
            Udx = UavgRings[ringCounter]; 
            U0x = 2 * Udx / (1 + sqrt(1 - Ct));

            if (ringrMedList_[ringCounter] - x < 1/resolution) {
                UdxPrev = UavgRings[ringCounter - 1]; 
            } else {
                UdxPrev = UavgRings[ringCounter]; 
            }
            U0xPrev = 2 * UdxPrev / (1 + sqrt(1 - Ct));

            g = gFunction(x, rootDistance_);
            F = FFunction(x, lambda_);

            gPrev = gFunction(x - 1/resolution, rootDistance_);
            FPrev = FFunction(x - 1/resolution, lambda_);

            func = (Udx/U0x) * g * F * x;
            funcPrev = (UdxPrev/U0xPrev) * gPrev * FPrev * (x - 1/resolution);
             
            integral += (1/resolution) * (func + funcPrev)/2;

            r = x * maxR_;
            // Info << "r = " << r << "; ringCounter = " << ringCounter << endl;

            if (
                (x + 1/resolution > ringrMedList_[ringCounter] / maxR_) and
                ( ringCounter < numberRings_ )
            ) {
                ringCounter += 1;
            }
        }
        Cp = 4 * lambda_ * q0 * integral;
        Info << "Cp = " << Cp << endl;
        
    }
    Info << "UrefYaw = " << UrefYaw << endl;

    //----- Calculate Thrust, Power and Torque------------------------------------------
    float density_ = 1.225;
    scalar upRho = 1; // Thrust
    float T = 0.50 * upRho * diskArea_ * pow(UrefYaw, 2) * Ct;
    Info << "Thrust fixed (with density): " << T * density_ << endl;

    // Power
    float P = 0.5 * upRho * diskArea_ * pow(UrefYaw, 3) * Cp;
    Info << "Power fixed [MW] (with density) " << P * density_ * 0.001 << endl;

    // Torque
    float Torque = 0;
    if (omega > 0)
    {
        Torque = P / omega;
    }
    //----- End of calculate Thrust, Power and Torque------------------------------------------

    //----- Initialization of parameters of loops of force calculation ------------------------------------------

    // For calculating node positions
    scalar tita_r = 0;              // angle between nodes in a certain ring
    scalar tita_n_Rad = 0;          // angle for the position of a certain node - rad
    scalar rMed_r = 0;              // radius for a certain ring
    scalar total_nodes_counter = 0; // Nodes counter

    // initalize WT vectors
    vector vector_n = vector(-1, 0, 0);
    vector vector_t = vector(0, -1, 0);
    vector vector_r = vector(0, 0, 1);

    // initalize blade and cell points
    vector Pi_ntr = vector(0, 0, 0);
    vector Bi_ntr = vector(0, 0, 0);

    // velocities for the profile BEM
    vector U_dPointCells_ntr = vector(0, 0, 0);
    scalar U_n = 0;
    scalar U_t = 0;
    scalar phi = 0;
    scalar radius = 0;
    scalar nblades_ = 3;

    // initalize distances from the blade point
    scalar dn = 0.0;
    scalar dt = 0.0;
    scalar dr = 0.0;

    // initialize blade section forces variables
    scalar F_n_Bi = 0;
    scalar F_tita_Bi = 0;
    vector F_tita_dir = vector(0, 0, 0);
    scalar V_point_F = 0.0;

    // for the interpolation of table 2
    float U_inf1 = GREAT;
    float U_inf2 = GREAT;
    float U_inf_point = GREAT;
    float fn1 = GREAT;
    float fn2 = GREAT;
    float fn_point = GREAT;
    float ft1 = GREAT;
    float ft2 = GREAT;
    float ft_point = GREAT;
    scalar rtable = 0;
    scalar posr = 0;
    scalar lastPos = 0;
    float x_point; // x_point = radius_point / maxR

    // count how many radius sections are
    int nR_ = 1;
    while (rList_[nR_] - rList_[nR_ - 1] > 0)
    {
        nR_ = nR_ + 1;
    }

    // angle of the rotation in this time
    scalar clockwiseSign = 1;
    bool clockwise = true;
    if (clockwise)
    {
        clockwiseSign = -1;
    }

    // Velocity field pointer
    // const volVectorField &U_ = mesh().lookupObject<volVectorField>("U");
    // volTensorField gradU = fvc::grad(U_);

    Info << "Starting loop through nodes" << endl;
    Info << " " << endl;
    total_nodes_counter = 0;
    //----- End of initialization of parameters of loops of force calculation ------------------------------------------

    //--- LOOP OVER RINGS FOR FORCE CALCULATION AND DISTTIBUTION -----------------------------------------------------------------
    // loop through rings and nodes for calculating forces and distributing them
    for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
    {
        tita_r = ringTitaList_[ring];
        rMed_r = ringrMedList_[ring];

        // for each node
        // -- LOOP OVER NODES IN RING FOR FORCE CALCULATION AND DISTRBUTION -------------------------------------------------
        for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
        {
            //----- Calculate node position ------------------------------------------
            tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

            scalar x_node = 0;
            scalar y_node = 0;
            scalar z_node = 0;

            // position of the node considering disk center = (0,0,0)
            if (ring != numberRings_)
            {
                x_node = -1 * rMed_r * sin(tita_n_Rad) * sin(yawRad);
                y_node = rMed_r * sin(tita_n_Rad) * cos(yawRad);
                z_node = rMed_r * cos(tita_n_Rad);
            }

            // move to turbine position
            x_node += diskPoint_[0];
            y_node += diskPoint_[1];
            z_node += diskPoint_[2];
            //----- End of calculate node position ------------------------------------------

            //----- Calculate radial and tangential vectors ------------------------------------------
            // blade vector
            vector bladeUniDir = vector(0, 0, 1); // we force this vector for the center node
            vector bladeDir = vector(0, 0, 1); // initialization 
            // if (ring == numberRings_)
            // {
            //     bladeUniDir = vector(0, 0, 1); // we force this vector for the center node
            // }
            // else
            if (ring != numberRings_)
            {
                bladeDir = vector(x_node - diskPoint_[0], y_node - diskPoint_[1], z_node - diskPoint_[2]);
                bladeUniDir = bladeDir / mag(bladeDir);
            }

            // calculate the tangential vector
            F_tita_dir = vector(uniDiskDir[1] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[1],
                                -1 * (uniDiskDir[0] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[0]),
                                uniDiskDir[0] * bladeUniDir[1] - uniDiskDir[1] * bladeUniDir[0]);

            F_tita_dir = F_tita_dir / mag(F_tita_dir);
            //----- End of calculate radial and tangential vectors ------------------------------------------

            //----- Calculate tensor to transfor cartesian coordinates to cylindrical coordinates -------------------------
            vector_n = -1 * uniDiskDir;
            vector_t = F_tita_dir;
            vector_r = bladeUniDir;

            tensor Transform(vector_n[0], vector_t[0], vector_r[0],
                             vector_n[1], vector_t[1], vector_r[1],
                             vector_n[2], vector_t[2], vector_r[2]);
            //----- End of calculate tensor to transfor cartesian coordinates to cylindrical coordinates -------------------------

            vector Bi = vector(x_node, y_node, z_node);

            radius = mag(diskPoint_ - Bi);

            // save for the output
            posList.append(radius);

            // change of coordinate system
            Bi_ntr = inv(Transform) & Bi;

            //----- Calculate velocity in node ------------------------------------------
            // calculate velocity in node
            vector U_dPointCells = vector(1000, 1000, 1000);
            // Info << "ring " << ring << endl;
            // Info << "nodeIterator " << nodeIterator << endl;
            // Info << "nodeCellID_[total_nodes_counter] = " << nodeCellID_[total_nodes_counter] << endl; 
            // Info << "total_nodes_counter = " << total_nodes_counter << endl; 
            
            if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
            {
                if (ring == numberRings_)
                {
                    U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
                }
                else
                {
                    U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                }

                // if (gradInterpolation_ == 1)
                if ( 
                    (gradInterpolation_ == 1) and
                    (ring != numberRings_)
                ){
                    vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
                    vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                    U_dPointCells += dU;
                }
            }

            // Info << "U_dPointCells = " << U_dPointCells << endl;

            // if (ring == numberRings_)
            // {
            //     if (nodeCellID_[nodesNumber_-1] != -1) //if the closer cell is in this procesor
            //     {
            //         U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
            //     }
            // }
            // else
            // {
            //     if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
            //     {
            //         U_dPointCells = U[nodeCellID_[total_nodes_counter]];
            //     }
            // }

            // // if (gradInterpolation_ == 1)
            // if ( 
            //     (gradInterpolation_ == 1) and
            //     (ring != numberRings_)
            // ){
            //     vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
            //     vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
            //     U_dPointCells += dU;
            // }
            
            reduce(U_dPointCells, minOp<vector>()); // take only normal values of U

            if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
            {
                U_dPointCells = vector(10, 0, 0);
                Info << "OpenFOAM cell Not found" << endl;
                Info << "ring: " << ring << endl;
                Info << "node: " << total_nodes_counter << endl;
                Info << "radius: " << radius << endl;
            }
            // Info << "U_dPointCells = " << U_dPointCells << endl;

            total_nodes_counter += 1;

            // change of coordinate system
            U_dPointCells_ntr = inv(Transform) & U_dPointCells;

            // velocities in the profile coordinates
            U_n = -1 * U_dPointCells_ntr[0];
            U_t = -1 * U_dPointCells_ntr[1];
            //----- End of calculate velocity in node ------------------------------------------

            // volume of the point cells weighted, for forces
            V_point_F = 0;

            //---- Loop over cells to weight in relation to distance to current node ---------------------------------- 
            forAll(cellsDisc, c)
            {
                // change of coordinate system
                // Pi_ntr = coordinate of cell in cylindrical cooridnate system
                Pi_ntr = inv(Transform) & mesh().cellCentres()[cellsDisc[c]];

                // calculate the distances from node to cell center in blade coordinate system
                dn = mag(Pi_ntr[0] - Bi_ntr[0]);
                dt = mag(Pi_ntr[1] - Bi_ntr[1]);
                dr = mag(Pi_ntr[2] - Bi_ntr[2]);

                //---for forces
                // calculate the wight dependig on the distances from the sectional point
                // very similar to convolution distribution form Mikkelsen PhD thesis in 2003
                // difference that here it decomposes in the directions.That part inside the exponential isn't exactly the same as in Mikkelsen proposal
                if (forceDistributionMethod_==0) { // no force distribution. Forces are applied directly in cell that contains the node

                }
                else if (forceDistributionMethod_==1) { // force distribution previously implented in this code
                    weightCells[cellsDisc[c]] = (1 / (En * Et * Er * pow(sqrt(M_PI), 3))) *
                                            exp(-1 * (pow(dn / En, 2) + pow(dt / Et, 2) + pow(dr / Er, 2)));
                }
                else if (forceDistributionMethod_==2) { // force distribution as proposed in Mikkelsen 2003
                    weightCells[cellsDisc[c]] = (1 / (pow(E, 3) * pow(sqrt(M_PI), 3))) *
                                            exp(-1 * ( (sqrt(pow(dn,2) + pow(dr,2) + pow(dt,2))) / (E) ) );
                }
                // Info << "weight of cell calculated" << endl;

                // distance of the cell center from sphere
                scalar dSphere = mag(mesh().cellCentres()[cellsDisc[c]] - diskPoint_);

                // calculate weight if its inside the sphere (cut exedent outside the disc)
                if (dSphere <= maxR * 1.15)
                {
                    weightCells[cellsDisc[c]] = weightCells[cellsDisc[c]] * 1;
                }
                else
                {
                    weightCells[cellsDisc[c]] = weightCells[cellsDisc[c]] * 0;
                }

                V_point_F += Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]];
            }
            // volume weighted
            reduce(V_point_F, sumOp<scalar>());
            // reduce(weightCells, sumOp<std::map<int, float>>());
            // Info << "ok calculating weight for cells for force distribution" << endl;
            //---- End of loop over cells to weight in relation to distance to current node ---------------------------------- 

            //---- force calculation with Numeric Actuator (Navarro Diaz 2019) -----------------------------------------------
            float nodeArea;
            if (forceCalculationMethod_ == 1) {
                //---- Calculate Uinf, fn and ft for each position of table 2 for interpolation ------------------------------ 
                posr = posrList_[ring + 1];
                rtable = rNodeList_[ring + 1];

                if (ring == numberRings_)
                {
                    posr = 0;
                    rtable = 0;
                }

                //---- Calculate Uinf, fn and ft for position 1 ---------------------------------- 
                scalar pos1table2 = posInTableUref2(
                    pos1,
                    Uref2List_,
                    UrefList_,
                    rList_,
                    rtable,
                    U_dPointCells_ntr,
                    UdiList_);

                U_inf1 = ((mag(U_dPointCells_ntr) - UdiList_[pos1table2]) * ((UinfList_[pos1table2 + nR_]) - (UinfList_[pos1table2])) / (UdiList_[pos1table2 + nR_] - UdiList_[pos1table2])) + (UinfList_[pos1table2]);
                fn1 = ((mag(U_dPointCells_ntr) - UdiList_[pos1table2]) * ((fnList_[pos1table2 + nR_]) - (fnList_[pos1table2])) / (UdiList_[pos1table2 + nR_] - UdiList_[pos1table2])) + (fnList_[pos1table2]);
                ft1 = ((mag(U_dPointCells_ntr) - UdiList_[pos1table2]) * ((ftList_[pos1table2 + nR_]) - (ftList_[pos1table2])) / (UdiList_[pos1table2 + nR_] - UdiList_[pos1table2])) + (ftList_[pos1table2]);
                //---- End of calculate Uinf, fn and ft for position 1 ---------------------------------- 

                //---- Calculate Uinf, fn and ft for position 2 ---------------------------------- 
                scalar pos2table2 = posInTableUref2(
                    pos2,
                    Uref2List_,
                    UrefList_,
                    rList_,
                    rtable,
                    U_dPointCells_ntr,
                    UdiList_);

                U_inf2 = ((mag(U_dPointCells_ntr) - UdiList_[pos2table2]) * ((UinfList_[pos2table2 + nR_]) - (UinfList_[pos2table2])) / (UdiList_[pos2table2 + nR_] - UdiList_[pos2table2])) + (UinfList_[pos2table2]);
                fn2 = ((mag(U_dPointCells_ntr) - UdiList_[pos2table2]) * ((fnList_[pos2table2 + nR_]) - (fnList_[pos2table2])) / (UdiList_[pos2table2 + nR_] - UdiList_[pos2table2])) + (fnList_[pos2table2]);
                ft2 = ((mag(U_dPointCells_ntr) - UdiList_[pos2table2]) * ((ftList_[pos2table2 + nR_]) - (ftList_[pos2table2])) / (UdiList_[pos2table2 + nR_] - UdiList_[pos2table2])) + (ftList_[pos2table2]);
                //---- End of calculate Uinf, fn and ft for position 2 ---------------------------------- 

                //---- Interpolate Uinf, fn and ft from values of position 1 and 2 ---------------------------------- 
                U_inf_point = ((UrefYaw - UrefList_[pos1]) * ((U_inf2) - (U_inf1)) / (UrefList_[pos2] - UrefList_[pos1])) + (U_inf1);
                fn_point = ((UrefYaw - UrefList_[pos1]) * ((fn2) - (fn1)) / (UrefList_[pos2] - UrefList_[pos1])) + (fn1);
                ft_point = ((UrefYaw - UrefList_[pos1]) * ((ft2) - (ft1)) / (UrefList_[pos2] - UrefList_[pos1])) + (ft1);
                //---- End of interpolate Uinf, fn and ft from values of position 1 and 2 ---------------------------------- 
                //---- End of calculate Uinf, fn and ft for each position of table 2 for interpolation ---------------------------------- 
                F_n_Bi = (fn_point * ringThickness_ * 3) / ringNodesList_[ring];
                F_n_Bi /= density_;
                F_tita_Bi = (ft_point * ringThickness_ * 3) / ringNodesList_[ring];
                F_tita_Bi /= density_;
            }
            else if (forceCalculationMethod_ == 2)
            // Analytical Actuator Disk by Sorensen 2020
            {
                // Info << "calculation method = Analytical by Sorensen 2020" << endl;
                U_inf_point = 2 * mag(U_dPointCells) / (1 + sqrt( 1 - Ct ));
                // Info << "U_inf_point = " << U_inf_point << endl;
                x_point = radius / maxR_;
                if (x_point == 0) {
                    fn_point = 0;
                    ft_point = 0;
                } else {
                    // Info << "g = " << gFunction(x_point, rootDistance_) << endl;
                    // Info << "F = " << FFunction(x_point, lambda_) << endl;
                    // Info << "U_inf_point = " << U_inf_point << endl;
                    // Info << "density = " << density_ << endl;
                    // Info << "q0 = " << q0 << endl;

                    fn_point = density_ * pow(U_inf_point, 2) * q0 * (gFunction(x_point, rootDistance_) * FFunction(x_point, lambda_) / x_point) * (lambda_ * x_point + q0 * (gFunction(x_point, rootDistance_) * FFunction(x_point, lambda_) / ( 2 * x_point)));
                    ft_point = density_ * pow(U_inf_point, 2) * q0 * (gFunction(x_point, rootDistance_) * FFunction(x_point, lambda_) / x_point) * (mag(U_dPointCells) / U_inf_point);

                    nodeArea = ringAreaList_[ring];
                    // fn_point = nodeArea * density_ * pow(U_inf_point, 2) * q0 * (gFunction(x_point, rootDistance_) * FFunction(x_point, lambda_) / x_point) * (lambda_ * x_point + q0 * (gFunction(x_point, rootDistance_) * FFunction(x_point, lambda_) / ( 2 * x_point))); 
                    // ft_point = nodeArea * density_ * pow(U_inf_point, 2) * q0 * (gFunction(x_point, rootDistance_) * FFunction(x_point, lambda_) / x_point) * (mag(U_dPointCells) / U_inf_point);

                    fn_point *= nodeArea;
                    ft_point *= nodeArea;

                    // Info << "fn = " << fn_point << endl;
                    // Info << "ft = " << ft_point << endl;
                }

                F_n_Bi = fn_point / density_; // without density to the solver
                F_tita_Bi = ft_point / density_; // without density to the solver
            }
            // Info << "fn = " << fn_point << endl;
            // Info << "ft = " << ft_point << endl;

            //---- Calculate total force of node from fn and ft ---------------------------------- 
            // Tangential and axial forces
            // divide the force depending on the amount of artificial blades (but because the force in table
            // is refered to 1 blade, i need to multiply by 3)

            // save for the output
            FnList.append(F_n_Bi * density_);
            FtList.append(F_tita_Bi * density_);
            UnList.append(U_n);
            UtList.append(U_t);

            // global from all the blades
            TorqueSects += mag(diskPoint_ - Bi) * F_tita_Bi * density_;
            Pcells += mag(diskPoint_ - Bi) * F_tita_Bi * density_ * omega;
            //---- End of calculate total force of node from fn and ft ---------------------------------- 

            // --- Force distribution in cells -----------------------------------------------------------------------------
            // loop over all the cells to apply source
            forAll(cellsDisc, c)
            {
                // apply the thrust force in the cell, volume weighed
                Usource[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_n_Bi) * uniDiskDir;

                // apply the tangential force in the cell, volume weighed
                Usource[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_tita_Bi) * F_tita_dir;

                // save the thrust force in the force field for paraview
                force[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_n_Bi) * uniDiskDir * -1 * density_ / Vcells[cellsDisc[c]];

                // save the tangential force in the force field for praview
                force[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_tita_Bi) * F_tita_dir * -1 * density_ / Vcells[cellsDisc[c]];

                // save the results of each cell, volume weighed
                Tcells += F_n_Bi * density_ * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);
                U_infCells += uniDiskDir * UrefYaw * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);

                Ftitacells += F_tita_Bi * density_ * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);

            } // close loop cells in sectional point
            // --- End of force distribution in cells -----------------------------------------------------------------------

        } // close loop in ring
        //---- End of loop through nodes in ring to calculate and distribute forces ----------------

    } // close loop rings
    //---- End of loop through rings to calculate and distribute forces ----------------

    //---colecting data from all the procesors
    reduce(Tcells, sumOp<scalar>());
    reduce(Ftitacells, sumOp<scalar>());
    reduce(U_infCells, sumOp<vector>());
    scalar U_infCellsMag = mag(U_infCells) / nodesNumber_;

    //-----print results------------------------------------------
    //---global turbine outputs
    if (Pstream::myProcNo() == 0) // if Im in the master proccesor
    {
        //"Turbine, time [s], Uref(calc) [m/s], Cp(calc), Ct(calc), Power_calc(Uref,Cp) [W], Thrust_calc(Uref,Ct) [N] , Center Ud [m/s], Average Ud [m/s], Uinf(disc cells) [m/s], Power(disc cells) [W], Thrust(disc cells) [N], Torque(disc cells) [Nm]"
        (*outTurbines) << this->name() << "," << t << "," << UrefYaw << "," << Cp << "," << Ct << "," << P * density_ << "," << T * density_ << "," << mag(U_dCenterCells) << "," << mag(U_dCells) << "," << U_infCellsMag << "," << Pcells << "," << Tcells << "," << TorqueSects << std::endl;
    }

    //---nodes outputs

    if (Pstream::myProcNo() == 0 and t > 0) // if Im in the master proccesor and from an specific time
    {
        // current time
        (*outRings) << "time [s] = " << t << std::endl;

        // turbine
        (*outRings) << "turbine = " << this->name() << std::endl;

        // rotational veloicty (fixed)
        (*outRings) << "rotational velocity [rad/s] (given) = " << omega << std::endl;

        // number of rings
        (*outRings) << "Number of rings (including ring 0) = " << numberRings_ + 1 << std::endl;

        // Ring thickness
        (*outRings) << "Ring thickness = " << ringThickness_ << std::endl;

        // AD ring parameters
        (*outRings) << " " << std::endl;
        (*outRings) << "AD ring parameters" << std::endl;
        for (int ring = 0; ring <= numberRings_; ring = ring + 1)
        {
            // Save rings and node per ring
            (*outRings) << "ring: " << ring;
            (*outRings) << " - nodes: " << ringNodesList_[ring];
            (*outRings) << " - tita (deg): " << ringTitaList_[ring];
            (*outRings) << " - r: " << ringrMedList_[ring];
            (*outRings) << " - node area: " << ringAreaList_[ring] << std::endl;
        }

        // AD ring parameters
        (*outRings) << " " << std::endl;
        (*outRings) << "Rings" << std::endl;

        // header
        (*outRings) << "node, ring, radius [m], node area [m2], normal velocity [m/s], tangential velocity [m/s], normal force [N], tangential force [N] " << std::endl;

        int i = 0;
        // nodes
        // for each ring
        for (int ring = 0; ring <= numberRings_; ring = ring + 1)
        {
            // for each node
            for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
            {

                (*outRings) << i + 1 << ", " << ring << ", " << posList[i] << ", " << ringAreaList_[ring] << ", "
                            << UnList[i] << ", " << UtList[i] << ", " << FnList[i] << ", " << FtList[i] << std::endl;

                i += 1;
            }
        }

        (*outRings) << std::endl;
    }

    return UrefYaw;
}
// ************************************************************************* //
