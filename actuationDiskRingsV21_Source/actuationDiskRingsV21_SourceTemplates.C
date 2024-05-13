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
    // scalar &UrefPrevious) const
    scalar &CtPrevious) const
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

    vectorField cellCentres = mesh().cellCentres();

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
        }
    }

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

    if (UdCorrection_ == 1) {
        // scalar axialInductionFactor = 0.5 - 0.5 * sqrt(1 - Ct);
        scalar axialInductionFactor = 0.5 - 0.5 * sqrt(1 - CtPrevious);
        float Ct_modified = CtPrevious / pow(1 - axialInductionFactor,2);
        scalar correctionScalar = pow( 1 + ( Ct_modified * E )/(2 * sqrt(2 * M_PI) * maxR_) ,-1);
        U_dCellsYaw *= correctionScalar;
    }

    Info << "U_dCenterCells not yawed: " << mag(U_dCenterCells) << endl;
    Info << "U_dCells not yawed: " << mag(U_dCells) << endl;
    Info << "U_dCells yawed: " << U_dCellsYaw << endl;
    Info << "U_dCells: " << U_dCells << endl;

    // --- save velocity in each node ----
    scalar t = mesh().time().value();
    const volVectorField &U_ = mesh().lookupObject<volVectorField>("U");
    volTensorField gradU = fvc::grad(U_);

    DynamicList<vector> U_dNodes;
    scalar tita_r = 0;              // angle between nodes in a certain ring
    scalar tita_n_Rad = 0;          // angle for the position of a certain node - rad
    scalar rMed_r = 0;              // radius for a certain ring
    scalar total_nodes_counter = 0; // Nodes counter
    for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
    {
        tita_r = ringTitaList_[ring];
        rMed_r = ringrMedList_[ring];

        // for each node
        // -- LOOP OVER NODES IN RING FOR FORCE CALCULATION AND DISTRBUTION -------------------------------------------------
        for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
        {
            //----- Calculate node position ------------------------------------------
            // Info << "node: " << total_nodes_counter << endl;
            tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

            vector Bi = getNodePosition(tita_n_Rad, rMed_r, yawRad, diskPoint_, ring, numberRings_);

            vector U_dPointCells = vector(1000, 1000, 1000);
            if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
            {
                // medir la velocidad en el nodo
                if (ring == numberRings_)
                {
                    U_dPointCells =  U[nodeCellID_[nodesNumber_-1]];
                }
                else
                {
                    U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                }
                if ( 
                    (gradInterpolation_ == 1) and
                    (ring != numberRings_)
                ){
                    vector dx = Bi - cellCentres[nodeCellID_[total_nodes_counter]];
                    vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                    U_dPointCells += dU;
                }
            }
            reduce(U_dPointCells, minOp<vector>()); // take only normal values of U
            if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
            {
                Info << "mag(U_dPointCells) > 100 for node " << total_nodes_counter << endl;
                if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
                {
                    vector dx = Bi - cellCentres[nodeCellID_[total_nodes_counter]];
                    vector dU = dx & gradU[nodeCellID_[total_nodes_counter]];
                    Pout << "U_dPointCells = " << U_dPointCells << endl;
                    Pout << "dx = " << dx << endl;
                    Pout << "dU = " << dU << endl;
                }
            }

            U_dNodes.append(U_dPointCells);

            total_nodes_counter += 1;
        }
    }
    // ----------------------------------

    float fn;
    float ft;
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
    float a3;
    float a4;
    float a5;
    float q0;
    float S0;
    float gNode;
    float FNode;
    DynamicList<float> UavgRings;
    DynamicList<float> phiavgRings;
    scalar scale_factor_n;
    scalar scale_factor_t;

    if (
        (ADmodel_ == 0) or
        (ADmodel_ == 5)
    ){
        // combine uniform and eliptic ADs since both do scale factoring
        if ( ADmodel_ == 0 ) {
            Info << "" << endl;
            Info << "Uniform AD model starting..." << endl;
        } else if ( ADmodel_ == 5 ) {
            Info << "" << endl;
            Info << "Eliptic AD model starting..." << endl;
        }

        scalar cosUinfAngle = (diskDir_[0] * uniDiskDir[0] + diskDir_[1] * uniDiskDir[1]) / (mag(diskDir_) * mag(uniDiskDir));
        UrefYaw = Uref_ * cosUinfAngle;
        // upRho = 1;
        omega = omega_;
        // scalar lambda = maxR * omega / UrefYaw;
        Ct = Ct_;
        Cp = Cp_;

        Info << "UrefYaw = " << UrefYaw << endl;
        Info << "Ct = " << Ct << endl;
        Info << "Cp = " << Cp << endl;
        Info << "omega = " << omega << endl;

        // thrust, power and torque without density (to then pass fn and ft without density to solver)
        float thrust = 0.5 * diskArea_ * pow(UrefYaw, 2) * Ct;
        float power = 0.5 * diskArea_ * pow(UrefYaw, 3) * Cp;
        float torque = 0;
        if (omega > 0)
        {
            torque = power / omega;
        }

        // uniform normal and tangential surface forces
        if ( ADmodel_ == 0 ){
            fn = thrust / diskArea_;
            ft = (3 * torque) / (2 * M_PI * pow(maxR_, 3));
            Info << "fn mean = " << fn << endl;
            Info << "ft mean = " << ft << endl;
        }

        // sums for scale factors (used if root or tip factor are in use) 
        scalar sumF_n_Bi = 0;
        scalar sumF_n_Bixfactor = 0;
        scalar sumTorque_Bi = 0;
        scalar sumTorque_Bixfactor = 0;

        if (
            (rootFactor_ == 0) and
            (tipFactor_ == 0)
        ){
            // if no root and tip factor, the scale factor will be 1, so:
            sumF_n_Bi = 1;
            sumF_n_Bixfactor = 1;
            sumTorque_Bi = 1;
            sumTorque_Bixfactor = 1;
        } else {
            // if root or tip factor are in use, scale factor needs to be calculated
            scalar tita_r = 0;
            scalar tita_n_Rad = 0;
            scalar rMed_r = 0;
            total_nodes_counter = 0;
            for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
            {
                tita_r = ringTitaList_[ring];
                rMed_r = ringrMedList_[ring];
                for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
                {
                    tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                    vector Bi = getNodePosition(tita_n_Rad, rMed_r, yawRad, diskPoint_, ring, numberRings_);
                    float radius = mag(diskPoint_ - Bi);

                    vector U_dPointCells = U_dNodes[total_nodes_counter];
                    tensor transform = getNodeTransformTensor(
                        Bi,
                        diskPoint_,
                        uniDiskDir_
                    );
                    float phi = getNodePhiAngle(
                        transform,
                        U_dPointCells,
                        radius,
                        omega
                    );
                    total_nodes_counter += 1;

                    float xNode = radius / maxR_;
                    float nodeArea = ringAreaList_[ring];

                    gNode = rootFactorFunction(rootFactor_, xNode, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);
                    FNode = tipFactorFunction(tipFactor_, xNode, lambda_, phi, U_dPointCells, UrefYaw);

                    if ( ADmodel_ == 0 ){
                        sumF_n_Bi += fn * nodeArea;
                        sumF_n_Bixfactor += gNode * FNode * fn * nodeArea;
                        if (xNode > 0) {
                            sumTorque_Bi += ft * nodeArea;
                            sumTorque_Bixfactor += gNode * FNode * ft * nodeArea;
                        }
                    } else if ( ADmodel_ == 5 ) {
                        float density = 1.225;
                        fn = 0.75 * density * pow(UrefYaw,2) * Ct * sqrt(1 - xNode);
                        sumF_n_Bixfactor += gNode * fn * nodeArea;
                    }
                }
            }
        }
        if ( ADmodel_ == 0 ){
            scale_factor_n = sumF_n_Bi / sumF_n_Bixfactor;
            scale_factor_t = sumTorque_Bi / sumTorque_Bixfactor;
        } else if ( ADmodel_ == 5 ) {
            scale_factor_n = thrust / sumF_n_Bixfactor;
        }

    } else if (
        (ADmodel_ == 1) or  // Numeric AD Navarro Diaz
        (ADmodel_ == 6)     // Numeric AD van der Laan
    ){
        //----- Numeric AD - Navarro Diaz 2019 -------------------------------------------------------------------
        //----- Find positions in table 1 (UdAvgList) for interpolation using UdCells ----------------------------
        Info << "" << endl;
        if (ADmodel_ == 1) {
            Info << "Numeric AD (Navarro Diaz) model starting..." << endl;
        } else if (ADmodel_ == 6) {
            Info << "Numeric AD (van der Laan) model starting..." << endl;
        }

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

        if (ADmodel_ == 1) {
            // pitch calibration only in Navarro Diaz AD
            pitch = ((UrefYaw - UrefList_[pos]) * ((pitchList_[pos + 1]) - (pitchList_[pos])) / (UrefList_[pos + 1] - UrefList_[pos])) + (pitchList_[pos]);
            Info << "pitch(deg) interpolated " << pitch << endl;
            pitch = pitch * 2 * M_PI / 360;
        }

        Ct = ((mag(U_dCells) - UdAvgList_[pos]) * (CtList_[pos + 1] - CtList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + CtList_[pos];
        Info << "Ct interpolated " << Ct << endl;

        Cp = ((mag(U_dCells) - UdAvgList_[pos]) * (CpList_[pos + 1] - CpList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + CpList_[pos];
        Info << "Cp interpolated " << Cp << endl;
        //---- End of interpolate UrefYaw, omega, pitch, Ct and Cp-----------------------------------------------
        
        // Calculate parameters used in analyitic AD for comparison
        float lambda = omega * maxR_ / UrefYaw;
        Info << "lambda = " << lambda << endl;
    } 
    else if (
        (ADmodel_ == 2) or
        (ADmodel_ == 3) or
        (ADmodel_ == 4)
    ) {
        Info << "" << endl;
        Info << "Analytic AD model starting..." << endl;

        // Initialize UrefYaw with UrefYaw from previous iteration
        // UrefYaw = UrefPrevious;

        // Initialize UrefYaw with Ct from previous iteration
        UrefYaw = 2 * mag(U_dCellsYaw) / (1 + sqrt(1 - CtPrevious));

        scalar UrefYawOld = UrefYaw;
        scalar UrefYawOld2 = UrefYawOld;
        float UrefDifference = GREAT;

        while (UrefDifference > 0.001 ) {
            // Get position in power curve table from UrefYaw
            float difference = GREAT;
            if (UrefYaw < UrefPowerCurveList_[0]) // if the U_d in the disc is out the table
            {
                pos = 0; // lower position
            }
            else
            {
                // search the value in the table
                for (int i = 0; i < (UrefPowerCurveList_.size() - 1); i = i + 1)
                {
                    if ((fabs(UrefYaw - UrefPowerCurveList_[i]) < difference) && ((UrefYaw - UrefPowerCurveList_[i]) >= 0))
                    {
                        difference = fabs(UrefYaw - UrefPowerCurveList_[i]);
                        pos = i; // the position of the lower value
                    }
                }
            }

            // Interpolate Ct from power curve table from UrefYaw
            Ct = ((UrefYaw - UrefPowerCurveList_[pos]) * ((CtPowerCurveList_[pos + 1]) - (CtPowerCurveList_ [pos])) / (UrefPowerCurveList_[pos + 1] - UrefPowerCurveList_[pos])) + (CtPowerCurveList_[pos]);

            // Calculate UrefYaw from interpolated Ct and formula used in Sorensen 2020
            UrefYaw = 2 * mag(U_dCellsYaw) / (1 + sqrt(1 - Ct));

            UrefDifference = fabs( UrefYaw - UrefYawOld );

            // Oscilating when calculating Ct and Uref, setting Uref as middle of oscilation
            if (UrefYaw == UrefYawOld2) {
                UrefYaw += (UrefYawOld - UrefYaw) * 0.5;
                Ct = 4 * (mag(U_dCells) / UrefYaw) * ( 1 - (mag(U_dCells) / UrefYaw) );
                UrefDifference = 0;
            }
            UrefYawOld2 = UrefYawOld;
            UrefYawOld = UrefYaw;
        }

        // Clipping of UrefYaw variation from power curve interpolation
        // (didn't prove to lower convergence iterations)
        // if (fabs( UrefYaw / UrefPrevious - 1 ) < 0.000001) {
        //     Info << "Using previous Uref, because difference < 0.001% for better convergence" << endl;
        //     UrefYaw = UrefPrevious;
        //     Ct = 4 * (mag(U_dCells) / UrefYaw) * ( 1 - (mag(U_dCells) / UrefYaw) );
        // }
        
        Info << "UrefYaw = " << UrefYaw << endl;
        Info << "Ct = " << Ct << endl;
        Info << "lambda_ = " << lambda_ << endl;
        omega = lambda_ * UrefYaw / maxR_;
        Info << "omega = " << omega << endl;
        pitch = 0;

        if (ADmodel_ == 2) {
            // q0 calculation considering uniform inflow
            // Cp calculation
            scalar tita_r = 0;
            scalar tita_n_Rad = 0;
            scalar rMed_r = 0;
            scalar total_nodes_counter = 0;
            for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
            {
                tita_r = ringTitaList_[ring];
                rMed_r = ringrMedList_[ring];

                // inicializar velocidad promedio del anillo
                float UavgRing = 0;
                float phiavgRing = 0;
                for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
                {
                    tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                    vector Bi = getNodePosition(tita_n_Rad, rMed_r, yawRad, diskPoint_, ring, numberRings_);
                    float radius = mag(diskPoint_ - Bi);
                    vector U_dPointCells = U_dNodes[total_nodes_counter];
                    tensor transform = getNodeTransformTensor(
                        Bi,
                        diskPoint_,
                        uniDiskDir_
                    );
                    float phi = getNodePhiAngle(
                        transform,
                        U_dPointCells,
                        radius,
                        omega
                    );
                    total_nodes_counter += 1;

                    UavgRing += mag(U_dPointCells); 
                    phiavgRing += phi;
                }
                UavgRing /= ringNodesList_[ring];
                phiavgRing /= ringNodesList_[ring];
                // dividir la suma de velocidades del anillo por la cantidad de nodos en el anillo
                // apendar a la lista de velocidades promedio de anillo
                UavgRings.append(UavgRing);
                phiavgRings.append(phiavgRing);
            }

            float resolution = 100;
            float integral = 0; 
            a1 = 0;
            a2 = 0;
            int ringCounter = 0;
            float gPrev = 0;
            float FPrev = 0;
            float funcPrev_a1 = 0;
            float funcPrev_a2 = 0;
            float funcPrev = 0;
            float g;
            float F;
            float func_a1;
            float func_a2;
            float func;
            float Udx;
            float U0x;
            float phi;
            for (float x = 1/resolution; x <= 1; x += 1/resolution) {
                Udx = UavgRings[ringCounter]; 
                phi = phiavgRings[ringCounter]; 
                U0x = 2 * Udx / (1 + sqrt(1 - Ct));
                vector UdxVector = vector(Udx,0,0);

                g = rootFactorFunction(rootFactor_, x, rootDistance_, phi, lambda_, UdxVector, UrefYaw);
                F = tipFactorFunction(tipFactor_, x, lambda_, phi, UdxVector, UrefYaw);

                func_a1 = pow(g * F,2) / x;   
                a1 += (1/resolution) * ( func_a1 + funcPrev_a1 )/2;          

                func_a2 = g * F * x;   
                a2 += (1/resolution) * ( func_a2 + funcPrev_a2 )/2;          

                func = (Udx/U0x) * g * F * x;
                 
                integral += (1/resolution) * (func + funcPrev)/2;

                if (
                    (x + 1/resolution > ringrMedList_[ringCounter] / maxR_) and
                    ( ringCounter < numberRings_ )
                ) {
                    ringCounter += 1;
                }

                funcPrev_a1 = func_a1;
                funcPrev_a2 = func_a2;
                funcPrev = func;
            }
            q0 = ( sqrt(16 * pow(lambda_,2) * pow(a2,2) + 8 * a1 * Ct) - 4 * lambda_ * a2 ) / ( 4 * a1 );
            Cp = 4 * lambda_ * q0 * integral;
        }
        else if (ADmodel_ == 3) {
            // q0 calculation for generalized inflow
            // loop through nodes to calculate the generalized a1 and a2
            a1 = 0;
            a2 = 0;
            float Cp_sum = 0;
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

                    vector Bi = getNodePosition(tita_n_Rad, rMed_r, yawRad, diskPoint_, ring, numberRings_);
                    float radius = mag(diskPoint_ - Bi);

                    vector U_dPointCells = U_dNodes[total_nodes_counter];
                    tensor transform = getNodeTransformTensor(
                        Bi,
                        diskPoint_,
                        uniDiskDir_
                    );
                    float phi = getNodePhiAngle(
                        transform,
                        U_dPointCells,
                        radius,
                        omega
                    );

                    total_nodes_counter += 1;

                    xNode = radius / maxR_;
                    gNode = rootFactorFunction(rootFactor_, xNode, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);  
                    FNode = tipFactorFunction(tipFactor_, xNode, lambda_, phi, U_dPointCells, UrefYaw);  
                    nodeArea = ringAreaList_[ring];

                    a1 += pow(UrefYaw, 2) * gNode * FNode * nodeArea;
                    if (xNode > 0) {
                        a2 += 0.5 * pow(UrefYaw, 2) * pow( gNode * FNode / xNode , 2) * nodeArea;
                    }

                    Cp_sum += UrefYaw * mag(U_dPointCells) * gNode * FNode * nodeArea;
                }
            }
            Info << "a1 = " << a1 << endl; 
            Info << "a2 = " << a2 << endl; 
            float density = 1.225;

            float nominalThrust = 0.5 * density * diskArea_ * pow(UrefYaw,2);
            float thrust = Ct * nominalThrust;
            q0 = ( - lambda_ * a1 + sqrt( pow(lambda_ * a1,2) + 4 * a2 * ( thrust / density ) ) )/(2 * a2);

            // Cp calculation
            float nominalPower = 0.5 * density * diskArea_ * pow(UrefYaw,3);
            float realPower = omega * maxR_ * density * q0 * Cp_sum ;
            Cp = realPower / nominalPower;
        }
        else if (ADmodel_ == 4) {
            a1 = 0;
            a2 = 0;
            a3 = 0;
            a4 = 0;
            a5 = 0;
            float b0;
            float b1;
            float b2;
            float Cp_sum = 0;
            float density = 1.225;
            float xNode;
            float UrefNode;
            float nodeArea;
            scalar tita_r = 0;
            scalar tita_n_Rad = 0;
            scalar rMed_r = 0;

            // q0 calculation 
            scalar total_nodes_counter = 0;
            for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
            {
                tita_r = ringTitaList_[ring];
                rMed_r = ringrMedList_[ring];
                for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
                {
                    tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                    vector Bi = getNodePosition(tita_n_Rad, rMed_r, yawRad, diskPoint_, ring, numberRings_);
                    float radius = mag(diskPoint_ - Bi);

                    vector U_dPointCells = U_dNodes[total_nodes_counter];
                    tensor transform = getNodeTransformTensor(
                        Bi,
                        diskPoint_,
                        uniDiskDir_
                    );
                    float phi = getNodePhiAngle(
                        transform,
                        U_dPointCells,
                        radius,
                        omega
                    );

                    total_nodes_counter += 1;
                    
                    xNode = radius / maxR_;
                    gNode = rootFactorFunction(rootFactor_, xNode, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);  
                    FNode = tipFactorFunction(tipFactor_, xNode, lambda_, phi, U_dPointCells, UrefYaw);  
                    nodeArea = ringAreaList_[ring];

                    a1 += nodeArea * pow(UrefYaw,2) * gNode * FNode;
                    a2 += nodeArea * pow(UrefYaw,2) * gNode * FNode * pow(xNode,2);
                    if (xNode > 0) {
                        a3 += nodeArea * pow(UrefYaw * gNode * FNode / xNode, 2);
                    }
                    a4 += nodeArea * pow(UrefYaw * gNode * FNode, 2);
                    a5 += nodeArea * pow(UrefYaw * gNode * FNode * xNode, 2);
                }
            }

            S0 = S0Function(Ct, Ct_rated_);

            float nominalThrust = 0.5 * density * diskArea_ * pow(UrefYaw,2);
            float thrust = Ct * nominalThrust;

            b0 = - S0 * lambda_ * a2 + 0.5 * pow(S0,2) * a5 - thrust / density;
            b1 = lambda_ * a1 - S0 * a4;
            b2 = 0.5 * a3;

            q0 = (-b1 + sqrt(pow(b1,2) - 4 * b0 * b2)) / (2 * b2);
            // ---

            // Cp calculation
            total_nodes_counter = 0;
            for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
            {
                tita_r = ringTitaList_[ring];
                rMed_r = ringrMedList_[ring];
                for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
                {
                    tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                    vector Bi = getNodePosition(tita_n_Rad, rMed_r, yawRad, diskPoint_, ring, numberRings_);
                    float radius = mag(diskPoint_ - Bi);

                    vector U_dPointCells = U_dNodes[total_nodes_counter];
                    tensor transform = getNodeTransformTensor(
                        Bi,
                        diskPoint_,
                        uniDiskDir_
                    );
                    float phi = getNodePhiAngle(
                        transform,
                        U_dPointCells,
                        radius,
                        omega
                    );

                    total_nodes_counter += 1;
                    
                    xNode = radius / maxR_;
                    gNode = rootFactorFunction(rootFactor_, xNode, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);  
                    FNode = tipFactorFunction(tipFactor_, xNode, lambda_, phi, U_dPointCells, UrefYaw);  
                    nodeArea = ringAreaList_[ring];

                    if (xNode > 0) {
                        Cp_sum += nodeArea * UrefYaw * mag(U_dPointCells) * (q0 / xNode - S0 * xNode) * gNode * FNode;
                    }
                }
            }
            float nominalPower = 0.5 * density * diskArea_ * pow(UrefYaw,3);
            float realPower = omega * maxR_ * density * Cp_sum;
            Cp = realPower / nominalPower;

            Info << "S0 = " << S0 << endl; 
        }

        Info << "q0 = " << q0 << endl; 
        Info << "Cp = " << Cp << endl;
    } 

    //----- Calculate Thrust, Power and Torque------------------------------------------
    float density_ = 1.225;
    scalar upRho = 1; 

    // Thrust
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
    total_nodes_counter = 0; // Nodes counter

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
    scalar rtable = 0;
    float ft_point = GREAT;
    scalar posr = 0;
    scalar lastPos = 0;
    float x_point; // x_point = radius_point / maxR

    // count how many radius sections are
    int nR_ = 1;
    if (
        (ADmodel_ == 1) or
        (ADmodel_ == 6)
    ){
        while (rList_[nR_] - rList_[nR_ - 1] > 0)
        {
            nR_ = nR_ + 1;
        }
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
    //----- End of initialization of parameters of loops of force calculation ------------------------------------------

    //--- LOOP OVER RINGS FOR FORCE CALCULATION AND DISTTIBUTION -----------------------------------------------------------------
    // loop through rings and nodes for calculating forces and distributing them
    total_nodes_counter = 0;
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

            vector Bi = getNodePosition(tita_n_Rad, rMed_r, yawRad, diskPoint_, ring, numberRings_);
            scalar x_node = Bi[0];
            scalar y_node = Bi[1];
            scalar z_node = Bi[2];
            //----- End of calculate node position ------------------------------------------

            tensor transform = getNodeTransformTensor(
                Bi,
                diskPoint_,
                uniDiskDir_
            );

            //----- End of calculate tensor to transfor cartesian coordinates to cylindrical coordinates -------------------------

            radius = mag(diskPoint_ - Bi);

            vector U_dPointCells = U_dNodes[total_nodes_counter];
            float phi = getNodePhiAngle(
                transform,
                U_dPointCells,
                radius,
                omega
            );

            // save for the output
            posList.append(radius);

            // change of coordinate system
            Bi_ntr = inv(transform) & Bi;

            //----- Calculate velocity in node ------------------------------------------
            // calculate velocity in node
            total_nodes_counter += 1;

            // change of coordinate system
            U_dPointCells_ntr = inv(transform) & U_dPointCells;

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
                Pi_ntr = inv(transform) & mesh().cellCentres()[cellsDisc[c]];

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
                else if (forceDistributionMethod_==1) { 
                    // force distribution previously implented in this code
                    weightCells[cellsDisc[c]] = (1 / (En * Et * Er * pow(sqrt(M_PI), 3))) *
                                            exp(-1 * (pow(dn / En, 2) + pow(dt / Et, 2) + pow(dr / Er, 2)));
                }
                else if (forceDistributionMethod_==2) { 
                    // force distribution as proposed in Mikkelsen 2003
                    weightCells[cellsDisc[c]] = (1 / (pow(E, 3) * pow(sqrt(M_PI), 3))) *
                                            exp(-1 * pow( (sqrt(pow(dn,2) + pow(dr,2) + pow(dt,2))) / (E) ,2) );
                }
                else if (forceDistributionMethod_==3) { 
                    // force distribution using optimization from Martinez et al. 2017
                    En = 0.17 * averageChordLength_;
                    Er = pow(10,-2) * averageChordLength_;
                    Et = pow(10,-2) * averageChordLength_;
                    float s_n = - 0.36 * averageChordLength_;
                    weightCells[cellsDisc[c]] = (1 / (En * Et * Er * pow(sqrt(M_PI), 3))) *
                                            exp(-1 * (pow((dn - s_n) / En, 2) + pow(dt / Et, 2) + pow(dr / Er, 2)));
                }
                else if (forceDistributionMethod_==4) {  
                    // force distribution with delta function of Li 2022
                    float diskCellSize = 0.5 * 2 * maxR_ / cellSize_;
                    float d_cellCentre_node = mag(Pi_ntr - Bi_ntr);
                    float d_adim = d_cellCentre_node / diskCellSize;
                    if (d_adim <= 0.5) {
                        weightCells[cellsDisc[c]] = 3/8 + M_PI/32 - pow(d_adim,2)/4;
                    } else if (d_adim <= 1.5) {
                        weightCells[cellsDisc[c]] = 1/4 + (1 - d_adim) * sqrt( -2 + 8 * d_adim - 4 * pow(d_adim,2) ) / 8 - asin(std::sqrt(2) * (d_adim - 1)) / 8;
                    } else if (d_adim <= 2.5) {
                        weightCells[cellsDisc[c]] = 17/16 + M_PI/64 - 3 * d_adim / 4 - pow(d_adim,2) / 8 + (d_adim - 2) * sqrt( -14 + 16 * d_adim - 4 * pow(d_adim,2)) / 16 + asin( std::sqrt(2) * (d_adim - 2) ) / 16;
                    } else {
                        weightCells[cellsDisc[c]] = 0;
                    }
                }

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
            reduce(V_point_F, sumOp<scalar>());
            //---- End of loop over cells to weight in relation to distance to current node ---------------------------------- 

            //---- force calculation with Numeric Actuator (Navarro Diaz 2019) -----------------------------------------------
            float nodeArea;
            if (ADmodel_ == 0) {
                float xNode = radius / maxR_;
                float nodeArea = ringAreaList_[ring];

                // Info << "xNode = " << xNode << endl;
                // Info << "phi = " << phi << endl;
                gNode = rootFactorFunction(rootFactor_, xNode, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);
                FNode = tipFactorFunction(tipFactor_, xNode, lambda_, phi, U_dPointCells, UrefYaw);
                // Info << "gNode = " << gNode << endl;
                // Info << "FNode = " << FNode << endl;

                F_n_Bi = gNode * FNode * scale_factor_n * fn * nodeArea;
                F_tita_Bi = gNode * FNode * scale_factor_t * ft * nodeArea;

            } else if (
                (ADmodel_ == 1) or
                (ADmodel_ == 6)
            ) {
                //---- Calculate Uinf, fn and ft for each position of table 2 for interpolation ------------------------------ 
                posr = posrList_[ring + 1];
                rtable = rNodeList_[ring + 1];

                if (ring == numberRings_)
                {
                    posr = 0;
                    rtable = 0;
                }

                if (ADmodel_ == 1) {
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
                    F_n_Bi = (fn_point * ringThickness_ * 3) / ringNodesList_[ring];
                    F_tita_Bi = (ft_point * ringThickness_ * 3) / ringNodesList_[ring];
                    F_n_Bi /= density_;
                    F_tita_Bi /= density_;
                } else if (ADmodel_ == 6) {
                    scalar pos1table2 = posInTableForcesVanDerLaan(
                        pos1,
                        UrefList_,
                        Uref2List_,
                        rList_,
                        rtable);
                    fn1 = fnList_[pos1table2];
                    ft1 = ftList_[pos1table2];
                    
                    scalar pos2table2 = posInTableForcesVanDerLaan(
                        pos2,
                        UrefList_,
                        Uref2List_,
                        rList_,
                        rtable);
                    fn2 = fnList_[pos2table2];
                    ft2 = ftList_[pos2table2];

                    fn_point = ((UrefYaw - UrefList_[pos1]) * ((fn2) - (fn1)) / (UrefList_[pos2] - UrefList_[pos1])) + (fn1);
                    ft_point = ((UrefYaw - UrefList_[pos1]) * ((ft2) - (ft1)) / (UrefList_[pos2] - UrefList_[pos1])) + (ft1);

                    float nodeArea = ringAreaList_[ring];
                    F_n_Bi = 3 * fn_point * density_ * pow(mag(U_dCells),2) * Ct * pow(UrefYaw/mag(U_dCells),2) * nodeArea / (6 * T);
                    Info << "F_n_Bi = " << F_n_Bi << endl;
                    F_tita_Bi = (ft_point * 3 * omega)  * density_ * pow(mag(U_dCells),3) * (Cp * pow(UrefYaw/mag(U_dCells),3) / omega) * nodeArea / (6 * P);
                    Info << "F_tita_Bi = " << F_tita_Bi << endl;
                    F_n_Bi /= density_;
                    F_tita_Bi /= density_;
                }
            }
            else if (
                (ADmodel_ == 2) or
                (ADmodel_ == 3)
            // Analytical Actuator Disk by Sorensen 2020
            ) {
                U_inf_point = 2 * mag(U_dPointCells) / (1 + sqrt( 1 - Ct ));
                x_point = radius / maxR_;
                if (x_point == 0) {
                    fn_point = 0;
                    ft_point = 0;
                } else {
                    float g_point = rootFactorFunction(rootFactor_, x_point, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);  
                    float F_point = tipFactorFunction(tipFactor_, x_point, lambda_, phi, U_dPointCells, UrefYaw);  
                    
                    fn_point = density_ * pow(UrefYaw, 2) * q0 * (g_point * F_point / x_point) * (lambda_ * x_point + q0 * (g_point * F_point / ( 2 * x_point)));
                    ft_point = density_ * pow(UrefYaw, 2) * q0 * (g_point * F_point / x_point) * (mag(U_dPointCells) / UrefYaw);

                    nodeArea = ringAreaList_[ring];

                    fn_point *= nodeArea;
                    ft_point *= nodeArea;
                }

                F_n_Bi = fn_point / density_; // without density to the solver
                F_tita_Bi = ft_point / density_; // without density to the solver
            }
            else if (ADmodel_ == 4)
            // Generalized Analytical Actuator Disk by Sorensen 2023
            {
                U_inf_point = 2 * mag(U_dPointCells) / (1 + sqrt( 1 - Ct ));
                x_point = radius / maxR_;
                if (x_point == 0) {
                    fn_point = 0;
                    ft_point = 0;
                } else {
                    float g_point = rootFactorFunction(rootFactor_, x_point, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);  
                    float F_point = tipFactorFunction(tipFactor_, x_point, lambda_, phi, U_dPointCells, UrefYaw);  
                    
                    fn_point = density_ * pow(UrefYaw, 2) * (q0 / x_point - S0 * x_point) * g_point * F_point * (lambda_ * x_point + 0.5 * (q0 / x_point - S0 * x_point) * g_point * F_point );
                    ft_point = density_ * UrefYaw * mag(U_dPointCells) * (q0 / x_point - S0 * x_point) * g_point * F_point;

                    nodeArea = ringAreaList_[ring];

                    fn_point *= nodeArea;
                    ft_point *= nodeArea;
                }

                F_n_Bi = fn_point / density_; // without density to the solver
                F_tita_Bi = ft_point / density_; // without density to the solver
            }
            else if (ADmodel_ == 5)
            // Eliptic Actuator Disk by Sorensen 1992 
            {
                x_point = radius / maxR_;
                gNode = rootFactorFunction(rootFactor_, x_point, rootDistance_, phi, lambda_, U_dPointCells, UrefYaw);

                fn_point = 0.75 * density_ * pow(UrefYaw,2) * Ct * sqrt(1 - x_point);
                ft_point = 0;

                nodeArea = ringAreaList_[ring];

                fn_point *= scale_factor_n * gNode * nodeArea;
                ft_point *= nodeArea;

                F_n_Bi = fn_point / density_; // without density to the solver
                F_tita_Bi = ft_point / density_; // without density to the solver
            }

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

    // return UrefYaw;
    if (
        (ADmodel_ == 0) or
        (ADmodel_ == 5)
    ){
        return CtPrevious;
    } else {
        return Ct;
    }
}
// ************************************************************************* //
