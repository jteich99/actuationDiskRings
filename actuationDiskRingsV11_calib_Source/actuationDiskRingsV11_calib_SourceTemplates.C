/*---------------------------------------------------------------------------*\
Author: Dimas Alejandro Barile
CSC - CONICET, Buenos Aires, 2022

-----------------------------------------------------------------------------
OpenFOAM 2.4

Actuation Disk

Actuation Disk made to represent wind turbines with a certain Thrust and Power
Cofficient (Ct, Cp). Nodes are arranged in the form of rings. The normal force
is distributed in a uniform way across AD plane, while the tangential force is
distributed in a way that the torque produced is uniform across de AD plane.

-----------------------------------------------------------------------------

Class
    Foam::fv::actuationDiskRingsV11_calib_Source

Description
    Actuation disk source

    Constant values for momentum source for actuation disk

    \f[
    T = 0.5 \rho A U_{inf}^2 Ct
    \f]

    \f[
    P = 0.5 \rho A U_{inf}^3 Cp
    \f]

    where:
    \vartable
        A   	| disk area
        U_inf 	| upstream velocity
        Ct	    | thrust coefficient
        Cp	    | power coefficient
        rho     | density
    \endvartable

    \heading Source usage

    Example usage:
    \verbatim
    actuationDiskRingsV11_calib_SourceCoeffs
    {
        fieldNames  (U);
        diskDir     (1 0 0);        // orientation of the flow
        Ct          $Ct;            //Ct
        Cp          $Cp;            //Cp
        Uinf        $Uinf;          //Inlet Velocity
        yaw         $yaw;           //yaw angle (positive is against clock wise)
        pitch       $pitch;         //Not used in this AD
        omega       $omega;         // rad/s
        cellSize    $cellSize;      //fixed min cell size in disc
        diskArea    $Area;          //disk area
        diskPoint   ($diskPoint_x $diskPoint_y $diskPoint_z); //disk center point
        rootFactor  1;              //1-On, 0-Off
        rootDistance    0.1;        //root as a fraction of Diameter [D]
        tipFactor   1;              //1-On, 0-Off
        nodesCellsRatio    2;       //Ratio Total amount of nodes / CellsInAD.
        rThicknessCellsizeRatio  0.5;//Ratio ring Thickness / cellSize
        density     1.225;          //kg/m3. 1.225 means air
        gradInterpolation  1;       //1-On, 0-Off
        #include        "./system/listfvOptions_actuationDiskV11_Source.txt"
    }
    \endverbatim


SourceFiles
    actuationDiskRingsV11_calib_Source.C
    actuationDiskRingsV11_calib_SourceTemplates.C

\*---------------------------------------------------------------------------*/

#include "actuationDiskRingsV11_calib_Source.H"
#include "volFields.H"
#include <math.h>
#include "fvc.H"
#include "fvCFD.H"
#include <map>
// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

/*
    Add resistance to the Ueqn

    Args:
        - Usource
        - cells
        - Vcells
        - rho
        - U
        - force


*/
template <class RhoFieldType>
void Foam::fv::actuationDiskRingsV11_calib_Source::addactuationDiskRingsV11_calib_AxialInertialResistance(
    vectorField &Usource,
    const labelList &cells,
    const scalarField &Vcells,
    const RhoFieldType &rho,
    const vectorField &U,
    vectorField &force) const
{

    //------------------//Information of the actuator line---------------------------

    // Info<< "Actuator disk N: " << this->name() << endl;
    // int tamano = CtList_.size();

    //------------------List definitions------------------------------------------

    // Scalar acumulation
    scalar Pcells = 0.0;      // Power in each line cell
    scalar Tcells = 0.0;      // Thrust in each line cell
    scalar Ftitacells = 0.0;  // Tangential force in each line cell
    scalar TorqueSects = 0.0; // Torque in each seaction
    // scalar weightCells[Vcells.size()]; //wight of All cell volumes for the AD
    // scalar weightCells_U[Vcells.size()]; //wight of All cell volumes for the AD for local velocity
    // scalar weightCellsAD[Vcells.size()]; //wight of All cell volumes (out and inside the sphere)
    // const Field<scalar> zoneCellVolumes(mesh().cellVolumes(), cells);	//line cell volumes

    std::map<int, float> weightCells; // we replace for maps so that we don't have a giant vector
    std::map<int, float> weightCells_U;
    std::map<int, float> weightCellsAD;

    // Vector acumulation
    vector U_dCenterCells_orient = vector(0, 0, 0);                   // Ud for the center cells of the sphere for orientation calculation
    vector U_dCenterCells = vector(0, 0, 0);                          // Ud for the center cells
    vector U_dCells = vector(0, 0, 0);                                // Ud in each cell
    vector U_infCells = vector(0, 0, 0);                              // Uinf in each cell
    const Field<vector> zoneCellCentres(mesh().cellCentres(), cells); // line cell centers

    // dinamic list (with append) of the AL blade points //Information
    DynamicList<scalar> posList;  // position
    DynamicList<scalar> UnList;   // local normal velocity
    DynamicList<scalar> UtList;   // local tangential velocity
    DynamicList<scalar> UrelList; // local relative velocity
    DynamicList<scalar> FnList;   // normal force
    DynamicList<scalar> FtList;   // tangential force

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
            if (mag(mesh().cellCentres()[cells[c]] - diskPoint_) < (centerRation_ * maxR))
            {
                Vcenter_orient += Vcells[cells[c]];
            }
        }
        reduce(Vcenter_orient, sumOp<scalar>());

        // Ud vector for the center cells of the sphere
        forAll(cells, c)
        {
            if (mag(mesh().cellCentres()[cells[c]] - diskPoint_) < (centerRation_ * maxR))
            {
                U_dCenterCells_orient += U[cells[c]] * (Vcells[cells[c]] / Vcenter_orient);
            }
        }
        reduce(U_dCenterCells_orient, sumOp<vector>());

        // no orientation z component
        U_dCenterCells_orient[2] = 0;
        // Info << "self orientation with no Z component "<< U_dCenterCells_orient/mag(U_dCenterCells_orient) << endl;

        // for self orientation depending on the velocity on the disc
        uniDiskDir = U_dCenterCells_orient / mag(U_dCenterCells_orient);
        ////Info<< "uniDiskDir " << uniDiskDir << endl;

        // angle between the self orientation vector and the inlet vector
        vector Vinlet = diskDir_ / mag(diskDir_);
        ////Info << "Vinlet " << Vinlet << endl;
        float cosAlpha = (Vinlet[0] * uniDiskDir[0] + Vinlet[1] * uniDiskDir[1]) / (mag(Vinlet) * mag(uniDiskDir));
        float Alpha = acos(cosAlpha) * 180.0 / M_PI;
        // Info<< "angle between the AD vector and the inlet vector "<< Alpha << endl;
    }
    else
    {
        // Info<<"Yaw specfication activated, uniDiskDir=diskYawed"<<endl;
        ////Info<< "yaw angle: " << yaw_ << endl;
        yawRad = yaw_ * 2 * M_PI / 360;
        // rotate the orginal diskDir with the yaw angle
        vector diskYawed = vector(diskDir_[0] * cos(yawRad) - diskDir_[1] * sin(yawRad), diskDir_[0] * sin(yawRad) + diskDir_[1] * cos(yawRad), diskDir_[2]);
        // Info << "new diskDir_ rotated with the yawed angle "<< diskYawed << endl;
        uniDiskDir = diskYawed;
    }

    //------------------Uref fixed yawed------------------------
    // calculate the component of Uinf for the yawed (or not) AD
    scalar cosUinfAngle = (diskDir_[0] * uniDiskDir[0] + diskDir_[1] * uniDiskDir[1]) / (mag(diskDir_) * mag(uniDiskDir));
    scalar UrefYaw = Uref_ * cosUinfAngle;

    // Info << "Uinf fixed (yawed): "<<UrefYaw<<endl;

    scalar upRho = 1;

    // Info<< "omega(rad/s) fixed: " << omega_ << endl;
    scalar omega = omega_;

    scalar tsr = maxR * omega / UrefYaw;
    // Info<< "TSR: " << tsr << endl;

    scalar Ct = Ct_;

    scalar Cp = Cp_;

    //-----Thrust fixed------------------------------------------

    // calculate the uniform thrust force without density
    float T = 0.50 * upRho * diskArea_ * pow(UrefYaw, 2) * Ct;
    // Info<< "Thrust fixed (with density): " << T *density_ << endl;

    // calculate the uniform thrust forces distribution on the disc [N/m2]
    float fn = T / diskArea_;
    // Info<< "Uniform thrust forces distribution on the disc [N/m2] (with density): " << fn *density_ << endl;

    //-----Torque fixed------------------------------------------
    // Power
    float P = 0.5 * upRho * diskArea_ * pow(UrefYaw, 3) * Cp;
    // Info << "Power fixed [MW] (with density) " << P*density_*0.001 << endl;

    // Torque
    float Torque = 0;
    if (omega > 0)
    {
        Torque = P / omega;
    }

    // Info << "Torque fixed (with density) " << Torque*density_ << endl;

    // calculate the uniform tangential forces distribution on the disc [N/m]
    float ft = (2 * Torque) / (pow(maxR, 2)); // It is multiplied by 2*pi*radius. We will divide it afterwards
    // Info<< "Uniform tangential forces distribution on the disc [N/m] (with density): " << ft *density_ << endl;

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
    // scalar E_U=0.5*minSep;

    DynamicList<label> cellsDisc;

    // Pout<<"total cells in the procesor: "<<cells.size()<<endl;
    forAll(cells, c)
    {
        // distance from AD plane
        scalar d = mag(uniDiskDir[0] * (mesh().cellCentres()[cells[c]][0] - diskPoint_[0]) + uniDiskDir[1] * (mesh().cellCentres()[cells[c]][1] - diskPoint_[1]) + uniDiskDir[2] * (mesh().cellCentres()[cells[c]][2] - diskPoint_[2])) / sqrt(pow(uniDiskDir[0], 2) + pow(uniDiskDir[1], 2) + pow(uniDiskDir[2], 2));

        // distance of the cell center from sphere
        scalar dSphere = mag(mesh().cellCentres()[cells[c]] - diskPoint_);

        if ((dSphere <= maxR * 1.15) and (d <= 3 * E))
        {
            cellsDisc.append(cells[c]);
        }
    } // Por ahora mantenemos el filtro igual a gonza porque está bueno que haya una buena
    // cantidad de celdas a un lado y al otro del AD para suavizar la gaussiana.

    // Pout<<"Disc cells in the procesor: "<<cellsDisc.size()<<endl;

    /*
    forAll(cellsDisc,c)
    {
        Pout<<"mesh().cellCentres()[cellsList[c]] "<<mesh().cellCentres()[cellsDisc[c]]<<endl;
    }
    */

    // Estas lineas tiran el listado de celdas que hay. Lo dejamos así porque seguramente
    //  va a venir muy bien para probar el código.

    //------------------blade section centers------------------------
    // get the time in this step
    scalar t = mesh().time().value();
    // Info << "current time "<< t << endl;

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

    scalar root = maxR * rootDistance_; // distance in meters
    Info << "root: " << root << endl;

    //----first loop over the points to calculate the scale factor due to new force distribution,
    //----only if root or tip factor are on

    // for normal force distribution
    scalar sumF_n_Bi = 0;
    scalar sumF_n_Bixfactor = 0;

    // for tangential force distribution
    scalar sumTorque_Bi = 0;
    scalar sumTorque_Bixfactor = 0;

    //- Velocity field pointer.
    const volVectorField &U_ = mesh().lookupObject<volVectorField>("U");
    volTensorField gradU = fvc::grad(U_);

    /*
    rootFactor_ y tipFactor_ are read in the .C file during initialization
    */
    if ((rootFactor_ == 0) and (tipFactor_ == 0))
    {
        sumF_n_Bi = 1;
        sumF_n_Bixfactor = 1;

        sumTorque_Bi = 1;
        sumTorque_Bixfactor = 1;
    }
    else
    {
        Info << "Starting first loop for calculating scale factor" << endl;
        // for each ring
        for (int ring = 0; ring <= (numberRings_ - 1); ring = ring + 1)
        // Not passing through the last ring to avoid errors with the center node
        {
            tita_r = ringTitaList_[ring];
            rMed_r = ringrMedList_[ring];

            // for each node
            for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
            {
                tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

                // position of the node
                scalar x_node = -1 * rMed_r * sin(tita_n_Rad) * sin(yawRad);
                scalar y_node = rMed_r * sin(tita_n_Rad) * cos(yawRad);
                scalar z_node = rMed_r * cos(tita_n_Rad);

                // move to turbine position
                x_node = x_node + diskPoint_[0];
                y_node = y_node + diskPoint_[1];
                z_node = z_node + diskPoint_[2];

                // blade vector
                vector bladeDir = vector(x_node - diskPoint_[0], y_node - diskPoint_[1], z_node - diskPoint_[2]);
                vector bladeUniDir = bladeDir / mag(bladeDir);

                // calculate the tangential vector
                F_tita_dir = vector(uniDiskDir[1] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[1], -1 * (uniDiskDir[0] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[0]), uniDiskDir[0] * bladeUniDir[1] - uniDiskDir[1] * bladeUniDir[0]); 
                F_tita_dir = F_tita_dir / mag(F_tita_dir);

                // calculate the tensor transformation of coordinates
                vector_n = -1 * uniDiskDir;
                vector_t = F_tita_dir;
                vector_r = bladeUniDir;
                tensor Transform(vector_n[0], vector_t[0], vector_r[0],
                                 vector_n[1], vector_t[1], vector_r[1],
                                 vector_n[2], vector_t[2], vector_r[2]); 
                vector Bi = vector(x_node, y_node, z_node);
                radius = mag(diskPoint_ - Bi);

                // change of coordinate system
                Bi_ntr = inv(Transform) & Bi;

                // calculate velocity in node from velocity in cell
                vector U_dPointCells = vector(1000, 1000, 1000);
                if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
                {
                    U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                    if (gradInterpolation_ == 1) // ?? condition
                    {
                        vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
                        vector dU = dx & gradU[nodeCellID_[total_nodes_counter]]; 
                        U_dPointCells += dU;
                    }
                }
                reduce(U_dPointCells, minOp<vector>()); 

                if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
                {
                    U_dPointCells = vector(10, 0, 0);
                    Info << "OpenFOAM cell Not found" << endl;
                    Info << "ring: " << ring << endl;
                    Info << "node: " << total_nodes_counter << endl;
                    Info << "radius: " << radius << endl;
                }

                total_nodes_counter += 1;

                // change of coordinate system
                U_dPointCells_ntr = inv(Transform) & U_dPointCells; 
                U_n = -1 * U_dPointCells_ntr[0];
                U_t = -1 * U_dPointCells_ntr[1];

                // phi angle, always positive the bottom part
                // PHI ANGLE CALCULATION -------------------------------------------------------------------------------
                if (omega > 0)
                {
                    phi = Foam::atan(U_n / (radius * omega - U_t));
                }
                if (((U_t - radius * omega) <= 0) and (U_n >= 0))
                {
                    ////Info << "case 1) (U_t - radius*omega)(-) and U_z(+) (most common case)"  << endl;
                }
                if (((U_t - radius * omega) <= 0) and (U_n < 0))
                {
                    ////Info << "case 3) (U_t - radius*omega)(-) and U_z(-)"  << endl;
                }
                if (((U_t - radius * omega) > 0) and (U_n >= 0))
                {
                    ////Info << "case 2) (U_t - radius*omega)(+) and U_z(+)"  << endl;
                    phi = phi + M_PI;
                }
                if (((U_t - radius * omega) > 0) and (U_n < 0))
                {
                    ////Info << "case 4) (U_t - radius*omega)(+) and U_z(-)"  << endl;
                    phi = phi - M_PI;
                }
                // PHI ANGLE CALCULATION END -------------------------------------------------------------------------------

                // Info << "radius: "<<radius<<endl;
                // TIP AND ROOT CORRECTION FACTORS ---------------------------------------------------------------
                scalar fcorr = 0.0;

                // TIP CORRECTION FACTORS ---------------------------------------------------------------
                scalar tipfactor = 1;
                scalar tipfactor_f = (nblades_ / 2) * (maxR - radius) / (radius * sin(phi));

                if (tipFactor_ == 1) // tip factor proposed by Shen 2005
                {
                    scalar c1 = 0.125;
                    scalar c2 = 27; // 27
                    scalar c3 = 0.1;
                    scalar g = 1;
                    // scalar tipfactor_f = (nblades_/2)*(maxR - radius)/(radius*sin(phi)); Movemos antes del if
                    // Info << "tipfactor_f " << tipfactor_f << endl;
                    g = exp(-c1 * (nblades_ * tsr - c2)) + c3;
                    if (tipfactor_f > 0)
                    {
                        if (((exp(-g * tipfactor_f)) > -1) and ((exp(-g * tipfactor_f)) < 1))
                        {
                            tipfactor = (2 / (M_PI)) * acos(exp(-g * tipfactor_f));
                        }
                    }
                }
                else if (tipFactor_ == 2) // tip factor proposed by Prandtl 
                {
                    Info << "calculating tipFactor" << endl;
                    scalar tipfactor_arg = exp(-(nblades_ * (1 - radius)) / (2 * radius * sin(phi)));
                    if ((tipfactor_arg > -1) and (tipfactor_arg < 1))
                    {
                        tipfactor = (2 / (M_PI)) * acos(tipfactor_arg);
                    }
                }

                // ROOT CORRECTION FACTORS ---------------------------------------------------------------
                scalar tipfactor = 1;
                scalar rootfactor = 1;
                if (rootFactor_ == 1) // root factor proposed by Glauert
                {
                    if (radius <= root)
                    {
                        rootfactor = 0;
                    }

                    if ((radius > root) and (tipfactor_f > 0) and (radius < maxR / 2.0))
                    {
                        scalar rootfactor_f = (nblades_ / 2) * (radius - 0.1 * maxR) / (radius * sin(phi));
                        scalar g = 1;

                        if (((exp(-g * rootfactor_f)) > -1) and ((exp(-g * rootfactor_f)) < 1))
                        {
                            rootfactor = (2 / (M_PI)) * acos(exp(-g * rootfactor_f));
                        }
                    }
                }
                else if (rootFactor_ == 2) // root factor proposed in Sorensen 2020 
                {
                    Info << "calculating rootFactor" << endl;
                    if (radius <= root)
                    {
                        rootfactor = 0;
                    }
                    if ((radius > root) and (radius < maxR / 2.0))
                    {
                        scalar constantA = 2.335;
                        scalar constantB = 4;
                        rootfactor = 1 - exp(-constantA * pow(radius / (root / maxR), constantB));
                    }
                }

                // Info << "rootfactor: "<<rootfactor<<endl;
                fcorr = tipfactor * rootfactor;
                // END OF TIP AND ROOT FACTORS CALCULATION -------------------------------------------

                //----FORCES CALCULATIONS-------------------
                // for normal force distribution
                sumF_n_Bi += fn * ringAreaList_[ring];
                sumF_n_Bixfactor += fcorr * fn * ringAreaList_[ring];

                // for tangential force distribution
                sumTorque_Bi += (ft / (2 * M_PI * radius)) * ringAreaList_[ring];
                sumTorque_Bixfactor += fcorr * (ft / (2 * M_PI * radius)) * ringAreaList_[ring];
            }
        }

        Info << "Calculating scale factor for central node" << endl;
        // Info << "ring: "<< numberRings_ << endl;
        // Info << "node: "<< nodesNumber_ << endl;
        // Calculation for de centrl node
        sumF_n_Bi += fn * ringAreaList_[numberRings_];
        // We only add normal force, not tangential
        // The fcorr is zero, that's why we don't add to sumF_n_Bixfactor
        // if there is root factor, it is zero. The tip factor is 1 in this point.

    } // close if root and tip factors
    // calculate the scale factors
    scalar scale_factor_F_n_Bi = sumF_n_Bi / sumF_n_Bixfactor;
    scalar scale_factor_Torque_Bi = sumTorque_Bi / sumTorque_Bixfactor;
    // Info<< "scale factor for the Fn distribution: "<< scale_factor_F_n_Bi <<endl;
    // Info<< "scale factor for the torque distribution: "<< scale_factor_Torque_Bi <<endl;

    //----TERMINADO EL PRIMER LOOP PARA EL CÁLCULO DEL SCALE FACTOR ARRANCAMOS EL SEGUNDO
    // LOOP PARA CALCULAR LAS FUERZAS.

    Info << "Starting loop for force calculation" << endl;
    Info << " " << endl;
    total_nodes_counter = 0;

    // for each ring

    //--- LOOP OVER RINGS FOR FORCE CALCULATION AND DISTTIBUTION -----------------------------------------------------------------
    // for (int ring = 0; ring <= (numberRings_ - 1); ring = ring + 1)
    // Not passing through the last ring to avoid errors with the center node
    for (int ring = 0; ring <= (numberRings_); ring = ring + 1)
    // Pass through last ring (center node)
    // FALTA COMENTAR Y DESPS ELIMINAR LA PARTE EN LA QUE HACE SOLO EL NODO DEL CENTRO
    {
        tita_r = ringTitaList_[ring];
        rMed_r = ringrMedList_[ring];
        // Info << "ring: "<< ring << endl;
        // Info << "nodes in ring: "<< ringNodesList_[ring] << endl;
        // Info << "rMed_r: "<< rMed_r << endl;

        // for each node
        // -- LOOP OVER NODES IN RING FOR FORCE CALCULATION AND DISTRBUTION -------------------------------------------------
        for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
        {
            tita_n_Rad = 2 * M_PI * (tita_r * (nodeIterator - 1)) / 360;

            scalar x_node = 0;
            scalar y_node = 0;
            scalar z_node = 0;

            // position of the node considering disk center = (0,0,0)
            if (ring == numberRings_)
            {
                scalar x_node = 0;
                scalar y_node = 0;
                scalar z_node = 0;
            }
            else
            {
                scalar x_node = -1 * rMed_r * sin(tita_n_Rad) * sin(yawRad);
                scalar y_node = rMed_r * sin(tita_n_Rad) * cos(yawRad);
                scalar z_node = rMed_r * cos(tita_n_Rad);
            }

            // move to turbine position
            x_node = x_node + diskPoint_[0];
            y_node = y_node + diskPoint_[1];
            z_node = z_node + diskPoint_[2];

            // blade vector
            vector bladeUniDir = vector(0, 0, 1); // we force this vector for the center node
            if (ring == numberRings_)
            {
                vector bladeUniDir = vector(0, 0, 1); // we force this vector for the center node
            }
            else
            {
                vector bladeDir = vector(x_node - diskPoint_[0], y_node - diskPoint_[1], z_node - diskPoint_[2]);
                vector bladeUniDir = bladeDir / mag(bladeDir);
            }
            ////Info << "blade uni vector : "<< bladeUniDir<< endl;

            // calculate the tangential vector
            F_tita_dir = vector(uniDiskDir[1] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[1],
                                -1 * (uniDiskDir[0] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[0]),
                                uniDiskDir[0] * bladeUniDir[1] - uniDiskDir[1] * bladeUniDir[0]);

            F_tita_dir = F_tita_dir / mag(F_tita_dir);
            ////Info << "F_tita_dir" << F_tita_dir << endl;

            // calculate the tensor transformation of coordinates
            vector_n = -1 * uniDiskDir;
            vector_t = F_tita_dir;
            vector_r = bladeUniDir;

            ////Info << "vector_n " << vector_n << endl;
            ////Info << "vector_t " << vector_t << endl;
            ////Info << "vector_r " << vector_r << endl;

            tensor Transform(vector_n[0], vector_t[0], vector_r[0],
                             vector_n[1], vector_t[1], vector_r[1],
                             vector_n[2], vector_t[2], vector_r[2]);

            vector Bi = vector(x_node, y_node, z_node);

            radius = mag(diskPoint_ - Bi);

            // save for the output
            posList.append(radius);

            // change of coordinate system
            Bi_ntr = inv(Transform) & Bi;

            vector U_dPointCells = vector(1000, 1000, 1000);
            if (nodeCellID_[total_nodes_counter] != -1) // if the closer cell is in this procesor
            {
                U_dPointCells = U[nodeCellID_[total_nodes_counter]];
                if (gradInterpolation_ == 1)
                {
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
            // Info << "node: "<< total_nodes_counter << endl;

            // change of coordinate system
            U_dPointCells_ntr = inv(Transform) & U_dPointCells;

            // velocities in the profile coordinates
            U_n = -1 * U_dPointCells_ntr[0];
            U_t = -1 * U_dPointCells_ntr[1];

            // phi angle, always positive the bottom part
            if (omega > 0)
            {
                phi = Foam::atan(U_n / (radius * omega - U_t));
            }

            if (((U_t - radius * omega) <= 0) and (U_n >= 0))
            {
                ////Info << "case 1) (U_t - radius*omega)(-) and U_z(+) (most common case)"  << endl;
            }

            if (((U_t - radius * omega) <= 0) and (U_n < 0))
            {
                ////Info << "case 3) (U_t - radius*omega)(-) and U_z(-)"  << endl;
            }

            if (((U_t - radius * omega) > 0) and (U_n >= 0))
            {
                ////Info << "case 2) (U_t - radius*omega)(+) and U_z(+)"  << endl;
                phi = phi + M_PI;
            }

            if (((U_t - radius * omega) > 0) and (U_n < 0))
            {
                ////Info << "case 4) (U_t - radius*omega)(+) and U_z(-)"  << endl;
                phi = phi - M_PI;
            }

            // volume of the point cells weighted, for forces
            V_point_F = 0;

            // loop over all the cells to weight for FORCES
            forAll(cellsDisc, c)
            {
                // change of coordinate system
                Pi_ntr = inv(Transform) & mesh().cellCentres()[cellsDisc[c]];

                // calculate the distances in blade coordinate system
                dn = mag(Pi_ntr[0] - Bi_ntr[0]);
                dt = mag(Pi_ntr[1] - Bi_ntr[1]);
                dr = mag(Pi_ntr[2] - Bi_ntr[2]);

                //---for forces
                // calculate the weight dependig on the distances from the sectional point
                // very similar to convolution distribution form Mikkelsen PhD thesis in 2003
                // difference that here it decomposes in the directions.That part inside the exponential isn't exactly the same as in Mikkelsen proposal
                if (forceDistributionMethod_==0) { // no force distribution. Forces are applied directly in cell that contains the node
                    // FINISH THIS PART
                    // !!!!!!
                    // !!!!!!
                    // !!!!!!
                }
                else if (forceDistributionMethod_==1) { // force distribution previously implented in this code
                    weightCells[cellsDisc[c]] = (1 / (En * Et * Er * pow(sqrt(M_PI), 3))) *
                                            exp(-1 * (pow(dn / En, 2) + pow(dt / Et, 2) + pow(dr / Er, 2)));
                }
                else if (forceDistributionMethod_==2) { // force distribution as proposed in Mikkelsen 2003
                    weightCells[cellsDisc[c]] = (1 / (pow(E, 3) * pow(sqrt(M_PI), 3))) *
                                            exp(-1 * ( (sqrt(pow(dn,2) + pow(dr,2) + pow(dt,2))) / (E) ) ));
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

            // volume weighted
            reduce(V_point_F, sumOp<scalar>());

            // Info << "radius: "<<radius<<endl;
            // TIP AND ROOT CORRECTION FACTORS ---------------------------------------------------------------
            scalar fcorr = 0.0;

            // TIP CORRECTION FACTORS ---------------------------------------------------------------
            scalar tipfactor = 1;
            scalar tipfactor_f = (nblades_ / 2) * (maxR - radius) / (radius * sin(phi));

            if (tipFactor_ == 1) // tip factor proposed by Shen 2005
            {
                scalar c1 = 0.125;
                scalar c2 = 27; // 27
                scalar c3 = 0.1;
                scalar g = 1;
                // scalar tipfactor_f = (nblades_/2)*(maxR - radius)/(radius*sin(phi)); Movemos antes del if
                // Info << "tipfactor_f " << tipfactor_f << endl;
                g = exp(-c1 * (nblades_ * tsr - c2)) + c3;
                if (tipfactor_f > 0)
                {
                    if (((exp(-g * tipfactor_f)) > -1) and ((exp(-g * tipfactor_f)) < 1))
                    {
                        tipfactor = (2 / (M_PI)) * acos(exp(-g * tipfactor_f));
                    }
                }
            }
            else if (tipFactor_ == 2) // tip factor proposed by Prandtl 
            {
                Info << "calculating tipFactor" << endl;
                scalar tipfactor_arg = exp(-(nblades_ * (1 - radius)) / (2 * radius * sin(phi)));
                if ((tipfactor_arg > -1) and (tipfactor_arg < 1))
                {
                    tipfactor = (2 / (M_PI)) * acos(tipfactor_arg);
                }
            }

            // ROOT CORRECTION FACTORS ---------------------------------------------------------------
            scalar tipfactor = 1;
            scalar rootfactor = 1;
            if (rootFactor_ == 1) // root factor proposed by Glauert
            {
                if (radius <= root)
                {
                    rootfactor = 0;
                }

                if ((radius > root) and (tipfactor_f > 0) and (radius < maxR / 2.0))
                {
                    scalar rootfactor_f = (nblades_ / 2) * (radius - 0.1 * maxR) / (radius * sin(phi));
                    scalar g = 1;

                    if (((exp(-g * rootfactor_f)) > -1) and ((exp(-g * rootfactor_f)) < 1))
                    {
                        rootfactor = (2 / (M_PI)) * acos(exp(-g * rootfactor_f));
                    }
                }
            }
            else if (rootFactor_ == 2) // root factor proposed in Sorensen 2020 
            {
                Info << "calculating rootFactor" << endl;
                if (radius <= root)
                {
                    rootfactor = 0;
                }
                if ((radius > root) and (radius < maxR / 2.0))
                {
                    scalar constantA = 2.335;
                    scalar constantB = 4;
                    rootfactor = 1 - exp(-constantA * pow(radius / (root / maxR), constantB));
                }
            }

            // Info << "rootfactor: "<<rootfactor<<endl;
            fcorr = tipfactor * rootfactor;
            // END OF TIP AND ROOT FACTORS CALCULATION -------------------------------------------

            // Tangential and axial forces
            // multiply the force by the corresponding area
            F_n_Bi = fcorr * scale_factor_F_n_Bi * fn * ringAreaList_[ring];
            // Info<< "normal force [N] in the point section " << F_n_Bi << endl;

            F_tita_Bi = fcorr * scale_factor_Torque_Bi * (ft / (2 * M_PI * radius)) * ringAreaList_[ring];
            // Info<< "tangential force [N] in the point section " << F_tita_Bi << endl;

            // save for the output
            FnList.append(F_n_Bi * density_);
            FtList.append(F_tita_Bi * density_);
            UnList.append(U_n);
            UtList.append(U_t);

            // global from all the blades
            TorqueSects += mag(diskPoint_ - Bi) * F_tita_Bi * density_;
            Pcells += mag(diskPoint_ - Bi) * F_tita_Bi * density_ * omega;

            // loop over all the cells to apply source
            forAll(cellsDisc, c)
            {
                // apply the thrust force in the cell, volume weighed
                Usource[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_n_Bi) * uniDiskDir;

                // apply the tangential force in the cell, volume weighed
                Usource[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_tita_Bi) * F_tita_dir;

                // save the thrust force in the force field for paraview
                force[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_n_Bi) * uniDiskDir * -1 * density_ / Vcells[cellsDisc[c]];
                // force[cellsDisc[c]]+= vector(0.001,0,0);

                // save the tangential force in the force field for praview
                force[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_tita_Bi) * F_tita_dir * -1 * density_ / Vcells[cellsDisc[c]];

                // save the results of each cell, volume weighed
                Tcells += F_n_Bi * density_ * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);

                U_infCells += uniDiskDir * UrefYaw * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);

                Ftitacells += F_tita_Bi * density_ * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);

            } // close loop cells in sectional point

        } // close loop nodes in this rings

    } // close loop rings

    Info << "Force calculation for the center node" << endl;
    Info << " " << endl;

    // Calculate for center node
    // Info << "ring: "<< numberRings_ << endl;
    // Info << "node: "<< nodesNumber_ << endl;

    // move to turbine position
    scalar x_node = diskPoint_[0];
    scalar y_node = diskPoint_[1];
    scalar z_node = diskPoint_[2];

    // blade vector
    vector bladeUniDir = vector(0, 0, 1); // we force this vector for the center node
    ////Info << "blade uni vector : "<< bladeUniDir<< endl;

    // calculate the tangential vector
    F_tita_dir = vector(uniDiskDir[1] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[1],
                        -1 * (uniDiskDir[0] * bladeUniDir[2] - uniDiskDir[2] * bladeUniDir[0]),
                        uniDiskDir[0] * bladeUniDir[1] - uniDiskDir[1] * bladeUniDir[0]);

    F_tita_dir = F_tita_dir / mag(F_tita_dir);
    ////Info << "F_tita_dir" << F_tita_dir << endl;

    // calculate the tensor transformation of coordinates
    vector_n = -1 * uniDiskDir;
    vector_t = F_tita_dir;
    vector_r = bladeUniDir;

    ////Info << "vector_n " << vector_n << endl;
    ////Info << "vector_t " << vector_t << endl;
    ////Info << "vector_r " << vector_r << endl;

    tensor Transform(vector_n[0], vector_t[0], vector_r[0],
                     vector_n[1], vector_t[1], vector_r[1],
                     vector_n[2], vector_t[2], vector_r[2]);

    vector Bi = vector(x_node, y_node, z_node);

    radius = mag(diskPoint_ - Bi);

    // save for the output
    posList.append(radius);

    // change of coordinate system
    Bi_ntr = inv(Transform) & Bi;

    vector U_dPointCells = vector(1000, 1000, 1000);
    if (nodeCellID_[nodesNumber_ - 1] != -1) // if the closer cell is in this procesor
    {
        U_dPointCells = U[nodeCellID_[nodesNumber_ - 1]];
        if (gradInterpolation_ == 1)
        {
            vector dx = Bi - mesh().cellCentres()[nodeCellID_[nodesNumber_ - 1]];
            vector dU = dx & gradU[nodeCellID_[nodesNumber_ - 1]];
            U_dPointCells += dU;
        }
    }
    reduce(U_dPointCells, minOp<vector>()); // take only normal values of U

    if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
    {
        U_dPointCells = vector(10, 0, 0);
        Info << "OpenFOAM cell Not found" << endl;
        Info << "ring: " << numberRings_ << endl;
        Info << "node: " << nodesNumber_ << endl;
        Info << "radius: " << radius << endl;
    }

    // change of coordinate system
    U_dPointCells_ntr = inv(Transform) & U_dPointCells;

    // velocities in the profile coordinates
    U_n = -1 * U_dPointCells_ntr[0];
    U_t = -1 * U_dPointCells_ntr[1];

    // volume of the point cells weighted, for forces
    V_point_F = 0;

    // loop over all the cells to weight for FORCES
    forAll(cellsDisc, c)
    {
        // change of coordinate system
        Pi_ntr = inv(Transform) & mesh().cellCentres()[cellsDisc[c]];

        // calculate the distances in blade coordinate system
        dn = mag(Pi_ntr[0] - Bi_ntr[0]);
        dt = mag(Pi_ntr[1] - Bi_ntr[1]);
        dr = mag(Pi_ntr[2] - Bi_ntr[2]);

        //---for forces
        // calculate the wight dependig on the distances from the sectional point
        weightCells[cellsDisc[c]] = (1 / (En * Et * Er * pow(sqrt(M_PI), 3))) *
                                    exp(-1 * (pow(dn / En, 2) + pow(dt / Et, 2) + pow(dr / Er, 2)));

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

    // We avoid calculating root and tip factor for the center node. Tip is always
    // 1 and root is zero in case we apply it.

    // Info << "radius: "<<radius<<endl;
    // Shen tip correction factor:
    scalar fcorr = 0.0;
    scalar tipfactor = 1;

    // Info << "tipfactor: "<<tipfactor<<endl;
    // Glauert root correction factor:
    scalar rootfactor = 1;
    if (rootFactor_ == 1 || rootFactor_ == 2) // If root factor is on
    {
        rootfactor = 0;
    }

    // Info << "rootfactor: "<<rootfactor<<endl;
    fcorr = tipfactor * rootfactor;
    // Tangential and axial forces

    // multiply the force by the corresponding area
    F_n_Bi = fcorr * scale_factor_F_n_Bi * fn * ringAreaList_[numberRings_];
    // Info<< "normal force [N] in center node: " << F_n_Bi << endl;

    F_tita_Bi = 0; // No tangential force at center node.
    // Info<< "tangential force [N] in the point section " << F_tita_Bi << endl;

    // save for the output
    FnList.append(F_n_Bi * density_);
    FtList.append(F_tita_Bi * density_);
    UnList.append(U_n);
    UtList.append(U_t);

    // global from all the blades
    TorqueSects += mag(diskPoint_ - Bi) * F_tita_Bi * density_;
    Pcells += mag(diskPoint_ - Bi) * F_tita_Bi * density_ * omega;

    // loop over all the cells to apply source
    forAll(cellsDisc, c)
    {
        // apply the thrust force in the cell, volume weighed
        Usource[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_n_Bi) * uniDiskDir;

        // apply the tangential force in the cell, volume weighed
        Usource[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_tita_Bi) * F_tita_dir;

        // save the thrust force in the force field for paraview
        force[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_n_Bi) * uniDiskDir * -1 * density_ / Vcells[cellsDisc[c]];
        // force[cellsDisc[c]]+= vector(0.001,0,0);

        // save the tangential force in the force field for praview
        force[cellsDisc[c]] += (((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F) * F_tita_Bi) * F_tita_dir * -1 * density_ / Vcells[cellsDisc[c]];

        // save the results of each cell, volume weighed
        Tcells += F_n_Bi * density_ * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);

        U_infCells += uniDiskDir * UrefYaw * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);

        Ftitacells += F_tita_Bi * density_ * ((Vcells[cellsDisc[c]] * weightCells[cellsDisc[c]]) / V_point_F);
    }

    //---colecting data from all the procesors
    // reduce(Pcells, sumOp<scalar>());
    ////Info<< "Pcells [W] (with density) from all the sections: " << Pcells/(3*numberSect) << endl;
    ////Info<< "Pcells [W] (with density) from all the sections: " << Pcells << endl;
    reduce(Tcells, sumOp<scalar>());
    // Info<< "Tcells (with density) from all the sections: " << Tcells << endl;
    reduce(Ftitacells, sumOp<scalar>());
    // Info<< "Ftitacells (with density) from all the sections: " << Ftitacells << endl;
    reduce(U_infCells, sumOp<vector>());
    // Info<< "U_infCells from all the sections: " << mag(U_infCells)/(3*numberSect) << endl;
    // reduce(TorqueSects, sumOp<scalar>());
    ////Info<< "TorqueSects (with density) from all the sections: " << TorqueSects << endl;

    //------------------Data from the disc-----------------------------------------------

    //------------------wight of the AD cells depending distance from plane and sphere-----

    // loop over the AD cells for weight calculation
    forAll(cellsDisc, c)
    {
        // distance from AD plane
        scalar d = mag(uniDiskDir[0] * (mesh().cellCentres()[cellsDisc[c]][0] - diskPoint_[0]) + uniDiskDir[1] * (mesh().cellCentres()[cellsDisc[c]][1] - diskPoint_[1]) + uniDiskDir[2] * (mesh().cellCentres()[cellsDisc[c]][2] - diskPoint_[2])) / sqrt(pow(uniDiskDir[0], 2) + pow(uniDiskDir[1], 2) + pow(uniDiskDir[2], 2));

        // evaluate weight gauss bell
        weightCellsAD[cellsDisc[c]] = (1 / (E * sqrt(M_PI))) * exp(-1 * pow((d / E), 2));

        // distance from sphere
        scalar dSphere = mag(mesh().cellCentres()[cellsDisc[c]] - diskPoint_);
        // evaluate weight
        if (dSphere <= maxR)
        {
            weightCellsAD[cellsDisc[c]] = weightCellsAD[cellsDisc[c]] * 1;
        }
        else
        {
            weightCellsAD[cellsDisc[c]] = weightCellsAD[cellsDisc[c]] * 0;
        }
    }

    //---volume of the center cells weighted
    scalar Vcenter = 0.0;

    forAll(cellsDisc, c)
    {

        if (mag(mesh().cellCentres()[cellsDisc[c]] - diskPoint_) < (0.30 * maxR))
        {
            Vcenter += Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]];
        }
    }

    reduce(Vcenter, sumOp<scalar>());

    //---volume of the AD cells weighted
    scalar V_AD = 0.0;

    forAll(cellsDisc, c)
    {
        V_AD += Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]];
    }

    reduce(V_AD, sumOp<scalar>());
    ////Info << "theorical AD volume: " << diskArea_ *minSep*2 << endl;
    ////Info << "total AD volume weighted: " << V_AD << endl;

    //---Ud for the center cells
    forAll(cellsDisc, c)
    {

        if (mag(mesh().cellCentres()[cellsDisc[c]] - diskPoint_) < (0.30 * maxR))
        {
            U_dCenterCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]]) / Vcenter);
        }
    }

    reduce(U_dCenterCells, sumOp<vector>());

    scalar cosU_dCenterCells = (U_dCenterCells[0] * uniDiskDir[0] + U_dCenterCells[1] * uniDiskDir[2]) / (mag(U_dCenterCells) * mag(uniDiskDir));
    scalar U_dCenterCellsYaw = mag(U_dCenterCells) * cosU_dCenterCells;

    // Info<< "U_dCenterCells not yawed: " << mag(U_dCenterCells) << endl;
    // Info<< "U_dCenterCells yawed: " << U_dCenterCellsYaw << endl;

    //---Ud for all the cells
    forAll(cellsDisc, c)
    {
        U_dCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]]) / V_AD);
    }

    reduce(U_dCells, sumOp<vector>());

    scalar cosU_dCells = (U_dCells[0] * uniDiskDir[0] + U_dCells[1] * uniDiskDir[2]) / (mag(U_dCells) * mag(uniDiskDir));
    scalar U_dCellsYaw = mag(U_dCells) * cosU_dCells;

    // Info<< "U_dCells not yawed: " << mag(U_dCells) << endl;
    // Info<< "U_dCells yawed: " << U_dCellsYaw << endl;

    //-----print results------------------------------------------
    //---global turbine outputs
    if (Pstream::myProcNo() == 0) // if Im in the master proccesor
    {
        //"Turbine, time [s], Uinf(fixed) [m/s], Cp(fixed), Ct(fixed), Power(Uinf,Cp) [W], Thrust(Uinf,Ct) [N] , Center Ud [m/s], Average Ud [m/s], Uinf(disc cells) [m/s], Power(disc cells) [W], Thrust(disc cells) [N], Torque(disc cells) [Nm]"
        (*outTurbines) << this->name() << "," << t << "," << Uref_ << "," << Cp_ << "," << Ct_ << "," << P * density_ << "," << T * density_ << "," << mag(U_dCenterCells) << "," << mag(U_dCells) << "," << mag(U_infCells) << "," << Pcells << "," << Tcells << "," << TorqueSects << std::endl;
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
}
// ************************************************************************* //
