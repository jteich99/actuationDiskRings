/*---------------------------------------------------------------------------*\
Author: Dimas Alejandro Barile
CSC - CONICET, Buenos Aires, 2022

-----------------------------------------------------------------------------
OpenFOAM 2.4

Actuation Disk

COMPLETAR

-----------------------------------------------------------------------------

Class
    Foam::fv::actuationDiskRingsV11_Source

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
    actuationDiskRingsV11_SourceCoeffs
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
    actuationDiskRingsV11_Source.C
    actuationDiskRingsV11_SourceTemplates.C

\*---------------------------------------------------------------------------*/

#include "actuationDiskRingsV11_Source.H"
#include "volFields.H"
#include <math.h>
#include "fvc.H"
#include "fvCFD.H"
#include <map>
// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template <class RhoFieldType>
void Foam::fv::actuationDiskRingsV11_Source::addactuationDiskRings_AxialInertialResistance(
    vectorField &Usource,
    const labelList &cells,
    const scalarField &Vcells,
    const RhoFieldType &rho,
    const vectorField &U,
    vectorField &force) const
{

    //------------------//Information of the actuator line---------------------------

    // Info<< "---Actuator disk N: " << this->name() << endl;
    // int tamano = CtList_.size();

    //------------------List definitions------------------------------------------

    // Scalar acumulation
    scalar Pcells = 0.0;      // Power in each line cell
    scalar Tcells = 0.0;      // Thrust in each line cell
    scalar Ftitacells = 0.0;  // Tangential force in each line cell
    scalar TorqueSects = 0.0; // Torque in each seaction
    // scalar weightCells[Vcells.size()]; //wight of All cell volumes for the AL
    // scalar weightCells_U[Vcells.size()]; //wight of All cell volumes for the AL for local velocity
    // scalar weightCellsAD[Vcells.size()]; //wight of All cell volumes (out and inside the sphere)
    // const Field<scalar> zoneCellVolumes(mesh().cellVolumes(), cells);	//line cell volumes

    std::map<int, float> weightCells; // we replace for maps so that we don't have a giant vector
    std::map<int, float> weightCells_U;
    std::map<int, float> weightCellsAD;
    std::map<int, float> weightCellsADCenter;

    // Vector acumulation
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

    //--------------------line orientation------------------------------------------------
    //-----Yaw rotation option
    vector uniDiskDir = vector(0, 0, 0);
    scalar yawRad;

    if (yaw_ == 360)
    {
        // Info<<"Because yaw = 360, no yaw restriction is applied, free orientation aviable"<< endl;

        // volume of the center cells of the sphere

        scalar Vcenter_orient = 0.0;

        forAll(cells, c)
        {

            if (mag(mesh().cellCentres()[cells[c]] - diskPoint_) < (0.30 * maxR))
            {
                Vcenter_orient += Vcells[cells[c]];
            }
        }

        reduce(Vcenter_orient, sumOp<scalar>());

        // Ud vector for the center cells of the sphere
        forAll(cells, c)
        {

            if (mag(mesh().cellCentres()[cells[c]] - diskPoint_) < (0.30 * maxR))
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
        // Info<< "yaw angle: " << yaw_ << endl;
        yawRad = yaw_ * 2 * M_PI / 360;

        // rotate the orginal diskDir with the yaw angle
        vector diskYawed = vector(diskDir_[0] * cos(yawRad) - diskDir_[1] * sin(yawRad), diskDir_[0] * sin(yawRad) + diskDir_[1] * cos(yawRad), diskDir_[2]);
        // Info << "new diskDir_ rotated with the yawed angle "<< diskYawed << endl;
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
    }
    Info << "new cell selection for the AD" << endl;

    //------------------wight of the AD cells depending distance from plane and sphere-----

    scalar Vcenter = 0.0;
    scalar centerRatio = centerRatio_; // in the future it has to be determined by user in fvOptions --> it will unlock to make a study of it's influence in the calculation of Ud
    scalar V_AD = 0.0;

    // loop over the AD cells for weight and volume calculation
    forAll(cellsDisc, c)
    {
        //---weight calculation
        // distance from AD plane
        scalar d = mag(uniDiskDir[0] * (mesh().cellCentres()[cellsDisc[c]][0] - diskPoint_[0]) + uniDiskDir[1] * (mesh().cellCentres()[cellsDisc[c]][1] - diskPoint_[1]) + uniDiskDir[2] * (mesh().cellCentres()[cellsDisc[c]][2] - diskPoint_[2])) / sqrt(pow(uniDiskDir[0], 2) + pow(uniDiskDir[1], 2) + pow(uniDiskDir[2], 2));

        // distance from sphere
        scalar dSphere = mag(mesh().cellCentres()[cellsDisc[c]] - diskPoint_);

        // weight calculation
        if (UdCellsMethod_ == 0)
        {
            // sum directly values in cells in AD
            scalar weightADplane = 1;
            scalar weightSphereCenter = 1;
            scalar weightSphereAD = 1;
        }
        else if (UdCellsMethod_ == 1)
        {
            // weight with Gaussian in distance to AD plane
            scalar weightADplane = (1 / (E * sqrt(M_PI))) * exp(-1 * pow((d / E), 2));
            scalar weightSphereCenter = 1;
            scalar weightSphereAD = 1;
        }
        else if (UdCellsMethod_ == 2)
        {
            // weight with Gaussian in distance to AD plane + distance to center of AD
            scalar weightADplane = (1 / (E * sqrt(M_PI))) * exp(-1 * pow((d / E), 2));
            scalar weightSphereCenter = (1 / ((centerRatio * maxR) * sqrt(M_PI))) * exp(-1 * pow((dSphere / (centerRatio * maxR)), 2));
            scalar weightSphereAD = (1 / (maxR * sqrt(M_PI))) * exp(-1 * pow((dSphere / maxR), 2));
        }

        // evaluate weight
        if (dSphere <= maxR)
        {
            weightCellsAD[cellsDisc[c]] = weightADplane * weightSphereAD;
            if (dSphere <= centerRatio_ * maxR)
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

    ////Info << "theorical AD volume: " << diskArea_ *minSep*2 << endl;
    ////Info << "total AD volume weighted: " << V_AD << endl;

    forAll(cellsDisc, c)
    {
        // dejo el cálculo de U_dCenterCells porque creo que se usa aguas abajo, pero todavia no se bien en qué y no quiero romper nada
        U_dCenterCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsADCenter[cellsDisc[c]]) / Vcenter);

        if (UdCenterToggle_ == 0)
        {
            U_dCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]]) / V_AD);
        }
        else if (UdCenterToggle == 1)
        {
            U_dCells += U[cellsDisc[c]] * ((Vcells[cellsDisc[c]] * weightCellsADCenter[cellsDisc[c]]) / Vcenter);
        }
    }
    reduce(U_dCells, sumOp<vector>());

    // U_dCells += U_dCenterCells;
    reduce(U_dCenterCells, sumOp<vector>());

    // scalar cosU_dCenterCells= (U_dCenterCells[0]*uniDiskDir[0]+U_dCenterCells[1]*uniDiskDir[2])/(mag(U_dCenterCells)*mag(uniDiskDir));
    // scalar U_dCenterCellsYaw=mag(U_dCenterCells)*cosU_dCenterCells;

    Info << "U_dCenterCells not yawed: " << mag(U_dCenterCells) << endl;
    // Info<< "U_dCenterCells yawed: " << U_dCenterCellsYaw << endl;

    // ineficiente que esté acá. lo mando arriba
    //---Ud for all the cells
    // forAll(cellsDisc, c)
    // {
    // 	U_dCells += U[cellsDisc[c]] * ( ( Vcells[cellsDisc[c]] * weightCellsAD[cellsDisc[c]] ) / V_AD );
    // }

    Info << "U_dCells not yawed: " << mag(U_dCells) << endl;

    scalar cosU_dCells = (U_dCells[0] * uniDiskDir[0] + U_dCells[1] * uniDiskDir[1]) / (mag(U_dCells) * mag(uniDiskDir));
    scalar U_dCellsYaw = mag(U_dCells) * cosU_dCells;

    Info << "U_dCells yawed: " << U_dCellsYaw << endl;

    //-----Using U_dCells, search in the UdAvgList_----------------------------

    float difference = GREAT;
    int pos = 0;
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

    int pos1 = pos;
    int pos2 = pos + 1;

    Info << "U_dCells: " << U_dCells << endl;
    // Info<< "postion of the bottom value in List: " << pos1 << endl;

    //-----U_up for the cell in the up stream point
    scalar upRho = 1;

    //----Interpolate omega and pitch-----------------------------------------------

    scalar UrefYaw = ((mag(U_dCells) - UdAvgList_[pos]) * ((UrefList_[pos + 1]) - (UrefList_[pos])) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + (UrefList_[pos]);
    Info << "UrefYaw(m/s) interpolated: " << UrefYaw << endl;

    scalar omega = ((mag(U_dCells) - UdAvgList_[pos]) * (omegaList_[pos + 1] - omegaList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + omegaList_[pos];
    Info << "omega(rad/s) interpolated: " << omega << endl;
    // scalar tsr = maxR*omega/UrefYaw;
    // Info<< "TSR: " << tsr << endl;

    scalar pitch = ((UrefYaw - UrefList_[pos]) * ((pitchList_[pos + 1]) - (pitchList_[pos])) / (UrefList_[pos + 1] - UrefList_[pos])) + (pitchList_[pos]);
    Info << "pitch(deg) interpolated " << pitch << endl;
    pitch = pitch * 2 * M_PI / 360;

    scalar Ct = ((mag(U_dCells) - UdAvgList_[pos]) * (CtList_[pos + 1] - CtList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + CtList_[pos];
    Info << "Ct interpolated " << Ct << endl;

    scalar Cp = ((mag(U_dCells) - UdAvgList_[pos]) * (CpList_[pos + 1] - CpList_[pos]) / (UdAvgList_[pos + 1] - UdAvgList_[pos])) + CpList_[pos];
    Info << "Cp interpolated " << Cp << endl;

    //-----Thrust fixed------------------------------------------
    float density_ = 1.225;

    // calculate the uniform thrust force with out density
    float T = 0.50 * upRho * diskArea_ * pow(UrefYaw, 2) * Ct;
    Info << "Thrust fixed (with density): " << T * density_ << endl;

    // calculate the uniform thrust forces distribution on the disc [N/m2]
    // float fn = T/diskArea_;
    // Info<< "Uniform thrust forces distribution on the disc [N/m2] (with density): " << fn *density_ << endl;

    //-----Torque fixed------------------------------------------
    // Power
    float P = 0.5 * upRho * diskArea_ * pow(UrefYaw, 3) * Cp;
    Info << "Power fixed [MW] (with density) " << P * density_ * 0.001 << endl;

    // Torque
    float Torque = 0;
    if (omega > 0)
    {
        Torque = P / omega;
    }

    // Info << "Torque fixed (with density) " << Torque*density_ << endl;

    // calculate the uniform tangential forces distribution on the disc [N/m]
    // float ft = (2*Torque) / (pow(maxR,2));
    // Info<< "Uniform tangential forces distribution on the disc [N/m] (with density): " << ft *density_ << endl;

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

    // count how many radius sections are
    int nR_ = 1;
    while (rList_[nR_] - rList_[nR_ - 1] > 0)
    {
        nR_ = nR_ + 1;
    }
    // Info<< "number of radius sections in table 2 " << nR_ << endl;

    scalar clockwiseSign = 1;
    // angle of the rotation in this time
    bool clockwise = true;
    if (clockwise)
    {
        clockwiseSign = -1;
    }

    //- Velocity field pointer.
    const volVectorField &U_ = mesh().lookupObject<volVectorField>("U");
    volTensorField gradU = fvc::grad(U_);

    Info << "Starting loop through nodes" << endl;
    Info << " " << endl;
    total_nodes_counter = 0;
    // for each ring

    // loop through rings and nodes for calculating forces and distributing them
    for (int ring = 0; ring <= (numberRings_ - 1); ring = ring + 1)
    // Not passing through the last ring to avoid errors with the center node
    {
        tita_r = ringTitaList_[ring];
        rMed_r = ringrMedList_[ring];
        // Info << "ring: "<< ring << endl;
        // Info << "nodes in ring: "<< ringNodesList_[ring] << endl;
        // Info << "rMed_r: "<< rMed_r << endl;

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

            //----FORCES INTERPOLATION-------------------

            // search for the most proximate radius in table

            /*
            scalar dist=VGREAT;
            rtable=0;
            posr=0;
            while( mag(radius-rList_[posr])<dist )
            {
                dist=mag(radius-rList_[posr]);
                rtable=rList_[posr];
                posr=posr+1;
            }
            */

            posr = posrList_[ring + 1];
            rtable = rNodeList_[ring + 1];

            // Info<< "Closer rtable: " << rtable << endl;

            //----interpolation procces-----------------
            // interpolate U_inf1 (bottom value)
            difference = 1000;
            pos = 0;
            // Info<< "mag(U_dPointCells_ntr) " << mag(U_dPointCells_ntr) << endl;
            // Info<<"Using UrefList_[pos1]" << endl;
            for (int i = 0; i < (Uref2List_.size() - 1); i = i + 1)
            {
                if ((Uref2List_[i] == UrefList_[pos1]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) >= 0))
                {
                    // Info<< "Uref2List_[i]: " << Uref2List_[i] << endl;
                    // Info<< "rList_[i]: " << rList_[i] << endl;
                    // Info<< "UdiList_[i]: " << UdiList_[i] << endl;

                    difference = fabs(mag(U_dPointCells_ntr) - UdiList_[i]);
                    // Info<< "difference " << difference << endl;

                    pos = i; // the position of the lower value
                }

                // just in case the value is outside the table
                if ((Uref2List_[i] == UrefList_[pos1]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) < 0))
                {
                    lastPos = i;
                }
            }

            if (difference == 1000)
            {
                // just in case the value is outside the table
                pos = lastPos;
            }

            U_inf1 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((UinfList_[pos + nR_]) - (UinfList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (UinfList_[pos]);
            fn1 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((fnList_[pos + nR_]) - (fnList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (fnList_[pos]);
            ft1 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((ftList_[pos + nR_]) - (ftList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (ftList_[pos]);

            // Info<< "U_inf1 : " << U_inf1  << endl;
            // Info<< "fn1 : " << fn1  << endl;
            // Info<< "ft1 : " << ft1  << endl;

            // interpolate U_inf2 (top value)
            difference = 1000;
            pos = 0;
            // Info<<"Using UrefList_[pos2]" << UrefList_[pos2]<< endl;
            for (int i = 0; i < (Uref2List_.size() - 1); i = i + 1)
            {
                if ((Uref2List_[i] == UrefList_[pos2]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) >= 0))
                {

                    // Info<< "Uref2List_[i]: " << Uref2List_[i] << endl;
                    // Info<< "UdiList_[i]: " << UdiList_[i] << endl;

                    difference = fabs(mag(U_dPointCells_ntr) - UdiList_[i]);
                    pos = i; // the position of the lower value
                }
                // Info<< "end of the loop" << endl;
                // just in case the value is outside the table
                if ((Uref2List_[i] == UrefList_[pos2]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) < 0))
                {
                    lastPos = i;
                    // Info << "lastPos = i" << endl;
                }
            }

            if (difference == 1000)
            {
                // just in case the value is outside the table
                pos = lastPos;
                // Info<< "the value is outside the table, pos=lastPos" << endl;
            }

            U_inf2 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((UinfList_[pos + nR_]) - (UinfList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (UinfList_[pos]);
            fn2 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((fnList_[pos + nR_]) - (fnList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (fnList_[pos]);
            ft2 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((ftList_[pos + nR_]) - (ftList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (ftList_[pos]);

            // Info<< "U_inf2 : " << U_inf2  << endl;
            // Info<< "fn2 : " << fn2  << endl;
            // Info<< "ft2 : " << ft2  << endl;

            // interpolate U_inf_point
            U_inf_point = ((UrefYaw - UrefList_[pos1]) * ((U_inf2) - (U_inf1)) / (UrefList_[pos2] - UrefList_[pos1])) + (U_inf1);
            fn_point = ((UrefYaw - UrefList_[pos1]) * ((fn2) - (fn1)) / (UrefList_[pos2] - UrefList_[pos1])) + (fn1);
            ft_point = ((UrefYaw - UrefList_[pos1]) * ((ft2) - (ft1)) / (UrefList_[pos2] - UrefList_[pos1])) + (ft1);

            // Info<< "U_inf_point : " << U_inf_point  << endl;
            // Info<< "fn_point [N/m] : " << fn_point  << endl;
            // Info<< "ft_point [N/m]: " << ft_point  << endl;

            //----FORCES CALCULATIONS-------------------
            // Info <<"----FORCES calculations for the disc----" << endl;

            // Tangential and axial forces
            // devide the force depending on the amount of artificial blades (but because the force in table
            // is refered to 1 blade, i need to multiply by 3)
            F_n_Bi = (fn_point * ringThickness_ * 3) / ringNodesList_[ring];
            F_n_Bi = F_n_Bi / density_; // without density to the solver
                                        // Info<< "normal force [N] in the point section (with density included) " << F_n_Bi*density_ << endl;

            F_tita_Bi = (ft_point * ringThickness_ * 3) / ringNodesList_[ring];
            F_tita_Bi = F_tita_Bi / density_; // without density to the solver
                                              // Info<< "tangential force [N] in the point section (with density included) " << F_tita_Bi*density_ << endl;

            // save for the output
            FnList.append(F_n_Bi * density_);
            FtList.append(F_tita_Bi * density_);
            UnList.append(U_n);
            UtList.append(U_t);

            // global from all the blades
            TorqueSects += mag(diskPoint_ - Bi) * F_tita_Bi * density_;
            Pcells += mag(diskPoint_ - Bi) * F_tita_Bi * density_ * omega;

            // Info<< "-----FORCES finished ----- " << endl;

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

        } // close loop in ring

    } // close loop rings

    Info << "Interpolation in the center node" << endl;
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

    //----FORCES INTERPOLATION-------------------

    // search for the most proximate radius in table
    rtable = 0;
    posr = 0;

    //----interpolation procces-----------------
    // interpolate U_inf1 (bottom value)
    difference = 1000;
    pos = 0;
    // Info<< "mag(U_dPointCells_ntr) " << mag(U_dPointCells_ntr) << endl;
    // Info<<"Using UrefList_[pos1]" << endl;
    for (int i = 0; i < (Uref2List_.size() - 1); i = i + 1)
    {
        if ((Uref2List_[i] == UrefList_[pos1]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) >= 0))
        {
            // Info<< "Uref2List_[i]: " << Uref2List_[i] << endl;
            // Info<< "rList_[i]: " << rList_[i] << endl;
            // Info<< "UdiList_[i]: " << UdiList_[i] << endl;

            difference = fabs(mag(U_dPointCells_ntr) - UdiList_[i]);
            // Info<< "difference " << difference << endl;

            pos = i; // the position of the lower value
        }

        // just in case the value is outside the table
        if ((Uref2List_[i] == UrefList_[pos1]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) < 0))
        {
            lastPos = i;
        }
    }

    if (difference == 1000)
    {
        // just in case the value is outside the table
        pos = lastPos;
    }

    U_inf1 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((UinfList_[pos + nR_]) - (UinfList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (UinfList_[pos]);
    fn1 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((fnList_[pos + nR_]) - (fnList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (fnList_[pos]);
    ft1 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((ftList_[pos + nR_]) - (ftList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (ftList_[pos]);

    // Info<< "U_inf1 : " << U_inf1  << endl;
    // Info<< "fn1 : " << fn1  << endl;
    // Info<< "ft1 : " << ft1  << endl;

    // interpolate U_inf2 (top value)
    difference = 1000;
    pos = 0;
    // Info<<"Using UrefList_[pos2]" << UrefList_[pos2]<< endl;
    for (int i = 0; i < (Uref2List_.size() - 1); i = i + 1)
    {
        if ((Uref2List_[i] == UrefList_[pos2]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) >= 0))
        {

            // Info<< "Uref2List_[i]: " << Uref2List_[i] << endl;
            // Info<< "UdiList_[i]: " << UdiList_[i] << endl;

            difference = fabs(mag(U_dPointCells_ntr) - UdiList_[i]);
            pos = i; // the position of the lower value
        }
        // Info<< "end of the loop" << endl;
        // just in case the value is outside the table
        if ((Uref2List_[i] == UrefList_[pos2]) and (rList_[i] == rtable) and (fabs(mag(U_dPointCells_ntr) - UdiList_[i]) < difference) && ((mag(U_dPointCells_ntr) - UdiList_[i]) < 0))
        {
            lastPos = i;
            // Info << "lastPos = i" << endl;
        }
    }

    if (difference == 1000)
    {
        // just in case the value is outside the table
        pos = lastPos;
        // Info<< "the value is outside the table, pos=lastPos" << endl;
    }

    U_inf2 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((UinfList_[pos + nR_]) - (UinfList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (UinfList_[pos]);
    fn2 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((fnList_[pos + nR_]) - (fnList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (fnList_[pos]);
    ft2 = ((mag(U_dPointCells_ntr) - UdiList_[pos]) * ((ftList_[pos + nR_]) - (ftList_[pos])) / (UdiList_[pos + nR_] - UdiList_[pos])) + (ftList_[pos]);

    // Info<< "U_inf2 : " << U_inf2  << endl;
    // Info<< "fn2 : " << fn2  << endl;
    // Info<< "ft2 : " << ft2  << endl;

    // interpolate U_inf_point
    U_inf_point = ((UrefYaw - UrefList_[pos1]) * ((U_inf2) - (U_inf1)) / (UrefList_[pos2] - UrefList_[pos1])) + (U_inf1);
    fn_point = ((UrefYaw - UrefList_[pos1]) * ((fn2) - (fn1)) / (UrefList_[pos2] - UrefList_[pos1])) + (fn1);
    ft_point = ((UrefYaw - UrefList_[pos1]) * ((ft2) - (ft1)) / (UrefList_[pos2] - UrefList_[pos1])) + (ft1);

    // Info<< "U_inf_point : " << U_inf_point  << endl;
    // Info<< "fn_point [N/m] : " << fn_point  << endl;
    // Info<< "ft_point [N/m]: " << ft_point  << endl;

    //----FORCES CALCULATIONS-------------------
    // Info <<"----FORCES calculations for the disc----" << endl;

    // Tangential and axial forces
    // devide the force depending on the amount of artificial blades (but because the force in table
    // is refered to 1 blade, i need to multiply by 3)
    F_n_Bi = (fn_point * rInt_ * 3);
    F_n_Bi = F_n_Bi / density_; // without density to the solver
                                // Info<< "normal force [N] in the point section (with density included) " << F_n_Bi*density_ << endl;

    F_tita_Bi = (ft_point * rInt_ * 3);
    F_tita_Bi = F_tita_Bi / density_; // without density to the solver
                                      // Info<< "tangential force [N] in the point section (with density included) " << F_tita_Bi*density_ << endl;

    // save for the output
    FnList.append(F_n_Bi * density_);
    FtList.append(F_tita_Bi * density_);
    UnList.append(U_n);
    UtList.append(U_t);

    // global from all the blades
    TorqueSects += mag(diskPoint_ - Bi) * F_tita_Bi * density_;
    Pcells += mag(diskPoint_ - Bi) * F_tita_Bi * density_ * omega;

    // Info<< "-----FORCES finished ----- " << endl;

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

    Info << "Done with force interpolation" << endl;

    //---colecting data from all the procesors
    reduce(Tcells, sumOp<scalar>());
    // Info<< "Tcells (with density) from all the sections: " << Tcells << endl;
    reduce(Ftitacells, sumOp<scalar>());
    // Info<< "Ftitacells (with density) from all the sections: " << Ftitacells << endl;
    reduce(U_infCells, sumOp<vector>());
    scalar U_infCellsMag = mag(U_infCells) / nodesNumber_;
    // Info<< "U_infCells from all the sections: " << U_infCellsMag << endl;
    // reduce(TorqueSects, sumOp<scalar>());
    ////Info<< "TorqueSects (with density) from all the sections: " << TorqueSects << endl;

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
}
// ************************************************************************* //
