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
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(actuationDiskRingsV21_Source, 0);
        addToRunTimeSelectionTable(
            option,
            actuationDiskRingsV21_Source,
            dictionary);
    }
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::actuationDiskRingsV21_Source::checkData() const
{

    if (magSqr(diskArea_) <= VSMALL)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV21_Source::checkData()")
            << "diskArea is approximately zero"
            << exit(FatalIOError);
    }
    if (Cp_ <= VSMALL || Ct_ <= VSMALL)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV21_Source::checkData()")
            << "Cp and Ct must be greater than zero"
            << exit(FatalIOError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV21_Source::checkData()")
            << "disk direction vector is approximately zero"
            << exit(FatalIOError);
    }

    if (returnReduce(diskCellId_, maxOp<label>()) == -1)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV21_Source::checkData()")
            << "disk center location " << diskPoint_ << " not found in mesh"
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuationDiskRingsV21_Source::actuationDiskRingsV21_Source(
    const word &name,
    const word &modelType,
    const dictionary &dict,
    const fvMesh &mesh)
    : cellSetOption(name, modelType, dict, mesh),
      // option(name, modelType, dict, mesh), //---new line
      diskDir_(coeffs_.lookup("diskDir")),
      Cp_(readScalar(coeffs_.lookup("Cp"))),
      Ct_(readScalar(coeffs_.lookup("Ct"))),
      Uref_(readScalar(coeffs_.lookup("Uinf"))),
      cellSize_(readScalar(coeffs_.lookup("cellSize", 1000))),
      yaw_(readScalar(coeffs_.lookup("yaw"))),
      omega_(readScalar(coeffs_.lookup("omega", 0.0))),
      lambda_(readScalar(coeffs_.lookup("lambda"))),
      diskArea_(readScalar(coeffs_.lookup("diskArea"))),
      centerRatio_(readScalar(coeffs_.lookup("centerRatio"))),
      rootDistance_(readScalar(coeffs_.lookup("rootDistance"))),
      diskPoint_(coeffs_.lookup("diskPoint")),
      UdCellsMethod_(readScalar(coeffs_.lookup("UdCellsMethod"))),
      UdCenterToggle_(readScalar(coeffs_.lookup("UdCenterToggle"))),
      forceDistributionMethod_(readScalar(coeffs_.lookup("forceDistributionMethod"))),
      ADmodel_(readScalar(coeffs_.lookup("ADmodel"))),
      rootFactor_(readScalar(coeffs_.lookup("rootFactor"))),
      tipFactor_(readScalar(coeffs_.lookup("tipFactor"))),
      Ct_rated_(readScalar(coeffs_.lookup("Ct_rated"))),
      averageChordLength_(readScalar(coeffs_.lookup("averageChordLength"))),
      nodesCellsRatio_(readScalar(coeffs_.lookup("nodesCellsRatio"))),
      rThicknessCellsizeRatio_(readScalar(coeffs_.lookup("rThicknessCellsizeRatio"))),
      powerCurve_table_(coeffs_.lookup("powerCurve_table")), 
      UdAvg_table_(coeffs_.lookup("UdAvg_table")), // for adaptation
      Udi_table_(coeffs_.lookup("Udi_table")),     // for adaptation
      gradInterpolation_(readScalar(coeffs_.lookup("gradInterpolation"))),
      diskCellId_(-1)
{
    //---define list using the table UdAvg_table
    UrefPowerCurveList_.setSize(powerCurve_table_.size());
    CtPowerCurveList_.setSize(powerCurve_table_.size());
    CpPowerCurveList_.setSize(powerCurve_table_.size());

    forAll(powerCurve_table_, i)
    {
        UrefPowerCurveList_[i] = powerCurve_table_[i][0];
        CtPowerCurveList_[i] = powerCurve_table_[i][1];
        CpPowerCurveList_[i] = powerCurve_table_[i][2];
    }

    if (ADmodel_ == 1)
    {
        //---define list using the table UdAvg_table
        UrefList_.setSize(UdAvg_table_.size());
        omegaList_.setSize(UdAvg_table_.size());
        pitchList_.setSize(UdAvg_table_.size());
        UdAvgList_.setSize(UdAvg_table_.size());
        CtList_.setSize(UdAvg_table_.size());
        CpList_.setSize(UdAvg_table_.size());

        forAll(UdAvg_table_, i)
        {
            UrefList_[i] = UdAvg_table_[i].first()[0];
            omegaList_[i] = UdAvg_table_[i].first()[1];
            pitchList_[i] = UdAvg_table_[i].first()[2];
            UdAvgList_[i] = UdAvg_table_[i].second()[0];
            CtList_[i] = UdAvg_table_[i].second()[1];
            CpList_[i] = UdAvg_table_[i].second()[2];
        }

        //---define list using the table Udi_table

        //- Original List
        Uref2List_orig.setSize(Udi_table_.size());
        UinfList_orig.setSize(Udi_table_.size());
        rList_orig.setSize(Udi_table_.size());
        UdiList_orig.setSize(Udi_table_.size());
        fnList_orig.setSize(Udi_table_.size());
        ftList_orig.setSize(Udi_table_.size());

        forAll(Udi_table_, i)
        {
            Uref2List_orig[i] = Udi_table_[i].first()[0];
            UinfList_orig[i] = Udi_table_[i].first()[1];
            rList_orig[i] = Udi_table_[i].first()[2];
            UdiList_orig[i] = Udi_table_[i].second()[0];
            fnList_orig[i] = Udi_table_[i].second()[1];
            ftList_orig[i] = Udi_table_[i].second()[2];
        }

        Info << "Print Avg Table: " << endl;
        Info << UdAvg_table_ << endl;
    }

    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info << "    - creating actuationDiskRingsV21_Source: "
         << this->name() << endl;
    diskCellId_ = mesh.findCell(diskPoint_);


    //---TEST--------------------------------------------------
    // all the variables should be declarade in .H in order to be seen here and i Template

    // calculate the radius
    maxR_ = sqrt(diskArea_ / M_PI);

    // calculate number of cells in the AD
    cellsInAd_ = round(diskArea_ / pow(cellSize_, 2));
    estimatedNodes_ = nodesCellsRatio_ * cellsInAd_,

    Info << "estimated number of nodes: " << estimatedNodes_ << endl;

    // calculate area corresponding to each node
    secArea_ = diskArea_ / estimatedNodes_;

    // calculate rInt_
    rInt_ = sqrt(secArea_ / M_PI);
    rExt_ = maxR_;

    // calculate estimate ring Thickness according to cellsize
    ringThickness_ = rThicknessCellsizeRatio_ * cellSize_;
    Info << "original ringThickness_: " << ringThickness_ << endl;

    // calculate number of rings
    numberRings_ = round((rExt_ - rInt_) / ringThickness_);
    Info << "number of rings: " << numberRings_ << endl;

    // correction of the ringThickness_
    ringThickness_ = (rExt_ - rInt_) / numberRings_;
    Info << "new ringThickness_: " << ringThickness_ << endl;

    // Total nodes in the AD
    nodesNumber_ = 0;

    // Nodes in each ring
    nodesI_ = 0;

    // Angle between nodes - deg
    titaI_ = 0;

    // Radius of each ring
    rMedI_ = rInt_ + ringThickness_ / 2;

    // Area of each ring
    areaI_ = 0;

    Info << "calculate number of rings" << endl;

    // for each ring of nodes
    for (int ring = 1; ring <= numberRings_; ring = ring + 1)
    {
        nodesI_ = estimatedNodes_ * (pow(rMedI_ + ringThickness_ / 2, 2) - pow(rMedI_ - ringThickness_ / 2, 2)) / (pow(rExt_, 2) - pow(rInt_, 2));
        ringNodesList_.append(round(nodesI_));
        nodesNumber_ += round(nodesI_);
        titaI_ = 360 / round(nodesI_); // deg
        ringTitaList_.append(titaI_);
        ringrMedList_.append(rMedI_);
        areaI_ = M_PI * (pow(rMedI_ + ringThickness_ / 2, 2) - pow(rMedI_ - ringThickness_ / 2, 2)) / round(nodesI_);
        ringAreaList_.append(areaI_);
        // Info << "ring: " << ring << endl;
        // Info << " - nodes: " << round(nodesI_) << endl;
        // Info << " - tita: " << titaI_ << endl;
        // Info << " - node area: " << areaI_ << endl;
        rMedI_ += ringThickness_;
    }

    int nR_orig = 1;
    if (ADmodel_ == 1)
    {
        // count how many radius sections are in the original table
        while (rList_orig[nR_orig] - rList_orig[nR_orig - 1] > 0)
        {
            nR_orig = nR_orig + 1;
        }
        Info << "nR_orig: " << nR_orig << endl;

        Info << "rList_orig: " << endl;
        for (int i = 0; i < nR_orig; i = i + 1)
        {
            Info << rList_orig[i] << endl;
        }
    }

    // Set up a list with node's radial positions
    rNodeList_.append(0); // First we add the center node
    for (int i = 0; i < numberRings_; i = i + 1)
    {
        rNodeList_.append(ringrMedList_(i));
    }

    Info << "rNodeList_: " << rNodeList_ << endl;

    if (ADmodel_ == 1)
    {
        // Set up a list with Uinf from original table
        for (int i = 0; i < (UinfList_orig.size() / UrefList_.size()); i = i + nR_orig)
        {
            UinfOnlyList_.append(UinfList_orig[i]);
        }

        Info << "Original UinfOnlyList_: " << UinfOnlyList_ << endl;

        // Set up new table lists with the corresponding size
        Uref2List_.setSize(rNodeList_.size() * UrefList_.size() * UinfOnlyList_.size());
        UinfList_.setSize(rNodeList_.size() * UrefList_.size() * UinfOnlyList_.size());
        rList_.setSize(rNodeList_.size() * UrefList_.size() * UinfOnlyList_.size());
        UdiList_.setSize(rNodeList_.size() * UrefList_.size() * UinfOnlyList_.size());
        fnList_.setSize(rNodeList_.size() * UrefList_.size() * UinfOnlyList_.size());
        ftList_.setSize(rNodeList_.size() * UrefList_.size() * UinfOnlyList_.size());

        // Set up new list with closer radial nodes from the old table.
        posrList_.setSize(rNodeList_.size());

        for (int i = 0; i < rNodeList_.size(); i = i + 1)
        {
            scalar dist = VGREAT;
            scalar posr = 0;
            while ((mag(rNodeList_[i] - rList_orig[posr]) < dist) && (posr < nR_orig - 1))
            {
                dist = mag(rNodeList_[i] - rList_orig[posr]);

                // if ((rNodeList_[i] - rList_orig[posr + 1]) >= 0)
                if (
                    mag(rNodeList_[i] - rList_orig[posr + 1]) < dist and
                    ((rNodeList_[i] - rList_orig[posr + 1]) >= 0)
                )
                {
                    posr += 1;
                }
            }

            posrList_[i] = posr;
        }

        Info << "posrList_ " << posrList_ << endl;

        // we fill up the new table
        for (int i = 0; i < UrefList_.size(); i = i + 1)
        {
            for (int j = 0; j < UinfOnlyList_.size(); j = j + 1)
            {
                // for (int k = 1; k < rNodeList_.size(); k = k + 1) // La posicion 0 la agregamos al final
                for (int k = 0; k < rNodeList_.size(); k = k + 1) 
                {
                    // Info << "List position: "<<i*(UinfOnlyList_.size()*rNodeList_.size())+j*(rNodeList_.size())+k << endl;
                    Uref2List_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = UrefList_[i];
                    UinfList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = UinfOnlyList_[j];
                    rList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = rNodeList_[k];
                    UdiList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = UdiList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]] + (rNodeList_[k] - rList_orig[posrList_[k]]) * (UdiList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k] + 1] - UdiList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]]) / (rList_orig[posrList_[k] + 1] - rList_orig[posrList_[k]]);
                    fnList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = fnList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]] + (rNodeList_[k] - rList_orig[posrList_[k]]) * (fnList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k] + 1] - fnList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]]) / (rList_orig[posrList_[k] + 1] - rList_orig[posrList_[k]]);
                    ftList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = ftList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]] + (rNodeList_[k] - rList_orig[posrList_[k]]) * (ftList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k] + 1] - ftList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]]) / (rList_orig[posrList_[k] + 1] - rList_orig[posrList_[k]]);
                }
            }
        }

        // we fill up the new table
        // for (int i = 0; i < UrefList_.size(); i = i + 1)
        // {
        //     for (int j = 0; j < UinfOnlyList_.size(); j = j + 1)
        //     {
        //         int k = 0;
        //         Uref2List_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = UrefList_[i];
        //         UinfList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = UinfOnlyList_[j];
        //         rList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = rNodeList_[k];

        //         UdiList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = UdiList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]] + (rNodeList_[k] - rList_orig[posrList_[k]]) * (UdiList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k] + 1] - UdiList_orig[i * (UinfOnlyList_.size() * nR_orig) + j * nR_orig + posrList_[k]]) / (rList_orig[posrList_[k] + 1] - rList_orig[posrList_[k]]);

        //         fnList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = fnList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k + rNodeList_.size() - 1] / ringNodesList_[ringNodesList_.size() - 1];
        //         ftList_[i * (UinfOnlyList_.size() * rNodeList_.size()) + j * (rNodeList_.size()) + k] = 0;
        //     }
        // }
    }

    // add the center node
    ringNodesList_.append(1);
    ringTitaList_.append(0);
    ringrMedList_.append(0);
    ringAreaList_.append(secArea_);
    nodesNumber_ = nodesNumber_ + 1;

    Info << "ringNodesList_: " << ringNodesList_ << endl;
    Info << "ringTitaList_: " << ringTitaList_ << endl;
    Info << "ringrMedList_: " << ringrMedList_ << endl;
    Info << "ringAreaList_: " << ringAreaList_ << endl;

    Info << "final nodes number: " << nodesNumber_ << endl;
    nodeCellID_.setSize(nodesNumber_);

    // direction of the disc
    yawRad_ = yaw_ * 2 * M_PI / 360;

    // rotate the orginal diskDir with the yaw angle
    diskYawed_ = vector(diskDir_[0] * cos(yawRad_) - diskDir_[1] * sin(yawRad_), diskDir_[0] * sin(yawRad_) + diskDir_[1] * cos(yawRad_), diskDir_[2]);
    uniDiskDir_ = diskYawed_; // REVISAMOS ESTA CUENTA DE GONZA Y ES CORRECTA

    tita_r_ = 0; // tita ring
    rMed_r_ = 0; // radius ring

    total_nodes_counter_ = 0;

    // Info << "Calculate nodes positions" << endl;

    // for each ring
    for (int ring = 0; ring <= numberRings_; ring = ring + 1)
    {
        // Info << "ring: "<< ring << endl;
        // Info << "nodes in ring: "<< ringNodesList_[ring] << endl;
        tita_r_ = ringTitaList_[ring];
        // Info << "tita_r_: "<< tita_r_ << endl;
        rMed_r_ = ringrMedList_[ring];
        // Info << "rMed_r_: "<< rMed_r_ << endl;

        for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator += 1)
        {
            // Info << "node in ring: "<< nodeIterator << endl;

            tita_n_Rad_ = 2 * M_PI * (tita_r_ * (nodeIterator - 1)) / 360;

            // Info << "tita_n_Rad_: "<< tita_n_Rad_<< endl;
            // Info << "yawRad_: "<< yawRad_<< endl;

            // position of the node
            scalar x_node_ = -1 * rMed_r_ * sin(tita_n_Rad_) * sin(yawRad_);
            scalar y_node_ = rMed_r_ * sin(tita_n_Rad_) * cos(yawRad_);
            scalar z_node_ = rMed_r_ * cos(tita_n_Rad_);

            // IMPORTANTE PARA EL CALCULO DE LA POSICION HICIMOS ALGO ANÁLOGO A COORDENADAS CILINDRICAS
            // Y PARA QUE ESTO VALGA SI O SI SE DEBE CUMPLIR QUE EL YAW ANGLE VAYA ENTRE 0 Y 180 GRADOS.
            // VAMOS CON ESA.

            // move to turbine position
            x_node_ = x_node_ + diskPoint_[0];
            y_node_ = y_node_ + diskPoint_[1];
            z_node_ = z_node_ + diskPoint_[2];

            vector Bi_ = vector(x_node_, y_node_, z_node_);

            // Info << "Bi_: "<< Bi_ << endl;

            // fin the closest cell in all the domain
            diskCellId_ = mesh.findCell(Bi_);
            // reduce(diskCellId_ , maxOp<label>());

            nodeCellID_[total_nodes_counter_] = diskCellId_;

            reduce(diskCellId_, maxOp<label>());
            // Info << "diskCellId_: "<< diskCellId_ << endl;
            total_nodes_counter_ += 1;

        } // close node loop
    }     // close ring loop
    Info << "total_nodes_counter: " << total_nodes_counter_ << endl;

    //------------------------------------------------------------

    checkData();

    outTurbines = new std::ofstream("outTurbines.csv", std::ios::app);
    if (Pstream::myProcNo() == 0) // In Master node
    {
        (*outTurbines) << "---ADRings.calibration---" << std::endl; // Head file
        (*outTurbines) << "0-Turbine, 1-time [s], 2-Uinf(fixed) [m/s], 3-Cp(fixed), 4-Ct(fixed), 5-Power(Uinf,Cp) [W], 6-Thrust(Uinf,Ct) [N] , 7-Center Ud [m/s], 8-Average Ud [m/s], 9-Uinf(disc cells) [m/s], 10-Power(disc cells) [W], 11-Thrust(disc cells) [N], 12-Torque(disc cells) [Nm]" << std::endl;
    }

    outRings = new std::ofstream("outRings.csv", std::ios::app);
    if (Pstream::myProcNo() == 0) // In Master node
    {
        (*outRings) << "---ADRings.calibration---" << std::endl; // Head file
    }

    if (ADmodel_ == 1)
    {
        outListFvOptions = new std::ofstream("outListFvOptions.txt", std::ios::app);
        if (Pstream::myProcNo() == 0) // In Master node
        {
            (*outListFvOptions) << "UdAvg_table" << std::endl;
            (*outListFvOptions) << "(" << std::endl;
            (*outListFvOptions) << "//(Uref[m/s], omega[rad/s], empty)(	UdAvg[m/s], 	Ct, 	Cp)" << std::endl;

            for (int i = 0; i < UdAvg_table_.size(); i = i + 1)
            {
                (*outListFvOptions) << "((" << UrefList_[i] << " " << omegaList_[i] << " 0) (" << UdAvgList_[i] << " " << CtList_[i] << " " << CpList_[i] << "))" << std::endl;
            }

            (*outListFvOptions) << ");" << std::endl;
            (*outListFvOptions) << "Udi_table" << std::endl;
            (*outListFvOptions) << "(" << std::endl;
            (*outListFvOptions) << "//(Uref[m/s], Uinf[m/s], r[m])(Udi(r)[m/s], fn(r)[N/m], ft(r)[N/m])" << std::endl;

            (*outListFvOptions) << ");" << std::endl;

            for (int i = 0; i < Uref2List_.size(); i = i + 1)
            {
                (*outListFvOptions) << "((" << Uref2List_[i] << " " << UinfList_[i] << " " << rList_[i] << ") (" << UdiList_[i] << " " << fnList_[i] << " " << ftList_[i] << "))" << std::endl;
            }

            (*outListFvOptions) << ");" << std::endl;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar UrefPrevious;
// void Foam::fv::actuationDiskRingsV21_Source::addSup(
void Foam::fv::actuationDiskRingsV21_Source::addSup(
    fvMatrix<vector> &eqn,
    const label fieldI)
{
    volVectorField force(
        IOobject(
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_),
        mesh_,
        dimensionedVector(
            "zero",
            eqn.dimensions() / dimVolume,
            // vector::zero //---new line
            Zero));

    const scalarField &cellsV = mesh_.V();
    vectorField &Usource = eqn.source();
    const vectorField &U = eqn.psi();
    Info << "time = " << mesh().time().value() << endl;
    // scalar UrefPrevious;
    if (mesh().time().value() == 1)
    {
        // UrefPrevious = Uinf() * 1.1;
        UrefPrevious = Uinf();
    }
    Info << "UrefPrevious = " << UrefPrevious << endl;
    // scalar UrefPrevious = Uinf();

    if (V() > VSMALL)
    {
        // addactuationDiskRings_AxialInertialResistance(
        // scalar Uref = addactuationDiskRings_AxialInertialResistance(
        UrefPrevious = addactuationDiskRings_AxialInertialResistance(
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U,
            // force);
            force,
            UrefPrevious);
        // UrefPrevious = Uref;
    }

    if (mesh_.time().writeTime())
    // if (mesh_.time().outputTime()) //---new line
    {
        force.write();
    }
}

void Foam::fv::actuationDiskRingsV21_Source::addSup(
    const volScalarField &rho,
    fvMatrix<vector> &eqn,
    const label fieldI)
{
    volVectorField force(
        IOobject(
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_),
        mesh_,
        dimensionedVector(
            "zero",
            eqn.dimensions() / dimVolume,
            // vector::zero //---new line
            Zero));

    const scalarField &cellsV = mesh_.V();
    vectorField &Usource = eqn.source();
    const vectorField &U = eqn.psi();
    Info << "time = " << mesh().time().value() << endl;
    // scalar UrefPrevious;
    if (mesh().time().value() == 1)
    {
        // UrefPrevious = Uinf() * 1.1;
        UrefPrevious = Uinf();
    }
    Info << "UrefPrevious = " << UrefPrevious << endl;
    // scalar UrefPrevious = Uinf();

    if (V() > VSMALL)
    {
        // addactuationDiskRings_AxialInertialResistance(
        // scalar Uref = addactuationDiskRings_AxialInertialResistance(
        UrefPrevious = addactuationDiskRings_AxialInertialResistance(
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U,
            // force);
            force,
            UrefPrevious);
        // UrefPrevious = Uref;
    }

    if (mesh_.time().writeTime())
    // if (mesh_.time().outputTime()) //---new line
    {
        force.write();
    }
}

bool Foam::fv::actuationDiskRingsV21_Source::read(const dictionary &dict)
{
    if (cellSetOption::read(dict))
    // if (option::read(dict)) //---new line
    {
        coeffs_.readIfPresent("diskDir", diskDir_);
        coeffs_.readIfPresent("Cp", Cp_);
        coeffs_.readIfPresent("Ct", Ct_);
        coeffs_.readIfPresent("Uinf", Uref_);
        coeffs_.readIfPresent("yaw", yaw_);
        coeffs_.readIfPresent("omega", omega_);
        coeffs_.readIfPresent("lambda", lambda_);
        coeffs_.readIfPresent("cellSize", cellSize_);
        coeffs_.readIfPresent("diskArea", diskArea_);
        coeffs_.readIfPresent("diskPoint", diskPoint_);
        coeffs_.readIfPresent("rootFactor", rootFactor_);
        coeffs_.readIfPresent("tipFactor", tipFactor_);
        coeffs_.readIfPresent("nodesCellsRatio", nodesCellsRatio_);
        coeffs_.readIfPresent("rThicknessCellsizeRatio", rThicknessCellsizeRatio_);
        coeffs_.readIfPresent("UdAvg_table", UdAvg_table_);
        coeffs_.readIfPresent("Udi_table", Udi_table_);
        coeffs_.readIfPresent("gradInterpolation", gradInterpolation_);
        coeffs_.readIfPresent("UdCellsMethod", UdCellsMethod_);
        coeffs_.readIfPresent("UdCenterToggle", UdCenterToggle_);
        coeffs_.readIfPresent("forceDistributionMethod", forceDistributionMethod_);
        coeffs_.readIfPresent("ADmodel", ADmodel_);
        coeffs_.readIfPresent("centerRatio", centerRatio_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}

Foam::fv::actuationDiskRingsV21_Source::~actuationDiskRingsV21_Source()
{
    (*outTurbines).close();
    (*outRings).close();
}

// ************************************************************************* //
