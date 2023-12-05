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
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuationDiskRingsV11_calib_Source, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuationDiskRingsV11_calib_Source,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::actuationDiskRingsV11_calib_Source::checkData() const
{

    if (magSqr(diskArea_) <= VSMALL)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV11_calib_Source::checkData()")
           << "diskArea is approximately zero"
           << exit(FatalIOError);
    }
    if (Cp_ <= VSMALL || Ct_ <= VSMALL)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV11_calib_Source::checkData()")
           << "Cp and Ct must be greater than zero"
           << exit(FatalIOError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV11_calib_Source::checkData()")
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }

    if (returnReduce(diskCellId_, maxOp<label>()) == -1)
    {
        FatalErrorIn("Foam::fv::actuationDiskRingsV11_calib_Source::checkData()")
           << "disk center location " << diskPoint_  << " not found in mesh"
           << exit(FatalIOError);
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuationDiskRingsV11_calib_Source::actuationDiskRingsV11_calib_Source
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    //option(name, modelType, dict, mesh), //---new line 
    diskDir_(coeffs_.lookup("diskDir")),
    Cp_(readScalar(coeffs_.lookup("Cp"))),
    Ct_(readScalar(coeffs_.lookup("Ct"))),
    Uref_(readScalar(coeffs_.lookup("Uinf"))),
    cellSize_(readScalar(coeffs_.lookup("cellSize",1000))),
    yaw_(readScalar(coeffs_.lookup("yaw"))),
    omega_(readScalar(coeffs_.lookup("omega",0.0))),
    diskArea_(readScalar(coeffs_.lookup("diskArea"))),
	diskPoint_(coeffs_.lookup("diskPoint")),
	rootFactor_(readScalar(coeffs_.lookup("rootFactor"))),
	rootDistance_(readScalar(coeffs_.lookup("rootDistance"))),
	tipFactor_(readScalar(coeffs_.lookup("tipFactor"))),
	nodesCellsRatio_(readScalar(coeffs_.lookup("nodesCellsRatio"))),
	rThicknessCellsizeRatio_(readScalar(coeffs_.lookup("rThicknessCellsizeRatio"))),
	density_(readScalar(coeffs_.lookup("density"))),
	gradInterpolation_(readScalar(coeffs_.lookup("gradInterpolation"))),
	diskCellId_(-1)
{
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuationDiskRingsV11_calib_Source: "
        << this->name() << endl;
	diskCellId_ = mesh.findCell(diskPoint_);


	//---TEST--------------------------------------------------
	//all the variables should be declarade in .H in order to be seen here and i Template

	//calculate the radius
	maxR_ = sqrt(diskArea_ / M_PI);

    //calculate number of cells in the AD
    cellsInAd_ = round(diskArea_/pow(cellSize_,2));
    estimatedNodes_ = nodesCellsRatio_*cellsInAd_,

    Info << "estimated number of nodes: "<<estimatedNodes_<<endl;

    //calculate area corresponding to each node
    secArea_ = diskArea_/estimatedNodes_;
    
    //calculate rInt_
    rInt_ = sqrt(secArea_ / M_PI);
    rExt_ = maxR_;

    //calculate estimate ring Thickness according to cellsize
    ringThickness_ = rThicknessCellsizeRatio_*cellSize_;
    Info << "original ringThickness_: "<<ringThickness_<<endl;

    //calculate number of rings
    numberRings_ = round((rExt_-rInt_)/ringThickness_);
    Info << "number of rings: "<<numberRings_<<endl;

    //correction of the ringThickness_
    ringThickness_ =(rExt_-rInt_)/numberRings_;
    Info << "new ringThickness_: "<<ringThickness_<<endl;

    //Total nodes in the AD
    nodesNumber_ = 0;

    //Nodes in each ring
    nodesI_ = 0;

    //Angle between nodes - deg
    titaI_ = 0;

    //Radius of each ring
    rMedI_ = rInt_ + ringThickness_/2;

    //Area of each ring
    areaI_ = 0;

    Info << "calculate number of rings" << endl;

    //for each ring of nodes
	for (int ring =1; ring<=numberRings_; ring=ring+1)
    {
        nodesI_ = estimatedNodes_*(pow(rMedI_+ringThickness_/2,2)-pow(rMedI_-ringThickness_/2,2))/(pow(rExt_,2)-pow(rInt_,2));
        ringNodesList_.append(round(nodesI_));
        nodesNumber_ += round(nodesI_);
        titaI_ = 360/round(nodesI_); //deg
        ringTitaList_.append(titaI_);
        ringrMedList_.append(rMedI_);
        areaI_ = M_PI*(pow(rMedI_+ringThickness_/2,2)-pow(rMedI_-ringThickness_/2,2))/round(nodesI_);
        ringAreaList_.append(areaI_);
        Info << "ring: "<<ring;
        Info << " - nodes: "<<round(nodesI_);
        Info << " - tita: "<<titaI_;
        Info << " - node area: "<<areaI_<<endl;
        rMedI_ += ringThickness_;
        }
    
    //add the center node
    ringNodesList_.append(1);
    ringTitaList_.append(0);
    ringrMedList_.append(0);
    ringAreaList_.append(secArea_);    
    nodesNumber_ = nodesNumber_ + 1;

    Info << "ringNodesList_: "<<ringNodesList_<<endl;
    Info << "ringTitaList_: "<<ringTitaList_<<endl;
    Info << "ringrMedList_: "<<ringrMedList_<<endl;
    Info << "ringAreaList_: "<<ringAreaList_<<endl;


    Info << "final nodes number: "<<nodesNumber_<<endl;
    nodeCellID_.setSize(nodesNumber_);

	//direction of the disc
	yawRad_=yaw_*2*M_PI/360;

	//rotate the orginal diskDir with the yaw angle
	diskYawed_ = vector(diskDir_[0]*cos(yawRad_)-diskDir_[1]*sin(yawRad_),diskDir_[0]*sin(yawRad_)+diskDir_[1]*cos(yawRad_),diskDir_[2]);
	uniDiskDir_=diskYawed_;
    
    tita_r_ = 0;//tita ring
    rMed_r_ = 0;//radius ring

    total_nodes_counter_ = 0;

    Info << "Calculate nodes positions" << endl;

	//for each ring
	for (int ring =0; ring<=numberRings_; ring=ring+1)
    {
        //Info << "ring: "<< ring << endl;
        //Info << "nodes in ring: "<< ringNodesList_[ring] << endl;
        tita_r_ = ringTitaList_[ring];
        //Info << "tita_r_: "<< tita_r_ << endl;
        rMed_r_ = ringrMedList_[ring];
        //Info << "rMed_r_: "<< rMed_r_ << endl;


        for (int nodeIterator = 1; nodeIterator <= ringNodesList_[ring]; nodeIterator+=1)
        {
		//Info << "node in ring: "<< nodeIterator << endl;
		tita_n_Rad_=2*M_PI*(tita_r_*(nodeIterator-1))/360;  
		
        //Info << "tita_n_Rad_: "<< tita_n_Rad_<< endl;
        //Info << "yawRad_: "<< yawRad_<< endl;

        //position of the node
        scalar x_node_ = -1*rMed_r_*sin(tita_n_Rad_)*sin(yawRad_);
        scalar y_node_ = rMed_r_*sin(tita_n_Rad_)*cos(yawRad_);
        scalar z_node_ = rMed_r_*cos(tita_n_Rad_);

        //IMPORTANTE PARA EL CALCULO DE LA POSICION HICIMOS ALGO ANÁLOGO A COORDENADAS CILINDRICAS
        //Y PARA QUE ESTO VALGA SI O SI SE DEBE CUMPLIR QUE EL YAW ANGLE VAYA ENTRE 0 Y 180 GRADOS. 
        //VAMOS CON ESA. 
        
        //move to turbine position
		x_node_=x_node_+diskPoint_[0];
		y_node_=y_node_+diskPoint_[1];
		z_node_=z_node_+diskPoint_[2];

		vector Bi_= vector(x_node_,y_node_,z_node_);

        //Info << "Bi_: "<< Bi_ << endl;
        
		//fin the closest cell in all the domain
		diskCellId_ = mesh.findCell(Bi_);
		//reduce(diskCellId_ , maxOp<label>());

		nodeCellID_[total_nodes_counter_]=diskCellId_;

        reduce(diskCellId_ , maxOp<label>());
        //Info << "diskCellId_: "<< diskCellId_ << endl;
        total_nodes_counter_ += 1;


		}//close node loop
	}//close ring loop
	Info<<"total_nodes_counter: "<<total_nodes_counter_<<endl;
	
	Info<<"nodeCellID: "<<nodeCellID_<<endl;

	//------------------------------------------------------------






    checkData();

	outTurbines = new std::ofstream("outTurbines.csv", std::ios::app);
	if (Pstream::myProcNo() == 0) //In Master node
	    {
	(*outTurbines) << "---ADRings.calibration---" << std::endl; //Head file
	(*outTurbines) << "0-Turbine, 1-time [s], 2-Uinf(fixed) [m/s], 3-Cp(fixed), 4-Ct(fixed), 5-Power(Uinf,Cp) [W], 6-Thrust(Uinf,Ct) [N] , 7-Center Ud [m/s], 8-Average Ud [m/s], 9-Uinf(disc cells) [m/s], 10-Power(disc cells) [W], 11-Thrust(disc cells) [N], 12-Torque(disc cells) [Nm]" << std::endl; //Head file
	}


	outRings = new std::ofstream("outRings.csv", std::ios::app);
	if (Pstream::myProcNo() == 0) //In Master node
	    {
	(*outRings) << "---ADRings.calibration---" << std::endl; //Head file
	}

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::actuationDiskRingsV11_calib_Source::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            Zero
            //vector::zero //---new line
        )
    );

    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        addactuationDiskRingsV11_calib_AxialInertialResistance
        (
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U,
            force
        );
    }

    if (mesh_.time().writeTime())
    //if (mesh_.time().outputTime()) //---new line
    {
        force.write();
    }

}


void Foam::fv::actuationDiskRingsV11_calib_Source::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            //vector::zero //---new line
            Zero
        )
    );


    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        addactuationDiskRingsV11_calib_AxialInertialResistance
        (
            Usource,
            cells_,
            cellsV,
            rho,
            U,
            force
        );
    }

    if (mesh_.time().writeTime())
    //if (mesh_.time().outputTime()) //---new line
    {
        force.write();
    }

}


bool Foam::fv::actuationDiskRingsV11_calib_Source::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    //if (option::read(dict)) //---new line
    {
        coeffs_.readIfPresent("diskDir", diskDir_);
        coeffs_.readIfPresent("Cp", Cp_);
        coeffs_.readIfPresent("Ct", Ct_);
	    coeffs_.readIfPresent("Uinf", Uref_);
        coeffs_.readIfPresent("yaw", yaw_);
        coeffs_.readIfPresent("omega", omega_);
        coeffs_.readIfPresent("pitch", pitch_);
        coeffs_.readIfPresent("cellSize", cellSize_);
        coeffs_.readIfPresent("diskArea", diskArea_);
        coeffs_.readIfPresent("diskPoint", diskPoint_);
        coeffs_.readIfPresent("rootFactor", rootFactor_);
        coeffs_.readIfPresent("rootDistance", rootDistance_);
        coeffs_.readIfPresent("tipFactor", tipFactor_);
        coeffs_.readIfPresent("nodesCellsRatio", nodesCellsRatio_);
        coeffs_.readIfPresent("rThicknessCellsizeRatio", rThicknessCellsizeRatio_);
        coeffs_.readIfPresent("density", density_);
        coeffs_.readIfPresent("gradInterpolation", gradInterpolation_);


        checkData();

        return true;
    }
    else
    {
        return false;
    }
}

Foam::fv::actuationDiskRingsV11_calib_Source::~actuationDiskRingsV11_calib_Source()
{
	(*outTurbines).close();
	(*outRings).close();
}


// ************************************************************************* //

