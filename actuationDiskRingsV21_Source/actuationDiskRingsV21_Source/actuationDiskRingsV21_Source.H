/*---------------------------------------------------------------------------*\
Authors: Juan Ignacio Teich & Dimas Alejandro Barile
CSC - CONICET, Buenos Aires

-----------------------------------------------------------------------------
OpenFOAM

Actuation Disk
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

#ifndef actuationDiskRingsV21_Source_H
#define actuationDiskRingsV21_Source_H

#include "cellSetOption.H"
// #include "fvOption.H" //---new line
#include "bladeModel.H"
#include "profileModelList.H"
#include <fstream>
#include <ostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
    namespace fv
    {
        /*---------------------------------------------------------------------------*\
                             Class actuationDiskRingsV21_Source Declaration
        \*---------------------------------------------------------------------------*/

        class actuationDiskRingsV21_Source
            // : public cellSetOption
            : public fv::cellSetOption
        // public option //---new line
        {

        protected:
            // Protected data

            //- Disk area normal
            vector diskDir_;

            //- Power coefficient
            scalar Cp_;

            //- Thrust coefficient
            scalar Ct_;

            //- Thrust coefficient for rated velocity
            scalar Ct_rated_;

            //- Rated velocity
            scalar U_rated_;

            //- Average chord length for optimum values in Gaussian force distribution 
            scalar averageChordLength_;

            //- Fixed Uinf for thrust calculation
            scalar Uref_;

            //- cell size in the disc
            scalar cellSize_;

            //- Fixed yaw angle compared to the inlet velocity vector
            scalar yaw_;

            //- Fixed rotational velocity [rad/s]
            scalar omega_;

            // - Tip speed ratio
            scalar lambda_;

            //- Fixed pitch angle [deg]
            scalar pitch_;

            //- Disk area
            scalar diskArea_;

            //- Disc center point sample
            point diskPoint_;

            //- Root factor on/off
            scalar rootFactor_;

            //- Tip factor on/off
            scalar tipFactor_;

            //- Ratio Total amount of nodes / CellsInAD
            scalar nodesCellsRatio_;

            // Ratio ring Thickness / cellSize
            scalar rThicknessCellsizeRatio_;

            // Interpolation by means of Grad(U).
            scalar gradInterpolation_;

            //- disk cell ID
            label diskCellId_;

            //--test--------------------------------------
            // the permanent variables have "_" in the end

            // for calculating number of rings and nodes

            // Root distance, % of R
            scalar rootDistance_;

            // AD radius
            scalar maxR_;

            // Cells in the AD
            scalar cellsInAd_;

            // Estimated nodes
            scalar estimatedNodes_;

            // Aproximate are associated to each node
            scalar secArea_;

            // Internal radius of the AD, inside there is only 1 node.
            scalar rInt_;

            // External radius of the AD, equal to maxR_
            scalar rExt_;

            // Thickness of all rings
            scalar ringThickness_;

            // Number of rings
            int numberRings_;

            // Total nodes in the AD
            scalar nodesNumber_;

            // Nodes in each ring
            scalar nodesI_;

            // Angle between nodes - deg
            scalar titaI_;

            // Radius of each ring
            scalar rMedI_;

            // sectional area of each ring
            scalar areaI_;

            // Lists of nodes, angle between nodes and radius of each ring.
            DynamicList<scalar> ringNodesList_;
            DynamicList<scalar> ringTitaList_;
            DynamicList<scalar> ringrMedList_;
            DynamicList<scalar> ringAreaList_;

            // direction of the disc
            scalar yawRad_;
            vector diskYawed_;
            vector uniDiskDir_;

            // For calculating node positions
            // angle between nodes in a certain ring
            scalar tita_r_;

            // angle for the position of a certain node
            scalar tita_n_Rad_;

            // radius for a certain ring
            scalar rMed_r_;

            // Nodes counter
            scalar total_nodes_counter_;

            // ratio to consider in center of AD
            scalar centerRatio_;

            // method to measure Velocity Ud
            //   - 0 = sum directly values in cells in AD
            //   - 1 = weight with Gaussian in distance to AD plane
            //   - 2 = weight with Gaussian in distance to AD plane + distance to center of AD
            int UdCellsMethod_;

            // toggle for measuring Ud with center of AD
            //   - 0 = use all of AD
            //   - 1 = use center of AD with centerRatio
            int UdCenterToggle_;

            int UdCorrection_;

            int forceDistributionMethod_;
            
            // AD model
            // - 0 : uniform
            // - 1 : numeric 
            // - 2 : analytic 
            // - 3 : analytic non-uniform inflow 
            // - 4 : generalized analytic non-uniform inflow 
            // - 5 : eliptic 
            int ADmodel_;

            // List of cell ID for each node
            List<label> nodeCellID_;

            //- power curve table. From turbine manufacturer data
            //
            List<vector> powerCurve_table_; //--- new line

            //- UrefPowerCurve
            List<scalar> UrefPowerCurveList_;

            //- CtPowerCurve
            List<scalar> CtPowerCurveList_;

            //- CpPowerCurve
            List<scalar> CpPowerCurveList_;

            //- UdAvg_table
            List<Tuple2<vector, vector>> UdAvg_table_;

            //- Uinf_vanDerLaan_table
            List<Tuple2<vector, vector>> Uinf_vanDerLaan_table_;

            //- UrefList
            List<scalar> UrefList_;

            //- omegaList
            List<scalar> omegaList_;

            //- pitchList
            List<scalar> pitchList_;

            //- UdAvgList
            List<scalar> UdAvgList_;

            //- CtList
            List<scalar> CtList_;

            //- CpList
            List<scalar> CpList_;

            //- Udi_table
            List<Tuple2<vector, vector>> Udi_table_; 

            //- forces_vanDerLaan_table
            List<Tuple2<vector, vector>> forces_vanDerLaan_table_; 

            //- Original lists
            //- Uref2List original
            List<scalar> Uref2List_orig;

            //- UinfList original
            List<scalar> UinfList_orig;

            //- rList original
            List<scalar> rList_orig;

            //- UdiList original
            List<scalar> UdiList_orig;

            //- fnList original
            List<scalar> fnList_orig;

            //- ftList original
            List<scalar> ftList_orig;

            //- Interpolated lists
            //- Uref2List
            List<scalar> Uref2List_;

            //- UinfList
            List<scalar> UinfList_;

            //- rList
            List<scalar> rList_;

            //- UdiList
            List<scalar> UdiList_;

            //- fnList
            List<scalar> fnList_;

            //- ftList
            List<scalar> ftList_;

            //- list for r positions in the original table
            List<scalar> posrList_;

            //- list por node radial positions
            DynamicList<scalar> UinfOnlyList_;

            //- Uinf legnth in table
            scalar UinfLength_;

            //- list por node radial positions
            DynamicList<scalar> rNodeList_;

            //--------------------------------------------

        private:
            // Private Member Functions

            //- Check data
            void checkData() const;
            
            // scalar UrefPrevious = 10;

            //- Add resistance to the UEqn
            template <class RhoFieldType>
            // void addactuationDiskRings_AxialInertialResistance(
            scalar addactuationDiskRings_AxialInertialResistance(
                vectorField &Usource,
                const labelList &cells,
                const scalarField &V,
                const RhoFieldType &rho,
                const vectorField &U,
                // vectorField &force) const;
                vectorField &force,
                // scalar &UrefPrevious) const;
                scalar &CtPrevious) const;

            //- Disallow default bitwise copy construct
            actuationDiskRingsV21_Source(const actuationDiskRingsV21_Source &);

            //- Disallow default bitwise assignment
            void operator=(const actuationDiskRingsV21_Source &);

            //- turbines output file
            std::ofstream *outTurbines;

            //- blades output file
            std::ofstream *outRings;

            //- New Calibration table output
            std::ofstream *outListFvOptions;

        public:

            //- Runtime type information
            TypeName("actuationDiskRingsV21_Source");

            // Constructors

            //- Construct from components
            actuationDiskRingsV21_Source(
                const word &name,
                const word &modelType,
                const dictionary &dict,
                const fvMesh &mesh);

            //- Destructor
            ~actuationDiskRingsV21_Source();

            // Member Functions

            // Access
            //- Return Cp
            scalar Cp() const
            {
                return Cp_;
            }

            //- Return Ct
            scalar Ct() const
            {
                return Ct_;
            }

            //- Return Ct_rated
            scalar Ct_rated() const
            {
                return Ct_rated_;
            }

            //- Return U_rated 
            scalar U_rated() const
            {
                return U_rated_;
            }

            //- Return averageChordLength 
            scalar averageChordLength() const
            {
                return averageChordLength_;
            }

            //- Return Uinf
            scalar Uinf() const
            {
                return Uref_;
            }

            //- Return cellSize
            scalar cellSize() const
            {
                return cellSize_;
            }

            //- Return yaw
            scalar yaw() const
            {
                return yaw_;
            }

            //- Return pitch
            scalar pitch() const
            {
                return pitch_;
            }

            //- Return omega
            scalar omega() const
            {
                return omega_;
            }

            //- Return TSR 
            scalar lambda() const
            {
                return lambda_;
            }

            //- Return rootFactor
            scalar rootFactor() const
            {
                return rootFactor_;
            }

            //- Return tipFactor
            scalar tipFactor() const
            {
                return tipFactor_;
            }

            //- Return rootDistance 
            scalar rootDistance() const
            {
                return rootDistance_;
            }

            //- table with UdAvg_table
            // List<Tuple2<vector, vector>> UdAvg_table() const
            List<Tuple2<vector, vector>> UdAvg_table() const //---new line
            {
                return UdAvg_table_;
            }

            //- table with Udi_table
            // List<Tuple2<vector, vector>> Udi_table() const
            List<Tuple2<vector, vector>> Udi_table() const //---new line
            {
                return Udi_table_;
            }

            //- Normal disk direction
            const vector &diskDir() const
            {
                return diskDir_;
            }

            //- Disk area
            scalar diskArea() const
            {
                return diskArea_;
            }

            //- method to measure Ud
            //  - 0 = sum directly values in cells in AD
            //  - 1 = weight with Gaussian in distance to AD plane
            //  - 2 = weight with Gaussian in distance to AD plane + distance to center of AD
            scalar UdCellsMethod() const
            {
                return UdCellsMethod_;
            }

            //- toggle for measuring Ud with center of AD
            //  - 0 = use all of AD
            //  - 1 = use center of AD with centerRatio
            scalar UdCenterToggle() const
            {
                return UdCenterToggle_;
            }

            scalar UdCorrection() const
            {
                return UdCorrection_;
            }

            scalar forceDistributionMethod() const
            {
                return forceDistributionMethod_;
            }

            scalar ADmodel() const
            {
                return ADmodel_;
            }

            //- ratio to consider in center of AD
            scalar centerRatio() const
            {
                return centerRatio_;
            }

            // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup(
                fvMatrix<vector> &eqn,
                const label fieldI);

            //- Source term to compressible momentum equation
            virtual void addSup(
                const volScalarField &rho,
                fvMatrix<vector> &eqn,
                const label fieldI);

            // IO

            //- Read dictionary
            virtual bool read(const dictionary &dict);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Foam::vector getNodePosition(
    Foam::scalar tita_n_Rad,
    Foam::scalar rMed_r,
    Foam::scalar yawRad,
    Foam::vector diskPoint_,
    int ring,
    int numberRings_);

Foam::vector getNodeVelocity(
    Foam::List<Foam::label> nodeCellID_,
    int total_nodes_counter,
    int ring,
    int numberRings_,
    Foam::vectorField U,
    Foam::volTensorField gradU,
    int nodesNumber_,
    int gradInterpolation_,
    Foam::vectorField cellCentres,
    Foam::vector Bi
);

Foam::vector getNodeFtitaDir(
    Foam::vector Bi,
    Foam::vector diskPoint,
    Foam::vector uniDiskDir
);

Foam::tensor getNodeTransformTensor(
    Foam::vector Bi,
    Foam::vector diskPoint,
    Foam::vector uniDiskDir
);

float getNodePhiAngle(
    Foam::tensor transform,
    Foam::vector U_dPointCells,
    float radius,
    float omega
);

float S0Function(
    float Ct,
    float Ct_rated
);

float posInTableUref2(
    int posI,
    Foam::List<Foam::scalar> Uref2List_,
    Foam::List<Foam::scalar> UrefList_,
    Foam::List<Foam::scalar> rList_,
    Foam::scalar rtable,
    Foam::vector U_dPointCells_ntr,
    Foam::List<Foam::scalar> UdiList_
);

float posInTableForcesVanDerLaan(
    int posUref,
    Foam::List<Foam::scalar> UrefList_,
    Foam::List<Foam::scalar> Uref2List_,
    Foam::List<Foam::scalar> rList_,
    Foam::scalar rtable
);

float sumForcesVanDerLaan(
    int force, // if 1 then sum fn, if 2 then sum ft
    int posUref,
    Foam::List<Foam::scalar> UrefList_,
    Foam::List<Foam::scalar> Uref2List_,
    Foam::List<Foam::scalar> rList_,
    Foam::List<Foam::scalar> fnList_,
    Foam::List<Foam::scalar> ftList_,
    Foam::scalar ringThickness,
    float maxR
);

float tipFactorFunction(
    int tipFactorType,
    float x,
    float lambda,
    // float phi
    float phi,
    Foam::vector UdPoint,
    float Uref
);

float rootFactorFunction(
    int rootFactorType,
    float x,
    float rootDistance,
    // float phi
    float phi,
    float lambda,
    Foam::vector UdPoint,
    float Uref
);

float gFunction(
    float x,
    float rootDistance
);

float FFunction(
    float x,
    float lambda
);

float a1Function(
    float rootDistance,
    float lambda
);

float a2Function(
    float rootDistance,
    float lambda
);

#ifdef NoRepository
#include "actuationDiskRingsV21_SourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
