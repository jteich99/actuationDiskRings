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

#ifndef actuationDiskRingsV11_calib_Source_H
#define actuationDiskRingsV11_calib_Source_H

#include "cellSetOption.H"
//#include "fvOption.H" //---new line
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
                     Class actuationDiskRingsV11_calib_Source Declaration
\*---------------------------------------------------------------------------*/

class actuationDiskRingsV11_calib_Source
:
    public cellSetOption
    //public option //---new line
{
    protected: // Protected data

        //- Disk area normal
        vector diskDir_;

        //- Power coefficient
        scalar Cp_;

        //- Thrust coefficient
        scalar Ct_;

        //- Fixed Uinf for thrust calculation
        scalar Uref_;

        //- cell size in the disc
        scalar cellSize_;

        //- Fixed yaw angle compared to the inlet velocity vector
        scalar yaw_;

        //- Fixed rotational velocity [rad/s]
        scalar omega_;

        //- Fixed pitch angle [deg]
        scalar pitch_;

        //- Disk area
        scalar diskArea_;

        //- Disc center point sample
        point diskPoint_;

        //- Root factor on/off
        scalar rootFactor_;
        
        //Root distance, % of R
        scalar rootDistance_;

        //- Tip factor on/off
        scalar tipFactor_;
        
        //- Ratio Total amount of nodes / CellsInAD
        scalar nodesCellsRatio_;
        
        //Ratio ring Thickness / cellSize   
        scalar rThicknessCellsizeRatio_;

        //Fluid density   
        scalar density_;
        
        //Interpolation by means of Grad(U).
        scalar gradInterpolation_; 

        //- disk cell ID
        label diskCellId_;

        //--test--------------------------------------
        //the permanent variables have "_" in the end

        //for calculating number of rings and nodes
        

        //AD radius
        scalar maxR_;

        //Cells in the AD
        scalar cellsInAd_;

        //Estimated nodes
        scalar estimatedNodes_;

        //Aproximate are associated to each node
        scalar secArea_;

        //Internal radius of the AD, inside there is only 1 node.
        scalar rInt_;        

        //External radius of the AD, equal to maxR_ 
        scalar rExt_;

        //Thickness of all rings
        scalar ringThickness_;

        //Number of rings
        int numberRings_;

        //Total nodes in the AD
        scalar nodesNumber_;
        
        //Nodes in each ring
        scalar nodesI_;

        //Angle between nodes - deg
        scalar titaI_;

        //Radius of each ring
        scalar rMedI_;

        //sectional area of each ring
        scalar areaI_;

        //Lists of nodes, angle between nodes and radius of each ring.
        DynamicList<scalar> ringNodesList_;
        DynamicList<scalar> ringTitaList_;
        DynamicList<scalar> ringrMedList_;
        DynamicList<scalar> ringAreaList_;

        //direction of the disc
        scalar yawRad_;
        vector diskYawed_;
        vector uniDiskDir_;
        
        //For calculating node positions
        //angle between nodes in a certain ring         
        scalar tita_r_;
        
        //angle for the position of a certain node
        scalar tita_n_Rad_;

        //radius for a certain ring
        scalar rMed_r_;

        //Nodes counter
        scalar total_nodes_counter_;

        //List of cell ID for each node
        List<label> nodeCellID_;

    //--------------------------------------------

private:

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        
        template<class RhoFieldType>
        void addactuationDiskRingsV11_calib_AxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U,
            vectorField& force
        ) const;


        //- Disallow default bitwise copy construct
        actuationDiskRingsV11_calib_Source(const actuationDiskRingsV11_calib_Source&);

        //- Disallow default bitwise assignment
        void operator=(const actuationDiskRingsV11_calib_Source&);
		
	//- turbines output file 
	std::ofstream * outTurbines;

	//- blades output file 
	std::ofstream * outRings;


public:

    //- Runtime type information
    TypeName("actuationDiskRingsV11_calib_Source");


    // Constructors

        //- Construct from components
        actuationDiskRingsV11_calib_Source
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    ~actuationDiskRingsV11_calib_Source(); 


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

			//- Return ratio Nodes/cellsInAD
            scalar nodesCellsRatio() const
            {
                return nodesCellsRatio_;
            }
            
            //- Return ratio ringThickness/cellSize
            scalar rThicknessCellsizeRatio() const
            {
                return rThicknessCellsizeRatio_;
            }
            
            //- Return density
            scalar density() const
            {
                return density_;
            }

            //- Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

            //- Disk area
            scalar diskArea() const
            {
                return diskArea_;
            }

         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "actuationDiskRingsV11_calib_SourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
