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

#ifndef actuationDiskRingsV11_Source_H
#define actuationDiskRingsV11_Source_H

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
                     Class actuationDiskRingsV11_Source Declaration
\*---------------------------------------------------------------------------*/

class actuationDiskRingsV11_Source
:
    public cellSetOption
    //public option //---new line
{

protected:

    // Protected data

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

        //- Tip factor on/off
        scalar tipFactor_;
        
        //- Ratio Total amount of nodes / CellsInAD
        scalar nodesCellsRatio_;
        
        //Ratio ring Thickness / cellSize   
        scalar rThicknessCellsizeRatio_;
        
        //Interpolation by means of Grad(U).
        scalar gradInterpolation_; 

	    //- disk cell ID
        label diskCellId_;

	    //--test--------------------------------------
	    //the permanent variables have "_" in the end

        //for calculating number of rings and nodes
        
        //Root distance, % of R
        scalar rootDistance_;

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

    //- UdAvg_table
        //List<Tuple2<vector, vector>> UdAvg_table_;
        List<Tuple2<vector, vector> > UdAvg_table_; //--- new line

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
        //List<Tuple2<vector, vector>> Udi_table_;
    	List<Tuple2<vector, vector> > Udi_table_; //--- new line

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

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addactuationDiskRings_AxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U,
            vectorField& force
        ) const;


        //- Disallow default bitwise copy construct
        actuationDiskRingsV11_Source(const actuationDiskRingsV11_Source&);

        //- Disallow default bitwise assignment
        void operator=(const actuationDiskRingsV11_Source&);
		
	//- turbines output file 
	std::ofstream * outTurbines;

	//- blades output file 
	std::ofstream * outRings;
	
	//- New Calibration table output
	std::ofstream * outListFvOptions;


public:

    //- Runtime type information
    TypeName("actuationDiskRingsV11_Source");


    // Constructors

        //- Construct from components
        actuationDiskRingsV11_Source
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    ~actuationDiskRingsV11_Source(); 


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


            //- table with UdAvg_table
            //List<Tuple2<vector, vector>> UdAvg_table() const
            List<Tuple2<vector, vector> > UdAvg_table() const //---new line
            {
                return UdAvg_table_;
            }

	    //- table with Udi_table
            //List<Tuple2<vector, vector>> Udi_table() const
            List<Tuple2<vector, vector> > Udi_table() const //---new line
            {
                return Udi_table_;
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
#   include "actuationDiskRingsV11_SourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
