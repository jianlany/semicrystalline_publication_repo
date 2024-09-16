#pragma once
#include <cmath>
#include <tuple>
#include "fvector.h"

class InputFile;

//! Stores and computes geometrical paramaters for the semicrystalline domain.
/*! 
  Layout of crystalline/amorphous layers
  i=nz-1             /////// last crystalline layer
  ...                ///////  
  i=hi_amorphous_z+1 /////// hi boundary/root atoms
  i=hi_amorphous_z   &&&&&&& 
  ...                &&&&&&&
  i=lo_amorphous_z   &&&&&&& 
  i=lo_amorphous_z-1 /////// lo boundary/root atoms
  ...                ///////
  i=0                /////// first crystalline later
 */
class Geometry {
public:
    //! Initializes geometry parameters from command line options.
    Geometry(const InputFile &opt);
    //! Returns the bounds of the amorphous region.
    std::pair<double,double> amorphous_region() const;

    //! Latice constants for crystalline PE in unrotated configuration.
    const double a, b, c;
    //! Number of (rotated) unit cells along each direction.
    const int nx, ny, nz;
    // z-indices of atoms in the amorphous domain, next to the crystalline boundary (roots).
    const int lo_amorphous_z, hi_amorphous_z;
    //! Size of unit cells in rotated configuration.
    double step_x, step_y, step_z;

    //! Returns z-coord of crystalline atoms that are bonded to the top of 
    //! the amorphous domain.
    double lower_root_coordinate() const { return (lo_amorphous_z-1)*step_z; }
    //! Returns z-coord of crystalline atoms that are bonded to the bottom of 
    //! the amorphous domain.
    double upper_root_coordinate() const { return (hi_amorphous_z+1)*step_z; }

    //! Unique vector between successive beads on chains within crystal region.
    FVector<3> bv;
    //! The shifted distance of 0.5*bv projected on the global x direction.
    double s;
};

