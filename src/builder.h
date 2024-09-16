#pragma once
#include <string>
#include <map>
#include <cmath>
#include "fvector.h"
#include <random>

class MDSystem;
class InputFile;
class Geometry;

/*
Structure of domain:
  --     ...      --  |-> crystalline
  --     ...      --  |-> crystalline                       |
  --    roots     --  |-> crystalline    hi_amorphous_z+1   | possible buffer
  --     ...      --  |-> amorphous      hi_amorphous_z
  --     ...      --  |-> amorphous
  --     ...      --  |-> amorphous      lo_amorphous_z
  --    roots     --  |-> crystalline    lo_amorphous_z-1   | possible buffer
  --     ...      --  |-> crystalline                       |
  --     ...      --  |-> crystalline
*/

class Builder {
public:
    //! Initializes the MD builder.
    Builder(const Geometry &g, MDSystem &s, const InputFile &opt);
    //! Generate an initial system.
    void generate();
    //! Return chain-segment-id together with the atoms on it which indicate
    //  to the part should be cut.
    std::map<int,std::pair<int,int>> determine_cut_in_segments();
    //! Writes a list of root atoms to a specified file.
    void write_root_atoms(std::string path) const;
    //! Writes a list of amourphous atoms to a specified file.
    void write_amorphous_atoms(std::string path) const;
    //! Writes a list of movable atoms to a specified file.
    void write_movable_atoms(std::string path) const;

    //! Reference to domain geometry parameters.
    const Geometry &geometry;
private:
    //! Amount to shear amorphous bond vectors so that the domain is periodic.
    double _compute_bond_shift() const;

    //! Reference to the MD system this builder modifies.
    MDSystem &_system;
    //! Stores the atom types along a chain.
    std::string _block;
    //! Percentage of the chain segments in the amorphous region which
    //  will have atoms removed.
    double _initial_chain_cut_ratio;
    //! Set of seed points in global coordinates from which chains are grown.
    std::vector<FVector<3>> _seed_pts;
    //! Ratio of removing atoms in the amorphous phase.
    double _removal_ratio;
    //! Number of crystall atom layers that are allowed to move in lammps.
    int _crystal_buffer_depth;
};
