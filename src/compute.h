#pragma once
#include <string>
#include <ctype.h>
#include "fvector.h"

class MDSystem;
class Geometry;
class InputFile;

//! Stores slice coordinate (first) and data (second).
typedef std::pair<std::vector<double>, std::vector<double>> Slice;
typedef std::vector<int> AtomSet;

//! Post processing of averaged structural and energy metrics for analysis.
/*!
    Implements routines for computing density, bond orientation parameter,
    pair, bond, and angle energies with window averaging along the z-direction.
 */
class Compute {
public:
    //! Sets a reference to the MDSystem.
    Compute(const Geometry &g, MDSystem &s, const InputFile &opt);
    //! Perform all computations.
    void equilibrium_check(std::string pairpath,
                           std::string anglepath,
                           std::string out,
                           std::vector<AtomSet> loops,
                           std::string chain_file,
                           std::string end_file);
    //! Perform only end2end computations.
    void perform_all_e2e(std::string out);
    //! Returns the density distributionalong z direction.
    std::vector<FVector<3>> slice_density(int n_slices) const;
    //! Returns the bond orientation distribution along z direction.
    Slice slice_bond_orientation(int n_slices) const;
    //! Returns the bond orientation distribution along z direction
    //  contributed only by loops.
    Slice slice_loop_bond_orientation(int n_slices,
                                      std::vector<AtomSet> loops) const;
    //! Computes the pair-energy density distribution along z direction.
    Slice slice_pair_energy(int n_slices, std::string pairpath) const;
    //! Computes the bond-energy density distribution along z direction.
    Slice slice_bond_energy(int n_slices) const;
    //! Computes the angle-energy density distribution along z direction.
    Slice slice_angle_energy(int n_slices, std::string anglepath) const;
    //! Compute the bond angle between atoms i1-i2-i3.
    double bond_angle(int i, int j, int k) const;
    //! Determine the chain ends of individual chains.
    std::vector<std::vector<int>> determin_chain_ends() const;
    //! Compute the average end2end distance squared <R.R> for system.
    double e2e_squared() const;
    //! Write individual chains in format "{chain_id} chain_type (ends)\n beads".
    void write_chains(std::string chain_file);
    //! Write end-atom coordinates.
    void write_end_atoms(std::string end_file);

private:
    //! Reference to domain geometry parameters.
    const Geometry &_geometry;
    //! This class will manipulate the bonds within the MDSystem.
    MDSystem &_system;
    //! CG system parameters.
    double _bead_mass;
    double _cutoff_distance;
    double _bond_stiffness, _bond_length;
    //! Thickness of each slice.
    double _dz;
    double _st;
    //! MC step index.
    int _step;
};

//! Read movable atoms.
AtomSet read_movable_atom(std::string movable_atom_path);
