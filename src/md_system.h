#pragma once
#include "fvector.h"

//! Stores per timestep data of the MD system.
struct StepData {
    //! Timestep number.
    int timestep = 0;
    //! Coordinates of all atoms.
    std::vector<FVector<3>> atom_coords; 
    //! Velocities of all atoms.
    std::vector<FVector<3>> atom_velcs; 
    //! The belonging molecule of each atom.
    std::vector<int> atom_mols;
    //! The type of each atom.
    std::vector<int> atom_types;
    //! The (averaged) von Mises stress of each atom.
    std::vector<double> atom_vms;
    //! The (averaged) centro-symmetry parameter of each atom.
    std::vector<double> atom_css;
    //! Stores min/max extent of simulation box in each dimension.
    FVector<2> box[3];
    //! Size of the simulation box along each dimension.
    FVector<3> DX;
};
// Shortcut defined for sets of atoms.
typedef std::vector<int> AtomSet;

//! Class that stores all MD simulation data.
/*! This class stores atom coords, types and bond information
 *  All indices are stored zero-based.
 */
class MDSystem {
    //! Writes the current system to a LAMMPS data file.
    friend void write_to_lammps(const MDSystem&, std::string);
    //! Read data can manipulate private data members.
    friend void read_data(std::string datapath, MDSystem &sys);
public:
    //! Adds an atom to the domain.
    void add_atom(const FVector<3> &x, int type);
    //! Update all atoms' belonging molecule from external source.
    void update_atom_mol(std::vector<int> mols);
    //! Update all atoms' type from external source.
    void update_atom_type(std::vector<int> types);
    //! Update all atoms' von Mises stress from external source.
    void update_atom_vm(std::vector<double> vms);
    //! Update all atoms' centro-symmetry parameter from external source.
    void update_atom_cs(std::vector<double> css);
    //! Update all the atom coordinates from external source.
    void update_atom_coordinate(std::vector<FVector<3>> coords);
    //! Return the molecule id of an atom.
    int get_atom_molecule(int atom_id) const;
    //! Adds a bond between atoms.
    void add_bond(int i, int j);
    //! Deletes a bond between i and j.
    void delete_bond(int i, int j);
    //! Returns the number of atoms in the system.
    size_t size() const { return _step.atom_coords.size(); }
    //! Returns the atoms at the ends of each chain.
    AtomSet end_atoms() const;
    //! Returns the index of the atom nearest to point x.
    size_t find_nearest_atom(const FVector<3> &x) const;
    //! Returns all atoms within a linear region.
    AtomSet atoms_in_slice(int dim, double lo, double hi) const;
    //! Filter atoms by the distance from the tail-end-atom position x.
    AtomSet find_atoms_within_cutoff(int id, double radius) const;
    //! Provide a shortcut way of accessing the coordinate of an atom.
    const FVector<3>& operator[](int i) const { return _step.atom_coords[i]; }
    //! Returns the belonging molecule of an atom.
    const int& mol(int id) const { return _step.atom_mols[id]; }
    //! Returns the type of an atom.
    const int& type(int id) const { return _step.atom_types[id]; }
    //! Returns the von Mises stress of an atom.
    const double& vm(int id) const { return _step.atom_vms[id]; }
    //! Returns the centro-symmetry parameter of an atom.
    const double& cs(int id) const { return _step.atom_css[id]; }
    //! Returns the Velocity of an atom.
    const FVector<3>& velocity(int id) const { return _step.atom_velcs[id]; }
    //! Assign zero velocities for every atom.
    void assign_initial_velocity(int natoms) { _step.atom_velcs.assign(natoms, zeros<3>()); };
    //! Returns all atoms connected to an atom by a bond.
    const AtomSet& bonded_neighbors(int id) const { return _bond_table[id]; }
    //! Returns the bounds along a single direcion.
    const FVector<2>& bounds(int i) const { return _step.box[i]; } 
    //! Returns all bounds information packed in a vector.
    std::vector<FVector<2>> bounds() const { 
        return {_step.box[0], _step.box[1], _step.box[2]}; 
    } 
    void set_bounds(double xlo, double xhi, 
                    double ylo, double yhi, 
                    double zlo, double zhi) {
        _step.box[0] = {xlo, xhi};
        _step.box[1] = {ylo, yhi};
        _step.box[2] = {zlo, zhi};
    }
    //! Returns the box length along a direction.
    const double box_length(int i) const {
        return _step.box[i](1) - _step.box[i](0);
    }
    //! Assign molecule id of all stoms in system.
    int assign_molecules_ids();
    //! Returns true if atoms i and j share a bond.
    bool are_bonded(int i, int j) const;
    //! Returns the minimum image distance between two atoms.
    double distance(int i, int j) const { return norm(bond_image(i, j)); }
    //! Returns the shortest image of a bond between atoms i and j.
    const FVector<3> bond_image(int i, int j) const {
        return bond_image((*this)[j] - (*this)[i]);
    }
    //! Computes the image distance vector given a periodic domain.
    const FVector<3> bond_image(FVector<3> dx) const {
        for (int i=0; i<3; ++i) {
            auto DX = box_length(i);
            while ( fabs(dx(i)-DX) < fabs(dx(i)) ) dx(i) -= DX;
            while ( fabs(dx(i)+DX) < fabs(dx(i)) ) dx(i) += DX;
        }
        return dx;
    }
    //! Wrap atom coordiante to be within the system box.
    void warp_atoms();

private:
    //! Stores type of each atom.
    std::vector<int> _atom_types;
    //! Stores which molecule (chain) each atom is located in.
    std::vector<int> _atom_molecule;
    //! Bond type is based on atom types.  Atom types are always stored
    //! symmetrically, i.e. {i,j} where i<j, so that the bond type does 
    //! not depend on the ordering of the atoms.
    std::vector<std::pair<int,int>> _bond_types;
    //! Connectivity of all atoms.
    std::vector<std::vector<int>> _bond_table;
    //! Stores all time-dependent data (atom coords, etc).
    StepData _step;
};

//! Reads coordinates and bonds from the specified lammps data file.
void read_data(std::string datapath, MDSystem &sys);

