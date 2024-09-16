#pragma once
#include <string>
#include <map>
#include <cmath>
#include "md_system.h"
#include "geometry.h"

class InputFile;
class Builder;

/*! Implements Monte Carlo operations for finding an equilibrium configuration
 *! of the amorphous region within a semi-crystalline polyethylene domain. */
class MonteCarlo {
public:
    //! Sets a reference to the MDSystem.
    MonteCarlo(const Geometry &g, MDSystem &s, const InputFile &opt);
    //! Return all the atom ids in the amorphous region
    AtomSet tail_atoms(const AtomSet &ends) const;
    //! Perform end-bridging Monte Carlo moves.
    void end_bridge();
    //! Returns atom ids where current segment enters the crystalline phase.
    // 1 point -> tail; 2 points -> bridge or loop.
    AtomSet segment_entry_points(int id) const; 
    //! Return all the atom ids on the current segment.
    AtomSet segment_points(int id) const;
    //! Return all the atoms on the same chain as the atom id.
    std::vector<int> determin_chain_atoms(int id) const;
    //! Compute the acceptance probablity of the virtual loop/bridge formed.
    double loop_bridge_acceptance_probability(int A, int T1, int T2) const;
    //! Returns the number of bonds between atoms id1 and id2.
    AtomSet shortest_path(int id1, int id2) const;
    //! Returns the number of bonds between two beads.
    int bonds_between_atoms(int id1, int id2) const;
    //! Returns all the roots as the boundary of crystalline-amorphous phase.
    AtomSet determine_roots() const;
    //! Returns the numbers of tails, bridges and loops in the amorphous region.
    std::vector<AtomSet> count_segments(bool report, std::string out) const;
    //! Returns the shift in x and y direction of the two frozen atoms
    //! connected to the two roots of a loop segment.
    FVector<2> loop_entry_shift(int id, AtomSet roots) const;
    //  Return the id of a neighbor of aim, the bond {aim, neighbor}
    //  is to be removed, if there is enough atoms.
    //  Otherwise return -1.
    int select_target_bond(int A, int T) const;
    //! Returns true if atom is an amorphous atom.
    bool is_amorphous(int i) const;
    //! Returns true if atom is an upper root atom.
    bool is_lower_root(int i) const;
    //! Returns true if atom is a lower root atom.
    bool is_upper_root(int i) const;
    //! Returns true if atom is either lower or upper root atom.
    bool is_root(int i) const { return is_lower_root(i) || is_upper_root(i); }
private:
    //! Reference to domain geometry parameters.
    const Geometry &_geometry;
    //! This class will manipulate the bonds within the MDSystem.
    MDSystem &_system;
    //! Distance used to search for available bridges starting from a tail-end.
    double _search_radius;
    //! Factor used to check if a chain segment is over-stretched.
    double _length_factor;
    //! Number of crystall atom layers that are allowed to move in lammps.
    int _crystal_buffer_depth;
    //! Sorted list of upper root atom ids.
    AtomSet _upper_root_atoms;
    //! Sorted list of lower root atom ids.
    AtomSet _lower_root_atoms;
    //! Sorted list of amorphous atom ids.
    AtomSet _amorphous_atoms;
};

