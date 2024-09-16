#include "md_system.h"
#include "string_tools.h"
#include <boost/format.hpp>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <tuple>

void MDSystem::add_atom(const FVector<3> &x, int type) {
    _step.atom_coords.push_back(x);
    _atom_types.push_back(type);
}

int MDSystem::get_atom_molecule(int atom_id) const {
    return _atom_molecule[atom_id];
}

// Update all atoms' belonging molecule from external source.
void MDSystem::update_atom_mol(std::vector<int> mols) {
    _step.atom_mols = mols;
}

// Update all atoms' type from external source.
void MDSystem::update_atom_type(std::vector<int> types) {
    _step.atom_types = types;
}

// Update all atoms' von Mises stress from external source.
void MDSystem::update_atom_vm(std::vector<double> vms) {
    _step.atom_vms = vms;
}

//! Update all atoms' centro-symmetry parameter from external source.
void MDSystem::update_atom_cs(std::vector<double> css) {
    _step.atom_css = css;
}

// Update all the atom coordinates from external source.
void MDSystem::update_atom_coordinate(std::vector<FVector<3>> coords) {
    _step.atom_coords = coords;
}

// Wrap atom coordiante to be within the system box.
void MDSystem::warp_atoms() {
    auto box = bounds();
    for (auto &x: _step.atom_coords) {
        for (int i=0; i<3; ++i) {
            while (x(i)<box[i](0)) x(i) += box_length(i);
            while (x(i)>box[i](1)) x(i) -= box_length(i);
        }
    }
}

void MDSystem::add_bond(int i, int j) {
    // Make sure that there are enough rows for max(i,j).
    while (std::max(i,j) >= _bond_table.size()) {
        _bond_table.emplace_back();
    }

    _bond_table[i].push_back(j);
    _bond_table[j].push_back(i);

    // Note that this works only if we already have _atom_types.
    // This requires to read coordinates info first, then read bond info.
    int itype = _atom_types[i];
    int jtype = _atom_types[j];

    if (itype > jtype) std::swap(itype, jtype);
    auto bond_type = std::make_pair(itype,jtype);

    auto iter = std::find(_bond_types.begin(), _bond_types.end(), bond_type);

    if (iter == _bond_types.end()) {
        _bond_types.push_back(bond_type);
    }
}

// Deletes a bond between i and j.  Since the bond table is symmetric, we need 
// to delete two entries.
void MDSystem::delete_bond(int i, int j) {
    auto iter = std::find(_bond_table[i].begin(), _bond_table[i].end(), j);
    if (iter == _bond_table[i].end()) {
        std::cerr << "Cannot remove bond (" << i << ", " << j << ")\n";
        exit(1);
    }
    _bond_table[i].erase(iter);
    iter = std::find(_bond_table[j].begin(), _bond_table[j].end(), i);
    if (iter == _bond_table[j].end()) {
        std::cerr << "Cannot remove bond (" << j << ", " << i << ")\n";
        exit(1);
    }
    _bond_table[j].erase(iter);
}

// Returns the atom ids for the ends of each chain.
AtomSet MDSystem::end_atoms() const {
    // Finds all atoms in the chain that are connected by only one bond.
    AtomSet chain_ends;
    for (int i=0; i<_bond_table.size(); ++i) {
        if (_bond_table[i].size() == 1) {
            chain_ends.push_back(i); 
        }
    }
    std::sort(chain_ends.begin(), chain_ends.end());
    return chain_ends;
}

// Returns the index of the atom nearest to point x.
// Not optimized - should use a bin table for this.
size_t MDSystem::find_nearest_atom(const FVector<3> &x) const {
    auto min_distance = 0.0;
    int nearest = -1;
    for (int i=0; i<size(); ++i) {
        auto distance = norm(bond_image((*this)[i] - x));
        if (i==0 || distance < min_distance) {
            nearest = i;
            min_distance = distance;
        }
    }
    return nearest;
}

// Returns all atoms within a linear region
AtomSet MDSystem::atoms_in_slice(int dim, double lo, double hi) const {
    AtomSet atoms;
    for (int i=0; i<size(); ++i) {
        auto x = (*this)[i](dim);
        if (x > lo && x < hi) {
            atoms.push_back(i);
        }
    }
    return atoms;
}

// Returns all atoms within a radius of a specific atom.
// Not optimized - should use a bin table for this.
AtomSet MDSystem::find_atoms_within_cutoff(int id, double radius) const {
    AtomSet within_radius;
    for (int i=0; i<size(); ++i) {
        if (i != id && distance(i, id) < radius) {
            within_radius.push_back(i);
        }
    }
    return within_radius;
}

// Assign molecule id of all atoms in system.
int MDSystem::assign_molecules_ids() {
    // All atoms start as unassigned to a molecule (-1).
    _atom_molecule.assign(size(), -1);

    int current_molecule = 0;
    for (int i=0; i<size(); ++i) {
        if (_atom_molecule[i] != -1) continue;
        _atom_molecule[i] = current_molecule;

        auto queue = bonded_neighbors(i);
        while (!queue.empty()) {
            auto next = queue.back(); queue.pop_back();
            _atom_molecule[next] = current_molecule;
            for (auto j: bonded_neighbors(next)) {
                if (_atom_molecule[j] == -1) {
                    queue.push_back(j);
                }
            }
        }
        // Found no more atoms successivly bonded to i, starting a new
        // molecule from the next unassigned molecule.
        current_molecule += 1;
    }
    std::cout << "Assigned " << current_molecule << " unique molecules.\n";
    return current_molecule;
}

// Returns true if atoms i and j share a bond.
bool MDSystem::are_bonded(int i, int j) const {
    return std::count(_bond_table[i].begin(), _bond_table[i].end(), j) > 0;
}

//! Returns the index of a member of a container.
template<class Container, typename T> 
inline size_t index_of(const Container &c, const T &t) {

    auto iter =  std::find(std::begin(c), std::end(c), 
                           decltype(*std::begin(c))(t));
    if (iter == c.end()) return -1;
    return std::distance(std::begin(c), iter);
}

// Reads coordinates and bonds from the specified lammps data file.
void read_data(std::string datapath, MDSystem &sys) {
    std::fstream fid(datapath, std::ios::in);
    if (!fid) {
        std::cerr << "Can't open LAMMPS data file " << datapath << "\n";
        exit(1);
    }
    auto line = read_line(fid);  // header
    line = read_line(fid);  // blank
    line = read_line(fid);  // # of atoms.
    int natoms = from_string<int>(split(line)[0]);
    
    FVector<2> box[3];
    while (fid) {
        line = read_line(fid);
        auto data = split(line);
        if (data.size() == 4 && data[3] == "xhi") {
            box[0](0) = from_string<double>(data[0]);
            box[0](1) = from_string<double>(data[1]);
            break;
        }
    }
    // Set box[1-2].  
    for (int i=1; i<=2; ++i) {
        auto data = split(read_line(fid));
        box[i](0) = from_string<double>(data[0]);
        box[i](1) = from_string<double>(data[1]);
    }   
    sys.set_bounds(box[0](0), box[0](1),
                   box[1](0), box[1](1),
                   box[2](0), box[2](1));

    while (fid) {
        if (startswith(read_line(fid), "Atoms")) break;
    }
    // Read atom coordinates and types.
    sys._atom_types.assign(natoms, -1);
    sys._step.atom_coords.assign(natoms, zeros<3>());
    while (fid) {
        line = read_line(fid);
        if (startswith(line, "Velocities")) break;
        auto data = split(line);
        if (data.size() != 9) continue;
        // Store atom-id as zero-indexed.
        int id = from_string<int>(data[0]) - 1;
        sys._atom_types[id] = from_string<int>(data[2]);
        sys._step.atom_coords[id] = {from_string<double>(data[3]),
                                     from_string<double>(data[4]),
                                     from_string<double>(data[5])};
    }
    // Read atom Velocities.
    sys._step.atom_velcs.assign(natoms, zeros<3>());
    while (fid) {
        line = read_line(fid);
        if (startswith(line, "Bonds")) break;
        auto data = split(line);
        if (data.size() != 4) continue;
        // Store atom-id as zero-indexed.
        int id = from_string<int>(data[0]) - 1;
        sys._step.atom_velcs[id] = {from_string<double>(data[1]),
                                    from_string<double>(data[2]),
                                    from_string<double>(data[3])};
    }
    // Read bonds.
    while (fid) {
        line = read_line(fid);
        if (startswith(line, "Angles")) break;
        auto bond = split(line);
        if (bond.size() != 4) continue;
        // Store bonds as zero-indexed.
        int i = from_string<int>(bond[2])-1;
        int j = from_string<int>(bond[3])-1;
        sys.add_bond(i,j);
    }
}

