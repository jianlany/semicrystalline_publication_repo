#include "builder.h"
#include "md_system.h"
#include "string_tools.h"
#include "input_file.h"
#include "geometry.h"
#include <set>

const double pi = acos(-1.0);

// Initializes the MD builder.
Builder::Builder(const Geometry &g, MDSystem &s, const InputFile &opt) 
  : geometry(g), _system(s) {
    // Sets simple constant values.
    _block = opt["block"];
    _initial_chain_cut_ratio = opt.get<double>("initial_chain_cut_ratio");
    _removal_ratio = opt.get<double>("removal_ratio");
    _crystal_buffer_depth = opt.get<int>("crystal_buffer_depth");
    // The global box has three layers: crystalline + amorphous + crystalline
    // The crystalline layers have equal thickness along z direction.
    _system.set_bounds(0, geometry.nx*geometry.step_x, 
                       0, geometry.ny*geometry.step_y, 
                       0, geometry.nz*geometry.step_z);

    // Template points.
    auto pa = FVector<3>{0.0, 0.0, 0.0};
    auto pb = FVector<3>{0.5*geometry.step_x, 0.5*geometry.step_y, 0.0}; // base-centered unit cell.
    // pb += 0.5 * geometry.bv;  // body-centered unit cell.
    // Generate seed points in the bottom plane of the global box.
    for (int i=0; i<geometry.nx; ++i) {
        for (int j=0; j<geometry.ny; ++j) {
            auto shift = FVector<3>{i*geometry.step_x, j*geometry.step_y, 0.0};
            _seed_pts.push_back(pa + shift);
            _seed_pts.push_back(pb + shift);
        }
    }
}

// Build a crystalline region which is in title angle theta with respect
// to the normal of the global basal plane (lamellar plane).
// The global coordinate system is:
// x - horizontal, point to right
//     the angle between lattice_a and x is the title ange theta
// y - point to inward, global rotation of the chains is around this axis
//     lattice_b in y direction
// z - vertical, point up
//     the angle between lattice_c and z is the title angle theta
void Builder::generate() {
    // Grow chains from seed points.
    // The growth along z is always upward (positive).
    // The growth along x is leftward (negative), but can have positive shift
    // due to periodicity, which corresponds to the chain's leaving from 
    // and entering into the global box.
    auto Lx = _system.box_length(0);
    auto Lz = _system.box_length(2);
    const auto segment_info = determine_cut_in_segments();
    std::set<int> occupied_seed_pts;

    for (int seed_id=0; seed_id<_seed_pts.size(); ++seed_id) {
        // If the seed point already has chain generated through it, continue.
        if (occupied_seed_pts.count(seed_id)) continue;
        // The new chain starts from the current seed point.
        auto x = _seed_pts[seed_id];
        _system.add_atom(x, _block[0]);
        occupied_seed_pts.insert(seed_id);

        // Grow a single chain with connectivity and periodicity.
        bool chain_complete = false;
        auto next_id = seed_id;
        while (!chain_complete) {
            // Vertically move through the global box from bottom to top, once.
            for (int i=0; i<geometry.nz; ++i) {
                x += geometry.bv;
                // In the amorphous region, the bond vector is always sheared.
                if (i>geometry.lo_amorphous_z && i<= geometry.hi_amorphous_z) {
                    x(0) += _compute_bond_shift();
                }

                // Check periodicity in x direction.
                if (x(0) < 0.0) x(0) += Lx;
                if (x(0) > Lx-1e-12) x(0) -= Lx;

                if (i < geometry.nz-1) {
                    auto iter = segment_info.find(next_id);
                    // Any atoms on the chain segments to be modified with
                    // its nz satisfying lo <= nz <= hi will be removed.
                    if (iter != segment_info.end()) {
                        auto lo = iter->second.first, hi = iter->second.second;
                        // The current x should be added to system.
                        if (i < lo || i > hi+1) {
                            _system.add_atom(x, _block[0]);
                            _system.add_bond(_system.size()-2, _system.size()-1);
                        }
                        // The beginning of a new chain, only add atom, not bond.
                        if (i == hi+1) {
                            _system.add_atom(x, _block[0]);
                        }
                    }
                    // Intact chain segment.
                    else {
                        _system.add_atom(x, _block[0]);
                        _system.add_bond(_system.size()-2, _system.size()-1);
                    }
                }
                // For the atom grown in the top plane (z==_nz*_step_z),
                // map it back to the bottom plane (z==0.0),
                // then check the mapped point overlaps with which seed point.
                else {
                    x(2) -= Lz;
                    auto iter = std::find_if(_seed_pts.begin(), _seed_pts.end(),
                                            [&x](const FVector<3> &y) { 
                                            return norm(x-y) < 1e-3; });
                    if (iter == _seed_pts.end()) {
                        std::cerr << "Wrapped around chain matches no seed points.\n";
                        exit(1);
                    }
                    next_id = std::distance(_seed_pts.begin(), iter);
                    // The current chain grows back to it starting seed point,
                    // chain grow complete.
                    if (next_id == seed_id) {
                        // Find the global id of the seed atom represented by
                        // the approximate position x.
                        auto id = _system.find_nearest_atom(x);
                        _system.add_bond(_system.size()-1, id);
                        chain_complete = true;
                    }
                    // Continue to grow the chain.
                    else {
                        occupied_seed_pts.insert(next_id);
                        _system.add_atom(_seed_pts[next_id], _block[0]);
                        _system.add_bond(_system.size()-2, _system.size()-1);
                    }
                }
            }
        }
    }
    _system.assign_initial_velocity(_system.size());
    _system.assign_molecules_ids();
    std::cout << "   Generated " << _system.size() << " atoms.\n";
    std::cout << "   Removed " << 2*geometry.nx*geometry.ny*geometry.nz - _system.size() << " atoms.\n";
}

// Return chain-segment-id together with the atoms on it which indicate
// to the part should be cut.
std::map<int,std::pair<int,int>> Builder::determine_cut_in_segments() {
    auto num_segment = _seed_pts.size();
    auto num_atom_on_segment = geometry.hi_amorphous_z - geometry.lo_amorphous_z;
    // Remove atoms to reduce density from 0.968g/cc to 0.784g/cc.
    auto num_atom_to_remove = int(num_segment * num_atom_on_segment * _removal_ratio);
    // The parameter _initial_chain_cut_ratio should be
    // less than 1/3 but higher than 1/5.
    auto num_segment_to_cut = int(num_segment * _initial_chain_cut_ratio);
    // Average number of atoms to remove on each cut segment.
    auto mu = num_atom_to_remove / num_segment_to_cut;
    std::cout << "Approximate # of atoms removed from amorphous region: " 
              << num_atom_to_remove << "\n";
    std::cout << "# of chain segments entering amorphous region: " 
              << num_segment << "\n";
    std::cout << "# of atoms on each chain segments: " 
              << num_atom_on_segment << "\n";
    std::cout << "# of chain segments to be modified in amorphous region: " 
              << num_segment_to_cut << "\n";
    std::cout << "# of bridges in amorphous region: " 
              << num_segment - num_segment_to_cut << "\n";
    std::cout << "Average # of atoms to remove on each modified segment: " 
              << mu << "\n";
    // Generate set of random chain segments id to be modified.
    std::set<int> seed_id_to_cut;
    while (seed_id_to_cut.size() < num_segment_to_cut) {
        seed_id_to_cut.insert(rand() % _seed_pts.size());
    }

    // Maximum number of atoms that can be removed on a chain segment.
    auto maximum = num_atom_on_segment - 2;
    // Generate random integers summing up close to num_atom_to_remove.
    int  sum = 0;
    std::vector<int> segment_atom_to_remove;
    for (int i=0; i<num_segment_to_cut; ++i) {
        auto k = 1 + rand() % maximum;
        sum += k;
        segment_atom_to_remove.push_back(k);
    }
    auto r = (double) num_atom_to_remove/sum;
    sum = 0;
    for (auto &k: segment_atom_to_remove) {
        k = (int) round(k*r);
        if (k < 1) k = 1;
        if (k > maximum) k = maximum;
        sum += k;
    }
    auto diff = num_atom_to_remove - sum;
    if (diff > 0) {
        for (auto &k: segment_atom_to_remove) {
            if (k == maximum) continue;
            auto n = rand() % (maximum-k);
            if (diff-n < 0 ) continue;
            k    += n;
            diff -= n;
        }
    }
    else if (diff < 0) {
        for (auto &k: segment_atom_to_remove) {
            if (k == 1) continue;
            auto n = rand() % k;
            if (diff+n > 0 ) continue;
            k    -= n;
            diff += n;
        }
    }
    sum = 0;
    for (auto &k: segment_atom_to_remove) sum += k;

    // Return mapped information.
    int j = 0;
    std::map<int,std::pair<int,int>> segment_info;
    for (int i=0; i<_seed_pts.size(); ++i) {
        if (seed_id_to_cut.count(i)) {
            auto r = ((double) rand() / (RAND_MAX));
            // Make sure nz_remove_lo >= geometry.lo_amorphous_z +1
            // and       nz_remove_hi <= geometry.hi_amorphous_z -1
            auto nz_remove_lo = int(
                          (num_atom_on_segment - segment_atom_to_remove[j]) * r
                          ) + geometry.lo_amorphous_z + 1;
            auto nz_remove_hi = nz_remove_lo + segment_atom_to_remove[j] - 2;
            if (nz_remove_hi > geometry.hi_amorphous_z -1) {
                nz_remove_hi = geometry.hi_amorphous_z -1;
                }
            segment_info[i] = std::make_pair(nz_remove_lo, nz_remove_hi);
            j += 1;
        }
    }
    return segment_info;
}

// Writes a list of *sorted* root atoms to a specified file.
// First line is upper root atoms, second line is lower root atoms.
void Builder::write_root_atoms(std::string path) const {
    std::fstream f(path, std::ios::out);
    auto zhi = geometry.upper_root_coordinate();
    auto zlo = geometry.lower_root_coordinate();
    double half_layer = 0.5*geometry.step_z;
    const double SEARCH_TOL = 1e-3;
    auto lo = _system.atoms_in_slice(2, zlo-SEARCH_TOL, zlo+half_layer+SEARCH_TOL);
    auto hi = _system.atoms_in_slice(2, zhi-SEARCH_TOL, zhi+half_layer+SEARCH_TOL);
    std::sort(lo.begin(), lo.end());
    std::sort(hi.begin(), hi.end());    
    for (auto i: hi) f << i+1 << " ";
    f << "\n";
    for (auto i: lo) f << i+1 << " ";
}

// Writes a list of amorphous atoms to a specified file.
void Builder::write_amorphous_atoms(std::string path) const {
    auto zhi = geometry.upper_root_coordinate();
    auto zlo = geometry.lower_root_coordinate();
    double half_layer = 0.5*geometry.step_z;
    const double SEARCH_TOL = 1e-3;
    auto amor_atoms = _system.atoms_in_slice(2, zlo+half_layer+SEARCH_TOL, zhi-SEARCH_TOL);
    std::fstream f(path, std::ios::out);
    for (auto &i: amor_atoms) {
        f << i+1 << " " << "\n";
    }
}

// Writes a list of movable atoms to a specified file.
// amorphous_atoms + buffer_atoms (including root_atoms).
void Builder::write_movable_atoms(std::string path) const {
    auto zhi = geometry.upper_root_coordinate();
    auto zlo = geometry.lower_root_coordinate();
    double half_layer = 0.5*geometry.step_z;
    const double SEARCH_TOL = 1e-3;
    // Determine the lower and upper bound of buffer atoms.
    auto lower_roots = _system.atoms_in_slice(2, zlo-SEARCH_TOL, zlo+half_layer+SEARCH_TOL);
    auto upper_roots = _system.atoms_in_slice(2, zhi-SEARCH_TOL, zhi+half_layer+SEARCH_TOL);
    std::vector<int> lower_buffers, upper_buffers;
    for (auto &i: lower_roots) {
        lower_buffers.push_back(i-_crystal_buffer_depth);
    }
    for (auto &i: upper_roots) {
        upper_buffers.push_back(i+_crystal_buffer_depth);
    }
    // Determine the smallest and largest z coordinates of buffer atoms.
    double hi = -1e100, lo = 1e100;
    for (auto &i: lower_buffers) {
        auto z = _system[i](2);
        if (z < lo) lo = z;
    }
    for (auto &i: upper_buffers) {
        auto z = _system[i](2);
        if (z > hi) hi = z;
    }
    // Determine and output all movable atoms.
    auto movable_atoms = _system.atoms_in_slice(2, lo-SEARCH_TOL, hi+SEARCH_TOL);
    std::fstream f(path, std::ios::out);
    std::string buffer = "group movable_atoms id ";
    for (auto i: movable_atoms) {
        auto s = std::to_string(i+1);

        // Flush buffer when it gets to long.
        if (buffer.size() + s.size() + 3 >= 80) {
            f << buffer << " &\n";
            buffer = "";
        }
        buffer += " " + s;
    }
    f << buffer << "\n";
}

/* Note that _step_x*bv(0) is not an integer, so beads at the top and bottom
 * will not line up properly to maintain periodic symmetry.  As a solution, we
 * shear all bonds in the amorphous region by a small value in the x-direction 
 * so that the bonds between the top and bottom planes are undeformed.
 */
double Builder::_compute_bond_shift() const {
    auto x_last = geometry.nz*geometry.bv(0);
    auto i = round(x_last / geometry.step_x);
    auto dz = geometry.hi_amorphous_z-geometry.lo_amorphous_z;
    return (i*geometry.step_x - x_last) / double(dz);
}
