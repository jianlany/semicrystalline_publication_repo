#include "end_bridge.h"
#include "md_system.h"
#include "string_tools.h"
#include "input_file.h"
#include <set>
#include <list>
#include <unordered_set>
#include <algorithm>

const double pi = acos(-1.0);

MonteCarlo::MonteCarlo(const Geometry &g, MDSystem &s, const InputFile &opt) 
  : _geometry(g), _system(s) {
    _length_factor = opt.get<double>("length_factor");
    _crystal_buffer_depth = opt.get<int>("crystal_buffer_depth");

    // Reads root atoms.
    std::fstream fid1("root_atoms.txt", std::ios::in);
    if (!fid1) {
        std::cerr << "Can't open root atom file\n";
        exit(1);
    }
    for (auto &s: split(read_line(fid1))) {
        _upper_root_atoms.push_back(from_string<int>(s)-1);
    }
    for (auto &s: split(read_line(fid1))) {
        _lower_root_atoms.push_back(from_string<int>(s)-1);
    }
    // Reads amorphous.
    std::fstream fid2("amorphous_atoms.txt", std::ios::in);
     if (!fid2) {
        std::cerr << "Can't open amorphous atom file\n";
        exit(1);
    }
    while (fid2) {
        int i;  fid2 >> i;
        _amorphous_atoms.push_back(i-1);
    }

    if (opt["current_stage"]=="melt")
        _search_radius = from_string<double>(split(opt["search_radius"])[0]);
    else
        _search_radius = from_string<double>(split(opt["search_radius"])[1]);
}

// Return all the amorphous or root atom ids
// which are connected to any chain ends (tail ends).
AtomSet MonteCarlo::tail_atoms(const AtomSet &ends) const {
    std::set<int> atoms;
    for (auto i: ends) {
        atoms.insert(i);
        // If the end is a root atom, then there are no more atoms to search.
        if (is_root(i)) continue;
        auto queue = _system.bonded_neighbors(i);
        // Search along bond connectivity until we reach a root atom.
        while (!queue.empty()) {
            auto next = queue.back(); queue.pop_back();
            atoms.insert(next);
            if (is_root(next)) {
                continue;
            }
            else {
                for (auto j: _system.bonded_neighbors(next)) {
                    if (atoms.count(j) == 0) {
                        queue.push_back(j);
                    }
                }
            }
        }
    }
    return AtomSet(atoms.begin(), atoms.end());
}

// Removes all elements of A that are in B (B must be sorted).
inline void exclude_atoms(std::vector<int> &A, const std::vector<int> &B) {
    A.erase(std::remove_if(A.begin(), A.end(), 
               [&](int i) { return std::binary_search(B.begin(), B.end(), i);}),
            A.end());
}

// Perform Monte Carlo moves.
void MonteCarlo::end_bridge() {
    // Determine the ids of chain ends which must be in the amorphous region.
    auto chain_ends = _system.end_atoms();
    // All amorphous or root atoms connected to chain (tail) ends.
    auto excludes = tail_atoms(chain_ends);
    // Removes atoms from set x that are not amorphous.
    auto remove_crystalline_atoms = [this](AtomSet &x) {
        auto iter = std::remove_if(x.begin(), x.end(), [this](int i) {
                return !(this->is_amorphous(i));
            });
        x.erase(iter, x.end());
    };

    const int max_attempt = 10000;
    for (int attempt=0; attempt<max_attempt; ++attempt) {
        // Find a valid attacking tail-end-bead A.
        int A;
        for (int i=1; i<max_attempt; ++i) {
            A = chain_ends[ rand() % chain_ends.size() ];
            auto test1 = !is_root(A);
            auto test2 = segment_entry_points(A).size() == 1;
            if (test1 && test2) break;
            if (i == max_attempt -1)
                std::cout << "Failed to find a valid attacking tail end\n";
        }
        // Filter atoms by the distance from the tail-end-atom.
        auto end_neighbors = _system.find_atoms_within_cutoff(A, _search_radius);
        // Return only neighbors that are amorphous atoms.
        remove_crystalline_atoms(end_neighbors);
        // Filter atoms by whether or not they are on tail segments.
        exclude_atoms(end_neighbors, excludes);
        // Filter atoms which are on the same chain of atom A.
        // This prevents a tail end attacking the chain it belongs to,
        // which may form an infinite ring chain.
        auto chain_atoms_of_A = determin_chain_atoms(A);
        exclude_atoms(end_neighbors, chain_atoms_of_A);
        // Handle special case.
        if (end_neighbors.size() < 1) continue;

        // Randomly select target atom T1 which is not a root atom.
        // If T1 is a root atom, then the shortest-path between T1 and
        // the roots of T1 will be {T1} alone, and this is not handled
        // by select_target_bond(A, T1).
        int T1;
        bool pass=true;
        for (int i=1; i<max_attempt; ++i) {
            T1 = end_neighbors[rand() % end_neighbors.size()];
            // Make sure the bonded neighbors of T1 are not root atoms.
            // So zero-length tail segment (single root atom as tail)
            // will never be generated.
            for (int j: _system.bonded_neighbors(T1)) {
                auto t1 = is_lower_root(j);
                auto t2 = is_upper_root(j);
                if (t1 || t2) {
                    pass = false;
                    break;
                }
            }
            if (!pass) break;
            auto test1 = !is_root(T1);
            auto test2 = segment_entry_points(T1).size() == 2;
            if (test1 && test2) break;
        }
        if (!pass) continue;
        // Select the atom T2 previously connected to T1.
        auto T2 = select_target_bond(A, T1);

        if (T2 >= 0) {
            std::cout << "Performed " << attempt+1 << " attempts to "
                << "find attacker and target:\n"
                << "  search radius is " << _search_radius << " A\n"
                << "  End atom " << A << " bonded with target " << T1 
                << ", bond length is " << _system.distance(A,T1) << " A"
                << "\n  Removed bond is " << T1 << " -- " << T2
                << "\n  In forward move, target atom T1 is picked from " 
                << end_neighbors.size() << " candidates.\n";
            // Perform end-bridging to alter the connectivity of system.
            _system.delete_bond(T1, T2);
            _system.add_bond(A, T1);

            // Starting from the connectivity-altered MDSystem (no alteration to coordinate),
            // determine the number of available targets of the new attacking bead T2,
            // for backward move.
            auto backward_chain_ends = _system.end_atoms();
            auto backward_excludes = tail_atoms(backward_chain_ends);
            auto backward_end_neighbors = _system.find_atoms_within_cutoff(T2, _search_radius);
            remove_crystalline_atoms(backward_end_neighbors);
            exclude_atoms(backward_end_neighbors, backward_excludes);
            auto chain_atoms_of_T2 = determin_chain_atoms(T2);
            exclude_atoms(backward_end_neighbors, chain_atoms_of_T2);
            if (backward_end_neighbors.empty()) {
                std::cout << "Backward move has no possible neighbors.\n";
                continue;
            }
            // Write number of available targets of the new attacking bead T2.
            std::fstream fid("num_candidates", std::ios::out);
            fid << end_neighbors.size() << " " << backward_end_neighbors.size();
            std::cout << "  In backward move, target atom T1 is picked from " 
                      << backward_end_neighbors.size() << " candidates.\n";

            // Since adding/deleting bonds changes molecular configuration, 
            // rebuild the molecule ids here.
            _system.assign_molecules_ids();
            break;
        }
        if (attempt == max_attempt-1) {
            std::cout << "!! Failed to find proper target in "
                      << attempt+1 << " tries\n"; 
        }
    }
}

// Returns root ids of segment containing atom id.
// 1 point -> tail; 2 points -> bridge or loop.
AtomSet MonteCarlo::segment_entry_points(int id) const {
    AtomSet entry_pts; 
    if (is_root(id)) {
        entry_pts.push_back(id);
    }
    AtomSet queue;
    for (int i: _system.bonded_neighbors(id)) {
        if (is_amorphous(i)) queue.push_back(i);
        if (is_root(i)) entry_pts.push_back(i);
    }
    std::unordered_set<int> visited_atoms = {id};
    while (!queue.empty()) {
        auto next = queue.back(); queue.pop_back();
        // If atom is amorphous, it cannot be a root.
        if (is_amorphous(next)) {
            visited_atoms.insert(next);
            for (auto i: _system.bonded_neighbors(next)) {
                if (visited_atoms.count(i) == 0) {
                    queue.push_back(i);
                }
            }
        }
        else if (is_root(next)) {
            entry_pts.push_back(next);
        }
        else {
            std::cerr << "Likely bug in MonteCarlo::segment_entry_points\n";
        }
    }
    return entry_pts;
}

// Return all the atom ids on the current segment in the amorphous phase,
// plus the connected root atom(s).
// Atom id should always be an amorphous atom, not root atom.
AtomSet MonteCarlo::segment_points(int id) const {
    AtomSet seg_pts = {id}; 
    std::unordered_set<int> visited_atoms = {id};
    auto queue = _system.bonded_neighbors(id);
    while (!queue.empty()) {
        auto next = queue.back(); queue.pop_back();
        if (is_amorphous(next)) {
            seg_pts.push_back(next);
            visited_atoms.insert(next);
            for (auto i: _system.bonded_neighbors(next)) {
                if (visited_atoms.count(i) == 0) {
                    queue.push_back(i);
                }
            }
        }
        else if (is_root(next)) seg_pts.push_back(next);
    }
    return seg_pts;
}

// Return all the atoms on the same chain as the atom id.
std::vector<int> MonteCarlo::determin_chain_atoms(int id) const {
    std::set<int> chain_atoms;
    chain_atoms.insert(id);
    auto queue = _system.bonded_neighbors(id);
    while (!queue.empty()) {
        auto next = queue.back(); queue.pop_back();
        if (chain_atoms.count(next) == 0) {
            chain_atoms.insert(next);
        }
        auto neighbors = _system.bonded_neighbors(next);
        for (auto j: neighbors) {
            if (chain_atoms.count(j) == 0) {
                queue.push_back(j);
            }
        }
    }
    return std::vector<int>(chain_atoms.begin(), chain_atoms.end());
}

// Compute the acceptance probablity of the virtual loop/bridge formed.
// A is the previously attacking tail end.
double MonteCarlo::loop_bridge_acceptance_probability(int A, int T1, int T2) const {
    const auto DX = std::max(_geometry.step_x, _geometry.step_y);
    // Entry point of the attacking tail.
    auto root = segment_entry_points(A)[0];    
    // Entry points of target atom.
    auto target_roots = segment_entry_points(T1);
    auto target_chain = shortest_path(target_roots[0], target_roots[1]);
    auto iter = std::find(target_chain.begin(), target_chain.end(), T1);

    auto end_root = target_chain.front();
    if (T2 == *(iter-1)) {
        // T1,T2 is the bond that will be cut, so the last element of the 
        // target chain will be the root on the other end of the bridge/loop.
        end_root = target_chain.back();
    }
    auto x1 = _system[root];
    auto x2 = _system[end_root];
    if (fabs(x1(2)-x2(2)) < _geometry.step_z) {
        // Compute the distance between the entry points of the loop.
        auto dx = norm(_system.bond_image(x1-x2));
        // Adjacent entry points should always have p==1.0;
        return std::min(1.0, exp(-(DX-dx)/DX));
    }    
    return 1.0;  // new segment was a bridge.
}

// Finds the shortest path between atoms id1 and id2.
AtomSet MonteCarlo::shortest_path(int id1, int id2) const {
    if (id1 == id2) return {id1};
    struct Node {
        Node(int i, Node* p) : id(i), parent(p) {} 
        int id;
        Node* parent;
        std::vector<Node> children;
    };
    Node head(id1, nullptr);

    std::list<Node*> queue = {&head};
    while (!queue.empty()) {
        auto n = queue.back(); queue.pop_back();
        const auto &neigh = _system.bonded_neighbors(n->id);

        for (auto j: neigh) {
            // We found the end of the path, now head back up
            // the tree to generate the path.
            if (j == id2) {
                AtomSet path = {id2};
                while (n) {
                    path.push_back(n->id);
                    n = n->parent;
                }
                return AtomSet(path.rbegin(), path.rend());
            }

            if (!n->parent || n->parent->id != j) {
                n->children.emplace_back(j, n);
            }
        }
        for (auto &c: n->children) queue.push_front(&c);
    }
    std::cerr << "No path between atoms " << id1 << " and " << id2 << "\n";
    exit(1);
}

// Returns the number of bonds between two beads.
int MonteCarlo::bonds_between_atoms(int id1, int id2) const {
    return shortest_path(id1, id2).size() - 1;
}

//  Return the id of a neighbor of aim, the bond {aim, neighbor}
//  is to be removed, if there is enough atoms.
//  Otherwise return -1.
int MonteCarlo::select_target_bond(int A, int T) const {
    const auto bond_length = 2.5;
    // Root A has only one entry point.
    auto root_A = segment_entry_points(A);
    // Root T has two entry points, 0,1.
    auto root_T = segment_entry_points(T);
    // Chain length from A to A entry point.
    auto len_A  = bonds_between_atoms(A, root_A[0]) * bond_length;
    // Chain length from T to either T entry point.
    auto len_T0 = bonds_between_atoms(T, root_T[0]) * bond_length;
    auto len_T1 = bonds_between_atoms(T, root_T[1]) * bond_length;
    // Available chain length if A connects through to T0.
    auto chainlen0 = len_A + len_T0;
    // Available chain length if A connects through to T1.
    auto chainlen1 = len_A + len_T1;

    const auto &A0 = _system[root_A[0]];
    const auto &T0 = _system[root_T[0]];
    const auto &T1 = _system[root_T[1]];
    // Distances between entry points on A and T.
    auto d0 = norm(_system.bond_image(A0-T0));
    auto d1 = norm(_system.bond_image(A0-T1));

    if (chainlen0 < _length_factor*d0 && chainlen1 < _length_factor*d1) {
        // Any choice will lead to an overly stretched chain.
        return -1;
    }
    if (chainlen0 < _length_factor*d0) {
        // The bond {T, shortest_path(T, root_T[1])[1]} should be kept,
        // The bond {T, shortest_path(T, root_T[0])[1]} should be removed.
        // Return the atom corresponding to the bond to be removed.
        return shortest_path(T, root_T[0])[1];
    }
    if (chainlen1 < _length_factor*d1) {
        // Return the atom corresponding to the bond to be removed.
        return shortest_path(T, root_T[1])[1];
    }
    return _system.bonded_neighbors(T)[rand()%2];
}

// Returns the numbers of tails, bridges and loops in the amorphous region.
std::vector<AtomSet> MonteCarlo::count_segments(bool report, std::string out) const {
    std::vector<AtomSet> loops, bridges, tails;
    std::vector<FVector<2>> loop_entry_shifts;
    std::set<int> visited_roots;

    AtomSet all_roots = _lower_root_atoms;
    all_roots.insert(all_roots.end(), _upper_root_atoms.begin(), 
                                      _upper_root_atoms.end());
    for (int i: all_roots) {
        if (visited_roots.count(i) == 1) continue;
        visited_roots.insert(i);
        // Finds the first neighbor of i that is amorphous atom.
        int a = -1;
        auto neigh = _system.bonded_neighbors(i);

        // Root is no longer connected to any amorphous atoms.
        // this is probably not reversible.
        if (neigh.size() < 2) {
            std::cout << "Warning, root not connected to amorphous atom.\n";
            tails.push_back({i});
            continue;
        }
        else {
            if (is_amorphous(neigh[0])) {
                a = neigh[0];
            }
            else if (is_amorphous(neigh[1])) {
                a = neigh[1];
            }
            else {
                std::cout << "Error: root is bounded to no amorphous atoms.\n";  
                exit(1);
            }
        }
        auto my_roots = segment_entry_points(a);
        if (my_roots.empty()) {
            std::cerr << "No roots found for atom " << a << ".\n";
            exit(1);
        }
        if (my_roots.size() == 1) {
            tails.push_back(segment_points(a));
        }
        else {
            // Determine types of segments.
            visited_roots.insert(my_roots[0]);
            visited_roots.insert(my_roots[1]);
            if (is_lower_root(my_roots[0]) ^ is_lower_root(my_roots[1])) { 
                bridges.push_back(segment_points(a));
            }
            else {
                loops.push_back(segment_points(a));
                auto shift = loop_entry_shift(a, my_roots);
                loop_entry_shifts.push_back(shift);
            }
        }
    }

    std::fstream f(out, std::ios::out);
    auto write_lammps_id = [&f] (const AtomSet &s) {
        f << "{";
        if (!s.empty()) f << s[0]+1;
        for (int i=1; i<s.size(); ++i) f << ", " << s[i]+1;
        f << "}\n";
    };

    f << "Tails (" << tails.size() << "):\n";
    for (auto &s: tails) {
        write_lammps_id(s);
    }
    f << "\nLoops (" << loops.size() << "):\n";
    f << "Format of loop statistics:\n";
    f << "num_half_cell_x: " << 0.5*_geometry.step_x << "  ";
    f << "num_half_cell_y: " << 0.5*_geometry.step_y << "  ";
    f << "small_shift_distance: " << _geometry.s << "\n";
    for (int i=0; i<loops.size(); ++i) {
        f << loop_entry_shifts[i](0) << " " << loop_entry_shifts[i](1) << " ";
        write_lammps_id(loops[i]);
    }
    f << "\nBridges (" << bridges.size() << "):\n";
    for (auto &s: bridges) {
        write_lammps_id(s);
    }
    if (report) {
        auto total = tails.size()/2 + loops.size() + bridges.size();
        std::cout << "   Number of roots: "     << all_roots.size()
                  << "\n   Number of tails: "   << tails.size()
                  << "\n   Number of bridges: " << bridges.size()
                  << "\n   Number of loops: "   << loops.size()
                  << "\n   Sum of segments: "   << total << "\n"; 
    }
    return loops;
}

// Returns the shift in x and y direction of the two frozen atoms
// connected to the two roots of a loop segment.
FVector<2> MonteCarlo::loop_entry_shift(int id, AtomSet roots) const {
    // Construct the complete loop atom list that can be used
    // to determine the shift.
    AtomSet latoms = segment_points(id);
    // As to loops, the two roots must both be bottom roots or top roots.
    int sign = is_lower_root(roots[0]) ? -1 : 1;
    // Plus 2 is used to ensure that we seek to frozen atoms.
    int num_buffer = _crystal_buffer_depth + 2;
    for (int i=1; i<=num_buffer; ++i) {
        latoms.push_back( roots[0]+sign*i );
        latoms.push_back( roots[1]+sign*i );
    }
    std::sort(latoms.begin(), latoms.end());
    // Sort the latoms due to connectivity.
    // Start from a frozen reference atom.
    auto first = roots[0]+sign*num_buffer;
    AtomSet loop_atoms = {first};
    std::unordered_set<int> visited_atoms = {first};
    auto queue = _system.bonded_neighbors(first);
    while (!queue.empty()) {
        auto next = queue.back(); queue.pop_back();
        if (std::binary_search(latoms.begin(), latoms.end(), next)) {
            loop_atoms.push_back(next);
            visited_atoms.insert(next);
            for (auto i: _system.bonded_neighbors(next)) {
                if (visited_atoms.count(i) == 0) {
                    queue.push_back(i);
                }
            }
        }
    }
    // Walk through the loop_atoms to compute the accumalated shift vector.
    FVector<3> shift = {0.0, 0.0, 0.0};
    for (int i=0; i<loop_atoms.size()-1; ++i) {
        shift += _system.bond_image(loop_atoms[i], loop_atoms[i+1]);
    }
    return {2.0*shift(0)/_geometry.step_x, 2.0*shift(1)/_geometry.step_y};
}

// Returns true if atom is an amorphous atom.
bool MonteCarlo::is_amorphous(int i) const {
    return std::binary_search(_amorphous_atoms.begin(), 
                              _amorphous_atoms.end(), i);
}
// Returns true if atom is an upper root atom.
bool MonteCarlo::is_lower_root(int i) const {
    return std::binary_search(_lower_root_atoms.begin(), 
                              _lower_root_atoms.end(), i);
}
// Returns true if atom is a lower root atom.
bool MonteCarlo::is_upper_root(int i) const {
    return std::binary_search(_upper_root_atoms.begin(), 
                              _upper_root_atoms.end(), i);
}
