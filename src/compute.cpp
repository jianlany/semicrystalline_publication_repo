#include "compute.h"
#include "end_bridge.h"
#include "md_system.h"
#include "string_tools.h"
#include "input_file.h"
#include <algorithm>
#include "table.h"

const double avogadros_number = 6.0221413e23;
const double pi = acos(-1.0);

Compute::Compute(const Geometry &g, MDSystem &s, const InputFile &opt) 
  : _geometry(g), 
    _system(s),
    _bead_mass(opt.get<double>("bead_mass")),
    _cutoff_distance(opt.get<double>("cutoff_distance")),
    _bond_stiffness(opt.get<double>("bond_stiffness")),
    _bond_length(opt.get<double>("bond_length")),
    _dz(opt.get<double>("slice_dz_factor")*g.step_z),
    _step(opt.get<int>("step"))
    {}

// Perform all computations.
void Compute::equilibrium_check(std::string pairpath, std::string anglepath,
                                std::string out, std::vector<AtomSet> loops,
                                std::string chain_file, std::string end_file) {
    /* This way of determining number of slices may cause problem
       for the computation for the last slice (the top crystalline slice),
       which may contain less number of atoms than its neighbor slice,
       thus the value for the top crystalline slice may show a decrease.
       However, this way is still better than computing _dz using an
       imput parameter for n_slices, which may cause high fluctuation
       at the slices around the crystalline/amorphous boundary. */
    int n_slices = int(floor(_system.box_length(2)/_dz)) + 1;

    std::fstream fid(out, std::ios::out);
    auto sd = slice_density(n_slices);
    fid << "Slice density, in (z, movable-atom-density, frozen-atom-density, combined-density):\n";
    for (auto &row: sd) {
        fid << std::setw(10) << row(0) << "  "
            << std::setw(10) << row(1) << "  "
            << std::setw(10) << row(2) << "  "
            << std::setw(10) << row(1)+row(2) << "\n";
    }

    auto sbo = slice_bond_orientation(n_slices);
    auto slbo = slice_loop_bond_orientation(n_slices, loops);
    fid << "\nSlice bond orientation, in (z, orientation-all, orientation-loop):\n";
    for (int i=0; i<sbo.first.size(); ++i) {
        fid << sbo.first[i] << "  " << sbo.second[i] 
            << "  " << slbo.second[i] <<"\n";
    }

    auto spe = slice_pair_energy(n_slices, pairpath);
    auto sbe = slice_bond_energy(n_slices);
    auto sae = slice_angle_energy(n_slices, anglepath);
    // Volume of a slice in nm^3
    auto v = _system.box_length(0)*_system.box_length(1)*_dz*1.0e-3;
    fid << "\nSlice energy density, slice volume: " << v << " nm^3\n";
    fid << "(z, pair-energy/nm^3, bond-energy/nm^3, angle-energy/nm^3, total/nm^3):\n";
    for (int i=0; i<spe.first.size(); ++i) {
        fid << spe.first[i] << "  " << spe.second[i] << " "
            << sbe.second[i] << "  " << sae.second[i] << " "
            << spe.second[i]+sbe.second[i]+sae.second[i] << "\n";
    }

    auto e2es = e2e_squared();
    fid << "\nAverage ene2end distance squared <R.R> for system in A^2:\n";
    fid << e2es << "\n";

    write_chains(chain_file);
    write_end_atoms(end_file);
}

// Perform only end2end computations.
void Compute::perform_all_e2e(std::string out) {
    std::fstream fid(out, std::ios_base::app | std::ios_base::out);
    if (_step == 0) fid << "mc-step  <R^2>, in A^2\n";
    auto e2es = e2e_squared();
    fid << _step << "  " << e2es << "\n";
}

// Returns the density distribution along z direction.
// dz is the shift of z for slice-bottom.
std::vector<FVector<3>> Compute::slice_density(int n_slices) const {
    auto mov_atom = read_movable_atom("movable_atoms");
    std::vector<FVector<3>> sd;
    auto box    = _system.bounds();
    auto zshift = box[2](0);
    // Volume of a slice in cm^3
    auto v = _system.box_length(0)*_system.box_length(1)*_dz*1.0e-24;
    for (int k=0; k<n_slices; ++k) {
        auto zlo=k*_dz, zhi=(k+1)*_dz;
        double m_mov=0.0, m_frz=0.0;
        for (int i=0; i<_system.size(); ++i) {
            // Find atoms in the slice {(k-1)*dz, k*dz} along z direction.
            auto zi = _system[i](2) - zshift; // ensure that the minimum zi in system starts from 0.
            if ((zi >= zlo && zi < zhi) ||
                (zi+_system.box_length(2) >= zlo && zi+_system.box_length(2) < zhi)) {
                // contribution from movable atoms.
                if (std::binary_search(mov_atom.begin(), mov_atom.end(), i))
                    m_mov += _bead_mass;
                // contribution from frozen atoms.
                else
                    m_frz += _bead_mass;
            }
        }
        m_mov /= avogadros_number; // convert mass to grams.
        m_frz /= avogadros_number;
        auto row = FVector<3>{0.5*(zlo + zhi), m_mov/v, m_frz/v};
        sd.push_back(row);
    }
    return sd;
}

// Returns the bond orientation distribution along z direction.
// dz is the shift of z for slice-bottom.
// st is the thickness of each slice along z direction.
Slice Compute::slice_bond_orientation(int n_slices) const {
    Slice sbo;
    auto box    = _system.bounds();
    auto zshift = box[2](0);
    auto unit = _geometry.bv/norm(_geometry.bv);
    for (int k=0; k<n_slices; ++k) { 
        double order_parameter = 0.0;
        int counter = 0;
        auto zlo=k*_dz, zhi=(k+1)*_dz;
        for (int i=0; i<_system.size(); ++i) {
            // Find atom i in slice.
            auto zi = _system[i](2) - zshift; // ensure that the minimum zi in system starts from 0.
            if ((zi >= zlo && zi < zhi) ||
                (zi+_system.box_length(2) >= zlo && zi+_system.box_length(2) < zhi)) {
                for (int j: _system.bonded_neighbors(i)) {
                    // Compute the bond second order parameter.
                    auto bond = _system.bond_image(i,j);
                    auto cos_fi = dot(bond,unit)/norm(bond);
                    order_parameter += 0.5 * (3.0*cos_fi*cos_fi - 1.0);
                    ++counter;
                }
            }
        }
        // Both order_parameter and counter are double-counted,
        // but the average below is invariant.
        if (counter > 0) order_parameter /= counter;
        sbo.first.push_back(0.5*(zlo + zhi));
        sbo.second.push_back(order_parameter);
    }
    return sbo;
}

// Returns the bond orientation distribution along z direction
// contributed only by loops.
Slice Compute::slice_loop_bond_orientation(int n_slices,
                             std::vector<AtomSet> loops) const {
    Slice slbo;
    auto unit = _geometry.bv/norm(_geometry.bv);
    //auto unit = FVector<3>{0.0, 0.0, 1.0};
    AtomSet loop_atoms;
    for (auto &loop: loops) {
        for (auto &i: loop) {
            loop_atoms.push_back(i);
        }
    }
    auto box = _system.bounds();
    auto zshift = box[2](0);
    for (int k=0; k<n_slices; ++k) { 
        double order_parameter = 0.0;
        int counter = 0;
        auto zlo=k*_dz, zhi=(k+1)*_dz;
        for (auto &i: loop_atoms) {
            // Find loop atom i in slice.
            auto zi = _system[i](2) - zshift; // ensure that the minimum zi in system starts from 0.
            if ((zi >= zlo && zi < zhi) ||
                (zi+_system.box_length(2) >= zlo && zi+_system.box_length(2) < zhi)) {
                for (int j: _system.bonded_neighbors(i)) {
                    // Compute the bond second order parameter.
                    auto bond = _system.bond_image(i,j);
                    auto cos_fi = dot(bond,unit)/norm(bond);
                    order_parameter += 0.5 * (3.0*cos_fi*cos_fi - 1.0);
                    ++counter;
                }
            }
        }
        // Both order_parameter and counter are double-counted,
        // but the average below is invariant.
        if (counter > 0) order_parameter /= counter;
        slbo.first.push_back(0.5*(zlo + zhi));
        slbo.second.push_back(order_parameter);
    }
    return slbo;
}

// Computes the pair-energy density distribution along z direction.
// dz is the shift of z for slice-bottom.
Slice Compute::slice_pair_energy(int n_slices, std::string pairpath) const {
    auto pair_table = Table(pairpath);
    Slice spe;
    // Volume of a slice in nm^3
    auto v = _system.box_length(0)*_system.box_length(1)*_dz*1.0e-3; 
    //! to do: Find an efficient way to incorporate the periodic slice atom
    //! identation, without looping over all atoms in system.
    auto box = _system.bounds();
    auto zshift = box[2](0);
    for (int k=0; k<n_slices; ++k) {
        int n_atom = 0;
        double slice_pe = 0.0;
        auto zlo=k*_dz, zhi=(k+1)*_dz;
        for (int i=0; i<_system.size(); ++i) {
            // Find atom i within slice.
            auto zi = _system[i](2) - zshift; // ensure that the minimum zi in system starts from 0.
            if ((zi >= zlo && zi < zhi) ||
                (zi+_system.box_length(2) >= zlo && zi+_system.box_length(2) < zhi)) {
                ++n_atom;
                // Note that cut_lo and cut_hi can lay outside the box.
                auto cut_lo = zi - _cutoff_distance;
                auto cut_hi = zi + _cutoff_distance;
                // Search image atoms of j within cutoff distance to atom i.
                double atom_pe = 0.0;
                for (int j=0; j<_system.size(); ++j) {
                    if (j==i || _system.are_bonded(i,j)) continue;

                    auto image_j = _system[i] + _system.bond_image(i, j);
                    if (image_j(2) >= cut_lo && image_j(2) < cut_hi) {
                        auto r = _system.distance(i, j);
                        // Allow half of the pair energy for (i,j) 
                        // contributes to atom i.
                        if (r<_cutoff_distance) {
                            atom_pe += 0.5 * pair_table.evaluate(r);
                        }
                    }
                }
                slice_pe += atom_pe;
            }
        }
        spe.first.push_back(0.5*(zlo + zhi));
        spe.second.push_back(slice_pe/v);
    }
    return spe;
}

// Computes the bond-energy density distribution along z direction.
// dz is the shift of z for slice-bottom.
// st is the thickness of each slice along z direction.
Slice Compute::slice_bond_energy(int n_slices) const {
    Slice sbe;
    // Volume of a slice in nm^3
    auto v = _system.box_length(0)*_system.box_length(1)*_dz*1.0e-3;
    auto box = _system.bounds();
    auto zshift = box[2](0);
    for (int k=0; k<n_slices; ++k) {
        int n_atom = 0;
        double slice_be = 0.0;
        auto zlo=k*_dz, zhi=(k+1)*_dz;
        for (int i=0; i<_system.size(); ++i) {
            // Find atom i within slice.
            auto zi = _system[i](2) - zshift; // ensure that the minimum zi in system starts from 0.
            if ((zi >= zlo && zi < zhi) ||
                (zi+_system.box_length(2) >= zlo && zi+_system.box_length(2) < zhi)) {
                ++n_atom;
                double atom_be = 0.0;
                for (int j: _system.bonded_neighbors(i)) {
                    // Compute the distance between (i,j).
                    auto dr = _system.distance(i, j) - _bond_length;
                    // Allow half of the bond energy for bond (i,j) 
                    // contributes to atom i, no matter if atom j
                    // is in the same slice or not.
                    atom_be += 0.5 * _bond_stiffness * dr * dr;
                }
                slice_be += atom_be;
            }
        }
        sbe.first.push_back(0.5*(zlo + zhi));
        sbe.second.push_back(slice_be/v);
    }
    return sbe;
}

// Computes the angle-energy density distribution along z direction.
// dz is the shift of z for slice-bottom.
// st is the thickness of each slice along z direction.
Slice Compute::slice_angle_energy(int n_slices, std::string anglepath) const {
    auto angle_table = Table(anglepath);
    Slice sae;
    // Volume of a slice in nm^3
    auto v = _system.box_length(0)*_system.box_length(1)*_dz*1.0e-3;
    auto box = _system.bounds();
    auto zshift = box[2](0);
    for (int k=0; k<n_slices; ++k) {
        int n_atom = 0;
        double slice_ae = 0.0;
        auto zlo=k*_dz, zhi=(k+1)*_dz;
        for (int i2=0; i2<_system.size(); ++i2) {
            // Find atom i2 within slice.
            auto zi2 = _system[i2](2) - zshift; // ensure that the minimum zi in system starts from 0.
            if ((zi2 >= zlo && zi2 < zhi) ||
                (zi2+_system.box_length(2) >= zlo && zi2+_system.box_length(2) < zhi)) {
                ++n_atom;
                double atom_ae = 0.0;
                // Find angle atoms i1-i2-i3 all within slice.
                //! \todo, consider branched case.
                if (_system.bonded_neighbors(i2).size() == 2) {
                    auto i1 = _system.bonded_neighbors(i2)[0];
                    auto i3 = _system.bonded_neighbors(i2)[1];
                    if (i1>i3) std::swap(i1, i3);
                    auto theta = bond_angle(i1, i2, i3);
                    // Fail if theta is not a number.
                    if (theta != theta) {
                        std::cerr << "Bond angle was NaN\n";
                        exit(1);
                    }
                    // Allow all the angle energy for (i1,i2,i3) 
                    // contributes to atom i2, no matter if atoms i1 and i3
                    // are in the same alice or not.
                    atom_ae += angle_table.evaluate(theta);
                }
                slice_ae += atom_ae;
            }
        }
        sae.first.push_back(0.5*(zlo + zhi));
        sae.second.push_back(slice_ae/v);
    }
    return sae;
}

// Compute the bond angle between atoms i,j,k.
double Compute::bond_angle(int i, int j, int k) const {
    auto r1 = _system.bond_image(i, j), r2 = _system.bond_image(k, j);
    auto c = dot(r1,r2)/(norm(r1)*norm(r2));
    c = std::min(std::max(c, -1.0), 1.0);
    return acos(c) * 180.0/pi;
}

// Determine the chain ends of individual chains.
std::vector<std::vector<int>> Compute::determin_chain_ends() const {
    // All atoms start as unassigned to a molecule (-1).
    std::vector<int> molecule;
    molecule.assign(_system.size(), -1);

    std::vector<std::vector<int>> chain_ends;
    int current_molecule = 0;
    int num_ends = 0;
    for (int i=0; i<_system.size(); ++i) {
        if (molecule[i] != -1) continue;
        molecule[i] = current_molecule;
        std::vector<int> ends;

        auto queue = _system.bonded_neighbors(i);
        if (queue.size() == 1) {
            ends.push_back(i);
            num_ends += 1;
        }
        while (!queue.empty()) {
            auto next = queue.back(); queue.pop_back();
            molecule[next] = current_molecule;
            auto neighbors = _system.bonded_neighbors(next);
            if (neighbors.size() == 1) {
                ends.push_back(next);
                num_ends += 1;
            }
            for (auto j: neighbors) {
                if (molecule[j] == -1) {
                    queue.push_back(j);
                }
            }
        }
        // Found no more atoms successivly bonded to i, starting a new
        // molecule from the next unassigned molecule.
        current_molecule += 1;
        chain_ends.push_back(ends);
    }
    return chain_ends;
}

// Compute the average end2end distance squared <R.R> for system.
double Compute::e2e_squared() const {
    auto ends = determin_chain_ends();
    double e2es = 0.0;
    int num_infinite_chain = 0;
    int num_linear_chain = 0;
    int num_branch_chain = 0;
    for (auto &e: ends) {
        if (e.size() < 1) {
            num_infinite_chain += 1;
            continue;
        }
        else if (e.size() == 2) {
            num_linear_chain += 1;
            auto x1 = _system[e[0]];
            auto x2 = _system[e[1]];
            // Distance between chain ends x1 and x2.
            auto dx = norm(_system.bond_image(x1-x2));
            e2es += dx*dx;
        }
        else {
            num_branch_chain += 1;
            continue;
        }
    }
    std::cout << "   Find infinite chains: " << num_infinite_chain << "\n";
    std::cout << "   Find linear chains: " << num_linear_chain << "\n";
    std::cout << "   Find branched chains: " << num_branch_chain << "\n\n";
    return e2es/num_linear_chain;
}

//! Read movable atoms.
AtomSet read_movable_atom(std::string movable_atom_path) {
    std::fstream fid(movable_atom_path, std::ios::in);
    if (!fid) {
        std::cerr << "Can't open LAMMPS data file " << movable_atom_path << "\n";
        exit(1);
    }
    AtomSet mov_atom;
    while (fid) {
        auto line = split(read_line(fid));
        for (auto &s: line) {
            int x = std::atoi(s.c_str());
            if (x==0) continue; 
            mov_atom.push_back(x-1);
        }
    }
    return mov_atom;
}

// Write individual chains in format "{chain_id} chain_type (ends)\n beads".
void Compute::write_chains(std::string chain_file) {
    auto num_chain = _system.assign_molecules_ids();
    auto ends = determin_chain_ends();
    // Determine chains which are finite.
    std::vector<int> finite_chains;
    std::fstream fid(chain_file, std::ios::out);
    fid << "All ends:\n";
    for (auto &e: ends) {
        if (e.size() < 1) continue;
        finite_chains.push_back( _system.get_atom_molecule(e[0]) );
        for (auto &i: e) fid << i+1 << " "; 
    }
    std::sort(finite_chains.begin(), finite_chains.end());
    fid << "\n\nChains:\n\n";
    // Chain loop.
    for (int c=0; c<num_chain; ++c) {
        // Write chain id.
        fid << "{" << c << "} ";
        // Write chain type.
        if (std::binary_search(finite_chains.begin(),
                               finite_chains.end(), c)) {
            fid << "finite ( ";
            // Write ends for chain c.
            for (auto &e: ends) {
                if (e.size() < 1) continue;
                if (_system.get_atom_molecule(e[0]) != c) continue;
                for (auto &i: e) fid << i+1 << " "; 
                fid << "):\n";
            }
        }
        else fid << "infinite:\n";
        // Write chain beads.
        for (int i=0; i<_system.size(); ++i) {
            auto chain = _system.get_atom_molecule(i);
            if (chain != c) continue;
            fid << i+1 << " ";
        }
        fid << "\n\n";
    }
}

// Write end-atom coordinates.
void Compute::write_end_atoms(std::string end_file) {
    auto ends = determin_chain_ends();
    std::fstream fid(end_file, std::ios::out);
    auto box    = _system.bounds();
    auto zshift = box[2](0);
    for (auto &e: ends) {
        if (e.size() < 1) continue;
        for (auto &i: e) {
            fid << i+1 << ":  "
                << _system[i](0) << "  "
                << _system[i](1) << "  " 
                << _system[i](2)-zshift << "\n";
        }
    }
}
