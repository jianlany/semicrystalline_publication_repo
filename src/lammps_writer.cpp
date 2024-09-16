#include "lammps_writer.h"
#include "md_system.h"
#include <boost/format.hpp>
#include <set>
#include <fstream>


// Writes the current system to a LAMMPS data file.
void write_to_lammps(const MDSystem &system, std::string lmpdatafile) {
    using boost::format;

    // Makes a vector of unique atom types.
    std::set<int> types;
    for (auto t: system._atom_types) {
        types.insert(t);
    }
    std::vector<int> lmp_atom_types(types.begin(), types.end());

    // Stores bonds as {i1, i2, type}.
    std::vector<std::tuple<int,int,int>> bonds;

    // Stores angles as {i1, i2, i3, type}.
    std::vector<std::tuple<int,int,int,int>> angles;
    std::set<std::tuple<int,int,int>> angle_types;
    for (int i=0; i<system._bond_table.size(); ++i) {
        for (int j: system._bond_table[i]) {
            if (i < j) {
                //! \todo Determine bond type.
                int type = 1;
                bonds.emplace_back(i,j,type);
            }
        }

        // If an atom is bonded to more than one atom it makes an angle.
        if (system._bond_table[i].size() == 2) {
            int i1 = system._bond_table[i][0];
            int i3 = system._bond_table[i][1];
            if (i1>i3) std::swap(i1, i3);
            //! \todo Determine angle type. 
            int type = 1;
            angles.emplace_back(i1, i, i3, type);

            int i1type = system._atom_types[i1];
            int i2type = system._atom_types[i];
            int i3type = system._atom_types[i3];

            // cba -> abc (symmetric angle types).
            if (i1type > i3type) std::swap(i1type, i3type);
            angle_types.insert(std::make_tuple(i1type, i2type, i3type));
        }
    }

    std::fstream fid(lmpdatafile, std::ios::out);
    fid << "LAMMPS 2005 data file built by pppp\n\n";
    fid << format(" %7d atoms\n") % system._atom_types.size();
    fid << format(" %7d bonds\n") % bonds.size();
    fid << format(" %7d angles\n") % angles.size();
    fid << format("%4d atom types\n") % lmp_atom_types.size();
    fid << format("%4d bond types\n") % system._bond_types.size();
    fid << format("%4d angle types\n") % angle_types.size();
    fid << "\n";

    auto &step = system._step;
    fid << format(" %15.9f %15.9f xlo xhi\n") % step.box[0](0) % step.box[0](1);
    fid << format(" %15.9f %15.9f ylo yhi\n") % step.box[1](0) % step.box[1](1);
    fid << format(" %15.9f %15.9f zlo zhi\n") % step.box[2](0) % step.box[2](1);

    fid << "\nMasses\n\n";
    for (int i=0; i<lmp_atom_types.size(); ++i) {
        //! \todo Get actual atom mass.
        double m = 28.0;
        fid << format("%d %f\n") %(i+1) %m;
        
    }
    fid << "\nAtoms\n\n";
    for (int i=0; i<step.atom_coords.size(); ++i) {
        auto &x = step.atom_coords[i];
        auto iter = std::find(lmp_atom_types.begin(), lmp_atom_types.end(), 
                              system._atom_types[i]);
        auto type = std::distance(lmp_atom_types.begin(), iter) + 1;
        // id mol type x y z
        fid << format("%6d %6d %3d %.12f %.12f %.12f\n")
                       %(i+1) %(system.get_atom_molecule(i)) %type %x(0) %x(1) % x(2);
    }

    fid << "\nVelocities\n\n";
    for (int i=0; i<step.atom_velcs.size(); ++i) {
        auto &v = system.velocity(i);
        // id vx vy vz
        fid << format("%6d %.12f %.12f %.12f\n")
                       %(i+1) %v(0) %v(1) % v(2);
    }

    fid << "\nBonds\n\n";
    for (int i=0; i<bonds.size(); ++i) {
        int i1,i2,type; std::tie(i1,i2,type) = bonds[i];
        fid << format("%6d %3d %6d %6d\n")
                       %(i+1) %type %(i1+1) %(i2+1);
    }

    fid << "\nAngles\n\n";
    for (int i=0; i<angles.size(); ++i) {
        int i1,i2,i3,type; std::tie(i1,i2,i3,type) = angles[i];
        fid << format("%6d %3d %6d %6d %6d\n")
                       %(i+1) %type %(i1+1) %(i2+1) %(i3+1);
    }
}

