#include "slip.h"
#include "end_bridge.h"
#include "md_system.h"
#include "string_tools.h"
#include "input_file.h"
#include <algorithm>
#include "table.h"
#include <unordered_set>
#include <unordered_map>
#include <set>

const double avogadros_number = 6.0221413e23;
const double pi = acos(-1.0);

Slip::Slip(const Geometry &g, MDSystem &s0, MDSystem &s1, const InputFile &opt) 
  : _geometry(g), 
    _sys0(s0),
    _sys1(s1),
    _trjfile("slip/step_"+opt["step"]+".lammpstrj"),
    _cutoff_distance(opt.get<double>("cutoff_distance")),
    _n_neighbor(opt.get<int>("n_neighbor")),
    _crystallinity(opt.get<double>("crystallinity")),
    _step(opt.get<int>("step"))
    {}


// Constructs a list to the nth bonded neighbor of atom i in the reference frame.
// e.g. if n = 2, then it will return the sequence { i-2, i-2, i, i+1, i+2 }.
AtomSet Slip::segment_around(int i, int n) {
    AtomSet segment = {i};
    int start=0, stop=1;
    for (int layer=0; layer<n; ++layer) {
        for (int j=start; j!=stop; ++j) {
            auto neigh = _sys0.bonded_neighbors(segment[j]);
            segment.insert(segment.end(), neigh.begin(), neigh.end());
        }
        start = stop;
        stop  = segment.size();
    }
    std::sort(segment.begin(), segment.end());
    segment.erase(std::unique(segment.begin(), segment.end()), segment.end());
    return segment;
}


// Perform slip computations.
//! _sys0 ~ reference, _sys1 ~ current.
void Slip::compute_slip(std::string out) {
    auto frames = time_frame(_trjfile);
    // Reference MDsystem.
    read_trajectory(_trjfile, frames, 0, _sys0);
    for (int frame_id=0; frame_id<frames.size()-2; frame_id+=1) {
        // Current MDsystem.
        read_trajectory(_trjfile, frames, frame_id+1, _sys1);
        std::cout << "Read frames 0, " << frame_id+1 << ".  ";
        // Compute slip vectors.
        auto svs = slip_vector(8.0, _n_neighbor);
        // Write/append the norm of slip vectors as the vx column
        // in a new trajectory file.
        if (frame_id==0) {
            // Set the slip vectors for the initial frame as zeros.
            CS svs0;
            svs0.assign(_sys0.size(), zeros<3>());
            output_trajectory(out, frames[0], _sys0, svs0, _crystallinity);
            output_trajectory(out, frames[1], _sys1, svs, _crystallinity);
            std::cout << "Output frames " << frame_id << ", " << frame_id+1 << ".\n";
        }
        else {
            output_trajectory(out, frames[frame_id+1], _sys1, svs, _crystallinity);
            std::cout << "Output frame " << frame_id+1 << ".\n";
        }
    }
}


//! Compute slip vector for all atoms at successsive time frames.
//! _sys0 ~ reference, _sys1 ~ current.
CS Slip::slip_vector(double cutoff, int n_neighbor) {
    CS svs;
    struct NearAtom { int id; double distance; };
    auto natom = _sys0.size();
    for (int i=0; i<natom; ++i) {
        // The beads bonded to bead i on the same chain.
        auto exclude = segment_around(i, 3);
        // Finds atoms within a cutoff distance to atom i
        // in the reference frame _sys0.
        std::vector<NearAtom> pair_neighbors;
        for (int j=0; j<natom; ++j) {
            if (j==i) continue;
            if (std::count(exclude.begin(), exclude.end(), j)) continue;
            // Compute the distance between atom i and image atom j.
            // Note that image_j = _sys[i] + _sys.bond_image(i, j)
            auto d = _sys0.distance(i, j);
            if (d<cutoff) {
                pair_neighbors.push_back(NearAtom{j, d});
            }
            if (pair_neighbors.size()==n_neighbor) break;
        }
        // Sorts i's neighbors by distance.
        std::sort(pair_neighbors.begin(), pair_neighbors.end(),
                  [](const NearAtom &a, const NearAtom &b) {
                  return a.distance < b.distance;});
        // Compute slip vector sv for atom i.
        auto s = pair_neighbors.size();
        auto n_loop = (n_neighbor < s) ? n_neighbor : s;
        auto sv = zeros<3>();
        for (int k=0; k<n_loop; ++k) {
            auto jj = pair_neighbors[k].id;
            auto Xij = _sys0.bond_image(i, jj);
            auto xij = _sys1.bond_image(i, jj);
            auto diff = Xij - xij;
            sv += diff;
        }
        svs.push_back(sv/n_loop);
    }
    return svs;
}


// Detect all the time frames in the trajectory file.
AtomSet time_frame(std::string trjpath) {
    std::fstream fid(trjpath, std::ios::in);
    if (!fid) {
        std::cerr << "Can't open LAMMPS trajectory file " << trjpath << "\n";
        exit(1);
    }
    std::vector<int> frames;
    int frame = 0;
    while (fid) {
        auto line = read_line(fid);
        if (line.find("ITEM: TIMESTEP") != std::string::npos) {
            line = read_line(fid);
            frame = from_string<int>(split(line)[0]);
            frames.push_back(frame);
        }
    }
    return frames;
}


// Reads coordinates of time frame i from the specified lammps trajectory file.
void read_trajectory(std::string trjpath, std::vector<int> frames,
                   int i, MDSystem &sys) {
    std::fstream fid(trjpath, std::ios::in);
    if (!fid) {
        std::cerr << "Can't open LAMMPS trajectory file " << trjpath << "\n";
        exit(1);
    }
    // Search for the i-th time frame.
    int frame = 0;
    while (fid) {
        auto line = read_line(fid);
        if (line.find("ITEM: TIMESTEP") != std::string::npos) {
            line = read_line(fid);
            frame = from_string<int>(split(line)[0]);
            if (frame == frames[i]) {
                for (int j=0; j<3; ++j)
                    read_line(fid);
                // Set box dimensions.  
                FVector<2> box[3];
                for (int i=0; i<=2; ++i) {
                    auto data = split(read_line(fid));
                    box[i](0) = from_string<double>(data[0]);
                    box[i](1) = from_string<double>(data[1]);
                }  
                sys.set_bounds(box[0](0), box[0](1),
                               box[1](0), box[1](1),
                               box[2](0), box[2](1));
                read_line(fid);
                break;
            }
        }
    }
    // Read coordinates in time frame i.
    AtomSet mols;
    AtomSet types;
    CS coords;
    SlipSet vms, css;
    auto natoms = sys.size();
    mols.assign(natoms, -1);
    types.assign(natoms, -1);
    coords.assign(natoms, zeros<3>());
    vms.assign(natoms, 0.0);
    css.assign(natoms, 0.0);
    while (fid) {
        auto line = read_line(fid);
        if (line.find("ITEM: TIMESTEP") != std::string::npos) break;
        // id mol type x y z
        auto data = split(line);
        // Store atom-id as zero-indexed.
        int id = from_string<int>(data[0]) - 1;
        mols[id]   = from_string<int>(data[1]);
        types[id]  = from_string<int>(data[2]);
        coords[id] = {from_string<double>(data[3]),
                      from_string<double>(data[4]),
                      from_string<double>(data[5])};
        //vms[id]    = from_string<double>(data[6]);
        //css[id]    = from_string<double>(data[7]);
    }
    sys.update_atom_mol(mols);
    sys.update_atom_type(types);
    sys.update_atom_coordinate(coords);
    sys.warp_atoms();
    sys.update_atom_vm(vms);
    sys.update_atom_cs(css);
}


// Write/append a single time frame to the output trajectory file.
void output_trajectory(std::string out, int frame, MDSystem &sys,
                       CS svs, double crystallinity) {
        std::fstream fid(out, std::ios_base::app | std::ios_base::out);
        fid << "ITEM: TIMESTEP\n";
        fid << frame << "\n";
        fid << "ITEM: NUMBER OF ATOMS\n";
        fid << sys.size() << "\n";
        fid << "ITEM: BOX BOUNDS pp pp pp\n";
        for (int i=0; i<=2; ++i) {
            fid << sys.bounds(i)(0) << " " << sys.bounds(i)(1) << "\n";
        }
        // vx ~ the magnitude of slip vector for each atom.
        // vy ~ the von Mises stress for each atom.
        // vz ~ the centro-symmetry parameter for each atom.
        fid << "ITEM: ATOMS id mol type x y z vx vy vz\n"; 
        for (int i=0; i<sys.size(); ++i) {
            // Approximately determine if the current atom is crystalline or amorphous.
            auto lo = sys.bounds(2)(0) + 0.5*sys.box_length(2)*crystallinity;
            auto hi = sys.bounds(2)(1) - 0.5*sys.box_length(2)*crystallinity;
            int type;
            if (sys.bonded_neighbors(i).size()==1) type = 0; // end atom
            else if (sys[i](2)>lo && sys[i](2)<hi) type = 1; // amorphous atom
            else                                   type = 2; // crystalline atom
            // Output
            fid << i+1 << " " 
                << sys.mol(i) << " "
                << type << " "
                << sys[i](0) << " "
                << sys[i](1) << " "
                << sys[i](2) << " "
                << norm(svs[i]) << " "
                << sys.vm(i) << " "
                << sys.cs(i) << "\n";
        }
}
