#pragma once
#include <string>
#include <ctype.h>
#include "fvector.h"

class MDSystem;
class Geometry;
class InputFile;

//! Stores coordinate set.
typedef std::vector<FVector<3>> CS;
typedef std::vector<int> AtomSet;
typedef std::vector<double> SlipSet;

//! Post processing of slip idemtation in semicrystalline system.
class Slip {
public:
    //! Sets a reference to the MDSystem.
    Slip(const Geometry &g, MDSystem &s0, MDSystem &s1, const InputFile &opt);
    //! Main function.
    void compute_slip(std::string out);
    //! Construct a list to the nth bonded neighbor of atom i in the reference frame.
    AtomSet segment_around(int i, int n);
    //! Compute slip vector for atoms at successsive time frames.
    CS slip_vector(double cutoff, int n_neighbor);


private:
    //! Reference to domain geometry parameters.
    const Geometry &_geometry;
    //! This class will manipulate the bonds within the MDSystem.
    MDSystem &_sys0;
    MDSystem &_sys1;
    //! Path of trajectory file.
    std::string _trjfile;
    //! CG system parameters.
    double _cutoff_distance;
    double _n_neighbor, _crystallinity;
    //! MC step index.
    int _step;
};

//! Detect all the time frames in the trajectory file.
AtomSet time_frame(std::string trjpath);
//! Reads coordinates of time frame i from the specified lammps trajectory file.
void read_trajectory(std::string trjpath, std::vector<int> frames,
                     int i, MDSystem &sys);
//! Write/append a single time frame to the output trajectory file.
void output_trajectory(std::string out, int frame, MDSystem &sys,
                       CS svs, double crystallinity);
