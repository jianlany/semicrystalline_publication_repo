#include "geometry.h"
#include "input_file.h"

// Initializes geometry parameters from command line options.
Geometry::Geometry(const InputFile &opt) 
 : a(opt.get<double>("a")),
   b(opt.get<double>("b")),
   c(opt.get<double>("c")),
   nx(opt.get<int>("nx")),
   ny(opt.get<int>("ny")),
   nz(opt.get<int>("nz")),
   lo_amorphous_z(int(nz*opt.get<double>("crystallinity"))/2),
   hi_amorphous_z(nz - lo_amorphous_z)
{
    // Seeds the random number generator.
    auto seed = from_string<unsigned>(opt["seed"]);
    seed     += from_string<unsigned>(opt["step"]);
    srand(seed);

    // Title angle between chain stems and the lamellar normal.
    // With Miller index (201).   
    auto theta = atan(2.0*c/a);
    step_x = a/cos(theta);  // horizontal
    step_y = b;             // inward
    step_z = c*cos(theta);  // vertical
    bv  = {-c*sin(theta), 0.0, c*cos(theta)};
    s   = norm(0.5*bv)*sin(theta);
}

// Returns the lower and upper z coordinates of the amorphous layer.
std::pair<double,double> Geometry::amorphous_region() const {
    return {lower_root_coordinate() + 1e-6*step_z,
            upper_root_coordinate() - 1e-6*step_z};
}
