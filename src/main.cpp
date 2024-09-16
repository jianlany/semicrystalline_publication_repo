#include <iostream>
#include "input_file.h"
#include "md_system.h"
#include "builder.h"
#include "end_bridge.h"
#include "lammps_writer.h"
#include "geometry.h"
#include "compute.h"
#include "slip.h"
#include "timer.h"

//! Performs a single MC step to build the system.
void build_system(const InputFile &opt);
//! Performs post processing for a step.
void post_process(const InputFile &opt);
//! Postprocessing, computes slip in system.
void post_slip(const InputFile &opt);

//! Main program reads input args and executes program mode.
int main(int argc, char **argv) {

    std::string modified_files = MODIFIED_FILES;
    if (modified_files.empty()) {
        std::cout << "Running with clean commit: " << GIT_REVISION << "\n";
    }
    else {
        std::cout << "Running with dirty commit " << GIT_REVISION << "\n";
        std::cout << "Files with modifications " 
                  << split(modified_files) << "\n";
    }

    if (argc < 2) {
        std::cerr << "Not enough input arguments.\n";
        std::cerr << "molbuilder (build|post|slip) [step] [inputfile] [current_stage]\n";
        exit(1);
    }
    // Input file is always arg 3.
    auto mode = std::string(argv[1]);
    auto opt  = InputFile(argv[3]);
    opt.add_option("step", argv[2]);
    opt.add_option("current_stage", argv[4]);

    if (mode == "build") {
        build_system(opt);
    }
    else if (mode == "post") {
        post_process(opt);
    }
    else if (mode == "slip") {
        post_slip(opt);
    }
    else {
        std::cerr << "Invalid mode.\n";
        exit(1);
    }
}

//! Generates trial.lammps (either first step or end bridge).
void build_system(const InputFile &opt) {
    auto system   = MDSystem();
    auto geometry = Geometry(opt);
    if (from_string<int>(opt["step"]) == 0) {
        auto builder = Builder(geometry, system, opt);
        builder.generate();
        builder.write_amorphous_atoms("amorphous_atoms.txt");
        builder.write_root_atoms("root_atoms.txt");
        builder.write_movable_atoms("movable_atoms");
    }
    else {
        read_data("trial.lammps", system);
        auto mc = MonteCarlo(geometry, system, opt);
        mc.end_bridge();
    }
    write_to_lammps(system, "trial_EB_start.lammps");
}

//! Computes post processing steps.
/*! 
 *! If e2e option is set, runs quick end-to-end calculation,
 *! otherwise, runs the whole set of post processing calcualtions.
 */
void post_process(const InputFile &opt) {
    const char* pair_path = "../pair.table.EE.44";
    const char* angle_path = "../angle.table.EEE.44";

    auto system   = MDSystem();
    auto geometry = Geometry(opt);
    auto lmp_data = "data/step_" + opt["step"] + ".lammps";
    read_data(lmp_data, system);
    auto compute = Compute(geometry, system, opt);
    if (opt.is_option_set("e2e")) {
        auto out_file = "postprocessing/e2e.txt";
        compute.perform_all_e2e(out_file);
    }
    else {
        auto mc = MonteCarlo(geometry, system, opt);
        auto out_file = "postprocessing/post_data_" + opt["step"] + ".txt";
        auto loops = mc.count_segments(true,
                        "postprocessing/lbt_atoms_" + opt["step"] + ".txt");
        auto chain_file = "postprocessing/chains_step_" + opt["step"] + ".txt";
        auto end_file   = "postprocessing/ends_step_" + opt["step"] + ".txt";
        compute.equilibrium_check(pair_path, angle_path, out_file,
                                  loops, chain_file, end_file);
    }
}

//! Postprocessing, computes slip in system.
void post_slip(const InputFile &opt) {
    auto geometry = Geometry(opt);
    auto lmp_data = "slip/step_" + opt["step"] + ".lammps";
    // MD systems at two successive time frames.
    // system0 ~ reference, system1 ~ current.
    auto system0   = MDSystem();
    auto system1   = MDSystem();
    // Set the same bond information to system0 and system1.
    // The coordinate and box informatoin will be modified later.
    read_data(lmp_data, system0);
    read_data(lmp_data, system1);
    auto slip = Slip(geometry, system0, system1, opt);
    auto out_file = "slip/step_" + opt["step"] + "_slip_vector.lammpstrj";
    slip.compute_slip(out_file);
}
