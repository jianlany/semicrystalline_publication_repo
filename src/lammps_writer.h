#pragma once
#include <iosfwd>

class MDSystem;

//! Writes the current system to a LAMMPS data file.
void write_to_lammps(const MDSystem &system, std::string lmpdatafile);
