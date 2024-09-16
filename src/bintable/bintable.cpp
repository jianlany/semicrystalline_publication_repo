#include "bintable.h"
#include "lammps_data.h"
#include "cgsystem.h"
#include <cmath>
#include <set>

namespace cg {

    BinTable::BinTable(double cutoff, int step, const CgSystem *cg) { 
        // Keep pointer reference to cgsystem so we can look up bins later.
        _cgsys       = cg;
        _cutoff      = cutoff;
        _step        = step;
        _granularity = 1;

        // Number of bins along each direction.
        // Granulity to make bins smaller.
        _xbins = int(cg->lmp_data.box_dx(step)/cutoff);
        _ybins = int(cg->lmp_data.box_dy(step)/cutoff);
        _zbins = int(cg->lmp_data.box_dz(step)/cutoff);
        _xbins = std::max(_xbins, 1) * _granularity;
        _ybins = std::max(_ybins, 1) * _granularity;
        _zbins = std::max(_zbins, 1) * _granularity;
        
        _bins.assign(_xbins*_ybins*_zbins, std::vector<int>());

        for (int i=0; i<cg->beads.size(); ++i) {
            _bins[_find_bin(i)].push_back(i);
        }
    }

    // Returns an array of all atoms that can be within cutoff of i.
    // This might not be the best implimentation because a lot of 
    // vectors need to be copied.
    std::vector<const std::vector<int>*> BinTable::neighbors(int i) const {
        // Scaled coordinate of the atom.
        int bi  = _find_bin(i);

        int bx = bi % _xbins;
        int by = (bi/_xbins)%_ybins;
        int bz = (bi/_xbins/_ybins);

        std::set<const std::vector<int>*> bins;

        for (int i=-_granularity; i<=_granularity; ++i) {
            int nx = (bx+i+_xbins) % _xbins;
            for (int j=-_granularity; j<=_granularity; ++j) {                
                int ny = (by+j+_ybins) % _ybins;
                for (int k=-_granularity; k<=_granularity; ++k) {
                    int nz = (bz+k+_zbins) % _zbins;
                    int nb = nx + ny*_xbins + nz*_xbins*_ybins;
                    bins.insert(&_bins[nb]);
                }
            }
        }
        return std::vector<const std::vector<int>*>(bins.begin(), bins.end()); 
    }

    // Returns the index of the bin.
    int BinTable::_find_bin(int i) const {
        double dx = _cgsys->lmp_data.box_dx(_step);
        double dy = _cgsys->lmp_data.box_dy(_step);
        double dz = _cgsys->lmp_data.box_dz(_step);

        auto r = _cgsys->lmp_data.get_bead_coordinate(_step, _cgsys->beads[i].center);
        int bx = std::min(int(double(_xbins)*r(0)/dx), _xbins-1);
        int by = std::min(int(double(_ybins)*r(1)/dy), _ybins-1);
        int bz = std::min(int(double(_zbins)*r(2)/dz), _zbins-1);

        return bx + by*_xbins + bz*_xbins*_ybins;
    }
}
