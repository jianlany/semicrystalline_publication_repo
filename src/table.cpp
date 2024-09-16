#include "table.h"
#include "string_tools.h"

Table::Table(std::string tablefile) {
    std::fstream fid(tablefile, std::ios::in);
    if (!fid) {
        std::cerr << "Can't open potential file " << tablefile << "\n";
        exit(1);
    }
    // Read until the first line with a keyword is found.
    while (fid) {
        auto line = trim(before(read_line(fid), "#"));
        if (line.size() > 0) {
            _keyword = line;
            break;
        }
    }
    // {N 1001 R min max}
    auto pieces = split(trim(read_line(fid)));
    int N = from_string<int>(pieces[1]);
    _x.assign(N, 0.0); 
    _y.assign(N, 0.0); 
    _dy.assign(N, 0.0); 
    for (int i=0; i<N; ++i) {
        int id;
        fid >> id >> _x[i] >> _y[i] >> _dy[i];
        if (fid.fail()) {
            std::cout << "Reading table " << tablefile << " failed.\n";
            exit(1);
        }
    } 
};
// Returns the f(x).
double Table::evaluate(double x) const {
    // Avoid crashing, if out of bounds, return first/last value.
    if (x <= _x.front()) return _y.front();
    if (x >= _x.back())  return _y.back();
    auto iter = std::lower_bound(_x.begin(), _x.end(), x);
    auto i = std::distance(_x.begin(), iter)-1;
    return _y[i] + (x-_x[i])*(_y[i+1]-_y[i])/(_x[i+1]-_x[i]);
}
// Returns the value of f'(x).
double Table::evaluate_derivative(double x) const {
    // Avoid crashing, if out of bounds, return first/last value.
    if (x <= _x.front()) return _dy.front();
    if (x >= _x.back())  return _dy.back();
    auto iter = std::lower_bound(_x.begin(), _x.end(), x);
    auto i = std::distance(_x.begin(), iter)-1;
    return _dy[i] + (x-_x[i])*(_dy[i+1]-_dy[i])/(_x[i+1]-_x[i]);
}

