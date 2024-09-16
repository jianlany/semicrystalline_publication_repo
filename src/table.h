#pragma once
#include <string>
#include <vector>

//! Implments a tabulated function from LAMMPS table style.
class Table {
public:
    Table(std::string tablefile);
    //! Returns the first value of x.
    double x_min() const { return _x.front(); }
    //! Returns the f(x).
    double evaluate(double x) const;
    //! Returns the value of f'(x).
    double evaluate_derivative(double x) const;
private:
    //! A keyword associated with this table value.
    std::string _keyword;
    //! The independent variable.
    std::vector<double> _x;
    //! The dependent variable.
    std::vector<double> _y;
    //! The derviative of the dependent variable.
    std::vector<double> _dy;
};
