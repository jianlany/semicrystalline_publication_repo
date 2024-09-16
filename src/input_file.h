#pragma once
#include <unordered_map>
#include "string_tools.h"

//! Reads an input file 
/*! Parses a human-readable input file for program options.  Each line in
 *! the input file defines an option. Valid formats are
 *!   key             # defines an option (true/false)
 *!   key value       # defines a key/value relationship (value my have spaces).
 *!   # defines a comment - anything after # is discarded.
 */   
class InputFile {
public:
    //! Parses input file and stores all user input.
    InputFile(std::string path);
    //! Returns option data represented as a string.
    const std::string& operator[](const std::string &key) const;
    //! Returns option data converted to a user-specified type.
    template<typename T> const T get(const std::string &key) const {
        return from_string<T>((*this)[key]);
    }
    //! Returns whether or not an option is set or defined.
    bool is_option_set(const std::string &key) const {
        return _options.count(key);
    }
    //! Allows options to be added (e.g. from command line arguments).
    void add_option(std::string key, std::string value) {
        _options[key] = value;
    }
private:
    std::unordered_map<std::string, std::string> _options;
};
