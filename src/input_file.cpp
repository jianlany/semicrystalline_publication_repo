#include "input_file.h"

// Parses input file and stores all user input.
InputFile::InputFile(std::string path) {
    std::fstream fid(path, std::ios::in);
    if (!fid) {
        std::cerr << "Unable to open input file " << path << "\n";
        exit(1);
    }
    while (fid) {
        // Reads line stripping off comments and split by whitespace. 
        auto line = split(before(read_line(fid), "#"));
        if (line.empty()) continue;
        else if (line.size() == 1) {
            _options[line[0]] = "true";
        }
        else {
            _options[line[0]] = join(line, 1);
        }
    }
}
// Returns option data represented as a string.
const std::string& InputFile::operator[](const std::string &key) const {
    auto iter = _options.find(key);
    if (iter == _options.end()) {
        std::cerr << "Required parameter {" << key << "} not found\n";
        exit(1);
    }
    else return iter->second;
}
