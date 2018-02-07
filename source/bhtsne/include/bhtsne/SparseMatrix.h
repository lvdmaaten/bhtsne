#pragma once

#include <vector>

namespace bhtsne {
struct SparseMatrix {
    std::vector<double> values;
    std::vector<unsigned int> columns;
    std::vector<unsigned int> rows;
};
}
