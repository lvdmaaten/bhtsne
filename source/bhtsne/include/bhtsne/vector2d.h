#pragma once

#include <vector>

namespace bhtsne {

    class Vector2D {
    public:
        Vector2D(size_t width, size_t height, double initValue = 0.0);

        double * operator[](size_t i);

        double & at(size_t i, size_t j);

    protected:
        std::vector<double> m_vector;
        size_t m_width;
    };

}