#pragma once

#include <vector>

namespace bhtsne {


template<typename T>
class Vector2D {
public:
    Vector2D(std::vector<T>& firstRow);
    Vector2D(size_t height, size_t width, T initValue = 0.0);

    void appendRow(std::vector<T>& row);

    T * operator[](size_t i);
    T & at(size_t i, size_t j);

protected:
    std::vector<T> m_vector;
    size_t m_width;
};


}