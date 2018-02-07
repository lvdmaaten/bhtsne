#pragma once

#include <vector>

namespace bhtsne {


template<typename T>
class Vector2D {
public:
    Vector2D() = default;
    Vector2D(size_t height, size_t width, T initValue = 0.0);
    Vector2D(Vector2D<T> && other);
    Vector2D(std::vector<std::vector<double>> vec);

    Vector2D<T> & operator=(Vector2D<T> && other);

    void initialize(size_t height, size_t width, T initValue = 0.0);
    void appendRow(std::vector<T>& row);

    size_t size() const;
    size_t width() const;
    size_t height() const;

    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();

    T * operator[](size_t i);
    const T * operator[](size_t i) const;
    T & at(size_t i, size_t j);

protected:
    std::vector<T> m_vector;
    size_t m_width;
};


}

#include <bhtsne/Vector2D.inl>