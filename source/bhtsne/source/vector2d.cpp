
#include <bhtsne/vector2d.h>

#include <cassert>


using namespace bhtsne;


template<typename T>
Vector2D<T>::Vector2D(std::vector<T>& firstRow)
{
    m_width = firstRow.size();
    appendRow(firstRow);
}

template<typename T>
Vector2D<T>::Vector2D(size_t height, size_t width, T initValue)
: m_width(width)
{
    m_vector.resize(height*width, initValue);
}

template<typename T>
void Vector2D<T>::appendRow(std::vector<T>& row)
{
    assert(m_width == row.size());
    m_vector.insert(m_vector.end(), row.begin(), row.end());
}

template<typename T>
T * Vector2D<T>::operator[](size_t i) {
    return m_vector.data() + i * m_width;
}

template<typename T>
T & Vector2D<T>::at(size_t i, size_t j) {
    return m_vector.at(i * m_width + j);
}