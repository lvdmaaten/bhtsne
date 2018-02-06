
#include <bhtsne/Vector2D.h>

#include <cassert>


namespace bhtsne 
{


template<typename T>
Vector2D<T>::Vector2D(size_t height, size_t width, T initValue)
{
    initialize(height, width, initValue);
}

template<typename T>
Vector2D<T>::Vector2D(Vector2D<T> && other)
{
    m_vector = std::move(other.m_vector);
    m_width = other.m_width;
}

template<typename T>
Vector2D<T>::Vector2D(std::vector<std::vector<double>> vec)
{
    if (vec.size() > 0)
    {
        m_width = vec.front().size();
        for (auto & each : vec)
        {
            appendRow(each);
        }
    }
}

template<typename T>
Vector2D<T> & Vector2D<T>::operator=(Vector2D<T> && other)
{
    m_vector = std::move(other.m_vector);
    m_width = other.m_width;
    return *this;
}

template<typename T>
void Vector2D<T>::initialize(size_t height, size_t width, T initValue)
{
    m_width = width;
    m_vector.resize(height*width, initValue);
}

template<typename T>
void Vector2D<T>::appendRow(std::vector<T>& row)
{
    assert(m_width == row.size());
    m_vector.insert(m_vector.end(), row.begin(), row.end());
}

template<typename T>
size_t Vector2D<T>::size() const
{
    return m_vector.size();
}

template<typename T>
size_t Vector2D<T>::width() const
{
    return m_width;
}

template<typename T>
size_t Vector2D<T>::height() const
{
    assert(size() != 0);
    return size() / width();
}

template<typename T>
typename std::vector<T>::iterator Vector2D<T>::begin()
{
    return m_vector.begin();
}

template<typename T>
typename std::vector<T>::iterator Vector2D<T>::end()
{
    return m_vector.end();
}

template<typename T>
T * Vector2D<T>::operator[](size_t i)
{
    assert(i < height());
    return m_vector.data() + (i * m_width);
}

template<typename T>
const T * Vector2D<T>::operator[](size_t i) const
{
    assert(i < height());
    return m_vector.data() + (i * m_width);
}

template<typename T>
T & Vector2D<T>::at(size_t i, size_t j)
{
    assert(i < height());
    return m_vector.at(i * m_width + j);
}


} //bhtsne