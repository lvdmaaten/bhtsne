#include <bhtsne/vector2d.h>

using namespace bhtsne;

Vector2D::Vector2D(size_t width, size_t height, double initValue)
        : m_width(width)
        , m_vector(width * height, initValue)
{

}


double * Vector2D::operator[](size_t i) {
    return m_vector.data() + i * m_width;
}

double & Vector2D::at(size_t i, size_t j) {
    return m_vector.at(i * m_width + j);
}