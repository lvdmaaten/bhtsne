
#include <algorithm>
#include <cmath>

#include "SpacePartitioningTreeTemplate.h"


using namespace bhtsne;


// Default constructor for SpacePartitioningTree -- build tree, too!
template<unsigned int D>
SpacePartitioningTree<D>::SpacePartitioningTree(const Vector2D<double> & data)
    : m_data(data)
    , m_pointIndex(0)
    , m_isLeaf(true)
    , m_cumulativeSize(1)
    , m_numberOfChildren(1u << D) // = 2^Dimensions
{
    auto numberOfPoints = static_cast<unsigned int>(data.height());
    assert(numberOfPoints > 0);
    // Compute mean, width, and height of current map (boundaries of SpacePartitioningTree)
    std::array<double, D> meanY;
    std::array<double, D> minY;
    std::array<double, D> maxY;
    meanY.fill(0);
    minY.fill(std::numeric_limits<double>::max());
    maxY.fill(std::numeric_limits<double>::min());

    for(unsigned int n = 0; n < numberOfPoints; ++n)
    {
        for(unsigned int d = 0; d < D; ++d)
        {
            auto value = data[n][d];
            meanY[d] += value;
            minY[d] = std::min(minY[d], value);
            maxY[d] = std::max(maxY[d], value);
        }
    }

    for (unsigned int d = 0; d < D; ++d)
    {
        meanY[d] /= numberOfPoints;
    }

    // set boundary
    m_centers = meanY;
    auto delta = 1e-5;
    for (unsigned int d = 0; d < D; ++d)
    {
        m_radii[d] = std::max(maxY[d] - meanY[d], meanY[d] - minY[d]) + delta;
    }

    // take first point for myself
    for (unsigned int d = 0; d < D; ++d)
    {
        m_centerOfMass[d] = data[0][d];
    }
    // insert the rest
    for (auto i = 1u; i < numberOfPoints; ++i)
    {
        insert(i);
    }
}


// Constructor for SpacePartitioningTree with particular size (do not fill the tree)
template<unsigned int D>
SpacePartitioningTree<D>::SpacePartitioningTree(const Vector2D<double> & data, const std::array<double, D> & centers,
                                                const std::array<double, D> & radii, unsigned int new_index)
    : m_centers(centers)
    , m_radii(radii)
    , m_data(data)
    , m_pointIndex(new_index)
    , m_isLeaf(true)
    , m_cumulativeSize(1)
    , m_numberOfChildren(1u << D) // = 2^Dimensions
{
    auto new_point = m_data[new_index];
    for (unsigned int d = 0; d < D; ++d)
    {
        m_centerOfMass[d] = new_point[d];
    }
}


// Insert a point into the SpacePartitioningTree
template<unsigned int D>
void SpacePartitioningTree<D>::insert(unsigned int new_index)
{
    auto new_point = m_data[new_index];

    // Online update of cumulative size and center-of-mass
    m_cumulativeSize++;
    auto avgAdjustment = (m_cumulativeSize - 1.0) / m_cumulativeSize;
    for (unsigned int d = 0; d < D; ++d)
    {
        m_centerOfMass[d] *= avgAdjustment;
        m_centerOfMass[d] += new_point[d] / m_cumulativeSize;
    }

    // leafs must be subdivided
    if (m_isLeaf)
    {
        // Don't add duplicates for now (this is not very nice)
        bool duplicate = true;
        for (unsigned int d = 0; d < D; d++)
        {
            if (new_point[d] != m_data[m_pointIndex][d])
            {
                duplicate = false;
                break;
            }
        }
        if (duplicate)
        {
            return;
        }

        // move my own point to child
        insertIntoChild(m_pointIndex);
        m_isLeaf = false;
    }

    // insert new point into correct child
    insertIntoChild(new_index);
}

template<unsigned int D>
void SpacePartitioningTree<D>::insertIntoChild(unsigned int new_index)
{
    auto childIndex = childIndexForPoint(m_data[new_index]);
    if (!m_children[childIndex])
    {
        std::array<double, D> child_center;
        std::array<double, D> halved_radius;
        for (unsigned int d = 0; d < D; ++d)
        {
            halved_radius[d] = m_radii[d] / 2.0;
            // if the d-th bit is set in the index, the child is below the center in the dimension d
            child_center[d] = (childIndex & (1 << d)) ? m_centers[d] - halved_radius[d] : m_centers[d] + halved_radius[d];
        }
        m_children[childIndex] = std::make_unique<SpacePartitioningTree>(m_data, child_center, halved_radius, new_index);
    }
    else
    {
        m_children[childIndex]->insert(new_index);
    }
}

template<unsigned int D>
unsigned int SpacePartitioningTree<D>::childIndexForPoint(const double * point)
{
    // if the child is below the center in the dimension d, the d-th bit is set in the index
    unsigned int childIndex = 0;
    for (unsigned int d = 0; d < D; ++d)
    {
        if (point[d] < m_centers[d])
        {
            childIndex |= (1 << d);
        }
    }
    return childIndex;
}

// Compute non-edge forces using Barnes-Hut algorithm
template<unsigned int D>
void SpacePartitioningTree<D>::computeNonEdgeForces(unsigned int pointIndex, double theta, double * forces,
                                                    double & forceSum)
{
    // Make sure that we spend no time on empty nodes or self-interactions
    if (m_isLeaf && m_pointIndex == pointIndex)
    {
        return;
    }

    std::array<double, D> distances;
    double sumOfSquaredDistances = 0.0;
    double maxRadius = 0.0;
    for (unsigned int d = 0; d < D; ++d)
    {
        // Compute distance between point and center-of-mass
        distances[d] = m_data[pointIndex][d] - m_centerOfMass[d];
        sumOfSquaredDistances += distances[d] * distances[d];
        maxRadius = std::max(m_radii[d], maxRadius);
    }

    // Check whether we can use this node as a "summary"
    if(m_isLeaf || maxRadius / std::sqrt(sumOfSquaredDistances) < theta)
    {
        // Compute and add t-SNE force between point and current node
        auto inverseDistSum = 1.0 / (1.0 + sumOfSquaredDistances);
        auto force = m_cumulativeSize * inverseDistSum;
        forceSum += force;
        force *= inverseDistSum;
        for (unsigned int d = 0; d < D; ++d)
        {
            forces[d] += force * distances[d];
        }
    }
    else
    {
        // Recursively apply Barnes-Hut to children
        for (auto & child : m_children)
        {
            // TODO: prevent invalid children in m_children
            if(child)
            {
                child->computeNonEdgeForces(pointIndex, theta, forces, forceSum);
            }
        }
    }
}


// Computes edge forces
template<unsigned int D>
void SpacePartitioningTree<D>::computeEdgeForces(const std::vector<unsigned int> & rows,
                                                 const std::vector<unsigned int> & columns,
                                                 const std::vector<double> & values,
                                                 Vector2D<double> & forces)
{
    // Loop over all edges in the graph
    std::array<double, D> distances;
    for(unsigned int n = 0; n < forces.height(); ++n)
    {
        for(auto i = rows[n]; i < rows[n + 1]; ++i)
        {
            // Compute pairwise distance and Q-value
            double sumOfSquaredDistances = 1.0;
            for (unsigned int d = 0; d < D; ++d)
            {
                distances[d] = m_data[n][d] - m_data[columns[i]][d];
                sumOfSquaredDistances += distances[d] * distances[d];
            }
            double force = values[i] / sumOfSquaredDistances;

            // Sum positive force
            for(unsigned int d = 0; d < D; ++d)
            {
                forces[n][d] += force * distances[d];
            }
        }
    }
}
