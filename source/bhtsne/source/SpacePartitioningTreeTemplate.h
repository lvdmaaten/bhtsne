/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

#pragma once

#include <array>
#include <memory>
#include <bhtsne/Vector2D.h>

namespace bhtsne {

    template<unsigned int D>
    class SpacePartitioningTree
    {
        // Axis-aligned bounding box stored as a center with half-dimensions to represent the boundaries of this quad tree
        std::array<double, D> m_centers;
        std::array<double, D> m_radii;

        // Indices in this space-partitioning tree node, corresponding center-of-mass, and list of all children
        const Vector2D<double> & m_data;
        std::array<double, D> m_centerOfMass;
        unsigned int m_pointIndex;

        // Properties of this node in the tree
        bool m_isLeaf;
        unsigned int m_cumulativeSize;

        // Children
        std::array<std::unique_ptr<SpacePartitioningTree>, 1u << D> m_children;
        unsigned int m_numberOfChildren;

    public:
        explicit SpacePartitioningTree(const Vector2D<double> & data);
        SpacePartitioningTree(const Vector2D<double> & data, const std::array<double, D> & centers,
                              const std::array<double, D> & radii, unsigned int new_index);
        SpacePartitioningTree(const SpacePartitioningTree & other) = delete;
        SpacePartitioningTree(SpacePartitioningTree && other) = default;

        void insert(unsigned int new_index);
        void insertIntoChild(unsigned int new_index);
        unsigned int childIndexForPoint(const double * point);

        // TODO make forces return value
        void computeNonEdgeForces(unsigned int pointIndex, double theta, double * forces, double & forceSum) const;
        void computeEdgeForces(const std::vector<unsigned int> & rows, const std::vector<unsigned int> & columns,
                               const std::vector<double> & values, Vector2D<double> & forces);
    };
}

#include "SpacePartitioningTreeTemplate.inl"