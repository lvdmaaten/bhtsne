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

#include <vector>
#include <memory>
#include <bhtsne/vector2d.h>

namespace bhtsne {
    class SPTree
    {
        struct Cell {
            unsigned int m_dimensions;
            std::vector<double> m_centers;
            std::vector<double> m_radii;

            Cell();
            Cell(unsigned int dimensions, std::vector<double> centers, std::vector<double> radii);
            bool containsPoint(const double* point);
        };

        // Fixed constants
        static const auto QT_NODE_CAPACITY = 1u;

        // Properties of this node in the tree
        unsigned int m_dimensions;
        bool m_isLeaf;
        unsigned int m_cumulativeSize;

        // Axis-aligned bounding box stored as a center with half-dimensions to represent the boundaries of this quad tree
        Cell m_boundary;

        // Indices in this space-partitioning tree node, corresponding center-of-mass, and list of all children
        const Vector2D<double> & m_data;
        std::vector<double> m_centerOfMass;
        std::vector<double> m_pointIndices;

        // Children
        std::vector<std::unique_ptr<SPTree>> m_children;
        unsigned int m_numberOfChildren;

    public:
        SPTree(const Vector2D<double> & data);
        SPTree(const Vector2D<double> & data, std::vector<double> centers, std::vector<double> radii);
        SPTree(const SPTree & other) = delete;
        SPTree(SPTree && other) = default;
        bool insert(unsigned int new_index);
        void computeNonEdgeForces(unsigned int pointIndex, double theta, double forces[], double& forceSum);
        void computeEdgeForces(const std::vector<unsigned int> & rows, const std::vector<unsigned int> & columns, const std::vector<double> & values, Vector2D<double> & forces);

    private:
        void init(std::vector<double> centers, std::vector<double> radii);
        void fill(unsigned int numberOfPoints);
        void subdivide();
    };
}