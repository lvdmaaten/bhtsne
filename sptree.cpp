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

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include "sptree.h"

struct SPTree::Node
{
    // For leaf nodes, the single data point it contains
    const double* point;

    // The total number of points in this subtree
    std::size_t size;

    // Axis-aligned bounding box stored as a center with half-dimension widths
    std::vector<double> center;
    const double* width;

    // The center of mass of the points in this subtree
    std::vector<double> center_of_mass;

    // For internal nodes, this node's children
    std::vector<Node*> children;
};

// Default constructor for SPTree -- build tree, too!
SPTree::SPTree(unsigned int D, double* inp_data, unsigned int N)
    : dimension(D),
      data(inp_data),
      buff(D)
{
    // Compute mean, width, and height of current map (boundaries of SPTree)
    double* point = data;
    std::vector<double> mean_Y(D, 0.0);
    std::vector<double> min_Y(D, DBL_MAX);
    std::vector<double> max_Y(D, -DBL_MAX);
    for(unsigned int n = 0; n < N; n++) {
        for(unsigned int d = 0; d < D; d++) {
            double value = point[d];
            mean_Y[d] += value;
            min_Y[d] = std::min(min_Y[d], value);
            max_Y[d] = std::max(max_Y[d], value);
        }
        point += D;
    }
    for(int d = 0; d < D; d++) mean_Y[d] /= N;

    // Construct the tree
    std::vector<double> width(D);
    max_width = 0.0;
    for(int d = 0; d < D; d++) {
        width[d] = std::max(max_Y[d] - mean_Y[d], mean_Y[d] - min_Y[d]) + 1e-5;
        max_width = std::max(max_width, width[d]);
    }
    widths.push_back(std::move(width));

    point = data;
    root = new_node(point, std::move(mean_Y), widths[0].data());

    for(unsigned int i = 1; i < N; i++) {
        point += D;
        insert(root, point);
    }

    for(Node& node : nodes) {
        double scale = 1.0 / node.size;
        for(unsigned int d = 0; d < D; ++d) {
            node.center_of_mass[d] *= scale;
        }
    }
}

// Destructor for SPTree
SPTree::~SPTree() = default;

// Create a new leaf node
SPTree::Node* SPTree::new_node(const double* point, std::vector<double> center, const double* width)
{
    nodes.emplace_back();
    Node* node = &nodes.back();
    node->point = point;
    node->size = 1;
    node->center = std::move(center);
    node->width = width;
    node->center_of_mass.assign(point, point + dimension);
    return node;
}

// Insert a point into the SPTree
void SPTree::insert(Node* node, const double* point)
{
    unsigned int depth = 0;

    while(node->point != point) {
        ++node->size;

        for(unsigned int d = 0; d < dimension; ++d) {
            node->center_of_mass[d] += point[d];
        }

        ++depth;

        if(node->point) {
            // If this is a leaf note, split it into an internal node
            insertChild(node, node->point, depth);
            node->point = nullptr;
        }

        node = insertChild(node, point, depth);
    }
}

// Find the right child node for a point, creating it if necessary
SPTree::Node* SPTree::insertChild(Node* node, const double* point, unsigned int depth) {
    // Find which child to insert into
    unsigned int i = 0;
    for(unsigned int d = 0; d < dimension; ++d) {
        if(point[d] > node->center[d]) {
            i |= 1 << d;
        }
    }

    if(i >= node->children.size()) {
        node->children.resize(1 << dimension, nullptr);
    }

    Node* child = node->children[i];

    if(!child) {
        if(depth >= widths.size()) {
            std::vector<double> width(dimension);
            for(unsigned int d = 0; d < dimension; ++d) {
                width[d] = 0.5 * node->width[d];
            }
            widths.push_back(std::move(width));
        }
        const double* width = widths[depth].data();

        std::vector<double> center(dimension);
        for(unsigned int d = 0; d < dimension; ++d) {
            if(i & (1 << d)) {
                center[d] = node->center[d] + width[d];
            }
            else {
                center[d] = node->center[d] - width[d];
            }
        }

        child = new_node(point, std::move(center), width);
        node->children[i] = child;
    }

    return child;
}

// Compute non-edge forces using Barnes-Hut algorithm
void SPTree::computeNonEdgeForces(unsigned int point_index, double theta, double neg_f[], double* sum_Q)
{
    double* point = data + point_index * dimension;
    computeNonEdgeForces(root, max_width * max_width, point, theta * theta, neg_f, sum_Q);
}

// Compute non-edge forces using Barnes-Hut algorithm
void SPTree::computeNonEdgeForces(Node* node, double max_width_sq, double* point, double theta_sq, double neg_f[], double* sum_Q)
{
    // Make sure that we spend no time on self-interactions
    if(node->point == point) return;

    // Compute distance between point and center-of-mass
    double D = 0.0;
    for(unsigned int d = 0; d < dimension; d++) buff[d] = point[d] - node->center_of_mass[d];
    for(unsigned int d = 0; d < dimension; d++) D += buff[d] * buff[d];

    // Optimize (max_width / sqrt(D) < theta) by squaring and multiplying through by D
    if(node->point || max_width_sq < theta_sq * D) {
        // Compute and add t-SNE force between point and current node
        D = 1.0 / (1.0 + D);
        double mult = node->size * D;
        *sum_Q += mult;
        mult *= D;
        for(unsigned int d = 0; d < dimension; d++) neg_f[d] += mult * buff[d];
    }
    else {
        // Recursively apply Barnes-Hut to children
        for(Node* child : node->children) {
            if(child) {
                computeNonEdgeForces(child, max_width_sq / 4.0, point, theta_sq, neg_f, sum_Q);
            }
        }
    }
}


// Computes edge forces
void SPTree::computeEdgeForces(unsigned int* row_P, unsigned int* col_P, double* val_P, int N, double* pos_f)
{
    
    // Loop over all edges in the graph
    unsigned int ind1 = 0;
    unsigned int ind2 = 0;
    double D;
    for(unsigned int n = 0; n < N; n++) {
        for(unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {
        
            // Compute pairwise distance and Q-value
            D = 1.0;
            ind2 = col_P[i] * dimension;
            for(unsigned int d = 0; d < dimension; d++) buff[d] = data[ind1 + d] - data[ind2 + d];
            for(unsigned int d = 0; d < dimension; d++) D += buff[d] * buff[d];
            D = val_P[i] / D;
            
            // Sum positive force
            for(unsigned int d = 0; d < dimension; d++) pos_f[ind1 + d] += D * buff[d];
        }
        ind1 += dimension;
    }
}

// Print out the tree
void SPTree::print()
{
    print(root);
}

void SPTree::print(Node* node)
{
    if(node->point) {
        printf("Leaf node; data = [");
        for(int d = 0; d < dimension; d++) {
            if(d > 0) {
                printf(", ");
            }
            printf("%f", node->point[d]);
        }
        printf("]\n");
    }
    else {
        printf("Intersection node with center-of-mass = [");
        for(int d = 0; d < dimension; d++) {
            if(d > 0) {
                printf(", ");
            }
            printf("%f", node->center_of_mass[d]);
        }
        printf("]; children are:\n");
        for(Node* child : node->children) {
            if(child) {
                print(child);
            }
        }
    }
}
