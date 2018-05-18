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



// Constructs cell
Cell::Cell(unsigned int inp_dimension)
    : corner(inp_dimension),
      width(inp_dimension)
{
}

Cell::Cell(unsigned int inp_dimension, double* inp_corner, double* inp_width)
    : corner(inp_dimension),
      width(inp_dimension)
{
    std::copy_n(inp_corner, inp_dimension, corner.begin());
    std::copy_n(inp_width, inp_dimension, width.begin());
}

double Cell::getCorner(unsigned int d) {
    return corner[d];
}

double Cell::getWidth(unsigned int d) {
    return width[d];
}

void Cell::setCorner(unsigned int d, double val) {
    corner[d] = val;
}

void Cell::setWidth(unsigned int d, double val) {
    width[d] = val;
}

// Checks whether a point lies in a cell
bool Cell::containsPoint(double point[])
{
    for(unsigned int d = 0; d < corner.size(); d++) {
        if(corner[d] - width[d] > point[d]) return false;
        if(corner[d] + width[d] < point[d]) return false;
    }
    return true;
}


// Default constructor for SPTree -- build tree, too!
SPTree::SPTree(unsigned int D, double* inp_data, unsigned int N)
{

    // Compute mean, width, and height of current map (boundaries of SPTree)
    int nD = 0;
    std::vector<double> mean_Y(D, 0.0);
    std::vector<double> min_Y(D, DBL_MAX);
    std::vector<double> max_Y(D, -DBL_MAX);
    for(unsigned int n = 0; n < N; n++) {
        for(unsigned int d = 0; d < D; d++) {
            double value = inp_data[nD + d];
            mean_Y[d] += value;
            min_Y[d] = std::min(min_Y[d], value);
            max_Y[d] = std::max(max_Y[d], value);
        }
        nD += D;
    }
    for(int d = 0; d < D; d++) mean_Y[d] /= (double) N;

    // Construct SPTree
    std::vector<double> width(D);
    for(int d = 0; d < D; d++) width[d] = std::max(max_Y[d] - mean_Y[d], mean_Y[d] - min_Y[d]) + 1e-5;
    init(NULL, D, inp_data, mean_Y.data(), width.data());
    fill(N);
}


// Constructor for SPTree with particular size and parent -- build the tree, too!
SPTree::SPTree(unsigned int D, double* inp_data, unsigned int N, double* inp_corner, double* inp_width)
{
    init(NULL, D, inp_data, inp_corner, inp_width);
    fill(N);
}


// Constructor for SPTree with particular size (do not fill the tree)
SPTree::SPTree(unsigned int D, double* inp_data, double* inp_corner, double* inp_width)
{
    init(NULL, D, inp_data, inp_corner, inp_width);
}


// Constructor for SPTree with particular size and parent (do not fill tree)
SPTree::SPTree(SPTree* inp_parent, unsigned int D, double* inp_data, double* inp_corner, double* inp_width) {
    init(inp_parent, D, inp_data, inp_corner, inp_width);
}


// Constructor for SPTree with particular size and parent -- build the tree, too!
SPTree::SPTree(SPTree* inp_parent, unsigned int D, double* inp_data, unsigned int N, double* inp_corner, double* inp_width)
{
    init(inp_parent, D, inp_data, inp_corner, inp_width);
    fill(N);
}


// Main initialization function
void SPTree::init(SPTree* inp_parent, unsigned int D, double* inp_data, double* inp_corner, double* inp_width)
{
    parent = inp_parent;
    dimension = D;
    data = inp_data;
    is_leaf = true;
    size = 0;
    cum_size = 0;

    boundary = Cell(dimension);
    for(unsigned int d = 0; d < D; d++) boundary.setCorner(d, inp_corner[d]);
    for(unsigned int d = 0; d < D; d++) boundary.setWidth( d, inp_width[d]);

    double max_width = 0.0;
    for(unsigned int d = 0; d < D; d++) {
        max_width = std::max(max_width, boundary.getWidth(d));
    }
    max_width_sq = max_width * max_width;

    unsigned int no_children = 1 << D;
    children.resize(no_children, nullptr);

    center_of_mass.resize(D, 0.0);

    buff.resize(D);
}


// Destructor for SPTree
SPTree::~SPTree()
{
    for(SPTree* child : children) {
        delete child;
    }
}


// Update the data underlying this tree
void SPTree::setData(double* inp_data)
{
    data = inp_data;
}


// Get the parent of the current tree
SPTree* SPTree::getParent()
{
    return parent;
}


// Insert a point into the SPTree
bool SPTree::insert(unsigned int new_index)
{
    // Ignore objects which do not belong in this quad tree
    double* point = data + new_index * dimension;
    if(!boundary.containsPoint(point))
        return false;
    
    // Online update of cumulative size and center-of-mass
    cum_size++;
    double mult1 = (double) (cum_size - 1) / (double) cum_size;
    double mult2 = 1.0 / (double) cum_size;
    for(unsigned int d = 0; d < dimension; d++) center_of_mass[d] *= mult1;
    for(unsigned int d = 0; d < dimension; d++) center_of_mass[d] += mult2 * point[d];
    
    // If there is space in this quad tree and it is a leaf, add the object here
    if(is_leaf && size < QT_NODE_CAPACITY) {
        index[size] = new_index;
        size++;
        return true;
    }
    
    // Don't add duplicates for now (this is not very nice)
    bool any_duplicate = false;
    for(unsigned int n = 0; n < size; n++) {
        bool duplicate = true;
        for(unsigned int d = 0; d < dimension; d++) {
            if(point[d] != data[index[n] * dimension + d]) { duplicate = false; break; }
        }
        any_duplicate = any_duplicate | duplicate;
    }
    if(any_duplicate) return true;
    
    // Otherwise, we need to subdivide the current cell
    if(is_leaf) subdivide();
    
    // Find out where the point can be inserted
    for(SPTree* child : children) {
        if(child->insert(new_index)) return true;
    }
    
    // Otherwise, the point cannot be inserted (this should never happen)
    return false;
}

    
// Create four children which fully divide this cell into four quads of equal area
void SPTree::subdivide() {

    // Create new children
    std::vector<double> new_corner(dimension);
    std::vector<double> new_width(dimension);
    for(unsigned int i = 0; i < children.size(); i++) {
        unsigned int div = 1;
        for(unsigned int d = 0; d < dimension; d++) {
            new_width[d] = .5 * boundary.getWidth(d);
            if((i / div) % 2 == 1) new_corner[d] = boundary.getCorner(d) - new_width[d];
            else                   new_corner[d] = boundary.getCorner(d) + new_width[d];
            div *= 2;
        }
        children[i] = new SPTree(this, dimension, data, new_corner.data(), new_width.data());
    }

    // Move existing points to correct children
    for(unsigned int i = 0; i < size; i++) {
        bool success = false;
        for(SPTree* child : children) {
            if(!success) success = child->insert(index[i]);
        }
        index[i] = -1;
    }

    // Empty parent node
    size = 0;
    is_leaf = false;
}


// Build SPTree on dataset
void SPTree::fill(unsigned int N)
{
    for(unsigned int i = 0; i < N; i++) insert(i);
}


// Checks whether the specified tree is correct
bool SPTree::isCorrect()
{
    for(unsigned int n = 0; n < size; n++) {
        double* point = data + index[n] * dimension;
        if(!boundary.containsPoint(point)) return false;
    }
    if(!is_leaf) {
        bool correct = true;
        for(SPTree* child : children) correct = correct && child->isCorrect();
        return correct;
    }
    else return true;
}



// Build a list of all indices in SPTree
void SPTree::getAllIndices(unsigned int* indices)
{
    getAllIndices(indices, 0);
}


// Build a list of all indices in SPTree
unsigned int SPTree::getAllIndices(unsigned int* indices, unsigned int loc)
{
    
    // Gather indices in current quadrant
    for(unsigned int i = 0; i < size; i++) indices[loc + i] = index[i];
    loc += size;
    
    // Gather indices in children
    if(!is_leaf) {
        for(SPTree* child : children) loc = child->getAllIndices(indices, loc);
    }
    return loc;
}


unsigned int SPTree::getDepth() {
    if(is_leaf) return 1;
    unsigned int depth = 0;
    for(SPTree* child : children) depth = std::max(depth, child->getDepth());
    return 1 + depth;
}


// Compute non-edge forces using Barnes-Hut algorithm
void SPTree::computeNonEdgeForces(unsigned int point_index, double theta, double neg_f[], double* sum_Q)
{
    
    // Make sure that we spend no time on empty nodes or self-interactions
    if(cum_size == 0 || (is_leaf && size == 1 && index[0] == point_index)) return;
    
    // Compute distance between point and center-of-mass
    double D = .0;
    unsigned int ind = point_index * dimension;
    for(unsigned int d = 0; d < dimension; d++) buff[d] = data[ind + d] - center_of_mass[d];
    for(unsigned int d = 0; d < dimension; d++) D += buff[d] * buff[d];

    // Optimize (max_width / sqrt(D) < theta) by squaring and multiplying through by D
    if(is_leaf || max_width_sq < theta * theta * D) {
    
        // Compute and add t-SNE force between point and current node
        D = 1.0 / (1.0 + D);
        double mult = cum_size * D;
        *sum_Q += mult;
        mult *= D;
        for(unsigned int d = 0; d < dimension; d++) neg_f[d] += mult * buff[d];
    }
    else {

        // Recursively apply Barnes-Hut to children
        for(SPTree* child : children) child->computeNonEdgeForces(point_index, theta, neg_f, sum_Q);
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


// Print out tree
void SPTree::print() 
{
    if(cum_size == 0) {
        printf("Empty node\n");
        return;
    }

    if(is_leaf) {
        printf("Leaf node; data = [");
        for(int i = 0; i < size; i++) {
            double* point = data + index[i] * dimension;
            for(int d = 0; d < dimension; d++) printf("%f, ", point[d]);
            printf(" (index = %d)", index[i]);
            if(i < size - 1) printf("\n");
            else printf("]\n");
        }        
    }
    else {
        printf("Intersection node with center-of-mass = [");
        for(int d = 0; d < dimension; d++) printf("%f, ", center_of_mass[d]);
        printf("]; children are:\n");
        for(SPTree* child : children) child->print();
    }
}

