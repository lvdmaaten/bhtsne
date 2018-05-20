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


#ifndef SPTREE_H
#define SPTREE_H

#include <deque>
#include <vector>

class SPTree
{
    // Internal node type
    struct Node;

    // The number of dimensions
    unsigned int dimension;

    // The actual points stored in the tree
    double* data;

    // The widths of the nodes at every level
    std::vector<std::vector<double>> widths;
    // The maximum width along any axis at the root level
    double max_width;

    // Holds all the child nodes (std::deque so that pointers are never invalidated)
    std::deque<Node> nodes;
    // Pointer to the root node
    Node* root;

    // A buffer we use when doing force computations
    std::vector<double> buff;

public:
    SPTree(unsigned int D, double* inp_data, unsigned int N);
    ~SPTree();

    void computeNonEdgeForces(unsigned int point_index, double theta, double neg_f[], double* sum_Q);
    void computeEdgeForces(unsigned int* row_P, unsigned int* col_P, double* val_P, int N, double* pos_f);
    void print();

private:
    Node* new_node(const double* point, std::vector<double> center, const double* width);
    void insert(Node* node, const double* point);
    Node* insertChild(Node* node, const double* point, unsigned int depth);
    void computeNonEdgeForces(Node* node, double max_width_sq, double* point, double theta_sq, double neg_f[], double* sum_Q);
    void print(Node* node);
};

#endif
