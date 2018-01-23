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


/* This code was adopted with minor modifications from Steve Hanov's great tutorial at http://stevehanov.ca/blog/index.php?id=130 */

#pragma once

#include <algorithm>
#include <vector>
#include <queue>
#include <memory>
#include <random>
#include <functional>

struct DataPoint
{
    unsigned int dimensions;
    unsigned int index;
    std::vector<double> data;

    DataPoint(unsigned int dimensions, unsigned int index, double * x);
};

using DistanceFunction = std::function<double(const DataPoint &, const DataPoint &)>;

class VantagePointTree
{
public:
    explicit VantagePointTree(DistanceFunction distanceFunction = euclideanDistance);

    // possible distance functions
    static double euclideanDistance(const DataPoint & a, const DataPoint & b);
    //TODO create some more common distance functions

    // Function to create a new VantagePointTree from data
    void create(const std::vector<DataPoint> & items);

    // Function that uses the tree to find the k nearest neighbors of target
    void search(const DataPoint & target, unsigned int k, std::vector<DataPoint> & results,
                std::vector<double> & distances);

private:
    std::vector<DataPoint> m_items;
    double m_maxDistance;
    std::mt19937 m_randomNumberGenerator;
    DistanceFunction m_distanceFunction;


    // Single node of a VP tree (has a point and radius; left children are closer to point than the radius)
    struct Node
    {
        unsigned int index; // index of point in node
        double threshold; // radius(?)
        std::unique_ptr<Node> leftChild; // points closer by than threshold
        std::unique_ptr<Node> rightChild; // points farther away than threshold

        explicit Node(unsigned int index = 0);
    };

    std::unique_ptr<Node> m_root;

    // An item on the intermediate result queue
    struct HeapItem {
        unsigned int index;
        double distance;

        bool operator<(const HeapItem & other) const;
    };

    // Function that (recursively) fills the tree
    std::unique_ptr<Node> buildFromPoints(unsigned int lower, unsigned int upper);

    // Helper function that searches the tree
    void search(const Node & node, const DataPoint & target, unsigned int k, std::priority_queue<HeapItem> & heap);
};
