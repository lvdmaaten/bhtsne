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


struct DataPoint
{
    unsigned int dimensions;
    unsigned int index;
    std::vector<double> data;

    DataPoint(unsigned int dimensions, unsigned int index, double * x)
        : dimensions(dimensions)
        , index(index)
        , data(x, x + dimensions)
    {}
};

double euclideanDistance(const DataPoint & a, const DataPoint & b) {
    assert(a.dimensions == b.dimensions);
    double squaredDistance = 0.0;
    for (size_t d = 0; d < a.dimensions; ++d)
    {
        double difference = a.data[d] - b.data[d];
        squaredDistance += difference * difference;
    }
    return sqrt(squaredDistance);
}


template<typename T, double calcDistance(const T &, const T &)>
class VantagePointTree
{
public:
    // Default constructor
    VantagePointTree()
        : m_maxDistance(0.0)
        , m_randomNumberGenerator(std::random_device()())
        , m_root(std::make_unique<Node>())
    {}

    // Function to create a new VantagePointTree from data
    void create(const std::vector<T> & items) {
        m_items = items;
        m_root = buildFromPoints(0, static_cast<unsigned int>(items.size()));
    }

    // Function that uses the tree to find the k nearest neighbors of target
    void search(const T & target, unsigned int k, std::vector<T> & results, std::vector<double> & distances)
    {

        // Use a priority queue to store intermediate results on
        std::priority_queue<HeapItem> heap;

        m_maxDistance = std::numeric_limits<double>::max();

        // Perform the search
        search(*m_root, target, k, heap);

        // Gather final results
        results.clear();
        distances.clear();
        while(!heap.empty())
        {
            results.push_back(m_items[heap.top().index]);
            distances.push_back(heap.top().distance);
            heap.pop();
        }

        // Results are in reverse order TODO is this really necessary?
        std::reverse(results.begin(), results.end());
        std::reverse(distances.begin(), distances.end());
    }

private:
    std::vector<T> m_items;
    double m_maxDistance;
    std::mt19937 m_randomNumberGenerator;


    // Single node of a VP tree (has a point and radius; left children are closer to point than the radius)
    struct Node
    {
        unsigned int index; // index of point in node
        double threshold; // radius(?)
        std::unique_ptr<Node> leftChild; // points closer by than threshold
        std::unique_ptr<Node> rightChild; // points farther away than threshold

        explicit Node(unsigned int index = 0)
            : index(index)
            , threshold(0.0)
        {}
    };

    std::unique_ptr<Node> m_root;


    // An item on the intermediate result queue
    struct HeapItem {
        HeapItem(unsigned int index, double distance)
            : index(index)
            , distance(distance)
        {}

        unsigned int index;
        double distance;

        bool operator<(const HeapItem & other) const
        {
            return distance < other.distance;
        }
    };

    // needed for std::nth_element
    struct DistanceComparator
    {
        const T & item;
        explicit DistanceComparator(const T & item)
            : item(item)
        {}
        bool operator()(const T & a, const T & b) {
            return calcDistance(item, a) < calcDistance(item, b);
        }
    };

    // Function that (recursively) fills the tree
    std::unique_ptr<Node> buildFromPoints(unsigned int lower, unsigned int upper)
    {
        if (upper == lower)
        {
            // we're done here, return null pointer (default constructed unique_ptr)
            return std::unique_ptr<Node>();
        }

        // Lower index is center of current node
        auto node = std::make_unique<Node>(lower);

        if (upper - lower > 1)
        {      // if we did not arrive at leaf yet

            // Choose an arbitrary point and move it to the start
            auto uniformDistribution = std::uniform_int_distribution<int>(lower, upper - 1);
            auto random_index = uniformDistribution(m_randomNumberGenerator);
            std::swap(m_items[lower], m_items[random_index]);

            // Partition around the median distance
            unsigned int median = (lower + upper) / 2;
            std::nth_element(m_items.begin() + lower + 1,
                             m_items.begin() + median,
                             m_items.begin() + upper,
                             DistanceComparator(m_items[lower]));

            // Threshold of the new node will be the distance to the median
            node->threshold = calcDistance(m_items[lower], m_items[median]);

            // Recursively build tree
            node->leftChild = buildFromPoints(lower + 1, median);
            node->rightChild = buildFromPoints(median, upper);
        }

        return node;
    }

    // Helper function that searches the tree
    void search(const Node & node, const T & target, unsigned int k, std::priority_queue<HeapItem> & heap)
    {
        // Compute distance between target and current node
        double distance = calcDistance(m_items[node.index], target);

        // If current node within radius tau
        if(distance < m_maxDistance)
        {
            if(heap.size() == k)
            {
                // remove furthest node from result list (if we already have k results)
                heap.pop();
            }
            heap.push(HeapItem(node.index, distance)); // add current node to result list
            if(heap.size() == k)
            {
                // update value of tau (farthest point in result list)
                m_maxDistance = heap.top().distance;
            }
        }

        // If the target lies within the radius of ball
        if(distance < node.threshold)
        {
            // if there can still be neighbors inside the ball, recursively search left child first
            if(distance - m_maxDistance <= node.threshold && node.leftChild)
            {
                search(*node.leftChild, target, k, heap);
            }

            // if there can still be neighbors outside the ball, recursively search right child
            if(distance + m_maxDistance >= node.threshold && node.rightChild)
            {
                search(*node.rightChild, target, k, heap);
            }

        // If the target lies outsize the radius of the ball
        }
        else
        {
            // if there can still be neighbors outside the ball, recursively search right child first
            if(distance + m_maxDistance >= node.threshold && node.rightChild)
            {
                search(*node.rightChild, target, k, heap);
            }

            // if there can still be neighbors inside the ball, recursively search left child
            if (distance - m_maxDistance <= node.threshold && node.leftChild)
            {
                search(*node.leftChild, target, k, heap);
            }
        }
    }
};
