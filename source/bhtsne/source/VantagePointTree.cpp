#include "VantagePointTree.h"

#include <cassert>

//using namespace bhtsne;

DataPoint::DataPoint(unsigned int dimensions, unsigned int index, double * x)
	: dimensions(dimensions)
	, index(index)
	, data(x, x + dimensions)
{}

double VantagePointTree::euclideanDistance(const DataPoint & a, const DataPoint & b)
{
    assert(a.dimensions == b.dimensions);
    double squaredDistance = 0.0;
    for (size_t d = 0; d < a.dimensions; ++d)
    {
        double difference = a.data[d] - b.data[d];
        squaredDistance += difference * difference;
    }
    return sqrt(squaredDistance);
}

VantagePointTree::VantagePointTree(DistanceFunction distanceFunction)
        : m_maxDistance(0.0)
        , m_randomNumberGenerator(std::random_device()())
        , m_distanceFunction(std::move(distanceFunction))
        , m_root(std::make_unique<Node>())
{}

void VantagePointTree::create(const std::vector<DataPoint> & items)
{
    m_items = items;
    m_root = buildFromPoints(0, static_cast<unsigned int>(items.size()));
}

void VantagePointTree::search(const DataPoint & target, unsigned int k, std::vector<DataPoint> & results,
                              std::vector<double> & distances)
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

std::unique_ptr<VantagePointTree::Node> VantagePointTree::buildFromPoints(unsigned int lower, unsigned int upper)
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
                         [this, lower](const DataPoint & a, const DataPoint & b){
                             return m_distanceFunction(m_items[lower], a) < m_distanceFunction(m_items[lower], b); });

        // Threshold of the new node will be the distance to the median
        node->threshold = m_distanceFunction(m_items[lower], m_items[median]);

        // Recursively build tree
        node->leftChild = buildFromPoints(lower + 1, median);
        node->rightChild = buildFromPoints(median, upper);
    }

    return node;
}

void VantagePointTree::search(const VantagePointTree::Node & node, const DataPoint & target, unsigned int k,
                              std::priority_queue<VantagePointTree::HeapItem> & heap)
{
    // Compute distance between target and current node
    double distance = m_distanceFunction(m_items[node.index], target);

    // If current node within radius tau
    if(distance < m_maxDistance)
    {
        if(heap.size() == k)
        {
            // remove furthest node from result list (if we already have k results)
            heap.pop();
        }
        heap.push(HeapItem{node.index, distance}); // add current node to result list
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

VantagePointTree::Node::Node(unsigned int index)
        : index(index)
        , threshold(0.0)
{}

bool VantagePointTree::HeapItem::operator<(const VantagePointTree::HeapItem &other) const {
    return distance < other.distance;
}
