#include <gmock/gmock.h>
#include "../../bhtsne/source/SpacePartitioningTreeTemplate.h"

using namespace bhtsne;

class SpacePartitioningTreeTest : public testing::Test
{
public:
};

TEST_F(SpacePartitioningTreeTest, OpenMPComputeNonEdgeForces)
{
    const unsigned int dataSize = 10;
    const unsigned int outputDimensions = 2;
    const double gradientAccuracy = 2.0;
    const auto result = Vector2D<double>{{
        { 0,56,19,80,58 },
        { 47,35,89,82,74 },
        { 17,85,71,51,30 },
        { 1,9,36,14,16 },
        { 98,44,11,0,0 },
        { 37,53,57,60,60 },
        { 16,66,45,35,5 },
        { 60,78,80,51,30 },
        { 87,72,95,92,53 },
        { 14,46,23,86,20 }
    }};

    auto tree = SpacePartitioningTree<3>(result);

    auto neg_f = Vector2D<double>(dataSize, outputDimensions, 0.0);
    double sum_Q = 0.0;
    for (unsigned int n = 0; n < dataSize; ++n)
    {
        tree.computeNonEdgeForces(n, gradientAccuracy, neg_f[n], sum_Q);
    }

    auto omp_neg_f = Vector2D<double>(dataSize, outputDimensions, 0.0);
    double omp_sum_Q = 0.0;
    // omp version on windows (2.0) does only support signed loop variables, should be unsigned
    #pragma omp parallel for reduction(+:omp_sum_Q)
    for (int n = 0; n < static_cast<int>(dataSize); ++n)
    {
        tree.computeNonEdgeForces(n, gradientAccuracy, omp_neg_f[n], omp_sum_Q);
    }

    ASSERT_DOUBLE_EQ(sum_Q, omp_sum_Q);

    for (unsigned int i = 0; i < neg_f.height(); ++i)
    {
        for (unsigned int j = 0; j < neg_f.width(); ++j)
        {
            ASSERT_FLOAT_EQ(neg_f[i][j], omp_neg_f[i][j]);
        }
    }

    // TODO this still sometimes (5% of the time max) fails...
}
