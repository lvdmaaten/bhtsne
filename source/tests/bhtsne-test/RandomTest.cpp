
#include <random>

#include <gmock/gmock.h>


class RandomTest : public testing::Test
{
public:
};

TEST_F(RandomTest, crossPlattformDeterministicRandom)
{
    auto expected = std::vector<unsigned long>{ 1125387415,2407456957,681542492,913057000,1194544295 };

    auto gen = std::mt19937(1337); //seed generator
    for (int i = 0; i < 5; ++i)
    {
        EXPECT_EQ(gen(), expected[i]);       
    }
}
