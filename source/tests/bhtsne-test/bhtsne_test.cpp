
#include <gmock/gmock.h>



class bhtsne_test: public testing::Test
{
public:
};

TEST_F(bhtsne_test, sanity_checks)
{
    EXPECT_EQ((unsigned int) 0, 0);
    EXPECT_EQ((unsigned int) 1, 1);
}
