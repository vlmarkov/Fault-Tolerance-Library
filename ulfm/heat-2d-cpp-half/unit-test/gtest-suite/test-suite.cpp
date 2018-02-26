#include "../../Grid.cpp"
#include "../../Task.cpp"
#include <gtest/gtest.h>

/**
 * @brief This test checks mpi rank assingment for the each Task instance.
 */
TEST(MpiRankAssingmentTest, Equals)
{
    const auto x  = 8192;
    const auto y  = 8192;
    const auto px = 12;
    const auto py = 12;
    const auto nx = x / px;
    const auto ny = y / py;

    Grid grid(x, y, nx, ny, px, py);

    for (auto i = 0; i < (px * py); ++i)
    {
        Task* task = grid.getTask(i);

        EXPECT_TRUE(task);

        ASSERT_EQ(task->getMpiRank(0), i);

        if (i < ((px * py) / 2))
        {
            ASSERT_EQ(task->getMpiRank(1), i + ((px * py) / 2));
        }
        else
        {
            ASSERT_EQ(task->getMpiRank(1), i - ((px * py) / 2));
        }
    }
}

/**
 * @brief This test checks redundancy layer numbers for the each Task instance.
 */
TEST(RedundancyLayerNumbersTest, Equals)
{
    const auto x  = 8192;
    const auto y  = 8192;
    const auto px = 12;
    const auto py = 12;
    const auto nx = x / px;
    const auto ny = y / py;

    Grid grid(x, y, nx, ny, px, py);

    for (auto i = 0; i < (px * py); ++i)
    {
        Task* task = grid.getTask(i);

        EXPECT_TRUE(task);

        ASSERT_EQ(task->getLayersNumber(), 2);
    }
}

/**
 * @brief This test checks repair status for the each Task instance.
 */
TEST(RepairStatusTest, Equals)
{
    const auto x  = 8192;
    const auto y  = 8192;
    const auto px = 12;
    const auto py = 12;
    const auto nx = x / px;
    const auto ny = y / py;

    Grid grid(x, y, nx, ny, px, py);

    for (auto i = 0; i < (px * py); ++i)
    {
        Task* task = grid.getTask(i);

        EXPECT_TRUE(task);

        ASSERT_EQ(task->getRepairStatus(), 1);
    }
}

/**
 * @brief This test checks local grid and local new grid for the each Task instance.
 */
TEST(LocalGridSetGetTest, ExpectTrue)
{
    const auto x  = 8192;
    const auto y  = 8192;
    const auto px = 12;
    const auto py = 12;
    const auto nx = x / px;
    const auto ny = y / py;

    Grid grid(x, y, nx, ny, px, py);

    for (auto i = 0; i < (px * py); ++i)
    {
        Task* task = grid.getTask(i);

        EXPECT_TRUE(task);

        for (auto l = 0; l < task->getLayersNumber(); ++l)
        {
            task->allocateLocalGrid(l);
            task->allocateLocalNewGrid(l);

            EXPECT_TRUE(task->getLocalGrid(l));
            EXPECT_TRUE(task->getLocalNewGrid(l));
        }
    }
}

/**
 * @brief This test checks neighbors at the real layer for the each Task instance.
 */
TEST(CheckNeighborsTest, ExpectTrue)
{
    const auto x  = 8192;
    const auto y  = 8192;
    const auto px = 12;
    const auto py = 12;
    const auto nx = x / px;
    const auto ny = y / py;

    Grid grid(x, y, nx, ny, px, py);

    for (auto i = 0; i < (px * py); ++i)
    {
        Task* task = grid.getTask(i);

        EXPECT_TRUE(task);

        if (task->getLeftNeighbor())
        {
            ASSERT_EQ(task->getLeftNeighbor()->getMpiRank(0),
                      task->getLeftNeighborRank(0));
        }

        if (task->getRightNeighbor())
        {
            ASSERT_EQ(task->getRightNeighbor()->getMpiRank(0),
                      task->getRightNeighborRank(0));
        }

        if (task->getUpNeighbor())
        {
            ASSERT_EQ(task->getUpNeighbor()->getMpiRank(0),
                      task->getUpNeighborRank(0));
        }


        if (task->getDownNeighbor())
        {
            ASSERT_EQ(task->getDownNeighbor()->getMpiRank(0),
                      task->getDownNeighborRank(0));
        }
    }
}

/**
 * @brief This test checks message tags for the each Task instance.
 */
TEST(CheckTagsTest, ExpectTrue)
{
    const auto x  = 8192;
    const auto y  = 8192;
    const auto px = 12;
    const auto py = 12;
    const auto nx = x / px;
    const auto ny = y / py;

    Grid grid(x, y, nx, ny, px, py);

    for (auto i = 0; i < (px * py); ++i)
    {
        Task* task = grid.getTask(i);

        EXPECT_TRUE(task);

        for (auto l = 0; l < task->getLayersNumber(); l++)
        {
            if (task->getLeftNeighbor())
            {
                ASSERT_EQ(task->getLeftNeighbor()->getRightTag(l),
                          task->getLeftTag(l));
            }

            if (task->getRightNeighbor())
            {
                ASSERT_EQ(task->getRightNeighbor()->getLeftTag(l),
                          task->getRightTag(l));
            }

            if (task->getUpNeighbor())
            {
                ASSERT_EQ(task->getUpNeighbor()->getDownTag(l),
                          task->getUpTag(l));
            }

            if (task->getDownNeighbor())
            {
                ASSERT_EQ(task->getDownNeighbor()->getUpTag(l),
                          task->getDownTag(l));
            }
        }
    }
}

/**
 * @brief This test checks redundancy tasks for the each Task instance.
 */
TEST(CheckRedundancyTasksTest, ExpectTrue)
{
    const auto x  = 8192;
    const auto y  = 8192;
    const auto px = 12;
    const auto py = 12;
    const auto nx = x / px;
    const auto ny = y / py;

    Grid grid(x, y, nx, ny, px, py);

    for (auto i = 0; i < (px * py); ++i)
    {
        Task* task = grid.getTask(i);

        EXPECT_TRUE(task);

        auto replacements = task->getRedundancyTasks();

        ASSERT_EQ(replacements.size(), 2);

        ASSERT_EQ(replacements[0]->getMpiRank(0), task->getMpiRank(0));
        ASSERT_EQ(replacements[0]->getMpiRank(1), task->getMpiRank(1));
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
