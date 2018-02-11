#include <cxxtest/TestSuite.h>

#include "../Grid.h"
#include "../Task.h"

class GridTaskTestSuite : public CxxTest::TestSuite
{
    public:
        void testGridCreation(void);
        void testGridKillRepairForward(void);
        void testGridKillRepairBackward(void);
};

void GridTaskTestSuite::testGridCreation(void)
{
    int nx  = 1024;
    int ny  = 1024;
    int npx = 256;
    int npy = 256;
    int px  = 2;
    int py  = 2;

    TS_TRACE("Test Grid creation: start");
/*
    std::cout << "Creating arguments: " << std::endl
              << nx << std::endl
              << ny << std::endl
              << npx << std::endl
              << npy << std::endl
              << px << std::endl
              << py << std::endl;
*/
    Grid grid(nx, ny, npx, npy, px, py);

    auto* task = grid.getTask(0);

    TS_ASSERT_EQUALS(task->getLayersNumber(), 2);

    TS_ASSERT_EQUALS(task->getUpTag(0),    0);
    TS_ASSERT_EQUALS(task->getLeftTag(0),  0);
    TS_ASSERT_EQUALS(task->getDownTag(0),  2);
    TS_ASSERT_EQUALS(task->getRightTag(0), 1);

    TS_ASSERT_EQUALS(task->getUpTag(1),    0);
    TS_ASSERT_EQUALS(task->getLeftTag(1),  0);
    TS_ASSERT_EQUALS(task->getDownTag(1),  6);
    TS_ASSERT_EQUALS(task->getRightTag(1), 5);

    TS_ASSERT_EQUALS(task->getUpNeighborRank(0),   -1);
    TS_ASSERT_EQUALS(task->getLeftNeighborRank(0), -1);
    TS_ASSERT_EQUALS(task->getDownNeighborRank(0),  2);
    TS_ASSERT_EQUALS(task->getRightNeighborRank(0), 1);

    TS_ASSERT_EQUALS(task->getUpNeighborRank(1),    2);
    TS_ASSERT_EQUALS(task->getLeftNeighborRank(1), -1);
    TS_ASSERT_EQUALS(task->getDownNeighborRank(1), -1);
    TS_ASSERT_EQUALS(task->getRightNeighborRank(1), 1);    

    TS_TRACE("Test Grid creation: done");
}

void GridTaskTestSuite::testGridKillRepairForward(void)
{
    int nx  = 1024;
    int ny  = 1024;
    int npx = 256;
    int npy = 256;
    int px  = 2;
    int py  = 2;

    TS_TRACE("Test Grid Kill&Reapir (Forward): start");
/*
    std::cout << "Creating arguments: " << std::endl
              << nx  << std::endl
              << ny  << std::endl
              << npx << std::endl
              << npy << std::endl
              << px  << std::endl
              << py  << std::endl;
*/

    Grid grid(nx, ny, npx, npy, px, py);

    for (auto i = 0; i < px; ++i)
    {
        TS_ASSERT_EQUALS(grid.kill(i), 0);
        TS_ASSERT_EQUALS(grid.repair(), 0);
    }

    TS_TRACE("Test  Grid Kill&Reapir (Forward): done");
}

void GridTaskTestSuite::testGridKillRepairBackward(void)
{
    int nx  = 1024;
    int ny  = 1024;
    int npx = 256;
    int npy = 256;
    int px  = 2;
    int py  = 2;

    TS_TRACE("Test Grid Kill&Reapir (Backward): start");
/*
    std::cout << "Creating arguments: " << std::endl
              << nx  << std::endl
              << ny  << std::endl
              << npx << std::endl
              << npy << std::endl
              << px  << std::endl
              << py  << std::endl;
*/

    Grid grid(nx, ny, npx, npy, px, py);

    for (auto i = (px * py) - 1; i > px; --i)
    {
        TS_ASSERT_EQUALS(grid.kill(i), 0);
        TS_ASSERT_EQUALS(grid.repair(), 0);
    }

    TS_TRACE("Test  Grid Kill&Reapir (Backward): done");
}