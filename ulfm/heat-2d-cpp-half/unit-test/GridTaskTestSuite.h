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

    auto* task_0 = grid.getTask(0);
    auto* task_1 = grid.getTask(1);
    auto* task_2 = grid.getTask(2);
    auto* task_3 = grid.getTask(3);

    TS_ASSERT(task_0);
    TS_ASSERT(task_1);
    TS_ASSERT(task_2);
    TS_ASSERT(task_3);

    // Checks 0 rank
    TS_ASSERT_EQUALS(task_0->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_0->getMpiRank(0),           0);
    TS_ASSERT_EQUALS(task_0->getMpiRank(1),           2);

    TS_ASSERT_EQUALS(task_0->getUpTag(0),             0);
    TS_ASSERT_EQUALS(task_0->getUpTag(1),             0);

    TS_ASSERT_EQUALS(task_0->getDownTag(0),           2);
    TS_ASSERT_EQUALS(task_0->getDownTag(1),           6);

    TS_ASSERT_EQUALS(task_0->getLeftTag(0),           0);
    TS_ASSERT_EQUALS(task_0->getLeftTag(1),           0);

    TS_ASSERT_EQUALS(task_0->getRightTag(0),          1);
    TS_ASSERT_EQUALS(task_0->getRightTag(1),          5);

    TS_ASSERT_EQUALS(task_0->getUpNeighborRank(0),   -1);
    TS_ASSERT_EQUALS(task_0->getUpNeighborRank(1),    2);

    TS_ASSERT_EQUALS(task_0->getDownNeighborRank(0),  2);
    TS_ASSERT_EQUALS(task_0->getDownNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_0->getLeftNeighborRank(0), -1);
    TS_ASSERT_EQUALS(task_0->getLeftNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_0->getRightNeighborRank(0), 1);
    TS_ASSERT_EQUALS(task_0->getRightNeighborRank(1), 1);

    // Checks 1 rank
    TS_ASSERT_EQUALS(task_1->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_1->getMpiRank(0),           1);
    TS_ASSERT_EQUALS(task_1->getMpiRank(1),           3);

    TS_ASSERT_EQUALS(task_1->getUpTag(0),             0);
    TS_ASSERT_EQUALS(task_1->getUpTag(1),             0);

    TS_ASSERT_EQUALS(task_1->getDownTag(0),           3);
    TS_ASSERT_EQUALS(task_1->getDownTag(1),           7);

    TS_ASSERT_EQUALS(task_1->getLeftTag(0),           1);
    TS_ASSERT_EQUALS(task_1->getLeftTag(1),           5);

    TS_ASSERT_EQUALS(task_1->getRightTag(0),          0);
    TS_ASSERT_EQUALS(task_1->getRightTag(1),          0);

    TS_ASSERT_EQUALS(task_1->getUpNeighborRank(0),   -1);
    TS_ASSERT_EQUALS(task_1->getUpNeighborRank(1),    3);

    TS_ASSERT_EQUALS(task_1->getDownNeighborRank(0),  3);
    TS_ASSERT_EQUALS(task_1->getDownNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_1->getLeftNeighborRank(0),  0);
    TS_ASSERT_EQUALS(task_1->getLeftNeighborRank(1),  0);

    TS_ASSERT_EQUALS(task_1->getRightNeighborRank(0),-1);
    TS_ASSERT_EQUALS(task_1->getRightNeighborRank(1),-1);

    // Checks 2 rank
    TS_ASSERT_EQUALS(task_2->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_2->getMpiRank(0),           2);
    TS_ASSERT_EQUALS(task_2->getMpiRank(1),           0);

    TS_ASSERT_EQUALS(task_2->getUpTag(0),             2);
    TS_ASSERT_EQUALS(task_2->getUpTag(1),             6);

    TS_ASSERT_EQUALS(task_2->getDownTag(0),           0);
    TS_ASSERT_EQUALS(task_2->getDownTag(1),           0);

    TS_ASSERT_EQUALS(task_2->getLeftTag(0),           0);
    TS_ASSERT_EQUALS(task_2->getLeftTag(1),           0);

    TS_ASSERT_EQUALS(task_2->getRightTag(0),          4);
    TS_ASSERT_EQUALS(task_2->getRightTag(1),          8);

    TS_ASSERT_EQUALS(task_2->getUpNeighborRank(0),    0);
    TS_ASSERT_EQUALS(task_2->getUpNeighborRank(1),   -1);

    TS_ASSERT_EQUALS(task_2->getDownNeighborRank(0), -1);
    TS_ASSERT_EQUALS(task_2->getDownNeighborRank(1),  0);

    TS_ASSERT_EQUALS(task_2->getLeftNeighborRank(0), -1);
    TS_ASSERT_EQUALS(task_2->getLeftNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_2->getRightNeighborRank(0), 3);
    TS_ASSERT_EQUALS(task_2->getRightNeighborRank(1), 3);

    // Checks 3 rank
    TS_ASSERT_EQUALS(task_3->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_3->getMpiRank(0),           3);
    TS_ASSERT_EQUALS(task_3->getMpiRank(1),           1);

    TS_ASSERT_EQUALS(task_3->getUpTag(0),             3);
    TS_ASSERT_EQUALS(task_3->getUpTag(1),             7);

    TS_ASSERT_EQUALS(task_3->getDownTag(0),           0);
    TS_ASSERT_EQUALS(task_3->getDownTag(1),           0);

    TS_ASSERT_EQUALS(task_3->getLeftTag(0),           4);
    TS_ASSERT_EQUALS(task_3->getLeftTag(1),           8);

    TS_ASSERT_EQUALS(task_3->getRightTag(0),          0);
    TS_ASSERT_EQUALS(task_3->getRightTag(1),          0);

    TS_ASSERT_EQUALS(task_3->getUpNeighborRank(0),    1);
    TS_ASSERT_EQUALS(task_3->getUpNeighborRank(1),   -1);

    TS_ASSERT_EQUALS(task_3->getDownNeighborRank(0), -1);
    TS_ASSERT_EQUALS(task_3->getDownNeighborRank(1),  1);

    TS_ASSERT_EQUALS(task_3->getLeftNeighborRank(0),  2);
    TS_ASSERT_EQUALS(task_3->getLeftNeighborRank(1),  2);

    TS_ASSERT_EQUALS(task_3->getRightNeighborRank(0),-1);
    TS_ASSERT_EQUALS(task_3->getRightNeighborRank(1),-1);

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

    auto* task_2 = grid.getTask(2);
    auto* task_3 = grid.getTask(3);

    TS_ASSERT(task_2);
    TS_ASSERT(task_3);

    // Checks 2 rank
    TS_ASSERT_EQUALS(task_2->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_2->getMpiRank(0),           2);
    TS_ASSERT_EQUALS(task_2->getMpiRank(1),           2);

    TS_ASSERT_EQUALS(task_2->getUpTag(0),             0);
    TS_ASSERT_EQUALS(task_2->getUpTag(1),             0);

    TS_ASSERT_EQUALS(task_2->getDownTag(0),           2);
    TS_ASSERT_EQUALS(task_2->getDownTag(1),           6);

    TS_ASSERT_EQUALS(task_2->getLeftTag(0),           0);
    TS_ASSERT_EQUALS(task_2->getLeftTag(1),           0);

    TS_ASSERT_EQUALS(task_2->getRightTag(0),          1);
    TS_ASSERT_EQUALS(task_2->getRightTag(1),          5);

    TS_ASSERT_EQUALS(task_2->getUpNeighborRank(0),   -1);
    TS_ASSERT_EQUALS(task_2->getUpNeighborRank(1),    2);

    TS_ASSERT_EQUALS(task_2->getDownNeighborRank(0),  2);
    TS_ASSERT_EQUALS(task_2->getDownNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_2->getLeftNeighborRank(0), -1);
    TS_ASSERT_EQUALS(task_2->getLeftNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_2->getRightNeighborRank(0), 3);
    TS_ASSERT_EQUALS(task_2->getRightNeighborRank(1), 3);

    // Checks 3 rank
    TS_ASSERT_EQUALS(task_3->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_3->getMpiRank(0),           3);
    TS_ASSERT_EQUALS(task_3->getMpiRank(1),           3);

    TS_ASSERT_EQUALS(task_3->getUpTag(0),             0);
    TS_ASSERT_EQUALS(task_3->getUpTag(1),             0);

    TS_ASSERT_EQUALS(task_3->getDownTag(0),           3);
    TS_ASSERT_EQUALS(task_3->getDownTag(1),           7);

    TS_ASSERT_EQUALS(task_3->getLeftTag(0),           1);
    TS_ASSERT_EQUALS(task_3->getLeftTag(1),           5);

    TS_ASSERT_EQUALS(task_3->getRightTag(0),          0);
    TS_ASSERT_EQUALS(task_3->getRightTag(1),          0);

    TS_ASSERT_EQUALS(task_3->getUpNeighborRank(0),   -1);
    TS_ASSERT_EQUALS(task_3->getUpNeighborRank(1),    3);

    TS_ASSERT_EQUALS(task_3->getDownNeighborRank(0),  3);
    TS_ASSERT_EQUALS(task_3->getDownNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_3->getLeftNeighborRank(0),  2);
    TS_ASSERT_EQUALS(task_3->getLeftNeighborRank(1),  2);

    TS_ASSERT_EQUALS(task_3->getRightNeighborRank(0),-1);
    TS_ASSERT_EQUALS(task_3->getRightNeighborRank(1),-1);

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

    for (auto i = (px * py); i > px; --i)
    {
        TS_ASSERT_EQUALS(grid.kill(i - 1), 0);
        TS_ASSERT_EQUALS(grid.repair(), 0);
    }

    auto* task_0 = grid.getTask(0);
    auto* task_1 = grid.getTask(1);

    TS_ASSERT(task_0);
    TS_ASSERT(task_1);

    // Checks 0 rank
    TS_ASSERT_EQUALS(task_0->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_0->getMpiRank(0),           0);
    TS_ASSERT_EQUALS(task_0->getMpiRank(1),           0);

    TS_ASSERT_EQUALS(task_0->getUpTag(0),             0);
    TS_ASSERT_EQUALS(task_0->getUpTag(1),             0);

    TS_ASSERT_EQUALS(task_0->getDownTag(0),           2);
    TS_ASSERT_EQUALS(task_0->getDownTag(1),           6);

    TS_ASSERT_EQUALS(task_0->getLeftTag(0),           0);
    TS_ASSERT_EQUALS(task_0->getLeftTag(1),           0);

    TS_ASSERT_EQUALS(task_0->getRightTag(0),          1);
    TS_ASSERT_EQUALS(task_0->getRightTag(1),          5);

    TS_ASSERT_EQUALS(task_0->getUpNeighborRank(0),   -1);
    TS_ASSERT_EQUALS(task_0->getUpNeighborRank(1),    0);

    TS_ASSERT_EQUALS(task_0->getDownNeighborRank(0),  0);
    TS_ASSERT_EQUALS(task_0->getDownNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_0->getLeftNeighborRank(0), -1);
    TS_ASSERT_EQUALS(task_0->getLeftNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_0->getRightNeighborRank(0), 1);
    TS_ASSERT_EQUALS(task_0->getRightNeighborRank(1), 1);

    // Checks 1 rank
    TS_ASSERT_EQUALS(task_1->getLayersNumber(),       2);

    TS_ASSERT_EQUALS(task_1->getMpiRank(0),           1);
    TS_ASSERT_EQUALS(task_1->getMpiRank(1),           1);

    TS_ASSERT_EQUALS(task_1->getUpTag(0),             0);
    TS_ASSERT_EQUALS(task_1->getUpTag(1),             0);

    TS_ASSERT_EQUALS(task_1->getDownTag(0),           3);
    TS_ASSERT_EQUALS(task_1->getDownTag(1),           7);

    TS_ASSERT_EQUALS(task_1->getLeftTag(0),           1);
    TS_ASSERT_EQUALS(task_1->getLeftTag(1),           5);

    TS_ASSERT_EQUALS(task_1->getRightTag(0),          0);
    TS_ASSERT_EQUALS(task_1->getRightTag(1),          0);

    TS_ASSERT_EQUALS(task_1->getUpNeighborRank(0),   -1);
    TS_ASSERT_EQUALS(task_1->getUpNeighborRank(1),    1);

    TS_ASSERT_EQUALS(task_1->getDownNeighborRank(0),  1);
    TS_ASSERT_EQUALS(task_1->getDownNeighborRank(1), -1);

    TS_ASSERT_EQUALS(task_1->getLeftNeighborRank(0),  0);
    TS_ASSERT_EQUALS(task_1->getLeftNeighborRank(1),  0);

    TS_ASSERT_EQUALS(task_1->getRightNeighborRank(0),-1);
    TS_ASSERT_EQUALS(task_1->getRightNeighborRank(1),-1);

    TS_TRACE("Test  Grid Kill&Reapir (Backward): done");
}