////////////////////////////////////////////////////////////////////////////////
#include <BVH.cpp>
#include <Eigen/Dense>
#include <catch.hpp>
////////////////////////////////////////////////////////////////////////////////

TEST_CASE("test_tree", "[tests]")
{
    Eigen::MatrixXd pts(4, 3);
    pts << 0, 0, 0,
        3, 0, 0,
        0, 3, 0,
        0, 3, 3;

    Eigen::MatrixXi tri(4, 3);
    tri << 0, 1, 2,
        0, 1, 3,
        0, 2, 3,
        1, 2, 3;

    BVH::BVH bvh;
    bvh.init(pts, tri, 1e-10);

    std::vector<unsigned int> pairs;
    bvh.intersect_box(pts.colwise().minCoeff(), pts.colwise().maxCoeff(), pairs);

    CHECK(pairs.size() == tri.rows());
}
