#include <SimpleBVH/BVH.hpp>
#include <catch2/catch_test_macros.hpp>

#include <iostream>

using namespace SimpleBVH;

TEST_CASE("test_tree", "[tests]")
{
    Eigen::MatrixXd pts(4, 3);
    // clang-format off
    pts << 0, 0, 0,
           3, 0, 0,
           0, 3, 0,
           0, 3, 3;
    // clang-format on

    Eigen::MatrixXi tri(4, 3);
    // clang-format off
    tri << 0, 1, 2,
           0, 1, 3,
           0, 2, 3,
           1, 2, 3;
    // clang-format on

    SimpleBVH::BVH bvh;
    bvh.init(pts, tri, 1e-10);

    std::vector<unsigned int> pairs;
    bvh.intersect_box(
        pts.colwise().minCoeff(), pts.colwise().maxCoeff(), pairs);

    CHECK(pairs.size() == tri.rows());
}

TEST_CASE("test_distance", "[tests]")
{
    Eigen::MatrixXd pts(4, 3);
    // clang-format off
    pts << 0, 0, 0,
           3, 0, 0,
           0, 3, 0,
           0, 3, 3;
    // clang-format on

    Eigen::MatrixXi edges(7, 2);
    // clang-format off
    edges << 0, 1,
           1, 2,
           2, 0,
           1, 3,
           3, 0,
           2, 3,
           1, 3;
    // clang-format on

    SimpleBVH::BVH bvh;
    VectorMax3d nearest_point;
    double sq_dist;
    bvh.init(pts, edges, 1e-10);
    bvh.set_leaf_distance_callback([&](const VectorMax3d& p, int f,
                                       VectorMax3d& cp, double& dist) {
        point_segment_squared_distance(
            p, { { pts.row(edges(f, 0)), pts.row(edges(f, 1)) } }, cp, dist);
    });
    bvh.set_get_point_callback([&](int f) { return pts.row(edges(f, 0)); });
    bvh.nearest_facet(Eigen::Vector3d(-1, -1, -1), nearest_point, sq_dist);

    CHECK(pts.row(0) == nearest_point.transpose());
}

TEST_CASE("test_distance_2d", "[tests]")
{
    Eigen::MatrixXd pts(4, 2);
    // clang-format off
    pts << 0, 0,
           3, 0,
           0, 3,
           0, 3;
    // clang-format on

    Eigen::MatrixXi edges(7, 2);
    // clang-format off
    edges << 0, 1,
           1, 2,
           2, 0,
           1, 3,
           3, 0,
           2, 3,
           1, 3;
    // clang-format on

    SimpleBVH::BVH bvh;
    VectorMax3d nearest_point;
    double sq_dist;
    bvh.init(pts, edges, 1e-10);
    bvh.set_leaf_distance_callback([&](const VectorMax3d& p, int f,
                                       VectorMax3d& cp, double& dist) {
        point_segment_squared_distance(
            p, { { pts.row(edges(f, 0)), pts.row(edges(f, 1)) } }, cp, dist);
    });
    bvh.set_get_point_callback([&](int f) { return pts.row(edges(f, 0)); });
    bvh.nearest_facet(Eigen::Vector2d(-1, -1), nearest_point, sq_dist);

    CHECK(pts.row(0) == nearest_point.transpose());
}
