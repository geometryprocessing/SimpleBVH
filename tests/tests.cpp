#include <SimpleBVH/BVH.hpp>
#include <catch2/catch_test_macros.hpp>

#include <iostream>
#include <fstream>

using namespace SimpleBVH;
using namespace Eigen;

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
    Vector3d nearest_point;
    double sq_dist;
    bvh.init(pts, edges, 1e-10);
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
    Vector2d nearest_point;
    double sq_dist;
    bvh.init(pts, edges, 1e-10);
    bvh.nearest_facet(Vector2d(-1, -1), nearest_point, sq_dist);

    CHECK(pts.row(0) == nearest_point.transpose());
}

TEST_CASE("test_distance_3d", "[tests]")
{
    MatrixXd vertices;
    MatrixXi facets;
    std::string root = BVH_ROOT_PATH;

    std::string mesh_filename = root + "/bunny.off";
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i) {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i) {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    Vector3d min = vertices.colwise().minCoeff();
    Vector3d max = vertices.colwise().maxCoeff();

    MatrixXd pts = MatrixXd::Random(100, 3);

    for (int d = 0; d < 3; ++d) {
        pts.col(d).array() += 1;
        pts.col(d) *= (max(d) - min(d)) / 2;
        pts.col(d).array() += min(d);
    }

    //     MatrixXd pts(1, 3);
    //     pts << 0.2, -0.2, 0.4;

    SimpleBVH::BVH bvh;
    bvh.init(vertices, facets, 1e-10);

    std::ofstream out("res.obj");

    Vector3d nearest_point;
    double sq_dist;
    for (int i = 0; i < pts.rows(); ++i) {
        bvh.nearest_facet(pts.row(i), nearest_point, sq_dist);

        out << "v " << pts(i, 0) << " " << pts(i, 1) << " " << pts(i, 2)
            << "\n";
        out << "v " << nearest_point(0) << " " << nearest_point(1) << " "
            << nearest_point(2) << "\n";
        out << "l " << 2 * i + 1 << " " << 2 * i + 2 << "\n\n";
    }
    //     std::ofstream out1("res1.obj");
    //     for (int i = 0; i < vertices.rows(); ++i) {
    //         out1 << "v " << vertices(i, 0) << " " << vertices(i, 1) << " "
    //              << vertices(i, 2) << "\n";
    //     }
    //     out1 << "f " << facets(66, 0) + 1 << " " << facets(66, 1) + 1 << " "
    //          << facets(66, 2) + 1 << "\n";
}
