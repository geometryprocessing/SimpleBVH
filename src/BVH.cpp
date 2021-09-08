#include <BVH.hpp>
#include <Morton.hpp>

#include <iostream>

namespace BVH {
namespace {
    bool box_box_intersection(
        const Eigen::Vector3d& min1,
        const Eigen::Vector3d& max1,
        const Eigen::Vector3d& min2,
        const Eigen::Vector3d& max2)
    {
        if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
            return 0;
        if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
            return 0;
        return 1;
    }
} // namespace

void BVH::init_boxes_recursive(
    const std::vector<std::array<Eigen::Vector3d, 2>>& cornerlist,
    int node_index,
    int b,
    int e)
{
    assert(b != e);
    assert(node_index < boxlist.size());

    if (b + 1 == e) {
        boxlist[node_index] = cornerlist[b];
        return;
    }
    int m = b + (e - b) / 2;
    int childl = 2 * node_index;
    int childr = 2 * node_index + 1;

    assert(childl < boxlist.size());
    assert(childr < boxlist.size());

    init_boxes_recursive(cornerlist, childl, b, m);
    init_boxes_recursive(cornerlist, childr, m, e);

    assert(childl < boxlist.size());
    assert(childr < boxlist.size());
    for (int c = 0; c < 3; ++c) {
        boxlist[node_index][0][c] =
            std::min(boxlist[childl][0][c], boxlist[childr][0][c]);
        boxlist[node_index][1][c] =
            std::max(boxlist[childl][1][c], boxlist[childr][1][c]);
    }
}

void BVH::box_search_recursive(
    const Eigen::Vector3d& bbd0,
    const Eigen::Vector3d& bbd1,
    std::vector<unsigned int>& list,
    int n,
    int b,
    int e) const
{
    assert(e != b);

    assert(n < boxlist.size());
    bool cut = box_intersects_box(bbd0, bbd1, n);

    if (cut == false)
        return;

    // Leaf case
    if (e == b + 1) {
        list.emplace_back(b);
        return;
    }

    int m = b + (e - b) / 2;
    int childl = 2 * n;
    int childr = 2 * n + 1;

    // assert(childl < boxlist.size());
    // assert(childr < boxlist.size());

    // Traverse the "nearest" child first, so that it has more chances
    // to prune the traversal of the other child.
    box_search_recursive(bbd0, bbd1, list, childl, b, m);
    box_search_recursive(bbd0, bbd1, list, childr, m, e);
}

int BVH::max_node_index(int node_index, int b, int e)
{
    assert(e > b);
    if (b + 1 == e) {
        return node_index;
    }
    int m = b + (e - b) / 2;
    int childl = 2 * node_index;
    int childr = 2 * node_index + 1;
    return std::max(max_node_index(childl, b, m), max_node_index(childr, m, e));
}

void BVH::init(const std::vector<std::array<Eigen::Vector3d, 2>>& cornerlist)
{
    n_corners = cornerlist.size();

    Eigen::MatrixXd box_centers(n_corners, 3);
    for (int i = 0; i < n_corners; ++i) {
        box_centers.row(i) = (cornerlist[i][0] + cornerlist[i][1]) / 2;
    }

    const Eigen::RowVector3d vmin = box_centers.colwise().minCoeff();
    const Eigen::RowVector3d vmax = box_centers.colwise().maxCoeff();
    const Eigen::RowVector3d center = (vmin + vmax) / 2;
    for (int i = 0; i < n_corners; i++) {
        // make box centered at origin
        box_centers.row(i) -= center;
    }

    // after placing box at origin, vmax and vmin are symetric.
    const Eigen::Vector3d scale_point = vmax - center;
    const double scale = scale_point.lpNorm<Eigen::Infinity>();
    // if the box is too big, resize it
    if (scale > 100) {
        box_centers /= scale;
    }

    struct sortstruct {
        int order;
        Resorting::MortonCode64 morton;
    };
    std::vector<sortstruct> list;
    const int multi = 1000;
    list.resize(n_corners);

    for (int i = 0; i < n_corners; i++) {
        const Eigen::MatrixXd tmp = box_centers.row(i) * multi;

        list[i].morton =
            Resorting::MortonCode64(int(tmp(0)), int(tmp(1)), int(tmp(2)));
        list[i].order = i;
    }

    const auto morton_compare = [](const sortstruct& a, const sortstruct& b) {
        return (a.morton < b.morton);
    };
    std::sort(list.begin(), list.end(), morton_compare);

    new2old.resize(n_corners);
    for (int i = 0; i < n_corners; i++) {
        new2old[i] = list[i].order;
    }

    std::vector<std::array<Eigen::Vector3d, 2>> sorted_cornerlist(n_corners);

    for (int i = 0; i < n_corners; i++) {
        sorted_cornerlist[i] = cornerlist[list[i].order];
    }

    boxlist.resize(
        max_node_index(1, 0, n_corners)
        + 1 // <-- this is because size == max_index + 1 !!!
    );

    init_boxes_recursive(sorted_cornerlist, 1, 0, n_corners);
}

void BVH::init(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double tol)
{
    assert(F.cols() == 3);
    assert(V.cols() == 3);

    std::vector<std::array<Eigen::Vector3d, 2>> cornerlist(F.rows());

    for (int i = 0; i < F.rows(); i++) {
        const Eigen::RowVector3i face = F.row(i);
        const Eigen::RowVector3d v0 = V.row(face(0));
        const Eigen::RowVector3d v1 = V.row(face(1));
        const Eigen::RowVector3d v2 = V.row(face(2));

        Eigen::Matrix3d tmp;
        tmp.row(0) = v0;
        tmp.row(1) = v1;
        tmp.row(2) = v2;

        const Eigen::RowVector3d min = tmp.colwise().minCoeff().array() - tol;
        const Eigen::RowVector3d max = tmp.colwise().maxCoeff().array() + tol;

        cornerlist[i][0] = min.transpose();
        cornerlist[i][1] = max.transpose();
    }

    init(cornerlist);
}

bool BVH::box_intersects_box(
    const Eigen::Vector3d& bbd0, const Eigen::Vector3d& bbd1, int index) const
{
    const auto& bmin = boxlist[index][0];
    const auto& bmax = boxlist[index][1];

    return box_box_intersection(bbd0, bbd1, bmin, bmax);
}
} // namespace BVH
