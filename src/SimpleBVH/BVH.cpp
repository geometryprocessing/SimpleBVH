#include "BVH.hpp"
#include <SimpleBVH/Morton.hpp>

#include <iostream>

namespace SimpleBVH {

void point_segment_squared_distance(
    const VectorMax3d& point,
    const std::array<VectorMax3d, 2>& f,
    VectorMax3d& closest_point,
    double& dist)
{
    const double l2 = (f[0] - f[1]).norm();
    const double t = (point - f[0]).dot(f[1] - f[0]);
    if (t <= 0.0 || l2 == 0.0) {
        closest_point = f[0];
    } else if (t > l2) {
        closest_point = f[1];
    } else {
        const double lambda1 = t / l2;
        const double lambda0 = 1.0 - lambda1;
        closest_point = lambda0 * f[0] + lambda1 * f[1];
    }
    dist = (point - closest_point).squaredNorm();
}

namespace {
    bool box_box_intersection(
        const VectorMax3d& min1,
        const VectorMax3d& max1,
        const VectorMax3d& min2,
        const VectorMax3d& max2)
    {
        if (min1.size() == 3) {
            if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
                return 0;
            if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
                return 0;
            return 1;
        } else {
            if (max1[0] < min2[0] || max1[1] < min2[1])
                return 0;
            if (max2[0] < min1[0] || max2[1] < min1[1])
                return 0;
            return 1;
        }
    }

    double point_box_center_squared_distance(
        const VectorMax3d& p, const std::array<VectorMax3d, 2>& B)
    {
        return (p - (B[0] + B[1]) / 2).norm();
    }

    double inner_point_box_squared_distance(
        const VectorMax3d& p, const std::array<VectorMax3d, 2>& B)
    {
        assert(p.size() == B[0].size());

        double result = std::pow(p[0] - B[0][0], 2);
        result = std::min(result, std::pow(p[0] - B[1][0], 2));
        for (int c = 1; c < p.size(); ++c) {
            result = std::min(result, std::pow(p[c] - B[0][c], 2));
            result = std::min(result, std::pow(p[c] - B[1][c], 2));
        }
        return result;
    }

    double point_box_signed_squared_distance(
        const VectorMax3d& p, const std::array<VectorMax3d, 2>& B)
    {
        assert(p.size() == B[0].size());

        bool inside = true;
        double result = 0.0;
        for (int c = 0; c < p.size(); c++) {
            if (p[c] < B[0][c]) {
                inside = false;
                result += std::pow(p[c] - B[0][c], 2);
            } else if (p[c] > B[1][c]) {
                inside = false;
                result += std::pow(p[c] - B[1][c], 2);
            }
        }
        if (inside) {
            result = -inner_point_box_squared_distance(p, B);
        }
        return result;
    }
} // namespace

void BVH::init_boxes_recursive(
    const std::vector<std::array<VectorMax3d, 2>>& cornerlist,
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
    for (int c = 0; c < cornerlist[0].size(); ++c) {
        boxlist[node_index][0][c] =
            std::min(boxlist[childl][0][c], boxlist[childr][0][c]);
        boxlist[node_index][1][c] =
            std::max(boxlist[childl][1][c], boxlist[childr][1][c]);
    }
}

void BVH::box_search_recursive(
    const VectorMax3d& bbd0,
    const VectorMax3d& bbd1,
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

void BVH::init(const std::vector<std::array<VectorMax3d, 2>>& cornerlist)
{
    n_corners = cornerlist.size();

    Eigen::MatrixXd box_centers(n_corners, cornerlist[0][0].size());
    for (int i = 0; i < n_corners; ++i) {
        box_centers.row(i) = (cornerlist[i][0] + cornerlist[i][1]) / 2;
    }

    const VectorMax3d vmin = box_centers.colwise().minCoeff();
    const VectorMax3d vmax = box_centers.colwise().maxCoeff();
    const VectorMax3d center = (vmin + vmax) / 2;
    for (int i = 0; i < n_corners; i++) {
        // make box centered at origin
        box_centers.row(i) -= center;
    }

    // after placing box at origin, vmax and vmin are symetric.
    const VectorMax3d scale_point = vmax - center;
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

        list[i].morton = Resorting::MortonCode64(
            int(tmp(0)), int(tmp(1)), tmp.size() == 3 ? int(tmp(2)) : 0);
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

    std::vector<std::array<VectorMax3d, 2>> sorted_cornerlist(n_corners);

    for (int i = 0; i < n_corners; i++) {
        sorted_cornerlist[i] = cornerlist[list[i].order];
    }

    boxlist.resize(
        max_node_index(1, 0, n_corners)
        + 1 // <-- this is because size == max_index + 1 !!!
    );
    for (auto& b : boxlist) {
        b[0].resize(cornerlist[0][0].size());
        b[1].resize(cornerlist[0][0].size());
    }

    init_boxes_recursive(sorted_cornerlist, 1, 0, n_corners);
}

void BVH::init(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double tol)
{
    assert(F.cols() == 3 || F.cols() == 2);
    assert(V.cols() == 3 || V.cols() == 2);

    std::vector<std::array<VectorMax3d, 2>> cornerlist(F.rows());
    if (F.cols() == 3) {
        for (int i = 0; i < F.rows(); i++) {
            const Eigen::RowVector3i face = F.row(i);
            const VectorMax3d v0 = V.row(face(0));
            const VectorMax3d v1 = V.row(face(1));
            const VectorMax3d v2 = V.row(face(2));

            Eigen::MatrixXd tmp(3, v0.size());
            tmp.row(0) = v0.transpose();
            tmp.row(1) = v1.transpose();
            tmp.row(2) = v2.transpose();

            const VectorMax3d min = tmp.colwise().minCoeff().array() - tol;
            const VectorMax3d max = tmp.colwise().maxCoeff().array() + tol;

            cornerlist[i][0] = min.transpose();
            cornerlist[i][1] = max.transpose();
        }
    } else if (F.cols() == 2) {
        for (int i = 0; i < F.rows(); i++) {
            const Eigen::RowVector2i face = F.row(i);
            const VectorMax3d v0 = V.row(face(0));
            const VectorMax3d v1 = V.row(face(1));

            Eigen::MatrixXd tmp(2, v0.size());
            tmp.row(0) = v0.transpose();
            tmp.row(1) = v1.transpose();

            const VectorMax3d min = tmp.colwise().minCoeff().array() - tol;
            const VectorMax3d max = tmp.colwise().maxCoeff().array() + tol;

            cornerlist[i][0] = min.transpose();
            cornerlist[i][1] = max.transpose();
        }
    }

    init(cornerlist);
}

bool BVH::box_intersects_box(
    const VectorMax3d& bbd0, const VectorMax3d& bbd1, int index) const
{
    const auto& bmin = boxlist[index][0];
    const auto& bmax = boxlist[index][1];

    return box_box_intersection(bbd0, bbd1, bmin, bmax);
}

void BVH::get_nearest_facet_hint(
    const VectorMax3d& p,
    int& nearest_f,
    VectorMax3d& nearest_point,
    double& sq_dist) const
{
    int b = 0;
    int e = n_corners;
    int n = 1;
    while (e != b + 1) {
        int m = b + (e - b) / 2;
        int childl = 2 * n;
        int childr = 2 * n + 1;
        if (point_box_center_squared_distance(p, boxlist[childl])
            < point_box_center_squared_distance(p, boxlist[childr])) {
            e = m;
            n = childl;
        } else {
            b = m;
            n = childr;
        }
    }
    nearest_f = b;

    nearest_point = getPoint(nearest_f);
    sq_dist = (p - nearest_point).norm();
}

void BVH::nearest_facet_recursive(
    const VectorMax3d& p,
    int& nearest_f,
    VectorMax3d& nearest_point,
    double& sq_dist,
    int n,
    int b,
    int e) const
{
    assert(e > b);

    // If node is a leaf: compute point-facet distance
    // and replace current if nearer
    if (b + 1 == e) {
        VectorMax3d cur_nearest_point;
        double cur_sq_dist;
        leafCallback(p, new2old[b], cur_nearest_point, cur_sq_dist);
        if (cur_sq_dist < sq_dist) {
            nearest_f = b;
            nearest_point = cur_nearest_point;
            sq_dist = cur_sq_dist;
        }
        return;
    }
    int m = b + (e - b) / 2;
    int childl = 2 * n;
    int childr = 2 * n + 1;

    assert(childl < boxlist.size());
    assert(childr < boxlist.size());

    double dl = point_box_signed_squared_distance(p, boxlist[childl]);
    double dr = point_box_signed_squared_distance(p, boxlist[childr]);

    // Traverse the "nearest" child first, so that it has more chances
    // to prune the traversal of the other child.
    if (dl < dr) {
        if (dl < sq_dist) {
            nearest_facet_recursive(
                p, nearest_f, nearest_point, sq_dist, childl, b, m);
        }
        if (dr < sq_dist) {
            nearest_facet_recursive(
                p, nearest_f, nearest_point, sq_dist, childr, m, e);
        }
    } else {
        if (dr < sq_dist) {
            nearest_facet_recursive(
                p, nearest_f, nearest_point, sq_dist, childr, m, e);
        }
        if (dl < sq_dist) {
            nearest_facet_recursive(
                p, nearest_f, nearest_point, sq_dist, childl, b, m);
        }
    }
}
} // namespace SimpleBVH
