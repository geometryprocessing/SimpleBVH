#pragma once

#include <Eigen/Core>

#include <vector>
#include <array>
#include <cassert>

namespace SimpleBVH {

using VectorMax3d =
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, 3, 1>;
using LeafCallback =
    std::function<void(const VectorMax3d&, int, VectorMax3d&, double&)>;
using GetPointCallback = std::function<VectorMax3d(int)>;

void point_segment_squared_distance(
    const VectorMax3d& point,
    const std::array<VectorMax3d, 2>& f,
    VectorMax3d& closest_point,
    double& dist);

void point_triangle_squared_distance(
    const VectorMax3d& point,
    const std::array<VectorMax3d, 3>& f,
    VectorMax3d& closest_point,
    double& dist);

class BVH {
public:
    void init(const std::vector<std::array<VectorMax3d, 2>>& cornerlist);

    void
    init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double tol);

    void clear()
    {
        boxlist.clear();
        new2old.clear();
        n_corners = -1;
    }

    void intersect_box(
        const VectorMax3d& bbd0,
        const VectorMax3d& bbd1,
        std::vector<unsigned int>& list) const
    {
        assert(bbd0.size() == 2 || bbd0.size() == 3);
        assert(bbd0.size() == bbd1.size());

        std::vector<unsigned int> tmp;
        assert(n_corners >= 0);
        box_search_recursive(bbd0, bbd1, tmp, 1, 0, n_corners);

        list.resize(tmp.size());
        for (int i = 0; i < tmp.size(); ++i)
            list[i] = new2old[tmp[i]];
    }

    int nearest_facet(
        const VectorMax3d& p, VectorMax3d& nearest_point, double& sq_dist) const
    {
        int nearest_facet;
        get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        nearest_facet_recursive(
            p, nearest_facet, nearest_point, sq_dist, 1, 0, n_corners);
        return nearest_facet;
    }

    void nearest_facet_with_hint(
        const VectorMax3d& p,
        int& nearest_facet,
        VectorMax3d& nearest_point,
        double& sq_dist) const
    {
        if (nearest_facet < 0) {
            get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        }
        nearest_facet_recursive(
            p, nearest_facet, nearest_point, sq_dist, 1, 0, n_corners);
    }

    int facet_in_envelope(
        const VectorMax3d& p,
        double sq_epsilon,
        VectorMax3d& nearest_point,
        double& sq_dist) const
    {
        int nearest_facet;
        get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        facet_in_envelope_recursive(
            p, sq_epsilon, nearest_facet, nearest_point, sq_dist, 1, 0,
            n_corners);
        return nearest_facet;
    }

    void facet_in_envelope_with_hint(
        const VectorMax3d& p,
        double sq_epsilon,
        int& nearest_facet,
        VectorMax3d& nearest_point,
        double& sq_dist) const
    {
        if (nearest_facet < 0) {
            get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        }
        facet_in_envelope_recursive(
            p, sq_epsilon, nearest_facet, nearest_point, sq_dist, 1, 0,
            n_corners);
    }

    void set_leaf_distance_callback(const LeafCallback& callback)
    {
        leafCallback = callback;
    }

    void set_get_point_callback(const GetPointCallback& callback)
    {
        getPoint = callback;
    }

private:
    void init_boxes_recursive(
        const std::vector<std::array<VectorMax3d, 2>>& cornerlist,
        int node_index,
        int b,
        int e);

    void box_search_recursive(
        const VectorMax3d& bbd0,
        const VectorMax3d& bbd1,
        std::vector<unsigned int>& list,
        int n,
        int b,
        int e) const;

    void get_nearest_facet_hint(
        const VectorMax3d& p,
        int& nearest_f,
        VectorMax3d& nearest_point,
        double& sq_dist) const;

    void nearest_facet_recursive(
        const VectorMax3d& p,
        int& nearest_f,
        VectorMax3d& nearest_point,
        double& sq_dist,
        int n,
        int b,
        int e) const;

    void facet_in_envelope_recursive(
        const VectorMax3d& p,
        double sq_epsilon,
        int& nearest_f,
        VectorMax3d& nearest_point,
        double& sq_dist,
        int n,
        int b,
        int e) const;

    bool box_intersects_box(
        const VectorMax3d& bbd0, const VectorMax3d& bbd1, int index) const;

    static int max_node_index(int node_index, int b, int e);

    void leaf_callback(
        const VectorMax3d& p, int f, VectorMax3d& np, double& sq_d) const;
    VectorMax3d point_callback(int f) const;

    std::vector<std::array<VectorMax3d, 2>> boxlist;
    std::vector<int> new2old;
    long n_corners = -1;
    LeafCallback leafCallback;
    GetPointCallback getPoint;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};
} // namespace SimpleBVH
