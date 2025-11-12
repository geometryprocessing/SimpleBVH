#pragma once

#include <Eigen/Core>

#include <vector>
#include <array>
#include <cassert>

namespace SimpleBVH {
using Vector3d = Eigen::Vector3d;
using Vector2d = Eigen::Vector2d;
using LeafCallback =
    std::function<void(const Vector3d&, int, Vector3d&, double&)>;
using GetPointCallback = std::function<Vector3d(int)>;

inline Vector2d to_2d(const Vector3d& x) { return Vector2d(x[0], x[1]); }
inline Vector3d to_3d(const Vector2d& x) { return Vector3d(x[0], x[1], 0); }

void point_segment_squared_distance(
    const Vector3d& point,
    const Vector3d& pa,
    const Vector3d& pb,
    Vector3d& closest_point,
    double& dist);

void point_triangle_squared_distance(
    const Vector3d& point,
    const Vector3d& pa,
    const Vector3d& pb,
    const Vector3d& pc,
    Vector3d& closest_point,
    double& dist);

class BVH {
public:
    void init(const std::vector<std::array<Vector3d, 2>>& cornerlist);

    void
    init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double tol);

    void clear()
    {
        boxlist.clear();
        new2old.clear();
        n_corners = -1;
    }

    void intersect_box(
        const Vector3d& bbd0,
        const Vector3d& bbd1,
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
        const Vector3d& p, Vector3d& nearest_point, double& sq_dist) const
    {
        int nearest_facet;
        get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        nearest_facet_recursive(
            p, nearest_facet, nearest_point, sq_dist, 1, 0, n_corners);
        return nearest_facet;
    }

    int nearest_facet(
        const Vector2d& p, Vector2d& nearest_point, double& sq_dist) const
    {
        Vector3d np3 = to_3d(nearest_point);
        Vector3d p3 = to_3d(p);
        int nf = nearest_facet(p3, np3, sq_dist);
        nearest_point = to_2d(np3);
        return nf;
    }

    void nearest_facet_with_hint(
        const Vector3d& p,
        int& nearest_facet,
        Vector3d& nearest_point,
        double& sq_dist) const
    {
        if (nearest_facet < 0) {
            get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        }
        nearest_facet_recursive(
            p, nearest_facet, nearest_point, sq_dist, 1, 0, n_corners);
    }

    void nearest_facet_with_hint(
        const Vector2d& p,
        int& nearest_facet,
        Vector2d& nearest_point,
        double& sq_dist) const
    {
        Vector3d np3 = to_3d(nearest_point);
        Vector3d p3 = to_3d(p);
        nearest_facet_with_hint(p3, nearest_facet, np3, sq_dist);
        nearest_point = to_2d(np3);
    }

    int facet_in_envelope(
        const Vector3d& p,
        double sq_epsilon,
        Vector3d& nearest_point,
        double& sq_dist) const
    {
        int nearest_facet;
        get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        facet_in_envelope_recursive(
            p, sq_epsilon, nearest_facet, nearest_point, sq_dist, 1, 0,
            n_corners);
        return nearest_facet;
    }

    int facet_in_envelope(
        const Vector2d& p,
        double sq_epsilon,
        Vector2d& nearest_point,
        double& sq_dist) const
    {
        Vector3d np3 = to_3d(nearest_point);
        Vector3d p3 = to_3d(p);
        int nearest_facet = facet_in_envelope(p3, sq_epsilon, np3, sq_dist);
        nearest_point = to_2d(np3);
        return nearest_facet;
    }

    void facet_in_envelope_with_hint(
        const Vector3d& p,
        double sq_epsilon,
        int& nearest_facet,
        Vector3d& nearest_point,
        double& sq_dist) const
    {
        if (nearest_facet < 0) {
            get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
        }

        facet_in_envelope_recursive(
            p, sq_epsilon, nearest_facet, nearest_point, sq_dist, 1, 0,
            n_corners);
    }

    void facet_in_envelope_with_hint(
        const Vector2d& p,
        double sq_epsilon,
        int& nearest_facet,
        Vector2d& nearest_point,
        double& sq_dist) const
    {
        Vector3d np3 = to_3d(nearest_point);
        Vector3d p3 = to_3d(p);
        facet_in_envelope_with_hint(
            p3, sq_epsilon, nearest_facet, np3, sq_dist);
        nearest_point = to_2d(np3);
    }

    void set_leaf_distance_callback(const LeafCallback& callback)
    {
        leafCallback = callback;
    }

    void set_get_point_callback(const GetPointCallback& callback)
    {
        getPoint = callback;
    }

    inline void point_facet_distance(
        const Vector3d& p,
        const int facet,
        Vector3d& closest_point,
        double& sq_d) const
    {
        const int f = new2old[facet];

        if (faces.cols() == 2) {
            point_segment_squared_distance(
                p, vertices.row(faces(f, 0)), vertices.row(faces(f, 1)),
                closest_point, sq_d);
            return;
        } else if (faces.cols() == 3) {
            point_triangle_squared_distance(
                p, vertices.row(faces(f, 0)), vertices.row(faces(f, 1)),
                vertices.row(faces(f, 2)), closest_point, sq_d);
            return;
        }
    }

    inline void point_facet_distance(
        const Vector2d& p,
        const int facet,
        Vector2d& closest_point,
        double& sq_d) const
    {
        Vector3d np3 = to_3d(closest_point);
        Vector3d p3 = to_3d(p);

        point_facet_distance(p3, facet, np3, sq_d);
        closest_point = to_2d(np3);
    }

private:
    void init_boxes_recursive(
        const std::vector<std::array<Vector3d, 2>>& cornerlist,
        int node_index,
        int b,
        int e);

    void box_search_recursive(
        const Vector3d& bbd0,
        const Vector3d& bbd1,
        std::vector<unsigned int>& list,
        int n,
        int b,
        int e) const;

    void get_nearest_facet_hint(
        const Vector3d& p,
        int& nearest_f,
        Vector3d& nearest_point,
        double& sq_dist) const;

    void nearest_facet_recursive(
        const Vector3d& p,
        int& nearest_f,
        Vector3d& nearest_point,
        double& sq_dist,
        int n,
        int b,
        int e) const;

    void facet_in_envelope_recursive(
        const Vector3d& p,
        double sq_epsilon,
        int& nearest_f,
        Vector3d& nearest_point,
        double& sq_dist,
        int n,
        int b,
        int e) const;

    bool box_intersects_box(
        const Vector3d& bbd0, const Vector3d& bbd1, int index) const;

    static int max_node_index(int node_index, int b, int e);

    void
    leaf_callback(const Vector3d& p, int f, Vector3d& np, double& sq_d) const;
    Vector3d point_callback(int f) const;

    std::vector<std::array<Vector3d, 2>> boxlist;
    std::vector<int> new2old;
    long n_corners = -1;
    LeafCallback leafCallback;
    GetPointCallback getPoint;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};
} // namespace SimpleBVH
