#pragma once

#include <Eigen/Core>

#include <vector>
#include <array>
#include <cassert>

namespace BVH {

using VectorMax3d =
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, 3, 1>;

class BVH {
private:
    std::vector<std::array<Eigen::Vector3d, 2>> boxlist;
    std::vector<int> new2old;
    size_t n_corners = -1;

    void init_boxes_recursive(
        const std::vector<std::array<Eigen::Vector3d, 2>>& cornerlist,
        int node_index,
        int b,
        int e);

    void box_search_recursive(
        const Eigen::Vector3d& bbd0,
        const Eigen::Vector3d& bbd1,
        std::vector<unsigned int>& list,
        int n,
        int b,
        int e) const;

    static int max_node_index(int node_index, int b, int e);

    bool box_intersects_box(
        const Eigen::Vector3d& bbd0,
        const Eigen::Vector3d& bbd1,
        int index) const;

    bool is_initialized = false;

public:
    void init(const std::vector<std::array<Eigen::Vector3d, 2>>& cornerlist);
    void
    init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double tol);

    inline void intersect_3D_box(
        const Eigen::Vector3d& bbd0,
        const Eigen::Vector3d& bbd1,
        std::vector<unsigned int>& list) const
    {
        std::vector<unsigned int> tmp;
        assert(n_corners >= 0);
        box_search_recursive(bbd0, bbd1, tmp, 1, 0, n_corners);

        list.resize(tmp.size());
        for (int i = 0; i < tmp.size(); ++i)
            list[i] = new2old[tmp[i]];
    }

    inline void intersect_box(
        const VectorMax3d& bbd0,
        const VectorMax3d& bbd1,
        std::vector<unsigned int>& list) const
    {
        Eigen::Vector3d bbd0_3D = Eigen::Vector3d::Zero();
        bbd0_3D.head(bbd0.size()) = bbd0.head(bbd0.size());

        Eigen::Vector3d bbd1_3D = Eigen::Vector3d::Zero();
        bbd1_3D.head(bbd1.size()) = bbd1.head(bbd1.size());

        intersect_3D_box(bbd0_3D, bbd1_3D, list);
    }

    inline void intersect_2D_box(
        const Eigen::Vector2d& bbd0,
        const Eigen::Vector2d& bbd1,
        std::vector<unsigned int>& list) const
    {
        Eigen::Vector3d bbd0_3D = Eigen::Vector3d::Zero();
        bbd0_3D.head<2>() = bbd0;

        Eigen::Vector3d bbd1_3D = Eigen::Vector3d::Zero();
        bbd1_3D.head<2>() = bbd1;

        intersect_3D_box(bbd0_3D, bbd1_3D, list);
    }
};
} // namespace BVH
