#pragma once

#include <Eigen/Core>

#include <vector>
#include <array>
#include <cassert>

namespace BVH {
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

    inline void intersect_box(
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
};
} // namespace BVH
