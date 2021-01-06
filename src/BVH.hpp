#pragma once

#include <Eigen/Core>

#include <vector>
#include <array>
#include <cassert>

namespace BVH
{
	class BVH
	{
	private:
		std::vector<std::array<Eigen::Vector3d, 2>> boxlist;
		std::vector<int> new2old;
		size_t n_corners = -1;

		void init_boxes_recursive(const std::vector<std::array<Eigen::Vector3d, 2>> &cornerlist,
								  int node_index,
								  int b, int e);

		// void triangle_search_recursive(const Eigen::Vector3d &triangle0, const Eigen::Vector3d &triangle1, const Eigen::Vector3d &triangle2,
		// 							   std::vector<unsigned int> &list,
		// 							   int n, int b, int e) const;

		// void point_search_recursive(
		// 	const Eigen::Vector3d &point,
		// 	std::vector<unsigned int> &list,
		// 	int n, int b, int e) const;
		// void segment_search_recursive(
		// 	const Eigen::Vector3d &seg0, const Eigen::Vector3d &seg1,
		// 	std::vector<unsigned int> &list,
		// 	int n, int b, int e) const;

		void box_search_recursive(
			const Eigen::Vector3d &bbd0, const Eigen::Vector3d &bbd1,
			std::vector<unsigned int> &list,
			int n, int b, int e) const;

		static int max_node_index(int node_index, int b, int e);

		// bool triangle_intersects_box(const Eigen::Vector3d &tri0, const Eigen::Vector3d &tri1, const Eigen::Vector3d &tri2, int index) const;
		// bool point_intersects_box(const Eigen::Vector3d &p, int index) const;
		// bool segment_intersects_box(const Eigen::Vector3d &seg0, const Eigen::Vector3d &seg1, int index) const;
		bool box_intersects_box(const Eigen::Vector3d &bbd0, const Eigen::Vector3d &bbd1, int index) const;

		bool is_initialized = false;
		void init(const std::vector<std::array<Eigen::Vector3d, 2>> &cornerlist);

	public:
		void init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const double tol);

		// inline void triangle_find_bbox(
		// 	const Eigen::Vector3d &triangle0, const Eigen::Vector3d &triangle1, const Eigen::Vector3d &triangle2,
		// 	std::vector<unsigned int> &list) const
		// {
		// 	assert(n_corners >= 0);
		// 	assert(boxlist.size() > 0);
		// 	// int de = algorithms::is_triangle_degenerated(triangle0, triangle1, triangle2);
		// 	// if (de == DEGENERATED_SEGMENT)
		// 	// {
		// 	// 	Eigen::Vector3d tmin, tmax;
		// 	// 	algorithms::get_tri_corners(triangle0, triangle1, triangle2, tmin, tmax);
		// 	// 	segment_find_bbox(tmin, tmax, list);
		// 	// 	return;
		// 	// }
		// 	// else if (de == DEGENERATED_POINT)
		// 	// {
		// 	// 	point_find_bbox(triangle0, list);
		// 	// 	return;
		// 	// }

		// 	triangle_search_recursive(triangle0, triangle1, triangle2, list, 1, 0, n_corners);
		// }

		// inline void intersect_point(const Eigen::Vector3d &p, std::vector<unsigned int> &list) const
		// {
		// 	assert(boxlist.size() > 0);
		// 	point_search_recursive(p, list, 1, 0, n_corners);
		// }
		// inline void intersect_segment(const Eigen::Vector3d &seg0, const Eigen::Vector3d &seg1, std::vector<unsigned int> &list) const
		// {
		// 	assert(boxlist.size() > 0);
		// 	segment_search_recursive(
		// 		seg0, seg1, list, 1, 0, n_corners);
		// }
		inline void intersect_box(const Eigen::Vector3d &bbd0, const Eigen::Vector3d &bbd1, std::vector<unsigned int> &list) const
		{
			list.clear();
			assert(n_corners >= 0);
			box_search_recursive(bbd0, bbd1, list, 1, 0, n_corners);
		}
	};
} // namespace BVH