#include <BVH.hpp>
#include <Morton.hpp>

#include <iostream>

namespace BVH
{
	namespace
	{
		bool box_box_intersection(const Eigen::Vector3d &min1, const Eigen::Vector3d &max1,
								  const Eigen::Vector3d &min2, const Eigen::Vector3d &max2)
		{
			if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
				return 0;
			if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
				return 0;
			return 1;
		}
	} // namespace

	void BVH::init_boxes_recursive(const std::vector<std::array<Eigen::Vector3d, 2>> &cornerlist,
								   int node_index, int b, int e)
	{
		assert(b != e);
		assert(node_index < boxlist.size());

		if (b + 1 == e)
		{
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
		for (int c = 0; c < 3; ++c)
		{
			boxlist[node_index][0][c] = std::min(boxlist[childl][0][c], boxlist[childr][0][c]);
			boxlist[node_index][1][c] = std::max(boxlist[childl][1][c], boxlist[childr][1][c]);
		}
	}

	// void BVH::triangle_search_recursive(const Eigen::Vector3d &triangle0, const Eigen::Vector3d &triangle1, const Eigen::Vector3d &triangle2, std::vector<unsigned int> &list,
	// 									int n, int b, int e) const
	// {
	// 	assert(e != b);

	// 	assert(n < boxlist.size());
	// 	bool cut = triangle_intersects_box(triangle0, triangle1, triangle2, n);

	// 	if (cut == false)
	// 		return;

	// 	// Leaf case
	// 	if (e == b + 1)
	// 	{
	// 		list.emplace_back(b);
	// 		return;
	// 	}

	// 	int m = b + (e - b) / 2;
	// 	int childl = 2 * n;
	// 	int childr = 2 * n + 1;

	// 	//assert(childl < boxlist.size());
	// 	//assert(childr < boxlist.size());

	// 	// Traverse the "nearest" child first, so that it has more chances
	// 	// to prune the traversal of the other child.
	// 	triangle_search_recursive(
	// 		triangle0, triangle1, triangle2, list,
	// 		childl, b, m);
	// 	triangle_search_recursive(
	// 		triangle0, triangle1, triangle2, list,
	// 		childr, m, e);
	// }

	// void BVH::point_search_recursive(const Eigen::Vector3d &point, std::vector<unsigned int> &list,
	// 								 int n, int b, int e) const
	// {
	// 	assert(e != b);

	// 	assert(n < boxlist.size());
	// 	bool cut = point_intersects_box(point, n);

	// 	if (cut == false)
	// 		return;

	// 	// Leaf case
	// 	if (e == b + 1)
	// 	{
	// 		list.emplace_back(b);
	// 		return;
	// 	}

	// 	int m = b + (e - b) / 2;
	// 	int childl = 2 * n;
	// 	int childr = 2 * n + 1;

	// 	//assert(childl < boxlist.size());
	// 	//assert(childr < boxlist.size());

	// 	// Traverse the "nearest" child first, so that it has more chances
	// 	// to prune the traversal of the other child.
	// 	point_search_recursive(
	// 		point, list,
	// 		childl, b, m);
	// 	point_search_recursive(
	// 		point, list,
	// 		childr, m, e);
	// }

	// void BVH::segment_search_recursive(const Eigen::Vector3d &seg0, const Eigen::Vector3d &seg1, std::vector<unsigned int> &list,
	// 								   int n, int b, int e) const
	// {
	// 	assert(e != b);

	// 	assert(n < boxlist.size());
	// 	bool cut = segment_intersects_box(seg0, seg1, n);

	// 	if (cut == false)
	// 		return;

	// 	// Leaf case
	// 	if (e == b + 1)
	// 	{
	// 		list.emplace_back(b);
	// 		return;
	// 	}

	// 	int m = b + (e - b) / 2;
	// 	int childl = 2 * n;
	// 	int childr = 2 * n + 1;

	// 	//assert(childl < boxlist.size());
	// 	//assert(childr < boxlist.size());

	// 	// Traverse the "nearest" child first, so that it has more chances
	// 	// to prune the traversal of the other child.
	// 	segment_search_recursive(
	// 		seg0, seg1, list,
	// 		childl, b, m);
	// 	segment_search_recursive(
	// 		seg0, seg1, list,
	// 		childr, m, e);
	// }

	void BVH::box_search_recursive(const Eigen::Vector3d &bbd0, const Eigen::Vector3d &bbd1, std::vector<unsigned int> &list,
								   int n, int b, int e) const
	{
		assert(e != b);

		assert(n < boxlist.size());
		bool cut = box_intersects_box(bbd0, bbd1, n);

		if (cut == false)
			return;

		// Leaf case
		if (e == b + 1)
		{
			list.emplace_back(b);
			return;
		}

		int m = b + (e - b) / 2;
		int childl = 2 * n;
		int childr = 2 * n + 1;

		//assert(childl < boxlist.size());
		//assert(childr < boxlist.size());

		// Traverse the "nearest" child first, so that it has more chances
		// to prune the traversal of the other child.
		box_search_recursive(
			bbd0, bbd1, list,
			childl, b, m);
		box_search_recursive(
			bbd0, bbd1, list,
			childr, m, e);
	}

	int BVH::max_node_index(int node_index, int b, int e)
	{
		assert(e > b);
		if (b + 1 == e)
		{
			return node_index;
		}
		int m = b + (e - b) / 2;
		int childl = 2 * node_index;
		int childr = 2 * node_index + 1;
		return std::max(
			max_node_index(childl, b, m),
			max_node_index(childr, m, e));
	}

	void BVH::init(const std::vector<std::array<Eigen::Vector3d, 2>> &cornerlist)
	{
		n_corners = cornerlist.size();

		boxlist.resize(
			max_node_index(
				1, 0, n_corners) +
			1 // <-- this is because size == max_index + 1 !!!
		);

		init_boxes_recursive(cornerlist, 1, 0, n_corners);
	}

	void BVH::init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const double tol)
	{
		assert(F.cols() == 3);
		assert(V.cols() == 3);

		struct sortstruct
		{
			int order;
			Resorting::MortonCode64 morton;
		};
		std::vector<sortstruct> list;
		const int multi = 1000;

		const Eigen::RowVector3d vmin = V.colwise().minCoeff();
		const Eigen::RowVector3d vmax = V.colwise().maxCoeff();
		const Eigen::RowVector3d center = (vmin + vmax) / 2;
		Eigen::MatrixXd v_shifted(V.rows(), V.cols());

		for (int i = 0; i < V.rows(); i++)
		{
			// make box centered at origin
			v_shifted.row(i) = V.row(i) - center;
		}

		// after placing box at origin, vmax and vmin are symetric.
		const Eigen::Vector3d scale_point = vmax - center;
		const double scale = scale_point.lpNorm<Eigen::Infinity>();
		// if the box is too big, resize it
		if (scale > 300)
		{
			v_shifted /= scale;
		}

		list.resize(F.rows());

		for (int i = 0; i < F.rows(); i++)
		{
			const Eigen::MatrixXd tmp = (v_shifted.row(F(i, 0)) + v_shifted.row(F(i, 1)) + v_shifted.row(F(i, 2))) * multi;

			list[i].morton = Resorting::MortonCode64(int(tmp(0)), int(tmp(1)), int(tmp(2)));
			list[i].order = i;
		}

		const auto morton_compare = [](const sortstruct &a, const sortstruct &b) {
			return (a.morton < b.morton);
		};
		std::sort(list.begin(), list.end(), morton_compare);

		new2old.resize(F.rows());
		for (int i = 0; i < F.rows(); i++)
		{
			new2old[i] = list[i].order;
		}

		std::vector<std::array<Eigen::Vector3d, 2>> cornerlist(F.rows());

		for (int i = 0; i < F.rows(); i++)
		{
			const Eigen::RowVector3i face = F.row(list[i].order);
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

	// bool BVH::triangle_intersects_box(const Eigen::Vector3d &tri0, const Eigen::Vector3d &tri1, const Eigen::Vector3d &tri2, int index) const
	// {
	// 	const auto &bmin = boxlist[index][0];
	// 	const auto &bmax = boxlist[index][1];
	// 	Eigen::Vector3d tmin, tmax;

	// 	algorithms::get_tri_corners(tri0, tri1, tri2, tmin, tmax);
	// 	bool cut = algorithms::box_box_intersection(tmin, tmax, bmin, bmax);
	// 	if (cut == false)
	// 		return false;

	// 	if (cut)
	// 	{

	// 		std::array<Vector2, 3> tri;
	// 		std::array<Vector2, 4> mp;
	// 		int o0, o1, o2, o3, ori;
	// 		for (int i = 0; i < 3; i++)
	// 		{
	// 			tri[0] = algorithms::to_2d(tri0, i);
	// 			tri[1] = algorithms::to_2d(tri1, i);
	// 			tri[2] = algorithms::to_2d(tri2, i);

	// 			mp[0] = algorithms::to_2d(bmin, i);
	// 			mp[1] = algorithms::to_2d(bmax, i);
	// 			mp[2][0] = mp[0][0];
	// 			mp[2][1] = mp[1][1];
	// 			mp[3][0] = mp[1][0];
	// 			mp[3][1] = mp[0][1];

	// 			for (int j = 0; j < 3; j++)
	// 			{
	// 				o0 = fastEnvelope::Predicates::orient_2d(mp[0], tri[j % 3], tri[(j + 1) % 3]);
	// 				o1 = fastEnvelope::Predicates::orient_2d(mp[1], tri[j % 3], tri[(j + 1) % 3]);
	// 				o2 = fastEnvelope::Predicates::orient_2d(mp[2], tri[j % 3], tri[(j + 1) % 3]);
	// 				o3 = fastEnvelope::Predicates::orient_2d(mp[3], tri[j % 3], tri[(j + 1) % 3]);
	// 				ori = fastEnvelope::Predicates::orient_2d(tri[(j + 2) % 3], tri[j % 3], tri[(j + 1) % 3]);
	// 				if (ori == 0)
	// 					continue;
	// 				if (ori * o0 <= 0 && ori * o1 <= 0 && ori * o2 <= 0 && ori * o3 <= 0)
	// 					return false;
	// 			}
	// 		}
	// 	}

	// 	return cut;
	// }

	// bool BVH::point_intersects_box(const Eigen::Vector3d &p, int index) const
	// {
	// 	const auto &bmin = boxlist[index][0];
	// 	const auto &bmax = boxlist[index][1];
	// 	if (p[0] < bmin[0] || p[1] < bmin[1] || p[2] < bmin[2])
	// 		return false;
	// 	if (p[0] > bmax[0] || p[1] > bmax[1] || p[2] > bmax[2])
	// 		return false;
	// 	return true;
	// }

	// bool BVH::segment_intersects_box(const Eigen::Vector3d &seg0, const Eigen::Vector3d &seg1, int index) const
	// {
	// 	const auto &bmin = boxlist[index][0];
	// 	const auto &bmax = boxlist[index][1];
	// 	Scalar min[3], max[3];
	// 	min[0] = std::min(seg0[0], seg1[0]);
	// 	min[1] = std::min(seg0[1], seg1[1]);
	// 	min[2] = std::min(seg0[2], seg1[2]);
	// 	max[0] = std::max(seg0[0], seg1[0]);
	// 	max[1] = std::max(seg0[1], seg1[1]);
	// 	max[2] = std::max(seg0[2], seg1[2]);
	// 	if (max[0] < bmin[0] || max[1] < bmin[1] || max[2] < bmin[2])
	// 		return false;
	// 	if (min[0] > bmax[0] || min[1] > bmax[1] || min[2] > bmax[2])
	// 		return false;
	// 	return true;
	// }

	bool BVH::box_intersects_box(const Eigen::Vector3d &bbd0, const Eigen::Vector3d &bbd1, int index) const
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];

		return box_box_intersection(bbd0, bbd1, bmin, bmax);
	}
} // namespace BVH
