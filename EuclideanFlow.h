#pragma once

#ifndef _EUCLIDEAN_FLOW_H_
#define _EUCLIDEAN_FLOW_H_

#include "BaseFlow.h"


#ifndef NEWTON_ERROR
#define NEWTON_ERROR 1e-4
#endif


namespace OpenMesh

{
	/*! \brief Class EuclideanFlow
	*
	*	Algorithm for computing Ricci flow
	*/

	class EuclideanFlow : public BaseFlow
	{
	public:
		/*! \brief CEucideanFlow constructor
		*  \param pMesh the input mesh
		*/
		EuclideanFlow(MeshFlow* pMesh);
		/*! \brief EuclideanFlow destructor
		*/
		~EuclideanFlow() {};

		virtual void CalculateMetric() = 0;

	protected:

		/*!
		*	Calculate each edge Length, has to be defined in the derivated classes
		*/
		virtual void Length(double u1, double u2, MeshFlow::EHandle e);

		/*!
		*	Cosine law, has to be defined in the derivated classes
		*/
		double CosineLaw(double a, double b, double c);

		/*!
		*	Calculate the edge weight
		*/
		virtual void CalculateEdgeWeight();
		/*embed a face to get height*/
		double HeightOnHalfedge(MeshFlow::HHandle h);
		/*!
		*	Set the target curvature on each vertex
		*/
		virtual void SetTargetCurvature();

		virtual double RToU(double r);
		virtual double UToR(double u);
		virtual void ComputeMetricWeight(MeshFlow::EHandle e);

		void Newton(double threshold, double stepLength);
		virtual void AcceleratedGradientDescent(double threshold, double stepLength);
		void CalculateHessianMatrix(Eigen::SparseMatrix<double>& pMatrix);
	};



	inline EuclideanFlow::EuclideanFlow(MeshFlow* pMesh) :BaseFlow(pMesh)
	{
	};

	//Compute the edge Length

	inline void EuclideanFlow::Length(double u1, double u2, MeshFlow::EHandle e)
	{
		double r1 = UToR(u1);
		//assert(abs(RToU(r1) - u1) < 1e-7);
		double r2 = UToR(u2);
		m_pMesh->Length(e) = sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * m_pMesh->MetricWeight(e));
	};


	//Calculate corner angle

	inline double EuclideanFlow::CosineLaw(double a, double b, double c)
	{
		double cs = (a * a + b * b - c * c) / (2.0 * a * b);
		if (cs <= 1.0 && cs >= -1.0) {
			return acos(cs);
		}

		else if (cs > 1.0) {
			return 0;
		}
		else {
			return PI;
		}
	};


	//Calculate edge weight


	inline void EuclideanFlow::CalculateEdgeWeight()
	{

		for (MeshFlow::EdgeIter eiter = m_pMesh->edges_begin(); eiter != m_pMesh->edges_end(); eiter++)
		{
			MeshFlow::EHandle e = *eiter;
			double Length = m_pMesh->Length(e);
			double weight = 0;
			MeshFlow::HalfedgeHandle h0 = m_pMesh->EdgeHalfedge(e, 0);
			MeshFlow::HalfedgeHandle h1 = m_pMesh->EdgeHalfedge(e, 1);
			if (!m_pMesh->is_boundary(h0))
				weight += HeightOnHalfedge(h0) / Length;
			if (!m_pMesh->is_boundary(h1))
				weight += HeightOnHalfedge(h1) / Length;
			m_pMesh->Weight(e) = weight;
		}
	}


	inline double EuclideanFlow::HeightOnHalfedge(MeshFlow::HHandle h)
	{
		double angle = m_pMesh->Angle(h);
		if (abs(angle - PI) < 1e-7 || abs(angle) < 1e-7) {
			return 0.0;
		}
		MeshFlow::EHandle e[3];
		MeshFlow::HHandle he[3];
		he[0] = h;
		e[0] = m_pMesh->edge_handle(he[0]);
		he[1] = m_pMesh->next_halfedge_handle(he[0]);
		e[1] = m_pMesh->edge_handle(he[1]);
		he[2] = m_pMesh->next_halfedge_handle(he[1]);
		e[2] = m_pMesh->edge_handle(he[2]);


		MeshLib::CPoint2 uv[3];

		uv[0] = MeshLib::CPoint2(0, 0);
		uv[1] = MeshLib::CPoint2(m_pMesh->Length(e[0]), 0);

		MeshLib::CPoint2 c1, c2;

		bool result = _circle_circle_intersection(MeshLib::CCircle(uv[0], m_pMesh->Length(e[2])),
			MeshLib::CCircle(uv[1], m_pMesh->Length(e[1])),
			c1, c2);
		assert(result);
		if (MeshLib::cross(uv[1] - uv[0], c1 - uv[0]) > 0)
		{
			uv[2] = c1;
		}
		else
		{
			uv[2] = c2;
		}

		std::vector<MeshFlow::VertexHandle> av;
		MeshLib::CCircle C[3];
		for (int i = 0; i < 3; ++i) {
			av.push_back(m_pMesh->from_vertex_handle(he[i]));
			C[i] = MeshLib::CCircle(uv[i], UToR(m_pMesh->U(av[i])));
		}
		MeshLib::CCircle otho = MeshLib::orthogonal(C);
		return otho.c()[1];
	}


	//set target curvature
	inline void EuclideanFlow::SetTargetCurvature()
	{

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); ++viter) {
			VertexHandle v = *viter;
			m_pMesh->TargetCurvature(v) = 0.0;
		}
		if (orbifold) return;
		m_pMesh->RequestBoundary();
		if (boundary_mode == Fixed && m_pMesh->boundaries.size() != 0) {
			std::vector<MeshFlow::HHandle> hb = m_pMesh->boundaries[0];
			int num = hb.size();
			for (int i = 0; i < 4; i++) {
				MeshFlow::HalfedgeHandle h = hb[num / 4 * i];
				MeshFlow::VertexHandle v = m_pMesh->to_vertex_handle(h);
				m_pMesh->TargetCurvature(v) = PI / 2;
			}
		}

	}


	inline double EuclideanFlow::RToU(double r)
	{
		return log(r);
	}


	inline double EuclideanFlow::UToR(double u)
	{
		return exp(u);
	}


	inline void EuclideanFlow::ComputeMetricWeight(MeshFlow::EHandle e)
	{
		MeshFlow::VertexHandle v1 = m_pMesh->EdgeVertex1(e);
		MeshFlow::VertexHandle v2 = m_pMesh->EdgeVertex2(e);
		double l = m_pMesh->OriginalLength(e);
		double r1 = UToR(m_pMesh->U(v1));
		double r2 = UToR(m_pMesh->U(v2));
		double weight = (l * l - r1 * r1 - r2 * r2) / (2 * r1 * r2);
		//if (weight < -1) {
		//	weight = -1;
		//}
		//if (weight > 1) {
		//	weight = 1;
		//}
		m_pMesh->MetricWeight(e) = weight;
	}
	;



	inline void EuclideanFlow::Newton(double threshold, double stepLength)
	{
		//reset_U();
		int num = m_pMesh->n_vertices();
		int max_iteration = 100;
		int max_error_num = 10;
		std::list<double> errors;
		for (int i = 0; i < max_iteration; ++i)
		{
			//the order of the following functions really matters
			CalculateEdgeLength();
			CalculateCornerAngle();
			CalculateVertexCurvature();
			CalculateEdgeWeight();
			if (boundary_mode == Circle) SetCircleCurvature();
			//double area = calculate_total_area();
			double new_error = CalculateCurvatureError();
			printf("Newton's Method: Current error is %f\r\n", new_error);

			if (errors.size() < 10) {
				errors.push_back(new_error);
			}
			else {
				errors.push_back(new_error);
				errors.erase(errors.begin());
			}

			double max_error = *(std::max_element(errors.begin(), errors.end()));
			double min_error = *(std::min_element(errors.begin(), errors.end()));
			if (new_error < threshold) {
				return;
			}
			if (errors.size() == 10) {
				if (max_error - min_error < 1e-15) {
					return;
				}
			}

			Eigen::SparseMatrix<double>  M(num, num);
			M.setZero();
			CalculateHessianMatrix(M);

			Eigen::VectorXd b(num);
			for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
			{
				MeshFlow::VertexHandle v = *viter;
				int idx = m_pMesh->Idx(v);
				b(idx) = m_pMesh->TargetCurvature(v) - m_pMesh->Curvature(v);
			}

			//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
			Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
			solver.compute(M);

			if (solver.info() != Eigen::Success)
			{
				std::cerr << "Waring: Eigen decomposition failed" << std::endl;
				continue;
			}
			Eigen::VectorXd x = solver.solve(b);
			Normalization(x, num);
			AdaptiveUpdateNewton(x, 1);
		}
	}


	inline void EuclideanFlow::AcceleratedGradientDescent(double threshold, double stepLength)
	{
		//printf("stepLength: %f", stepLength);

		int num = m_pMesh->n_vertices();
		int max_iteration = 50000;

		std::list<double> errors;
		int max_error_num = 10;

		Eigen::VectorXd x(num);
		Eigen::VectorXd xb(num);
		Eigen::VectorXd y;
		double c;
		//Eigen::VectorXd U(num);

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VertexHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				x(idx) = 0;
				xb(idx) = 0;
				continue;
			}
			x(idx) = m_pMesh->U(v);
			xb(idx) = m_pMesh->U(v);
		}
		y = x;
		for (int i = 0; i < max_iteration; ++i)
		{
			//the order of the following functions really matters
			CalculateEdgeLength();
			CalculateCornerAngle();
			CalculateVertexCurvature();
			//CalculateEdgeWeight();
			//double area = calculate_total_area();
			double new_error = CalculateCurvatureError();
			if (i % 10 == 0)
				printf("Accelerated Gradient Method: Current error is %f\r\n", new_error);

			if (errors.size() < 10) {
				errors.push_back(new_error);
			}
			else {
				errors.push_back(new_error);
				errors.erase(errors.begin());
			}

			if (boundary_mode == Circle && i % 10 == 0) {
				SetCircleCurvature();
			}

			double max_error = *(std::max_element(errors.begin(), errors.end()));
			double min_error = *(std::min_element(errors.begin(), errors.end()));
			if (new_error < threshold) {
				return;
			}
			if (errors.size() > 2) {
				if (max_error - min_error < 1e-15) {
					return;
				}
			}

			if (i == 0) {
				x = UpdateAgRicci(y, 0.5, true);
			}
			else {
				xb = x;
				x = UpdateAgRicci(y, 0.5, true);
				Eigen::VectorXd diff = x - xb;
				y = x + double(i - 1) / double(i + 2) * diff;
			}

			if (boundary_mode == Free) {
				double c = 0;
				Eigen::VectorXd diff = x - xb;
				for (int i = 0; i < num; i++) {
					c += diff(i);
				}
				c = c / m_pMesh->n_vertices();
				for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
				{
					MeshFlow::VHandle v = *viter;
					//		//if (v->boundary()) continue;
					double u = m_pMesh->U(v);
					m_pMesh->U(v) = u - c;
				}
			}
		}
	}


	inline void EuclideanFlow::CalculateHessianMatrix(Eigen::SparseMatrix<double>& M)
	{
		std::vector<Eigen::Triplet<double> > M_coefficients;

		//set A
		for (MeshFlow::EdgeIter eiter = m_pMesh->edges_begin(); eiter != m_pMesh->edges_end(); eiter++)
		{
			MeshFlow::EHandle e = *eiter;
			MeshFlow::VHandle v1 = m_pMesh->EdgeVertex1(e);
			MeshFlow::VHandle v2 = m_pMesh->EdgeVertex2(e);
			if (boundary_mode == Free && (m_pMesh->is_boundary(v1) || m_pMesh->is_boundary(v2))) continue;
			int idx1 = m_pMesh->Idx(v1);
			int idx2 = m_pMesh->Idx(v2);
			M_coefficients.push_back(Eigen::Triplet<double>(idx1, idx2, -m_pMesh->Weight(e)));
			M_coefficients.push_back(Eigen::Triplet<double>(idx2, idx1, -m_pMesh->Weight(e)));

		}

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VertexHandle  v = *viter;
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) continue;
			int idx = m_pMesh->Idx(v);
			double w = 0;
			for (MeshFlow::VertexEdgeIter veiter = m_pMesh->ve_iter(v); veiter.is_valid(); veiter++)
			{
				MeshFlow::EHandle pE = *veiter;
				w += m_pMesh->Weight(pE);
			}
			M_coefficients.push_back(Eigen::Triplet<double>(idx, idx, w));
		}

		M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());

	}

}

#endif  _EUCLIDEAN_FLOW_H_