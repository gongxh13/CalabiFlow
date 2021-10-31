
#ifndef _BASE_FLOW_H_
#define _BASE_FLOW_H_

#include <map>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "MeshFlow.h"
#include "Meshlib/core/Geometry/Circle.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif


namespace OpenMesh
{

	enum BoundaryMode { Fixed, Circle, Free };
	/*! \brief BaseClass BaseFlow
	*
	*	Algorithm for computing general Ricci flow
	*/

	class BaseFlow
	{
	public:
		/*! \brief BaseFlow constructor
		 *  \param pMesh the input mesh
		 */
		BaseFlow(MeshFlow* pMesh);
		/*! \brief BaseFlow destructor
		 */
		~BaseFlow();
		/*!	Computing the metric
		 */
		virtual void CalculateMetric() = 0;

		/*Flags*/
		/*Fixed: fix boundary cuvature (if no boundary, set all target curvarue zero)
		  Circle: map to circle,
		  Free: free boundary setting*/
		BoundaryMode boundary_mode;

		// whether use orbifold mode
		char orbifold;

	protected:
		/*!
		 *	the input mesh
		 */
		MeshFlow* m_pMesh;

		/* Adaptive step Length*/
		void AdaptiveUpdateNewton(Eigen::VectorXd& x, double beta);
		Eigen::VectorXd UpdateAgCalabi(Eigen::VectorXd y, double beta, bool normalize);
		Eigen::VectorXd UpdateAgRicci(Eigen::VectorXd y, double beta, bool normalize);
		/* for Calabi Flow */
		Eigen::VectorXd CalculateGradient();

		/*!
		 *	Calculate each edge Length, has to be defined in the derivated classes
		 */
		virtual void Length(double u1, double u2, EdgeHandle e) = 0;
		/*!
		 *	Calculate each edge Length
		 */
		void CalculateEdgeLength();

		/*!
		 *	Cosine law, has to be defined in the derivated classes
		 */
		virtual double CosineLaw(double a, double b, double c) = 0;

		/*!
		 *	Calculate corner angle
		 */
		virtual void CalculateCornerAngle();

		/*!
		 *	Calculate vertex curvature
		 */
		void CalculateVertexCurvature();

		/*calcualte original edge Length in R^3*/
		void CalculateOriginalEdgeLength();

		/*calculate initial circle metric*/
		virtual void CalculateInitialMetric();

		/*!
		 *	Calculate vertex curvature error
		 */
		double   CalculateCurvatureError();

		/*!
		 *	Calculate the edge weight
		 */
		virtual void CalculateEdgeWeight() = 0;

		/*!
		 *	Set the target curvature on each vertex
		 */
		virtual void  SetTargetCurvature() = 0;

		/*!
		 *	Newton's method to optimize the entropy energy
		 * \param threshold err bound
		 * \param stepLength step Length
		 */
		virtual void   Newton(double threshold, double stepLength) = 0;
		/*!
		 *	Normalization
		 * \param du the du vector
		 * \param n dimension of the du vector
		 */
		virtual void Normalization(Eigen::VectorXd& du, int n);
		/*!
		 *	calculate hessian matrix Hessain
		 * \param SparseMatrix
		 */
		virtual void CalculateHessianMatrix(Eigen::SparseMatrix<double>& pMatrix) = 0;

		virtual double RToU(double r) = 0;

		virtual double UToR(double u) = 0;

		virtual void ComputeMetricWeight(EdgeHandle e) = 0;

		double SetCircleCurvature();

	};


	//Constructor

	inline BaseFlow::BaseFlow(MeshFlow* pMesh) : m_pMesh(pMesh)
	{
		pMesh->RequestBoundary();
		orbifold = false;
		boundary_mode = Fixed;

		int idx = 0;
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VertexHandle v = *viter;
			m_pMesh->Idx(v) = idx;
			idx++;
		}
	};


	//Destructor

	inline BaseFlow::~BaseFlow()
	{
	};


	inline void BaseFlow::AdaptiveUpdateNewton(Eigen::VectorXd& direction, double beta)
	{
		int num = m_pMesh->n_vertices();
		double init_energy = CalculateCurvatureError();
		if (init_energy < 1e-5) {
			return;
		}
		Eigen::VectorXd ori_U(num);

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			ori_U(idx) = m_pMesh->U(v);
		}

		bool stop = true;
		while (true) {

			for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
			{
				MeshFlow::VertexHandle v = *viter;
				int idx = m_pMesh->Idx(v);
				m_pMesh->U(v) = ori_U(idx) + beta * direction(idx);
			}

			CalculateEdgeLength();
			CalculateCornerAngle();
			CalculateVertexCurvature();

			double new_energy = CalculateCurvatureError();
			if (new_energy > init_energy) {
				beta = beta * 0.5;
				if (beta < 0.0001) {
					for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
					{
						MeshFlow::VertexHandle v = *viter;
						int idx = m_pMesh->Idx(v);
						m_pMesh->U(v) = ori_U(idx);
					}
					break;
				}
				continue;
			}
			break;
		}
	}

	inline Eigen::VectorXd BaseFlow::UpdateAgCalabi(Eigen::VectorXd y, double beta, bool normalize)
	{
		int num = m_pMesh->n_vertices();
		Eigen::VectorXd ori_U(num);
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			ori_U(idx) = m_pMesh->U(v);
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				y(idx) = 0;
				continue;
			}
			m_pMesh->U(v) = y(idx);
		}
		CalculateEdgeLength();
		CalculateCornerAngle();
		CalculateVertexCurvature();
		CalculateEdgeWeight();
		Eigen::VectorXd init_grad = CalculateGradient();
		if (normalize) {
			Normalization(init_grad, num);
		}
		bool stop = true;
		Eigen::VectorXd xm = y + 0.01 * init_grad;
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VertexHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				xm(idx) = 0;
				continue;
			}
			m_pMesh->U(v) = xm(idx);
		}
		return xm;
	}

	inline Eigen::VectorXd BaseFlow::UpdateAgRicci(Eigen::VectorXd y, double beta, bool normalize)
	{
		int num = m_pMesh->n_vertices();
		Eigen::VectorXd ori_U(num);
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			ori_U(idx) = m_pMesh->U(v);
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				y(idx) = 0;
				continue;
			}
			m_pMesh->U(v) = y(idx);
		}
		CalculateEdgeLength();
		CalculateCornerAngle();
		CalculateVertexCurvature();
		Eigen::VectorXd init_grad(num);
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				init_grad(idx) = 0;
			}
			else {
				init_grad(idx) = m_pMesh->TargetCurvature(v) - m_pMesh->Curvature(v);
			}
		}
		if (normalize) {
			Normalization(init_grad, num);
		}
		bool stop = true;
		Eigen::VectorXd xm = y + 0.01 * init_grad;
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VertexHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				xm(idx) = 0;
				continue;
			}
			m_pMesh->U(v) = xm(idx);
		}
		return xm;
	}

	inline Eigen::VectorXd BaseFlow::CalculateGradient()
	{
		int num = m_pMesh->n_vertices();
		Eigen::SparseMatrix<double>  M(num, num);
		M.setZero();
		CalculateHessianMatrix(M);
		Eigen::VectorXd b(num);
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VertexHandle v = *viter;
			int idx = m_pMesh->Idx(v);
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				b(idx) = 0;
			}
			else {
				b(idx) = m_pMesh->TargetCurvature(v) - m_pMesh->Curvature(v);
			}
		}
		Eigen::VectorXd grad = M.transpose() * b;
		if (boundary_mode == Free) {
			for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
			{
				MeshFlow::VertexHandle v = *viter;
				int idx = m_pMesh->Idx(v);
				if (m_pMesh->is_boundary(v)) {
					grad(idx) = 0;
				}
			}
		}
		return grad;
	}

	//Compute the edge Length

	inline void BaseFlow::CalculateEdgeLength()
	{

		for (MeshFlow::EdgeIter eiter = m_pMesh->edges_begin(); eiter != m_pMesh->edges_end(); eiter++)
		{
			MeshFlow::EHandle e = *eiter;

			MeshFlow::VertexHandle v1 = m_pMesh->EdgeVertex1(e);
			MeshFlow::VertexHandle v2 = m_pMesh->EdgeVertex2(e);

			double u1 = m_pMesh->U(v1);
			double u2 = m_pMesh->U(v2);
			Length(u1, u2, e);
		}

	};


	//Calculate corner angle

	inline void BaseFlow::CalculateCornerAngle()
	{

		for (MeshFlow::FaceIter fiter = m_pMesh->faces_begin(); fiter != m_pMesh->faces_end(); fiter++)
		{
			MeshFlow::FaceHandle f = *fiter;

			MeshFlow::HalfedgeHandle he[3];

			he[0] = m_pMesh->halfedge_handle(f);

			for (int i = 0; i < 3; i++)
			{
				he[(i + 1) % 3] = m_pMesh->next_halfedge_handle(he[i]);
			}

			double l[3];
			for (int i = 0; i < 3; i++)
			{
				MeshFlow::EdgeHandle e = m_pMesh->edge_handle(he[i]);
				l[i] = m_pMesh->Length(e);
			}

			for (int i = 0; i < 3; i++)
			{
				m_pMesh->Angle(he[(i + 1) % 3]) = CosineLaw(l[(i + 1) % 3], l[(i + 2) % 3], l[i]);
			}
		}

	};

	//Calculate vertex curvature

	inline void BaseFlow::CalculateVertexCurvature()
	{

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VertexHandle v = *viter;
			double k;
			if (orbifold) {
				k = PI * 2 / m_pMesh->Order(v);
			}
			else {
				k = (m_pMesh->is_boundary(v)) ? PI : PI * 2;
			}

			for (MeshFlow::VertexIHalfedgeIter vh = m_pMesh->vih_iter(v); vh.is_valid(); ++vh)
			{
				MeshFlow::HalfedgeHandle he = *vh;
				k -= m_pMesh->Angle(he);
			}
			m_pMesh->Curvature(v) = k;

		}

	}

	inline void BaseFlow::CalculateOriginalEdgeLength()
	{
		for (MeshFlow::EdgeIter eiter = m_pMesh->edges_begin(); eiter != m_pMesh->edges_end(); ++eiter) {
			MeshFlow::EdgeHandle e = *eiter;
			MeshFlow::VertexHandle v1 = m_pMesh->EdgeVertex1(e);
			MeshFlow::VertexHandle v2 = m_pMesh->EdgeVertex2(e);
			double Length = (m_pMesh->point(v1) - m_pMesh->point(v2)).norm();
			m_pMesh->OriginalLength(e) = Length;
		}
	}

	inline void BaseFlow::CalculateInitialMetric()
	{
		CalculateOriginalEdgeLength();
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); ++viter) {
			MeshFlow::VertexHandle v = *viter;
			int m = 0;
			double r = 1e30;
			//double r = 0.0;
			for (MeshFlow::VertexFaceIter vfiter = m_pMesh->vf_iter(v); vfiter.is_valid(); ++vfiter) {
				MeshFlow::FaceHandle f = *vfiter;
				++m;
				double tem_r = 0.0;
				for (MeshFlow::FaceEdgeIter feiter = m_pMesh->fe_iter(f); feiter.is_valid(); ++feiter) {
					MeshFlow::EdgeHandle e = *feiter;
					if (m_pMesh->EdgeVertex1(e) != v && m_pMesh->EdgeVertex2(e) != v) {
						tem_r -= m_pMesh->OriginalLength(e);
					}
					else {
						tem_r += m_pMesh->OriginalLength(e);
					}
				}
				r = tem_r / 2.0 < r ? tem_r / 2.0 : r;
				//r += tem_r/2.0;
			}
			//r = r / m;
			m_pMesh->U(v) = RToU(r);
		}

		for (MeshFlow::EdgeIter eiter = m_pMesh->edges_begin(); eiter != m_pMesh->edges_end(); ++eiter) {
			MeshFlow::EdgeHandle e = *eiter;
			ComputeMetricWeight(e);
		}

	}
	;


	//compute curvature error


	inline double BaseFlow::CalculateCurvatureError()
	{
		double max_error = -1;

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			MeshFlow::VHandle v = *viter;
			if (boundary_mode == Free && m_pMesh->is_boundary(v)) {
				continue;
			}
			double k = m_pMesh->TargetCurvature(v) - m_pMesh->Curvature(v);
			k = fabs(k);
			if (k > max_error)
			{
				max_error = k;
			}
		}
		//printf("Vertex id is %d\n", vert->Id() );
		return max_error;
	};

	//compute metric



	inline void BaseFlow::Normalization(Eigen::VectorXd& x, int num)
	{
		double s = 0;
		for (int i = 0; i < num; i++)
		{
			s += x(i);
		}
		s /= num;

		for (int i = 0; i < num; i++)
		{
			x(i) -= s;
		}
	}

	inline double BaseFlow::SetCircleCurvature()
	{
		double max_dif = 0;
		std::vector<std::vector<MeshFlow::HHandle>> pLs = m_pMesh->boundaries;
		int id = 0;
		for (std::vector<std::vector<MeshFlow::HHandle>>::iterator liter = pLs.begin(); liter != pLs.end(); liter++)
		{
			std::vector<MeshFlow::HHandle> pHes = *liter;

			double sum = 0;
			double inv_sum = 0;

			for (std::vector<MeshFlow::HHandle>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter++)
			{
				MeshFlow::HHandle  he = *hiter;
				MeshFlow::EHandle  pe = m_pMesh->edge_handle(he);
				double Length = m_pMesh->Length(pe);
				sum += Length;
				inv_sum += 1 / Length;
			}

			for (std::vector<MeshFlow::HHandle>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter++)
			{
				MeshFlow::HHandle ce = *hiter;
				MeshFlow::VHandle pv = m_pMesh->to_vertex_handle(ce);
				MeshFlow::HHandle he = ce;
				MeshFlow::HHandle te = m_pMesh->next_halfedge_handle(he);

				double L = (m_pMesh->Length(m_pMesh->edge_handle(he))
					+ m_pMesh->Length(m_pMesh->edge_handle(te))) / 2.0;

				// map all boundaries to circular holes
				if (id == 0)
				{
					double tk = 2 * PI * L / sum;
					double ndif = abs(m_pMesh->TargetCurvature(pv) - tk);
					ndif = ndif > 100 ? 100 : ndif;
					max_dif = ndif > max_dif ? ndif : max_dif;
					m_pMesh->TargetCurvature(pv) = tk;
				}
				else
				{
					double tk = -2 * PI * L / sum;
					double ndif = abs(m_pMesh->TargetCurvature(pv) - tk);
					ndif = ndif > 100 ? 100 : ndif;
					max_dif = ndif > max_dif ? ndif : max_dif;
					m_pMesh->TargetCurvature(pv) = tk;
				}
			}
			id++;
		}
		return max_dif;
	}

}
#endif  #pragma once
