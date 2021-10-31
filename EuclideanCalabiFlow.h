#pragma once
#pragma once
/*! \file TangentialRicciFlow.h
*  \brief General Euclidean Ricci flow algorithm
*  \author David Gu
*  \date   documented on 10/17/2010
*
*	Algorithm for general Ricci Flow
*/

#ifndef _CALABI_EUC
#define _CALABI_EUC

#include <map>
#include <vector>
#include <Eigen/Sparse>
#include "EuclideanFlow.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif


namespace OpenMesh
{
	/*! \brief Class EuclideanCalabiFlow
	*
	*	Algorithm for computing Ricci flow
	*/

	class EuclideanCalabiFlow : public EuclideanFlow
	{
	public:
		/*! \brief EuclideanCalabiFlow constructor
		*  \param p_mesh the input mesh
		*/
		EuclideanCalabiFlow(MeshFlow* p_mesh);
		/*! \brief EuclideanCalabiFlow destructor
		*/
		~EuclideanCalabiFlow() {};
		/*!	Computing the metric
		*/
		void CalculateMetric();
		/*!
		*	Curvature flow, override
		*/
	protected:

		void AcceleratedGradientDescent(double threshold, double stepLength);

	};



	inline EuclideanCalabiFlow::EuclideanCalabiFlow(MeshFlow* p_mesh) :EuclideanFlow(p_mesh)
	{
	};



	//compute metric


	inline void EuclideanCalabiFlow::CalculateMetric()
	{
		double error = NEWTON_ERROR;
		if (m_pMesh->boundaries.size() == 0) {
			error = error * 0.1;
		}
		CalculateOriginalEdgeLength();
		//if (m_pMesh->freebound) reset_U();
		CalculateInitialMetric();
		CalculateEdgeLength();
		SetTargetCurvature();
		AcceleratedGradientDescent(error, 1);
	}

	inline void EuclideanCalabiFlow::AcceleratedGradientDescent(double threshold, double stepLength)
	{
		//printf("stepLength: %f", stepLength);
		int num = m_pMesh->n_vertices();
		int max_iteration = 20000;
		std::list<double> errors;
		int max_error_num = 10;

		Eigen::VectorXd x(num);
		Eigen::VectorXd xb(num);
		Eigen::VectorXd y;
		//Eigen::VectorXd U(num);
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			VertexHandle  v = *viter;
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
			CalculateEdgeWeight();
			if (boundary_mode == Circle && i % 10 == 0) {
				SetCircleCurvature();
			}
			double new_error = CalculateCurvatureError();
			if (i % 10 == 0) {
				printf("Accelerated Gradient Method: Current error is %f\r\n", new_error);
			}


			if (errors.size() < max_error_num) {
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
			if (errors.size() == max_error_num) {
				if (max_error - min_error < 1e-15) {
					return;
				}
			}

			if (i == 0) {
				x = UpdateAgCalabi(y, 0.5, true);
			}
			else {
				xb = x;
				x = UpdateAgCalabi(y, 0.5, true);
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
}


#endif  _CALABI_EUC