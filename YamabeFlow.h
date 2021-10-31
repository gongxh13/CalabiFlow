#pragma once
/*! \file TangentialRicciFlow.h
 *  \brief General Euclidean Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _EQUIVALENCE_RICCI_FLOW_H_
#define _EQUIVALENCE_RICCI_FLOW_H_

#include <map>
#include <vector>
#include <Eigen/Sparse>
#include "EuclideanFlow.h"

namespace OpenMesh
{
	/*! \brief Class CYamabeFlow
	*
	*	Algorithm for computing Ricci flow
	*/

	class CYamabeFlow : public EuclideanFlow
	{
	public:
		/*! \brief CYamabeFlow constructor
		 *  \param p_mesh the input mesh
		 */
		CYamabeFlow(MeshFlow* p_mesh);
		/*! \brief CYamabeFlow destructor
		 */
		~CYamabeFlow() {};
		/*!	Computing the metric
		 */
		void CalculateMetric();

	protected:

		/*!
		 *	Calculate each edge Length, has to be defined in the derivated classes
		 */
		void Length(double u1, double u2, EdgeHandle e);


		//void Normalization();
		/*!
		 *	Calculate the edge weight
		 */
		void CalculateEdgeWeight();

		void CalculateInitialMetric();


	};



	inline CYamabeFlow::CYamabeFlow(MeshFlow* p_mesh) :EuclideanFlow(p_mesh)
	{
	};

	//Compute the edge Length

	inline void CYamabeFlow::Length(double u1, double u2, EdgeHandle  e)
	{
		m_pMesh->Length(e) = exp(u1) * exp(u2) * m_pMesh->OriginalLength(e);
	};


	//Calculate edge weight


	inline void CYamabeFlow::CalculateEdgeWeight()
	{
		for (MeshFlow::EdgeIter eiter = m_pMesh->edges_begin(); eiter != m_pMesh->edges_end(); eiter++) {
			double weight = 0.0;
			EdgeHandle e = *eiter;
			HalfedgeHandle h0 = m_pMesh->EdgeHalfedge(e, 0);
			HalfedgeHandle h1 = m_pMesh->EdgeHalfedge(e, 1);
			if (!m_pMesh->is_boundary(h0)) {
				double angle = m_pMesh->Angle(m_pMesh->next_halfedge_handle(h0));
				if (abs(angle) > 1e-5 && abs(angle - PI) > 1e-5) {
					double cotan = 1 / tan(angle);
					if (!isnan(cotan))
						weight += cotan;
				}
			}
			if (!m_pMesh->is_boundary(h1)) {
				double angle = m_pMesh->Angle(m_pMesh->next_halfedge_handle(h1));
				if (abs(angle) > 1e-5 && abs(angle - PI) > 1e-5) {
					double cotan = 1 / tan(angle);
					weight += cotan;
				}
			}
			m_pMesh->Weight(e) = weight;
		}
	}


	inline void CYamabeFlow::CalculateInitialMetric()
	{
		CalculateOriginalEdgeLength();
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); ++viter) {
			VertexHandle v = *viter;
			m_pMesh->U(v) = 0.0;
		}
	}



	inline void CYamabeFlow::CalculateMetric()
	{
		double error = NEWTON_ERROR;
		if (m_pMesh->boundaries.size() == 0) {
			error = error * 0.1;
		}
		CalculateInitialMetric();
		SetTargetCurvature();
		CalculateEdgeLength();
		if (boundary_mode == Circle) {
			double result = SetCircleCurvature();
		}
		if (boundary_mode != Free) {
			Newton(error, 1.0);
		}
		AcceleratedGradientDescent(error, 1.0);
	}






}

#endif  _EQUIVALENCE_RICCI_FLOW_H_