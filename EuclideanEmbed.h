/*! \file EuclideanEmbed.h
 *  \brief Isometrically embed a mesh with flat metric onto the plane
 *  \author David Gu
 *  \date  documented on 10/17/2010
 *
 *  Embed a mesh with flat metric on canonical domain
 */
#ifndef  _EUCLIDEAN_EMBED_H_
#define  _EUCLIDEAN_EMBED_H_

#include <vector>
#include <queue>

#include "BaseEmbed.h"


namespace OpenMesh
{
	/*! \brief CEmbed class
	 *
	 *   embed a mesh with flat metric onto the plane
	 */

	class EuclideanEmbed : public CBaseEmbed
	{
	public:
		/*! \brief EuclideanEmbed constructor */
		EuclideanEmbed(MeshFlow* pMesh);
		/*! \brief EuclideanEmbed destructor */
		~EuclideanEmbed() {};

		// store the generators of deck transition
		//std::vector<EuclideanTransformation> transition_generators;

	protected:
		/*! embed the first face
		 * \param head the first face
		 */
		void EmbedFirstFace(FaceHandle head);
		/*!
		 *	embed one face
		 */
		bool EmbedFace(FaceHandle head);

		void NormalizeCoord2D();

	};


	//constructor

	inline EuclideanEmbed::EuclideanEmbed(MeshFlow* pMesh) :CBaseEmbed(pMesh)
	{
	}

	//embed the first face

	inline void EuclideanEmbed::EmbedFirstFace(FaceHandle head)
	{
		HalfedgeHandle he[3];
		EdgeHandle edges[3];
		he[0] = m_pMesh->halfedge_handle(head);
		he[1] = m_pMesh->next_halfedge_handle(he[0]);
		he[2] = m_pMesh->next_halfedge_handle(he[1]);

		for (int i = 0; i < 3; i++) {
			edges[i] = m_pMesh->edge_handle(he[i]);
		}

		std::vector<VertexHandle> av;

		for (int i = 0; i < 3; i++)
		{
			av.push_back(m_pMesh->to_vertex_handle(he[(i + 2) % 3]));
		}
		MeshLib::CPoint2 uv_p[3];
		Vec2d uv_v[3];
		EdgeHandle e = m_pMesh->edge_handle(he[0]);
		uv_v[0] = Vec2d(0, 0);
		uv_v[1] = Vec2d(m_pMesh->Length(e), 0);
		uv_p[0] = Convert<Vec2d, MeshLib::CPoint2>(uv_v[0]);
		uv_p[1] = Convert<Vec2d, MeshLib::CPoint2>(uv_v[1]);
		m_pMesh->set_texcoord2D(av[0], uv_v[0]);
		m_pMesh->set_texcoord2D(av[1], uv_v[1]);

		MeshLib::CPoint2 c1, c2;

		_circle_circle_intersection(MeshLib::CCircle(uv_p[0], m_pMesh->Length(edges[2])),
			MeshLib::CCircle(uv_p[1], m_pMesh->Length(edges[1])),
			c1, c2);

		if (MeshLib::cross(uv_p[1] - uv_p[0], c1 - uv_p[0]) > 0)
		{
			uv_p[2] = c1;
		}
		else
		{
			uv_p[2] = c2;
		}

		m_pMesh->set_texcoord2D(av[2], Convert<MeshLib::CPoint2, Vec2d>(uv_p[2]));

		for (int i = 0; i < 3; i++)
		{
			m_pMesh->Touched(m_pMesh->to_vertex_handle(he[i])) = true;
		}

	};


	//vertex A, vertex B are known, vertex C is unknown


	inline bool EuclideanEmbed::EmbedFace(FaceHandle head)
	{
		std::vector<VertexHandle> av;

		for (MeshFlow::FaceVertexIter fviter = m_pMesh->fv_iter(head); fviter.is_valid(); fviter++)
		{
			VertexHandle pV = *fviter;
			av.push_back(pV);
		}

		VertexHandle A;
		A.invalidate();
		VertexHandle B;
		B.invalidate();
		VertexHandle C;
		C.invalidate();

		for (int i = 0; i < 3; i++)
		{
			if (m_pMesh->Touched(av[i])) continue;
			C = av[(i + 0) % 3];
			A = av[(i + 1) % 3];
			B = av[(i + 2) % 3];
			break;
		}

		if (!C.is_valid()) return true;

		//radius of the first circle
		double r1 = m_pMesh->Length(m_pMesh->VertexEdge(A, C));
		//radius of the second circle
		double r2 = m_pMesh->Length(m_pMesh->VertexEdge(B, C));
		//center of the first circle
		Vec2d c1 = m_pMesh->texcoord2D(A);
		//center of the second circle
		Vec2d c2 = m_pMesh->texcoord2D(B);

		MeshLib::CPoint2 c1_p = Convert<Vec2d, MeshLib::CPoint2>(c1);
		MeshLib::CPoint2 c2_p = Convert<Vec2d, MeshLib::CPoint2>(c2);

		MeshLib::CPoint2 i1;
		MeshLib::CPoint2 i2;

		int result = _circle_circle_intersection(MeshLib::CCircle(c1_p, r1), MeshLib::CCircle(c2_p, r2), i1, i2);

		Vec2d i1_v = Convert<MeshLib::CPoint2, Vec2d>(i1);
		Vec2d i2_v = Convert<MeshLib::CPoint2, Vec2d>(i2);

		if (result == 1)
		{
			if (MeshLib::cross(c2_p - c1_p, i1 - c1_p) > 0)
				m_pMesh->set_texcoord2D(C, i1_v);
			else
				m_pMesh->set_texcoord2D(C, i2_v);

			m_pMesh->Touched(C) = true;
			return true;
		}

		else
		{
			return false;
		}
	}

	inline void EuclideanEmbed::NormalizeCoord2D()
	{
		double u_min, u_max, v_min, v_max;
		u_min = v_min = 1e30;
		u_max = v_max = -1e30;

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			VertexHandle pV = *viter;
			Vec2d uv = m_pMesh->texcoord2D(pV);
			double u_param = uv[0];
			double v_param = uv[1];
			if (u_min > u_param)
				u_min = u_param;
			if (u_max < u_param)
				u_max = u_param;
			if (v_min > v_param)
				v_min = v_param;
			if (v_max < v_param)
				v_max = v_param;
		}

		double range = (u_max - u_min) > (v_max - v_min) ? u_max - u_min : v_max - v_min;

		printf("range = %lf\n", range);

		Vec2d s(0, 0);
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			VertexHandle pV = *viter;
			s += m_pMesh->texcoord2D(pV);
		}
		s = s / m_pMesh->n_vertices();
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			VertexHandle pV = *viter;
			Vec2d uv = m_pMesh->texcoord2D(pV);
			m_pMesh->set_texcoord2D(pV, uv - s);
		}


		double max_dist = 0;

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			VertexHandle pV = *viter;
			Vec2d uv = m_pMesh->texcoord2D(pV);
			max_dist = uv.norm() > max_dist ? uv.norm() : max_dist;
		}

		range = max_dist;

		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++)
		{
			VertexHandle pV = *viter;
			Vec2d uv = m_pMesh->texcoord2D(pV);
			m_pMesh->set_texcoord2D(pV, uv / range);
		}

	}
} // end of namespace OpenMesh


#endif #pragma once
