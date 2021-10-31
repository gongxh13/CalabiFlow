#ifndef  _BASE_EMBED_H_
#define  _BASE_EMBED_H_

#include <vector>
#include <queue>

#include "MeshFlow.h"

namespace OpenMesh
{
	/*! \brief CBaseEmbed class
	 *
	 *   embed a mesh with canonical metric onto the canonical domain
	 */

	class CBaseEmbed
	{
	public:
		/*! \brief CEmbed constructor */
		CBaseEmbed(MeshFlow* pMesh);
		/*! \brief CEmbed destructor */
		~CBaseEmbed() {};
		/*! _embed the mesh */
		void embed();


	protected:

		/*! embed the first face
		 * \param head the first face
		 */
		virtual void EmbedFirstFace(FaceHandle head) = 0;
		/*!
		 *
		 */
		virtual bool EmbedFace(FaceHandle head) = 0;

		/*choose root face*/
		FaceHandle InitRootFace();

		/*! initialization */
		void Initialize();

		/*normalize uv*/
		virtual void NormalizeCoord2D() {}

		template <class T1, class T2>
		T2 Convert(T1 t);


	protected:

		/*! the mesh to be embedded */
		MeshFlow* m_pMesh;

		/*! queue of faces */
		std::queue<FaceHandle> m_queue;
	};


	//constructor

	inline CBaseEmbed::CBaseEmbed(MeshFlow* pMesh)
	{
		m_pMesh = pMesh;
	};


	inline FaceHandle CBaseEmbed::InitRootFace()
	{
		/*******************************\
		find the closest face, the eye position is
		set as (0,0,5), so the distance is z
		\*******************************/

		double z = 0;
		FaceHandle root_face;
		for (MeshFlow::FaceIter fiter = m_pMesh->faces_begin(); fiter != m_pMesh->faces_end(); ++fiter) {
			FaceHandle f = *fiter;
			m_pMesh->Touched(f) = false;
			MeshFlow::Point s(0, 0, 0);
			for (MeshFlow::FaceVertexIter fviter = m_pMesh->fv_iter(f); fviter.is_valid(); ++fviter) {
				VertexHandle v = *fviter;
				s += m_pMesh->point(v);
			}
			s /= 3;
			if (s[2] > z) {
				z = s[2];
				root_face = f;
			}
		}
		return root_face;
	}

	inline void CBaseEmbed::Initialize()
	{
		for (MeshFlow::VertexIter viter = m_pMesh->vertices_begin(); viter != m_pMesh->vertices_end(); viter++) {
			VertexHandle v = *viter;
			m_pMesh->Touched(v) = false;
			m_pMesh->set_texcoord2D(v, Vec2d(0, 0));
		}

		FaceHandle root_face = InitRootFace();
		EmbedFirstFace(root_face);
		m_pMesh->Touched(root_face) = true;

		for (MeshFlow::FaceFaceIter ffiter = m_pMesh->ff_iter(root_face); ffiter.is_valid(); ++ffiter) {
			FaceHandle df = *ffiter;
			m_queue.push(df);
			m_pMesh->Touched(df) = true;
		}

	}


	inline void CBaseEmbed::embed()
	{

		Initialize();
		while (!m_queue.empty())
		{
			FaceHandle head = m_queue.front();
			m_queue.pop();
			assert(m_pMesh->Touched(head));

			//CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge(head);
			for (MeshFlow::FaceFaceIter ffiter = m_pMesh->ff_iter(head); ffiter.is_valid(); ++ffiter) {
				FaceHandle df = *ffiter;
				if (!(m_pMesh->Touched(df)))
				{
					m_queue.push(df);
					m_pMesh->Touched(df) = true;
				}
			}

			bool result = EmbedFace(head);
			if (!result) {
				return;
			}
		}

		NormalizeCoord2D();

	};

	template<class T1, class T2>
	inline T2 CBaseEmbed::Convert(T1 t)
	{
		return T2(t[0], t[1]);
	}

}


#endif #pragma once
