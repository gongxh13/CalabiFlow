#ifndef _UTILITIES_FLOW_H_
#define _UTILITIES_FLOW_H_

#include <Eigen\Dense>
#include "MeshFlow.h"

// Convert Triangle openmesh data structure to Eigen matrices
namespace OpenMesh {

	inline void LoadMeshDataToMatrix(OpenMesh::MeshFlow* Mesh, Eigen::MatrixXd& V, Eigen::MatrixXd& VN, Eigen::MatrixXi& F, Eigen::MatrixXd& FN)
	{
		int nv = Mesh->n_vertices();
		int nf = Mesh->n_faces();
		V.resize(nv, 3);
		VN.resize(nv, 3);
		F.resize(nf, 3);
		FN.resize(nf, 3);

		Mesh->update_normals();

		/*Reindex the idx*/
		int idx = 0;
		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter) {
			OpenMesh::VertexHandle v = *viter;
			Mesh->Idx(v) = idx;
			idx++;
		}

		/*load vertex data*/
		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter) {
			OpenMesh::VertexHandle v = *viter;

			int idx = Mesh->Idx(v);
			OpenMesh::Vec3d point = Mesh->point(v);
			OpenMesh::Vec3d normal = Mesh->normal(v);

			for (int i = 0; i < 3; i++) {
				V(idx, i) = point[i];
				VN(idx, i) = normal[i];
			}
		}

		int fidx = 0;
		for (OpenMesh::MeshFlow::FaceIter fiter = Mesh->faces_begin(); fiter != Mesh->faces_end(); ++fiter) {
			OpenMesh::FaceHandle f = *fiter;
			int i = 0;
			for (OpenMesh::MeshFlow::FaceVertexIter fviter = Mesh->fv_iter(f); fviter.is_valid(); ++fviter) {
				OpenMesh::VertexHandle v = *fviter;
				F(fidx, i) = Mesh->Idx(v);
				i++;
			}
			OpenMesh::Vec3d normal = Mesh->normal(f);
			for (int j = 0; j < 3; j++) {
				FN(fidx, j) = normal[j];
			}
			fidx++;
		}
	}

	inline void LoadMeshFlatMappingResultToMatrix(OpenMesh::MeshFlow* Mesh, Eigen::MatrixXd& UV)
	{
		int nv = Mesh->n_vertices();
		UV.resize(nv, 2);

		/*load vertex data*/
		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter) {
			OpenMesh::VertexHandle v = *viter;

			int idx = Mesh->Idx(v);
			OpenMesh::Vec2d uv = Mesh->texcoord2D(v);

			for (int i = 0; i < 2; i++) {
				UV(idx, i) = uv[i];
			}
		}
	}

	inline void LoadMeshSphericalMappingResultToMatrix(OpenMesh::MeshFlow* Mesh, Eigen::MatrixXd& V)
	{
		int nv = Mesh->n_vertices();
		V.resize(nv, 3);

		/*load vertex data*/
		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter) {
			OpenMesh::VertexHandle v = *viter;

			int idx = Mesh->Idx(v);
			OpenMesh::Vec3d point = Mesh->texcoord3D(v);
			for (int i = 0; i < 3; i++) {
				V(idx, i) = point[i];
			}
		}
	}

	inline void NormalizeMesh(OpenMesh::MeshFlow* Mesh)
	{
		OpenMesh::MeshFlow::Point s(0, 0, 0);
		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter)
		{
			OpenMesh::VertexHandle  v = *viter;
			s = s + Mesh->point(v);
		}
		s = s / Mesh->n_vertices();

		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter)
		{
			OpenMesh::VertexHandle v = *viter;
			OpenMesh::MeshFlow::Point p = Mesh->point(v);
			p = p - s;
			Mesh->set_point(v, p);
		}

		double d = 0;
		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter)
		{
			OpenMesh::VertexHandle v = *viter;
			OpenMesh::MeshFlow::Point p = Mesh->point(v);
			for (int k = 0; k < 3; k++)
			{
				d = (d > fabs(p[k])) ? d : fabs(p[k]);
			}
		}

		for (OpenMesh::MeshFlow::VertexIter viter = Mesh->vertices_begin(); viter != Mesh->vertices_end(); ++viter)
		{
			OpenMesh::VertexHandle v = *viter;
			OpenMesh::MeshFlow::Point p = Mesh->point(v);
			p = p / d;
			Mesh->set_point(v, p);
		}

	}
}


#endif // !_Mesh_IO_H_

#pragma once
