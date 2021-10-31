#include "BaseFlow.h"
#include "YamabeFlow.h"
#include "EuclideanEmbed.h"
#include "UtilitiesFlow.h"
#include "EuclideanCalabiFlow.h"
#ifndef _FLOW_CONTROLER_H_
#define _FLOW_CONTROLER_H_


namespace OpenMesh {

	enum FlowType { EuclideanCalabi, EuclideanRicci, CETM, HyperbolicCalabi, HyperbolicRicci, SphericalCalabi, SphericalRicci };

	class  ControlerFlow
	{
	public:
		ControlerFlow() { m_mesh = NULL; boundary_mode = Fixed; orbifold = false;  flow_type = EuclideanRicci; }
		~ControlerFlow() { delete m_mesh; }

		/*IO*/
		bool LoadMeshFromFile(const char* input);
		bool SaveMeshToFile(const char* input);
		void MeshConvertToMatrix(Eigen::MatrixXd& V, Eigen::MatrixXd& VN, Eigen::MatrixXd& UV, Eigen::MatrixXd& UV3, Eigen::MatrixXi& F, Eigen::MatrixXd& FN);

		/*algorithms*/

		void RunFlow();




		//Flags, thses flags can be access outside the class
		FlowType flow_type;
		BoundaryMode boundary_mode;
		bool orbifold;
		int orbifold_mode;


	protected:
		MeshFlow* m_mesh;

		inline  void RunEuclideanRicciFlow();
		inline  void RunCETM();
		inline  void RunEuclideanCalabiFlow();
		inline  void RunHyperbolicRicciFlow();
		inline  void RunHyperbolicCalabiFlow();
		inline  void RunSphericalRicciFlow();
		inline  void RunSphericalCalabiFlow();

	};

	void ControlerFlow::RunCETM()
	{
		// flow
		CYamabeFlow CETM(m_mesh);
		switch (boundary_mode) {
		case 0:
			CETM.boundary_mode = Fixed;
			break;
		case 1:
			CETM.boundary_mode = Circle;
			break;
		case 2:
			CETM.boundary_mode = Free;
		}
		CETM.CalculateMetric();

		m_mesh->RequestBoundary();
		if (m_mesh->boundaries.size() == 0) {
			//CMeshReconstructor Slice(m_mesh);
			//MeshFlow* new_mesh = Slice.SliceToDisk();
			//delete m_mesh;
			//m_mesh = new_mesh;
		}

		//embed
		EuclideanEmbed euc_ebd(m_mesh);
		euc_ebd.embed();
	}


	inline bool ControlerFlow::LoadMeshFromFile(const char* input)
	{
		if (m_mesh) {
			delete m_mesh;
		}
		m_mesh = new OpenMesh::MeshFlow;

		/*Check extension of file*/
		std::string mesh_file_name_string = std::string(input);
		size_t last_dot = mesh_file_name_string.rfind('.');
		if (last_dot == std::string::npos)
		{
			printf("Error: No file extension found in %s\n", input);
			return false;
		}
		std::string extension = mesh_file_name_string.substr(last_dot + 1);

		/*Choose load function according to extension*/
		if (extension == "obj" || extension == "OBJ")
		{
			m_mesh->ReadObjFile(input);
		}
		if (extension == "m" || extension == "M")
		{
			m_mesh->ReadMFile(input);
		}
		NormalizeMesh(m_mesh);
	}

	inline bool ControlerFlow::SaveMeshToFile(const char* output)
	{
		if (!m_mesh) {
			return false;
		}
		/*Check extension of file*/
		std::string mesh_file_name_string = std::string(output);
		size_t last_dot = mesh_file_name_string.rfind('.');
		if (last_dot == std::string::npos)
		{
			printf("Error: No file extension found in %s\n", output);
			return false;
		}
		std::string extension = mesh_file_name_string.substr(last_dot + 1);

		/*Choose load function according to extension*/
		if (extension == "obj" || extension == "OBJ")
		{
			m_mesh->WriteObjFile(output);
		}
		if (extension == "m" || extension == "M")
		{
			m_mesh->WriteMFile(output);
		}
	}

	inline void ControlerFlow::MeshConvertToMatrix(Eigen::MatrixXd& V, Eigen::MatrixXd& VN, Eigen::MatrixXd& UV, Eigen::MatrixXd& UV3, Eigen::MatrixXi& F, Eigen::MatrixXd& FN)
	{
		LoadMeshDataToMatrix(m_mesh, V, VN, F, FN);
		LoadMeshFlatMappingResultToMatrix(m_mesh, UV);
		LoadMeshSphericalMappingResultToMatrix(m_mesh, UV3);
	}


	inline void ControlerFlow::RunFlow()
	{
		if (!m_mesh) {
			return;
		}
		switch (flow_type)
		{
		case OpenMesh::EuclideanCalabi:
			if (!orbifold) RunEuclideanCalabiFlow();
			//else RunEuclideanCalabiFlowWithOrbifold();
			break;
		//case OpenMesh::EuclideanRicci:
		//	if (!orbifold) RunEuclideanRicciFlow();
		//	else RunEuclideanRicciFlowWithOrbifold();
		//	break;
		//case OpenMesh::CETM:
		//	RunCETM();
		//	break;
		//case OpenMesh::HyperbolicCalabi:
		//	if (!orbifold) RunHyperbolicCalabiFlow();
		//	else RunHyperbolicCalabiFlowWithOrbifold();
		//	break;
		//case OpenMesh::HyperbolicRicci:
		//	if (!orbifold) RunHyperbolicRicciFlow();
		//	else RunHyperbolicRicciFlowWithOrbifold();
		//	break;
		//case OpenMesh::SphericalCalabi:
		//	if (!orbifold) RunSphericalCalabiFlow();
		//	else RunSphericalCalabiFlowWithOrbifold();
		//	break;
		//case OpenMesh::SphericalRicci:
		//	if (!orbifold) RunSphericalRicciFlow();
		//	else RunSphericalRicciFlowWithOrbifold();
		//	break;
		default:
			break;
		}

	}

	void ControlerFlow::RunEuclideanRicciFlow()
	{
		//EuclideanRicciFlow euc_ricci(m_mesh);
		//switch (boundary_mode) {
		//case 0:
		//	euc_ricci.boundary_mode = Fixed;
		//	break;
		//case 1:
		//	euc_ricci.boundary_mode = Circle;
		//	break;
		//case 2:
		//	euc_ricci.boundary_mode = Free;
		//}
		//euc_ricci.CalculateMetric();

		//m_mesh->RequestBoundary();
		//if (m_mesh->boundaries.size() == 0) {
		//	CMeshReconstructor Slice(m_mesh);
		//	MeshFlow* new_mesh = Slice.SliceToDisk();
		//	delete m_mesh;
		//	m_mesh = new_mesh;
		//}

		////embed
		//EuclideanEmbed euc_ebd(m_mesh);
		//euc_ebd.embed();
	}

	void ControlerFlow::RunEuclideanCalabiFlow()
	{
		EuclideanCalabiFlow euc_calabi(m_mesh);
		switch (boundary_mode) {
		case 0:
			euc_calabi.boundary_mode = Fixed;
			break;
		case 1:
			euc_calabi.boundary_mode = Circle;
			break;
		case 2:
			euc_calabi.boundary_mode = Free;
		}
		euc_calabi.CalculateMetric();

		m_mesh->RequestBoundary();
		//if (m_mesh->boundaries.size() == 0) {
		//	CMeshReconstructor Slice(m_mesh);
		//	MeshFlow* new_mesh = Slice.SliceToDisk();
		//	delete m_mesh;
		//	m_mesh = new_mesh;
		//}

		//embed
		EuclideanEmbed euc_ebd(m_mesh);
		euc_ebd.embed();
	}

	void ControlerFlow::RunHyperbolicRicciFlow()
	{
		//HyperbolicRicciFlow hyper_ricci(m_mesh);
		//hyper_ricci.CalculateMetric();

		//CMeshReconstructor Slice(m_mesh);
		//MeshFlow* new_mesh = Slice.SliceToDisk();
		//delete m_mesh;
		//m_mesh = new_mesh;

		//HyperbolicEmbed hyper_ebd(m_mesh);
		//hyper_ebd.embed();
	}

	void ControlerFlow::RunHyperbolicCalabiFlow()
	{
		//HyperbolicCalabiFlow hyper_ricci(m_mesh);
		//hyper_ricci.CalculateMetric();

		//CMeshReconstructor Slice(m_mesh);
		//MeshFlow* new_mesh = Slice.SliceToDisk();
		//delete m_mesh;
		//m_mesh = new_mesh;

		//HyperbolicEmbed hyper_ebd(m_mesh);
		//hyper_ebd.embed();
	}

	void ControlerFlow::RunSphericalRicciFlow()
	{
		//SphericalRicciFlow sph_ricci(m_mesh);
		//sph_ricci.CalculateMetric();
		//SphericalEmbed sph_ebd(m_mesh);
		//sph_ebd.embed();
	}

	void ControlerFlow::RunSphericalCalabiFlow()
	{
		//SphericalRicciFlow sph_calabi(m_mesh);
		//sph_calabi.CalculateMetric();
		//SphericalEmbed sph_ebd(m_mesh);
		//sph_ebd.embed();
	}


}// end namespace OpenMesh
#endif // !_FLOW_CONTROLER_H_

#pragma once
