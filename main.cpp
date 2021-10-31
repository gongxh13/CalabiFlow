#include "ControlerFlow.h"
#include <string.h>

int main(int argc, char* argv[])
{
	// Load a mesh in OFF format
	OpenMesh::ControlerFlow controler;
	controler.boundary_mode = OpenMesh::Fixed;
	controler.orbifold = false;
	controler.flow_type = OpenMesh::EuclideanCalabi;

	printf("===========Start======================\n");

	std::string fileName = "C:\\Users\\dell\\Desktop\\PolycubeSpaceCode\\PolycubeModels\\airplane1_input.obj";
	const char* p = fileName.c_str();
	bool result = controler.LoadMeshFromFile(p);

	controler.RunFlow();

	std::string fileNameSae = "C:\\Users\\dell\\Desktop\\gongxhfile\\project\\CalabiExtract\\result\\enbed\\airplane1_input.obj";
	const char* pSave = fileNameSae.c_str();
	result = controler.SaveMeshToFile(pSave);
	printf("===========End======================\n");
	// Plot the mesh

	return 0;
}
