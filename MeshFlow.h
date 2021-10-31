#ifndef _FLOW_MESH_H_
#define _FLOW_MESH_H_



//==============OPENMESH=============

#define _USE_MATH_DEFINES
#define NOMINMAX
#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS

//=============OPENMESH===============

#include <map>
#include <vector>
#include <queue>
#include "MeshLib\core\Parser\parser.h"
#include "MeshLib\core\Parser\strutil.h"
#include <queue>
#include <algorithm>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef MAX_LINE
#define MAX_LINE 4096
#endif // !

/*
Customize the mesh data structure
*/

namespace OpenMesh {
	struct MyTraits :public DefaultTraits
	{
		typedef Vec3d Point;
		typedef Vec3d Normal;
		typedef Vec2d TexCoord2D;
		typedef Vec3d TexCoord3D;
		typedef Vec3d Color;
		VertexAttributes(OpenMesh::Attributes::TexCoord2D | OpenMesh::Attributes::TexCoord3D | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
		FaceAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::Status);
		HalfedgeAttributes(OpenMesh::Attributes::Status);
		EdgeAttributes(OpenMesh::Attributes::Status);
	};

	typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> BaseMesh;


	class MeshFlow : public BaseMesh
	{
	public:
		/*constructor*/
		MeshFlow();

		/*create new face*/
		FaceHandle CreateFace(std::vector<VertexHandle>, int id);
		/*create new vertex*/
		VertexHandle CreateVertex(Point p, int id);
		/*map from id to vertex*/
		VHandle IdVertex(int id);
		/*map from id to face*/
		FHandle IdFace(int id);


		/*IO mesh file*/
		/*read from .m file*/
		void ReadMFile(const char* input);

		/*read from .obj file*/
		void ReadObjFile(const char* input);

		/*read .m file*/
		void WriteMFile(const char* input);

		/*read .obj file*/
		void WriteObjFile(const char* input);

		/*get boundaries*/
		void RequestBoundary();
		std::vector<std::vector<HalfedgeHandle>> boundaries;

		VertexHandle base_point;
		std::vector<std::vector<HalfedgeHandle>> homotopy_basis;

		// commen functions to get data
		template<class ItemHandle>
		char& Touched(ItemHandle item) { return GetProperty<char>("Touched", item); }
		template<class ItemHandle>
		std::string& String(ItemHandle item) { return GetProperty<std::string>("string", item); }
		template<class ItemHandle>
		int& Id(ItemHandle item) { return GetProperty<int>("id", item); }


		// some functions to get data of vertex;
		int& Parent(VertexHandle v) { return GetProperty<int>("parent", v); }
		//char& root(VertexHandle v) { return GetProperty<char>("root", v); }
		int& Dist(VertexHandle v) { return GetProperty<int>("dist", v); }
		int& Idb(VertexHandle v) { return GetProperty<int>("idb", v); }
		int& Idx(VertexHandle v) { return GetProperty<int>("idx", v); }
		double& U(VHandle v) { return GetProperty<double>("u", v); }
		double& Curvature(VHandle v) { return GetProperty<double>("k", v); }
		double& TargetCurvature(VHandle v) { return GetProperty<double>("target_k", v); }
		int& Order(VHandle v) { return GetProperty<int>("order", v); }
		char& IsBranch(VHandle v) { return GetProperty<char>("branch", v); }
		int& MaxWedge(VHandle v) { return GetProperty<int>("max_wedge", v); }

		// functions to get data of edge;
		char& IsOnSlice(EdgeHandle e) { return GetProperty<char>("on_Slice", e); }
		double& Length(EdgeHandle e) { return GetProperty<double>("Length", e); }
		double& OriginalLength(EdgeHandle e) { return GetProperty<double>("OriginalLength", e); }
		double& Weight(EHandle e) { return GetProperty<double>("weight", e); }
		double& MetricWeight(EHandle e) { return GetProperty<double>("metric_weight", e); }


		// function to get data of halfeges;
		double& Angle(HHandle h) { return GetProperty<double>("angle", h); }
		int& Wedge(HHandle h) { return GetProperty<int>("wedge", h); }
		double& Dtheta1Du1(HHandle h) { return GetProperty<double>("Dtheta1Du1", h); }
		double& Dtheta1Du2(HHandle h) { return GetProperty<double>("Dtheta1Du2", h); }

		// function to get data of face;



		// connectivity 
		/*return the edge between two vertices*/
		EHandle VertexEdge(VertexHandle v0, VertexHandle v1);

		/*return one of the vertices of one edge*/
		VertexHandle EdgeVertex1(EdgeHandle e);

		/*return another one of the vertices of one edge*/
		VertexHandle EdgeVertex2(EdgeHandle e);

		/*return halfedge of edge*/
		HHandle EdgeHalfedge(EdgeHandle e, unsigned i);




		/*covering space*/
		template<class P>
		void PushCovering(VertexHandle v, P p);
		std::map<VertexHandle, std::vector<Vec2d>> flat_covering;
		std::map<VertexHandle, std::vector<Vec3d>> sphere_covering;

	protected:

		/*handle string property*/
		void ToString(VertexHandle);
		void ToString(EdgeHandle);
		void ToString(HalfedgeHandle);
		void ToString(FaceHandle);
		void FromString(VertexHandle);
		void FromString(EdgeHandle);
		void FromString(HalfedgeHandle);
		void FromString(FaceHandle);

		/*initial vertex properties*/
		void InitializeMeshProperties();
		void InitializeVertexProperties();
		void InitializeEdgeProperties();
		void InitializeFaceProperties();
		void InitializeHalfedgeProperties();
		template<class PropHandle>
		void RegisteProperty(PropHandle handle, std::string name);

		/*get property*/
		template<class T>
		T& GetProperty(std::string name, VertexHandle item_handle);
		template<class T>
		T& GetProperty(std::string name, EdgeHandle item_handle);
		template<class T>
		T& GetProperty(std::string name, FaceHandle item_handle);
		template<class T>
		T& GetProperty(std::string name, HalfedgeHandle item_handle);

		/*! map between vetex and its id*/
		std::map<int, VertexHandle> m_map_IdVertex;
		/*! map between face and its id*/
		std::map<int, FaceHandle> m_map_IdFace;

	};

	inline MeshFlow::MeshFlow()
	{
		InitializeMeshProperties();
	}

	inline FaceHandle MeshFlow::CreateFace(std::vector<VertexHandle> verts, int id)
	{
		FaceHandle f = add_face(verts);
		m_map_IdFace[id] = f;
		this->Id(f) = id;
		return f;
	}

	inline VertexHandle MeshFlow::CreateVertex(Point p, int id)
	{
		assert(id != 0);
		VertexHandle v = this->add_vertex(p);
		m_map_IdVertex[id] = v;
		this->Id(v) = id;
		return v;
	}

	inline VertexHandle MeshFlow::IdVertex(int id)
	{
		return m_map_IdVertex[id];
	}

	inline FaceHandle MeshFlow::IdFace(int id)
	{
		return m_map_IdFace[id];
	}

	inline void MeshFlow::ReadMFile(const char* input) {
		std::fstream is(input, std::fstream::in);

		if (is.fail())
		{
			fprintf(stderr, "Error in opening file %s\n", input);
			return;
		}

		char buffer[MAX_LINE];
		int id;

		while (is.getline(buffer, MAX_LINE))
		{

			std::string line(buffer);
			line = strutil::trim(line);

			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token == "Vertex")
			{
				stokenizer.nextToken();
				token = stokenizer.getToken();
				id = strutil::parseString<int>(token);

				Point p;
				for (int i = 0; i < 3; i++)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					p[i] = strutil::parseString<float>(token);
				}

				VertexHandle v = this->CreateVertex(p, id);

				if (!stokenizer.nextToken("\t\r\n")) continue;
				token = stokenizer.getToken();

				int sp = (int)token.find("{");
				int ep = (int)token.find("}");

				if (sp >= 0 && ep >= 0)
				{
					String(v) = token.substr(sp + 1, ep - sp - 1);
					FromString(v);
				}
				continue;
			}


			if (token == "Face")
			{

				stokenizer.nextToken();
				token = stokenizer.getToken();
				id = strutil::parseString<int>(token);

				std::vector<VertexHandle> verts;
				while (stokenizer.nextToken())
				{
					token = stokenizer.getToken();
					if (strutil::startsWith(token, "{")) break;
					int vid = strutil::parseString<int>(token);
					verts.push_back(IdVertex(vid));
				}

				FaceHandle f = CreateFace(verts, id);

				if (strutil::startsWith(token, "{"))
				{
					String(f) = strutil::trim(token, "{}");
				}
				continue;
			}

			//read in edge attributes
			if (token == "Edge")
			{
				stokenizer.nextToken();
				token = stokenizer.getToken();
				int id0 = strutil::parseString<int>(token);

				stokenizer.nextToken();
				token = stokenizer.getToken();
				int id1 = strutil::parseString<int>(token);


				VertexHandle v0 = IdVertex(id0);
				VertexHandle v1 = IdVertex(id1);

				EdgeHandle edge = VertexEdge(v0, v1);

				if (!stokenizer.nextToken("\t\r\n")) continue;
				token = stokenizer.getToken();

				int sp = (int)token.find("{");
				int ep = (int)token.find("}");

				if (sp >= 0 && ep >= 0)
				{
					String(edge) = token.substr(sp + 1, ep - sp - 1);
				}
				continue;
			}
		}
	}

	inline void MeshFlow::ReadObjFile(const char* input)
	{
		std::fstream f(input, std::fstream::in);
		if (f.fail()) return;

		char cmd[1024];

		int  vid = 1;
		int  fid = 1;

		std::vector<Vec2d> uvs;
		std::vector<Vec3d> normals;


		while (f.getline(cmd, 1024))
		{
			std::string line(cmd);
			line = strutil::trim(line);

			strutil::Tokenizer stokenizer(line, " \t\r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token == "v")
			{
				Vec3d p;
				for (int i = 0; i < 3; i++)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					p[i] = strutil::parseString<float>(token);
				}

				VertexHandle v = CreateVertex(p, vid);
				vid++;
				continue;
			}


			if (token == "vt")
			{
				Vec2d uv;
				for (int i = 0; i < 2; i++)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					uv[i] = strutil::parseString<float>(token);
				}
				uvs.push_back(uv);
				continue;
			}


			if (token == "vn")
			{
				Vec3d n;
				for (int i = 0; i < 3; i++)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					n[i] = strutil::parseString<float>(token);
				}
				normals.push_back(n);
				continue;
			}




			if (token == "f")
			{
				std::vector<VertexHandle> v;
				for (int i = 0; i < 3; i++)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();


					strutil::Tokenizer tokenizer(token, " /\t\r\n");

					int ids[3];
					int k = 0;
					while (tokenizer.nextToken())
					{
						std::string token = tokenizer.getToken();
						ids[k] = strutil::parseString<int>(token);
						k++;
					}


					v.push_back(m_map_IdVertex[ids[0]]);
					if (uvs.size() > 0 && normals.size() > 0) {
						set_texcoord2D(v[i], uvs[ids[1] - 1]);
						set_normal(v[i], normals[ids[2] - 1]);
					}
					else if (normals.size() > 0)
						set_normal(v[i], normals[ids[1] - 1]);
					else if (uvs.size() > 0) {
						set_texcoord2D(v[i], uvs[ids[1] - 1]);
					}

				}
				CreateFace(v, fid++);
			}
		}

		f.close();
	}

	inline void MeshFlow::WriteMFile(const char* output)
	{
		//write traits to string
		for (VertexIter viter = vertices_begin(); viter != vertices_end(); viter++)
		{
			VertexHandle pV = *viter;
			String(pV) = "";
			ToString(pV);
		}

		for (EdgeIter eiter = edges_begin(); eiter != edges_end(); eiter++)
		{
			EdgeHandle pE = *eiter;
			String(pE) = "";
			ToString(pE);
		}

		for (FaceIter fiter = faces_begin(); fiter != faces_end(); fiter++)
		{
			FaceHandle pF = *fiter;
			String(pF) = "";
			ToString(pF);
		}

		for (HalfedgeIter hiter = halfedges_begin(); hiter != halfedges_end(); ++hiter)
		{
			HalfedgeHandle h = *hiter;
			String(h) = "";
			ToString(h);
		}


		std::fstream _os(output, std::fstream::out);
		if (_os.fail())
		{
			fprintf(stderr, "Error is opening file %s\n", output);
			return;
		}


		//remove vertices
		for (VertexIter viter = vertices_begin(); viter != vertices_end(); viter++)
		{
			VertexHandle v = *viter;

			_os << "Vertex " << Id(v);

			for (int i = 0; i < 3; i++)
			{
				_os << " " << point(v)[i];
			}
			if (String(v).size() > 0)
			{
				_os << " " << "{" << String(v) << "}";
			}
			_os << std::endl;
		}

		for (FaceIter fiter = faces_begin(); fiter != faces_end(); fiter++)
		{
			FaceHandle f = *fiter;

			_os << "Face " << Id(f);
			for (FaceVertexIter fviter = fv_iter(f); fviter.is_valid(); ++fviter) {
				_os << " " << Id(*fviter);
			}

			if (String(f).size() > 0)
			{
				_os << " " << "{" << String(f) << "}";
			}
			_os << std::endl;
		}

		_os.close();
	}

	inline void MeshFlow::WriteObjFile(const char* output)
	{
		std::fstream _os(output, std::fstream::out);
		if (_os.fail())
		{
			fprintf(stderr, "Error is opening file %s\n", output);
			return;
		}

		int vid = 1;
		for (VertexIter viter = vertices_begin(); viter != vertices_end(); viter++)
		{
			VertexHandle v = *viter;
			Idx(v) = vid++;
		}

		for (VertexIter viter = vertices_begin(); viter != vertices_end(); viter++)
		{
			VertexHandle v = *viter;

			_os << "v";

			for (int i = 0; i < 3; i++)
			{
				_os << " " << point(v)[i];
			}
			_os << std::endl;

			_os << "vt";

			for (int i = 0; i < 2; i++)
			{
				_os << " " << texcoord2D(v)[i];
			}
			_os << std::endl;

			_os << "vn";

			for (int i = 0; i < 3; i++)
			{
				_os << " " << normal(v)[i];
			}
			_os << std::endl;

		}

		for (FaceIter fiter = faces_begin(); fiter != faces_end(); fiter++)
		{
			FaceHandle f = *fiter;

			_os << "f";

			for (FaceVertexIter fviter = fv_iter(f); fviter.is_valid(); ++fviter) {
				_os << " " << Idx(*fviter) << "/" << Idx(*fviter) << '/' << Idx(*fviter);
			}

			_os << std::endl;
		}

		_os.close();
	}

	inline void MeshFlow::RequestBoundary()
	{
		boundaries.clear();
		bool has_boundary = false;
		for (HalfedgeIter hiter = halfedges_begin(); hiter != halfedges_end(); ++hiter) {
			HalfedgeHandle h = *hiter;
			Touched(h) = false;
			if (is_boundary(h)) {
				has_boundary = true;
			}
		}
		if (!has_boundary) return;
		char stop;
		do {
			stop = true;
			HalfedgeHandle hs, hc;
			for (HalfedgeIter hiter = halfedges_begin(); hiter != halfedges_end(); ++hiter) {
				HalfedgeHandle h = *hiter;
				char Touched = GetProperty<char>("Touched", h);
				if (is_boundary(h) && !Touched) {

					stop = false;
					hs = h;
					hc = next_halfedge_handle(h);
					break;
				}
			}
			if (hs.is_valid()) {
				std::vector<HalfedgeHandle> boundary;
				boundary.push_back(hs);
				Touched(hs) = true;
				while (hc != hs) {
					boundary.push_back(hc);
					Touched(hc) = true;
					hc = next_halfedge_handle(hc);
				}
				boundaries.push_back(boundary);
			}
		} while (!stop);
	}

	inline EdgeHandle MeshFlow::VertexEdge(VertexHandle v0, VertexHandle v1)
	{
		HalfedgeHandle h01 = find_halfedge(v0, v1);
		return edge_handle(h01);
	}

	inline VertexHandle MeshFlow::EdgeVertex1(EdgeHandle e)
	{
		return from_vertex_handle(halfedge_handle(e, 0));
	}

	inline VertexHandle MeshFlow::EdgeVertex2(EdgeHandle e)
	{
		return to_vertex_handle(halfedge_handle(e, 0));
	}

	inline HalfedgeHandle MeshFlow::EdgeHalfedge(EdgeHandle e, unsigned i)
	{
		return halfedge_handle(e, i);
	}

	inline void MeshFlow::InitializeMeshProperties()
	{
		InitializeVertexProperties();
		InitializeEdgeProperties();
		InitializeFaceProperties();
		InitializeHalfedgeProperties();
	}

	inline void MeshFlow::InitializeVertexProperties()
	{
		VPropHandleT<std::string> v_string;
		RegisteProperty(v_string, "string");
		VPropHandleT<int> id;
		RegisteProperty(id, "id");
		VPropHandleT<int> idx;
		RegisteProperty(idx, "idx");
		VPropHandleT<char> Touched;
		RegisteProperty(Touched, "Touched");
		VPropHandleT<double> u;
		RegisteProperty(u, "u");
		VPropHandleT<double> k;
		RegisteProperty(k, "k");
		VPropHandleT<double> target_k;
		RegisteProperty(target_k, "target_k");
		VPropHandleT<int> idb;
		RegisteProperty(idb, "idb");
		VPropHandleT<int> order;
		RegisteProperty(order, "order");
		VPropHandleT<char> branch;
		RegisteProperty(branch, "branch");
		VPropHandleT<int> parent;
		RegisteProperty(parent, "parent");
		VPropHandleT<int> dist;
		RegisteProperty(dist, "dist");
		VPropHandleT<int> max_wedge;
		RegisteProperty(max_wedge, "max_wedge");

	}

	inline void MeshFlow::InitializeEdgeProperties()
	{
		EPropHandleT<std::string>  string;
		RegisteProperty(string, "string");
		EPropHandleT<double> Length;
		RegisteProperty(Length, "Length");
		EPropHandleT<double> OriginalLength;
		RegisteProperty(OriginalLength, "OriginalLength");
		EPropHandleT<double> metric_weight;
		RegisteProperty(metric_weight, "metric_weight");
		EPropHandleT<double> weight;
		RegisteProperty(weight, "weight");
		EPropHandleT<char> Touched;
		RegisteProperty(Touched, "Touched");
		EPropHandleT<char> on_Slice;
		RegisteProperty(on_Slice, "on_Slice");

	}

	inline void MeshFlow::InitializeFaceProperties()
	{
		FPropHandleT<std::string> f_string;
		RegisteProperty(f_string, "string");
		FPropHandleT<int> id;
		RegisteProperty(id, "id");
		FPropHandleT<int> idx;
		RegisteProperty(idx, "idx");
		FPropHandleT<char> Touched;
		RegisteProperty(Touched, "Touched");
		FPropHandleT<Vec3d> ergb;
		RegisteProperty(Touched, "ergb");
	}

	inline void MeshFlow::InitializeHalfedgeProperties()
	{
		HPropHandleT<std::string> h_string;
		RegisteProperty(h_string, "string");
		HPropHandleT<char> Touched;
		RegisteProperty(Touched, "Touched");
		HPropHandleT<double> angle;
		RegisteProperty(angle, "angle");
		HPropHandleT<int> wedge;
		RegisteProperty(wedge, "wedge");
		HPropHandleT<double> Dtheta1Du1;
		RegisteProperty(Dtheta1Du1, "Dtheta1Du1");
		HPropHandleT<double> Dtheta1Du2;
		RegisteProperty(Dtheta1Du2, "Dtheta1Du2");
	}

	template<class T>
	inline T& MeshFlow::GetProperty(std::string name, VertexHandle item_handle)
	{
		VPropHandleT<T> prop;
		this->get_property_handle(prop, name);
		return property(prop, item_handle);
	}

	template<class T>
	inline T& MeshFlow::GetProperty(std::string name, EdgeHandle item_handle)
	{
		EPropHandleT<T> prop;
		this->get_property_handle(prop, name);
		return property(prop, item_handle);
	}

	template<class T>
	inline T& MeshFlow::GetProperty(std::string name, FaceHandle item_handle)
	{
		FPropHandleT<T> prop;
		this->get_property_handle(prop, name);
		return property(prop, item_handle);
	}

	template<class T>
	inline T& MeshFlow::GetProperty(std::string name, HalfedgeHandle item_handle)
	{
		HPropHandleT<T> prop;
		this->get_property_handle(prop, name);
		return property(prop, item_handle);
	}

	template<class P>
	inline void MeshFlow::PushCovering(VertexHandle v, P p)
	{
		if (typeid(P) == typeid(Vec2d)) {
			Vec2d vp;
			vp[0] = p[0];
			vp[1] = p[1];
			flat_covering[v].push_back(vp);
		}
		else if (typeid(P) == typeid(Vec3d)) {
			Vec3d vp;
			vp[0] = p[0];
			vp[1] = p[1];
			vp[2] = p[2];
			sphere_covering[v].push_back(vp);
		}
	}

	template<class PropHandle>
	inline void MeshFlow::RegisteProperty(PropHandle handle, std::string name)
	{
		add_property(handle, name);
		property(handle).set_persistent(true);
	}

	inline void OpenMesh::MeshFlow::ToString(VertexHandle v)
	{
		MeshLib::CParser parser(String(v));
		parser._removeToken("uv");
		parser._removeToken("rgb");
		parser._toString(String(v));
		std::stringstream iss;
		iss << "uv=(" << texcoord2D(v)[0] << " " << texcoord2D(v)[1] << ") ";
		iss << "uv3=(" << texcoord3D(v)[0] << " " << texcoord3D(v)[1] << " " << texcoord3D(v)[2] << ") ";
		iss << "normal=(" << normal(v)[0] << " " << normal(v)[1] << " " << normal(v)[2] << ") ";
		iss << "idb=(" << Idb(v) << ") ";
		iss << "branch=(" << int(IsBranch(v)) << ") ";
		iss << "order=(" << this->Order(v) << ") ";

		iss << "flat_covering=(" << texcoord2D(v)[0] << " " << texcoord2D(v)[1];
		for (int i = 0; i < flat_covering[v].size(); i++) {
			iss << ";" << flat_covering[v][i][0] << " " << flat_covering[v][i][1];
		}
		iss << ") ";

		iss << "sphere_covering=(" << texcoord3D(v)[0] << " " << texcoord3D(v)[1] << " " << texcoord3D(v)[2];
		for (int i = 0; i < sphere_covering[v].size(); i++) {
			iss << ";" << sphere_covering[v][i][0] << " " << sphere_covering[v][i][1] << " " << sphere_covering[v][i][2];
		}
		iss << ") ";

		if (String(v).length() > 0)
		{
			String(v) += " ";
		}
		String(v) += iss.str();
		//std::cout<< String(v)<<std::endl;
	}

	inline void OpenMesh::MeshFlow::ToString(EdgeHandle e)
	{
	}

	inline void OpenMesh::MeshFlow::ToString(HalfedgeHandle h)
	{
	}

	inline void OpenMesh::MeshFlow::ToString(FaceHandle f)
	{
	}

	inline void OpenMesh::MeshFlow::FromString(VertexHandle v)
	{
		std::string m_string = GetProperty<std::string>("string", v);
		MeshLib::CParser parser(m_string);

		for (std::list<MeshLib::CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			MeshLib::CToken* token = *iter;

			if (token->m_key == "rgb")
			{
				Color p;
				std::string line = strutil::trim(token->m_value, "()");
				line.erase(0, line.find_first_not_of("()"));
				line.erase(line.find_last_not_of("()") + 1);
				std::istringstream iss(line);
				iss >> p[0] >> p[1] >> p[2];
				set_color(v, p);
				continue;
			}


			if (token->m_key == "uv")
			{
				TexCoord2D uv;
				std::string line = strutil::trim(token->m_value, "()");
				line.erase(0, line.find_first_not_of("()"));
				line.erase(line.find_last_not_of("()") + 1);
				std::istringstream iss(line);
				iss >> uv[0] >> uv[1];
				set_texcoord2D(v, uv);
				//std::cout << uv[0] << uv[1] << std::endl;
				continue;
			}

			if (token->m_key == "order")
			{
				std::string line = strutil::trim(token->m_value, "()");
				std::istringstream iss(line);
				iss >> Order(v);
				continue;
			}

			if (token->m_key == "idb")
			{
				std::string line = strutil::trim(token->m_value, "()");
				std::istringstream iss(line);
				iss >> Idb(v);
				continue;
			}

			if (token->m_key == "branch")
			{
				std::string line = strutil::trim(token->m_value, "()");
				std::istringstream iss(line);
				int _branch;
				iss >> _branch;
				if (_branch)
					IsBranch(v) = true;
				else
					IsBranch(v) = false;
				continue;
			}
		}
	}

	inline void OpenMesh::MeshFlow::FromString(EdgeHandle e)
	{
	}

	inline void OpenMesh::MeshFlow::FromString(HalfedgeHandle h)
	{
	}

	inline void OpenMesh::MeshFlow::FromString(FaceHandle f)
	{
	}


} //end namespace OpenMesh
#endif#pragma once
