/*!
*      \file DynamicMesh.h
*      \brief Dynamic Mesh Class for edge swap, face split, edge split
*
*      Dynamic Mesh Class for edge swap, face split, edge split
*	   \author David Gu
*      \date 10/07/2010
*
*/

#ifndef  _DYNAMIC_MESH_H_
#define  _DYNAMIC_MESH_H_

#include <map>
#include <vector>
#include <queue>

#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"


namespace MeshLib
{
  
/*-------------------------------------------------------------------------------------------------------------------------------------

	Dynamic Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CDynamicMesh class : Dynamic mesh
* 
*  Mesh supports FaceSlit, EdgeSlit, EdgeSwap operations
*/
template<typename V, typename E, typename F, typename H>
class CDynamicMesh : public CBaseMesh<V,E,F,H>
{
public:
	/*! CDynamicMesh constructor */
	CDynamicMesh(){ m_vertex_id = 0; m_face_id = 0; };
	/*! CDynamicMesh destructor */
	~CDynamicMesh();

	//dynamic mesh editing
	/*! Split a Face to three small faces 
	 * \param pFace the face to be split
	 */
	V *  splitFace( F * pFace );
	/*! Split one edge to two edges 
	 * \param pEdge the edge to be split
	 */
	V *  splitEdge( E * pEdge );
	/*! Swap an edge
	* \param edge the edge to be swapped
	*/
	void swapEdge( E * edge );
	
protected:
	/*! attach halfeges to an edge
	* \param he0, he1 the halfedges
	* \param e edge
	*/
	void __attach_halfedge_to_edge( H * he0, H * he1, E * e );
	/*! next vertex id */
	int  m_vertex_id;
	/*! next face id */
	int  m_face_id;
};

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
CDynamicMesh<V,E,F,H>::~CDynamicMesh()
{
}


/*---------------------------------------------------------------------------*/

//insert a vertex in the center of a face, split the face to 3 faces

template<typename CVertex, typename CEdge, typename CFace, typename CHalfEdge>
CVertex * CDynamicMesh<CVertex,CEdge,CFace,CHalfEdge>::splitFace( CFace * pFace )
{

	CVertex * pV = createVertex( ++m_vertex_id );
	
	CVertex   *  v[3];
    CHalfEdge *  h[3];
	CHalfEdge *  hs[3];
	
	CEdge     *  eg[3];

	h[0] =  faceHalfedge( pFace );
	h[1] =  faceNextCcwHalfEdge( h[0] );
	h[2] =  faceNextCcwHalfEdge( h[1] );

	for( int i = 0; i < 3; i ++ )
	{
		v[i] = halfedgeTarget( h[i] );
		eg[i] = halfedgeEdge( h[i] );
		hs[i] = halfedgeSym( h[i] );
	}
	

	CFace * f = new CFace();
	assert( f != NULL );
	f->id() = ++m_face_id;
	m_faces.push_back( f );

	//create halfedges
	tHalfEdge hes[3];
	for(int i = 0; i < 3; i ++ )
	{
		hes[i] = new CHalfEdge;
		assert( hes[i] );
	}

	//linking to each other
	for(int i = 0; i < 3; i ++ )
	{
		hes[i]->he_next() = hes[(i+1)%3];
		hes[i]->he_prev() = hes[(i+2)%3];
	}

	//linking to face
	for(int i = 0; i < 3; i ++ )
	{
		hes[i]->face()   = f;
		f->halfedge()    = hes[i];
	}


	f = new CFace();
	assert( f != NULL );
	f->id() = ++m_face_id;
	m_faces.push_back( f );

	//create halfedges
	tHalfEdge hes2[3];

	for(int i = 0; i < 3; i ++ )
	{
		hes2[i] = new CHalfEdge;
		assert( hes2[i] );
	}

	//linking to each other
	for(int i = 0; i < 3; i ++ )
	{
		hes2[i]->he_next() = hes2[(i+1)%3];
		hes2[i]->he_prev() = hes2[(i+2)%3];
	}

	//linking to face
	for(int i = 0; i < 3; i ++ )
	{
		hes2[i]->face()   = f;
		f->halfedge()    = hes2[i];
	}

	CEdge * e[3];
	for( int i = 0; i < 3; i ++ )
	{
		e[i] = new CEdge();
		assert( e[i] );
		m_edges.push_back( e[i] );
	}

	__attach_halfedge_to_edge( h[1], hes[0], e[0] );
	__attach_halfedge_to_edge( hes[2], hes2[1], e[1] );
	__attach_halfedge_to_edge( h[2], hes2[0], e[2] );
	__attach_halfedge_to_edge( h[0], hs[0], eg[0] );
	__attach_halfedge_to_edge( hes[1], hs[1], eg[1] );
	__attach_halfedge_to_edge( hes2[2], hs[2], eg[2] );
	


	pV->halfedge() = h[1];


	h[1]->vertex() = pV;
	h[2]->vertex() = v[2];

	hes[0]->vertex() = v[0];
	hes[1]->vertex() = v[1];
	hes[2]->vertex() = pV;

	hes2[0]->vertex() = pV;
	hes2[1]->vertex() = v[1];
	hes2[2]->vertex() = v[2];

	for( int i = 0; i < 3; i ++ )
	{
		v[i]->halfedge() = hs[i]->he_sym();
	}
	return pV;

};
/*---------------------------------------------------------------------------*/

template<typename CVertex, typename CEdge, typename CFace, typename CHalfEdge>
void CDynamicMesh<CVertex,CEdge,CFace,CHalfEdge>::swapEdge( CEdge * edge )
{

  CHalfEdge * he_left   = edgeHalfedge( edge, 0 );
  CHalfEdge * he_right  = edgeHalfedge( edge, 1 );

  if( he_right == NULL )  return;
  
  CHalfEdge * ph[6];

  ph[0] = he_left;
  ph[1] = faceNextCcwHalfEdge( ph[0] );
  ph[2] = faceNextCcwHalfEdge( ph[1] );


  ph[3] = he_right;
  ph[4] = faceNextCcwHalfEdge( ph[3] );
  ph[5] = faceNextCcwHalfEdge( ph[4] );

  CVertex * pv[4];

  pv[0] = halfedgeTarget( ph[0]);
  pv[1] = halfedgeTarget( ph[1]);
  pv[2] = halfedgeTarget( ph[2]);
  pv[3] = halfedgeTarget( ph[4]);

  bool boundary[4];

  for( int i = 0; i < 4; i ++ )
  {
    boundary[i] = pv[0]->boundary();
  }

  int pi[6];
  CEdge * pe[6];
  for( int i = 0; i < 6; i ++ )
  {
    CHalfEdge * he = ph[i];
    CEdge        *   e= halfedgeEdge( he );
    pe[i] = e;

    if( e->halfedge(0) == he )
    {
      pi[i] = 0;
    }
    else
    {
      assert( e->halfedge(1) == he );
      pi[i] = 1;
    }
  }

  //relink the vertices

  ph[0]->target() = pv[1];
  if( !boundary[1] )
  {
    pv[1]->halfedge() = ph[0];
  }
  ph[1]->target() = pv[2];
  if( !boundary[2] )
  {
    pv[2]->halfedge() = ph[1];
  }
  
  ph[2]->target() = pv[3];
  if( !boundary[3] )
  {
    pv[3]->halfedge() = ph[2];
  }
  
  ph[3]->target() = pv[3];
  if( !boundary[3] )
  {
    pv[3]->halfedge() = ph[3];
  }
  
  ph[4]->target() = pv[0];
  if( !boundary[0] )
  {
    pv[0]->halfedge() = ph[4];
  }
  
  ph[5]->target() = pv[1];
  if( !boundary[1] )
  {
    pv[1]->halfedge() = ph[5];
  }
  
  //relink the edge-halfedge pointers

  ph[1]->edge() = pe[2];
  pe[2]->halfedge( pi[2] ) = ph[1];

  ph[2]->edge() = pe[4];
  pe[4]->halfedge( pi[4] ) = ph[2];

  ph[4]->edge() = pe[5];
  pe[5]->halfedge( pi[5] ) = ph[4];

  ph[5]->edge() = pe[1];
  pe[1]->halfedge( pi[1] ) = ph[5];
  

  for(int i = 0; i < 6; i ++ )
  {
      CHalfEdge * he = ph[i];
      CHalfEdge * sh = halfedgeSym( he );
      assert( he->target() == sh->he_prev()->target() );
      assert( he->he_prev()->target() == sh->target() );
  }

};

/*---------------------------------------------------------------------------*/


//insert a vertex in the center of an edge, split each neighboring face into 2 faces

template<typename CVertex, typename CEdge, typename CFace, typename CHalfEdge>
CVertex * CDynamicMesh<CVertex,CEdge,CFace,CHalfEdge>::splitEdge( CEdge * pEdge )
{

	CVertex * pV = createVertex( ++ m_vertex_id );


	CHalfEdge *  h[12];
	CHalfEdge *  s[6];
	CVertex   *  v[6];
	CEdge     * eg[6];
	
	h[0] = edgeHalfedge( pEdge, 0 );
	h[1] = faceNextCcwHalfEdge( h[0] );
	h[2] = faceNextCcwHalfEdge( h[1] );

	h[3] = edgeHalfedge( pEdge, 1 );
	h[4] = faceNextCcwHalfEdge( h[3] );
	h[5] = faceNextCcwHalfEdge( h[4] );


	CFace * f[4];
	f[0] = halfedgeFace( h[0] );
	f[1] = halfedgeFace( h[3] );

	for( int i = 0; i < 6; i ++ )
	{
		eg[i] = halfedgeEdge( h[i] );
		v[i] = halfedgeVertex(h[i]);
		s[i] = halfedgeSym( h[i] );
	}

	f[2] = new CFace();
	assert( f[2] != NULL );
	f[2]->id() = ++ m_face_id;
	m_faces.push_back( f[2] );

	//create halfedges
	for(int i = 6; i < 9; i ++ )
	{
		h[i] = new CHalfEdge;
		assert( h[i] );
	}

	//linking to each other
	for(int i = 0; i < 3; i ++ )
	{
		h[i+6]->he_next() = h[6+(i+1)%3];
		h[i+6]->he_prev() = h[6+(i+2)%3];
	}

	//linking to face
	for(int i = 6; i < 9; i ++ )
	{
		h[i]->face()   = f[2];
		f[2]->halfedge()= h[i];
	}


	f[3] = new CFace();
	assert( f[3] != NULL );
	f[3]->id() = ++m_face_id;
	m_faces.push_back( f[3] );

	//create halfedges
	for(int i = 9; i < 12; i ++ )
	{
		h[i] = new CHalfEdge;
		assert( h[i] );
	}

	//linking to each other
	for(int i = 0; i < 3; i ++ )
	{
		h[i+9]->he_next() = h[9+(i+1)%3];
		h[i+9]->he_prev() = h[9+(i+2)%3];
	}

	//linking to face
	for(int i = 9; i < 12; i ++ )
	{
		h[i]->face()   = f[3];
		f[3]->halfedge()    = h[i];
	}

	CEdge * e[3];

	for( int i = 0; i < 3; i ++ )
	{
		e[i] = new CEdge();
		m_edges.push_back( e[i] );
		assert( e[i] );
	}
	
	__attach_halfedge_to_edge(h[2], h[6] , e[0]);
	__attach_halfedge_to_edge(h[8], h[9] , e[1]);
	__attach_halfedge_to_edge(h[4], h[11], e[2]);

	__attach_halfedge_to_edge(h[0], h[3], eg[0]);
	__attach_halfedge_to_edge(h[1], s[1], eg[1]);
	__attach_halfedge_to_edge(h[5], s[5], eg[5]);

	__attach_halfedge_to_edge(h[7],  s[2], eg[2]);
	__attach_halfedge_to_edge(h[10], s[4], eg[4]);


	h[0]->vertex() = v[0];
	h[1]->vertex() = v[1];
	h[2]->vertex() = pV;
	h[3]->vertex() = pV;
	h[4]->vertex() = v[4];
	h[5]->vertex() = v[5];
	h[6]->vertex() = v[1];
	h[7]->vertex() = v[2];
	h[8]->vertex() = pV;
	h[9]->vertex() = v[2];
	h[10]->vertex()= v[4];
	h[11]->vertex()= pV;

	

	v[0]->halfedge() = h[0];
	v[1]->halfedge() = h[1];
	v[2]->halfedge() = h[7];
	v[4]->halfedge() = h[4];
	pV->halfedge()   = h[3];

	for( int k = 0; k < 4; k ++ )
	{
		CHalfEdge * pH = faceHalfedge( f[k] );
		for( int i = 0; i < 3; i ++ )
		{
			assert( pH->vertex() == pH->he_sym()->he_prev()->vertex() );
			pH = faceNextCcwHalfEdge( pH );
		}
	}
	return pV;	

};

/*---------------------------------------------------------------------------*/

template<typename CVertex, typename CEdge, typename CFace, typename CHalfEdge>
void CDynamicMesh<CVertex,CEdge,CFace,CHalfEdge>::__attach_halfedge_to_edge( CHalfEdge * he0, CHalfEdge * he1, CEdge * e )
{
	e->halfedge(0 ) = he0;
	e->halfedge(1 ) = he1;

	he0->edge() = e;
	he1->edge() = e;

};



}
#endif  _DYNAMIC_MESH_H_