/* 
* Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom
* the Software is furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
* OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _BVH_H_
#define _BVH_H_

#include <vector>
#include <algorithm>

#include "bunnykiller.h"
#include "Geometry/Aggregate/Aggregate.h"
#include "Geometry/Primitives/3D/Triangle.h"

template<class T, unsigned D>
class BVH: public Aggregate<T,D>
{
public:
	// Although Wald [2007] suggests that using just 8 bins is enough for
	// BVH construction, we opt for a more conservative approach, since we
	// do not aim for super-efficient construction.
	static int nb_bins_construction;
	// This is a rather standard number of primitives.
	static int min_nb_primitives_leaf;
	static int max_nb_primitives_leaf;
private:
	/** Class modelling the nodes in the hierarchy, including leaves */
	struct Node
	{
		std::vector<int>::iterator begin, end;
		AABB<D> bb;
		Node *child1, *child2;

		Node(): child1(0), child2(0){}
	};

	std::vector<Node> nodes;
	std::vector<int> indices;

	// Class Comp, used for comparison when partitioning primitives in 
	// the building of the tree.
	template<int axis>
	struct Comp
	{
		Comp(BVH<T,D> *bvh): b_bvh(bvh){}
		BVH<T,D> *b_bvh;
		bool operator ()(const int obj1, const int obj2)const
		{
			if( b_bvh->primitives[obj1].get_center()[axis] <= b_bvh->primitives[obj2].get_center()[axis])
				return false;

			return true;
		}
		
	};
	friend struct Comp<0>; friend struct Comp<1>; friend struct Comp<2>;
	Comp<0> comp0; Comp<1> comp1; Comp<2> comp2;

	void sort_objects(const std::vector<int>::iterator &begin, const std::vector<int>::iterator &end, const int axis);

	void build_bounding_box(const std::vector<int>::iterator &begin, const std::vector<int>::iterator &end);
	
	/** Divide the primitives contained in 'objects' into two
		groups, using a SAH. The binning technique used cames from:
		Wald, I. 2007. "On fast Construction of SAH-based Bounding Volume 
		Hierarchies". In Symposium on Interactive Ray Tracing 2007. 
		http://www.sci.utah.edu/~wald/Publications/2007/FastBuild/download/fastbuild.pdf */ 
	void split_objects_sah(const Node &node, Node **node1, Node **node2);
	void split_objects_median(const Node &node, Node **node1, Node **node2);
	void split_objects(const Node &node, Node **node1, Node **node2);

	void build_tree(Node &node);
	
	/** Recursive traversal of the tree to find an intersection */
	bool intersect(const Node &node, Ray<D>& r, Intersection<D> &it, float max_distance ) const;
	bool intersect(const Node &node, const Ray<D>& r, float max_distance ) const;
	
	/** Recursive traversal of tree to find all primitives intersecting a ray along its way
	 */
	bool intersect_all(const Node &node, const Ray<D>& r, float max_distance, std::vector<const T*>& intersected_primitives) const;
	bool intersect_all(const Node &node, const Ray<D>& r, float max_distance, std::vector<Intersection<D>*>& intersections) const;
	bool intersect_all(const Node &node, const VectorN<D>& p, std::vector<const T*>& intersected_primitives) const;
	
public:
	BVH();
	BVH(std::vector<T> &objects);
	~BVH();

	void freeze();
	void clear();

	bool intersect(Ray<D>& r, Intersection<D> &it, float max_distance ) const;
	bool intersect(const Ray<D>& r, float max_distance ) const;

	bool intersect_all(const Ray<D>& r, float max_distance, std::vector<const T*>& intersected_primitives) const;
	bool intersect_all(const Ray<D>& r, float max_distance, std::vector<Intersection<D>*>& intersections) const;
	bool intersect_all(const VectorN<D>& p, std::vector<const T*>& intersected_primitives) const;

	Real get_intersection_cost()const{return 0.;}
}; //BVH

//-----------------------------------------------------------------------
// Constructors & Desctructor
template<class T, unsigned D>
BVH<T,D>::BVH():Aggregate<T,D>(),comp0(this),comp1(this),comp2(this)
{

}

template<class T, unsigned D>
BVH<T,D>::BVH(std::vector<T> &objects):Aggregate<T,D>(objects),comp0(this),comp1(this),comp2(this)
{
	freeze();
}

template<class T, unsigned D>
BVH<T,D>::~BVH()
{
	clear();
}

//------------------------------------------------------------------------
// Functions to build the tree...
template<class T, unsigned D>
void BVH<T,D>::freeze()
{
	if( this->frozen )
		return;

	// Clear tree data to rebuild it
	indices.clear();
	nodes.clear();

	// Create primitives index
	indices.reserve(this->primitives.size());
	for(int i=0; i<this->primitives.size(); ++i)
		indices.push_back(i);

	// Since each primitive is referenced only once, the number of nodes
	// is bounded by 2N-1, with N the number of primitives.
	nodes.reserve(2*this->primitives.size()-1);

	// Create root node
	Node n;
	nodes.push_back(n);
	nodes.back().begin = indices.begin(); nodes.back().end = indices.end();
	
	nodes.back().bb = this->bb;
	
	// Recursively create BVH based on root node
	build_tree(nodes.back());
	
	// Pack together the primitives that belong to the same leaf
	std::vector<T> aux_primitives = this->primitives;
	for(int i=0; i<this->primitives.size(); ++i)
	{
		this->primitives[i] = aux_primitives[indices[i]];
		indices[i] = i;
	}

	this->frozen = true;
}

template<>
void BVH<Triangle, 3>::freeze()
{
	if (this->frozen)
		return;

	// Clear tree data to rebuild it
	indices.clear();
	nodes.clear();

	
	AABB<3> bb = AABB<3>();
	for (unsigned int i = 0; i < primitives.size(); ++i)
		bb.add(this->primitives[i].get_bounding_box());
	Real max_diagonal = bb.get_diagonal()/16;

	std::list<Triangle> triangles;
	for (unsigned int i = 0; i < this->primitives.size(); ++i)
		primitives[i].tessellate(triangles, max_diagonal);
	
	
	// Create primitives index
	primitives.clear(); 
	primitives.reserve(triangles.size());
	indices.reserve(triangles.size());
	
	int i = 0;
	while (!triangles.empty())
	{
		primitives.push_back(triangles.front()); 
		triangles.pop_front();
		indices.push_back(i);

		++i;
	}

	// Since each primitive is referenced only once, the number of nodes
	// is bounded by 2N-1, with N the number of primitives.
	nodes.reserve(2 * this->primitives.size() - 1);

	// Create root node
	Node n;
	nodes.push_back(n);
	nodes.back().begin = indices.begin(); nodes.back().end = indices.end();

	nodes.back().bb = this->bb;

	// Recursively create BVH based on root node
	build_tree(nodes.back());

	// Pack together the primitives that belong to the same leaf
	std::vector<Triangle> aux_primitives = this->primitives;
	for (unsigned int i = 0; i<this->primitives.size(); ++i) {
		this->primitives[i] = aux_primitives[indices[i]];
		indices[i] = i;
	}

	this->frozen = true;
}

template<class T, unsigned D>
void BVH<T,D>::build_tree( Node &node )
{
	Node *n1 = 0, *n2 = 0;
	split_objects(node, &n1, &n2);
	if( n1 && n2 ) //Note that it's actually not needed comparing both, but that way the code is more cute.
	{	
		node.child1 = n1;
		node.child2 = n2;

		build_tree(*n1);
		build_tree(*n2);
	}
}

template<class T, unsigned D>
void BVH<T,D>::build_bounding_box(const std::vector<int>::iterator &begin, const std::vector<int>::iterator &end)
{
	this->bb = AABB<D>();
	for (std::vector<int>::iterator i=begin; i!=end; ++i)
		this->bb.add(this->primitives[*i].get_bounding_box());
}

template<class T, unsigned D>
void BVH<T,D>::sort_objects(const std::vector<int>::iterator &begin, const std::vector<int>::iterator &end, const int axis)
{
	switch(axis)
	{
	case 0:
		std::stable_sort(begin, end, comp0);
		break;
	case 1:
		std::stable_sort(begin, end, comp1);
		break;
	case 2:
		std::stable_sort(begin, end, comp2);
		break;
	}
}

template<class T, unsigned D>
void BVH<T, D>::split_objects(const Node &node, Node **node1, Node **node2)
{
#ifdef _SAH_
	split_objects_sah(node, node1, node2);
#else
	split_objects_median(node, node1, node2);
#endif
}

template<class T, unsigned D>
void BVH<T, D>::split_objects_sah(const Node &node, Node **node1, Node **node2)
{
	if( (node.end-node.begin) <= min_nb_primitives_leaf )
	{	(*node1) = (*node2) = 0; return;	}

	int axis = node.bb.get_max_axis();
	int dominant_axis = 0;
	Real dim = node.bb.get_center()[axis];

	// We first compute the cost of not-traversing the node
	Real C_nt = 0;
	for( std::vector<int>::iterator i = node.begin; i != node.end; ++i )
		C_nt += this->primitives[*i].get_intersection_cost();
	
	// We compute the individual intersection cost as the average cost;
	Real C_i = 8.;//C_nt / static_cast<Real>(node.end-node.begin);
	
	// Sort the primitives
	sort_objects(node.begin, node.end, axis);

	// We now compute the SAH heuristic for the defined bins, so we keep the 
	// best possible partition.
	Real delta_position = node.bb.dimensions()[axis]/static_cast<Real>(nb_bins_construction+1);
	Real current_position = delta_position+node.bb._min[axis];
	
	// Precompute all fixed values to be used
	Real C_t = 1.;//11471.59804;//2*bb.get_intersection_cost();
	Real iSAV = 1./node.bb.area();
	
	// ... and use as best current case the end of traversal
	Real best_SA = (C_t > C_nt)?C_nt:std::numeric_limits<Real>::infinity();
	std::vector<int>::iterator best_switchpoint = node.end;
	AABB<D> best_bb[2];
	int N_c[2]; N_c[0]=N_c[1]=0;
	while(1) {
		int min_nb_bin=0;
		Real min_nb_position;
		bool to_switch_nb_bin = true;

		// ... then compute SAH of each bin
		for( int b=0; b<nb_bins_construction; ++b) {
			AABB<D> bb_c[2]; int c=0;
			N_c[0]=N_c[1]=0;
			bool toswitch = true; std::vector<int>::iterator switchpoint = node.end;
			
			// For each primitive in the node evaluate at wich side of the partition
			// it is (its centroid), and add it to each respective bounding box
			for( std::vector<int>::iterator i = node.begin; i != node.end; ++i )
			{
				//printf("%f vs %f: %s\n",primitives[*i].get_center()[axis],current_position,toswitch?"true":"false");
				if( toswitch && (this->primitives[*i].get_center()[axis] < current_position) )
				{	c=1; toswitch=false; switchpoint=i; }
				
				

				bb_c[c].add(this->primitives[*i].get_bounding_box());
				++N_c[c];
			}
			
			Real p1 = bb_c[0].area()*iSAV, p2 = bb_c[1].area()*iSAV;
			Real sum_p = p1+p2;

			// Then, compute this partition's cost...
			Real SA = C_t + ((N_c[0])?(bb_c[0].area()*N_c[0]):0. + (N_c[1])?(bb_c[1].area()*N_c[1]):0.) * iSAV*C_i;
			
			// Set if a change between the number of primitives on one side or in the other, to
			// keep track the minimum in case the sampling is too coarse to spot it.
			if(to_switch_nb_bin && (N_c[0] < N_c[1]))
			{	to_switch_nb_bin = false; min_nb_bin = b-1; 
				min_nb_position=current_position-delta_position;	}
			
			// ... and if it's lower than the best current partition, keep it!
			// Note that, for each leaf, we continue until the nb_primitives_leaf is reached
			if( SA < best_SA && !(N_c[0]<min_nb_primitives_leaf || N_c[1]<min_nb_primitives_leaf ) )
			{
				best_SA = SA;
				best_switchpoint = switchpoint;
				best_bb[0] = bb_c[0]; 
				best_bb[1] = bb_c[1];
			}

			// ... finally, update the partition position.
			current_position += delta_position;
		}
		
		// If found a good partition, end recursion
		if( best_switchpoint != node.end )
			break;
		
		// If no good partition reached, and there's a small enough number 
		// of primitives, just end and store the current node as a leaf
		if( N_c[0]<max_nb_primitives_leaf || N_c[1]<max_nb_primitives_leaf )
		{	(*node1) = (*node2) = 0; return;	}	

		// Else, iterate in the change interval, looking at it
		// with finer sampling
		current_position = min_nb_position;
		delta_position /= static_cast<Real>(nb_bins_construction+1);
		// ...or, if delta is too small, change the partition axis
		/*if( delta_position < 1e-5 )
		{
			axis = node.bb.get_n_dominant_axis(dominant_axis);
			sort_objects(node.begin, node.end, axis);

			// We now compute the SAH heuristic for the defined bins, so we keep the 
			// best possible partition.
			delta_position = node.bb.dimensions()[axis]/static_cast<Real>(nb_bins_construction+1);
			current_position = delta_position+node.bb._min[axis];
		}*/
//		{	(*node1) = (*node2) = 0; return;	}
	}

	// ... and return the partitioned nodes
	Node n_1;
	n_1.begin = node.begin; 
	n_1.end = best_switchpoint; 
	n_1.bb = best_bb[0];
  	nodes.push_back(n_1); *node1 = &(nodes.back());

	Node n_2;
	n_2.begin = best_switchpoint; 
	n_2.end = node.end; 
	n_2.bb = best_bb[1];
  	nodes.push_back(n_2); *node2 = &nodes.back();
}

template<class T, unsigned D>
void BVH<T, D>::split_objects_median(const Node &node, Node **node1, Node **node2)
{
	if ((node.end - node.begin) <= min_nb_primitives_leaf) {
		(*node1) = (*node2) = 0; return;
	}

	int axis = node.bb.get_max_axis();
//	int dominant_axis = 0;
//	Real dim = node.bb.get_center()[axis];

	// Sort the primitives
	sort_objects(node.begin, node.end, axis);

	// ... and return the partitioned nodes
	Node n_1;
	
	int div_point = (node.end - node.begin) / 2;
	n_1.begin = node.begin;
	n_1.end = node.begin + div_point;

	n_1.bb = AABB<D>();
	for (std::vector<int>::iterator i = n_1.begin; i != n_1.end; ++i)
		n_1.bb.add(this->primitives[*i].get_bounding_box());
	
	nodes.push_back(n_1); *node1 = &(nodes.back());

	Node n_2;
	n_2.begin = node.begin + div_point;
	n_2.end = node.end;
	
	n_2.bb = AABB<D>();
	for (std::vector<int>::iterator i = n_2.begin; i != n_2.end; ++i)
		n_2.bb.add(this->primitives[*i].get_bounding_box());
	nodes.push_back(n_2); *node2 = &nodes.back();
}

//------------------------------------------------------------------------
// Clearing...
template<class T, unsigned D>
void BVH<T,D>::clear()
{
	Aggregate<T,D>::clear();
	nodes.clear();
	indices.clear();
}

//------------------------------------------------------------------------
// Tree traversal...
// Single Intersection
template<class T, unsigned D>
bool BVH<T,D>::intersect(Ray<D>& r, Intersection<D> &it, float max_distance ) const
{
	if (!this->frozen) {
		fprintf(stderr,"Warning! Intersecting over a non-build BVH! Using Aggregate's method...\n");
		return Aggregate<T,D>::intersect(r,max_distance);
	}

	Real t;
	if( nodes.front().bb.intersect(r,t) )
		return (t<max_distance) && intersect(nodes.front(),r,it,max_distance);

	return false;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect(const Ray<D>& r, float max_distance ) const
{
	if (!this->frozen) {
		fprintf(stderr,"Warning! Intersecting over a non-build BVH! Using Aggregate's method...\n");
		return Aggregate<T,D>::intersect(r,max_distance);
	}

	Real t;
	if( nodes.front().bb.intersect(r,t) )
		return (t<max_distance) && (intersect(nodes.front(),r,max_distance));

	return false;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect(const Node &node, Ray<D>& r, Intersection<D> &it, float max_distance ) const
{
	// If an internal node of the tree...
	if( node.child1 && node.child2 )
	{
		// Test if intersection against children's BB, and check ranges!
		Real t1, t2;
		bool i1 = node.child1->bb.intersect(r, t1) && t1<max_distance;
		bool i2 = node.child2->bb.intersect(r, t2) && t2<max_distance;

		// If it not, just return no intersection...
		if( !(i1 || i2) )
			return false;
		
		t1 = (i1)?t1:std::numeric_limits<Real>::infinity();
		t2 = (i2)?t2:std::numeric_limits<Real>::infinity();

		// Then, test which is the closest children check that it is 
		// closest than the current intersection, and traverse down 
		// the tree. Afterwards, traverse the second child if needed.
		if( t1 < t2 )
		{
			if( t1 < r.get_parameter())
			{
				bool hit = intersect(*node.child1, r, it, max_distance);
				if( t2 < r.get_parameter() )
					hit = hit | intersect(*node.child2, r, it, max_distance);
				
				return hit;
			}
		}
		else
		{
			if( t2 < r.get_parameter())
			{
				bool hit = intersect(*node.child2, r, it, max_distance);
				if( t1 < r.get_parameter() )
					hit = hit | intersect(*node.child1, r, it, max_distance);

				return hit;
			}
		}
		return false;
	}
	
	//Real t; 
	//return bb.intersect(r,t);

	// ... if a leave.
	bool hit(false);
	// Test intersection on all primitives
	for( std::vector<int>::iterator i = node.begin; i != node.end; ++i )
		hit = this->primitives[*i].intersect(r,it,max_distance) | hit;

	return hit;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect(const Node &node, const Ray<D>& r, float max_distance ) const
{
	// If an internal node of the tree...
	if( node.child1 && node.child2 )
	{
		// Test if intersection against children's BB, and check ranges!
		Real t1, t2;
		bool i1 = node.child1->bb.intersect(r, t1) && t1<max_distance;
		bool i2 = node.child2->bb.intersect(r, t2) && t2<max_distance;

		// If it not, just return no intersection...
		if( !(i1 || i2) )
			return false;
		

		if( i1 )
			if( intersect(*node.child1, r, max_distance ) )
				return true;

		if( i2 )
			if( intersect(*node.child2, r, max_distance ) )
				return true;
		
		return false;
	}
	
	// ... if a leave.
	bool hit(false);
	// Test intersection on all primitives
	for( std::vector<int>::iterator i = node.begin; !hit && i != node.end; ++i )
		hit = this->primitives[*i].intersect(r, max_distance) | hit;

	return hit;
}

// Multiple intersections...
template<class T, unsigned D>
bool BVH<T,D>::intersect_all(const Ray<D>& r, float max_distance, std::vector<const T*>& intersected_primitives) const
{
	if( !this->frozen )
	{
		//THIS SHOULD BE CHANGED BY AN intersect_all IMPLEMENTATION IN Geometry/Aggregate/Aggregate.h 
		fprintf(stderr,"Warning! Intersecting over a non-build BVH! Returning error...\n");
		return false;
	}
	
	Real t;
	if( nodes.front().bb.intersect(r,t) )
	{
		return (t<max_distance) && intersect_all(nodes.front(),r,max_distance,intersected_primitives);
	}
	return false;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect_all(const Ray<D>& r, float max_distance, std::vector<Intersection<D>*>& intersections) const
{
	if( !this->frozen )
	{
		//THIS SHOULD BE CHANGED BY AN intersect_all IMPLEMENTATION IN Geometry/Aggregate/Aggregate.h 
		fprintf(stderr,"Warning! Intersecting over a non-build BVH! Returning error...\n");
		return false;
	}
	
	Real t;
	if( nodes.front().bb.intersect(r,t) )
	{
		return (t<max_distance) && intersect_all(nodes.front(),r,max_distance,intersections);
	}
	return false;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect_all(const VectorN<D>& p, std::vector<const T*>& intersected_primitives) const
{
	if( !this->frozen )
	{
		//THIS SHOULD BE CHANGED BY AN intersect_all IMPLEMENTATION IN Geometry/Aggregate/Aggregate.h 
		fprintf(stderr,"Warning! Intersecting over a non-build BVH! Returning error...\n");
		return false;
	}
	
	if( nodes.front().bb.intersect(p) )
	{
		return intersect_all(nodes.front(),p,intersected_primitives);
	}
	return false;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect_all(const Node &node, const Ray<D>& r, float max_distance, std::vector<const T*>& intersected_primitives) const
{
	// If an internal node of the tree...
	if( node.child1 && node.child2 )
	{
		// Test if intersection against children's BB, and check ranges!
		Real t1, t2;
		bool i1 = node.child1->bb.intersect(r, t1) && t1<max_distance;
		bool i2 = node.child2->bb.intersect(r, t2) && t2<max_distance;
		bool hit1(false), hit2(false);
		// If it not, just return no intersection...
		if( !(i1 || i2) )
			return false;
		
		if( i1 )
			hit1 = intersect_all(*node.child1, r, max_distance, intersected_primitives);
		
		if( i2 )
			hit2 = intersect_all(*node.child2, r, max_distance, intersected_primitives);
		
		return hit1 || hit2;
	}
	
	// ... if a leave.
	bool hit(false), new_hit(false);
	// Test intersection on all primitives
	for( std::vector<int>::iterator i = node.begin; i != node.end; ++i )
	{
		new_hit = this->primitives[*i].intersect(r, max_distance);
		hit = new_hit | hit;
		
		if (new_hit) 
			intersected_primitives.push_back(&(this->primitives[*i]));
	}
	
	return hit;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect_all(const Node &node, const Ray<D>& r, float max_distance, std::vector<Intersection<D>*>& intersections) const
{
	// If an internal node of the tree...
	if( node.child1 && node.child2 )
	{
		// Test if intersection against children's BB, and check ranges!
		Real t1, t2;
		bool i1 = node.child1->bb.intersect(r, t1) && t1<max_distance;
		bool i2 = node.child2->bb.intersect(r, t2) && t2<max_distance;
		bool hit1(false), hit2(false);
		// If it not, just return no intersection...
		if( !(i1 || i2) )
			return false;
		
		if( i1 )
			hit1 = intersect_all(*node.child1, r, max_distance, intersections);
		
		if( i2 )
			hit2 = intersect_all(*node.child2, r, max_distance, intersections);
		
		return hit1 || hit2;
	}
	
	// ... if a leave.
	bool hit(false), new_hit(false);
	// Test intersection on all primitives
	for( std::vector<int>::iterator i = node.begin; i != node.end; ++i )
	{
		Intersection<D> *it = new Intersection<D>(r, NULL,VectorN<D>(), Vector2());
		new_hit = this->primitives[*i].intersect(r, *it, max_distance);
		hit = new_hit | hit;
		if (new_hit) 
			intersections.push_back(it);
		else
			delete[] it;
	}
	
	return hit;
}

template<class T, unsigned D>
bool BVH<T,D>::intersect_all(const Node &node, const VectorN<D>& p, std::vector<const T*>& intersected_primitives) const
{
	// If an internal node of the tree...
	if ( node.child1 && node.child2 ) {
		// Test if intersection against children's BB, and check ranges!
		Real t1, t2;
		bool i1 = node.child1->bb.intersect(p);
		bool i2 = node.child2->bb.intersect(p);
		bool hit1(false), hit2(false);
		// If it not, just return no intersection...
		if( !(i1 || i2) )
			return false;
		
		if( i1 )
			hit1 = intersect_all(*node.child1, p, intersected_primitives);
		
		if( i2 )
			hit2 = intersect_all(*node.child2, p, intersected_primitives);
		
		return hit1 || hit2;
	}
	
	// ... if a leave.
	bool hit(false), new_hit(false);
	// Test intersection on all primitives
	for ( std::vector<int>::iterator i = node.begin; i != node.end; ++i ) {
		new_hit = this->primitives[*i].intersect(p);
		hit = new_hit | hit;
		
		if (new_hit) 
			intersected_primitives.push_back(&(this->primitives[*i]));
	}
	
	return hit;
}

//------------------------------------------------------------------------
// Static data initialization...
template<class T, unsigned D>
int BVH<T,D>::nb_bins_construction = 6;
template<class T, unsigned D>
int BVH<T,D>::min_nb_primitives_leaf = 4;
template<class T, unsigned D>
int BVH<T,D>::max_nb_primitives_leaf = 8;

#endif //_BVH_H_
