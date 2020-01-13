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

#ifndef _KDTREE_H_
#define _KDTREE_H_

// This code has been adapted from Adolfo Munoz's Mjolnir-RT, developed
//	at Universidad de Zaragoza (Spain).
// Particularly, the functions to find the N-closest nearest neighbors
//	have been added to the original code.

#include <vector>
#include <list>
#include <cmath>
#include <algorithm>

#include "bunnykiller.h"

template <class T, unsigned N>
class KDTree
{
public:
	class Node
	{
	friend class KDTree;
	private:
		T d;
		std::vector<Real> p;
		int axis;
	public:
		Node() : p(N), axis(-1)
		{}

		Node(const std::vector<Real>& _p, const T& _data) :
			d(_data), p(_p), axis(-1)
		{
			if (_p.size() != N) throw("Wrong dimension number in kd-tree");
		}

		const std::vector<Real>& point() const
		{
			return p;
		}

		const T& data() const
		{
			return d;
		}
	};
private:
	std::vector<Node> nodes;
	std::vector<Node> balanced;
public:
	KDTree() {}
	
	void clear() {
		nodes.clear();
		balanced.clear();
	}

	void store(const std::vector<Real>& point, const T& data)
	{
		nodes.push_back(Node(point,data));
	}
private:
	static void median_split(std::vector<Node>& p, size_t start, size_t end, size_t median,
	                         size_t axis);

	static void balance_segment(std::vector<Node>& pbal, std::vector<Node>& porg, size_t index,
	                           size_t start, size_t end, const std::vector<Real>& bbmin,
							   const std::vector<Real>& bbmax);
public:
	void balance();

	static Real distance(const std::vector<Real>& p1, const std::vector<Real>& p2)
	{
		Real sol2 = 0.0;
		for (size_t i = 0; i<N; i++)
			sol2 += (p1[i] - p2[i])*(p1[i] - p2[i]);
		return std::sqrt(sol2);
	}
private:
	size_t closest(const std::vector<Real>& p, size_t index, size_t best) const;

	void find(const std::vector<Real>& p, size_t index, Real radius,
	          std::vector<const Node*> &nodes) const;

	void find(const std::vector<Real>& p, size_t index, size_t nb_elements,
	          Real &dist_worst, std::vector<const Node*> &nodes,
			  std::vector<std::pair<unsigned, Real> > &dist) const;

	//I've removed static for compiling problems 
	//static class HeapComparison
	class HeapComparison
	{
	public:
		bool operator()(const std::pair<unsigned, Real> &val1,
		                const std::pair<unsigned, Real> &val2) const
		{
			return val1.second < val2.second;
		}
	};
	void update_heap_nodes(const Node &node, const Real distance, size_t nb_elements,
						   std::vector<const Node*> &nodes,
						   std::vector<std::pair<unsigned, Real> > &dist)const;
public:
	//========================================================================================================
	// Fixed Radius
	size_t find(const std::vector<Real>& p, Real radius) const
	{
		std::vector<const Node*> local_nodes;
		find(p, 1, radius, local_nodes);
		return local_nodes.size();
	}

	size_t find(const std::vector<Real>& p, Real radius, std::vector<const Node*> &nodes) const
	{
		find(p, 1, radius, nodes);
		return nodes.size();
	}
	//========================================================================================================
	// Nearest Neighbor search
	void find(const std::vector<Real>& p, size_t nb_elements, std::vector<const Node*> &nodes,
	          Real &max_distance) const
	{
		nodes.clear(); 
		max_distance = std::numeric_limits<Real>::infinity();
		
		if (balanced.empty())
			return;

		Real md = max_distance;

		nodes.reserve(nb_elements);
		std::vector<std::pair<unsigned, Real> > dist;
		dist.reserve(nb_elements);

		find(p, 1, nb_elements, md, nodes, dist);
		max_distance = md;
	}

	const Node& find(const std::vector<Real>& p) const
	{
		return balanced[closest(p, 1, 1)];
	}

	inline const Node& operator[](const size_t idx)const
	{
		#ifdef _SAFE_CHECK_		
		if ( idx > balanced.size()-1) throw("Out-of-range");
		#endif 
		
		return balanced[idx];
	}

	inline size_t nb_elements() const {
		return balanced.size();
	}

	inline bool is_empty() const
	{
		return balanced.empty();
	}
};
//--------------------------------------------------------------------------------------------------
//Private Find(radius)
template <class T, unsigned N >
void KDTree<T, N>::find(const std::vector<Real>& p, size_t index, Real radius,
                       std::vector<const Node*> &nodes) const
{
	// We check if our node enters
	if (distance(balanced[index].point(),p) < radius) {
		nodes.push_back(&balanced[index]);
	}
	// Now we check that this is not a leaf node
	if (index < ((balanced.size()-1)/2)) {
		Real distaxis=p[balanced[index].axis] - balanced[index].point()[balanced[index].axis];
		if (distaxis < 0.0) { // left node first
			find(p, 2*index, radius,nodes);
			if (radius > std::abs(distaxis)) // Maybe we can find more nodes on the other child
				find(p,2*index + 1,radius,nodes);
		} else { //right node first
			find(p, 2*index + 1, radius,nodes);
			if (radius > std::abs(distaxis)) // Maybe we can find more nodes on the other child
				find(p,2*index,radius,nodes);
		}
	}
}

//
//--------------------------------------------------------------------------------------------------
//Private Find(N-Nearest Neighbors)
template <class T, unsigned N>
void KDTree<T, N>::update_heap_nodes(const Node &node, const Real distance, size_t nb_elements,
		std::vector<const Node*> &nodes, std::vector<std::pair<unsigned, Real> > &dist)const
{
	// If there's still buffer for  more, don't bother with heaps...
	if (nodes.size() < nb_elements) {
		dist.push_back(std::pair<unsigned, Real>(nodes.size(), distance));
		nodes.push_back(&node);

		//...unless you've reach max size, then prepare the heap...
		if (nodes.size() == nb_elements)
			std::make_heap(dist.begin(), dist.end(), HeapComparison());
	} else {
		size_t idx = dist.front().first;
		nodes[idx] = &node;
		// Pop removed element
		pop_heap(dist.begin(), dist.end(), HeapComparison()); dist.pop_back();
		// Push new one
		dist.push_back(std::pair<unsigned, Real>(idx, distance));
		push_heap(dist.begin(), dist.end(), HeapComparison());
	}
}
template <class T, unsigned N>
void KDTree<T,N>::find(const std::vector<Real>& p, size_t index, size_t nb_elements, Real &dist_worst,
		std::vector<const Node*> &nodes, std::vector<std::pair<unsigned, Real> > &dist) const
{
	Real aux;
	//We check if our node is better
	if ((aux = distance(balanced[index].point(),p)) < dist_worst) {
		update_heap_nodes(balanced[index], aux, nb_elements, nodes, dist);
		dist_worst = (nodes.size() < nb_elements) ?
				std::numeric_limits<Real>::infinity() :
				dist.front().second;
	}

	//Now we check that this is not a leaf node
	if (index < ((balanced.size() - 1)/2)) {
		Real distaxis = p[balanced[index].axis] - balanced[index].point()[balanced[index].axis];
		if (distaxis < 0.0) { // left node first
			find(p,2*index,nb_elements,dist_worst,nodes, dist);
			if (dist_worst > std::abs(distaxis)) // Maybe we can find more nodes on the other child
				find(p, 2*index + 1, nb_elements, dist_worst, nodes, dist);
		} else { //right node first
			find(p,2*index + 1, nb_elements, dist_worst, nodes, dist);
			if (dist_worst > std::abs(distaxis)) // Maybe we can find more nodes on the other child
				find(p, 2*index, nb_elements, dist_worst, nodes, dist);
		}
	}
}

//--------------------------------------------------------------------------------------------------
// Closest
template <class T, unsigned N>
size_t KDTree<T,N>::closest(const std::vector<Real>& p, size_t index, size_t best) const
{
	size_t sol = best;
	Real distbest=distance(p,balanced[best].point());
	Real aux;
	//We check if our node is better
	if ((aux=distance(balanced[index].point(),p))<distbest) { sol=index; distbest=aux; }
	// Now we check that this is not a leaf node
	if (index<((balanced.size()-1)/2)) {
		Real distaxis=p[balanced[index].axis] - balanced[index].point()[balanced[index].axis];
		if (distaxis < 0.0) { // left node first
			size_t candidate = closest(p,2*index,sol);
			if ((aux = distance(balanced[candidate].point(),p))<distbest) { sol=candidate; distbest=aux; }
			if (distbest > std::abs(distaxis)) // Maybe the best solution is on the other side
			{
				candidate=closest(p,2*index + 1,sol);
				if ((aux = distance(balanced[candidate].point(),p))<distbest) {
					sol=candidate; distbest=aux;
				}
			}
		} else { // right node first
			int candidate = closest(p, 2*index+1, sol);
			if ((aux = distance(balanced[candidate].point(),p))<distbest) {
				sol=candidate;
				distbest=aux;
			}
			if (distbest > std::abs(distaxis)) { // Maybe the best solution is on the other side
				candidate = closest(p,2*index,sol);
				if ((aux = distance(balanced[candidate].point(),p))<distbest) {
					sol=candidate;
					distbest=aux;
				}
			}
		}
	}
	return sol;
}

#define myswap(array,a,b) { Node aux=(array)[(a)]; (array)[(a)]=(array)[(b)]; (array)[(b)]=aux;}

//--------------------------------------------------------------------------------------------------
//Balance Tree
template <class T, unsigned N >
void KDTree<T, N>::median_split(std::vector<Node>& p, size_t start, size_t end, size_t median, size_t axis)
{
	size_t left = start;
	size_t right = end;

	while (right > left) {
		Real v = p[right].point()[axis];
		size_t i=left-1;
		size_t j=right;
		for (;;) {
			while(v > p[++i].point()[axis]);
			while(v < p[--j].point()[axis] && j>left);
			if (i >= j)
				break;
			myswap(p, i, j);
		}

		myswap(p,i,right);
		if (i >= median)
			right = i-1;
		if (i <= median)
			left = i+1;
	}
}

template <class T, unsigned N >
void KDTree<T, N>::balance_segment(std::vector<Node>& pbal, std::vector<Node>& porg, size_t index,
		size_t start, size_t end, const std::vector<Real>& bbmin, const std::vector<Real>& bbmax)
{
	size_t median=1;
	while((4*median)<=(end-start+1))
		median += median;

	if ((3*median)<=(end-start+1)) {
		median += median;
		median += start-1;
	} else
		median = end-median+1;

	// elegimos el eje más apropiado...
	unsigned axis=0;
	for (unsigned i=1;i<N;i++)
		if ((bbmax[i] - bbmin[i]) > (bbmax[axis] - bbmin[axis]))
			axis = i;

	// partimos el bloque de fotones por la mediana
	median_split(porg, start, end, median, axis);

	pbal[index]=porg[median];
	pbal[index].axis=axis;

	// y por último balanceamos recursivamente los bloques izquierdo y derecho
	if(median>start) {
		// balancear el segmento izquierdo
		if(start<median-1) {
			std::vector<Real> newbbmax=bbmax;
			newbbmax[axis]=pbal[index].point()[axis];
			balance_segment(pbal, porg, 2*index, start, median-1,bbmin,newbbmax);
		} else {
			pbal[2*index]=porg[start];
		}
	}

	if(median<end) {
		// balancear el segmento derecho
		if(median+1<end) {
			std::vector<Real> newbbmin=bbmin;
			newbbmin[axis]=pbal[index].point()[axis];
			balance_segment(pbal, porg, 2*index+1, median+1, end,newbbmin,bbmax);
		} else {
			pbal[2*index+1]=porg[end];
		}
	}
}

template <class T, unsigned N >
void KDTree<T,N>::balance()
{
	if (nodes.empty())
		return;

	std::vector<Node> aux(nodes.size()+1);
	balanced.resize(nodes.size()+1);
	size_t i;
	std::vector<Real> bbmax = nodes.front().point();
	std::vector<Real> bbmin = nodes.front().point();
	//balanced[0] and aux[0] do not contain any useful information
	for (i = 1; !(nodes.empty()); i++, nodes.erase(nodes.begin())) {
		aux[i] = nodes.front();
		for (size_t j = 0; j<N; j++) {
			if (aux[i].point()[j] < bbmin[j])
				bbmin[j] = aux[i].point()[j];
			if (aux[i].point()[j] > bbmax[j])
				bbmax[j] = aux[i].point()[j];
		}
	}
	nodes.clear();

	KDTree<T, N>::balance_segment(balanced, aux, 1, 1, balanced.size()-1,bbmin,bbmax);
}

#endif
