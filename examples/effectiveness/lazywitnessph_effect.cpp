/*
 
 Snap-Ripser: C++ code for computing of $\epsilon$-net induced Lazy witness persistence barcodes
  
 Copyright (c) National University of Singapore
 
 Author: Naheed Anjum Arafat on 6/10/19.
 
 Acknowledgement: Adapted from Ripser c++ (copyright: Ulrich Bauer, License: MIT)

*/

//#define ASSEMBLE_REDUCTION_MATRIX
//#define USE_COEFFICIENTS

//#define INDICATE_PROGRESS
#define PRINT_PERSISTENCE_PAIRS

#define USE_GOOGLE_HASHMAP
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include "stdafx.h"
#include <cstdlib>
#include <set>
#include <typeinfo>
#include <random>
#include <map>
#include <stack>
#include <tuple>
#include <functional>
#include <ctime>

#ifdef USE_GOOGLE_HASHMAP
#include <sparsehash/sparse_hash_map>
template <class Key, class T> class hash_map : public google::sparse_hash_map<Key, T> {
public:
    inline void reserve(size_t hint) { this->resize(hint); }
};
#else
template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};
#endif

typedef float value_t;
typedef int64_t index_t;
typedef int16_t coefficient_t;

std::vector<std::string> splitpath(
                                   const std::string& str
                                   , const std::set<char> delimiters)
{
    std::vector<std::string> result;
    
    char const* pch = str.c_str();
    char const* start = pch;
    for(; *pch; ++pch)
    {
        if (delimiters.find(*pch) != delimiters.end())
        {
            if (start != pch)
            {
                std::string str(start, pch);
                result.push_back(str);
            }
            else
            {
                result.push_back("");
            }
            start = pch + 1;
        }
    }
    result.push_back(start);
    
    return result;
}

class binomial_coeff_table {
    std::vector<std::vector<index_t>> B;
    
public:
    binomial_coeff_table(index_t n, index_t k) : B(n + 1) {
        for (index_t i = 0; i <= n; i++) {
            B[i].resize(k + 1);
            for (index_t j = 0; j <= std::min(i, k); j++)
                if (j == 0 || j == i)
                    B[i][j] = 1;
                else
                    B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
        }
    }
    
    index_t operator()(index_t n, index_t k) const {
        assert(n < B.size() && k < B[n].size());
        return B[n][k];
    }
};

bool is_prime(const coefficient_t n) {
    if (!(n & 1) || n < 2) return n == 2;
    for (coefficient_t p = 3, q = n / p, r = n % p; p <= q; p += 2, q = n / p, r = n % p)
        if (!r) return false;
    return true;
}

coefficient_t normalize(const coefficient_t n, const coefficient_t modulus) {
    return n > modulus/2 ? n - modulus : n;
}

std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m) {
    std::vector<coefficient_t> inverse(m);
    inverse[1] = 1;
    // m = a * (m / a) + m % a
    // Multipying with inverse(a) * inverse(m % a):
    // 0 = inverse(m % a) * (m / a)  + inverse(a)  (mod m)
    for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
    return inverse;
}

#ifdef USE_COEFFICIENTS
struct __attribute__((packed)) entry_t {
    index_t index : 8 * (sizeof(index_t) - sizeof(coefficient_t));
    coefficient_t coefficient;
    entry_t(index_t _index, coefficient_t _coefficient)
    : index(_index), coefficient(_coefficient) {}
    entry_t(index_t _index) : index(_index), coefficient(1) {}
    entry_t() : index(0), coefficient(1) {}
};

static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

entry_t make_entry(index_t _index, coefficient_t _coefficient) {
    return entry_t(_index, _coefficient);
}
index_t get_index(const entry_t& e) { return e.index; }
index_t get_coefficient(const entry_t& e) { return e.coefficient; }
void set_coefficient(entry_t& e, const coefficient_t c) { e.coefficient = c; }

bool operator==(const entry_t& e1, const entry_t& e2) {
    return get_index(e1) == get_index(e2) && get_coefficient(e1) == get_coefficient(e2);
}

std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
    stream << get_index(e) << ":" << get_coefficient(e);
    return stream;
}

#else

typedef index_t entry_t;
const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }
void set_coefficient(entry_t& e, const coefficient_t c) {}

#endif

const entry_t& get_entry(const entry_t& e) { return e; }

typedef std::pair<value_t, index_t> diameter_index_t;
value_t get_diameter(const diameter_index_t& i) { return i.first; }
index_t get_index(const diameter_index_t& i) { return i.second; }

typedef std::pair<index_t, value_t> index_diameter_t;
index_t get_index(const index_diameter_t& i) { return i.first; }
value_t get_diameter(const index_diameter_t& i) { return i.second; }

class diameter_entry_t : public std::pair<value_t, entry_t> {
public:
    diameter_entry_t() {}
    diameter_entry_t(const entry_t& e) : std::pair<value_t, entry_t>(0, e) {}
    diameter_entry_t(value_t _diameter, index_t _index, coefficient_t _coefficient)
    : std::pair<value_t, entry_t>(_diameter, make_entry(_index, _coefficient)) {}
    diameter_entry_t(const diameter_index_t& _diameter_index, coefficient_t _coefficient)
    : std::pair<value_t, entry_t>(get_diameter(_diameter_index),
                                  make_entry(get_index(_diameter_index), _coefficient)) {}
    diameter_entry_t(const diameter_index_t& _diameter_index)
    : diameter_entry_t(_diameter_index, 1) {}
};

const entry_t& get_entry(const diameter_entry_t& p) { return p.second; }
entry_t& get_entry(diameter_entry_t& p) { return p.second; }
const index_t get_index(const diameter_entry_t& p) { return get_index(get_entry(p)); }
const coefficient_t get_coefficient(const diameter_entry_t& p) {
    return get_coefficient(get_entry(p));
}
const value_t& get_diameter(const diameter_entry_t& p) { return p.first; }
void set_coefficient(diameter_entry_t& p, const coefficient_t c) {
    set_coefficient(get_entry(p), c);
}

template <typename Entry> struct greater_diameter_or_smaller_index {
    bool operator()(const Entry& a, const Entry& b) {
        return (get_diameter(a) > get_diameter(b)) ||
        ((get_diameter(a) == get_diameter(b)) && (get_index(a) < get_index(b)));
    }
};

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> class compressed_distance_matrix {
public:
    std::vector<value_t> distances;
    std::vector<value_t*> rows;
    
    void init_rows();
    
    // Called when passed either upper/lower triangular matrix.
    compressed_distance_matrix(std::vector<value_t>&& _distances)
    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
        assert(distances.size() == size() * (size() - 1) / 2);
        init_rows();
    }
    
    // Gets called when passed square dist matrix
    template <typename DistanceMatrix>
    compressed_distance_matrix(const DistanceMatrix& mat)
    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
        init_rows();
        
        for (index_t i = 1; i < size(); ++i)
            for (index_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
    }
    
    value_t operator()(const index_t i, const index_t j) const;
    
    size_t size() const { return rows.size(); }
    
    
    void print() const ;
};

class sparse_distance_matrix {
public:
    std::vector<std::vector<index_diameter_t>> neighbors;
    
    index_t num_edges;
    
    sparse_distance_matrix(std::vector<std::vector<index_diameter_t>>&& _neighbors,
                           index_t _num_edges)
    : neighbors(std::move(_neighbors)), num_edges(_num_edges) {}
    
    template <typename DistanceMatrix>
    sparse_distance_matrix(const DistanceMatrix& mat, value_t threshold)
    : neighbors(mat.size()), num_edges(0) {
        
        for (index_t i = 0; i < size(); ++i)
            for (index_t j = 0; j < size(); ++j)
                if (i != j && mat(i, j) <= threshold) {
                    ++num_edges;
                    neighbors[i].push_back(std::make_pair(j, mat(i, j)));
                }
    }
    
    size_t size() const { return neighbors.size(); }
    void print() const {
        for (int i = 0; i< size(); i++){
            //            neighbors[i].sort();
            for (int j = 0; j< neighbors[i].size(); j++){
                std::cout << neighbors[i][j].second << " ";
            }
            std::cout << std::endl;
        }
    }
};


template <> void compressed_distance_matrix<LOWER_TRIANGULAR>::init_rows() {
    //    std::cout << " init lower tri \n";
    value_t* pointer = &distances[0];
    //    std::cout << *pointer <<" ";
    for (index_t i = 1; i < size(); ++i) {
        rows[i] = pointer;
        //        std::cout << *pointer <<" ";
        pointer += i;
    }
    
    //    // printing
    //    std::cout <<"\n * * * * * * * * * \n";
    //    for (index_t i = 1; i< size(); ++i){
    //        for(index_t j = 0; j< i; ++j)
    //            if (i==j) {std::cout << 0 <<" "; continue;}
    //            else{
    //                std::cout << *(rows[i]+j)<<" ";
    //            }
    //    std::cout <<"\n";
    //    }
    //    std::cout <<" ** * ** * ** * *\n";
}
// Called from compressed_distance_matrix(std::vector<value_t>&& _distances) constructor when template typename is UPPER_TRIANGULAR
template <> void compressed_distance_matrix<UPPER_TRIANGULAR>::init_rows() {
    std::cout << " init upper tri \n";
    value_t* pointer = &distances[0];
    std::cout << *pointer <<" ";
    for (index_t i = 0; i < size() - 1; ++i) {
        //        std::cout <<i<<" ";
        rows[i] = pointer;
        std::cout << *pointer <<" ";
        pointer += size() - i - 1;
        //        std::cout <<size() - i - 1<<"\n";
    }
    std::cout <<"\n";
}

template <>
value_t compressed_distance_matrix<UPPER_TRIANGULAR>::operator()(const index_t i,
                                                                 const index_t j) const {
    return i == j ? 0 : i > j ? rows[j][i] : rows[i][j];
}


template <>
value_t compressed_distance_matrix<LOWER_TRIANGULAR>::operator()(const index_t i,
                                                                 const index_t j) const {
    return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

template<>
void compressed_distance_matrix <LOWER_TRIANGULAR>::print() const{
    // printing
    std::cout <<"\n * * * * * * * * * \n";
    for (index_t i = 0; i< size(); ++i){
        for(index_t j = 0; j<= i; ++j)
            if (i==j) {std::cout << 0 <<" "; continue;}
            else{
                std::cout << *(rows[i]+j)<<" ";
            }
        std::cout <<"\n";
    }
    std::cout <<" ** * ** * ** * *\n";
}

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;

class euclidean_distance_matrix {
public:
    std::vector<std::vector<value_t>> points;
    
    euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points)
    : points(std::move(_points)) {
        for (auto p : points) { assert(p.size() == points.front().size()); }
    }
    
    value_t operator()(const index_t i, const index_t j) const {
        assert(i < points.size());
        assert(j < points.size());
        return std::sqrt(std::inner_product(
                                            points[i].begin(), points[i].end(), points[j].begin(), value_t(), std::plus<value_t>(),
                                            [](value_t u, value_t v) { return (u - v) * (u - v); }));
    }
    
    size_t size() const { return points.size(); }
};

class union_find {
    std::vector<index_t> parent;
    std::vector<uint8_t> rank;
    
public:
    union_find(index_t n) : parent(n), rank(n, 0) {
        for (index_t i = 0; i < n; ++i) parent[i] = i;
    }
    
    index_t find(index_t x) {
        index_t y = x, z;
        while ((z = parent[y]) != y) y = z;
        while ((z = parent[x]) != y) {
            parent[x] = y;
            x = z;
        }
        return z;
    }
    void link(index_t x, index_t y) {
        if ((x = find(x)) == (y = find(y))) return;
        if (rank[x] > rank[y])
            parent[y] = x;
        else {
            parent[x] = y;
            if (rank[x] == rank[y]) ++rank[y];
        }
    }
};

template <typename Heap> diameter_entry_t pop_pivot(Heap& column, coefficient_t modulus) {
    if (column.empty())
        return diameter_entry_t(-1);
    else {
        auto pivot = column.top();
        
#ifdef USE_COEFFICIENTS
        coefficient_t coefficient = 0;
        do {
            coefficient = (coefficient + get_coefficient(column.top())) % modulus;
            column.pop();
            
            if (coefficient == 0) {
                if (column.empty())
                    return diameter_entry_t(-1);
                else
                    pivot = column.top();
            }
        } while (!column.empty() && get_index(column.top()) == get_index(pivot));
        if (get_index(pivot) != -1) { set_coefficient(pivot, coefficient); }
#else
        column.pop();
        while (!column.empty() && get_index(column.top()) == get_index(pivot)) {
            column.pop();
            if (column.empty())
                return diameter_entry_t(-1);
            else {
                pivot = column.top();
                column.pop();
            }
        }
#endif
        return pivot;
    }
}

template <typename Heap> diameter_entry_t get_pivot(Heap& column, coefficient_t modulus) {
    diameter_entry_t result = pop_pivot(column, modulus);
    if (get_index(result) != -1) column.push(result);
    return result;
}

template <typename ValueType> class compressed_sparse_matrix {
    std::vector<size_t> bounds;
    std::vector<ValueType> entries;
    
public:
    size_t size() const { return bounds.size(); }
    
    typename std::vector<ValueType>::const_iterator cbegin(size_t index) const {
        assert(index < size());
        return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
    }
    
    typename std::vector<ValueType>::const_iterator cend(size_t index) const {
        assert(index < size());
        return entries.cbegin() + bounds[index];
    }
    
    template <typename Iterator> void append_column(Iterator begin, Iterator end) {
        for (Iterator it = begin; it != end; ++it) { entries.push_back(*it); }
        bounds.push_back(entries.size());
    }
    
    void append_column() { bounds.push_back(entries.size()); }
    
    void push_back(ValueType e) {
        assert(0 < size());
        entries.push_back(e);
        ++bounds.back();
    }
    
    void pop_back() {
        assert(0 < size());
        entries.pop_back();
        --bounds.back();
    }
    
    template <typename Collection> void append_column(const Collection collection) {
        append_column(collection.cbegin(), collection.cend());
    }
};

template <typename DistanceMatrix> class ripser {
    DistanceMatrix dist;
    index_t n, dim_max;
    value_t threshold;
    float ratio;
    coefficient_t modulus;
    const binomial_coeff_table binomial_coeff;
    std::vector<coefficient_t> multiplicative_inverse;
    mutable std::vector<index_t> vertices;
    mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator> neighbor_it;
    mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator> neighbor_end;
    mutable std::vector<diameter_entry_t> coface_entries;
    std::vector<std::vector<std::pair<value_t,value_t>>> pintervals;
public:
    mutable std::vector< std::vector< std::vector<int> > > cocycles_by_dim;
    
    ripser(DistanceMatrix&& _dist, index_t _dim_max, value_t _threshold, float _ratio,
           coefficient_t _modulus)
    : dist(std::move(_dist)), n(dist.size()),
    dim_max(std::min(_dim_max, index_t(dist.size() - 2))), threshold(_threshold),
    ratio(_ratio), modulus(_modulus), binomial_coeff(n, dim_max + 2),
    multiplicative_inverse(multiplicative_inverse_vector(_modulus)) {}
    
    index_t get_next_vertex(index_t& v, const index_t idx, const index_t k) const {
        if (binomial_coeff(v, k) > idx) {
            index_t count = v;
            while (count > 0) {
                index_t i = v;
                index_t step = count >> 1;
                i -= step;
                if (binomial_coeff(i, k) > idx) {
                    v = --i;
                    count -= step + 1;
                } else
                    count = step;
            }
        }
        assert(binomial_coeff(v, k) <= idx && binomial_coeff(v + 1, k) > idx);
        return v;
    }
    
    index_t get_edge_index(const index_t i, const index_t j) const {
        return binomial_coeff(i, 2) + j;
    }
    
    template <typename OutputIterator>
    OutputIterator get_simplex_vertices(index_t idx, const index_t dim, index_t v,
                                        OutputIterator out) const {
        --v;
        for (index_t k = dim + 1; k > 0; --k) {
            get_next_vertex(v, idx, k);
            *out++ = v;
            idx -= binomial_coeff(v, k);
        }
        return out;
    }
    
    // Given dimension, and index (sort of id) of a k-simplex, compute the maximum of the pairwise distance of the vertices in the simplex.
    value_t compute_diameter(const index_t index, index_t dim) const {
        value_t diam = -std::numeric_limits<value_t>::infinity();
        
        vertices.clear();
        get_simplex_vertices(index, dim, dist.size(), std::back_inserter(vertices));
        
        for (index_t i = 0; i <= dim; ++i)
            for (index_t j = 0; j < i; ++j) {
                diam = std::max(diam, dist(vertices[i], vertices[j]));
            }
        return diam;
    }
    
    class simplex_coboundary_enumerator;
    
    void assemble_columns_to_reduce(std::vector<diameter_index_t>& simplices,
                                    std::vector<diameter_index_t>& columns_to_reduce,
                                    hash_map<index_t, index_t>& pivot_column_index, index_t dim);
    
    void compute_dim_0_pairs(std::vector<diameter_index_t>& edges,
                             std::vector<diameter_index_t>& columns_to_reduce,std::ofstream& myfile) {
        
        pintervals.push_back(std::vector<std::pair<value_t,value_t>>());
        union_find dset(n);
        
        edges = get_edges(); // Gather all the edges less than threshold
        std::cout <<" num edges: "<<edges.size()<<"\n";
        std::sort(edges.rbegin(), edges.rend(),
                  greater_diameter_or_smaller_index<diameter_index_t>()); // sort the edges increasing order of their weights (if weight equal smaller index first)
        
#ifdef PRINT_PERSISTENCE_PAIRS
        std::cout << "persistence intervals in dim 0:" << std::endl;
#endif
        
        std::vector<index_t> vertices_of_edge(2);
        for (auto e : edges) {
            vertices_of_edge.clear();
            get_simplex_vertices(get_index(e), 1, n, std::back_inserter(vertices_of_edge));
            index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);
            //            std::cout << "edge: "<<vertices_of_edge[0] <<" "<<vertices_of_edge[1]<< "diameter "<< get_diameter(e)<<"\n";
            if (u != v) {
                //#ifdef PRINT_PERSISTENCE_PAIRS
                if (get_diameter(e) != 0){
                    std::cout << " [0," << get_diameter(e) << ")" << std::endl;
                    //                    myfile << " 0," << get_diameter(e) << "\n" <<std::flush;
                    pintervals[0].push_back(std::make_pair(0,get_diameter(e)));
                }
                //#endif
                dset.link(u, v);
                //                std::cout << u <<" "<<v<<" link \n";
            } else
                columns_to_reduce.push_back(e); // Don't understand
            //            std::cout << " - - - - - -";
        }
        std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());
        
        //#ifdef PRINT_PERSISTENCE_PAIRS
        for (index_t i = 0; i < n; ++i){
            if (dset.find(i) == i){
                std::cout << " 0,Inf\n";
                value_t inf = std::numeric_limits<value_t>::max();
                pintervals[0].push_back(std::make_pair(0,inf));
            }
        }
        //#endif
        //        std::cout << " \n columns to reduce => \n";
        //        for (auto cr:columns_to_reduce){
        //            vertices_of_edge.clear();
        //            get_simplex_vertices(get_index(cr), 1, n, std::back_inserter(vertices_of_edge));
        //            std::cout << vertices_of_edge[0]<<" "<<vertices_of_edge[1]<<"\n";
        //        }
    }
    
    template <typename Column, typename Iterator>
    diameter_entry_t add_coboundary_and_get_pivot(Iterator column_begin, Iterator column_end,
                                                  coefficient_t factor_column_to_add,
#ifdef ASSEMBLE_REDUCTION_MATRIX
                                                  Column& working_reduction_column,
#endif
                                                  Column& working_coboundary, const index_t& dim,
                                                  hash_map<index_t, index_t>& pivot_column_index,
                                                  bool& might_be_apparent_pair);
    
    void compute_pairs(std::vector<diameter_index_t>& columns_to_reduce,
                       hash_map<index_t, index_t>& pivot_column_index, index_t dim, std::ofstream& myfile) {
        pintervals.push_back(std::vector<std::pair<value_t,value_t>>());
#ifdef PRINT_PERSISTENCE_PAIRS
        std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
#endif
        
#ifdef ASSEMBLE_REDUCTION_MATRIX
        compressed_sparse_matrix<diameter_entry_t> reduction_matrix;
#else
#ifdef USE_COEFFICIENTS
        std::vector<diameter_entry_t> reduction_matrix;
#endif
#endif
        
        std::vector<diameter_entry_t> coface_entries;
        
        for (index_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
             ++index_column_to_reduce) {
            auto column_to_reduce = columns_to_reduce[index_column_to_reduce];
            
#ifdef ASSEMBLE_REDUCTION_MATRIX
            std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
            greater_diameter_or_smaller_index<diameter_entry_t>>
            working_reduction_column;
#endif
            
            std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
            greater_diameter_or_smaller_index<diameter_entry_t>>
            working_coboundary;
            
            value_t diameter = get_diameter(column_to_reduce);
            
#ifdef INDICATE_PROGRESS
            if ((index_column_to_reduce + 1) % 1000000 == 0)
                std::cout << "\033[K"
                << "reducing column " << index_column_to_reduce + 1 << "/"
                << columns_to_reduce.size() << " (diameter " << diameter << ")"
                << std::flush << "\r";
#endif
            
            index_t index_column_to_add = index_column_to_reduce;
            
            diameter_entry_t pivot;
            
            // start with factor 1 in order to initialize working_coboundary
            // with the coboundary of the simplex with index column_to_reduce
            coefficient_t factor_column_to_add = 1;
            
#ifdef ASSEMBLE_REDUCTION_MATRIX
            // initialize reduction_matrix as identity matrix
            reduction_matrix.append_column();
#endif
#ifdef USE_COEFFICIENTS
            reduction_matrix.push_back(diameter_entry_t(column_to_reduce, 1));
#endif
            
            bool might_be_apparent_pair = (index_column_to_reduce == index_column_to_add);
            
            while (true) {
#ifdef ASSEMBLE_REDUCTION_MATRIX
#ifdef USE_COEFFICIENTS
                auto reduction_column_begin = reduction_matrix.cbegin(index_column_to_add),
                reduction_column_end = reduction_matrix.cend(index_column_to_add);
#else
                std::vector<diameter_entry_t> coeffs;
                coeffs.push_back(columns_to_reduce[index_column_to_add]);
                for (auto it = reduction_matrix.cbegin(index_column_to_add);
                     it != reduction_matrix.cend(index_column_to_add); ++it)
                    coeffs.push_back(*it);
                auto reduction_column_begin = coeffs.begin(), reduction_column_end = coeffs.end();
#endif
#else
#ifdef USE_COEFFICIENTS
                auto reduction_column_begin = &reduction_matrix[index_column_to_add],
                reduction_column_end = &reduction_matrix[index_column_to_add] + 1;
#else
                auto reduction_column_begin = &columns_to_reduce[index_column_to_add],
                reduction_column_end = &columns_to_reduce[index_column_to_add] + 1;
#endif
#endif
                
                pivot = add_coboundary_and_get_pivot(
                                                     reduction_column_begin, reduction_column_end, factor_column_to_add,
#ifdef ASSEMBLE_REDUCTION_MATRIX
                                                     working_reduction_column,
#endif
                                                     working_coboundary, dim, pivot_column_index, might_be_apparent_pair);
                
                if (get_index(pivot) != -1) {
                    auto pair = pivot_column_index.find(get_index(pivot));
                    
                    if (pair != pivot_column_index.end()) {
                        index_column_to_add = pair->second;
                        factor_column_to_add = modulus - get_coefficient(pivot);
                    } else {
#ifdef PRINT_PERSISTENCE_PAIRS
                        value_t death = get_diameter(pivot);
                        if (death > diameter * ratio) {
#ifdef INDICATE_PROGRESS
                            std::cout << "\033[K";
#endif
                            //	std::cout << " [" << diameter << "," << death << ")" << std::endl
                            //	          << std::flush;
                            //myfile  << diameter << "," << death << std::endl << std::flush;
                            pintervals[dim].push_back(std::make_pair(diameter,death));
#ifdef ASSEMBLE_REDUCTION_MATRIX
                            
                            //Representative cocycle
                            auto cocycle = working_reduction_column;
                            diameter_entry_t e;
                            std::vector<index_t> simplex;
                            std::vector<int> thiscocycle;
                            while (get_index(e = get_pivot(cocycle, modulus)) != -1) {
                                simplex.clear();
                                get_simplex_vertices(get_index(e), dim, n, std::back_inserter(simplex));
                                for (size_t k = 0; k < simplex.size(); k++) {
                                    thiscocycle.push_back((int)simplex[k]);
                                }
                                thiscocycle.push_back(normalize(get_coefficient(e), modulus));
                                cocycle.pop();
                            }
                            //                            std::cout << " basis: ";
                            //                            for(auto v: thiscocycle){
                            //                                std::cout<<v<<" ";
                            //                            }
                            //                            std::cout <<std::endl << std::flush;
                            cocycles_by_dim[dim].push_back(thiscocycle);
#endif
                        }
                        
#endif
                        
                        
                        pivot_column_index.insert(
                                                  std::make_pair(get_index(pivot), index_column_to_reduce));
                        
#ifdef USE_COEFFICIENTS
                        const coefficient_t inverse =
                        multiplicative_inverse[get_coefficient(pivot)];
#endif
                        
#ifdef ASSEMBLE_REDUCTION_MATRIX
                        // replace current column of reduction_matrix (with a single diagonal 1
                        // entry) by reduction_column (possibly with a different entry on the
                        // diagonal)
#ifdef USE_COEFFICIENTS
                        reduction_matrix.pop_back();
#else
                        pop_pivot(working_reduction_column, modulus);
#endif
                        
                        while (true) {
                            diameter_entry_t e = pop_pivot(working_reduction_column, modulus);
                            if (get_index(e) == -1) break;
#ifdef USE_COEFFICIENTS
                            set_coefficient(e, inverse * get_coefficient(e) % modulus);
                            assert(get_coefficient(e) > 0);
#endif
                            reduction_matrix.push_back(e);
                        }
#else
#ifdef USE_COEFFICIENTS
                        reduction_matrix.pop_back();
                        reduction_matrix.push_back(diameter_entry_t(column_to_reduce, inverse));
#endif
#endif
                        break;
                    }
                } else {
#ifdef PRINT_PERSISTENCE_PAIRS
                    //					std::cout << " [" << diameter << ", )" << std::endl << std::flush;
                    //myfile << diameter << ", Inf" << std::endl << std::flush;
#endif
                    value_t inf = std::numeric_limits<value_t>::max();
                    pintervals[dim].push_back(std::make_pair(diameter,inf));
#ifdef ASSEMBLE_REDUCTION_MATRIX
                    
                    //Representative cocycle
                    auto cocycle = working_reduction_column;
                    diameter_entry_t e;
                    std::vector<index_t> simplex;
                    std::vector<int> thiscocycle;
                    while (get_index(e = get_pivot(cocycle, modulus)) != -1) {
                        simplex.clear();
                        get_simplex_vertices(get_index(e), dim, n, std::back_inserter(simplex));
                        for (size_t k = 0; k < simplex.size(); k++) {
                            thiscocycle.push_back((int)simplex[k]);
                        }
                        thiscocycle.push_back(normalize(get_coefficient(e), modulus));
                        cocycle.pop();
                    }
                    //                    std::cout << " basis: ";
                    //                    for(auto v:thiscocycle){
                    //                        std::cout<<v<<" ";
                    //                    }
                    //                    std::cout <<std::endl << std::flush;
                    cocycles_by_dim[dim].push_back(thiscocycle);
#endif
                    break;
                }
            }
        }
        
#ifdef INDICATE_PROGRESS
        std::cout << "\033[K";
#endif
    }
    
    std::vector<diameter_index_t> get_edges();
    
    std::vector<std::vector<std::pair<value_t,value_t>>>
    compute_barcodes(std::string filename) {
        // Interval Files will be stored as distancefilename_0.csv or _1.csv
        //        std::cout << "Class type = "<<typeid(dist).name()<<"\n";
        //        dist.print();
        std::vector<diameter_index_t> simplices, columns_to_reduce;
        cocycles_by_dim.push_back(std::vector< std::vector<int> >());
        std::ofstream myfile;
        std::cout<<filename+"_0.csv"<<"\n";
        //        std::cout << " dd "<< simplices.size() << " " << columns_to_reduce.size()<< " \n";
        // myfile.open(filename+cluster_id+"_0.csv"); // Write 0-th persistent intervals
        {
            compute_dim_0_pairs(simplices, columns_to_reduce,myfile); // This function populates columns_to_reduce with edges for whom (0,inf) created in Dim 0
        }
        //        std::cout << "\n dd"<< simplices.size()<< " "<< columns_to_reduce.size()<< " sz\n";
        //        std::cout << "dimmax = "<<dim_max<<"\n";
        // myfile.close();
        for (index_t dim = 1; dim <= dim_max; ++dim) {
            //            std::cout << "ph in dim: "<<dim<<"\n";
            hash_map<index_t, index_t> pivot_column_index;
            cocycles_by_dim.push_back(std::vector< std::vector<int> >());
            pivot_column_index.reserve(columns_to_reduce.size());
            std::cout<<filename+"_"+std::to_string(dim)+".csv"<<"\n";
            //   myfile.open(filename+cluster_id+"_"+std::to_string(dim)+".csv");
            {
                compute_pairs(columns_to_reduce, pivot_column_index, dim,myfile);
            }
            // myfile.close();
            if (dim < dim_max) {
                assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
                                           dim + 1);
            }
        }
        return pintervals;
    }
};

template <> class ripser<compressed_lower_distance_matrix>::simplex_coboundary_enumerator {
private:
    index_t idx_below, idx_above, v, k;
    std::vector<index_t> vertices;
    const diameter_entry_t simplex;
    const coefficient_t modulus;
    const compressed_lower_distance_matrix& dist;
    const binomial_coeff_table& binomial_coeff;
    
public:
    simplex_coboundary_enumerator(const diameter_entry_t _simplex, index_t _dim,
                                  const ripser& parent)
    : idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1), k(_dim + 1),
    vertices(_dim + 1), simplex(_simplex), modulus(parent.modulus), dist(parent.dist),
    binomial_coeff(parent.binomial_coeff) {
        parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n, vertices.begin());
    }
    
    bool has_next() {
        while ((v != -1) && (binomial_coeff(v, k) <= idx_below)) {
            idx_below -= binomial_coeff(v, k);
            idx_above += binomial_coeff(v, k + 1);
            --v;
            --k;
            assert(k != -1);
        }
        return v != -1;
    }
    
    diameter_entry_t next() {
        value_t coface_diameter = get_diameter(simplex);
        for (index_t w : vertices) coface_diameter = std::max(coface_diameter, dist(v, w));
        index_t coface_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
        coefficient_t coface_coefficient =
        (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;
        return diameter_entry_t(coface_diameter, coface_index, coface_coefficient);
    }
};

template <typename DistanceMatrix>
template <typename Column, typename Iterator>
diameter_entry_t ripser<DistanceMatrix>::add_coboundary_and_get_pivot(
                                                                      Iterator column_begin, Iterator column_end, coefficient_t factor_column_to_add,
#ifdef ASSEMBLE_REDUCTION_MATRIX
                                                                      Column& working_reduction_column,
#endif
                                                                      Column& working_coboundary, const index_t& dim, hash_map<index_t, index_t>& pivot_column_index,
                                                                      bool& might_be_apparent_pair) {
    for (auto it = column_begin; it != column_end; ++it) {
        diameter_entry_t simplex = *it;
        set_coefficient(simplex, get_coefficient(simplex) * factor_column_to_add % modulus);
        
#ifdef ASSEMBLE_REDUCTION_MATRIX
        working_reduction_column.push(simplex);
#endif
        
        coface_entries.clear();
        simplex_coboundary_enumerator cofaces(simplex, dim, *this);
        while (cofaces.has_next()) {
            diameter_entry_t coface = cofaces.next();
            if (get_diameter(coface) <= threshold) {
                coface_entries.push_back(coface);
                if (might_be_apparent_pair && (get_diameter(simplex) == get_diameter(coface))) {
                    if (pivot_column_index.find(get_index(coface)) == pivot_column_index.end()) {
                        return coface;
                    }
                    might_be_apparent_pair = false;
                }
            }
        }
        for (auto coface : coface_entries) working_coboundary.push(coface);
    }
    
    return get_pivot(working_coboundary, modulus);
}

template <> class ripser<sparse_distance_matrix>::simplex_coboundary_enumerator {
private:
    const ripser& parent;
    
    index_t idx_below, idx_above, v, k, max_vertex_below;
    const diameter_entry_t simplex;
    const coefficient_t modulus;
    const sparse_distance_matrix& dist;
    const binomial_coeff_table& binomial_coeff;
    
    std::vector<index_t>& vertices;
    std::vector<std::vector<index_diameter_t>::const_reverse_iterator>& neighbor_it;
    std::vector<std::vector<index_diameter_t>::const_reverse_iterator>& neighbor_end;
    index_diameter_t x;
    
public:
    simplex_coboundary_enumerator(const diameter_entry_t _simplex, index_t _dim,
                                  const ripser& _parent)
    : parent(_parent), idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1),
    k(_dim + 1), max_vertex_below(parent.n - 1), simplex(_simplex), modulus(parent.modulus),
    dist(parent.dist), binomial_coeff(parent.binomial_coeff), vertices(parent.vertices),
    neighbor_it(parent.neighbor_it), neighbor_end(parent.neighbor_end) {
        
        neighbor_it.clear();
        neighbor_end.clear();
        vertices.clear();
        
        parent.get_simplex_vertices(idx_below, _dim, parent.n, std::back_inserter(vertices));
        
        for (auto v : vertices) {
            neighbor_it.push_back(dist.neighbors[v].rbegin());
            neighbor_end.push_back(dist.neighbors[v].rend());
        }
    }
    
    bool has_next(bool all_cofaces = true) {
        for (auto &it0 = neighbor_it[0], &end0 = neighbor_end[0]; ; ++it0) {
            // != operator was overloaded in snap. Hence I changed it.
            if(it0 == end0) break;
            x = *it0;
            for (size_t idx = 1; idx < neighbor_it.size(); ++idx) {
                auto &it = neighbor_it[idx], end = neighbor_end[idx];
                while (get_index(*it) > get_index(x))
                    if (++it == end) return false;
                auto y = *it;
                if (get_index(y) != get_index(x))
                    goto continue_outer;
                else
                    x = std::max(x, y);
            }
            return all_cofaces || !(k > 0 && parent.get_next_vertex(max_vertex_below, idx_below,
                                                                    k) > get_index(x));
        continue_outer:;
        }
        return false;
    }
    
    diameter_entry_t next() {
        ++neighbor_it[0];
        
        while (k > 0 && parent.get_next_vertex(max_vertex_below, idx_below, k) > get_index(x)) {
            idx_below -= binomial_coeff(max_vertex_below, k);
            idx_above += binomial_coeff(max_vertex_below, k + 1);
            --k;
        }
        
        value_t coface_diameter = std::max(get_diameter(simplex), get_diameter(x));
        
        coefficient_t coface_coefficient =
        (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;
        
        return diameter_entry_t(coface_diameter,
                                idx_above + binomial_coeff(get_index(x), k + 1) + idx_below,
                                coface_coefficient);
    }
};

// return set of edges whose weight/(end-point distances) less than threshold
template <> std::vector<diameter_index_t> ripser<compressed_lower_distance_matrix>::get_edges() {
    std::vector<diameter_index_t> edges;
    for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
        value_t diameter = compute_diameter(index, 1);
        if (diameter <= threshold) edges.push_back(std::make_pair(diameter, index));
    }
    return edges;
}

// return set of edges whose weight/(end-point distances) less than threshold
template <> std::vector<diameter_index_t> ripser<sparse_distance_matrix>::get_edges() {
    std::vector<diameter_index_t> edges;
    for (index_t i = 0; i < n; ++i)
        for (auto n : dist.neighbors[i]) {
            index_t j = get_index(n);
            if (i > j) edges.push_back(std::make_pair(get_diameter(n), get_edge_index(i, j)));
        }
    return edges;
}

template <>
void ripser<compressed_lower_distance_matrix>::assemble_columns_to_reduce(
                                                                          std::vector<diameter_index_t>& simplices, std::vector<diameter_index_t>& columns_to_reduce,
                                                                          hash_map<index_t, index_t>& pivot_column_index, index_t dim) {
    index_t num_simplices = binomial_coeff(n, dim + 1);
    
    columns_to_reduce.clear();
    
#ifdef INDICATE_PROGRESS
    std::cout << "\033[K"
    << "assembling " << num_simplices << " columns" << std::flush << "\r";
#endif
    
    for (index_t index = 0; index < num_simplices; ++index) {
        if (pivot_column_index.find(index) == pivot_column_index.end()) {
            value_t diameter = compute_diameter(index, dim);
            if (diameter <= threshold) columns_to_reduce.push_back(std::make_pair(diameter, index));
#ifdef INDICATE_PROGRESS
            if ((index + 1) % 1000000 == 0)
                std::cout << "\033[K"
                << "assembled " << columns_to_reduce.size() << " out of " << (index + 1)
                << "/" << num_simplices << " columns" << std::flush << "\r";
#endif
        }
    }
    
#ifdef INDICATE_PROGRESS
    std::cout << "\033[K"
    << "sorting " << num_simplices << " columns" << std::flush << "\r";
#endif
    
    std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              greater_diameter_or_smaller_index<diameter_index_t>());
#ifdef INDICATE_PROGRESS
    std::cout << "\033[K";
#endif
}

template <>
void ripser<sparse_distance_matrix>::assemble_columns_to_reduce(
                                                                std::vector<diameter_index_t>& simplices, std::vector<diameter_index_t>& columns_to_reduce,
                                                                hash_map<index_t, index_t>& pivot_column_index, index_t dim) {
    
#ifdef INDICATE_PROGRESS
    std::cout << "\033[K"
    << "assembling columns" << std::flush << "\r";
#endif
    
    --dim;
    columns_to_reduce.clear();
    
    std::vector<diameter_index_t> next_simplices;
    
    for (diameter_index_t simplex : simplices) {
        simplex_coboundary_enumerator cofaces(simplex, dim, *this);
        
        while (cofaces.has_next(false)) {
            auto coface = cofaces.next();
            
            next_simplices.push_back(std::make_pair(get_diameter(coface), get_index(coface)));
            
            if (pivot_column_index.find(get_index(coface)) == pivot_column_index.end())
                columns_to_reduce.push_back(
                                            std::make_pair(get_diameter(coface), get_index(coface)));
        }
    }
    
    simplices.swap(next_simplices);
    
#ifdef INDICATE_PROGRESS
    std::cout << "\033[K"
    << "sorting " << columns_to_reduce.size() << " columns" << std::flush << "\r";
#endif
    
    std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              greater_diameter_or_smaller_index<diameter_index_t>());
#ifdef INDICATE_PROGRESS
    std::cout << "\033[K";
#endif
}

enum file_format {
    LOWER_DISTANCE_MATRIX,
    UPPER_DISTANCE_MATRIX,
    DISTANCE_MATRIX,
    POINT_CLOUD,
    DIPHA,
    SPARSE,
    RIPSER,
    GRAPH,
    WGRAPH
};

template <typename T> T read(std::istream& s) {
    T result;
    s.read(reinterpret_cast<char*>(&result), sizeof(T));
    return result; // on little endian: boost::endian::little_to_native(result);
}

compressed_lower_distance_matrix read_point_cloud(std::istream& input_stream) {
    std::vector<std::vector<value_t>> points;
    
    std::string line;
    value_t value;
    while (std::getline(input_stream, line)) {
        std::vector<value_t> point;
        std::istringstream s(line);
        while (s >> value) {
            point.push_back(value);
            s.ignore();
        }
        if (!point.empty()) points.push_back(point);
        assert(point.size() == points.front().size());
    }
    
    euclidean_distance_matrix eucl_dist(std::move(points));
    
    index_t n = eucl_dist.size();
    
    std::cout << "point cloud with " << n << " points in dimension "
    << eucl_dist.points.front().size() << std::endl;
    
    std::vector<value_t> distances;
    
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));
    
    return compressed_lower_distance_matrix(std::move(distances));
}

sparse_distance_matrix read_sparse_distance_matrix(std::istream& input_stream) {
    
    std::vector<std::vector<index_diameter_t>> neighbors;
    
    index_t num_edges = 0;
    
    std::string line;
    while (std::getline(input_stream, line)) {
        std::istringstream s(line);
        size_t i, j;
        value_t value;
        s >> i;
        s >> j;
        s >> value;
        if (i != j) {
            neighbors.resize(std::max({neighbors.size(), i + 1, j + 1}));
            neighbors[i].push_back(std::make_pair(j, value));
            neighbors[j].push_back(std::make_pair(i, value));
            ++num_edges;
        }
    }
    
    for (index_t i = 0; i < neighbors.size(); ++i)
        std::sort(neighbors[i].begin(), neighbors[i].end());
    
    return sparse_distance_matrix(std::move(neighbors), num_edges);
}

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
    std::vector<value_t> distances;
    value_t value;
    while (input_stream >> value) {
        distances.push_back(value);
        input_stream.ignore();
    }
    
    return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_upper_distance_matrix(std::istream& input_stream) {
    std::vector<value_t> distances;
    value_t value;
    while (input_stream >> value) {
        distances.push_back(value);
        input_stream.ignore();
    }
    
    return compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
}

compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
    std::vector<value_t> distances;
    
    std::string line;
    value_t value;
    for (int i = 0; std::getline(input_stream, line); ++i) {
        std::istringstream s(line);
        for (int j = 0; j < i && s >> value; ++j) {
            distances.push_back(value);
            s.ignore();
        }
    }
    
    return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_dipha(std::istream& input_stream) {
    if (read<int64_t>(input_stream) != 8067171840) {
        std::cerr << "input is not a Dipha file (magic number: 8067171840)" << std::endl;
        exit(-1);
    }
    
    if (read<int64_t>(input_stream) != 7) {
        std::cerr << "input is not a Dipha distance matrix (file type: 7)" << std::endl;
        exit(-1);
    }
    
    index_t n = read<int64_t>(input_stream);
    
    std::vector<value_t> distances;
    
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i > j)
                distances.push_back(read<double>(input_stream));
            else
                read<double>(input_stream);
    
    return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_ripser(std::istream& input_stream) {
    std::vector<value_t> distances;
    while (!input_stream.eof()) distances.push_back(read<value_t>(input_stream));
    return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream, file_format format) {
    switch (format) {
        case LOWER_DISTANCE_MATRIX:
            return read_lower_distance_matrix(input_stream);
        case UPPER_DISTANCE_MATRIX:
            return read_upper_distance_matrix(input_stream);
        case DISTANCE_MATRIX:
            return read_distance_matrix(input_stream);
        case POINT_CLOUD:
            return read_point_cloud(input_stream);
        case DIPHA:
            return read_dipha(input_stream);
        default:
            return read_ripser(input_stream);
    }
}

// Read weighted undirected graph, and store it as weighted directed network.
template <class nodewType,class edgewType> // for now nodewType=int , edgewType = double for compatibility
TPt <TNodeEDatNet<nodewType, edgewType> > readGraph(std::istream& is)
{
    
    //size_t vertices, edges;
    
    //if (is >> vertices >> edges)
    //{
    //is.ignore(1024, '\n'); // assume not overly long lines
    std::string line;
    TPt <TNodeEDatNet<nodewType, edgewType> > Net = TNodeEDatNet<nodewType, edgewType>::New();
    //while (edges--)
    //{
    //
    value_t min_cost = 10000000;
    for (int i = 0; std::getline(is, line); ++i) {
        std::istringstream s(line);
        //if (is && std::getline(is, line))
        //{
        int from, to;
        value_t cost;
        s>>from;
        s.ignore();
        s>>to;
        s.ignore();
        s>>cost;
        s.ignore();
        //std::istringstream line_stream(line);
        
        
        //std::cout << from << " "<<to<<" "<<cost<<std::endl;
        //if (!(Net.insert({ { from, to }, cost }).second))
        //  throw std::runtime_error("Duplicate edge");
        if(!Net->IsNode(from))
            Net->AddNode(from);
        if(!Net->IsNode(to))
            Net->AddNode(to);
        // add bidirectional edges
        Net->AddEdge(from,to);
        Net->AddEdge(to,from);
        Net->SetEDat(from,to,cost);
        Net->SetEDat(to,from,cost);
        min_cost = std::min(min_cost,cost);
        //}
    }
    std::cout << "Min edge weight= "<<min_cost<<"\n";
    
    return Net;
}


void print_usage_and_exit(int exit_code) {
    std::cerr
    << "Usage: "
    << "ripser "
    << "[options] [filename]" << std::endl
    << std::endl
    << "Options:" << std::endl
    << std::endl
    << "  --help           print this screen" << std::endl
    << "  --format         use the specified file format for the input. Options are:"
    << std::endl
    << "                     lower-distance (lower triangular distance matrix; default)"
    << std::endl
    << "                     upper-distance (upper triangular distance matrix)" << std::endl
    << "                     distance       (full distance matrix)" << std::endl
    << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
    << "                     dipha          (distance matrix in DIPHA file format)" << std::endl
    << "                     sparse         (sparse distance matrix in Sparse Triplet format)"
    << std::endl
    << "                     ripser         (distance matrix in Ripser binary file format)"
    << std::endl
    << "  --dim <k>        compute persistent homology up to dimension <k>" << std::endl
    << "  --threshold <t>  compute Rips complexes up to diameter <t>" << std::endl
#ifdef USE_COEFFICIENTS
    << "  --modulus <p>    compute homology with coefficients in the prime field Z/<p>Z"
#endif
    << std::endl;
    
    exit(exit_code);
}


template <class nodewType,class edgewType>
void Compute_EigenVectorCentr(TPt <TNodeEDatNet<nodewType, edgewType> > Graph, TIntFltH& NIdEigenH, const double& Eps, const int& MaxIter) {
    const int NNodes = Graph->GetNodes();
    NIdEigenH.Gen(NNodes);
    // initialize vector values
    for (auto NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
        NIdEigenH.AddDat(NI.GetId(), 1.0/NNodes);
        IAssert(NI.GetId() == NIdEigenH.GetKey(NIdEigenH.Len()-1));
    }
    TFltV TmpV(NNodes);
    for (int iter = 0; iter < MaxIter; iter++) {
        int j = 0;
        // add neighbor values
        for (auto NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
            TmpV[j] = 0;
            for (int e = 0; e < NI.GetOutDeg(); e++) {
                // Unweighted Eigenvector
                //TmpV[j] += NIdEigenH.GetDat(NI.GetOutNId(e)) ;
                
                // Weighted eigenvector
                TmpV[j] += (NIdEigenH.GetDat(NI.GetOutNId(e)) * Graph->GetEDat(NI.GetId(),NI.GetOutNId(e)));
            }
        }
        
        // normalize
        double sum = 0;
        for (int i = 0; i < TmpV.Len(); i++) {
            sum += (TmpV[i]*TmpV[i]);
        }
        sum = sqrt(sum);
        for (int i = 0; i < TmpV.Len(); i++) {
            TmpV[i] /= sum;
        }
        
        // compute difference
        double diff = 0.0;
        j = 0;
        for (auto NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
            diff += fabs(NIdEigenH.GetDat(NI.GetId())-TmpV[j]);
        }
        
        
        if (diff < Eps) {
            break;
        }
        else{
            // set new values
            j = 0;
            for (auto NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
                NIdEigenH.AddDat(NI.GetId(), TmpV[j]);
            }
        }
    }
}

void run_partialBFS_wgraph(const TPt<TNodeEDatNet<TInt, TFlt>>& Net, bool InitBigQ, int StartNId, float epsilon, std::set <TInt> &candidates, std::set <TInt> &candidates_gteps_lt2eps, hash_map <int, bool> &marked , int &num_markedv, hash_map <int, double> &dist_to_cover){
    // Append a set of new candidates for landmark selection in the next iteration.
    int v;
    TSnapQueue<int> Queue(InitBigQ?Net->GetNodes():1024);
    TIntFltH NIdDistH(InitBigQ?Net->GetNodes():1024);
    
    // Initializations
    IAssert(Net->IsNode(StartNId));
    dist_to_cover[StartNId] = 0;
    NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
    Queue.Push(StartNId); // BFS traversal queue
    candidates.erase(StartNId);
    //    std::cout << "Partial BFS from node "<<StartNId<<"\n";
    while (! Queue.Empty()){
        const int NId = Queue.Top();  Queue.Pop();
        const int Dist = NIdDistH.GetDat(NId);
        auto NodeI = Net->GetNI(NId);
        //        std::cout <<NId<<" ";
        for (v = 0; v < NodeI.GetOutDeg(); v++) {
            const int DstNId = NodeI.GetOutNId(v);
            bool DstNId_marked = marked[DstNId];
            
            if (!NIdDistH.IsKey(DstNId)){
                auto weight = Net->GetEDat(NId,DstNId);
                NIdDistH.AddDat(DstNId, Dist+weight);
                if (Dist + weight <= 2*epsilon){
                    Queue.Push(DstNId);
                    if (Dist + weight <= epsilon){
                        candidates.erase(DstNId);
                        if(!DstNId_marked){
                            marked[DstNId] = true;
                            num_markedv++;
                            candidates_gteps_lt2eps.erase(DstNId);
                        }
                        //                        //Compare! If (v'(key) -> d_l1) was shortest before, now if v' -> d_l2 is shorter, need to update.
                        if (dist_to_cover.find(DstNId)==dist_to_cover.end()) // does not exists
                            dist_to_cover[DstNId] = Dist + weight;
                        else{
                            //                            std::cout << "collision\n";
                            dist_to_cover[DstNId] = min(Dist + weight ,dist_to_cover[DstNId]);
                        }
                        
                        //                        std::cout << " Marked "<<DstNId << "\n";
                    }
                    else{
                        if(!DstNId_marked && candidates_gteps_lt2eps.find(DstNId) == candidates_gteps_lt2eps.end()) // unmarked and not added already
                            candidates_gteps_lt2eps.insert(DstNId);
                    }
                }
                else{ // Dist + 1> 2*epsilon. Add DstNId as candidate
                    if(!DstNId_marked)
                        candidates.insert(DstNId);
                }
            }
        }
    }
}

template <class nodewType,class edgewType>
TIntV weighted_select_landmarks(const TPt<TNodeEDatNet<nodewType, edgewType>>& Net, const std::string criteria, hash_map <int, double> &dist_to_cover, int num_landmarks = 1, float epsilon = 0.01){
    string weight_means = "cost"; // distance/cost/length = reciprocal of similarity/closeness
    TIntV selected_vertices;
    selected_vertices.Reserve(num_landmarks);
    if (criteria.compare("degree")==0){
        std:: cout << "degree selection" <<std::endl;
        std:: vector <std::pair <double,int>> vec_d_intpair;
        for (auto NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
            const int NId = NI.GetId();
            double DegCentr = 0;
            for (int v = 0; v < NI.GetOutDeg(); v++) {
                const int nbr_id = NI.GetOutNId(v);
                if (weight_means == "cost")
                    DegCentr += Net->GetEDat(NId,nbr_id);
                else
                    DegCentr += (1.0/Net->GetEDat(NId,nbr_id)); // closeness
            }
            
            vec_d_intpair.push_back(std::make_pair(DegCentr,NId));
        }
        // weight_means = "cost". Nodes that have higher total cost=> most similar
        std::sort(vec_d_intpair.rbegin(), vec_d_intpair.rend()); // descending order, matlab version
        // in my opinion it should be ->
        //     sort(vec_d_intpair.begin(), vec_d_intpair.end());
        
        int N_ = Net->GetNodes();
        int max_K = TMath::Mn(num_landmarks, N_);
        
        for (int tries = 0; tries < max_K; tries++){
            selected_vertices.Add(vec_d_intpair[tries].second);
            //      std:: cout << vec_d_intpair[tries].first << " "<< vec_d_intpair[tries].second<<std::endl;
        }
    }
    if (criteria.compare("eigen")==0){
        
        TIntFltH EigH;
        double epsilon = 0.1;
        int max_num_iteration = 1000;
        //        Compute_EigenVectorCentr<TInt,TFlt>(Net,EigH,epsilon,max_num_iteration);
        std:: vector <std::pair <double,int>> vec_d_intpair;
        for (auto NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
            const TInt NId = NI.GetId();
            double EigCentr = EigH.GetDat(NId);
            vec_d_intpair.push_back(std::make_pair(EigCentr,NId));
        }
        sort(vec_d_intpair.rbegin(), vec_d_intpair.rend());
        int N_ = Net->GetNodes();
        int max_K = TMath::Mn(num_landmarks, N_);
        
        for (int tries = 0; tries < max_K; tries++){
            selected_vertices.Add(vec_d_intpair[tries].second);
        }
    }
    if (criteria.compare("maxmin")==0){
        // The code for point-cloud needs to be modified
        selected_vertices.Reserve(num_landmarks);
        int N_ =  Net->GetNodes();
        int max_K = TMath::Mn(num_landmarks,N_);
        int StartNId;
        int NextStartNId;
        TIntV NodeIdV;
        Net->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        dist_to_cover.reserve(N_);
        
        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++)
            dist_to_cover[NI.GetId()] = std::numeric_limits<value_t>::infinity();
        
        for (int tries = 0; tries < max_K; tries++) {
            if (tries == 0){
                StartNId = NodeIdV[tries];
            }
            TBreathFSW<TPt<TNodeEDatNet<nodewType, edgewType>>> BFS(Net,true,false);
            BFS.DoBfs(StartNId, true, false, -1, TInt::Mx);
            value_t cur_max = -std::numeric_limits<value_t>::infinity();
            for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){ // Update min dist to landmark set
                int DstNId = NI.GetId();
                double dist = BFS.GetHops(StartNId,DstNId);
                dist_to_cover[DstNId] = min ( dist_to_cover[DstNId], dist );
                if (cur_max < dist_to_cover[DstNId]){
                    NextStartNId = DstNId;
                    cur_max = dist_to_cover[DstNId];
                }
            }
            selected_vertices.Add(StartNId);
            StartNId = NextStartNId;
        }
        
    }
    
    if (criteria.compare("epsmaxmin")==0){
        // The code for point-cloud needs to be modified
        int N_ =  Net->GetNodes();
        int StartNId;
        int NextStartNId;
        TIntV NodeIdV;
        Net->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        dist_to_cover.reserve(N_);
        hash_map <int,bool> marked;
        marked.reserve(N_);
        
        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++)
            dist_to_cover[NI.GetId()] = std::numeric_limits<int>::max();
        
        int num_markedv = 0;
        for (int tries = 0; num_markedv<N_; tries++) {
            if (tries == 0){
                StartNId = NodeIdV[tries];
                // Run eps-BFS and mark its cover
                
            }
            //            TBreathFS<PUNGraph> BFS(g);
            //            BFS.DoBfs(StartNId, true, false, -1, TInt::Mx);
            
            bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
            int v;
            TSnapQueue<int> Queue(InitBigQ?Net->GetNodes():1024);
            TIntFltH NIdDistH(InitBigQ?Net->GetNodes():1024);
            
            // run bfs on whole graph
            // Initializations
            IAssert(Net->IsNode(StartNId));
            NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
            Queue.Push(StartNId); // BFS traversal queue
            marked[StartNId] = true;
            num_markedv++;
            // Run epsilon-BFS
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                auto Dist = NIdDistH.GetDat(NId);
                auto NodeI = Net->GetNI(NId);
                for (v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    if (!NIdDistH.IsKey(DstNId)){
                        if (Dist + Net->GetEDat(NId,DstNId) <= epsilon && !marked[DstNId]){
                            num_markedv++;
                            marked[DstNId] = true;
                        }
                        Queue.Push(DstNId);
                        NIdDistH.AddDat(DstNId, Dist + Net->GetEDat(NId,DstNId) );
                    }
                }
            }
            
            value_t cur_max = std::numeric_limits<int>::min();
            for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){ // Update min dist to landmark set
                int DstNId = NI.GetId();
                double dist = NIdDistH.GetDat(DstNId);
                dist_to_cover[DstNId] = min ( dist_to_cover[DstNId], dist );
                if (cur_max < dist_to_cover[DstNId]){
                    NextStartNId = DstNId;
                    cur_max = dist_to_cover[DstNId];
                }
            }
            selected_vertices.Add(StartNId);
            StartNId = NextStartNId;
        }
        
    }
    
    if (criteria.compare("epsrand")==0){
        std::cout << "epsrand\n";
        vector <TInt> landmarks;
        std::default_random_engine generator; // random generator
        
        // Select first landmark
        TIntV NodeIdV;
        int N_ = Net->GetNodes();
        Net->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        
        int num_markedv = 0;
        int StartNId = NodeIdV[0];
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        TIntV vertices(N_);
        Net->GetNIdV(vertices);
        std::set <TInt> candidates;
        for (int i=0; i<vertices.Len();i++)   candidates.insert(vertices[i]);
        while (num_markedv< N_){
            //            std::cout <<" Before marking: "<<num_markedv<<" ";
            //            int temp_ = num_markedv;
            if (!marked[StartNId]){
                marked[StartNId] = true;
                num_markedv++;
            }
            //            std::cout << "marked # "<<num_markedv<<" total # "<<N_<<"\n";
            // run partial bfs and mark
            // BFS house-keeping data structures
            bool InitBigQ = N_>10*1024? true: false;; // Whether to initialize large sized queue or hashtable.
            
            
            int v;
            TSnapQueue<int> Queue(InitBigQ?Net->GetNodes():1024);
            TIntFltH NIdDistH(InitBigQ?Net->GetNodes():1024);
            //            hash_map <int,int> NIdDistH;
            
            // Initializations
            IAssert(Net->IsNode(StartNId));
            //            NIdDistH[StartNId] = 0;
            candidates.erase(StartNId);
            landmarks.push_back(StartNId);
            NIdDistH.AddDat(StartNId, 0.0); // keep track of shortest path distances
            dist_to_cover.insert({StartNId,0.0});
            Queue.Push(StartNId); // BFS traversal queue
            //    marked[StartNId] = true; // marked flag
            //    num_markedv++;
            
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                auto Dist = NIdDistH.GetDat(NId);
                //                const int Dist = NIdDistH[NId];
                auto NodeI = Net->GetNI(NId);
                //        std::cout <<NId<<" ";
                for (v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    //            if (! marked[DstNId]) {
                    if (!NIdDistH.IsKey(DstNId)){
                        //                    if (NIdDistH.find(DstNId)==NIdDistH.end()){ // DstNId key does not exists
                        NIdDistH.AddDat(DstNId, Dist+ Net->GetEDat(NId,DstNId) );
                        //                          NIdDistH[DstNId] = Dist + 1;
                        if (Dist + Net->GetEDat(NId,DstNId) <= epsilon){
                            Queue.Push(DstNId);
                            if (dist_to_cover.find(DstNId)==dist_to_cover.end()) // does not exists
                                dist_to_cover[DstNId] = Dist + Net->GetEDat(NId,DstNId);
                            else
                                dist_to_cover[DstNId] = min(Dist + Net->GetEDat(NId,DstNId),dist_to_cover[DstNId]);
                            if(!marked[DstNId]){
                                marked[DstNId] = true;
                                num_markedv++;
                                int n_removed = candidates.erase(DstNId);
                                assert(n_removed == 1);
                            }
                            //                    std::cout << " Marked "<<DstNId << "\n";
                        }
                    }
                }
                
            }
            //            std::cout <<" After : "<<num_markedv<<" Diff: "<< num_markedv - temp_<< " \n";
            // Choose from candidates
            if (num_markedv < N_){ // if candidates exists.
                // Choose the next startnid for Bfs/ next landmark
                std::uniform_int_distribution<int> distribution(0,candidates.size()-1);
                set<TInt>::const_iterator it(candidates.begin());
                std::advance(it,distribution(generator));
                StartNId = *it;
                //                landmarks.push_back(StartNId);
            }
            std::cout << "Landmark size "<<landmarks.size()<<" Candidates remaining "<<candidates.size()<<"\n";
        }
        std::cout <<" landmarks: ";
        for(auto i:landmarks) std::cout << i<<" ";
        std::cout<<"\n";
        
        selected_vertices.Reserve(landmarks.size());
        for (auto i:landmarks)  selected_vertices.Add(i);
    }
    
    //    if (criteria.compare("eps_densebfs")==0){
    //        std:: vector <std::pair <double,int>> vec_d_intpair;
    //        for (auto NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
    //            vec_d_intpair.push_back(std::make_pair(NI.GetOutDeg(),NId));
    //        }
    //        // weight_means = "cost". Nodes that have higher total cost=> most similar
    //        std::sort(vec_d_intpair.rbegin(), vec_d_intpair.rend()); // descending order, matlab version
    //
    //    }
    if (criteria.compare("eps_sptree_prune")==0){
        std::cout << "eps_sptree_prune\n";
        int N_ = Net->GetNodes();
        int E_ = Net->GetEdges();
        // Disjoint set
        //        typedef boost::associative_property_map<std::map<int,int>> rankType;
        
        typedef std::map<int,int> rank_t;
        typedef std::map<int,int> parent_t;
        rank_t rankmap;
        parent_t parentmap;
        boost::associative_property_map<rank_t> rank_pmap(rankmap);
        boost::associative_property_map<parent_t> parent_pmap(parentmap);
        
        typedef boost::associative_property_map<rank_t> Rank;
        typedef boost::associative_property_map<parent_t> Parent;
        boost::disjoint_sets<Rank,Parent> dset(rank_pmap,parent_pmap);
        
        //        // Marked map
        //        hash_map <int, bool> marked; // keep track of marked vertices
        //        marked.reserve(N_); // reserve memory for marked hash map
        //        // Distance map -> use `dist_to_cover' map
        
        //        int num_markedv = 0;
        int set_cnt = Net->GetNodes();
        hash_map< int,int> parent;
        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
            int n_id = NI.GetId();
            dist_to_cover[n_id] = 0.0;
            dset.make_set(n_id);
            parent[n_id] = n_id;
        }
        
        // Sort Edges. **************
        std::vector <std::tuple<value_t,int,int>> edgew_idmap;
        
        for (auto EI = Net->BegEI(); EI < Net->EndEI(); EI++){
            int u = EI.GetSrcNId();
            int v = EI.GetDstNId();
            edgew_idmap.push_back( std::make_tuple(Net->GetEDat(u,v),u,v) );
        }
        std::sort(edgew_idmap.rbegin(),edgew_idmap.rend());
        
        // Tree is represented as an weighted bidirectional network
        TPt<TNodeEDatNet<TInt, TFlt>> spt = TNodeEDatNet<TInt, TFlt>::New();
        
        for (auto edge_t = edgew_idmap.begin(); set_cnt>1 ;edge_t++){
            int u,v;
            double weight;
            
            weight = get<0>(*edge_t);
            u = get<1>(*edge_t);
            v = get<2>(*edge_t);
            if( u > v )     continue;
            
            std::cout << "Edge "<< "("<< u <<","<< v << " = "<<weight<<") ; set-count" <<set_cnt<<" \n";
            
            int rep_u = dset.find_set(u);
            int rep_v = dset.find_set(v);
            std::cout << "rep_u "<<rep_u<< " rep_v "<<rep_v<<"\n";
            if (rep_u!=rep_v){
                if (!spt->IsNode(u))
                    spt->AddNode(u);
                if (!spt->IsNode(v))
                    spt->AddNode(v);
                // Add edge in both direction for parent <-> child traversal
                
                if (parent[u] == u && parent[v] == v){
                    if (u<v) { parent[v] = u;   spt->AddEdge(u,v,Net->GetEDat(u,v)); std::cout << u<<" -> "<<v<<" " << u<< " Parent "<<parent[u]<< " " << v <<" parent "<<parent[v]<<"\n";}
                    else { parent[u] = v;   spt->AddEdge(v,u,Net->GetEDat(v,u));  std::cout << u<<" -> "<<v<<" " << u<< " Parent "<<parent[u]<< " " << v <<" parent "<<parent[v]<<"\n";}
                }
                else if (parent[u] != u && parent[v] == v) { parent[v] = u; spt->AddEdge(u,v,Net->GetEDat(u,v)); std::cout << u<<" -> "<<v<<" " << u<< " Parent "<<parent[u]<< " " << v <<" parent "<<parent[v]<<"\n";}
                else if (parent[u] == u && parent[v] != v) { parent[u] = v; spt->AddEdge(v,u,Net->GetEDat(v,u)); std::cout << u<<" -> "<<v<<" " << u<< " Parent "<<parent[u]<< " " << v <<" parent "<<parent[v]<<"\n";}
                else{  // (parent[u] != u && parent[v] != v){
                    //                    std::cout << u<< " rep " <<rep_u<< " "<<v <<" rep "<< rep_v<<"\n";
                    //                    if (rep_u < rep_v)  { parent[rep_v] = rep_u; spt->AddEdge(rep_u,rep_v,Net->GetEDat(rep_u,rep_v)); }
                    //                    else { parent[rep_u] = rep_v; spt->AddEdge(rep_v,rep_u,Net->GetEDat(rep_v,rep_u)); }
                    if (u<v)    { parent[v] = u;   spt->AddEdge(u,v,Net->GetEDat(u,v)); std::cout << u<<" -> "<<v<<" " << u<< " Parent "<<parent[u]<< " " << v <<" parent "<<parent[v]<<"\n";}
                    else
                    { parent[u] = v;   spt->AddEdge(v,u,Net->GetEDat(v,u));  std::cout << u<<" -> "<<v<<" " << u<< " Parent "<<parent[u]<< " " << v <<" parent "<<parent[v]<<"\n";}
                }
                dset.union_set(u,v);
                set_cnt--;
            }
            
        }
        
        //        for(auto EI= spt->BegEI(); EI<spt->EndEI();EI++){
        //            int u = EI.GetSrcNId();
        //            int v = EI.GetDstNId();
        //            std::cout << u<<" -> "<<v<<" " << u<< " Parent "<<parent[u]<< " " << v <<" parent "<<parent[v]<<"\n";
        //
        //        }
        // Mark the nodes in its eps-cover in spt. the distance between d_spt(root,*) > d_G(root,*) , so if d_spt < epsilon => d_G < epsilon as well. can not be larger.
        
        TIntFltH NIdDistH(spt->GetNodes());
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        
        // Traverse the spt and gather landmark. Start from root.
        int root;
        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
            if (parent[NI.GetId()]==NI.GetId()){    root = NI.GetId(); break;}
        }
        int num_markedv = 0;
        TSnapQueue<int> candidates(spt->GetNodes());
        std::set <int> nodes_inside_ring;
        
        candidates.Push(root);
        while(num_markedv<N_){
            std::set<int> candidates_gt2pes;
            std::set <int> Leafvec;
            boost::disjoint_sets<Rank,Parent> ringdset(rank_pmap,parent_pmap);
            // Initializations
            int StartNId = candidates.Top(); candidates.Pop();
            dist_to_cover[StartNId] = 0;
            NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
            marked[StartNId] = true;
            num_markedv++;
            selected_vertices.Add(StartNId);
            
            //            std::cout << "Partial BFS from node "<<StartNId<<"\n";
            TSnapQueue<int> Queue(spt->GetNodes());
            Queue.Push(StartNId); // BFS traversal queue
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                double Dist = NIdDistH.GetDat(NId);
                auto NodeI = spt->GetNI(NId);
                std::cout <<NId<<"\'s Children :\n";
                for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    
                    bool DstNId_marked = marked[DstNId];
                    
                    std::cout << DstNId<<"  "<< "Marked? "<<DstNId_marked << " weight = "<<spt->GetEDat(NId,DstNId)<<" \n";
                    if (!NIdDistH.IsKey(DstNId)){
                        NIdDistH.AddDat(DstNId, Dist + spt->GetEDat(NId,DstNId));
                        std::cout << "Dist "<<Dist + spt->GetEDat(NId,DstNId)<<"\n";
                        if (Dist + spt->GetEDat(NId,DstNId) <= 2*epsilon){
                            Queue.Push(DstNId);
                            if(Dist + spt->GetEDat(NId,DstNId) <= epsilon){
                                std::cout << "Marking "<< DstNId<<"\n";
                                num_markedv++;
                                marked[DstNId] = true;
                                if (dist_to_cover.find(DstNId)==dist_to_cover.end()) // does not exists
                                    dist_to_cover[DstNId] = Dist + Net->GetEDat(NId,DstNId); // Do not use spt's weight here. **** IMportant
                                else{
                                    //                            std::cout << "collision\n";
                                    dist_to_cover[DstNId] = std::min(Dist + Net->GetEDat(NId,DstNId),dist_to_cover[DstNId]);
                                }
                            }
                            else{
                                // we look for leaf node inside (eps,2eps) ring, since some nodes in its eps-ancestor may not be marked by later landmarks
                                if (spt->GetNI(DstNId).GetOutDeg() == 0){ // leaf
                                    std::cout << "Leaf "<<DstNId<<"\n";
                                    Leafvec.insert(DstNId);
                                    //                                    ringdset.make_set(DstNId);
                                    selected_vertices.Add(DstNId); // For now, just add those leaves as landmarks.
                                    marked[DstNId] = true;
                                    num_markedv++;
                                    // Traverse parent and mark ancestors
                                }
                                
                                //                                nodes_inside_ring.push_back(DstNId);
                            }
                            //                        std::cout << " Marked "<<DstNId << "\n";
                            
                        }
                        else{
                            std::cout << " Candidate add: "<<DstNId<<"\n";
                            candidates_gt2pes.insert(DstNId);
                        }
                    }
                }
            }
            std::cout << "- -- - -\n BFS QUEUE END \n ---------\n";
            
            
            //            bool visited[N_] = {false};
            ////            if (Leafvec.size()==1)  selected_vertices.Add(Leafvec[0]);
            //            hash_map<int,int> ancestor;
            //            for(int t_p = 0; tp < Leafvec.size(); t_p++){
            //                int u = Leafvec[t_p];
            //                if (dset.find_set(u) != u)  continue;
            //                int lca=-1;
            //                while (true){ // back-track till you find a marked parent meaning <= epsilon dist.
            //                    if(marked[u]){
            //                        ancestor[Leafvec[t_p]] = u;
            //                        break;
            //                    }
            //                    if (!visited[u]){
            //                        visited[u] = true;
            //                        u = p[u];
            //                    }
            //                    else{
            //                        lca = u;
            //                        dset.union_set(Leafvec[t_p],)
            //                        break;
            //                    }
            //                }
            //                if (lca==-1){
            //                    // lca not found
            //
            //                }
            //            }
            
            //                for(int s_p = t_p+1; sp< Leafvec.size(); s_p++){
            //                    int v = Leafvec[s_p];
            //                    //                                      selected_vertices.Add(DstNId);
            //
            //
            //                    while(true){
            //                        if (visited[v]) { lca = v; break; }
            //                        v = p[v];
            //                    }
            //                    if (lca != -1){
            //                        dset.union_set(Leafvec[t_p],Leafvec[s_p]);
            //                    }
            //                    else{
            //
            //                    }
            //                }
            //            }
            
            // Printing Candidates
            std::cout << " Printing candidates\n";
            int len_ = candidates.Len();
            for(;len_>0;len_--){
                auto val = candidates.Top(); candidates.Pop();
                std::cout << val << " ";
                candidates.Push(val);
            }
            std::cout <<"\n";
            
            // Printing candidates gt 2eps
            std::cout << " Printing candidates gt 2eps\n";
            for(auto it = candidates_gt2pes.begin(); it != candidates_gt2pes.end(); it++){
                std::cout << *it << " ";
            }
            std::cout <<"\n";
            
            // Prune candidates
            // Task 1: Run bfs on the spt from each of the candidates. but this time only in reverse direction i.e. look at indegrees/parents
            // Task 2: Remove vertices that are covered by the candidates from nodes_inside_ring.
            //            for(auto it = candidates_gt2pes.begin(); it != candidates_gt2pes.end(); it++){
            //                std::cout << " Pruning \n";
            //                int u = *it;
            //                TIntFltH u_NidDisth(100);
            //                TSnapQueue<int> Queue(spt->GetNodes());
            //
            //                if (!marked[u]){
            //                    std::cout << "unmarked "<<u<<"\n";
            //                    Queue.Push(u);
            //                    u_NidDisth.AddDat(u,0);
            ////                    marked[u] = true;
            ////                    num_markedv++;
            //                    candidates.Push(u); // subset of candidates> 2eps that are not marked by anyone.
            //                    while (! Queue.Empty() ){
            //                        const int NId = Queue.Top();  Queue.Pop();
            //                        marked[NId] = true;
            //                        auto Dist = u_NidDisth.GetDat(NId);
            //                        // Handle its parent
            //                        if (!u_NidDisth.IsKey(parent[NId]) && !marked[parent[NId]]){
            //                            u_NidDisth.AddDat(parent[NId], Dist + spt->GetEDat(parent[NId],NId));
            //                            Queue.Push(parent[NId]);
            //                        }
            ////                            marked[NId] = true;
            ////                            num_markedv++;
            ////                            candidates.Push(NId);
            ////                            if (!marked[parent[NId]]){
            ////                                u_NidDisth.AddDat(parent[NId], Dist + spt->GetEDat(parent[NId],NId));
            ////                                Queue.Push(parent[NId]);
            ////                            }
            //                            auto NodeI = spt->GetNI(parent[NId]);
            //                            for (int v = 0; v < NodeI.GetOutDeg(); v++) {
            //                                auto sibling = NodeI.GetOutNId(v);
            //                                if (sibling == NId) continue;
            //                                if (!u_NidDisth.IsKey(sibling) && !marked[sibling]){
            //                                    if ( Dist + spt->GetEDat(parent[NId],NId) + spt->GetEDat(parent[NId],sibling) <= epsilon){
            ////                                        marked[sibling] = true;
            //                                        u_NidDisth.AddDat(sibling, Dist + spt->GetEDat(parent[NId],NId) + spt->GetEDat(parent[NId],sibling) );
            //                                        Queue.Push(sibling);
            //                                    }
            //                                }
            //                            }
            ////                        }
            //                    }
            //                }
            //            }
            for(auto u: candidates_gt2pes){
                std::cout << " Mark from > 2eps points back ward. \n";
                //                TIntFltH u_NidDisth(100);
                //                TSnapQueue<int> Queue(spt->GetNodes());
                candidates.Push(u);
                while (!marked[parent[u]]){
                    std::cout << "unmarked "<<u<<"\n";
                    u = parent[u];
                    marked[u] = true;
                    num_markedv++;
                }
            }
            
            for (auto u: Leafvec){
                std::cout << " Mark from leaves";
                while (!marked[parent[u]]){
                    std::cout << "unmarked "<<u<<"\n";
                    u = parent[u];
                    marked[u] = true;
                    num_markedv++;
                }
            }
            
        }
    }
    if (criteria.compare("eps_kruskal")==0){
        int N_ = Net->GetNodes();
        int E_ = Net->GetEdges();
        // Disjoint set
        //        typedef boost::associative_property_map<std::map<int,int>> rankType;
        
        // Disjoint set
        //        std::map<int,int> rank;
        //        std::map<int,int> parent;
        typedef std::map<int,int> rank_t;
        typedef std::map<int,int> parent_t;
        rank_t rankmap;
        parent_t parentmap;
        boost::associative_property_map<rank_t> rank_pmap(rankmap);
        boost::associative_property_map<parent_t> parent_pmap(parentmap);
        
        typedef boost::associative_property_map<rank_t> Rank;
        typedef boost::associative_property_map<parent_t> Parent;
        boost::disjoint_sets<Rank,Parent> dset(rank_pmap,parent_pmap);
        
        // Marked map
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        // Distance map -> use `dist_to_cover' map
        
        int num_markedv = 0;
        
        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
            int n_id = NI.GetId();
            dist_to_cover[n_id] = 0.0;
            dset.make_set(n_id);
            
        }
        //        std::cout << " - - - - - - - - - - - - - - \n";
        //        std::cout << " Starting \n";
        //        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
        //            int n_id = NI.GetId();
        //            std::cout << n_id <<" => "<< dset.find_set(n_id)<<"\n";
        //        }
        //        std::cout << " - - - - - - - - - - - - - - \n";
        
        // Improvement on the size of landmarks:- If both u, v are marked check d(rep_u)+d(rep_v) + w(u,v) <= epsilon. If so merge u,v.
        // To do so, maintain d(rep_u) as the max of distances from root in that tree to the farthest node.
        //        hash_map <int, value_t> farthest_dist_tree;
        
        // Sort!!!
        //
        for (auto EI = Net->BegEI(); EI < Net->EndEI() && num_markedv < N_;EI++){
            int u = EI.GetSrcNId();
            int v = EI.GetDstNId();
            if( u > v )     continue;
            value_t weight = Net->GetEDat(u,v);
            std::cout << "Edge "<< "("<< u <<","<< v << " = "<<Net->GetEDat(u,v)<<")\n";
            std::cout << " Marked status = "<< marked[u]<< ","<<marked[v]<<"\n";
            //            if (num_markedv>= N_) break;
            int rep_u = dset.find_set(u);
            int rep_v = dset.find_set(v);
            
            if (rep_u!=rep_v){
                if (!marked[u] && marked[v]){
                    if ( dist_to_cover[v] + weight <= epsilon){
                        dset.link(rep_u,rep_v);
                        marked[u] = true;
                        dist_to_cover[u] = dist_to_cover[v] + weight;
                        num_markedv++;
                        //                        farthest_dist_tree[dset.find_set(u)] = max( farthest_dist_tree[dset.find_set(u)], dist_to_cover[u]);
                    }
                }
                else if(marked[u] && !marked[v]){
                    if ( dist_to_cover[u] + weight <= epsilon){
                        dset.link(rep_u,rep_v);
                        marked[v] = true;
                        dist_to_cover[v] = dist_to_cover[u] + weight;
                        num_markedv++;
                        //                        farthest_dist_tree[dset.find_set(u)] = max( farthest_dist_tree[dset.find_set(u)], dist_to_cover[v]);
                        //                                                std::cout << " - - - - - - - - - - - - - - \n";
                        std::cout << "Merge "<< u << " "<<v << " -> "<< dset.find_set(u)<<"\n";
                        std::cout << rankmap[u]<<" "<<rankmap[v]<<"\n";
                        //                                                for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
                        //                                                    int n_id = NI.GetId();
                        //                                                    std::cout << n_id <<" => "<< dset.find_set(n_id)<<"\n";
                        //                                                }
                        //                                                std::cout << " - - - - - - - - - - - - - - \n";
                    }
                }
                else if(!marked[u] && !marked[v]){
                    if (weight <= epsilon){
                        dset.union_set(u,v);
                        marked[u] = true;
                        marked[v] = true;
                        if (dset.find_set(u) == u){
                            dist_to_cover[v] = weight;
                            //                            farthest_dist_tree[u] = weight; // Improvement
                        }
                        else{
                            dist_to_cover[u] = weight;
                            //                            farthest_dist_tree[v] = weight;
                            
                        }
                        num_markedv+=2;
                        //                                                std::cout << " - - - - - - - - - - - - - - \n";
                        std::cout << "W* Merge "<< u << " "<<v << " -> "<< dset.find_set(u)<<"\n";
                        std::cout << rankmap[u]<<" "<<rankmap[v]<<"\n";
                        //                                                for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
                        //                                                    int n_id = NI.GetId();
                        //                                                    std::cout << n_id <<" => "<< dset.find_set(n_id)<<"\n";
                        //                                                }
                        //                                                std::cout << " - - - - - - - - - - - - - - \n";
                    }
                }
                else{
                    //                    if ( weight + farthest_dist_tree[rep_v] > epsilon || weight + farthest_dist_tree[rep_u] > epsilon) // it is hopeless to merge
                    //                        continue;
                    //                    else{
                    //                        std::cout << " Both marked => \n";
                    //                        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
                    //                            int n_id = NI.GetId();
                    //                            std::cout << n_id <<" => "<< dset.find_set(n_id)<<" dist to cover= "<<dist_to_cover[n_id]<<"\n";
                    //                            std::cout << "Now. "<< weight + farthest_dist_tree[rep_v]<<"\n";
                    //                            std::cout << "And. "<< weight + farthest_dist_tree[rep_u]<<"\n";
                    //                        }
                    //
                    //                        // Both u, v marked
                    //                        // Do nothing
                    //                        if( rep_u == u && rep_v == v && (weight + farthest_dist_tree[rep_v] <= epsilon) ){ // u is a representative. then farthest_dist_tree[u] contains the ball radius.
                    //                                        // Check whether w[u,v] + farthest_dist_tree[rep_v] <= epsilon. Then merge
                    //                            std::cout << "Special case1 "<<rankmap[u]<<" "<<rankmap[v]<<"\n";
                    //                        }
                    //                        if( rep_u == u && rep_v == v && weight + farthest_dist_tree[rep_u] <= epsilon){ // u is a representative. then farthest_dist_tree[u] contains the ball radius.
                    //                            // Check whether w[u,v] + farthest_dist_tree[rep_v] <= epsilon. Then merge
                    //                            std::cout << "Special case2 "<<rankmap[u]<<" "<<rankmap[v]<<"\n";
                    //                        }
                    //                    }
                }
            }
            
        }
        
        for(auto NI = Net->BegNI(); NI < Net->EndNI(); NI++){
            int n_id = NI.GetId();
            if (dset.find_set(n_id) == n_id)
                selected_vertices.Add(n_id);
        }
    }
    if (criteria.compare("eps_filtering")==0){
        // 1. Construct BFS tree.
        // 2. Select eps-Sample from tree traversal
        // 3. eps-sparsity Pruning
        
        // 1.
        int N_ = Net->GetNodes();
        std::cout << "Nodes: "<<N_<<"\n";
        PNGraph sptree = TNGraph::New( N_, N_ - 1 ); // Spanning tree as an undirected graph
        const TInt StartNId = Net->BegNI().GetId();
        
        bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
        TSnapQueue<int> Queue(InitBigQ?Net->GetNodes():1024);
        TIntFltH NIdDistH(InitBigQ?Net->GetNodes():1024);
        
        // Initializations
        IAssert(Net->IsNode(StartNId));
        NIdDistH.AddDat(StartNId, 0.0); // keep track of shortest path distances
        Queue.Push(StartNId); // BFS traversal queue
        // Run epsilon-BFS
        while (! Queue.Empty()){
            const int NId = Queue.Top();  Queue.Pop();
            if (!sptree->IsNode(NId))  sptree->AddNode(NId);
            auto Dist = NIdDistH.GetDat(NId);
            auto NodeI = Net->GetNI(NId);
            for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                const int DstNId = NodeI.GetOutNId(v);
                if (!NIdDistH.IsKey(DstNId)){
                    Queue.Push(DstNId);
                    NIdDistH.AddDat(DstNId, Dist + Net->GetEDat(NId,DstNId) );
                    if (!sptree->IsNode(DstNId))   sptree->AddNode(DstNId);
                    sptree->AddEdge(NId,DstNId);
                    
                }
            }
        }
        
        // 2.
        std::cout << " Tree Nodes: "<<sptree->GetNodes()<<" Edges: "<< sptree->GetEdges()<<"\n";
        // eps-BFS in the bfs tree.
        // traverse the tree level order. Mark nodes. And take any unmarked successor node down the tree.
        std::vector <int> eps_sample;
        
        TSnapQueue<int> candidates(InitBigQ?Net->GetNodes():1024);
        candidates.Push(Net->BegNI().GetId());
        while(!candidates.Empty()){
            auto StartNId = candidates.Top(); candidates.Pop();
            eps_sample.push_back(StartNId);
            std::cout << " BFS:- "<<StartNId<<"\n";
            //            marked[StartNId] = false; // mark the candidates as false. But its eps-cover as true.
            TSnapQueue<int> Queue2(InitBigQ?sptree->GetNodes():1024);
            TIntFltH NIdDistH2(InitBigQ?sptree->GetNodes():1024);
            NIdDistH2.AddDat(StartNId, 0.0); // keep track of shortest path distances
            Queue2.Push(StartNId); // BFS traversal queue
            while ( !Queue2.Empty() ){
                const int NId = Queue2.Top();  Queue2.Pop();
                auto Dist = NIdDistH2.GetDat(NId);
                auto NodeI = sptree->GetNI(NId);
                std::cout << " Visited: \n";
                for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    if (!NIdDistH2.IsKey(DstNId)){
                        if ( Dist + Net->GetEDat(NId,DstNId) > epsilon){
                            candidates.Push(DstNId);
                            std::cout << " Pushing in candidates: "<<DstNId<<"\n";
                        }
                        else{
                            Queue2.Push(DstNId);
                            NIdDistH2.AddDat(DstNId, Dist + Net->GetEDat(NId,DstNId) );
                            std::cout << DstNId<<"\n";
                        }
                    }
                }
            }
            std::cout << " Next\n";
        }
        
        std::cout << " Size of epsilon-sample: "<< eps_sample.size()<<"\n";
        // 3.
        
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        for(int i = 0; i<eps_sample.size(); i++){
            auto StartNId = eps_sample[i];
            TSnapQueue<int> Queue2(InitBigQ?sptree->GetNodes():1024);
            TIntFltH NIdDistH2(InitBigQ?sptree->GetNodes():1024);
            if(!marked[StartNId]){
                selected_vertices.Add(StartNId);
                Queue2.Push(StartNId);
                NIdDistH2.AddDat(StartNId, 0.0);
                dist_to_cover[StartNId] = 0.0;
                marked[StartNId] = true;
                while(!Queue2.Empty()){
                    // eps-traverse the graph.
                    const int NId = Queue2.Top();  Queue2.Pop();
                    auto Dist = NIdDistH2.GetDat(NId);
                    auto NodeI = sptree->GetNI(NId);
                    for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                        const int DstNId = NodeI.GetOutNId(v);
                        if (!NIdDistH2.IsKey(DstNId) && Dist + Net->GetEDat(NId,DstNId) <= epsilon ){
                            Queue2.Push(DstNId);
                            NIdDistH2.AddDat(DstNId, Dist + Net->GetEDat(NId,DstNId) );
                            dist_to_cover[DstNId] = Dist + Net->GetEDat(NId,DstNId);
                            marked[DstNId] = true;
                        }
                    }
                }
            }
        }
    }
    
    if (criteria.compare("epsnetring")==0){
        std::cout << "epsring\n";
        vector <TInt> landmarks;
        std::default_random_engine generator; // random generator
        
        // Select first landmark
        TIntV NodeIdV;
        int N_ = Net->GetNodes();
        Net->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        //        landmarks.push_back(NodeIdV[0]);
        int num_markedv = 0;
        int StartNId = NodeIdV[0];
        //        int StartNId = 0;
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        std::set <TInt> candidates_gt2pes; // Always maintain the (eps,2eps) unmarked ring
        candidates_gt2pes.insert(StartNId); // Candidates always maintain a set of potential landmarks to select for each iteration.
        // Except inside Partial BFS candidates contain all unmarked potential landmarks.
        std::set <TInt> candidates_gteps;
        //        std::cout << " First landmark = "<<StartNId<<"\n";
        bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
        
        while (num_markedv < N_){
            
            //            std::cout <<" Before marking: "<<num_markedv<<" ";
            int temp_ = num_markedv;
            if (!marked[StartNId]){
                marked[StartNId] = true;
                num_markedv++;
                landmarks.push_back(StartNId);
            }
            else{
                //                std::cout << "Marked "<<num_markedv<<"\n";
                //                std::cout << "invariant violeted. check why!! \n";
            }
            //            std::cout << "marked # "<<num_markedv<<" total # "<<N_<<"\n";
            // run partial bfs and mark
            // BFS house-keeping data structures
            
            
            run_partialBFS_wgraph(Net,InitBigQ,StartNId,epsilon,candidates_gt2pes,candidates_gteps,marked,num_markedv,dist_to_cover);
            
            // TO DO: Need to verify that previously existing key in dist_to_cover is not overwritten here.
            //            std::cout <<" After : "<<num_markedv<<" Diff: "<< num_markedv - temp_<< " \n";
            // Choose from candidates
            if (num_markedv < N_){ // if candidates exists.
                if(candidates_gt2pes.size()>0){
                    // Choose the next startnid for Bfs/ next landmark
                    std::uniform_int_distribution<int> distribution(0,candidates_gt2pes.size()-1);
                    set<TInt>::const_iterator it(candidates_gt2pes.begin());
                    std::advance(it,distribution(generator));
                    StartNId = *it;
                }
                else{ // there are some vertices > eps <2eps . So the bfs couldn't find any candidates. Leaving holes in that range
                    //                    std::cout << " Holes in the selection\n";
                    std::uniform_int_distribution<int> distribution(0,candidates_gteps.size()-1);
                    set<TInt>::const_iterator it(candidates_gteps.begin());
                    std::advance(it,distribution(generator));
                    StartNId = *it;
                }
            }
            
            
            //            else{ // increment by delta and run bfs again
            ////                std::cout << "\n delta inc \n";
            //                for(int delta = 1; candidates.size()==0 && num_markedv!=N_ ;delta++){
            //
            //                    float ep_ = pow(2.0,delta);
            ////                    std::cout << "delta "<<delta<<" ep_ "<<ep_<<"\n";
            //                    run_partialBFS(g,InitBigQ,StartNId,ep_,candidates,marked,num_markedv);
            ////                    std::cout << "inside delta "<< num_markedv<<"\n";
            //                    if (candidates.size()>0){
            //                        std::uniform_int_distribution<int> distribution(0,candidates.size()-1);
            //                        StartNId = candidates[distribution(generator)];
            //                        landmarks.push_back(StartNId);
            //                    }
            //                }
            
            //            }
            //                        std::cout << "Landmark size "<<landmarks.size()<<" \n";
            //            std::cout <<" landmarks: ";
            //            for(auto i:landmarks) std::cout << i<<" ";
            //            std::cout<<"\n";
        }
        //        std::cout << "Total marked = "<< num_markedv<<"\n";
        //        run_partialBFS(g,false,StartNId,epsilon,candidates,marked,num_markedv,dist_to_cover);
        
        selected_vertices.Reserve(landmarks.size());
        for (auto i:landmarks)  selected_vertices.Add(i);
        std::cout << "Total marked = "<< num_markedv<<"\n";
    }
    
    if (criteria.compare("eps_baseline2")==0){
        // Store eps-bfs tree for each node and store for each vertex -> size of the tree (i.e. # of vertices eps-covered)
        int N_ = Net->GetNodes();
        std::cout << "Eps- baseline\n";
        // Step 1. For each node run epsilon-BFS from it, and store the eps-cover size
        //        TIntH eps_bfs_sz(N_); // key: Nodeid (TINt) value: number of vertices eps-covered.
        std:: vector <std::pair <int,int>> eps_bfs_sz;
        //        std::vector <int> eps_bfs_sz;
        hash_map < int, TIntFltH > covers;
        //        hash_map < int, std::map<int,int> > covers;
        bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
        for (auto NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
            const TInt StartNId = NI.GetId();
            
            int v;
            TSnapQueue<int> Queue(InitBigQ?Net->GetNodes():1024);
            TIntFltH NIdDistH(InitBigQ?Net->GetNodes():1024);
            
            // Initializations
            IAssert(Net->IsNode(StartNId));
            NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
            Queue.Push(StartNId); // BFS traversal queue
            // Run epsilon-BFS
            //            std::cout << " Start: -> \n";
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                auto Dist = NIdDistH.GetDat(NId);
                auto NodeI = Net->GetNI(NId);
                //                std::cout << "NID : "<<NId<<" "<< Dist<<"\n";
                for (v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    if (!NIdDistH.IsKey(DstNId)){
                        auto weight = Net->GetEDat(NId,DstNId);
                        if (Dist + weight <= epsilon){
                            //                            std::cout << " DstNId: "<<DstNId<<"\n";
                            Queue.Push(DstNId);
                            NIdDistH.AddDat(DstNId, Dist+ weight);
                        }
                    }
                }
            }
            eps_bfs_sz.push_back(std::make_pair(NIdDistH.Len(),StartNId));
            //            eps_bfs_sz[StartNId] = NIdDistH.Len();
            covers[StartNId] = NIdDistH;
            //            eps_bfs_sz[StartNId] = NIdDistH.Len();
            //
        }
        
        std::cout << "Sort by dat\n";
        // Step 2. Sort the hashtable eps_bfs_sz in descending order of value
        //        std::sort(eps_bfs_sz.rbegin(),eps_bfs_sz.rend());
        //        std::cout << "PRint: \n";
        //        for (int i = 0; i< N_; i++){
        //            std::cout << eps_bfs_sz[i].second<< " cover size: "<< eps_bfs_sz[i].first<<"\n";
        //        }
        
        vector <bool> marked(N_,false);
        int total_marked = 0;
        // Difference from eps_baseline:-
        //                          0. After taking ith landmark and its eps-cover -
        //              1. Update the cover-size of other unmarked vertices.
        //              2. Sort Again.
        //              3. Select the vertex that has the largest cover size.
        for (int i = 0; total_marked<N_; i++){
            std::sort (eps_bfs_sz.begin(), eps_bfs_sz.end(),[ ]( const std::pair<int,int> &lhs, const std::pair<int,int> &rhs )
                       {
                           return lhs.first > rhs.first;
                       });
            
            
            
            int node_id = eps_bfs_sz[0].second;
            int cov_sz = eps_bfs_sz[0].first;
            //            std::cout << " node_id= "<<node_id<<" Sz: "<<cov_sz<<"\n";
            // 0.
            //                std::cout << "PRint before: \n";
            //                for (int x = 0; x< eps_bfs_sz.size(); x++){
            //                    std::cout << eps_bfs_sz[x].second<< " cover size: "<< eps_bfs_sz[x].first<<"\n";
            //                }
            int num_markedv = 0;
            // mark its covered vertices
            for (auto it = covers[node_id].BegI(); it<covers[node_id].EndI(); it++){
                const TInt cov_node = it.GetKey();
                if (!marked[cov_node]){
                    num_markedv++;
                    total_marked++;
                    marked[cov_node] = true;
                }
                // 1. Update
                for(int j = 1; j<eps_bfs_sz.size(); j++){
                    int DstNId = eps_bfs_sz[j].second;
                    //                        for (auto it2 = covers.begin(); it2 != covers.end(); it2++){
                    //                            if (it2->first == cov_node) continue;
                    //                            auto toupdate = covers[j];
                    if (covers[DstNId].IsKey(cov_node)){
                        covers[DstNId].DelKey(cov_node);
                        eps_bfs_sz[j].first--;
                    }
                    
                }
                auto distance = it.GetDat();
                if (dist_to_cover.find(cov_node) == dist_to_cover.end())
                    dist_to_cover[cov_node] = distance;
                else
                    dist_to_cover[cov_node] = std::fmin(TFlt(distance),dist_to_cover[cov_node]);
            }
            covers.erase(node_id);
            eps_bfs_sz.erase(eps_bfs_sz.begin());
            
            selected_vertices.Add(node_id);
            //                std::cout << "Added: "<<node_id<< " Marked# "<<num_markedv<<"\n";
        }
        
    }
    return selected_vertices;
}

void run_partialBFS(const PUNGraph &g, bool InitBigQ, int StartNId, float epsilon, std::set <TInt> &candidates, std::set <TInt> &candidates_gteps_lt2eps, hash_map <int, bool> &marked , int &num_markedv, hash_map <int, int> &dist_to_cover){
    // Append a set of new candidates for landmark selection in the next iteration.
    int v;
    TSnapQueue<int> Queue(InitBigQ?g->GetNodes():1024);
    TIntH NIdDistH(InitBigQ?g->GetNodes():1024);
    
    // Initializations
    IAssert(g->IsNode(StartNId));
    dist_to_cover[StartNId] = 0;
    NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
    Queue.Push(StartNId); // BFS traversal queue
    candidates.erase(StartNId);
    //    std::cout << "Partial BFS from node "<<StartNId<<"\n";
    while (! Queue.Empty()){
        const int NId = Queue.Top();  Queue.Pop();
        const int Dist = NIdDistH.GetDat(NId);
        auto NodeI = g->GetNI(NId);
        //        std::cout <<NId<<" ";
        for (v = 0; v < NodeI.GetOutDeg(); v++) {
            const int DstNId = NodeI.GetOutNId(v);
            bool DstNId_marked = marked[DstNId];
            
            if (!NIdDistH.IsKey(DstNId)){
                NIdDistH.AddDat(DstNId, Dist+1);
                if (Dist + 1 <= 2*epsilon){
                    Queue.Push(DstNId);
                    if (Dist + 1 <= epsilon){
                        candidates.erase(DstNId);
                        if(!DstNId_marked){
                            marked[DstNId] = true;
                            num_markedv++;
                            candidates_gteps_lt2eps.erase(DstNId);
                        }
                        //                        //Compare! If (v'(key) -> d_l1) was shortest before, now if v' -> d_l2 is shorter, need to update.
                        if (dist_to_cover.find(DstNId)==dist_to_cover.end()) // does not exists
                            dist_to_cover[DstNId] = Dist + 1;
                        else{
                            //                            std::cout << "collision\n";
                            dist_to_cover[DstNId] = min(Dist + 1,dist_to_cover[DstNId]);
                        }
                        
                        //                        std::cout << " Marked "<<DstNId << "\n";
                    }
                    else{
                        if(!DstNId_marked && candidates_gteps_lt2eps.find(DstNId) == candidates_gteps_lt2eps.end()) // unmarked and not added already
                            candidates_gteps_lt2eps.insert(DstNId);
                    }
                }
                else{ // Dist + 1> 2*epsilon. Add DstNId as candidate
                    if(!DstNId_marked)
                        candidates.insert(DstNId);
                }
            }
        }
    }
}

TIntV select_landmarks(const PUNGraph &g, const std::string criteria, hash_map <int, int> &dist_to_cover, int num_landmarks = 1, float epsilon = 0.01){
    TIntV selected_vertices;
    
    
    if (criteria.compare("degree")==0){
        selected_vertices.Reserve(num_landmarks);
        std:: cout << "degree selection" <<std::endl;
        std:: vector <std::pair <double,int>> vec_d_intpair;
        for (TUNGraph::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
            const int NId = NI.GetId();
            const double DegCentr = g->GetNI(NId).GetDeg();
            vec_d_intpair.push_back(std::make_pair(DegCentr,NId));
        }
        std::sort(vec_d_intpair.rbegin(), vec_d_intpair.rend());
        int N_ = g->GetNodes();
        int max_K = TMath::Mn(num_landmarks, N_);
        for (int tries = 0; tries < max_K; tries++){
            selected_vertices.Add(vec_d_intpair[tries].second);
            // std:: cout << vec_d_intpair[tries].first << " "<< vec_d_intpair[tries].second<<std::endl;
        }
    }
    if (criteria.compare("random")==0){
        selected_vertices.Reserve(num_landmarks);
        std:: cout << "random selection " << TMath::Mn(num_landmarks, g->GetNodes()) <<std::endl;
        TIntV NodeIdV;
        
        g->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        int N_ =  g->GetNodes();
        int max_K = TMath::Mn(num_landmarks,N_);
        for (int tries = 0; tries < max_K; tries++) {
            selected_vertices.Add(NodeIdV[tries]);
        }
    }
    if (criteria.compare("eigen")==0){
        selected_vertices.Reserve(num_landmarks);
        std:: cout << "eigenvector selection" <<std::endl;
        TIntFltH EigH;
        TSnap::GetEigenVectorCentr(g, EigH);
        std:: vector <std::pair <double,int>> vec_d_intpair;
        for (TUNGraph::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
            const TInt NId = NI.GetId();
            double EigCentr = EigH.GetDat(NId);
            vec_d_intpair.push_back(std::make_pair(EigCentr,NId));
        }
        sort(vec_d_intpair.begin(), vec_d_intpair.end());
        int N_ = g->GetNodes();
        int max_K = TMath::Mn(num_landmarks, N_);
        
        for (int tries = 0; tries < max_K; tries++){
            selected_vertices.Add(vec_d_intpair[(N_-1)-tries].second);
            //            std:: cout << vec_d_intpair[(N_-1)-tries].first << " "<< vec_d_intpair[(N_-1)-tries].second<<std::endl;
        }
    }
    if (criteria.compare("maxmin")==0){
        // The code for point-cloud needs to be modified
        selected_vertices.Reserve(num_landmarks);
        int N_ =  g->GetNodes();
        int max_K = TMath::Mn(num_landmarks,N_);
        int StartNId;
        int NextStartNId;
        TIntV NodeIdV;
        g->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        hash_map <int,int> min_dist_landmarkset;
        min_dist_landmarkset.reserve(N_);
        
        for(auto NI = g->BegNI(); NI < g->EndNI(); NI++)
            min_dist_landmarkset[NI.GetId()] = std::numeric_limits<int>::max();
        
        for (int tries = 0; tries < max_K; tries++) {
            if (tries == 0){
                StartNId = NodeIdV[tries];
            }
            TBreathFS<PUNGraph> BFS(g); // run bfs on whole graph instead of BFS(subG)
            BFS.DoBfs(StartNId, true, false, -1, TInt::Mx);
            int cur_max = TInt::Mn;
            for(auto NI = g->BegNI(); NI < g->EndNI(); NI++){ // Update min dist to landmark set
                int DstNId = NI.GetId();
                int dist = BFS.GetHops(StartNId,DstNId);
                min_dist_landmarkset[DstNId] = min ( min_dist_landmarkset[DstNId], dist );
                if (cur_max < min_dist_landmarkset[DstNId]){
                    NextStartNId = DstNId;
                    cur_max = min_dist_landmarkset[DstNId];
                }
            }
            selected_vertices.Add(StartNId);
            StartNId = NextStartNId;
        }
        // Actually min_dist_landmarkset is your dist_to_cover.
        dist_to_cover = std::move(min_dist_landmarkset);
    }
    if (criteria.compare("epsmaxmin")==0){
        // The code for point-cloud needs to be modified
        selected_vertices.Reserve(num_landmarks);
        int N_ =  g->GetNodes();
        int max_K = TMath::Mn(num_landmarks,N_);
        int StartNId;
        int NextStartNId;
        TIntV NodeIdV;
        g->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        hash_map <int,int> min_dist_landmarkset;
        min_dist_landmarkset.reserve(N_);
        hash_map <int,bool> marked;
        marked.reserve(N_);
        
        for(auto NI = g->BegNI(); NI < g->EndNI(); NI++)
            min_dist_landmarkset[NI.GetId()] = std::numeric_limits<int>::max();
        
        int num_markedv = 0;
        for (int tries = 0; num_markedv<N_; tries++) {
            if (tries == 0){
                StartNId = NodeIdV[tries];
                // Run eps-BFS and mark its cover
                
            }
            //            TBreathFS<PUNGraph> BFS(g);
            //            BFS.DoBfs(StartNId, true, false, -1, TInt::Mx);
            
            bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
            int v;
            TSnapQueue<int> Queue(InitBigQ?g->GetNodes():1024);
            TIntH NIdDistH(InitBigQ?g->GetNodes():1024);
            
            // run bfs on whole graph
            // Initializations
            IAssert(g->IsNode(StartNId));
            NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
            Queue.Push(StartNId); // BFS traversal queue
            marked[StartNId] = true;
            num_markedv++;
            // Run epsilon-BFS
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                const int Dist = NIdDistH.GetDat(NId);
                auto NodeI = g->GetNI(NId);
                for (v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    if (!NIdDistH.IsKey(DstNId)){
                        if (Dist + 1 <= epsilon && !marked[DstNId]){
                            num_markedv++;
                            marked[DstNId] = true;
                        }
                        Queue.Push(DstNId);
                        NIdDistH.AddDat(DstNId, Dist+1);
                    }
                }
            }
            
            int cur_max = TInt::Mn;
            for(auto NI = g->BegNI(); NI < g->EndNI(); NI++){ // Update min dist to landmark set
                int DstNId = NI.GetId();
                int dist = NIdDistH.GetDat(DstNId);
                min_dist_landmarkset[DstNId] = min ( min_dist_landmarkset[DstNId], dist );
                if (cur_max < min_dist_landmarkset[DstNId]){
                    NextStartNId = DstNId;
                    cur_max = min_dist_landmarkset[DstNId];
                }
            }
            selected_vertices.Add(StartNId);
            StartNId = NextStartNId;
        }
        // Actually min_dist_landmarkset is your dist_to_cover.
        dist_to_cover = std::move(min_dist_landmarkset);
    }
    
    if (criteria.compare("eps_filtering")==0){
        // 1. Construct BFS tree.
        // 2. Select eps-Sample from tree traversal
        // 3. eps-sparsity Pruning
        
        // 1.
        int N_ = g->GetNodes();
        std::cout << "Nodes: "<<N_<<"\n";
        PNGraph sptree = TNGraph::New( N_, N_ - 1 ); // Spanning tree as an undirected graph
        
        // Select first landmark
        std::default_random_engine generator; // random generator
        TIntV NodeIdV;
        g->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        const TInt StartNId = NodeIdV[0]; // First vertex is chosen u.a.r
        std::cout <<" StartNode: "<<StartNId<<"\n";
        
        bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
        TSnapQueue<int> Queue(InitBigQ?g->GetNodes():1024);
        TIntH NIdDistH(InitBigQ?g->GetNodes():1024);
        
        // Initializations
        IAssert(g->IsNode(StartNId));
        NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
        Queue.Push(StartNId); // BFS traversal queue
        // Run epsilon-BFS
        while (! Queue.Empty()){
            const int NId = Queue.Top();  Queue.Pop();
            if (!sptree->IsNode(NId))  sptree->AddNode(NId);
            const int Dist = NIdDistH.GetDat(NId);
            auto NodeI = g->GetNI(NId);
            for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                const int DstNId = NodeI.GetOutNId(v);
                if (!NIdDistH.IsKey(DstNId)){
                    Queue.Push(DstNId);
                    NIdDistH.AddDat(DstNId, Dist + 1 );
                    if (!sptree->IsNode(DstNId))   sptree->AddNode(DstNId);
                    sptree->AddEdge(NId,DstNId);
                    
                }
            }
        }
        
        // 2.
        std::cout << " Tree Nodes: "<<sptree->GetNodes()<<" Edges: "<< sptree->GetEdges()<<"\n";
        // eps-BFS in the bfs tree.
        // traverse the tree level order. Mark nodes. And take any unmarked successor node down the tree.
        std::vector <int> eps_sample;
        
        TSnapQueue<int> candidates(InitBigQ?g->GetNodes():1024);
        candidates.Push(StartNId);
        while(!candidates.Empty()){
            auto StartNId = candidates.Top(); candidates.Pop();
            eps_sample.push_back(StartNId);
            //            std::cout << " BFS:- "<<StartNId<<"\n";
            //            marked[StartNId] = false; // mark the candidates as false. But its eps-cover as true.
            TSnapQueue<int> Queue2(InitBigQ?sptree->GetNodes():1024);
            TIntH NIdDistH2(InitBigQ?sptree->GetNodes():1024);
            NIdDistH2.AddDat(StartNId, 0); // keep track of shortest path distances
            Queue2.Push(StartNId); // BFS traversal queue
            while ( !Queue2.Empty() ){
                const int NId = Queue2.Top();  Queue2.Pop();
                const int Dist = NIdDistH2.GetDat(NId);
                auto NodeI = sptree->GetNI(NId);
                //                std::cout << " Visited: \n";
                for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    if (!NIdDistH2.IsKey(DstNId)){
                        if ( Dist + 1 > epsilon){
                            candidates.Push(DstNId);
                            //                            std::cout << " Pushing in candidates: "<<DstNId<<"\n";
                        }
                        else{
                            Queue2.Push(DstNId);
                            NIdDistH2.AddDat(DstNId, Dist + 1 );
                            //                            std::cout << DstNId<<"\n";
                        }
                    }
                }
            }
            //            std::cout << " Next\n";
        }
        
        std::cout << " Size of epsilon-sample: "<< eps_sample.size()<<"\n";
        // 3.
        
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        for(int i = 0; i<eps_sample.size(); i++){
            auto StartNId = eps_sample[i];
            TSnapQueue<int> Queue2(InitBigQ?sptree->GetNodes():1024);
            TIntH NIdDistH2(InitBigQ?sptree->GetNodes():1024);
            if(!marked[StartNId]){
                selected_vertices.Add(StartNId);
                Queue2.Push(StartNId);
                NIdDistH2.AddDat(StartNId, 0);
                dist_to_cover[StartNId] = 0.0;
                marked[StartNId] = true;
                while(!Queue2.Empty()){
                    // eps-traverse the graph.
                    const int NId = Queue2.Top();  Queue2.Pop();
                    const int Dist = NIdDistH2.GetDat(NId);
                    auto NodeI = sptree->GetNI(NId);
                    for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                        const int DstNId = NodeI.GetOutNId(v);
                        if (!NIdDistH2.IsKey(DstNId) && Dist + 1 <= epsilon ){
                            Queue2.Push(DstNId);
                            NIdDistH2.AddDat(DstNId, Dist + 1 );
                            dist_to_cover[DstNId] = Dist + 1;
                            marked[DstNId] = true;
                        }
                    }
                }
            }
        }
    }
    
    if (criteria.compare("eps_densedfs")==0){
        std::cout << "Eps densebfs\n";
        int N_ = g->GetNodes();
        std:: vector <std::pair <double,int>> vec_d_intpair;
        for (TUNGraph::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
            const int NId = NI.GetId();
            const double DegCentr = g->GetNI(NId).GetDeg();
            vec_d_intpair.push_back(std::make_pair(DegCentr,NId));
        }
        std::sort(vec_d_intpair.rbegin(), vec_d_intpair.rend());
        
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        
        int total_marked = 0;
        // Traverse and Run DFS
        for (auto it : vec_d_intpair){
            if(total_marked>=N_) break;
            int cur_id = it.second;
            int Deg = it.first;
            
            if (!marked[cur_id]){
                //                std::cout << "Degree: "<<Deg<<"\n";
                int num_markedv = 0;
                dist_to_cover[cur_id] = 0;
                marked[cur_id] = true;
                total_marked++;
                std::stack <int> st;
                st.push(cur_id);
                
                while(!st.empty()){
                    int SrcNId = st.top(); st.pop();
                    auto NodeI = g->GetNI(SrcNId);
                    const int Dist = dist_to_cover[SrcNId];
                    int nextSrcNId=-1;
                    int maxdeg = -1;
                    for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                        int DstNId = NodeI.GetOutNId(v);
                        if(!marked[DstNId] && Dist+1 <= epsilon){
                            marked[DstNId] = true;
                            num_markedv++;
                            total_marked++;
                            dist_to_cover[DstNId] = Dist + 1;
                            st.push(DstNId);
                            //                            if (g->GetNI(DstNId).GetOutDeg() > maxdeg){
                            //                                nextSrcNId = DstNId;
                            //                                maxdeg = g->GetNI(DstNId).GetOutDeg();
                            //                            }
                        }
                    }
                    //                    if (nextSrcNId!=-1){
                    //                        st.push(nextSrcNId);
                    //                    }
                }
                selected_vertices.Add(cur_id);
                //                std::cout << "Marked# "<<num_markedv<<"\n";
            }
            
        }
    }
    if (criteria.compare("eps_kruskal")==0){
        typedef boost::associative_property_map<std::map<int,int>> rankType;
        
        std::cout << " Eps kruskal\n";
        int N_ = g->GetNodes();
        int E_ = g->GetEdges();
        
        // Disjoint set
        std::map<int,int> rank;
        std::map<int,int> parent;
        boost::disjoint_sets<rankType,rankType> dset(boost::make_assoc_property_map(rank),boost::make_assoc_property_map(parent));
        
        // Marked map
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        // Distance map -> use `dist_to_cover' map
        int num_markedv = 0;
        
        for(auto NI = g->BegNI(); NI < g->EndNI(); NI++){
            int n_id = NI.GetId();
            //            std::cout << " id = "<<n_id<<"\n";
            dist_to_cover[n_id] = 0;
            dset.make_set(n_id);
        }
        // All edge weights are 1. No sorting needs to be done.
        for (auto EI = g->BegEI(); EI < g->EndEI() ;EI++){
            int u = EI.GetSrcNId();
            int v = EI.GetDstNId();
            if (num_markedv>= N_) break;
            int rep_u = dset.find_set(u);
            int rep_v = dset.find_set(v);
            if (rep_u!=rep_v){
                if (!marked[u] && marked[v]){
                    if ( dist_to_cover[v] + 1 <= epsilon){
                        dset.link(rep_u,rep_v);
                        marked[u] = true;
                        dist_to_cover[u] = dist_to_cover[v] + 1;
                        num_markedv++;
                        //                        std::cout << " - - - - - - - - - - - - - - \n";
                        //                        for(auto NI = g->BegNI(); NI < g->EndNI(); NI++){
                        //                            int n_id = NI.GetId();
                        //                            std::cout << n_id <<" => "<< dset.find_set(n_id)<<"\n";
                        //                        }
                        //                        std::cout << " - - - - - - - - - - - - - - \n";
                    }
                }
                else if(marked[u] && !marked[v]){
                    if ( dist_to_cover[u] + 1 <= epsilon){
                        dset.link(rep_u,rep_v);
                        marked[v] = true;
                        dist_to_cover[v] = dist_to_cover[u] + 1;
                        num_markedv++;
                        //                        std::cout << " - - - - - - - - - - - - - - \n";
                        //                        for(auto NI = g->BegNI(); NI < g->EndNI(); NI++){
                        //                            int n_id = NI.GetId();
                        //                            std::cout << n_id <<" => "<< dset.find_set(n_id)<<"\n";
                        //                        }
                        //                        std::cout << " - - - - - - - - - - - - - - \n";
                    }
                }
                else if(!marked[u] && !marked[v]){
                    if (1 <= epsilon){
                        dset.union_set(u,v);
                        marked[u] = true;
                        marked[v] = true;
                        if (dset.find_set(u) == u)
                            dist_to_cover[v] = 1;
                        else
                            dist_to_cover[u] = 1;
                        num_markedv+=2;
                        //                        std::cout << " - - - - - - - - - - - - - - \n";
                        //                        for(auto NI = g->BegNI(); NI < g->EndNI(); NI++){
                        //                            int n_id = NI.GetId();
                        //                            std::cout << n_id <<" => "<< dset.find_set(n_id)<<"\n";
                        //                        }
                        //                        std::cout << " - - - - - - - - - - - - - - \n";
                    }
                }
                else{
                    // Do nothing
                }
            }
            
        }
        
        for(auto NI = g->BegNI(); NI < g->EndNI(); NI++){
            int n_id = NI.GetId();
            if (dset.find_set(n_id) == n_id)
                selected_vertices.Add(n_id);
        }
    }
    if (criteria.compare("eps_baseline")==0){
        // Store eps-bfs tree for each node and store for each vertex -> size of the tree (i.e. # of vertices eps-covered)
        int N_ = g->GetNodes();
        std::cout << "Eps- baseline\n";
        // Step 1. For each node run epsilon-BFS from it, and store the eps-cover size
        //        TIntH eps_bfs_sz(N_); // key: Nodeid (TINt) value: number of vertices eps-covered.
        std:: vector <std::pair <int,int>> eps_bfs_sz;
        hash_map < int, TIntH > covers;
        
        for (TUNGraph::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
            const TInt StartNId = NI.GetId();
            
            bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
            int v;
            TSnapQueue<int> Queue(InitBigQ?g->GetNodes():1024);
            TIntH NIdDistH(InitBigQ?g->GetNodes():1024);
            
            // Initializations
            IAssert(g->IsNode(StartNId));
            NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
            Queue.Push(StartNId); // BFS traversal queue
            // Run epsilon-BFS
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                const int Dist = NIdDistH.GetDat(NId);
                auto NodeI = g->GetNI(NId);
                for (v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    if (!NIdDistH.IsKey(DstNId)){
                        if (Dist + 1 <= epsilon){
                            Queue.Push(DstNId);
                            NIdDistH.AddDat(DstNId, Dist+1);
                        }
                    }
                }
            }
            eps_bfs_sz.push_back(std::make_pair(NIdDistH.Len(),StartNId));
            covers[StartNId] = NIdDistH;
            //            eps_bfs_sz[StartNId] = NIdDistH.Len();
            //
        }
        
        std::cout << "Sort by dat\n";
        // Step 2. Sort the hashtable eps_bfs_sz in descending order of value
        std::sort(eps_bfs_sz.rbegin(),eps_bfs_sz.rend());
        
        
        vector <bool> marked(N_,false);
        for (int i = 0; i< N_; i++){
            int node_id = eps_bfs_sz[i].second;
            int cov_sz = eps_bfs_sz[i].first;
            if(!marked[node_id]){
                int num_markedv = 0;
                // mark its covered vertices
                for (auto it = covers[node_id].BegI();it<covers[node_id].EndI(); it++){
                    const TInt cov_node = it.GetKey();
                    if (!marked[cov_node])
                        num_markedv++;
                    marked[cov_node] = true;
                    const int distance = it.GetDat();
                    if (dist_to_cover.find(cov_node) == dist_to_cover.end())
                        dist_to_cover[cov_node] = distance;
                    else
                        dist_to_cover[cov_node] = min(distance,dist_to_cover[cov_node]);
                }
                selected_vertices.Add(node_id);
                //                std::cout << "Node: "<< node_id<< " Marked# "<<num_markedv<<"\n";
                //                        std::cout << "PRint: \n";
                //                        for (int i = 0; i< N_; i++){
                //                            std::cout << eps_bfs_sz[i].second<< " cover size: "<< eps_bfs_sz[i].first<<"\n";
                //                        }
            }
        }
        
    }
    
    if (criteria.compare("eps_baseline2")==0){
        // Store eps-bfs tree for each node and store for each vertex -> size of the tree (i.e. # of vertices eps-covered)
        int N_ = g->GetNodes();
        std::cout << "Eps- baseline\n";
        // Step 1. For each node run epsilon-BFS from it, and store the eps-cover size
        //        TIntH eps_bfs_sz(N_); // key: Nodeid (TINt) value: number of vertices eps-covered.
        std:: vector <std::pair <int,int>> eps_bfs_sz;
        //        std::vector <int> eps_bfs_sz;
        hash_map < int, TIntH > covers;
        //        hash_map < int, std::map<int,int> > covers;
        bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
        for (TUNGraph::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
            const TInt StartNId = NI.GetId();
            
            
            int v;
            TSnapQueue<int> Queue(InitBigQ?g->GetNodes():1024);
            TIntH NIdDistH(InitBigQ?g->GetNodes():1024);
            
            // Initializations
            IAssert(g->IsNode(StartNId));
            NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
            Queue.Push(StartNId); // BFS traversal queue
            // Run epsilon-BFS
            //            std::cout << " Start: -> \n";
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                const int Dist = NIdDistH.GetDat(NId);
                auto NodeI = g->GetNI(NId);
                //                std::cout << "NID : "<<NId<<" "<< Dist<<"\n";
                for (v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    if (!NIdDistH.IsKey(DstNId)){
                        if (Dist + 1 <= epsilon){
                            //                            std::cout << " DstNId: "<<DstNId<<"\n";
                            Queue.Push(DstNId);
                            NIdDistH.AddDat(DstNId, Dist+1);
                        }
                    }
                }
            }
            eps_bfs_sz.push_back(std::make_pair(NIdDistH.Len(),StartNId));
            //            eps_bfs_sz[StartNId] = NIdDistH.Len();
            covers[StartNId] = NIdDistH;
            //            eps_bfs_sz[StartNId] = NIdDistH.Len();
            //
        }
        
        std::cout << "Sort by dat\n";
        // Step 2. Sort the hashtable eps_bfs_sz in descending order of value
        //        std::sort(eps_bfs_sz.rbegin(),eps_bfs_sz.rend());
        //        std::cout << "PRint: \n";
        //        for (int i = 0; i< N_; i++){
        //            std::cout << eps_bfs_sz[i].second<< " cover size: "<< eps_bfs_sz[i].first<<"\n";
        //        }
        
        vector <bool> marked(N_,false);
        int total_marked = 0;
        // Difference from eps_baseline:-
        //                          0. After taking ith landmark and its eps-cover -
        //              1. Update the cover-size of other unmarked vertices.
        //              2. Sort Again.
        //              3. Select the vertex that has the largest cover size.
        for (int i = 0; total_marked<N_; i++){
            std::sort (eps_bfs_sz.begin(), eps_bfs_sz.end(),[ ]( const std::pair<int,int> &lhs, const std::pair<int,int> &rhs )
                       {
                           return lhs.first > rhs.first;
                       });
            
            
            
            int node_id = eps_bfs_sz[0].second;
            int cov_sz = eps_bfs_sz[0].first;
            //            std::cout << " node_id= "<<node_id<<" Sz: "<<cov_sz<<"\n";
            // 0.
            //                std::cout << "PRint before: \n";
            //                for (int x = 0; x< eps_bfs_sz.size(); x++){
            //                    std::cout << eps_bfs_sz[x].second<< " cover size: "<< eps_bfs_sz[x].first<<"\n";
            //                }
            int num_markedv = 0;
            // mark its covered vertices
            for (auto it = covers[node_id].BegI(); it<covers[node_id].EndI(); it++){
                const TInt cov_node = it.GetKey();
                if (!marked[cov_node]){
                    num_markedv++;
                    total_marked++;
                    marked[cov_node] = true;
                }
                // 1. Update
                for(int j = 1; j<eps_bfs_sz.size(); j++){
                    int DstNId = eps_bfs_sz[j].second;
                    //                        for (auto it2 = covers.begin(); it2 != covers.end(); it2++){
                    //                            if (it2->first == cov_node) continue;
                    //                            auto toupdate = covers[j];
                    if (covers[DstNId].IsKey(cov_node)){
                        covers[DstNId].DelKey(cov_node);
                        eps_bfs_sz[j].first--;
                    }
                    
                }
                const int distance = it.GetDat();
                if (dist_to_cover.find(cov_node) == dist_to_cover.end())
                    dist_to_cover[cov_node] = distance;
                else
                    dist_to_cover[cov_node] = min(distance,dist_to_cover[cov_node]);
            }
            covers.erase(node_id);
            eps_bfs_sz.erase(eps_bfs_sz.begin());
            
            selected_vertices.Add(node_id);
            //                std::cout << "Added: "<<node_id<< " Marked# "<<num_markedv<<"\n";
        }
        
    }
    if (criteria.compare("epsnetring")==0){
        std::cout << "epsring\n";
        vector <TInt> landmarks;
        std::default_random_engine generator; // random generator
        
        // Select first landmark
        TIntV NodeIdV;
        int N_ = g->GetNodes();
        g->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        //        landmarks.push_back(NodeIdV[0]);
        int num_markedv = 0;
        int StartNId = NodeIdV[0];
        //        int StartNId = 0;
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        std::set <TInt> candidates_gt2pes; // Always maintain the (eps,2eps) unmarked ring
        candidates_gt2pes.insert(StartNId); // Candidates always maintain a set of potential landmarks to select for each iteration.
        // Except inside Partial BFS candidates contain all unmarked potential landmarks.
        std::set <TInt> candidates_gteps;
        //        std::cout << " First landmark = "<<StartNId<<"\n";
        while (num_markedv < N_){
            
            //            std::cout <<" Before marking: "<<num_markedv<<" ";
            int temp_ = num_markedv;
            if (!marked[StartNId]){
                marked[StartNId] = true;
                num_markedv++;
                landmarks.push_back(StartNId);
            }
            else{
                //                std::cout << "Marked "<<num_markedv<<"\n";
                //                std::cout << "invariant violeted. check why!! \n";
            }
            //            std::cout << "marked # "<<num_markedv<<" total # "<<N_<<"\n";
            // run partial bfs and mark
            // BFS house-keeping data structures
            bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
            
            run_partialBFS(g,InitBigQ,StartNId,epsilon,candidates_gt2pes,candidates_gteps,marked,num_markedv,dist_to_cover);
            
            // TO DO: Need to verify that previously existing key in dist_to_cover is not overwritten here.
            //            std::cout <<" After : "<<num_markedv<<" Diff: "<< num_markedv - temp_<< " \n";
            // Choose from candidates
            if (num_markedv < N_){ // if candidates exists.
                if(candidates_gt2pes.size()>0){
                    // Choose the next startnid for Bfs/ next landmark
                    std::uniform_int_distribution<int> distribution(0,candidates_gt2pes.size()-1);
                    set<TInt>::const_iterator it(candidates_gt2pes.begin());
                    std::advance(it,distribution(generator));
                    StartNId = *it;
                }
                else{ // there are some vertices > eps <2eps . So the bfs couldn't find any candidates. Leaving holes in that range
                    //                    std::cout << " Holes in the selection\n";
                    std::uniform_int_distribution<int> distribution(0,candidates_gteps.size()-1);
                    set<TInt>::const_iterator it(candidates_gteps.begin());
                    std::advance(it,distribution(generator));
                    StartNId = *it;
                }
            }
            
            
            //            else{ // increment by delta and run bfs again
            ////                std::cout << "\n delta inc \n";
            //                for(int delta = 1; candidates.size()==0 && num_markedv!=N_ ;delta++){
            //
            //                    float ep_ = pow(2.0,delta);
            ////                    std::cout << "delta "<<delta<<" ep_ "<<ep_<<"\n";
            //                    run_partialBFS(g,InitBigQ,StartNId,ep_,candidates,marked,num_markedv);
            ////                    std::cout << "inside delta "<< num_markedv<<"\n";
            //                    if (candidates.size()>0){
            //                        std::uniform_int_distribution<int> distribution(0,candidates.size()-1);
            //                        StartNId = candidates[distribution(generator)];
            //                        landmarks.push_back(StartNId);
            //                    }
            //                }
            
            //            }
            //            std::cout << "Landmark size "<<landmarks.size()<<" \n";
            //            std::cout <<" landmarks: ";
            //            for(auto i:landmarks) std::cout << i<<" ";
            //            std::cout<<"\n";
        }
        //        std::cout << "Total marked = "<< num_markedv<<"\n";
        //        run_partialBFS(g,false,StartNId,epsilon,candidates,marked,num_markedv,dist_to_cover);
        
        selected_vertices.Reserve(landmarks.size());
        for (auto i:landmarks)  selected_vertices.Add(i);
        std::cout << "Total marked = "<< num_markedv<<"\n";
    }
    
    if (criteria.compare("epsrand")==0){
        std::cout << "epsrand\n";
        vector <TInt> landmarks;
        std::default_random_engine generator; // random generator
        
        // Select first landmark
        TIntV NodeIdV;
        int N_ = g->GetNodes();
        g->GetNIdV(NodeIdV);
        TRnd r;
        r.Randomize(); // If i do not call Randomize(), same random sequence is generated.
        NodeIdV.Shuffle(r);
        
        int num_markedv = 0;
        int StartNId = NodeIdV[0];
        hash_map <int, bool> marked; // keep track of marked vertices
        marked.reserve(N_); // reserve memory for marked hash map
        TIntV vertices(N_);
        g->GetNIdV(vertices);
        std::set <TInt> candidates;
        for (int i=0; i<vertices.Len();i++)   candidates.insert(vertices[i]);
        while (num_markedv< N_){
            //            std::cout <<" Before marking: "<<num_markedv<<" ";
            //            int temp_ = num_markedv;
            if (!marked[StartNId]){
                marked[StartNId] = true;
                num_markedv++;
            }
            //            std::cout << "marked # "<<num_markedv<<" total # "<<N_<<"\n";
            // run partial bfs and mark
            // BFS house-keeping data structures
            bool InitBigQ = N_>10*1024? true: false; // Whether to initialize large sized queue or hashtable.
            
            
            int v;
            TSnapQueue<int> Queue(InitBigQ?g->GetNodes():1024);
            TIntH NIdDistH(InitBigQ?g->GetNodes():1024);
            //            hash_map <int,int> NIdDistH;
            
            // Initializations
            IAssert(g->IsNode(StartNId));
            //            NIdDistH[StartNId] = 0;
            candidates.erase(StartNId);
            landmarks.push_back(StartNId);
            NIdDistH.AddDat(StartNId, 0); // keep track of shortest path distances
            dist_to_cover.insert({StartNId,0});
            Queue.Push(StartNId); // BFS traversal queue
            //    marked[StartNId] = true; // marked flag
            //    num_markedv++;
            
            while (! Queue.Empty()){
                const int NId = Queue.Top();  Queue.Pop();
                const int Dist = NIdDistH.GetDat(NId);
                //                const int Dist = NIdDistH[NId];
                auto NodeI = g->GetNI(NId);
                //        std::cout <<NId<<" ";
                for (v = 0; v < NodeI.GetOutDeg(); v++) {
                    const int DstNId = NodeI.GetOutNId(v);
                    //            if (! marked[DstNId]) {
                    if (!NIdDistH.IsKey(DstNId)){
                        //                    if (NIdDistH.find(DstNId)==NIdDistH.end()){ // DstNId key does not exists
                        NIdDistH.AddDat(DstNId, Dist+1);
                        //                          NIdDistH[DstNId] = Dist + 1;
                        if (Dist + 1 <= epsilon){
                            Queue.Push(DstNId);
                            if (dist_to_cover.find(DstNId)==dist_to_cover.end()) // does not exists
                                dist_to_cover[DstNId] = Dist + 1;
                            else
                                dist_to_cover[DstNId] = min(Dist + 1,dist_to_cover[DstNId]);
                            if(!marked[DstNId]){
                                marked[DstNId] = true;
                                num_markedv++;
                                int n_removed = candidates.erase(DstNId);
                                assert(n_removed == 1);
                            }
                            //                    std::cout << " Marked "<<DstNId << "\n";
                        }
                    }
                }
                
            }
            //            std::cout <<" After : "<<num_markedv<<" Diff: "<< num_markedv - temp_<< " \n";
            // Choose from candidates
            if (num_markedv < N_){ // if candidates exists.
                // Choose the next startnid for Bfs/ next landmark
                std::uniform_int_distribution<int> distribution(0,candidates.size()-1);
                set<TInt>::const_iterator it(candidates.begin());
                std::advance(it,distribution(generator));
                StartNId = *it;
                //                landmarks.push_back(StartNId);
            }
            //            std::cout << "Landmark size "<<landmarks.size()<<" Candidates remaining "<<candidates.size()<<"\n";
        }
        //        std::cout <<" landmarks: ";
        //        for(auto i:landmarks) std::cout << i<<" ";
        //        std::cout<<"\n";
        
        selected_vertices.Reserve(landmarks.size());
        for (auto i:landmarks)  selected_vertices.Add(i);
    }
    
    
    return selected_vertices;
}

//void run_starclustering(const PUNGraph &g,const TIntV &centers, std::vector <TIntV> &clusters,TBreathFS<PUNGraph>&BFS){
//    // I am going to populate clusters variable. To avoid returning it
//    //    g.DoBfs_concurrent(centers,clusters)
//
//    BFS.DoBfs_concurrent(centers,clusters);
//    //    for(int i=0;i < BFS.cl_boundary.size();i++){
//    //        std::cout <<"cluster bf "<<i<<"\n";
//    //        for(int it=0; it<BFS.cl_boundary[i].Len();it++) std::cout<<BFS.cl_boundary[i][it]<<" ";
//    //        std::cout<<"\n";
//    //    }
//}

//void construct_flatdist_matrix(const PUNGraph& graph, TIntV& subg_vertices,std::vector <float>& flat_matrix){
//    // First construct subgraph
//    PUNGraph SubG = TSnap::GetSubGraph(graph, subg_vertices);
////    TBreathFS<PUNGraph> BFS(SubG, true, false);
////     BFS.DoBfs(SrcNId, true, ! IsDir, DstNId, TInt::Mx);
//    int subg_length = subg_vertices.Len();
//    for (int i = 0; i < subg_length; i++){
//        const int SrcNId = subg_vertices[i];
//        for (int j = i; j< subg_length; j++) {
//            const int DstNId = subg_vertices[j];
//            int dist = TSnap::GetShortPath(SubG,SrcNId,DstNId);
//            flat_matrix.push_back(dist);
//        }
//    }
//}

// improvement 1: do bfs for only subg_length times.
//template <class nodewType,class edgewType>
//void weighted_construct_flatdist_matrix(const TPt<TNodeEDatNet<nodewType, edgewType>>& Net, const TIntV& subg_vertices,std::vector <value_t>& flat_matrix){
//    // First construct subgraph
//    //TNodeEDatNet<nodewType, edgewType>explicit_subgraph(subg_vertices.Len(), 0);
//    TPt<TNodeEDatNet<nodewType, edgewType>> pSubG = TNodeEDatNet<nodewType, edgewType>::New();
//    TIntSet NIdSet(subg_vertices.Len());
//    for (int i = 0; i < subg_vertices.Len(); i++) {
//        if(Net->IsNode(subg_vertices[i])){
//            pSubG->AddNode(subg_vertices[i]);
//            NIdSet.AddKey(subg_vertices[i]);
//        }
//    }
//    //  std::cout << "Adding edges\n\n";
//    for (int n = 0; n < NIdSet.Len(); n++) {
//        const int SrcNId = NIdSet[n];
//        // std::cout << "xxxx ---- "<<SrcNId << " ";
//        auto NI = Net->GetNI(SrcNId);
//        for (int edge = 0; edge < NI.GetOutDeg(); edge++) {
//            const int OutNId = NI.GetOutNId(edge);
//            if (NIdSet.IsKey(OutNId)) {
//                pSubG->AddEdge(SrcNId, OutNId);
//                //pSubG->AddEdge(OutNId,SrcNId);
//            }
//        }
//    }
//
//    for (auto EI = pSubG->BegEI(); EI < pSubG->EndEI(); EI++) {
//        pSubG->SetEDat(EI.GetSrcNId(),EI.GetDstNId(),Net->GetEDat(EI.GetSrcNId(),EI.GetDstNId()));
//    }
//
//    // print the subgraph
//    //    for (TUNGraph::TNodeI NI = SubG->BegNI(); NI < SubG->EndNI(); NI++) {
//    //        printf("node id %d with out-degree %d and in-degree %d\n",
//    //               NI.GetId(), NI.GetOutDeg(), NI.GetInDeg());
//    //    }
//    // std::cout << "printing subgraph "<<std::endl;
//    //  for (auto EI = pSubG->BegEI(); EI < pSubG->EndEI(); EI++) {
//    //      std:: cout << EI.GetSrcNId() << " "<< EI.GetDstNId() << " = "<<pSubG->GetEDat(EI.GetSrcNId(),EI.GetDstNId()) << std::endl;
//    //  }
//    // Run Bfs and construct flat matrix
//    int subg_length = subg_vertices.Len();
//    for (int i = 0; i < subg_length; i++){
//        const int SrcNId = subg_vertices[i];
//        std::cout << SrcNId<<" ";
//        //TBreathFSW<PIntNEDNet> BFS(Net,true,false);
//        //BFS.DoBfs_concurrent(selected_vertices,clusters);
//        TBreathFSW<TPt<TNodeEDatNet<nodewType, edgewType>>> BFS(Net,true,false);
//        //TBreathFS<PUNGraph> BFS(SubG);
//        BFS.DoBfs(SrcNId, true, false, -1, TInt::Mx);
//        // always compute lower distance matrix as flat matrix
//        for (int j = 0; j< i; j++) {
//            const int DstNId = subg_vertices[j];
//            value_t dist = BFS.GetHops(SrcNId,DstNId);
//            flat_matrix.push_back(dist);
//        }
//    }
//    std::cout <<"\n";
//    std::cout << "flat matrix length: "<<flat_matrix.size()<< "\n";
//}

//void construct_flatdist_matrix(const PUNGraph& graph, const TIntV& subg_vertices,std::vector <value_t>& flat_matrix){
//    // First construct subgraph
//    PUNGraph SubG = TSnap::GetSubGraph(graph, subg_vertices);
//    // print the subgraph
//    //    for (TUNGraph::TNodeI NI = SubG->BegNI(); NI < SubG->EndNI(); NI++) {
//    //        printf("node id %d with out-degree %d and in-degree %d\n",
//    //               NI.GetId(), NI.GetOutDeg(), NI.GetInDeg());
//    //    }
//    //    for (TUNGraph::TEdgeI EI = SubG->BegEI(); EI < SubG->EndEI(); EI++) {
//    //        printf("edge (%d, %d)\n", EI.GetSrcNId(), EI.GetDstNId());
//    //    }
//
//    int subg_length = subg_vertices.Len();
//    if(subg_length==1) return;
//    for (int i = 0; i < subg_length; i++){
//        const int SrcNId = subg_vertices[i];
//        TBreathFS<PUNGraph> BFS(graph); // run bfs on whole graph instead of BFS(subG)
//        BFS.DoBfs(SrcNId, true, false, -1, TInt::Mx);
//        for (int j = 0; j< i; j++) { // compute lower distance matrix as flat matrix
//            const int DstNId = subg_vertices[j];
//            value_t dist = BFS.GetHops(SrcNId,DstNId);
//            flat_matrix.push_back(dist);
//        }
//    }
//}

template <class nodewType,class edgewType>
void weighted_construct_landmark_distmat(const TPt<TNodeEDatNet<nodewType, edgewType>>& Net, std::vector <value_t>& flat_matrix_lower_tri, hash_map <int, double> &m, const TIntV &landmarks){
    /*
     Mimics unweighted case. Except the graph type and BFS traversal call.
     */
    int L = landmarks.Len();
    int N = Net->GetNodes();
    vector < TIntFltH > dists;
    int maxint = std::numeric_limits<int>::max();
    int minint = std::numeric_limits<int>::min();
    
    for (int i=0; i<L; i++){ // Pre-compute BFS distances from landmarks to vertices.
        TBreathFSW<TPt<TNodeEDatNet<nodewType, edgewType>>> BFS(Net,true,false);
        BFS.DoBfs(landmarks[i], true, false, -1, TInt::Mx);
        
        dists.push_back(std::move(BFS.NIdDistH));
    }
    std::cout << "dists 0-len: "<<dists[0].Len()<<"\n";
    
    // The rest is the same as unweighted graph
    //    // eps-net-ring Rips code
    //        for (int i=0; i<L; i++){
    //            for(int j=0; j<i; j++){
    //                TFlt Dist;
    //                value_t e_ij = (dists[i].IsKeyGetDat(j,Dist)?Dist.Val:std::numeric_limits<int>::max());
    //                flat_matrix_lower_tri.push_back( e_ij );
    //
    //            }
    //        }
    
    // Lazy witness code
    for (int i=0; i<L; i++){
        for(int j=0; j<i; j++){
            value_t e_ij = std::numeric_limits<value_t>::infinity(); // distance to witness
            for (auto NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
                //                //        std::cout << "bfs "<<xx<<" ";
                //                xx++;
                auto n = NI.GetId();
                if( landmarks[i] == n || landmarks[j] == n){ continue;}
                
                // Find the nearest witness and distance to that witness. Javaplex implementation O(k^2*N)
                TFlt Dist;
                value_t dist_i = (dists[i].IsKeyGetDat(n, Dist)?Dist.Val:std::numeric_limits<int>::max());
                value_t dist_j = (dists[j].IsKeyGetDat(n, Dist)?Dist.Val:std::numeric_limits<int>::max());
                value_t mx = std::max(dist_i,dist_j);
                //                std::cout << " i "<<i <<" j "<<j<<" mx: "<<mx<<"\n";
                if ( mx < m[n] )
                    mx = 0.0;
                else
                    mx = mx - m[n];
                //                if (mx < e_ij)
                //                    e_ij = mx;
                e_ij = std::min(mx,e_ij);
            }
            flat_matrix_lower_tri.push_back( e_ij );
            
        }
    }
    
}

void construct_landmark_distmat(PUNGraph &g, std::vector <value_t>& flat_matrix_lower_tri, hash_map <int, int> &m, const TIntV &landmarks){
    // L = number of landmarks
    int L = landmarks.Len();
    int N = g->GetNodes();
    
    //    for (int i=0; i<L; i++){
    //        for(int j=0; j<i; j++){
    //            value_t e_ij = std::numeric_limits<value_t>::infinity(); // distance to witness
    //            int witness; // the witness
    //            // Find the nearest witness and distance to that witness. Javaplex implementation O(k^2*N)
    //            for (auto NI = g->BegNI(); NI < g->EndNI(); NI++) {
    //                auto n = NI.GetId();
    //                if( landmarks[i] == n || landmarks[j] == n) continue;
    //                double d_max = -std::numeric_limits<value_t>::infinity();
    //
    //            }
    //        }
    //    }
    /** // Will take a long long time (And its wrong code)
     int maxint = std::numeric_limits<int>::max();
     int minint = std::numeric_limits<int>::min();
     for (int i=0; i<L; i++)
     for(int j=0; j<i; j++)
     flat_matrix_lower_tri.push_back(maxint);
     
     int xx = 1;
     for (auto NI = g->BegNI(); NI < g->EndNI(); NI++) {
     std::cout << "bfs "<<xx<<"\n";
     xx++;
     auto n = NI.GetId();
     TBreathFS<PUNGraph>BFS(g);
     BFS.DoBfs(n,true, false);
     
     //        double d_max = -std::numeric_limits<value_t>::infinity();
     int idx_mat = 0;
     for (int i=0; i<L; i++){
     for(int j=0; j<i; j++){
     if( landmarks[i] == n || landmarks[j] == n){ continue;}
     //                value_t e_ij = std::numeric_limits<value_t>::infinity(); // distance to witness
     //                int witness; // the witness
     // Find the nearest witness and distance to that witness. Javaplex implementation O(k^2*N)
     value_t dist_i = BFS.GetHops(n,i);
     value_t dist_j = BFS.GetHops(n,j);
     value_t mx = std::max(std::max(dist_i,dist_j) - m[n],0.0f);
     if ( mx < flat_matrix_lower_tri[idx_mat]){
     flat_matrix_lower_tri[idx_mat] = mx;
     }
     idx_mat++;
     }
     }
     }
     **/
    
    // Alternate way
    vector < TIntH > dists;
    int maxint = std::numeric_limits<int>::max();
    int minint = std::numeric_limits<int>::min();
    //    for (int i=0; i<L; i++)
    //        for(int j=0; j<=i; j++){
    //            if (i==j) flat_matrix_lower_tri.push_back(0.0);
    //            else
    //                flat_matrix_lower_tri.push_back(maxint);
    //        }
    //    std::cout << "flat mat "<<flat_matrix_lower_tri.size()<<"\n";
    for (int i=0; i<L; i++){ // Pre-compute BFS distances from landmarks to vertices.
        TBreathFS<PUNGraph>BFS(g);
        BFS.DoBfs(landmarks[i],true, false);
        dists.push_back(std::move(BFS.NIdDistH));
    }
    std::cout << "dists 0-len: "<<dists[0].Len()<<"\n";
    //
    //    // eps-net-ring Rips code
    //        for (int i=0; i<L; i++){
    //            for(int j=0; j<i; j++){
    //                TInt Dist;
    //                value_t e_ij = (dists[i].IsKeyGetDat(j,Dist)?Dist.Val:std::numeric_limits<int>::max());
    //                flat_matrix_lower_tri.push_back( e_ij );
    //
    //            }
    //        }
    
    // Lazy witness code
    for (int i=0; i<L; i++){
        for(int j=0; j<i; j++){
            value_t e_ij = std::numeric_limits<value_t>::infinity(); // distance to witness
            for (auto NI = g->BegNI(); NI < g->EndNI(); NI++) {
                //                //        std::cout << "bfs "<<xx<<" ";
                //                xx++;
                auto n = NI.GetId();
                if( landmarks[i] == n || landmarks[j] == n){ continue;}
                
                //                int witness; // the witness
                // Find the nearest witness and distance to that witness. Javaplex implementation O(k^2*N)
                TInt Dist;
                value_t dist_i = (dists[i].IsKeyGetDat(n, Dist)?Dist.Val:std::numeric_limits<int>::max());
                value_t dist_j = (dists[j].IsKeyGetDat(n, Dist)?Dist.Val:std::numeric_limits<int>::max());
                value_t mx = std::max(dist_i,dist_j);
                //                std::cout << " i "<<i <<" j "<<j<<" mx: "<<mx<<"\n";
                if ( mx < m[n] )
                    mx = 0.0;
                else
                    mx = mx - m[n];
                if (mx < e_ij)
                    e_ij = mx;
            }
            flat_matrix_lower_tri.push_back( e_ij );
            
        }
    }
    
    //    for (int i=0; i<L; i++){
    //        for(int j=0; j<i; j++){
    //            std::cout << flat_matrix_lower_tri[i*L+j]<<" ";
    //        }
    //        std::cout <<"\n";
    //    }
    //    std::cout <<"\n";
    std::cout << "construct_landmark_distmat Ends\n";
    dists.clear();
}


//void print_merged_barcode(std::vector<std::vector<std::vector<std::pair<value_t,value_t>>>>& all_intervals, int max_dim, std::string outfilename){
//    for(int cluster_id=0;cluster_id<all_intervals.size();cluster_id++){
//        std::cout <<"cluster id: "<<cluster_id<<"\n";
//        if(all_intervals[cluster_id].size()==0) continue;
//        for (int d = 0; d<= max_dim; d++){
//            std::cout << "dimension "<<d<<": "<< all_intervals[cluster_id][d].size()<<"\n";
//            std::ofstream myfile;
//            myfile.open(outfilename+std::to_string(cluster_id)+"_"+std::to_string(d)+".csv");
//            {
//                for(int k=0;k<all_intervals[cluster_id][d].size();k++)
//                    if(all_intervals[cluster_id][d][k].second != std::numeric_limits<value_t>::max())
//                        myfile<< all_intervals[cluster_id][d][k].first << ", "<< all_intervals[cluster_id][d][k].second<<"\n";
//                    else
//                        myfile<< all_intervals[cluster_id][d][k].first << ", Inf"<<"\n";
//            }
//            myfile.close();
//        }
//    }
//}

void print_barcode(std::vector<std::vector<std::pair<value_t,value_t>>>& all_intervals, int max_dim, std::string outfilename){
    std::cout << " Printing barcodes \n"<<"max_dim = "<<max_dim<<"\n";
    for (int d = 0; d<= max_dim; d++){
        std::cout << "dimension "<<d<<": "<< all_intervals[d].size()<<"\n";
        std::ofstream myfile;
        myfile.open(outfilename+"_"+std::to_string(d)+".csv");
        {
            for(int k=0;k<all_intervals[d].size();k++)
                if(all_intervals[d][k].second != std::numeric_limits<value_t>::max()){
                    myfile<< all_intervals[d][k].first << ", "<< all_intervals[d][k].second<<"\n";
                    std::cout << all_intervals[d][k].first << ", "<< all_intervals[d][k].second<<"\n";
                }
                else{
                    myfile<< all_intervals[d][k].first << ", Inf"<<"\n";
                    std::cout << all_intervals[d][k].first << ", Inf"<<"\n";
                }
        }
        myfile.close();
    }
    
}

int main(int argc, char** argv) {
    
    char* filename = nullptr;
    int num_landmarks = -1;
    float epsilon = 0;
    file_format format = GRAPH;
    index_t dim_max = 1;
    value_t threshold = std::numeric_limits<value_t>::max();
    bool connected_;
    float ratio = 1;
    std::string heuristic;
    int iterations = 1; // optional parameter used for running each heuristic several times
#ifdef USE_COEFFICIENTS
    coefficient_t modulus = 2;
#else
    const coefficient_t modulus = 2;
#endif
    
    for (index_t i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help") {
            print_usage_and_exit(0);
        }
        else if (arg == "--iterations"){
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            iterations = std::stol(parameter, &next_pos);
            
        }
        else if (arg == "--heuristic"){
            heuristic = std::string(argv[++i]);
        }
        else if (arg == "--dim") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            dim_max = std::stol(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        } else if (arg == "--num-landmarks") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            num_landmarks = std::stol(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        }else if (arg == "--epsilon") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            epsilon = std::stof(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        }else if (arg == "--threshold") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            threshold = std::stof(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        } else if (arg == "--ratio") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            ratio = std::stof(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        } else if (arg == "--format") {
            std::string parameter = std::string(argv[++i]);
            if (parameter == "lower-distance")
                format = LOWER_DISTANCE_MATRIX;
            else if (parameter == "upper-distance")
                format = UPPER_DISTANCE_MATRIX;
            else if (parameter == "distance")
                format = DISTANCE_MATRIX;
            else if (parameter == "point-cloud"){
                format = POINT_CLOUD;
                std::cout << "point cloud"<<std::endl;
            }
            else if (parameter == "wgraph"){
                format = WGRAPH;
                std::cout <<"weighted graph"<<std::endl;
            }
            else if (parameter == "graph"){
                format = GRAPH;
                std::cout <<"graph"<<std::endl;
            }
            else if (parameter == "dipha")
                format = DIPHA;
            else if (parameter == "sparse")
                format = SPARSE;
            else if (parameter == "ripser")
                format = RIPSER;
            else
                print_usage_and_exit(-1);
#ifdef USE_COEFFICIENTS
        } else if (arg == "--modulus") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            modulus = std::stol(parameter, &next_pos);
            if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
        } else {
            if (filename) { print_usage_and_exit(-1); }
            filename = argv[i];
        }
    }
    
    std::ifstream file_stream(filename);
    if (filename && file_stream.fail()) {
        std::cerr << "couldn't open file " << filename << std::endl;
        exit(-1);
    }
    
    if (format == SPARSE) {
        string cluster_id = "";
        sparse_distance_matrix dist =
        read_sparse_distance_matrix(filename ? file_stream : std::cin);
        std::cout << "sparse distance matrix with " << dist.size() << " points and "
        << dist.num_edges << "/" << (dist.size() * (dist.size() - 1)) / 2 << " entries"
        << std::endl;
        
        ripser<sparse_distance_matrix>(std::move(dist), dim_max, threshold, ratio, modulus)
        .compute_barcodes(filename);
    }
    else if(format == GRAPH){
        //        std::string heuristic = "epsnetring";
        //        std::string heuristic = "epsrand";
        //        std::string heuristic = "epsmaxmin";
        //        std::string heuristic = "random";
        //        std::string heuristic = "eps_baseline";
        //        std::string heuristic = "eps_baseline2";
        //        std::string heuristic = "eps_kruskal";
        //        std::string heuristic = "eps_densedfs";
        //          std::string heuristic = "eps_sptree_prune";
        //          std::string heuristic = "eps_filtering";
        
        
        std::cout <<"fname (unweighted)"<< filename << std::endl;
        std::set<char> delims{'/'};
        std::set<char> delims2{'.'};
        std::vector<std::string> path = splitpath(filename, delims);
        std::vector<std::string> path2 = splitpath(path.back(), delims2);
        
        //for(auto it:path) std::cout<<it<<"-";
        //  std::cout<<"\n";
        //std::cout<<path2[path2.size()-2]<<"\n";
        
        // Need to handle weighted graphs.
        PUNGraph inputg = TSnap::LoadEdgeList<PUNGraph>(filename,0,1);
        TIntPrV SccSzCnt;
        //        TSnap::GetSccSzCnt(inputg, SccSzCnt);
        connected_ = TSnap::IsConnected(inputg);
        std::cout << " * * * * * \n";
        std::cout << "Graph with "<<inputg->GetNodes()<< " Nodes and "<< inputg->GetEdges() << " Edges\n" << std::endl;
        std::cout << (connected_?"Connected":"Not connected") << " \n* * * * * * \n";
        //        // Doing Graph Statistics.
        //  TSnap::PrintInfo(inputg, filename, "", inputg->GetNodes()>Kilo(1));
        ////        std::cout << TSnap::GetBfsFullDiam(inputg,100)<<std::endl;
        //        std::cout << inputg->GetNodes()<< " "<< inputg->GetEdges() << std::endl;
        //TVec < TPair<TInt, TInt> > CntV;
        //TSnap::GetWccSzCnt(inputg, CntV);
        //        std::cout << CntV.Len() <<std::endl;
        //        int sum = 0;
        //        for (int i = 0; i< CntV.Len(); i++){
        //            auto pair = CntV[i];
        //            std:: cout<< pair.GetVal1() << " "<< pair.GetVal2() << std::endl;
        ////            sum += pair.GetVal1().Val * pair.GetVal2().Val;
        //        }
        ////        std::cout << "sum : "<<sum; // checking the distribution
        //
        // Select starclusters
        //        std::vector <int> centers = select_cluster_centers(inputg, N_starclusters,"degree");
        
        for (int _iter = 1; _iter <= iterations; _iter++){
            
            // Measure Time
            std::clock_t c_start = std::clock();
            std::clock_t tot_start = std::clock();
            
            TIntV landmarks(1024);
            hash_map <int, int> dist_nearest_nbr_inL;
            
            
            
            if (connected_)
                landmarks = select_landmarks(inputg, heuristic, dist_nearest_nbr_inL, num_landmarks, epsilon); // it will populate partial dist_nearest_nbr array
            else{
                TCnComV components;
                TSnap::GetSccs(inputg,components);
                for(int i=0; i< components.Len(); i++){
                    auto comp = components[i];
                    PUNGraph SubG = TSnap::GetSubGraph(inputg, comp.NIdV);
                    //                std::cout << SubG->GetNodes() << " size\n";
                    TIntV land_ = select_landmarks(SubG, heuristic, dist_nearest_nbr_inL, num_landmarks, epsilon);
                    landmarks.AddV(land_);
                    //                std::cout << i+1 <<"-th component\n";
                }
            }
            std::clock_t c_end = std::clock();
            
            
            // print for each v, its distance to nearest neighbor
            //                    for (auto NI = inputg->BegNI(); NI < inputg->EndNI(); NI++) {
            //                        const int NId = NI.GetId();
            //                        std::cout << NId << " -> "<<dist_nearest_nbr_inL[NId]<<"\n";
            //                    }
            std::cout << " Size of dist_nearest_nbr hashmap => "<<dist_nearest_nbr_inL.size()<< " Supposed to be -> "<<inputg->GetNodes()<<"\n";
            
            std::cout << "Landmark size = "<<landmarks.Len()<<" / "<< inputg->GetNodes()<<" \n";
            //        for (int i=0; i<landmarks.Len(); i++ ){ std::cout << landmarks[i] << "\n"<< std::endl;}
            std::string output_directory;
            if (heuristic.rfind("eps", 0) == 0){
                stringstream stream;
                stream << fixed << setprecision(2) << epsilon;
                output_directory = heuristic+"_"+path2[path2.size()-2]+"_"+stream.str();
            }
            else    output_directory = heuristic+"_"+path2[path2.size()-2]+"_"+std::to_string(num_landmarks);
            
                    mkdir(output_directory.c_str(), S_IRWXU);
            std::string outfilename;
            //for(int j = 0; j < path.size()-1; j++)  outfilename+=(path[j]+"/");
            outfilename+=output_directory+"/";
            outfilename+=path2[0];
            outfilename+=("_"+std::to_string(_iter));
            std::cout << outfilename<<"\n";
            
            
            // compute the k x N matrix
            int cluster_length = landmarks.Len();
            int flat_matrix_length = ((cluster_length%2==0?cluster_length/2:cluster_length)*((cluster_length-1)%2==0?(cluster_length-1)/2:(cluster_length-1)));
            if (flat_matrix_length>0){
                std::vector <value_t> flat_matrix_lower_tri;
                flat_matrix_lower_tri.reserve(flat_matrix_length);
                construct_landmark_distmat(inputg,flat_matrix_lower_tri,dist_nearest_nbr_inL,landmarks);
                //        for(auto i:flat_matrix_lower_tri){
                //            std::cout <<i<<" ";
                //        }
                //        std::cout <<"\n";
                
                threshold = *std::max_element(std::begin(flat_matrix_lower_tri), std::end(flat_matrix_lower_tri)) + 1;
                // Pass the distance matrix to ripser to expand and compute PH
                compressed_lower_distance_matrix dist = compressed_lower_distance_matrix(std::move(flat_matrix_lower_tri));
                //            dist.print();
                int d_sz = dist.size();
                std::cout << "threshold "<<threshold<<" max dim = "<<dim_max<<"\n";
                std::vector<std::vector<std::pair<value_t,value_t>>> pIntervals;
                pIntervals.reserve(dim_max);
                pIntervals = ripser<compressed_lower_distance_matrix>(std::move(dist), dim_max, threshold, ratio,
                                                                      modulus).compute_barcodes(outfilename);
                print_barcode(pIntervals,std::min(dim_max, index_t(d_sz- 2)),outfilename);
            }
            else{
                // Only have 1 vertex. Then only have (0,Inf] in H_0.
                // TO DO: Write that in a file for H_0.
            }
            
            
            std::clock_t tot_end = std::clock();
            long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
            
            std::cout << "(Landmark selection) CPU time used: " << time_elapsed_ms << " ms\n";
            long double totaltime_elapsed_ms = 1000.0 * (tot_end-tot_start) / CLOCKS_PER_SEC;
            std::cout << "(Total Execution) CPU time used: " << totaltime_elapsed_ms << " ms\n";
            
            ofstream myfile (output_directory+"/"+output_directory+"_"+std::to_string(_iter)+".time");
            
            if (myfile.is_open())
            {
                //            for (index_t i = 1; i < argc; ++i) {
                //                myfile << argv[i]<< " ";
                //            }
                myfile << "Land_time "<< time_elapsed_ms<<" ";
                myfile << "Tot_time "<< totaltime_elapsed_ms<<" ";
                myfile << "#Land "<<landmarks.Len()<<" ";
                myfile << "graphV "<<inputg->GetNodes()<<" ";
                myfile << "graphE "<<inputg->GetEdges()<<" ";
                //            myfile << "Approx_Diam "<< TSnap::GetBfsFullDiam(inputg,inputg->GetNodes()/log(inputg->GetNodes()))<<"\n";
            }
            myfile.close();
            
        }
        
        //        }
        
        //        //        TIntVV
        //        std::vector <TIntV> clusters(N_starclusters);
        //        TBreathFS<PUNGraph> BFS(inputg, true, false);
        //        run_starclustering(inputg,centers,clusters,BFS);
        //
        //        // Cluster statistics
        //        //        int sum = 0;
        //        //        for (int i = 0; i< N_starclusters; i++){
        //        //            sum += clusters[i].Len();
        //        //            std::cout << "Cluster "<< i << " : " << clusters[i].Len() << std::endl;
        //        //        }
        //        //        std:: cout << sum << std::endl;
        //        std::vector<std::vector<std::vector<std::pair<value_t,value_t>>>> all_intervals(N_starclusters);
        //        //TSnap::PrintInfo(inputg, filename, "", inputg->GetNodes()>Kilo(1));
        //        for (int i = 0; i< N_starclusters; i++){
        //            int cluster_length = clusters[i].Len();
        //            int flat_matrix_length = ((cluster_length%2==0?cluster_length/2:cluster_length)*((cluster_length-1)%2==0?(cluster_length-1)/2:(cluster_length-1)));
        //            //            std::cout << flat_matrix_length << " flat length"<<std::endl;
        //            std::vector <value_t> flat_matrix;
        //            flat_matrix.reserve(flat_matrix_length);
        //            if (cluster_length!=0){ // Don't process empty clusters
        //                std::cout << "cluster "<<i <<" "<< " constructing distance matrix"<< std::endl;
        //                if(cluster_length==1){
        //                    std::vector<std::pair<value_t,value_t>> vec_pair({std::make_pair<value_t,value_t>(0,std::numeric_limits<value_t>::max())});
        //                    all_intervals[0]=std::vector<std::vector<std::pair<value_t,value_t>>>({vec_pair});
        //                    continue;
        //                }
        //                // Construct flat matrix and threshold
        //                value_t threshold_computed;
        //                construct_flatdist_matrix(inputg, clusters[i], flat_matrix);
        //                //                for (auto it:flat_matrix){
        //                //                    std::cout << it << " ";
        //                //                }
        //                //                std::cout << std::endl;
        //                std::cout << "threshold_computed: "<<threshold_computed<<" threshold: "<<threshold<<std::endl;
        //                //                threshold = std::min(threshold_computed, threshold);
        //
        //                value_t max_dist = *std::max_element(std::begin(flat_matrix), std::end(flat_matrix));
        //
        //                // std::cout << "threshold : "<<threshold << " max/largest diameter: "<< max_dist<<std::endl;
        //                compressed_lower_distance_matrix dist = compressed_lower_distance_matrix(std::move(flat_matrix));
        //
        //                //                for (auto it : dist.distances)
        //                //                    std::cout <<it << std::endl;
        //                // Compute threshold
        //                if (threshold == std::numeric_limits<value_t>::max()) {
        //                    value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
        //                    for (index_t i = 0; i < dist.size(); ++i) {
        //                        value_t r_i = -std::numeric_limits<value_t>::infinity();
        //                        for (index_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
        //                        enclosing_radius = std::min(enclosing_radius, r_i);
        //                    }
        //
        //                    threshold = enclosing_radius;
        //                }
        //                std::cout << "threshold : "<<threshold << " max/largest diameter: "<< max_dist<< std::endl;
        //                std::cout << "subgraph size: "<< cluster_length<<std::endl;
        //                // Compute barcode
        //                all_intervals[i] = ripser<compressed_lower_distance_matrix>(std::move(dist), dim_max, threshold, ratio,
        //                                                                            modulus).compute_barcodes(outfilename,std::to_string(i));
        //
        //            }
        //            //            for(int j =0; j< clusters[i].Len();j++){
        //            //                std::cout << clusters[i][j] << " ";
        //            //            }
        //            //            std::cout << std::endl;
        //
        //        }
        
        
        
        
    }
    else if(format == WGRAPH){
        //        std::string heuristic = "epsnetring";
        //                std::string heuristic = "epsrand";
        //        std::string heuristic = "random";
        //                std::string heuristic = "eps_baseline";
        //        std::string heuristic = "eps_kruskal";
        //                std::string heuristic = "eps_densedfs";
        //        std::string heuristic = "eps_sptree_prune";
        //        std::string heuristic = "eps_filtering";
        //         std::string heuristic = "epsmaxmin";
        
        std::cout <<"fname "<< filename << std::endl;
        // Read weighted graph.
        std::ifstream file_stream(filename);
        TPt <TNodeEDatNet<TInt, TFlt> > Net = readGraph<TInt,TFlt>(filename ? file_stream : std::cin);
        connected_ = TSnap::IsConnected(Net);
        std::cout << " * * * * * \n";
        std::cout << (connected_?"Connected":"Not connected") << " \n* * * * * * \n";
        
        for (int _iter = 1; _iter <= iterations; _iter++){
            std::clock_t c_start = std::clock();
            std::clock_t tot_start = std::clock();
            TIntV landmarks(1024);
            hash_map <int, double> dist_nearest_nbr_inL;
            dist_nearest_nbr_inL.reserve(Net->GetNodes());
            
            
            /* Select the landmarks (connected/not connected network) */
            if (connected_)
                landmarks = weighted_select_landmarks(Net, heuristic, dist_nearest_nbr_inL, num_landmarks, epsilon); // it will populate partial dist_nearest_nbr array
            else{
                TCnComV components;
                TSnap::GetSccs(Net,components);
                for(int i=0; i< components.Len(); i++){
                    auto comp = components[i];
                    auto subg_vertices = comp.NIdV;
                    /* Construct weighted subgraph */
                    TPt<TNodeEDatNet<TInt, TFlt>> pSubG = TNodeEDatNet<TInt, TFlt>::New();
                    TIntSet NIdSet(subg_vertices.Len());
                    for (int i = 0; i < subg_vertices.Len(); i++) {
                        if(Net->IsNode(subg_vertices[i])){
                            pSubG->AddNode(subg_vertices[i]);
                            NIdSet.AddKey(subg_vertices[i]);
                        }
                    }
                    //  std::cout << "Adding edges\n\n";
                    for (int n = 0; n < NIdSet.Len(); n++) {
                        const int SrcNId = NIdSet[n];
                        // std::cout << "xxxx ---- "<<SrcNId << " ";
                        auto NI = Net->GetNI(SrcNId);
                        for (int edge = 0; edge < NI.GetOutDeg(); edge++) {
                            const int OutNId = NI.GetOutNId(edge);
                            if (NIdSet.IsKey(OutNId)) {
                                pSubG->AddEdge(SrcNId, OutNId);
                                //pSubG->AddEdge(OutNId,SrcNId);
                            }
                        }
                    }
                    
                    for (auto EI = pSubG->BegEI(); EI < pSubG->EndEI(); EI++) {
                        pSubG->SetEDat(EI.GetSrcNId(),EI.GetDstNId(),Net->GetEDat(EI.GetSrcNId(),EI.GetDstNId()));
                    }
                    
                    std::cout << pSubG->GetNodes() << " size\n";
                    TIntV land_ = weighted_select_landmarks(pSubG, heuristic, dist_nearest_nbr_inL, num_landmarks, epsilon);
                    landmarks.AddV(land_);
                    std::cout << i+1 <<"-th component\n";
                }
            }
            std::clock_t c_end = std::clock();
            
            //        for (auto NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
            //            const int NId = NI.GetId();
            //            std::cout << NId << " -> "<<dist_nearest_nbr_inL[NId]<<"\n";
            //        }
            std::cout << " Size of dist_nearest_nbr hashmap => "<<dist_nearest_nbr_inL.size()<< " Supposed to be -> "<<Net->GetNodes()<<"\n";
            
            std::cout << "Landmark size = "<<landmarks.Len()<<" / "<< Net->GetNodes()<<" \n";
            //        for (int i=0; i<landmarks.Len(); i++ ){ std::cout << landmarks[i] << "\n"<< std::endl;}
            
            
            std::set<char> delims{'/'};
            std::set<char> delims2{'.'};
            std::vector<std::string> path = splitpath(filename, delims);
            std::vector<std::string> path2 = splitpath(path.back(), delims2);
            //for(auto it:path) std::cout<<it<<"-";
            //  std::cout<<"\n";
            //std::cout<<path2[path2.size()-2]<<"\n";
            //std::string output_directory = heuristic+path2[path2.size()-2];
            std::string output_directory;
            if (heuristic.rfind("eps", 0) == 0){
                stringstream stream;
                stream << fixed << setprecision(2) << epsilon;
                output_directory = heuristic+"_"+path2[path2.size()-2]+"_"+stream.str();
            }
            else    output_directory = heuristic+"_"+path2[path2.size()-2]+"_"+std::to_string(num_landmarks);
            std::cout<< "making directory (if not exists): "<<output_directory<<"\n";
            mkdir(output_directory.c_str(), S_IRWXU);
            std::string outfilename;
            //for(int j = 0; j < path.size()-1; j++)  outfilename+=(path[j]+"/");
            outfilename+=output_directory+"/";
            outfilename+=path2[0];
            outfilename+=("_"+std::to_string(_iter));
            std::cout << "Writing barcodes to file: "<<outfilename<<"\n";
            
            
            // Dist matrix, Simplicial complex and Homology computation
            // compute the k x N matrix
            value_t diametr = -1; // stores largest shortest path.
            int cluster_length = landmarks.Len();
            int flat_matrix_length = ((cluster_length%2==0?cluster_length/2:cluster_length)*((cluster_length-1)%2==0?(cluster_length-1)/2:(cluster_length-1)));
            if (flat_matrix_length>0){
                std::vector <value_t> flat_matrix_lower_tri;
                flat_matrix_lower_tri.reserve(flat_matrix_length);
                weighted_construct_landmark_distmat(Net,flat_matrix_lower_tri,dist_nearest_nbr_inL,landmarks);
                //        for(auto i:flat_matrix_lower_tri){
                //            std::cout <<i<<" ";
                //        }
                //        std::cout <<"\n";
                
                threshold = *std::max_element(std::begin(flat_matrix_lower_tri), std::end(flat_matrix_lower_tri)) + 1;
                diametr = std::max(diametr,(threshold-1));
                // Pass the distance matrix to ripser to expand and compute PH
                compressed_lower_distance_matrix dist = compressed_lower_distance_matrix(std::move(flat_matrix_lower_tri));
                //            dist.print();
                int d_sz = dist.size();
                std::cout << "threshold "<<threshold<<" max dim = "<<dim_max<<"\n";
                std::vector<std::vector<std::pair<value_t,value_t>>> pIntervals;
                pIntervals.reserve(dim_max);
                pIntervals = ripser<compressed_lower_distance_matrix>(std::move(dist), dim_max, threshold, ratio,
                                                                      modulus).compute_barcodes(outfilename);
                print_barcode(pIntervals,std::min(dim_max, index_t(d_sz- 2)),outfilename);
            }
            else{
                // Only have 1 vertex. Then only have (0,Inf] in H_0.
                // TO DO: Write that in a file for H_0.
            }
            
            
            std::clock_t tot_end = std::clock();
            long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
            
            std::cout << "(Landmark selection) CPU time used: " << time_elapsed_ms << " ms\n";
            long double totaltime_elapsed_ms = 1000.0 * (tot_end-tot_start) / CLOCKS_PER_SEC;
            std::cout << "(Total Execution) CPU time used: " << totaltime_elapsed_ms << " ms\n";
            
            ofstream myfile (output_directory+"/"+output_directory+"_"+std::to_string(_iter)+".time");
            
            if (myfile.is_open())
            {
                //            for (index_t i = 1; i < argc; ++i) {
                //                myfile << argv[i]<< " ";
                //            }
                myfile << "Land_time "<< time_elapsed_ms<<" ";
                myfile << "Tot_time "<< totaltime_elapsed_ms<<" ";
                myfile << "#Land "<<landmarks.Len()<<" ";
                myfile << "graphV "<<Net->GetNodes()<<" ";
                myfile << "graphE "<<Net->GetEdges()<<" ";
                myfile << "Approx_Diam "<< diametr <<"\n";
            }
            myfile.close();
        }
        
    }
    else {
        string cluster_id = "";
        compressed_lower_distance_matrix dist =
        read_file(filename ? file_stream : std::cin, format);
        
        value_t min = std::numeric_limits<value_t>::infinity(),
        max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
        int num_edges = 0;
        
        if (threshold == std::numeric_limits<value_t>::max()) {
            value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
            for (index_t i = 0; i < dist.size(); ++i) {
                value_t r_i = -std::numeric_limits<value_t>::infinity();
                for (index_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
                enclosing_radius = std::min(enclosing_radius, r_i);
            }
            
            threshold = enclosing_radius;
        }
        
        for (auto d : dist.distances) {
            min = std::min(min, d);
            max = std::max(max, d);
            max_finite =
            d != std::numeric_limits<value_t>::infinity() ? std::max(max, d) : max_finite;
            if (d <= threshold) ++num_edges;
        }
        
        std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;
        
        if (threshold >= max) {
            std::cout << "distance matrix with " << dist.size() << " points" << std::endl;
            ripser<compressed_lower_distance_matrix>(std::move(dist), dim_max, threshold, ratio,
                                                     modulus)
            .compute_barcodes(filename);
        } else {
            std::cout << "sparse distance matrix with " << dist.size() << " points and "
            << num_edges << "/" << (dist.size() * dist.size() - 1) / 2 << " entries"
            << std::endl;
            
            ripser<sparse_distance_matrix>(sparse_distance_matrix(std::move(dist), threshold),
                                           dim_max, threshold, ratio, modulus)
            .compute_barcodes(filename);
        }
    }
}
