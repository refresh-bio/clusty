#pragma once 

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include <iostream>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <random>
#include <limits>
#include <cmath>
#include <memory>
#include <vector>
#include <sstream>

#include "clustering.h"
#include "distances.h"
#include "clock.h"
#include "utils.h" 
#include "memory_monotonic.h"
#include "log.h"

namespace linkage_algorithm_heaptrix
{
   /** Group is composed of two other groups.
    *  It is a part of a dengrogram. */
   struct group 
   {
      std::size_t id, left, right; // indices of left and right groups
      double distance; // distances between left and right group
   };

   /// Dendrogram 
   struct dendrogram
   {
      std::vector<group> groups;
   };

  inline std::ostream & operator << (std::ostream & sos, const group & g)
   {
      return sos << g.id << ": d(" << g.left << ", " << g.right << ") == " << g.distance;
   }

   inline std::ostream & operator << (std::ostream & sos, const dendrogram & d)
   {
      sos << "--------------------" << std::endl;
      for (const auto g : d.groups)
      {
         sos << g << std::endl;
      }
      sos << "--------------------" << std::endl;
      return sos;
   }


   template <class Distance, typename AggregationRule>
      class linkage : public HierarchicalClustering<Distance> 
   {
      const std::size_t MAX_SIZE_T = std::numeric_limits<std::size_t>::max(); // Index of the location in the heap. MAX means no index set.
      const double MAX_DOUBLE = std::numeric_limits<double>::max(); 
      const double INF_DOUBLE = std::numeric_limits<double>::infinity(); 
      const std::string MAX_LABEL {"--"};

      /*** An element of a matrix and a heap simultaneously.  */
      #pragma pack(push, 8)
      struct element
      {
          double _value;              ///< stored value
          int32_t _row;           ///< index of a row in the matrix
          int32_t _column;        ///< index of a column in the matrix
          std::size_t _index_heap;    ///< index of an element in the heap

         element (const int32_t row, const int32_t column, const std::size_t index_heap, const double value)
         {
            _row = row;
            _column = column;
            _index_heap = index_heap;
            _value = value;
         }

         /*
            auto operator <=> (const element & p)
            {
            return this->_value <=> p._value;
            }
            */

         auto operator < (const element & p) const
         {
             if(this->_value != p._value)
                 return this->_value < p._value;

             if(this->_row != p._row)
                 return this->_row < p._row;

             return this->_column < p._column;
         }

         auto operator ==(const element& p) const
         {
            return (this->_value == p._value) && (this->_row == p._row) && (this->_column == p._column);
         }

         void print (std::ostream & sos)
         {
            sos << "{(" << _row << ", " << _column << ") [" << _index_heap << "] " << _value << "}";
         }
      };
      #pragma pack(pop)

      class matrix_row_ht
      {
      public:
          using ht_elem_t = std::pair<std::size_t, element*>;

      private:
          static constexpr size_t empty = std::numeric_limits<std::size_t>::max();
          static constexpr size_t removed = std::numeric_limits<std::size_t>::max() - 1;

          std::vector<ht_elem_t> ht;

          double max_fill_factor;
          size_t ht_size;
          size_t ht_filled;
          size_t ht_mask;
          size_t when_restruct;

          std::size_t determine_size(std::size_t req_size)
          {
              std::size_t s = (std::size_t)(req_size / max_fill_factor);
              while (s & (s - 1))
                  s &= s - 1;

              return 2 * s;
          }

          void restruct()
          {
              std::vector<ht_elem_t> ht_old(std::move(ht));

              ht.clear();
              ht_size *= 2;
              ht_mask = ht_size - 1;
              when_restruct = (std::size_t)(ht_size * max_fill_factor);
              ht_filled = 0;

              ht.resize(ht_size, std::make_pair(empty, nullptr));

              for (const auto& [idx, ptr] : ht_old)
                  if (ptr)
                      insert(idx, ptr);
          }

          // MurMur3 hash
          std::size_t hash_mm(std::size_t x) const
          {
              x ^= x >> 33;
              x *= 0xff51afd7ed558ccdLL;
              x ^= x >> 33;
              x *= 0xc4ceb9fe1a85ec53LL;
              x ^= x >> 33;

              return x;
          }

      public:
          matrix_row_ht(std::size_t init_size = 8, double max_fill_factor = 0.8) :
              max_fill_factor(max_fill_factor)
          {
              ht_size = determine_size(init_size) / 2;
              ht_mask = ht_size - 1;

              // Enforce init at first insert
              when_restruct = 0;
              ht_filled = 0;
          }

          matrix_row_ht(matrix_row_ht&& rhs) noexcept
          {
              ht = move(rhs.ht);

              max_fill_factor = rhs.max_fill_factor;
              ht_size = rhs.ht_size;
              ht_filled = rhs.ht_filled;
              ht_mask = rhs.ht_mask;
              when_restruct = rhs.when_restruct;
          }

          matrix_row_ht& operator=(matrix_row_ht&& rhs) 
          {
              if (this == &rhs)
                  return *this;

              ht = move(rhs.ht);

              max_fill_factor = rhs.max_fill_factor;
              ht_size = rhs.ht_size;
              ht_filled = rhs.ht_filled;
              ht_mask = rhs.ht_mask;
              when_restruct = rhs.when_restruct;

              return *this;
          }

          matrix_row_ht& operator=(const matrix_row_ht& rhs) 
          {
              if (this == &rhs)
                  return *this;

              ht = rhs.ht;

              max_fill_factor = rhs.max_fill_factor;
              ht_size = rhs.ht_size;
              ht_filled = rhs.ht_filled;
              ht_mask = rhs.ht_mask;
              when_restruct = rhs.when_restruct;

              return *this;
          }

          void insert(size_t idx, element* ptr)
          {
              if (ht_filled == when_restruct)
                  restruct();

              auto pos = hash_mm(idx) & ht_mask;

              while (true)
              {
                  if (ht[pos].second == nullptr)
                      break;

                  pos = (pos + 1) & ht_mask;
              }        

              ht[pos].first = idx;
              ht[pos].second = ptr;

              ++ht_filled;
          }

          typename std::vector<matrix_row_ht::ht_elem_t>::iterator find(size_t idx)
          {
              if (ht_filled == 0)
                  return ht.end();

              auto pos = hash_mm(idx) & ht_mask;

              while (true)
              {
                  if (ht[pos].first == idx)
                      return ht.begin() + pos;

                  if (ht[pos].first == empty)
                      return ht.end();

                  pos = (pos + 1) & ht_mask;
              }
          }

          void erase(typename std::vector<matrix_row_ht::ht_elem_t>::iterator it)
          {
              it->first = removed;
              it->second = nullptr;
          }

          void erase(size_t idx)
          {
              if (ht_filled == 0)
                  return;

              auto pos = hash_mm(idx) & ht_mask;

              while (true)
              {
                  if (ht[pos].first == idx)
                  {
                      ht[pos].first = removed;
                      ht[pos].second = nullptr;
                      return;
                  }

                  if (ht[pos].first == empty)
                      return;

                  pos = (pos + 1) & ht_mask;
              }
          }

          typename std::vector<matrix_row_ht::ht_elem_t>::iterator begin()
          {
              return ht.begin();
          }

          typename std::vector<matrix_row_ht::ht_elem_t>::iterator end()
          {
              return ht.end();
          }

          void clear()
          {
              ht_size = 8;
              when_restruct = 0;
              ht_filled = 0;
              ht.clear();
              ht.shrink_to_fit();
          }

          void prefetch()
          {
              _my_prefetch(ht.data());
          }
      };

      /*** row of the matrix */
      struct matrix_row
      {
         std::size_t id; // index
         matrix_row_ht _values; 

         matrix_row() : id(std::numeric_limits<std::size_t>::max()) {  }
         matrix_row(std::size_t id) : id(id) {  }

         void clear() { _values.clear(); id = std::numeric_limits<std::size_t>::max(); }

         matrix_row(const matrix_row& rhs)
         {
             id = rhs.id;
             _values = rhs._values;
         }

         matrix_row(matrix_row&& rhs) noexcept
         {
             id = std::move(rhs.id);
             _values = std::move(rhs._values);
         };

         matrix_row& operator=(matrix_row&&) noexcept = default;
         matrix_row& operator=(const matrix_row&) = default;
      };

            //------------------- 
      struct matrix
      {
          const std::size_t MAX_SIZE_T = std::numeric_limits<std::size_t>::max(); 
          
          std::vector<matrix_row> _rows;

         /** @return true, if the w-th row exists in the matrix.
          */ 
         bool exists_row(const std::size_t w)
         {
             if (w >= _rows.size())
                 return false;

             return _rows[w].id != MAX_SIZE_T;
         }

         /** @return true, if the cell in the row-th row and col-th column exists. */
         bool exists_cell (const std::size_t row, const std::size_t col)
         {
            if (not exists_row(row))
               return false;
            else
               return _rows[row]._values.find(col) != _rows[row]._values.end();
         }

         bool add_row(const std::size_t w)
         {
             if (w < _rows.size())
             {
                 if (_rows[w].id != MAX_SIZE_T)
                 {
                     std::cerr << "Row " << w << " already exists" << std::endl;
                     return false;
                 }
                 _rows[w].id = w;
                 return true;
             }
                 
             _rows.resize((std::size_t)(std::max<size_t>(w, 16) * 1.2), MAX_SIZE_T);
             _rows[w].id = w;

             return true;
         }

         void prefetch_row(const std::size_t w)
         {
             if (w >= _rows.size())
                 return;

             _my_prefetch(_rows.data() + w);
         }

         void prefetch_row_data(const std::size_t w)
         {
             if (w >= _rows.size())
                 return;

             _rows[w]._values.prefetch();
         }

         bool remove_row(const std::size_t w)
         {
             if (w >= _rows.size())
                 return false;

             _rows[w].clear();

             return true;
         }

         size_t get_max_row() const
         {
             for (auto p = _rows.rbegin(); p != _rows.rend(); ++p)
                 if (p->id != MAX_SIZE_T)
                     return p->id;

             return MAX_SIZE_T;
         }
      };
      //------------------- 

      class heap
      {
         /// Heap is index from 0.
         /// parent  --> children
         /// i       --> 2 * i + 1  and  2 * i + 2
         /// child   --> parent
         /// i       --> floor((i - 1) / 2)
         std::vector<element*> _heap;

         bool less (const element* l, const element* p)
         {
            return *l < *p;
         }

         public:
         std::size_t size() const 
         {
            return _heap.size();
         }

         void reserve(std::size_t req_size)
         {
             _heap.reserve(req_size);
         }

         public:
         void clear()
         {
            _heap.clear();
         }

         public:
         /// The method ONLY pushes back an item into the _heap vector WITHOUT sifting. 
         void push_back(element* p)
         {
            _heap.push_back(p);
         }

         protected:
         REFRESH_FORCE_INLINE void sift_up (const std::size_t last)
         {
            auto child = last;

            while (child > 0)
            {
               auto parent = (child - 1) / 2;

               if(child > 0)
                    _my_prefetch(_heap.data() + (child - 1) / 2);

               if (less(_heap[child], _heap[parent]))
               {
                  std::swap(_heap[child], _heap[parent]);
                  std::swap(_heap[child]->_index_heap, _heap[parent]->_index_heap);
                  child = parent;
               }
               else 
                  return;
            }
         }

         protected:
         REFRESH_FORCE_INLINE void sift_down (const std::size_t first)
         {
            auto parent = first;
            auto last = _heap.size() - 1;

            if (parent * 2 + 1 <= last)
                _my_prefetch(_heap.data() + parent * 2 + 1);

            while (2 * parent + 1 <= last) // While parent has a child.
            {
               auto child1 = 2 * parent + 1;
               auto child2 = 2 * parent + 2;
               // Which child is smaller?
               auto smaller_child = child1;

               if(child1 * 2 + 1 <= last)
                   _my_prefetch(_heap.data() + child1 * 2 + 1);

               _my_prefetch(_heap[child1]);
               if (child2 <= last)
                   _my_prefetch(_heap[child2]);

               if (child2 <= last)
               {
                  smaller_child = less(_heap[child1], _heap[child2]) ? child1 : child2;
               }
               
               if (2 * smaller_child + 1 <= last)
                   _my_prefetch(&_heap[2 * smaller_child + 1]);

               if (less (_heap[smaller_child], _heap[parent]))
               {
                  std::swap(_heap[smaller_child], _heap[parent]);
                  std::swap(_heap[smaller_child]->_index_heap, _heap[parent]->_index_heap);
                  parent = smaller_child;
               }
               else
                  return ;
            }
         }

         public:
         void make_heap ()
         {
            std::size_t start = (_heap.size() + 1) / 2;
            bool last = false;
            while (not last)
            {
               sift_down(start);
               if (start == 0)
                  last = true;
               else
                  start--;
            }
            // wpisanie indeksow:
            for (std::size_t i = 0; i < _heap.size(); ++i)
            {
               _heap[i]->_index_heap = i;
            }
         }

         public:
         void print_heap() const 
         {
            std::cout << "-------------" << std::endl;
            for (const auto & p : _heap)
            {
               p->print(std::cout); 
               std::cout << " ";
            }
            std::cout << std::endl;
            std::cout << "-------------" << std::endl;
         }

         public:
         void print_heap_indices() const 
         {
            std::cout << "-------------" << std::endl;
            for (const auto & p : _heap)
               std::cout << p->_index_heap << " ";
            std::cout << std::endl;
            std::cout << "-------------" << std::endl;
         }

         public:
         bool empty () const 
         {
            return _heap.empty();
         }

         public:
         /// @return element on the top (does not remove it from the heap).
         element* top()
         {
            return _heap.front();
         }

         public:
         void insert(element* p)
         {
            _heap.push_back(p);
            p->_index_heap = _heap.size() - 1;
            sift_up(p->_index_heap);
         }

         public:
         void remove(const element* p)
         {
            if (p == nullptr)
               return;

            std::size_t index = p->_index_heap;

            if (index >= _heap.size())
            {
               std::stringstream sos;
               sos << "incorrect index: index (" << index << ") out of range [0, " << _heap.size() - 1 << "]"; 
               throw std::string (sos.str());
            }

            if (_heap.size() > 1)
            {
             //   auto old_value = _heap[index]->_value;

                bool move_up = less(_heap.back(), p);

                _heap[index] = _heap.back();
                _heap[index]->_index_heap = index;
                _heap.pop_back();

                if(move_up)
                    sift_up(index);
                else
                    sift_down(index);

               _my_prefetch(_heap.back());
            }
            else if (_heap.size() == 1)
            {
               _heap.pop_back();
            }
         }

         void replace(const element* p_old, element* p_new)
         {
             auto index = p_old->_index_heap;

             _heap[index] = p_new;
             p_new->_index_heap = index;

             if (less(p_new, p_old))
                 sift_up(index);
             else
             {
                 sift_down(index);
                 _my_prefetch(_heap.back());
             }
         }

         public:
         /// @return element on the top and remove it from the heap
         element* pop()
         {
            element* returned_value = nullptr;
            if (_heap.size() > 1)
            {
               returned_value = _heap.front();
               std::swap(_heap[0], _heap[_heap.size() - 1]);
               std::swap(_heap.front()->_index_heap, _heap.back()->_index_heap);
               _heap.pop_back();
               sift_down(0);
            }
            else if (_heap.size() == 1)
            {
               returned_value = _heap.front();
               _heap.pop_back();

            }
            return returned_value;
         }
      };

      matrix _matrix;
      dendrogram _dendrogram;
      heap _heap;
      refresh::memory_monotonic_unsafe *mma;
      std::vector<element*> mma_buf;

      element* element_allocate(const std::size_t row, const std::size_t column, const std::size_t index_heap, const double value)
      {
          element* ptr;

          if(mma_buf.empty())
              ptr = (element *) mma->allocate(sizeof(element));
          else
          {
              ptr = mma_buf.back();
              mma_buf.pop_back();
          }

          ptr->_row = row;
          ptr->_column = column;
          ptr->_index_heap = index_heap;
          ptr->_value = value;

          return ptr;
      }

      void element_free(element* ptr)
      {
          if (mma_buf.size() == mma_buf.capacity())
              mma_buf.reserve(2 * mma_buf.size());

          mma_buf.emplace_back(ptr);
      }

      public:
      void print_indices_of_matrix_rows()
      {
         std::cout << "---" << std::endl;

         for (const auto & [row_id, row] : _matrix._rows)
         {
            std::cout << row_id << " ";
         }
         std::cout << std::endl;
         std::cout << "---" << std::endl;
      }

      public:
      void print_matrix()
      {
         std::cout << "---" << std::endl;

         for (const auto & [row_id, row] : _matrix._rows)
         {
            std::cout << row_id << ": ";
            auto it  = row._value.begin();
            auto end = row._value.end();

            std::size_t k = 0;
            while (it != end)
            {
               while (it->second->_column != k)
               {
                  std::cout << MAX_LABEL << " ";
                  ++k;
               }
               std::cout << it->second->_value << " ";
               ++k;
               ++it;
            }
            std::cout << std::endl;
         }
         std::cout << "---" << std::endl;
      }

      public:
      void add_value (std::size_t w, std::size_t k, const double _value)
      {
          if (w > k)
              std::swap(w, k);

         if (not _matrix.exists_cell(w, k))
         {
            if (not _matrix.exists_row(w))
                _matrix.add_row(w);

            if (not _matrix.exists_row(k))
                _matrix.add_row(k);

            element* pEl = element_allocate(w, k, _heap.size(), _value);
            _matrix._rows[w]._values.insert(k, pEl);
            _matrix._rows[k]._values.insert(w, pEl);
            _heap.push_back(pEl);
         }
      }

      public:
      void read_matrix (SparseMatrix<Distance> & m)
      {
         _matrix._rows.clear();
       
         _heap.reserve((std::size_t)(1.1 * m.num_elements()));

         for (size_t i = 0; i < m.num_objects(); ++i) {
             
             for (const Distance* p = m.begin(i); p < m.end(i); ++p) {
                 const Distance& edge = *p;
                 if (edge.get_d() != MAX_DOUBLE and edge.get_d() != INF_DOUBLE and i != edge.get_id())
                     add_value(i, edge.get_id(), edge.get_d());
             }

             m.clear_row(i);
         }

         _heap.make_heap();
      };

      protected:
      AggregationRule aggregation;

      protected:
      void do_clustering()
      {
         std::size_t id_of_the_next_group = _matrix.get_max_row() + 1; // id of aggregated group
         std::size_t number_of_objects = id_of_the_next_group;

         //std::cerr << std::endl;

         std::vector<element*> merged_column;
         std::vector<element*> heap_insert_buffer;

         while (number_of_objects > 1 and not _heap.empty())
         {
          //   if(number_of_objects % 1000 == 0)
          //       std::cerr << std::to_string(number_of_objects) + "     \r";

            // find the minimal distance
            auto pMinimal = _heap.top(); // I am not removing it from the heap. I am only reading the value. It will be removed in row merger.
            std::size_t r_min = pMinimal->_row;    // index of a row to merge
            std::size_t c_min = pMinimal->_column; // index of a row to merge
            double minimal_distance = pMinimal->_value;     // minimal distance

            // r_min should be less than k_min – for easier merging
            if (r_min > c_min)
               std::swap(r_min, c_min);

            // Create a new group in the dendrogram.
            group g {id_of_the_next_group, _matrix._rows[r_min].id, _matrix._rows[c_min].id, minimal_distance };
            _dendrogram.groups.push_back(g); 

            ///////////////////////
            // merge matrix rows
            matrix_row merged_row(id_of_the_next_group);  

            merged_column.clear();

            heap_insert_buffer.clear();

            // Iterate all values in r_min row and merge:
            for (const auto & [column_id, p] : _matrix._rows[r_min]._values)
            {
                if (!p)                 // need for custom HT
                    continue;

                if (column_id == c_min)
                    continue;

               // Is the a counterpart in c_min row?
               auto it = _matrix._rows[c_min]._values.find(column_id);
               if (it != _matrix._rows[c_min]._values.end())
               {
                  // Yes, there is – merge!
                  double merged = aggregation (p->_value, (*it).second->_value);
                  auto pNew = element_allocate(column_id, id_of_the_next_group, MAX_SIZE_T, merged); // new element for the merged column
                  heap_insert_buffer.emplace_back(pNew);

                  merged_row._values.insert(column_id, pNew); // insert into the matrix
                  merged_column.push_back(pNew);
               }
               else // it == end
               {
                  // There is no finite counterpart.
                  double merged = aggregation (p->_value, MAX_DOUBLE);
                  if (merged != MAX_DOUBLE)
                  {
                     auto pNew = element_allocate(column_id, id_of_the_next_group, MAX_SIZE_T, merged);
                     heap_insert_buffer.emplace_back(pNew);

                     merged_row._values.insert(column_id, pNew); // insert into the matrix
                     merged_column.push_back(pNew);
                  }
               }
            }
            // Iterate all values in c_min row and merge: 
            for (const auto & [column_id, p] : _matrix._rows[c_min]._values)
            {
                if (!p)                 // need for custom HT
                    continue;

                if (column_id == r_min)
                    continue;
                
                // Is the a counterpart in c_min row?
               auto it = _matrix._rows[r_min]._values.find(column_id);
               if (it == _matrix._rows[r_min]._values.end())
               {
                  // If there is a counterpart, it has already been merged in the loop above.
                  // If there is no counterpart, merge it!

                  // There is no finite counterpart.
                  double merged = aggregation (MAX_DOUBLE, p->_value);
                  if (merged != MAX_DOUBLE)
                  {
                     auto pNew = element_allocate(column_id, id_of_the_next_group, MAX_SIZE_T, merged);
                     heap_insert_buffer.emplace_back(pNew);

                     merged_row._values.insert(column_id, pNew); // insert into the matrix
                     merged_column.push_back(pNew);
                  }
               }
            }
            // add merged row
            _matrix.add_row(id_of_the_next_group);
            _matrix._rows[id_of_the_next_group] = std::move(merged_row);

            // add values in the merged column
            for (const auto p : merged_column)
               _matrix._rows[p->_row]._values.insert(p->_column, p);

            _matrix.prefetch_row(r_min);
            _matrix.prefetch_row(c_min);

            // remove rows r_min and c_min
            // first remove from the heap
            for (const auto & [_, p] : _matrix._rows[r_min]._values)
            {
                if (!p)                 // need for custom HT
                    continue;

                auto c_row = r_min != p->_column ? p->_column : p->_row;

                if (!heap_insert_buffer.empty())
                {
                    _heap.replace(p, heap_insert_buffer.back());
                    heap_insert_buffer.pop_back();
                }
                else
                    _heap.remove(p);

               _matrix._rows[c_row]._values.erase(r_min);

               element_free(p);
            }

            _matrix.prefetch_row_data(r_min);
            _matrix.prefetch_row_data(c_min);

            for (const auto & [_, p] : _matrix._rows[c_min]._values)
            {
                if (!p)                 // need for custom HT
                    continue;
                
                auto c_row = c_min != p->_column ? p->_column : p->_row;

                if (!heap_insert_buffer.empty())
                {
                    _heap.replace(p, heap_insert_buffer.back());
                    heap_insert_buffer.pop_back();
                }
                else
                    _heap.remove(p);

                _matrix._rows[c_row]._values.erase(c_min);

                element_free(p);
            }

            // then remove from the matrix
            _matrix.remove_row(r_min);
            _matrix.remove_row(c_min);

            ////////////////////////////////////
            number_of_objects--; // Two merged groups are removed and one result group added.
            id_of_the_next_group++; 

            for (const auto e : heap_insert_buffer)
                _heap.insert(e);
            heap_insert_buffer.clear();
         }

         mma->release();
         mma_buf.clear();
         mma_buf.shrink_to_fit();

         //std::cerr << std::to_string(number_of_objects) + "     \n";
      }

      public:
      /*
      dendrogram do_clustering (const std::vector<std::vector<double>> & m)
      {
         ksi::clock stopwatch;
         std::cout << "Loading data into data structure ";
         stopwatch.start();
         read_matrix(m);
         stopwatch.stop(); 
         double elapsed_time = (double) stopwatch.elapsed_milliseconds() / 1000;
         std::cout << "done in " << elapsed_time << " s. ";
         std::cout << "Clustering ";
         stopwatch.start();
         do_clustering();
         stopwatch.stop();
         elapsed_time = (double) stopwatch.elapsed_milliseconds() / 1000;
         std::cout << "done in " << elapsed_time << " s. ";
         return _dendrogram;
      }
      */

      public:
      dendrogram do_clustering (SparseMatrix<Distance> & m)
      {
         ksi::clock stopwatch;
         LOG_VERBOSE << "Loading data into heap ";
         stopwatch.start();
         read_matrix(m);
         stopwatch.stop(); 
         double elapsed_time = (double) stopwatch.elapsed_milliseconds() / 1000;
         LOG_VERBOSE << "done in " << elapsed_time << " s. ";
         LOG_VERBOSE << "Performing linkage ";
         stopwatch.start();
         do_clustering();
         stopwatch.stop();
         elapsed_time = (double) stopwatch.elapsed_milliseconds() / 1000;
         LOG_VERBOSE << "done in " << elapsed_time << " s. ";
         return _dendrogram;
      }

      protected:
      std::vector<node_t> makeDendrogram(const dendrogram & d, const std::size_t nObjects)
      {
         std::vector<node_t> dendro (nObjects);

         for (const auto & g : d.groups) // for each group
         {
            dendro.emplace_back(g.left, g.right, g.distance);
         }

         return dendro;
      }

      public:
      int run (
            SparseMatrix<Distance>& matrix,
            const std::vector<int>& objects,
            double threshold,
            std::vector<int>& assignments
            )
      {
         auto dendrogram = do_clustering(matrix);
         std::vector<node_t> node_t_dendrogram = makeDendrogram(dendrogram, objects.size());
         assignments.resize(objects.size());
         int n_clusters = this->dendrogramToAssignments(node_t_dendrogram, threshold, assignments);

         return n_clusters;
      }

      linkage()
      {
          mma = new refresh::memory_monotonic_unsafe(16 << 20, std::max<size_t>(16, alignof(element)));
      }

      ~linkage()
      {
          delete mma;
      }
      ////////////////
   };

   template<typename T>
   struct my_max {
       T operator()(const T &x, const T &y) const
       {
           return std::max<T>(x, y);
       }
   };
   
   template<typename T>
   struct my_min {
       T operator()(const T &x, const T &y) const
       {
           return std::min<T>(x, y);
       }
   };
   

   template <class Distance, typename AggregationRule = my_max<double>>
      class complete_linkage : public linkage<Distance, AggregationRule>
   {
   };

   template <class Distance>
      class single_linkage : public linkage<Distance, my_min<double>> 
   {
   };
}


inline std::ostream & operator << (std::ostream & sos, const std::vector<std::tuple<std::size_t, std::size_t, double>> & dist)
{
   for (const auto& [w, k, d] : dist)
      sos << "[" << w << ", " << k << ", " << d << "] ";
   sos << std::endl;
   return sos;
}

template <class Distance>
class SingleLinkage : public linkage_algorithm_heaptrix::single_linkage<Distance> 
{
   int operator()
      (
       SparseMatrix<Distance>& distances,
       const std::vector<int>& objects,
       double threshold,
       std::vector<int>& assignments
      ) override 
      {
   //      std::cout << "started" << std::endl;
         return this->run (distances, objects, threshold, assignments);
      }
};

template <class Distance>
class CompleteLinkage : public linkage_algorithm_heaptrix::complete_linkage<Distance> 
{
   int operator()
      (
       SparseMatrix<Distance>& distances,
       const std::vector<int>& objects,
       double threshold,
       std::vector<int>& assignments
      ) override 
      {
     //    std::cout << "started" << std::endl;
         return this->run (distances, objects, threshold, assignments);
      }
};
