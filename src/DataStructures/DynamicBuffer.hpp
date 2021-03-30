// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <vector>

#include "DataStructures/DataVector.hpp"

/**
 * \ingroup DataStructuresGroup
 *
 * \brief A dynamically sized vector of DataVectors.
 *
 * \details This class is useful when one wants to create a
 * `std::vector<DataVector>` where the size of the vector is not known at
 * compile-time. It allocates all `DataVector`s in a single memory chunk rather
 * than allocating each `DataVector` individually. If the size of the vector is
 * known at compile time, a `TempBuffer` object should be used instead.
 *
 * If needed in the future, this class can be templated on the two containers.
 *
 */
class DynamicBuffer {
 public:
  DynamicBuffer() = default;
  DynamicBuffer(size_t outer_size, size_t inner_size);
  ~DynamicBuffer() = default;
  DynamicBuffer(const DynamicBuffer& other) noexcept;
  DynamicBuffer(DynamicBuffer&& other) noexcept;
  DynamicBuffer& operator=(DynamicBuffer other) noexcept;
  // copy assignment can also move assign, explicit deletion silences warning
  DynamicBuffer& operator=(DynamicBuffer&& other) = delete;

  inline DataVector& operator[](size_t index) { return data_vectors_[index]; }
  inline DataVector& at(size_t index) { return data_vectors_.at(index); }
  inline const DataVector& operator[](size_t index) const {
    return data_vectors_[index];
  }
  inline const DataVector& at(size_t index) const {
    return data_vectors_.at(index);
  }

  inline std::vector<DataVector>::iterator begin() {
    return data_vectors_.begin();
  }
  inline std::vector<DataVector>::iterator end() { return data_vectors_.end(); }
  inline const std::vector<DataVector>::const_iterator begin() const {
    return data_vectors_.begin();
  }
  inline const std::vector<DataVector>::const_iterator end() const {
    return data_vectors_.end();
  }

  size_t size() const noexcept { return data_vectors_.size(); }

  friend void swap(DynamicBuffer& first, DynamicBuffer& second) {
    using std::swap;
    swap(first.inner_size_, second.inner_size_);
    swap(first.buffer_, second.buffer_);
    swap(first.data_vectors_, second.data_vectors_);
  }

 private:
  size_t inner_size_;
  // vector of non-owning DataVectors pointing into `buffer_`
  std::vector<DataVector> data_vectors_;
  // memory buffer for all `DataVector`s
  std::vector<double> buffer_;

  // sets data references for all `data_vectors_` into `buffer_`
  void set_references_();
};
