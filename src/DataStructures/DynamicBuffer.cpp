// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/DynamicBuffer.hpp"

DynamicBuffer::DynamicBuffer(size_t outer_size, size_t inner_size)
    : inner_size_(inner_size),
      data_vectors_(outer_size),
      buffer_(outer_size * inner_size) {
  set_references_();
}

DynamicBuffer::DynamicBuffer(const DynamicBuffer& other) noexcept
    : inner_size_(other.inner_size_),
      data_vectors_(other.size()),
      buffer_(other.buffer_) {
  set_references_();
}

DynamicBuffer::DynamicBuffer(DynamicBuffer&& other) noexcept : DynamicBuffer() {
  swap(*this, other);
}

DynamicBuffer& DynamicBuffer::operator=(DynamicBuffer other) noexcept {
  swap(*this, other);
  return *this;
}

void DynamicBuffer::set_references_() {
  for (size_t i = 0; i < size(); ++i) {
    data_vectors_[i].set_data_ref(&buffer_[inner_size_ * i], inner_size_);
  }
}
