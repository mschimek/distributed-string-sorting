#pragma once
#include <vector>
#include <iostream>
#include <bitset>

#include "golomb_bit_stream.hpp"
#include "delta_stream.hpp"

struct BlockWriter {
  std::vector<size_t> data;

  void PutRaw(const size_t value) { data.push_back(value); }
  size_t block_size() { return sizeof(size_t); };
};

template<typename T>
struct PrintMe;

template <typename OutputIterator>
struct BlockWriterIt {
  using outItValueType = typename OutputIterator::container_type::value_type;
  using type = typename std::enable_if<std::is_same<outItValueType, size_t>::value, size_t>::type;

  OutputIterator outIt;
  BlockWriterIt(OutputIterator outIt) : outIt(outIt) {}

  void PutRaw(const type value) { outIt = value; }
  size_t block_size() { return sizeof(size_t); };
};

struct BlockReader {
  using iterator = std::vector<size_t>::const_iterator;

  iterator curPos;
  iterator end;
  std::vector<size_t> decoded_data;

  BlockReader(iterator begin, iterator end) :  curPos(begin), end(end) {}
  bool HasNext() { return curPos != end; }
  template<typename T>
  T GetRaw() { return *(curPos++); }

};

template <typename T>
void printBits(const T& value) {
  std::cout << std::bitset<sizeof(T) * 8>(value) << std::endl;
}

template <typename InputIterator, typename OutputIterator>
inline void getDeltaEncoding(InputIterator begin, InputIterator end, BlockWriterIt<OutputIterator>& blockWriter) {
  using BlockWriter = BlockWriterIt<OutputIterator>;
  using GolombWriter = thrill::core::GolombBitStreamWriter<BlockWriter>;
  using DeltaStreamWriter = thrill::core::DeltaStreamWriter<GolombWriter, size_t>;

  GolombWriter golombWriter(blockWriter, 8);
  DeltaStreamWriter deltaStreamWriter(golombWriter);
  InputIterator it = begin;
  while (it != end) {
    deltaStreamWriter.Put(*it);
    std::cout << *it << std::endl;
    ++it;
  }
}


template <typename InputIterator, typename OutputIterator>
inline void getDeltaEncoding(InputIterator begin, InputIterator end, OutputIterator out) {
  using StreamWriter = thrill::core::GolombBitStreamWriter<BlockWriter>;
  BlockWriterIt<OutputIterator> blockWriter(out);
  getDeltaEncoding(begin, end, blockWriter);
}

std::vector<size_t> getDecoding(const std::vector<size_t>& values) {
  using GolombReader = thrill::core::GolombBitStreamReader<BlockReader>;

  BlockReader reader(values.begin(), values.end());
  std::vector<size_t> decodedValues;
  const size_t b = 8;
  GolombReader golombReader(reader, 8);

  while(golombReader.HasNext())
    decodedValues.push_back(golombReader.Next<size_t>());

  return decodedValues;
}
std::vector<size_t> getDeltaDecoding(const std::vector<size_t>& values) {
  using GolombReader = thrill::core::GolombBitStreamReader<BlockReader>;
  using DeltaReader = thrill::core::DeltaStreamReader<GolombReader, size_t>;

  BlockReader reader(values.begin(), values.end());
  std::vector<size_t> decodedValues;
  const size_t b = 8;
  GolombReader golombReader(reader, 8);
  DeltaReader deltaReader(golombReader);


  while(deltaReader.HasNext())
    decodedValues.push_back(deltaReader.Next<size_t>());

  return decodedValues;
}
