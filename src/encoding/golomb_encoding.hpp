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
inline void getDeltaEncoding(InputIterator begin, InputIterator end, BlockWriterIt<OutputIterator>& blockWriter, size_t b) {
  using BlockWriter = BlockWriterIt<OutputIterator>;
  using GolombWriter = thrill::core::GolombBitStreamWriter<BlockWriter>;
  using DeltaStreamWriter = thrill::core::DeltaStreamWriter<GolombWriter, size_t>;

  GolombWriter golombWriter(blockWriter, b);
  DeltaStreamWriter deltaStreamWriter(golombWriter);
  InputIterator it = begin;
  while (it != end) {
    deltaStreamWriter.Put(*it);
    ++it;
  }
}


template <typename InputIterator, typename OutputIterator>
inline void getDeltaEncoding(InputIterator begin, InputIterator end, OutputIterator out, size_t b) {
  using StreamWriter = thrill::core::GolombBitStreamWriter<BlockWriter>;
  BlockWriterIt<OutputIterator> blockWriter(out);
  getDeltaEncoding(begin, end, blockWriter, b);
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

template <typename InputIterator, typename OutputIterator>
void getDeltaDecoding(const InputIterator begin, const InputIterator end, OutputIterator out, size_t b) {
  using GolombReader = thrill::core::GolombBitStreamReader<BlockReader>;
  using DeltaReader = thrill::core::DeltaStreamReader<GolombReader, size_t>;

  BlockReader reader(begin, end);
  GolombReader golombReader(reader, b);
  DeltaReader deltaReader(golombReader);


  while(golombReader.HasNext()) {
    out = deltaReader.Next<size_t>();
  }
}
