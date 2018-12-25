#pragma once
#include <vector>
#include <string>
#include <cassert>
#include <cstring>

namespace dss_schimek {

  class EmptyByteEncoder{
    public:
      static std::string getName() {
        return "EmptyByteEncoder";
      }
  };

  class SequentialDelayedByteEncoder {
    /*
     * Encoding: 
     * (|Number Chars : size_t | Number Strings: size_t | [char : unsigned char] | [number : size_t] |)*
     *
     */
    public:
      static std::string getName() {
        return "SequentialDelayedByteEncoder";
      }
      size_t computeNumberOfSendBytes(size_t charsToSend, size_t numbersToSend) const {
        return charsToSend + sizeof(size_t) * numbersToSend + 2 * sizeof(size_t);
      }

      size_t computeNumberOfSendBytes(const std::vector<size_t>& charsToSend,
          const std::vector<size_t>& numbersToSend) const {
        assert(charsToSend.size() == numbersToSend.size());
        size_t numberOfSendBytes = 0;
        for (size_t i = 0; i < charsToSend.size(); ++i) {
          numberOfSendBytes += computeNumberOfSendBytes(charsToSend[i], numbersToSend[i]); 
        }
        return numberOfSendBytes;
      }

      std::pair<size_t, size_t> computeNumberOfRecvData(const char unsigned* buffer, size_t size) const {
        size_t recvChars = 0, recvNumbers  = 0;
        size_t i = 0; 
        while(i < size) {
          size_t curRecvChars = 0;
          size_t curRecvNumbers = 0;
          memcpy(&curRecvChars, buffer + i, sizeof(size_t));
          recvChars += curRecvChars;
          i += sizeof(size_t);
          memcpy(&curRecvNumbers, buffer + i, sizeof(size_t));
          recvNumbers += curRecvNumbers;
          i += sizeof(size_t);
          i += curRecvChars + sizeof(size_t) * curRecvNumbers;
        }
        return std::make_pair(recvChars, recvNumbers); 
      }


      // think about common interface
      unsigned char*  write(unsigned char* buffer,
          const unsigned char* charsToWrite,
          size_t numChars,
          const size_t* numbersToWrite,
          size_t numNumbers) const {
        const size_t alignmentSizeT = alignof(size_t);
        const size_t sizeOfSizeT = sizeof(size_t);

        memcpy(buffer, &numChars, sizeOfSizeT);
        buffer += sizeOfSizeT;
        memcpy(buffer, &numNumbers, sizeOfSizeT);
        buffer += sizeOfSizeT;
        memcpy(buffer, charsToWrite, numChars);
        buffer += numChars;
        const size_t bytesToWrite = sizeof(size_t) * numNumbers;
        memcpy(buffer, numbersToWrite, bytesToWrite);
        buffer += bytesToWrite;
        return buffer;
      }

      void read_(unsigned char* buffer,
          std::vector<unsigned char>& charsToRead,
          std::vector<size_t>& numbersToRead) const {

        size_t numCharsToRead = 0;
        size_t numNumbersToRead = 0;
        memcpy(&numCharsToRead, buffer, sizeof(size_t));
        buffer += sizeof(size_t);
        memcpy(&numNumbersToRead, buffer, sizeof(size_t));
        buffer += sizeof(size_t);
        charsToRead.clear();
        charsToRead.reserve(numCharsToRead);
        numbersToRead.clear();
        numbersToRead.reserve(numNumbersToRead);

        for (size_t i = 0; i < numCharsToRead; ++i) 
          charsToRead.emplace_back(buffer[i]);
        for (size_t i = 0; i < numNumbersToRead; ++i) {
          size_t tmp;
          std::memcpy(&tmp, buffer + numCharsToRead + i * sizeof(size_t), sizeof(size_t));
          numbersToRead.emplace_back(tmp);
        }  
      }
      std::pair<std::vector<unsigned char>, std::vector<size_t>> read(unsigned char* buffer, size_t size) const {

        auto recvData = computeNumberOfRecvData(buffer, size);
        size_t numCharsToRead = recvData.first;
        size_t numNumbersToRead = recvData.second;
        std::vector<unsigned char> charsToRead;
        std::vector<size_t> numbersToRead;
        charsToRead.reserve(numCharsToRead);
        numbersToRead.reserve(numNumbersToRead);

        size_t curPos = 0;
        while(curPos < size) {
          memcpy(&numCharsToRead, buffer + curPos, sizeof(size_t));
          curPos += sizeof(size_t);
          memcpy(&numNumbersToRead, buffer + curPos, sizeof(size_t));
          curPos += sizeof(size_t);

          for(size_t i = 0; i < numCharsToRead; ++i) 
            charsToRead.emplace_back(*(buffer + curPos + i));
          curPos += numCharsToRead;
          for (size_t i = 0; i < numNumbersToRead; ++i) {
            size_t tmp;
            std::memcpy(&tmp, buffer + curPos + i * sizeof(size_t), sizeof(size_t));
            numbersToRead.emplace_back(tmp);
          }
          curPos += numNumbersToRead * sizeof(size_t);
        }
        return make_pair(std::move(charsToRead), std::move(numbersToRead));
      }
  };

  class SequentialByteEncoder {
    /*
     * Encoding: 
     * (|Number Chars : size_t | Number Strings: size_t | [char : unsigned char] | [number : size_t] |)*
     *
     */
    public:
      static std::string getName() {
        return "SequentialByteEncoder";
      }
      size_t computeNumberOfSendBytes(size_t charsToSend, size_t numbersToSend) const {
        return charsToSend + sizeof(size_t) * numbersToSend + 2 * sizeof(size_t);
      }

      size_t computeNumberOfSendBytes(const std::vector<size_t>& charsToSend,
          const std::vector<size_t>& numbersToSend) const {
        assert(charsToSend.size() == numbersToSend.size());
        size_t numberOfSendBytes = 0;
        for (size_t i = 0; i < charsToSend.size(); ++i) {
          numberOfSendBytes += computeNumberOfSendBytes(charsToSend[i], numbersToSend[i]); 
        }
        return numberOfSendBytes;
      }

      std::pair<size_t, size_t> computeNumberOfRecvData(const char unsigned* buffer, size_t size) const {
        size_t recvChars = 0, recvNumbers  = 0;
        size_t i = 0; 
        while(i < size) {
          size_t curRecvChars = 0;
          size_t curRecvNumbers = 0;
          memcpy(&curRecvChars, buffer + i, sizeof(size_t));
          recvChars += curRecvChars;
          i += sizeof(size_t);
          memcpy(&curRecvNumbers, buffer + i, sizeof(size_t));
          recvNumbers += curRecvNumbers;
          i += sizeof(size_t);
          i += curRecvChars + sizeof(size_t) * curRecvNumbers;
        }
        return std::make_pair(recvChars, recvNumbers); 
      }

      // start work here
      template<typename StringSet>
        unsigned char* write(unsigned char* buffer,
            const StringSet ss,
            const size_t* numbersToWrite) const {

          using String = typename StringSet::String;
          using CharIt = typename StringSet::CharIterator;

          unsigned char* const startOfBuffer = buffer;
          size_t numChars = 0;
          const size_t size = ss.size();

          buffer += sizeof(size_t);
          memcpy(buffer, &size, sizeof(size_t));
          buffer += sizeof(size_t);
          auto beginOfSet = ss.begin();
          for (size_t i = 0; i < size; ++i) {
            String str = ss[beginOfSet + i];
            size_t stringLength = ss.get_length(str) + 1;
            numChars += stringLength;
            memcpy(buffer, ss.get_chars(str, 0), stringLength);
            buffer += stringLength;
          }
          memcpy(startOfBuffer, &numChars, sizeof(size_t));
          memcpy(buffer, numbersToWrite, size * sizeof(size_t));
          buffer += size * sizeof(size_t);

          return buffer;
        }




      std::pair<std::vector<unsigned char>, std::vector<size_t>> read(unsigned char* buffer, size_t size) const {

        auto recvData = computeNumberOfRecvData(buffer, size);
        size_t numCharsToRead = recvData.first;
        size_t numNumbersToRead = recvData.second;
        std::vector<unsigned char> charsToRead;
        std::vector<size_t> numbersToRead;
        charsToRead.reserve(numCharsToRead);
        numbersToRead.reserve(numNumbersToRead);

        size_t curPos = 0;
        while(curPos < size) {
          memcpy(&numCharsToRead, buffer + curPos, sizeof(size_t));
          curPos += sizeof(size_t);
          memcpy(&numNumbersToRead, buffer + curPos, sizeof(size_t));
          curPos += sizeof(size_t);

          for(size_t i = 0; i < numCharsToRead; ++i) 
            charsToRead.emplace_back(*(buffer + curPos + i));
          curPos += numCharsToRead;
          for (size_t i = 0; i < numNumbersToRead; ++i) {
            size_t tmp;
            std::memcpy(&tmp, buffer + curPos + i * sizeof(size_t), sizeof(size_t));
            numbersToRead.emplace_back(tmp);
          }
          curPos += numNumbersToRead * sizeof(size_t);
        }
        return make_pair(std::move(charsToRead), std::move(numbersToRead));
      }
  };

  class InterleavedByteEncoder {
    /*
     * Encoding: 
     * (|Number Chars : size_t | Number Strings : size_t | (char : unsigned char | number : size_t)^NumberStrings )*
     */
    public:
      static std::string getName() {
        return "InterleavedByteEncoder";
      }
      size_t computeNumberOfSendBytes(size_t charsToSend, size_t numbersToSend) const {
        return charsToSend + sizeof(size_t) * numbersToSend + 2 * sizeof(size_t);
      }

      size_t computeNumberOfSendBytes(const std::vector<size_t>& charsToSend,
          const std::vector<size_t>& numbersToSend) const {
        assert(charsToSend.size() == numbersToSend.size());
        size_t numberOfSendBytes = 0;
        for (size_t i = 0; i < charsToSend.size(); ++i) {
          numberOfSendBytes += computeNumberOfSendBytes(charsToSend[i], numbersToSend[i]); 
        }
        return numberOfSendBytes;
      }

      std::pair<size_t, size_t> computeNumberOfRecvData(
          const char unsigned* buffer, size_t size) const {
        size_t recvChars = 0, recvNumbers  = 0;
        size_t i = 0; 
        while(i < size) {
          size_t curRecvChars = 0;
          size_t curRecvNumbers = 0;
          memcpy(&curRecvChars, buffer + i, sizeof(size_t));
          recvChars += curRecvChars;
          i += sizeof(size_t);
          memcpy(&curRecvNumbers, buffer + i, sizeof(size_t));
          recvNumbers += curRecvNumbers;
          i += sizeof(size_t);
          i += curRecvChars + sizeof(size_t) * curRecvNumbers;
        }
        return std::make_pair(recvChars, recvNumbers); 
      }

      template<typename StringSet>
        unsigned char* write(unsigned char* buffer,
            const StringSet ss,
            const size_t* numbersToWrite) const {

          using String = typename StringSet::String;
          using CharIt = typename StringSet::CharIterator;

          unsigned char* const startOfBuffer = buffer;
          size_t numChars = 0;
          const size_t size = ss.size();

          buffer += sizeof(size_t);
          memcpy(buffer, &size, sizeof(size_t));
          buffer += sizeof(size_t);
          auto beginOfSet = ss.begin();
          for (size_t i = 0; i < size; ++i) {
            String str = ss[beginOfSet + i];
            size_t stringLength = ss.get_length(str) + 1;
            numChars += stringLength;
            memcpy(buffer, ss.get_chars(str, 0), stringLength);
            buffer += stringLength;
            memcpy(buffer, numbersToWrite + i, sizeof(size_t));
            buffer += sizeof(size_t);
          }
          memcpy(startOfBuffer, &numChars, sizeof(size_t));

          return buffer;
        }

      std::pair<std::vector<unsigned char>, std::vector<size_t>> read(
          unsigned char* buffer, size_t size) const {

        const auto [numTotalCharsToRead, numTotalNumbersToRead] = 
          computeNumberOfRecvData(buffer, size);
        std::vector<unsigned char> charsToRead;
        std::vector<size_t> numbersToRead;
        charsToRead.reserve(numTotalCharsToRead);
        numbersToRead.reserve(numTotalNumbersToRead);

        size_t curPos = 0;
        while(curPos < size) { // loop over all interval (= pieces received from other PEs)
          size_t numCharsToRead = 0, numNumbersToRead = 0;
          memcpy(&numCharsToRead, buffer + curPos, sizeof(size_t));
          curPos += sizeof(size_t);
          memcpy(&numNumbersToRead, buffer + curPos, sizeof(size_t));
          curPos += sizeof(size_t);

          for(size_t i = 0; i < numNumbersToRead; ++i) {
            while(true) {
              const unsigned char curChar = *(buffer + curPos);
              ++curPos;
              charsToRead.emplace_back(curChar);
              if (curChar == 0)
                break;
            }
            size_t curNumber = 0;
            memcpy(&curNumber, buffer + curPos, sizeof(size_t));
            curPos += sizeof(size_t);
            numbersToRead.emplace_back(curNumber);
          }
        }
        return make_pair(std::move(charsToRead), std::move(numbersToRead));
      }
  };
}
