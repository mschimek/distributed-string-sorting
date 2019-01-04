/*******************************************************************************
 * src/sequential/bingmann-lcp_losertree.hpp
 *
 * Implementation of a LCP aware multiway losertree.
 *
 *******************************************************************************
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once

#include <iostream>
#include <algorithm>
#include <utility>

#include "merge/stringtools.hpp"
#include "strings/stringptr.hpp"
#include <tlx/define/likely.hpp>

namespace bingmann {

using namespace stringtools;

typedef unsigned char* string;

/******************************************************************************/
// LcpStringLoserTree

template <size_t K>
class LcpStringLoserTree
{
    typedef LcpStringPtr Stream;

    struct Node
    {
        size_t idx;
        lcp_t  lcp;
    };

private:
    Stream streams[K + 1];
    Node nodes[K + 1];

    //! play one comparison edge game: contender is the node below
    //! defender. After the game, defender contains the lower index, contender
    //! the winning index, and defender.lcp = lcp(s_loser,s_winner).
    void updateNode(Node& contender, Node& defender)
    {
        const Stream& defenderStream = streams[defender.idx];

        if (TLX_UNLIKELY(defenderStream.empty()))
            return;

        const Stream& contenderStream = streams[contender.idx];

        if (TLX_UNLIKELY(contenderStream.empty()))
        {
            std::swap(defender, contender);
            return;
        }
#if 1
        if (defender.lcp > contender.lcp)
        {
            // CASE 2: curr->lcp > contender->lcp => curr < contender
            std::swap(defender, contender);
        }
        else if (defender.lcp == contender.lcp)
        {
            // CASE 1: compare more characters
            lcp_t lcp = defender.lcp;

            string s1 = defenderStream.firstString() + lcp;
            string s2 = contenderStream.firstString() + lcp;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != 0 && *s1 == *s2)
                s1++, s2++, lcp++;

            if (*s1 < *s2) // CASE 1.1: curr < contender
                std::swap(defender, contender);

            // update inner node with lcp(s_1,s_2)
            defender.lcp = lcp;
        }
        else {
            // CASE 3: curr->lcp < contender->lcp => contender < curr  => nothing to do
        }
#else
        lcp_compare(contender.idx, contenderStream.firstString(), contender.lcp,
                    defender.idx, defenderStream.firstString(), defender.lcp,
                    contender.idx, contender.lcp, defender.idx, defender.lcp);
#endif
        assert(scmp(streams[contender.idx].firstString(),
                    streams[defender.idx].firstString()) <= 0);

        assert(calc_lcp(streams[contender.idx].firstString(),
                        streams[defender.idx].firstString()) == defender.lcp);
    }

    void initTree(lcp_t knownCommonLcp)
    {
        //std::cout << "inittree start\n";
        for (size_t k = 1; k <= K; k++)
        {
            Node contender;
            contender.idx = k;
            contender.lcp = knownCommonLcp;

            size_t nodeIdx = K + k;

            //std::cout << "nodeIdx " << nodeIdx << "\n";

            while (nodeIdx % 2 == 0 && nodeIdx > 2)
            {
                nodeIdx >>= 1;
                //std::cout << "play against " << nodeIdx << "\n";
                updateNode(contender, nodes[nodeIdx]);
            }
            nodeIdx = (nodeIdx + 1) / 2;
            //std::cout << "save as " << nodeIdx << "\n";
            nodes[nodeIdx] = contender;
        }
        //std::cout << "inittree done\n";
    }

public:
    LcpStringLoserTree(const LcpStringPtr& input, const std::pair<size_t, size_t>* ranges,
                       lcp_t knownCommonLcp = 0)
    {
        for (size_t i = 1; i <= K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i - 1];

            streams[i] = input.sub(currRange.first, currRange.second);
        }

        initTree(knownCommonLcp);
    }



    void writeElementsToStream(LcpStringPtr outStream, const size_t length)
    {
        const LcpStringPtr end = outStream.sub(length, 0);
        while (outStream < end)
        {
            // take winner and put into output

            size_t winnerIdx = nodes[1].idx;
            //std::cout << "winnerIdx " << winnerIdx << std::endl;
            outStream.setFirst(streams[winnerIdx].firstString(), nodes[1].lcp);
            ++outStream;

            // advance winner stream

            Stream& stream = streams[winnerIdx];
            ++stream;

            // run new items from winner stream up the tree

            Node& contender = nodes[1];

            if (!stream.empty())
                contender.lcp = streams[winnerIdx].firstLcp();

            size_t nodeIdx = winnerIdx + K;
            //std::cout << "nodeIdx " << nodeIdx << "\n";

            while (nodeIdx > 2) {
                nodeIdx = (nodeIdx + 1) / 2;
                //std::cout << "play against " << nodeIdx << "\n";
                updateNode(contender, nodes[nodeIdx]);
            }
            //std::cout << "play against " << nodeIdx << "\n";

            // for (size_t nodeIdx = (K + winnerIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
            // {
            //     updateNode(contender, nodes[nodeIdx]);
            // }
        }
    }
};




} // namespace bingmann

namespace dss_schimek {
  template <size_t K, typename StringSet>
class LcpStringLoserTree_
{
    typedef dss_schimek::StringLcpPtrMergeAdapter<StringSet> Stream;
    using CharIt = typename StringSet::CharIterator;

    struct Node
    {
        size_t idx;
        lcp_t  lcp;
    };

private:
    Stream streams[K + 1];
    Node nodes[K + 1];

    //! play one comparison edge game: contender is the node below
    //! defender. After the game, defender contains the lower index, contender
    //! the winning index, and defender.lcp = lcp(s_loser,s_winner).
    void updateNode(Node& contender, Node& defender)
    {
        //std::cout << "\t\t\t contender.idx: " << contender.idx  << " contender.lcp:" << contender.lcp  << " contender: " << streams[contender.idx].firstStringChars() << " defender.idx:" << defender.idx << " defender.lcp:" << defender.lcp <<  " defender " << streams[defender.idx].firstStringChars() << std::endl;
        const Stream& defenderStream = streams[defender.idx];

        if (TLX_UNLIKELY(defenderStream.empty()))
            return;

        const Stream& contenderStream = streams[contender.idx];

        if (TLX_UNLIKELY(contenderStream.empty()))
        {
            std::swap(defender, contender);
            return;
        }
#if 1
        if (defender.lcp > contender.lcp)
        {
            // CASE 2: curr->lcp > contender->lcp => curr < contender
            std::swap(defender, contender);
        }
        else if (defender.lcp == contender.lcp)
        {
            // CASE 1: compare more characters
            lcp_t lcp = defender.lcp;

            CharIt s1 = defenderStream.firstStringChars() + lcp;
            CharIt s2 = contenderStream.firstStringChars() + lcp;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != 0 && *s1 == *s2)
                s1++, s2++, lcp++;

            if (*s1 < *s2) // CASE 1.1: curr < contender
                std::swap(defender, contender);

            // update inner node with lcp(s_1,s_2)
            defender.lcp = lcp;
        }
        else {
            // CASE 3: curr->lcp < contender->lcp => contender < curr  => nothing to do
        }
#else
        lcp_compare(contender.idx, contenderStream.firstString(), contender.lcp,
                    defender.idx, defenderStream.firstString(), defender.lcp,
                    contender.idx, contender.lcp, defender.idx, defender.lcp);
#endif
        assert(scmp(streams[contender.idx].firstStringChars(),
                    streams[defender.idx].firstStringChars()) <= 0);

        assert(calc_lcp(streams[contender.idx].firstStringChars(),
                        streams[defender.idx].firstStringChars()) == defender.lcp);
    }

    void initTree(lcp_t knownCommonLcp)
    {
        //std::cout << "inittree start\n";
        for (size_t k = 1; k <= K; k++)
        {
            Node contender;
            contender.idx = k;
            contender.lcp = knownCommonLcp;

            size_t nodeIdx = K + k;

            //std::cout << "nodeIdx " << nodeIdx << "\n";

            while (nodeIdx % 2 == 0 && nodeIdx > 2)
            {
                nodeIdx >>= 1;
                //std::cout << "play against " << nodeIdx << "\n";
                updateNode(contender, nodes[nodeIdx]);
            }
            nodeIdx = (nodeIdx + 1) / 2;
            //std::cout << "save as " << nodeIdx << "\n";
            nodes[nodeIdx] = contender;
        }
        //std::cout << "inittree done\n";
    }

public:
    LcpStringLoserTree_(const dss_schimek::StringLcpPtrMergeAdapter<StringSet>& input, const std::pair<size_t, size_t>* ranges,
                       lcp_t knownCommonLcp = 0)
    {
        for (size_t i = 1; i <= K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i - 1];

            streams[i] = input.sub(currRange.first, currRange.second);
        }

        initTree(knownCommonLcp);
    }



    void writeElementsToStream(dss_schimek::StringLcpPtrMergeAdapter<StringSet> outStream, const size_t length)
    {
        const dss_schimek::StringLcpPtrMergeAdapter<StringSet> end = outStream.sub(length, 0);
        size_t counter = 0;
        for (size_t i = 1; i < K + 1; ++i) {
                std::cout << "\t\t\t" << i << " idx:" << nodes[i].idx << " lcp:" << nodes[i].lcp << std::endl;
              }
        std::cout << "length: " << length << std::endl;
        while (outStream < end)
        {
            // take winner and put into output

            size_t winnerIdx = nodes[1].idx;
            std::cout << "winnerIdx " << winnerIdx << std::endl;
            outStream.setFirst(streams[winnerIdx].firstString(), nodes[1].lcp);
            ++outStream;
            std::cout << "counter " << counter << " streams[winnerIdx]" << streams[winnerIdx].firstStringChars() << std::endl;

            counter++;
            // advance winner stream

            Stream& stream = streams[winnerIdx];
            ++stream;
            // run new items from winner stream up the tree

            Node& contender = nodes[1];

            if (!stream.empty()) {
              contender.lcp = streams[winnerIdx].firstLcp();
              std::cout << "winnerIdx: " << winnerIdx << " contender.lcp: " << contender.lcp << std::endl;

            }
            else {
              std::cout << "empty" << std::endl;
            }

            size_t nodeIdx = winnerIdx + K;
            //std::cout << "nodeIdx " << nodeIdx << "\n";

              std::cout << "\t\t\trecalculate tree" << std::endl;
              while (nodeIdx > 2) {
                nodeIdx = (nodeIdx + 1) / 2;
                //std::cout << "play against " << nodeIdx << "\n";
                updateNode(contender, nodes[nodeIdx]);
            }
              for (size_t i = 1; i < K + 1; ++i) {
                std::cout << "\t\t\t" << i << " idx:" << nodes[i].idx << " lcp:" << nodes[i].lcp << std::endl;
              }
            //std::cout << "play against " << nodeIdx << "\n";

            // for (size_t nodeIdx = (K + winnerIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
            // {
            //     updateNode(contender, nodes[nodeIdx]);
            // }
        }
    }
};
}

/******************************************************************************/
