//
// This file is part of the libWetCloth open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing
// officer.
//
// We would like to hear from you if you appreciate this work.
//
// The Clear BSD License
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors
// may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
// THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef THREAD_UTILS
#define THREAD_UTILS

#include <vector>
template <typename T, typename A>
using vsize_type = typename std::vector<T, A>::size_type;

#define DEBUG_PARALLEL 1
// #define NO_PARALLEL 1

#ifndef NO_PARALLEL
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
// #include <tbb/parallel_for_each.h> // NOTE could use this instead of manual
// for_each
// #include <tbb/parallel_reduce.h>
#include <thread>
#include <time.h>
#include <iostream>
#endif

namespace threadutils {

#ifdef NO_PARALLEL

inline unsigned get_num_threads() { return 1U; }

template <typename Index, typename Callable>
static void parallel_for(Index start, Index end, Callable func) {
  for (Index i = start; i < end; ++i) {
    func(i);
  }
}

template <typename Data, typename A, typename Callable>
static void parallel_for_each(std::vector<Data, A> &vec, Callable func) {
  for (auto &data : vec) {
    func(data);
  }
}

template <typename Data, typename A, typename Callable>
static void parallel_for_each(const std::vector<Data, A> &vec, Callable func) {
  for (const auto &data : vec) {
    func(data);
  }
}

#else  // PARALLEL =====================================

inline unsigned get_num_threads() {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
  return std::thread::hardware_concurrency();
#else
  return 1U;
#endif
}

template <typename Index, typename Callable>
static void parallel_for(Index start, Index end, Callable func) {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
  tbb::parallel_for(start, end, static_cast<Index>(1), func);
#else
  for (Index i = start; i < end; ++i) {
    func(i);
  }
#endif
}

template <typename Data, typename A, typename Callable>
static void parallel_for_each(std::vector<Data, A> &vec, Callable func) {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
  parallel_for(static_cast<vsize_type<Data, A> >(0), vec.size(),
               [&](vsize_type<Data, A> i) { func(vec[i]); });
#else
  for (auto &data : vec) {
    func(data);
  }
#endif
}

template <typename Data, typename A, typename Callable>
static void parallel_for_each(const std::vector<Data, A> &vec, Callable func) {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
  parallel_for(static_cast<vsize_type<Data, A> >(0), vec.size(),
               [&](vsize_type<Data, A> i) { func(vec[i]); });
#else
  for (const auto &data : vec) {
    func(data);
  }
#endif
}

#endif
}     // namespace threadutils
#endif /* THREAD_UTILS */
