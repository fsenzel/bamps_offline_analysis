
/** @file 
 *  Providing type_index based on C++0x proposition from gcc 4.6
 */

#ifndef TYPEINDEX_H
#define TYPEINDEX_H

#include <typeinfo>
#include <cstddef>
#include <boost/functional/hash.hpp>

using std::type_info;
using std::size_t;

/** @brief The class type_index provides a simple wrapper for type_info which can be used as an index type in associative containers (23.6) and in unordered associative containers (23.7). */
struct type_index
{
  type_index(const type_info& __rhs)
  : _M_target(&__rhs) { }
  
  bool
  operator==(const type_index& __rhs) const
  { return *_M_target == *__rhs._M_target; }
  
  bool
  operator!=(const type_index& __rhs) const
  { return *_M_target != *__rhs._M_target; }
  
  bool
  operator<(const type_index& __rhs) const
  { return _M_target->before(*__rhs._M_target); }
  
  bool
  operator<=(const type_index& __rhs) const
  { return !__rhs._M_target->before(*_M_target); }
  
  bool
  operator>(const type_index& __rhs) const
  { return __rhs._M_target->before(*_M_target); }
  
  bool
  operator>=(const type_index& __rhs) const
  { return !_M_target->before(*__rhs._M_target); }
  
  size_t hash_code() const
  {
    boost::hash<std::string> string_hash;
    return string_hash( _M_target->name() );
  }
  
  const char*
  name() const
  { return _M_target->name(); }
  
  private:
    const type_info* _M_target;
};

template<typename _Tp> struct hash;

/// std::hash specialization for type_index.
template<>
struct hash<type_index>
{
  typedef size_t        result_type;
  typedef type_index  argument_type;
  
  size_t
  operator()(const type_index& __ti) const
  { return __ti.hash_code(); }
};

#endif 


// The coypright note from the original gcc version

// C++0x typeindex -*- C++ -*-

// Copyright (C) 2010 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.