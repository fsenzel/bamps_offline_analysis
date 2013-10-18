#=============================================================================
# Find properties of the host system:
# * SSE/AVX support
# * RDTSCP opcode support (needed for TSC timing)
#
#
# Parts are taken from "Vc" by Matthias Kretz:
# (copyright see below)
# 1) The macros "_my_find", "AutodetectHostArchitecture" 
#    are copied one-to-one (and enhanced by some
#    additional hosts)
# 2) "SetHostVectorSupport" is a boiled down version of
#    "OptimizeForArchitecture"
#=============================================================================


# Taken from "Vc":
#=============================================================================
# Copyright 2010-2013 Matthias Kretz <kretz@kde.org>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * The names of Kitware, Inc., the Insight Consortium, or the names of
#    any consortium members, or of any contributors, may not be used to
#    endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

macro(_my_find _list _value _ret)
   list(FIND ${_list} "${_value}" _found)
   if(_found EQUAL -1)
      set(${_ret} FALSE)
   else(_found EQUAL -1)
      set(${_ret} TRUE)
   endif(_found EQUAL -1)
endmacro(_my_find)

macro(AutodetectHostArchitecture)
   set(TARGET_ARCHITECTURE "generic")
   set(_vendor_id)
   set(_cpu_family)
   set(_cpu_model)
   if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      file(READ "/proc/cpuinfo" _cpuinfo)
      string(REGEX REPLACE ".*vendor_id[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _vendor_id "${_cpuinfo}")
      string(REGEX REPLACE ".*cpu family[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_family "${_cpuinfo}")
      string(REGEX REPLACE ".*model[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_model "${_cpuinfo}")
      string(REGEX REPLACE ".*flags[ \t]*:[ \t]+([^\n]+).*" "\\1" _cpu_flags "${_cpuinfo}")
   elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      exec_program("/usr/sbin/sysctl -n machdep.cpu.vendor" OUTPUT_VARIABLE _vendor_id)
      exec_program("/usr/sbin/sysctl -n machdep.cpu.model"  OUTPUT_VARIABLE _cpu_model)
      exec_program("/usr/sbin/sysctl -n machdep.cpu.family" OUTPUT_VARIABLE _cpu_family)
      exec_program("/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE _cpu_flags)
      string(TOLOWER "${_cpu_flags}" _cpu_flags)
      string(REPLACE "." "_" _cpu_flags "${_cpu_flags}")
   elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      get_filename_component(_vendor_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;VendorIdentifier]" NAME CACHE)
      get_filename_component(_cpu_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;Identifier]" NAME CACHE)
      mark_as_advanced(_vendor_id _cpu_id)
      string(REGEX REPLACE ".* Family ([0-9]+) .*" "\\1" _cpu_family "${_cpu_id}")
      string(REGEX REPLACE ".* Model ([0-9]+) .*" "\\1" _cpu_model "${_cpu_id}")
   endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   if(_vendor_id STREQUAL "GenuineIntel")
      if(_cpu_family EQUAL 6)
         # Any recent Intel CPU except NetBurst
         if(_cpu_model EQUAL 58)
            set(TARGET_ARCHITECTURE "ivy-bridge")
         elseif(_cpu_model EQUAL 46) # Xeon 7500 series
            set(TARGET_ARCHITECTURE "westmere")
         elseif(_cpu_model EQUAL 45) # Xeon TNG
            set(TARGET_ARCHITECTURE "sandy-bridge")
         elseif(_cpu_model EQUAL 44) # Xeon 5600 series
            set(TARGET_ARCHITECTURE "westmere")
         elseif(_cpu_model EQUAL 42) # Core TNG
            set(TARGET_ARCHITECTURE "sandy-bridge")
         elseif(_cpu_model EQUAL 37) # Core i7/i5/i3
            set(TARGET_ARCHITECTURE "westmere")
         elseif(_cpu_model EQUAL 31) # Core i7/i5
            set(TARGET_ARCHITECTURE "westmere")
         elseif(_cpu_model EQUAL 30) # Core i7/i5
            set(TARGET_ARCHITECTURE "westmere")
         elseif(_cpu_model EQUAL 29)
            set(TARGET_ARCHITECTURE "penryn")
         elseif(_cpu_model EQUAL 28)
            set(TARGET_ARCHITECTURE "atom")
         elseif(_cpu_model EQUAL 26)
            set(TARGET_ARCHITECTURE "nehalem")
         elseif(_cpu_model EQUAL 23)
            set(TARGET_ARCHITECTURE "penryn")
         elseif(_cpu_model EQUAL 15)
            set(TARGET_ARCHITECTURE "merom")
         elseif(_cpu_model EQUAL 14)
            set(TARGET_ARCHITECTURE "core")
         elseif(_cpu_model LESS 14)
            message(WARNING "Your CPU (family ${_cpu_family}, model ${_cpu_model}) is not known. Auto-detection of optimization flags failed and will use the generic CPU settings with SSE2.")
            set(TARGET_ARCHITECTURE "generic")
         else()
            message(WARNING "Your CPU (family ${_cpu_family}, model ${_cpu_model}) is not known. Auto-detection of optimization flags failed and will use the 65nm Core 2 CPU settings.")
            set(TARGET_ARCHITECTURE "merom")
         endif()
      elseif(_cpu_family EQUAL 7) # Itanium (not supported)
         message(WARNING "Your CPU (Itanium: family ${_cpu_family}, model ${_cpu_model}) is not supported by OptimizeForArchitecture.cmake.")
      elseif(_cpu_family EQUAL 15) # NetBurst
         list(APPEND _available_vector_units_list "sse" "sse2")
         if(_cpu_model GREATER 2) # Not sure whether this must be 3 or even 4 instead
            list(APPEND _available_vector_units_list "sse" "sse2" "sse3")
         endif(_cpu_model GREATER 2)
      endif(_cpu_family EQUAL 6)
   elseif(_vendor_id STREQUAL "AuthenticAMD")
      if(_cpu_family EQUAL 21) # 15h
         set(TARGET_ARCHITECTURE "bulldozer")
      elseif(_cpu_family EQUAL 20) # 14h
      elseif(_cpu_family EQUAL 18) # 12h
         if(_cpu_model EQUAL 1) # Llano
	    set(TARGET_ARCHITECTURE "llano")
         endif()
      elseif(_cpu_family EQUAL 16) # 10h
         set(TARGET_ARCHITECTURE "barcelona")
      elseif(_cpu_family EQUAL 15)
         set(TARGET_ARCHITECTURE "k8")
         if(_cpu_model GREATER 64) # I don't know the right number to put here. This is just a guess from the hardware I have access to
            set(TARGET_ARCHITECTURE "k8-sse3")
         endif(_cpu_model GREATER 64)
      endif()
   endif(_vendor_id STREQUAL "GenuineIntel")
endmacro()

macro(SetHostVectorSupport)
   set(TARGET_ARCHITECTURE "auto" CACHE STRING "CPU architecture to optimize for. Using an incorrect setting here can result in crashes of the resulting binary because of invalid instructions used.\nSetting the value to \"auto\" will try to optimize for the architecture where cmake is called.\nOther supported values are: \"none\", \"generic\", \"core\", \"merom\" (65nm Core2), \"penryn\" (45nm Core2), \"nehalem\", \"westmere\", \"sandy-bridge\", \"ivy-bridge\", \"atom\", \"k8\", \"k8-sse3\", \"barcelona\", \"istanbul\", \"magny-cours\", \"bulldozer\", \"llano\", \"interlagos\".")
   set(_force)
   if(NOT _last_target_arch STREQUAL "${TARGET_ARCHITECTURE}")
#      message(STATUS "target changed from \"${_last_target_arch}\" to \"${TARGET_ARCHITECTURE}\"")
      set(_force FORCE)
   endif()
   set(_last_target_arch "${TARGET_ARCHITECTURE}" CACHE STRING "" FORCE)
   mark_as_advanced(_last_target_arch)
   string(TOLOWER "${TARGET_ARCHITECTURE}" TARGET_ARCHITECTURE)

   set(_available_vector_units_list)

   if(TARGET_ARCHITECTURE STREQUAL "auto")
      AutodetectHostArchitecture()
      message(STATUS "Detected CPU: ${TARGET_ARCHITECTURE}")
   endif(TARGET_ARCHITECTURE STREQUAL "auto")

   if(TARGET_ARCHITECTURE STREQUAL "core")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "merom")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "penryn")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3")
      message(STATUS "Sadly the Penryn architecture exists in variants with SSE4.1 and without SSE4.1.")
      if(_cpu_flags MATCHES "sse4_1")
         message(STATUS "SSE4.1: enabled (auto-detected from this computer's CPU flags)")
         list(APPEND _available_vector_units_list "sse4.1")
      else()
         message(STATUS "SSE4.1: disabled (auto-detected from this computer's CPU flags)")
      endif()
   elseif(TARGET_ARCHITECTURE STREQUAL "nehalem")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2")
   elseif(TARGET_ARCHITECTURE STREQUAL "westmere")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2")
   elseif(TARGET_ARCHITECTURE STREQUAL "ivy-bridge")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2" "avx" "rdrnd" "f16c")
   elseif(TARGET_ARCHITECTURE STREQUAL "sandy-bridge")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2" "avx")
   elseif(TARGET_ARCHITECTURE STREQUAL "atom")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "k8")
      list(APPEND _available_vector_units_list "sse" "sse2")
   elseif(TARGET_ARCHITECTURE STREQUAL "k8-sse3")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "interlagos")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4a" "sse4.1" "sse4.2" "avx" "xop" "fma4")
   elseif(TARGET_ARCHITECTURE STREQUAL "bulldozer")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4a" "sse4.1" "sse4.2" "avx" "xop" "fma4")
   elseif(TARGET_ARCHITECTURE STREQUAL "barcelona")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "istanbul")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "magny-cours")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "llano")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "generic")
      list(APPEND _march_flag_list "generic")
   elseif(TARGET_ARCHITECTURE STREQUAL "none")
      # add this clause to remove it from the else clause
   else(TARGET_ARCHITECTURE STREQUAL "core")
      message(FATAL_ERROR "Unknown target architecture: \"${TARGET_ARCHITECTURE}\". Please set TARGET_ARCHITECTURE to a supported value.")
   endif(TARGET_ARCHITECTURE STREQUAL "core")

   if(NOT TARGET_ARCHITECTURE STREQUAL "none")
      _my_find(_available_vector_units_list "sse2" SSE2_FOUND)
      _my_find(_available_vector_units_list "sse3" SSE3_FOUND)
      _my_find(_available_vector_units_list "ssse3" SSSE3_FOUND)
      _my_find(_available_vector_units_list "sse4.1" SSE4_1_FOUND)
      _my_find(_available_vector_units_list "sse4.2" SSE4_2_FOUND)
      _my_find(_available_vector_units_list "sse4a" SSE4a_FOUND)
      if(DEFINED Vc_AVX_INTRINSICS_BROKEN AND Vc_AVX_INTRINSICS_BROKEN)
         UserWarning("AVX disabled per default because of old/broken compiler")
         set(AVX_FOUND false)
         set(XOP_FOUND false)
         set(FMA4_FOUND false)
      else()
         _my_find(_available_vector_units_list "avx" AVX_FOUND)
         if(DEFINED Vc_FMA4_INTRINSICS_BROKEN AND Vc_FMA4_INTRINSICS_BROKEN)
            UserWarning("FMA4 disabled per default because of old/broken compiler")
            set(FMA4_FOUND false)
         else()
            _my_find(_available_vector_units_list "fma4" FMA4_FOUND)
         endif()
         if(DEFINED Vc_XOP_INTRINSICS_BROKEN AND Vc_XOP_INTRINSICS_BROKEN)
            UserWarning("XOP disabled per default because of old/broken compiler")
            set(XOP_FOUND false)
         else()
            _my_find(_available_vector_units_list "xop" XOP_FOUND)
         endif()
      endif()

      set(HAVE_VECTOR_SSE2   ${SSE2_FOUND}   CACHE BOOL "We have support of SSE2. If SSE2 instructions are not enabled the SSE implementation will be disabled." ${_force})
      set(HAVE_VECTOR_SSE3   ${SSE3_FOUND}   CACHE BOOL "We have support of SSE3. If SSE3 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSSE3  ${SSSE3_FOUND}  CACHE BOOL "We have support of SSSE3. If SSSE3 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSE4_1 ${SSE4_1_FOUND} CACHE BOOL "We have support of SSE4.1. If SSE4.1 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSE4_2 ${SSE4_2_FOUND} CACHE BOOL "We have support of SSE4.2. If SSE4.2 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSE4a  ${SSE4a_FOUND}  CACHE BOOL "We have support of SSE4a. If SSE4a instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_AVX    ${AVX_FOUND}    CACHE BOOL "We have support of AVX. This will double some of the vector sizes relative to SSE." ${_force})
      set(HAVE_VECTOR_XOP    ${XOP_FOUND}    CACHE BOOL "We have support of XOP." ${_force})
      set(HAVE_VECTOR_FMA4   ${FMA4_FOUND}   CACHE BOOL "We have support of FMA4." ${_force})

  endif ()

endmacro()

macro(SetHostRDTSCP)
  execute_process(COMMAND cat /proc/cpuinfo
    COMMAND grep flags
    COMMAND grep -c rdtscp
    OUTPUT_VARIABLE HAVE_RDTSCP)

  set(HAVE_RDTSCP ${HAVE_RDTSCP} CACHE INTEGER "CPU has RDTSCP opcode (indicates also number of cores)" FORCE)
endmacro()

### 
# now DO the macros:
###

SetHostVectorSupport()
SetHostRDTSCP()

###
