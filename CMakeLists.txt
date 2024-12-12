cmake_minimum_required(VERSION 3.27)
project(heptagon_24autumn_recruitment CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

include_directories(include)

file(GLOB SOURCE_FILES "src/*.cc" "main.cc")

add_executable (SmithWaterman ${SOURCE_FILES})

enable_testing()

set(dataset ${CMAKE_SOURCE_DIR}/data)

add_test(NAME "tiny case" COMMAND ${CMAKE_BINARY_DIR}/SmithWaterman ${dataset}/tiny_q.fna ${dataset}/tiny_t.fna ${dataset}/tiny_ref.txt)
add_test(NAME "small case" COMMAND ${CMAKE_BINARY_DIR}/SmithWaterman ${dataset}/small_q.fna ${dataset}/small_t.fna ${dataset}/small_ref.txt)
add_test(NAME "medium case" COMMAND ${CMAKE_BINARY_DIR}/SmithWaterman ${dataset}/medium_q.fna ${dataset}/medium_t.fna ${dataset}/medium_ref.txt)

if (CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
  target_compile_options(SmithWaterman PRIVATE -O3 -acc -Minfo=acc)
  find_package(OpenACC) # 若成功，执行cmake -B 会输出 "-- Found OpenACC_CXX: -acc"  
  if(OpenACC_CXX_FOUND)
    target_link_libraries(SmithWaterman OpenACC::OpenACC_CXX) # 否则无法正确链接acc的一系列库
  endif()
endif()