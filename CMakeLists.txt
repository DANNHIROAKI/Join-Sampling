cmake_minimum_required(VERSION 3.10)

project(RectSamplerExperiments LANGUAGES CXX)

# 要用到 <filesystem> 可以改成 20，看你需求
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 让编译器多啰嗦一点，方便调 bug（可按需删掉或改成 Release）
if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 头文件所在目录
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/algorithms
)

# 主程序（目前所有实现都在 .hpp 里，是 header-only）
add_executable(rect_sampler_experiments
    main.cpp
)

# 如果将来你把 structures/algorithms 拆成 .cpp，可以在这里加上：
# add_executable(rect_sampler_experiments
#     main.cpp
#     src/geometry.cpp
#     src/structures.cpp
#     src/baseline1.cpp
#     src/baseline2.cpp
#     src/ours.cpp
# )
