cmake_minimum_required(VERSION 3.21)
project(GamePhysicsProject)

set(CMAKE_CXX_STANDARD 14)

add_executable(GamePhysicsProject
        main.cpp
        StableSolver.cpp
        Visualizer.cpp
        wrapper/imgui_impl_glut.cpp
        wrapper/imgui_impl_opengl2.cpp
        externals/imgui/imgui.cpp
        externals/imgui/imgui_draw.cpp
        externals/imgui/imgui_tables.cpp
        externals/imgui/imgui_widgets.cpp)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

target_include_directories(GamePhysicsProject PUBLIC ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
target_link_libraries(GamePhysicsProject ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
