cmake_minimum_required(VERSION 3.21)
project(GamePhysicsProject)

set(CMAKE_CXX_STANDARD 14)

add_executable(GamePhysicsProject main.cpp StableSolver.cpp Visualizer.cpp Visualizer.h externals/imgui/imgui.cpp wrapper/imgui_impl_glut.cpp wrapper/imgui_impl_opengl2.cpp externals/imgui/imgui_draw.cpp externals/imgui/imgui_tables.cpp externals/imgui/imgui_widgets.cpp)

# IMPORTANT: Remove -framework if you are on Windows or Linux
# You need OpenGL and GLUT installed globally on your system
target_link_libraries(GamePhysicsProject "-framework OpenGL" "-framework GLUT")
