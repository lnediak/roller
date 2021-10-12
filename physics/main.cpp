#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <chrono>
#include <iostream>

#define GLAD_GL_IMPLEMENTATION
#include "gl_program.hpp"

#include "prism.hpp"
#include "solver.hpp"

namespace {

void errCallback(int, const char *err) { std::cerr << err << std::endl; }

} // namespace

int main() {
  glfwSetErrorCallback(&errCallback);
  if (!glfwInit()) {
    std::cout << "GLFW Initialization Failed." << std::endl;
    return 1;
  }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  const std::size_t width = 1000, height = 1000;
  GLFWwindow *window =
      glfwCreateWindow(width, height, "test lel", nullptr, nullptr);
  if (!window) {
    std::cout << "GLFW window creation failed." << std::endl;
    return 1;
  }
  glfwMakeContextCurrent(window);
  if (!gladLoadGL((GLADloadfunc)glfwGetProcAddress)) {
    std::cout << "GLAD loading failed." << std::endl;
  }
  glClearColor(0.f, 0.f, 0.f, 1.f);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glFrontFace(GL_CCW);

  roller::SliceDirs sd = {{0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0},
                          1.5,       1.5,       128};
  roller::OneObjWorld<roller::Prism<GLMesh>> world;
  world.objs.emplace_back(v::DVec<3>{10, 10, 1}, v::DVec<3>{0, 0, -5}, 0);
  world.objs.emplace_back(v::DVec<3>{1, 2, 1}, v::DVec<3>{0, 8, 5}, 1);
  // world.objs.emplace_back(v::DVec<3>{1, 2, 1}, v::DVec<3>{0, 8, 2}, 1);
  roller::Solver<decltype(world) &> solver(world, 1);

  GLProgram prog;
  prog.compileProgram();
  GLTrianglesRecorder fun(prog);

  std::size_t fpsCount = 0;
  auto beg = std::chrono::high_resolution_clock::now();
  while (!glfwWindowShouldClose(window)) {
    glfwSwapBuffers(window);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glfwPollEvents();

    solver.step(0.01);
    for (std::size_t i = 0; i < world.numObjs(); i++) {
      world.getObj(i).exportAllTriangles(sd, fun);
    }
    fun.renderAll();
    fun.clear();

    GLint err = glGetError();
    if (err) {
      std::cout << "OpenGL Error: " << err << std::endl;
    }

    auto dur = std::chrono::high_resolution_clock::now() - beg;
    double secs =
        std::chrono::duration_cast<std::chrono::duration<double>>(dur).count();
    fpsCount++;
    if (secs >= 1) {
      std::cout << fpsCount << "fps" << std::endl;
      fpsCount = 0;
      beg = std::chrono::high_resolution_clock::now();
    }
  }
}

