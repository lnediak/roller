#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <chrono>
#include <iostream>

#define GLAD_GL_IMPLEMENTATION
#include "gl_program.hpp"

#include "prism.hpp"
// #include "solver.hpp"

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

  roller::SliceDirs sd = {
      {15, 10, -20}, {-1, 0, 0}, {0, 0.7071, 0.7071}, {0, -0.7071, 0.7071}, 1.5,
      1.5,           128};
  std::vector<roller::Prism<GLMesh>> world;
  /*world.objs.emplace_back(v::DVec<3>{20, 30, 1}, v::DVec<3>{0, 33, -10}, 0);
  world.objs.emplace_back(v::DVec<3>{1, 2, 1}, v::DVec<3>{0, 8, -2}, 1);
  world.objs.emplace_back(v::DVec<3>{1, 1, 1}, v::DVec<3>{1.9, 8, 5}, 1);
  world.objs[2].pi.lm[2] = -1;
  world.objs[2].pi.am = {0, 1, 1};
  world.objs.emplace_back(v::DVec<3>{20, 30, 1}, v::DVec<3>{0, 33, 10}, 0);
  world.objs.emplace_back(v::DVec<3>{20, 1, 10}, v::DVec<3>{0, 33, 0}, 0);
  world.objs.emplace_back(v::DVec<3>{1, 30, 10}, v::DVec<3>{-20, 33, 0}, 0);
  world.objs.emplace_back(v::DVec<3>{1, 30, 10}, v::DVec<3>{20, 33, 0}, 0);
  world.objs.emplace_back(v::DVec<3>{20, 1, 10}, v::DVec<3>{0, 0, 0}, 0);
  world.objs.back().doRender = false;
  // world.objs.emplace_back(v::DVec<3>{1, 2, 1}, v::DVec<3>{0, 8, 2}, 1);
  */
  v::DVec<3> sl1 = {9.67501, 6.7767, 8.67646};
  roller::Pose pose1;
  pose1.p = {1.10539, -11.8452, -12.344};
  pose1.q = {0.643683, -0.473474, 0.57681, 0.16966};
  world.emplace_back(sl1 / 2, pose1.toWorldCoords(sl1 / 2), 0);
  world.back().pi.pose.q = pose1.q;
  v::DVec<3> sl2 = {0.145758, 7.73647, 2.88868};
  roller::Pose pose2;
  pose2.p = {-2.82045, 0.561555, -9.22222};
  pose2.q = {0.189117, -0.711096, -0.556603, 0.385708};
  world.emplace_back(sl2 / 2, pose2.toWorldCoords(sl2 / 2), 0);
  world.back().pi.pose.q = pose2.q;

  // roller::Solver<decltype(world) &> solver(world, 5);

  GLProgram prog;
  prog.compileProgram();
  GLTrianglesRecorder fun(prog);

  std::size_t fpsCount = 0;
  auto beg = std::chrono::high_resolution_clock::now();
  while (!glfwWindowShouldClose(window)) {
    glfwSwapBuffers(window);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glfwPollEvents();

    for (int spam = 0; spam < 3; spam++) {
      // solver.step(0.003);
    }
    sd.c[0] -= 0.01;
    for (std::size_t i = 0; i < world.size(); i++) {
      world[i].exportAllTriangles(sd, fun);
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

