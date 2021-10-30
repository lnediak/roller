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

  roller::SliceDirs sd = {{100, 80, -120},
                          {-1, 0, 0},
                          {0, 0.7071, 0.7071},
                          {0, -0.7071, 0.7071},
                          1.5,
                          1.5,
                          256};
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
  v::DVec<3> b1 = {99.09796747704624, -94.23041213362805, 32.9722090258958};
  v::DVec<3> sl1 = {13.309027721583904, 72.250546882134714, 78.944365750799378};
  v::DVec<3> x1 = {0.56412915188846058, 0.58984828315316462,
                   0.57778655475086138};
  v::DVec<3> y1 = {0.68630241956766769, 0.054092773712392336,
                   -0.72530197899058901};
  v::DVec<3> z1 = {-0.45907220443539276, 0.80570030079016597,
                   -0.37429899334818539};
  roller::Pose pose1;
  pose1.p = b1;
  pose1.q = roller::rotMatToQuaternion(roller::DMat3x3{x1, y1, z1}.transpose());
  v::DVec<4> qtest = roller::rotMatToQuaternion(pose1.toRotationMatrix());
  if (v::norm2(pose1.q - qtest) >= 1e-8) {
    std::cerr << "CHECK YOUR rotMatToQuaternion!!!!!!" << std::endl;
    return 1;
  }
  world.emplace_back(sl1 / 2, pose1.toWorldCoords(sl1 / 2), 0);
  world.back().pi.pose.q = pose1.q;
  v::DVec<3> b2 = {69.661213826329458, -83.031677045612597,
                   -62.664839198712251};
  v::DVec<3> sl2 = {76.28469153294175, 19.901083146143105, 78.650146791614162};
  v::DVec<3> x2 = {0.45711062741863417, 0.35964249065534648,
                   -0.81345384205076765};
  v::DVec<3> y2 = {0.00023988627748436353, 0.91454966018202322,
                   0.4044735609463963};
  v::DVec<3> z2 = {0.8894098136891172, -0.18508429961858811,
                   0.41796409567593762};
  roller::Pose pose2;
  pose2.p = b2;
  pose2.q = roller::rotMatToQuaternion(roller::DMat3x3{x2, y2, z2}.transpose());
  qtest = roller::rotMatToQuaternion(pose2.toRotationMatrix());
  if (v::norm2(pose2.q - qtest) >= 1e-8) {
    std::cerr << "CHECK YOUR rotMatToQuaternion!!!!!!" << std::endl;
    return 1;
  }
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
    sd.c[0] -= 0.1;
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

