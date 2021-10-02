#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <chrono>
#include <iostream>

#define GLAD_GL_IMPLEMENTATION
#include "gl_program.hpp"

#include "bool_terrain_generator.hpp"
#include "generator_perlin.hpp"
#include "terrain_manager.hpp"

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

  roller::SliceDirs sd = {{0, 0, 0},
                          {0.2673, -0.5246, 0.8083},
                          {-0.5345, 0.6172, 0.5774},
                          {0.8018, 0.5864, 0.1155},
                          1.5,
                          1.5,
                          96};
  const std::size_t numGradVecs = 4096;
  roller::TerrainManager<
      roller::BoolTerrainGenerator<roller::GeneratorPerlin<3>>, GLMesh>
      terrainManager({{{{64, 64, 64},
                        roller::getGradVecs(numGradVecs, 3, 1),
                        numGradVecs - 1,
                        3,
                        0.5}}},
                     16383);

  GLProgram prog;
  prog.compileProgram();

  std::size_t fpsCount = 0;
  auto beg = std::chrono::high_resolution_clock::now();
  while (!glfwWindowShouldClose(window)) {
    glfwSwapBuffers(window);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glfwPollEvents();

    // std::cout << "NEW FRAME" << std::endl << std::endl << std::endl <<
    // std::endl << std::endl;

    GLTrianglesRecorder fun = prog.makeTriangleRecorder();
    terrainManager.exportAllTriangles(sd, fun);
    v::DVec<16> projMat = roller::genProjMat(sd);
    prog.setProjMat(&projMat[0]);
    fun.renderAll();

    GLint err = glGetError();
    if (err) {
      std::cout << "OpenGL Error: " << err << std::endl;
    }

    sd.c += 0.2;
    // sd.c[2] -= 0.01;

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

