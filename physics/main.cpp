#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <chrono>
#include <iostream>

#define GLAD_GL_IMPLEMENTATION
#include "gl_program.hpp"

#include "generic_obj.hpp"
#include "prism.hpp"
#include "world.hpp"

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

  roller::SliceDirs sd = {{18, 8, -10}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0},
                          1.5,          1.5,       256};
  typedef roller::GenericPrimWrapper<roller::Prism<GLMesh>> PrimType;
  roller::World<PrimType, roller::GenericObj> world(9.81, 0.02);
  PrimType groundPrims[] = {
      {{v::DVec<3>{20, 1, 20}, v::DVec<3>{0, -10, 0}, 0}},
      {{v::DVec<3>{20, 1, 20}, v::DVec<3>{0, 10, 0}, 0}},
      {{v::DVec<3>{20, 10, 1}, v::DVec<3>{0, 0, -20}, 0}},
      {{v::DVec<3>{20, 10, 1}, v::DVec<3>{0, 0, 20}, 0}},
      {{v::DVec<3>{1, 10, 20}, v::DVec<3>{-20, 0, 0}, 0}},
      {{v::DVec<3>{1, 10, 20}, v::DVec<3>{20, 0, 0}, 0}},
  };
  groundPrims[0].doRender = false;
  const std::size_t groundPrimsLen = sizeof(groundPrims) / sizeof(PrimType);
  roller::CPhysInfo cpiTmp;
  roller::GenericObj groundObjTmp(groundPrims + 0, groundPrims + groundPrimsLen,
                                  cpiTmp);
  world.addObj(cpiTmp, std::move(groundObjTmp), groundPrims + 0,
               groundPrims + groundPrimsLen);
  PrimType prims[] = {
      {{v::DVec<3>{1, 1, 2}, v::DVec<3>{0, 0, 5}, 1}},
      {{v::DVec<3>{1, 1, 2}, v::DVec<3>{0, 0, -5}, 1}},
  };
  const std::size_t primsLen = sizeof(prims) / sizeof(PrimType);
  for (std::size_t i = 0; i < primsLen; i++) {
    roller::GenericObj tmpObj(prims + i, prims + i + 1, cpiTmp);
    cpiTmp.pi.lm[2] = i * 5;
    cpiTmp.updateAux();
    world.addObj(cpiTmp, std::move(tmpObj), prims + i, prims + i + 1);
  }

  GLProgram prog;
  prog.compileProgram();
  GLTrianglesRecorder trirec(prog);

  std::size_t fpsCount = 0;
  auto beg = std::chrono::high_resolution_clock::now();
  while (!glfwWindowShouldClose(window)) {
    glfwSwapBuffers(window);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glfwPollEvents();

    world.step(0.01, roller::genericImpulseCalc<PrimType>);

    auto renderPrim = [&sd, &trirec](PrimType &prim) -> void {
      prim.exportAllTriangles(sd, trirec);
    };
    world.applyToPrims(renderPrim);
    trirec.renderAll();
    trirec.clear();

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
