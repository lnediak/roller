#ifndef GL_PROGRAM_HPP_
#define GL_PROGRAM_HPP_

#include "glad/gl.h"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct GLProgram {

  GLuint prog;
  GLuint posLoc;
  GLuint colLoc;
  GLuint projMatLoc;

  GLProgram() : prog(0) {}
  GLProgram(const GLProgram &) = delete;
  GLProgram(GLProgram &&other) : prog(other.prog) { other.prog = 0; }
  GLProgram &operator=(const GLProgram &) = delete;
  GLProgram &operator=(GLProgram &&other) {
    if (prog != other.prog) {
      this->~GLProgram();
      prog = other.prog;
    }
    other.prog = 0;
    return *this;
  }
  ~GLProgram() {
    if (prog) {
      glDeleteProgram(prog);
    }
  }

  void compileProgram() {
    std::string vsrc = "#version 130\n"
                       "in vec3 pos;"
                       "in vec4 col;"
                       "out vec4 outCol;"
                       "uniform mat4 projMat;"
                       "void main() {"
                       "  gl_Position = projMat * vec4(pos, 1.0);"
                       "  outCol = col;"
                       "}";
    std::string fsrc = "#version 130\n"
                       "in vec4 outCol;"
                       "out vec4 fragCol;"
                       "void main() {"
                       "  fragCol = outCol;"
                       "}";
    GLuint vs = createShader(vsrc, GL_VERTEX_SHADER);
    GLuint fs = createShader(fsrc, GL_FRAGMENT_SHADER);

    GLuint pro = glCreateProgram();
    if (!pro) {
      throw std::runtime_error("Failed to create OpenGL program.");
    }
    glAttachShader(pro, vs);
    glAttachShader(pro, fs);
    glLinkProgram(pro);
    glDetachShader(pro, vs);
    glDetachShader(pro, fs);
    glDeleteShader(vs);
    glDeleteShader(fs);

    GLint llen = 0;
    glGetProgramiv(pro, GL_INFO_LOG_LENGTH, &llen);
    if (llen) {
      std::string infoLog(llen, ' ');
      glGetProgramInfoLog(pro, llen, NULL, &infoLog[0]);
      std::cerr << infoLog << std::endl;
    }

    GLint success = 0;
    glGetProgramiv(pro, GL_LINK_STATUS, &success);
    if (success == GL_FALSE) {
      glDeleteProgram(pro);
      throw std::runtime_error("Failed to link program.");
    }
    if (prog) {
      glDeleteProgram(prog);
    }
    prog = pro;
    glUseProgram(prog);

    posLoc = glGetAttribLocation(prog, "pos");
    colLoc = glGetAttribLocation(prog, "col");
    projMatLoc = glGetUniformLocation(prog, "projMat");
  }

  /// row-major please
  void setProjMat(const double *p) {
    GLfloat mat[16];
    for (int i = 0; i < 16; i++) {
      mat[i] = p[i];
    }
    glUniformMatrix4fv(projMatLoc, 1, GL_TRUE, mat);
  }

private:
  GLuint createShader(const std::string &src, GLenum stype) const {
    const char *csrc = src.c_str();
    const GLint slen = src.length();

    GLuint shad = glCreateShader(stype);
    if (!shad) {
      throw std::runtime_error("Failed to create OpenGL shader.");
    }
    glShaderSource(shad, 1, &csrc, &slen);
    glCompileShader(shad);

    GLint llen = 0;
    glGetShaderiv(shad, GL_INFO_LOG_LENGTH, &llen);
    if (llen) {
      std::string infoLog(llen, ' ');
      glGetShaderInfoLog(shad, llen, NULL, &infoLog[0]);
      std::cerr << infoLog << std::endl;
    }

    GLint success = 0;
    glGetShaderiv(shad, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
      glDeleteShader(shad);
      throw std::runtime_error(std::string() +
                               "Failed to compile shader. Source:\n" + src);
    }
    return shad;
  }
};

struct GLMesh {

  GLuint vao;
  GLuint buf;
  std::size_t nm = 0;

  GLMesh() : vao(0), buf(0) {}
  explicit GLMesh(int) {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &buf);
  }
  GLMesh(const GLMesh &) = delete;
  GLMesh(GLMesh &&other) : vao(other.vao), buf(other.buf), nm(other.nm) {
    other.vao = 0;
    other.buf = 0;
  }
  GLMesh &operator=(const GLMesh &) = delete;
  GLMesh &operator=(GLMesh &&other) {
    if (vao != other.vao) {
      this->~GLMesh();
      vao = other.vao;
      buf = other.buf;
      nm = other.nm;
    }
    other.vao = 0;
    other.buf = 0;
    return *this;
  }

  operator bool() const { return !!vao; }

  ~GLMesh() {
    if (buf) {
      glDeleteBuffers(1, &buf);
    }
    if (vao) {
      glDeleteVertexArrays(1, &vao);
    }
  }
};

struct GLTrianglesRecorder {

  struct Action {
    bool isMesh;
    union {
      GLMesh *mesh;
      double projMat[16];
    };
  };
  std::vector<Action> meshes;
  GLProgram &prog;

  GLTrianglesRecorder(GLProgram &prog) : prog(prog) {}

  void clear() { meshes.clear(); }

  void operator()(GLMesh &tag, const float *triangles, std::size_t sz) {
    if (!sz) {
      tag.~GLMesh();
      new (&tag) GLMesh();
      return;
    }
    if (!tag) {
      tag.~GLMesh();
      new (&tag) GLMesh(0);
    }
    Action act;
    act.isMesh = true;
    act.mesh = &tag;
    meshes.push_back(act);
    tag.nm = sz / 4;
    glBindBuffer(GL_ARRAY_BUFFER, tag.buf);
    glBufferData(GL_ARRAY_BUFFER, sz * sizeof(float), triangles,
                 GL_STATIC_DRAW);
    glVertexAttribPointer(prog.posLoc, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float),
                          0);
    glVertexAttribPointer(prog.colLoc, 4, GL_UNSIGNED_BYTE, GL_TRUE,
                          4 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(prog.posLoc);
    glEnableVertexAttribArray(prog.colLoc);
  }
  void operator()(GLMesh &tag) {
    if (tag) {
      Action act;
      act.isMesh = true;
      act.mesh = &tag;
      meshes.push_back(act);
    }
  }
  void operator()(const double *projMat) {
    Action act;
    act.isMesh = false;
    std::memcpy(act.projMat, projMat, 16 * sizeof(double));
    meshes.push_back(act);
  }

  void renderAll() const {
    for (const Action &a : meshes) {
      if (a.isMesh) {
        glBindVertexArray(a.mesh->vao);
        glDrawArrays(GL_TRIANGLES, 0, a.mesh->nm);
      } else {
        prog.setProjMat(a.projMat);
      }
    }
  }
};

#endif // GL_PROGRAM_HPP_

