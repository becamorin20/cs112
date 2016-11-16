#ifndef INPUTMODULE_H
#define INPUTMODULE_H

#ifdef WIN32
#include <windows.h>
#endif

#include <stdlib.h>
#ifdef __APPLE__
#include <glut/glut.h>
#else
#include <GL/glut.h>
#endif


extern int lightposmodifier;
extern GLfloat current_pos[];
extern GLfloat light_position[];
#ifdef __cplusplus
extern "C" {
#endif
  
void readKeyboard(unsigned char key, int x, int y);
void readSpecialKeys(int key, int x, int y);
void mouseButtHandler(int button, int state, int x, int y);
void mouseMoveHandler(int x, int y);
void setUserView();
void transformLight();

#ifdef __cplusplus
}
#endif


#endif









