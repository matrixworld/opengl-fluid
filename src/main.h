#ifndef __MAIN_H__
#define __MAIN_H__

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "fluid.h"

//moveMode:
#define VIEW_TRANSLATE 0
#define VIEW_ZOOM 1

int main(int argc, char **argv);
void disp(void); //called each time the screen is displayed
void reshape(int w, int h); //called each time window is reshaped
void keyb(unsigned char key, int x, int y); //called when key pressed
void mouseClick(int button, int state, int x, int y);
void mouseMove(int x, int y);
void redisp(int value);
void queuenext();

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))

#endif
