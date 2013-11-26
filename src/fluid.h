#ifndef __FLUID_H__
#define __FLUID_H__

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SIGN(A) (((A)>0)?(1):(((A)<0)?(-1):(0)))
#define ABS(A) ((A)<0?(-(A)):(A))

typedef double real;

typedef struct {
    int nx, ny;
    real *fdiv;
    real xmin,xmax;
    real ymin,ymax;
    real *ubar, *vbar, *vort;
    real *fx, *fy;
    real *u, *u1;
    real *v, *v1;
    real *T, *T1;
    real *rho, *rho1;
    real *div, *p;
    real dt, hx, hy;
    real dx,dy;
    real nu,eps,buoy;
    
} fluid;

fluid* alloc_fluid(int pnx, int pny);
void free_fluid(fluid* f);
void fluid_calc_spacing(fluid* f);
inline real fluid_interpolate_u(fluid* f, real x, real y);
inline real fluid_interpolate_v(fluid* f, real x, real y);
inline real fluid_interpolate_T(fluid* f, real x, real y);
void fluid_advect_u(fluid* f);
void fluid_advect_v(fluid* f);
void fluid_advect_T(fluid* f);
void fluid_project(fluid* f,int stepmax,int eastwall, int westwall, int northwall, int southwall);
real fluid_calc_div(fluid* f);
void fluid_draw(fluid* f,real scale);
void fluid_get_cell(fluid* f, real x, real y, int* i, int* j);
void fluid_zero_temperature(fluid* f);
void fluid_lower_temperature(fluid* f);
void fluid_calc_vbar(fluid* f);
void fluid_calc_vorticity(fluid* f);
void fluid_calc_forces(fluid* f,int vortOn,int buoyOn);
void fluid_add_forces(fluid* f);
void fluid_zero_fdiv(fluid* f);
void fluid_read_image(fluid* f);
void fluid_zero_velocity(fluid* f);

#endif /*__FLUID_H__*/
