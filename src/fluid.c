#include "fluid.h" 

//#define LINEAR_INTERPOLATION
//#define LINEAR_INTERPOLATION_T

//Obvious enough... allocate the arrays and such
fluid* alloc_fluid(int pnx, int pny) {
    fprintf(stderr,"Allocating fluid...");
    int i,j;
    fluid* f = malloc(sizeof(fluid));

    f->nx = pnx;
    f->ny = pny;
    f->u = malloc(sizeof(real)*(f->nx+1)*f->ny);
    f->u1 = malloc(sizeof(real)*(f->nx+1)*f->ny);
    f->v = malloc(sizeof(real)*f->nx*(f->ny+1));
    f->v1 = malloc(sizeof(real)*f->nx*(f->ny+1));
    f->T = malloc(sizeof(real)*f->nx*f->ny);
    f->T1 = malloc(sizeof(real)*f->nx*f->ny);
    //f->rho = malloc(sizeof(real)*f->nx*f->ny);
    //f->rho1 = malloc(sizeof(real)*f->nx*f->ny);
    f->p = malloc(sizeof(real)*f->nx*f->ny);
    f->div = malloc(sizeof(real)*f->nx*f->ny);
    f->fdiv = malloc(sizeof(real)*f->nx*f->ny);
    f->ubar = malloc(sizeof(real)*f->nx*f->ny);
    f->vbar = malloc(sizeof(real)*f->nx*f->ny);
    f->vort = malloc(sizeof(real)*f->nx*f->ny);
    f->fx = malloc(sizeof(real)*f->nx*f->ny);
    f->fy = malloc(sizeof(real)*f->nx*f->ny);

    for(i=0;i<f->nx+1;i++) {
        for(j=0;j<f->ny;j++) {
            f->u[i+(f->nx+1)*j] = (drand48()-0.5)*0.0001;
        }
    }
    for(i=0;i<f->nx;i++) {
        for(j=0;j<f->ny+1;j++) {
            f->v[i+f->nx*j] = (drand48()-0.5)*0.0001;
        }
    }
    for(i=0;i<f->nx;i++) {
        for(j=0;j<f->ny;j++) {
            f->T[i+f->nx*j] = 0.0;
        }
    }
    fprintf(stderr,"done.\n");
    return f;
}

void fluid_zero_velocity(fluid* f) {
    int i;
    for(i=0;i<(f->nx+1)*f->ny;i++) {
        f->u[i] = 0.0;
    }
    for(i=0;i<f->nx*(f->ny+1);i++) {
        f->v[i] = 0.0;
    }
    for(i=0;i<f->nx*f->ny;i++) {
        f->p[i] = 0.0;
    }
}

void fluid_zero_temperature(fluid* f) {
    int i,j;
    for(j=0;j<f->ny;j++) {
        for(i=0;i<f->nx;i++) {
            //if(j<f->ny/2) {
                //f->T[i+f->nx*j] = 1.0;
            //} else {
                f->T[i+f->nx*j] = 0.0;
            //}
        }
    }
}

void fluid_read_image(fluid* f) {
    int i,j,k;
    FILE* image=fopen("cy.dat","r");
    for(j=f->ny-1;j>=0;j--) {
        for(i=0;i<f->nx;i++) {
            fscanf(image,"%i",&k);
            f->T[i+f->nx*j]=1.0*k/255.0;
        }
    }
    fclose(image);
}

//Free the arrays when done...
void free_fluid(fluid* f) {
    if(f) {
        fprintf(stderr,"Freeing fluid...");
        free(f->u);
        free(f->u1);
        free(f->v);
        free(f->v1);
        free(f->T);
        free(f->T1);
        //free(f->rho);
        //free(f->rho1);
        free(f->p);
        free(f->div);
        free(f->ubar);
        free(f->vbar);
        free(f->vort);
        free(f->fx);
        free(f->fy);
        free(f->fdiv);
        free(f);
        fprintf(stderr,"done.\n");
    } else {
        fprintf(stderr,"Error: fluid already freed\n");
    }
}

//Calculate the cell-centered average velocity from face-centered velocities
void fluid_calc_vbar(fluid* f) {
    int i,j;
    for(j=0;j<f->ny;j++) {
        for(i=0;i<f->nx;i++) {
            f->ubar[i+f->nx*j] = 0.5*(f->u[i+(f->nx+1)*j] + f->u[i+1+(f->nx+1)*j]);
            f->vbar[i+f->nx*j] = 0.5*(f->v[i+f->nx*j] + f->v[i+f->nx*(j+1)]);
        }
    }
}

//Calculate vorticity from the cell-centered velocities
void fluid_calc_vorticity(fluid* f) {
    int i,j;
    for(j=1;j<f->ny-1;j++) {
        for(i=1;i<f->nx-1;i++) {
            f->vort[i+f->nx*j] = 
                (f->vbar[(i+1)+f->nx*j]-f->vbar[(i-1)+f->nx*j])/f->hx +
                (f->ubar[i+f->nx*(j-1)]-f->ubar[i+f->nx*(j+1)])/f->hy;
        }
    }
}

//Calculate vorticity confinement and buoyancy forces and add
//them to the fx and fy arrays.
void fluid_calc_forces(fluid* f,int vortOn,int buoyOn) {
    int i,j;
    real gomx,gomy,gomm,Nx,Ny;
    memset(f->fx,0,sizeof(real)*f->nx*f->ny);
    memset(f->fy,0,sizeof(real)*f->nx*f->ny);
    if(vortOn) {
        for(i=2;i<f->nx-2;i++) {
            for(j=2;j<f->ny-2;j++) {
                gomx = (ABS(f->vort[(i+1)+f->nx*j])-ABS(f->vort[(i-1)+f->nx*j]))/f->hx*0.5;
                gomy = (ABS(f->vort[i+f->nx*(j+1)])-ABS(f->vort[i+f->nx*(j-1)]))/f->hy*0.5;
                gomm = 1.0/sqrt(gomx*gomx+gomy*gomy+1.0e-7)*f->vort[i+f->nx*j];
                Nx = gomx*gomm;
                Ny = gomy*gomm;
                f->fx[i+f->nx*j] = Ny*f->nx*f->eps;
                f->fy[i+f->nx*j] = -Nx*f->ny*f->eps;
            }
        }
    }
    if(buoyOn) {
        for(j=0;j<f->ny;j++) {
            for(i=0;i<f->nx;i++) {
                f->fy[i+f->nx*j] += pow(ABS(f->T[i+f->nx*j]),1.00)*f->buoy/MAX(f->hx,f->hy);
            }
        }
    }
}

//Add forces to the fluid.  It should have an additional factor of
// (1/2)dt, but it doesn't really make a difference since the
// coefficients are arbitrary anyway.
void fluid_add_forces(fluid* f) {
    int i,j;
    for(i=1;i<f->nx;i++) {
        for(j=0;j<f->ny;j++) {
            f->u[i+(f->nx+1)*j] += f->fx[i-1+f->nx*j]+f->fx[i+f->nx*j];
        }
    }
    for(i=0;i<f->nx;i++) {
        for(j=1;j<f->ny;j++) {
            f->v[i+f->nx*j] += f->fy[i+f->nx*(j-1)]+f->fy[i+f->nx*j];
        }
    }
}

//Calculate the grid spacing.  Self explanatory.
void fluid_calc_spacing(fluid* f) {
    f->dx = (f->xmax-f->xmin);
    f->dy = (f->ymax-f->ymin);
    f->hx = f->dx/f->nx;
    f->hy = f->dy/f->ny;
}

//Figure out which cell a given (x,y) pair is inside.  Return
//the result in (i,j)
void fluid_get_cell(fluid* f, real x, real y, int* i, int* j) {
    real nx,ny;
    int ix,iy;
    nx=(x-f->xmin)/f->dx*f->nx;
    ny=(y-f->ymin)/f->dy*f->ny;
    ny=(ny<0.0)?(0.0):(ny);
    ny=(ny>f->ny-1)?(f->ny-1):(ny);
    nx=(nx<0.0)?(0.0):(nx);
    nx=(nx>f->nx-1)?(f->nx-1):(nx);
    ix=(int)nx;
    iy=(int)ny;
    iy=(iy<0)?(0):(iy);
    iy=(iy>f->ny-1)?(f->ny-1):(iy);
    ix=(ix<0)?(0):(ix);
    ix=(ix>f->nx-1)?(f->nx-1):(ix);
    *i=ix;
    *j=iy;
}

//Cubic interpolation in the range from p(0) to p(1) with exterior points p(-1), p(2)
//This version is just a four-point cubic polynomial fit, but it doesn't preserve
//monotonicity
 real cubic_interpolation(real pn1, real p0, real p1, real p2, real t) {
    return (-0.166666666666666666*pn1 + 0.5*p0 - 0.5*p1 + 0.166666666666666666*p2)*t*t*t +
           (0.5*pn1 - p0 + 0.5*p1)*t*t +
           (-0.333333333333333333*pn1 - 0.5*p0 + p1 - 0.166666666666666666*p2)*t +
           p0;
}

//Hermite Monotonic Cubic interpolation from p(0) to p(1).  Fit through two points using
//the exterior points to calculate slopes.  Copied from "Visual Simulation of Smoke",
//by Fedkiw, Stam, and Jensen.
 real cubic_monotonic_interpolation(real pn1, real p0, real p1, real p2, real t) {
    real dk,dk1,delk,a3,a2;
    dk=(p1-pn1)*0.5;
    dk1=(p2-p0)*0.5;
    delk=p1-p0;
    if(SIGN(dk)!=SIGN(delk))
        dk=0.0;
    if(SIGN(dk1)!=SIGN(delk))
        dk1=0.0;
    a3=dk+dk1-2.0*delk;
    a2=3.0*delk-2.0*dk-dk1;
    return a3*t*t*t+a2*t*t+dk*t+p0;
}

//Interpolate the u-component of velocity subject to du/dn=0 on boundaries and with
//cubic monotonic interpolation
 real fluid_interpolate_cubic_u(fluid* f, real x, real y) {
    int ix,iy,ixm,ixp,iym,iyp;
    real nx,ny,tx,ty,f1,f2,f3,f4;
    nx=(x-f->xmin)/f->dx*f->nx;
    ny=(y-f->ymin)/f->dy*f->ny-0.5;
    ny=(ny<0.0)?(0.0):(ny);
    ny=(ny>f->ny-1)?(f->ny-1):(ny);
    nx=(nx<0.0)?(0.0):(nx);
    nx=(nx>f->nx)?(f->nx):(nx);
    ix=(int)nx;
    iy=(int)ny;
    iy=(iy<0)?(0):(iy);
    iy=(iy>f->ny-2)?(f->ny-2):(iy);
    ix=(ix<0)?(0):(ix);
    ix=(ix>f->nx-1)?(f->nx-1):(ix);
    tx=nx-ix;
    ty=ny-iy;
    ixm=MAX(ix-1,0);
    ixp=MIN(ix+2,f->nx);
    iym=MAX(iy-1,0);
    iyp=MIN(iy+2,f->ny-1);
    f1=cubic_monotonic_interpolation(f->u[ixm+(f->nx+1)*iym],f->u[ix+(f->nx+1)*iym],f->u[ix+1+(f->nx+1)*iym],f->u[ixp+(f->nx+1)*iym],tx);
    f2=cubic_monotonic_interpolation(f->u[ixm+(f->nx+1)*iy],f->u[ix+(f->nx+1)*iy],f->u[ix+1+(f->nx+1)*iy],f->u[ixp+(f->nx+1)*iy],tx);
    f3=cubic_monotonic_interpolation(f->u[ixm+(f->nx+1)*(iy+1)],f->u[ix+(f->nx+1)*(iy+1)],f->u[ix+1+(f->nx+1)*(iy+1)],f->u[ixp+(f->nx+1)*(iy+1)],tx);
    f4=cubic_monotonic_interpolation(f->u[ixm+(f->nx+1)*iyp],f->u[ix+(f->nx+1)*iyp],f->u[ix+1+(f->nx+1)*iyp],f->u[ixp+(f->nx+1)*iyp],tx);
    return cubic_monotonic_interpolation(f1,f2,f3,f4,ty);
}

//Interpolate the u-component of velocity subject to du/dn=0 on boundaries
//and with linear interpolation
 real fluid_interpolate_u(fluid* f, real x, real y) {
    real nx,ny,tx,ty,f11,f12,f21,f22,f1,f2;
    int ix,iy;
    nx=(x-f->xmin)/f->dx*f->nx;
    ny=(y-f->ymin)/f->dy*f->ny-0.5;
    ny=(ny<0.0)?(0.0):(ny);
    ny=(ny>f->ny-1)?(f->ny-1):(ny);
    nx=(nx<0.0)?(0.0):(nx);
    nx=(nx>f->nx)?(f->nx):(nx);
    ix=(int)nx;
    iy=(int)ny;
    iy=(iy<0)?(0):(iy);
    iy=(iy>f->ny-2)?(f->ny-2):(iy);
    ix=(ix<0)?(0):(ix);
    ix=(ix>f->nx-1)?(f->nx-1):(ix);
    tx=nx-ix;
    ty=ny-iy;
    f11=f->u[ix+(f->nx+1)*iy];
    f21=f->u[ix+1+(f->nx+1)*iy];
    f12=f->u[ix+(f->nx+1)*(iy+1)];
    f22=f->u[ix+1+(f->nx+1)*(iy+1)];
    f1=f11+(f21-f11)*tx;
    f2=f12+(f22-f12)*tx;
    return f1+(f2-f1)*ty;
}

//Cubic interpolation for v
 real fluid_interpolate_cubic_v(fluid* f, real x, real y) {
    int ix,iy,ixm,ixp,iym,iyp;
    real nx,ny,tx,ty,f1,f2,f3,f4;
    nx=(x-f->xmin)/f->dx*f->nx-0.5;
    ny=(y-f->ymin)/f->dy*f->ny;
    ny=(ny<0.0)?(0.0):(ny);
    ny=(ny>f->ny)?(f->ny):(ny);
    nx=(nx<0.0)?(0.0):(nx);
    nx=(nx>f->nx-1)?(f->nx-1):(nx);
    ix=(int)nx;
    iy=(int)ny;
    ix=(ix<0)?(0):(ix);
    ix=(ix>f->nx-2)?(f->nx-2):(ix);
    iy=(iy<0)?(0):(iy);
    iy=(iy>f->ny-1)?(f->ny-1):(iy);
    tx=nx-ix;
    ty=ny-iy;
    ixm=MAX(ix-1,0);
    ixp=MIN(ix+2,f->nx-1);
    iym=MAX(iy-1,0);
    iyp=MIN(iy+2,f->ny);
    f1=cubic_monotonic_interpolation(f->v[ixm+f->nx*iym],f->v[ix+f->nx*iym],f->v[ix+1+f->nx*iym],f->v[ixp+f->nx*iym],tx);
    f2=cubic_monotonic_interpolation(f->v[ixm+f->nx*iy],f->v[ix+f->nx*iy],f->v[ix+1+f->nx*iy],f->v[ixp+f->nx*iy],tx);
    f3=cubic_monotonic_interpolation(f->v[ixm+f->nx*(iy+1)],f->v[ix+f->nx*(iy+1)],f->v[ix+1+f->nx*(iy+1)],f->v[ixp+f->nx*(iy+1)],tx);
    f4=cubic_monotonic_interpolation(f->v[ixm+f->nx*iyp],f->v[ix+f->nx*iyp],f->v[ix+1+f->nx*iyp],f->v[ixp+f->nx*iyp],tx);
    return cubic_monotonic_interpolation(f1,f2,f3,f4,ty);
}

//Linear interpolation for v
 real fluid_interpolate_v(fluid* f, real x, real y) {
    real nx,ny,tx,ty,f11,f12,f21,f22,f1,f2;
    int ix,iy;
    nx=(x-f->xmin)/f->dx*f->nx-0.5;
    ny=(y-f->ymin)/f->dy*f->ny;
    ny=(ny<0.0)?(0.0):(ny);
    ny=(ny>f->ny)?(f->ny):(ny);
    nx=(nx<0.0)?(0.0):(nx);
    nx=(nx>f->nx-1)?(f->nx-1):(nx);
    ix=(int)nx;
    iy=(int)ny;
    ix=(ix<0)?(0):(ix);
    ix=(ix>f->nx-2)?(f->nx-2):(ix);
    iy=(iy<0)?(0):(iy);
    iy=(iy>f->ny-1)?(f->ny-1):(iy);
    tx=nx-ix;
    ty=ny-iy;
    f11=f->v[ix+f->nx*iy];
    f21=f->v[ix+1+f->nx*iy];
    f12=f->v[ix+f->nx*(iy+1)];
    f22=f->v[ix+1+f->nx*(iy+1)];
    f1=f11+(f21-f11)*tx;
    f2=f12+(f22-f12)*tx;
    return f1+(f2-f1)*ty;
}

//Cubic interpolation for T
 real fluid_interpolate_cubic_T(fluid* f, real x, real y) {
    real nx,ny,tx,ty,f1,f2,f3,f4;
    int ix,iy,ixm,ixp,iym,iyp;
    nx=(x-f->xmin)/f->dx*f->nx-0.5;
    ny=(y-f->ymin)/f->dy*f->ny-0.5;
    ny=(ny<0.0)?(0.0):(ny);
    ny=(ny>f->ny-1)?(f->ny-1):(ny);
    nx=(nx<0.0)?(0.0):(nx);
    nx=(nx>f->nx-1)?(f->nx-1):(nx);
    ix=(int)nx;
    iy=(int)ny;
    ix=(ix<0)?(0):(ix);
    ix=(ix>f->nx-2)?(f->nx-2):(ix);
    iy=(iy<0)?(0):(iy);
    iy=(iy>f->ny-2)?(f->ny-2):(iy);
    tx=nx-ix;
    ty=ny-iy;
    ixm=MAX(ix-1,0);
    ixp=MIN(ix+2,f->nx-1);
    iym=MAX(iy-1,0);
    iyp=MIN(iy+2,f->ny-1);
    f1=cubic_monotonic_interpolation(f->T[ixm+f->nx*iym],f->T[ix+f->nx*iym],f->T[ix+1+f->nx*iym],f->T[ixp+f->nx*iym],tx);
    f2=cubic_monotonic_interpolation(f->T[ixm+f->nx*iy],f->T[ix+f->nx*iy],f->T[ix+1+f->nx*iy],f->T[ixp+f->nx*iy],tx);
    f3=cubic_monotonic_interpolation(f->T[ixm+f->nx*(iy+1)],f->T[ix+f->nx*(iy+1)],f->T[ix+1+f->nx*(iy+1)],f->T[ixp+f->nx*(iy+1)],tx);
    f4=cubic_monotonic_interpolation(f->T[ixm+f->nx*iyp],f->T[ix+f->nx*iyp],f->T[ix+1+f->nx*iyp],f->T[ixp+f->nx*iyp],tx);
    return cubic_interpolation(f1,f2,f3,f4,ty);
}

//Linear interpolation for T
 real fluid_interpolate_T(fluid* f, real x, real y) {
    real nx,ny,tx,ty,f11,f12,f21,f22,f1,f2;
    int ix,iy;
    nx=(x-f->xmin)/f->dx*f->nx-0.5;
    ny=(y-f->ymin)/f->dy*f->ny-0.5;
    ny=(ny<0.0)?(0.0):(ny);
    ny=(ny>f->ny-1)?(f->ny-1):(ny);
    nx=(nx<0.0)?(0.0):(nx);
    nx=(nx>f->nx-1)?(f->nx-1):(nx);
    ix=(int)nx;
    iy=(int)ny;
    ix=(ix<0)?(0):(ix);
    ix=(ix>f->nx-2)?(f->nx-2):(ix);
    iy=(iy<0)?(0):(iy);
    iy=(iy>f->ny-2)?(f->ny-2):(iy);
    tx=nx-ix;
    ty=ny-iy;
    f11=f->T[ix+f->nx*iy];
    f21=f->T[ix+1+f->nx*iy];
    f12=f->T[ix+f->nx*(iy+1)];
    f22=f->T[ix+1+f->nx*(iy+1)];
    f1=f11+(f21-f11)*tx;
    f2=f12+(f22-f12)*tx;
    return f1+(f2-f1)*ty;
}

//Advect the u-component of velocity and put into u1.
void fluid_advect_u(fluid* f) {
    int i,j;
    real x0,y0,u0,v0,x1,y1,u1,v1,x2,y2;
    for(j=0;j<f->ny;j++) {
        y0=f->ymin+f->hy*(j+0.5);
        for(i=0;i<f->nx+1;i++) {
            x0=f->xmin+f->hx*i;
            u0=f->u[i+(f->nx+1)*j];
#ifdef LINEAR_INTERPOLATION
            //Linear interpolation, Euler integration:
            v0=fluid_interpolate_v(f,x0,y0);
            //x1=x0-f->dt*u0*0.5;
            //y1=y0-f->dt*v0*0.5;
            //u1=fluid_interpolate_u(f,x1,y1);
            //v1=fluid_interpolate_v(f,x1,y1);
            x2=x0-f->dt*u0;
            y2=y0-f->dt*v0;
            f->u1[i+(f->nx+1)*j]=fluid_interpolate_u(f,x2,y2);
#else
            //Cubic interpolation, RK2 integration
            v0=fluid_interpolate_cubic_v(f,x0,y0);
            x1=x0-f->dt*u0*0.5;
            y1=y0-f->dt*v0*0.5;
            u1=fluid_interpolate_cubic_u(f,x1,y1);
            v1=fluid_interpolate_cubic_v(f,x1,y1);
            x2=x0-f->dt*u1;
            y2=y0-f->dt*v1;
            f->u1[i+(f->nx+1)*j]=fluid_interpolate_cubic_u(f,x2,y2);
#endif
        }
    }
}

//Advect the v-component of velocity and put it into v1
void fluid_advect_v(fluid* f) {
    int i,j;
    real x0,y0,u0,v0,x1,y1,u1,v1,x2,y2;
    for(j=0;j<f->ny+1;j++) {
        y0=f->ymin+f->hy*j;
        for(i=0;i<f->nx;i++) {
            x0=f->xmin+f->hx*(i+0.5);
#ifdef LINEAR_INTERPOLATION
            //Linear interpolation, Euler integration
            u0=fluid_interpolate_u(f,x0,y0);
            v0=f->v[i+f->nx*j];
            //x1=x0-f->dt*u0*0.5;
            //y1=y0-f->dt*v0*0.5;
            //u1=fluid_interpolate_u(f,x1,y1);
            //v1=fluid_interpolate_v(f,x1,y1);
            x2=x0-f->dt*u0;
            y2=y0-f->dt*v0;
            f->v1[i+f->nx*j]=fluid_interpolate_v(f,x2,y2);
#else
            //Cubic interpolation, RK2 integration
            u0=fluid_interpolate_cubic_u(f,x0,y0);
            v0=f->v[i+f->nx*j];
            x1=x0-f->dt*u0*0.5;
            y1=y0-f->dt*v0*0.5;
            u1=fluid_interpolate_cubic_u(f,x1,y1);
            v1=fluid_interpolate_cubic_v(f,x1,y1);
            x2=x0-f->dt*u1;
            y2=y0-f->dt*v1;
            f->v1[i+f->nx*j]=fluid_interpolate_cubic_v(f,x2,y2);
#endif
        }
    }
}

void fluid_advect_T(fluid* f) {
    int i,j;
    real x0,y0,u0,v0,x1,y1,u1,v1,x2,y2;
    for(j=0;j<f->ny;j++) {
        y0=f->ymin+f->hy*(j+0.5);
        for(i=0;i<f->nx;i++) {
            x0=f->xmin+f->hx*(i+0.5);
#ifdef LINEAR_INTERPOLATION
            //Linear interpolation, Euler integration
            u0=fluid_interpolate_u(f,x0,y0);
            v0=f->v[i+f->nx*j];
            //x1=x0-f->dt*u0*0.5;
            //y1=y0-f->dt*v0*0.5;
            //u1=fluid_interpolate_u(f,x1,y1);
            //v1=fluid_interpolate_v(f,x1,y1);
            x2=x0-f->dt*u0;
            y2=y0-f->dt*v0;
            f->T1[i+f->nx*j]=fluid_interpolate_T(f,x2,y2);
#else
            //Cubic interpolation, RK2 integration
            u0=fluid_interpolate_cubic_u(f,x0,y0);
            v0=f->v[i+f->nx*j];
            x1=x0-f->dt*u0*0.5;
            y1=y0-f->dt*v0*0.5;
            u1=fluid_interpolate_cubic_u(f,x1,y1);
            v1=fluid_interpolate_cubic_v(f,x1,y1);
            x2=x0-f->dt*u1;
            y2=y0-f->dt*v1;
#ifdef LINEAR_INTERPOLATION_T
            f->T1[i+f->nx*j]=fluid_interpolate_T(f,x2,y2);
#else
            f->T1[i+f->nx*j]=fluid_interpolate_cubic_T(f,x2,y2);
#endif
#endif
        }
    }
}

//Calculate the divergence and add in the manually-specified divergence field, fdiv
real fluid_calc_div(fluid* f) {
    int i,j;
    real sum=0.0,div=0.0;
    for(j=0;j<f->ny;j++) {
        for(i=0;i<f->nx;i++) {
            div=
                (f->u1[i+1+(f->nx+1)*j]-f->u1[i+(f->nx+1)*j])/f->hx +
                (f->v1[i+f->nx*(j+1)]-f->v1[i+f->nx*j])/f->hy;
            f->div[i+f->nx*j] = div + f->fdiv[i+f->nx*j];
            sum += div;
        }
    }
    return div;
}

//Zero out the forced divergence field.
void fluid_zero_fdiv(fluid* f) {
    memset(f->fdiv,0,sizeof(real)*f->nx*f->ny);
}

//Project the velocity onto its divergence-free part.  Use the Helholtz theorem which says
//every vector field can be decomposed into a divergence-free and a curl-free part.  This is
//accomplished with a simple Poisson equation for the pressure.  
//Step 1: calculate divergence
//Step 2: calculate pressure
//Step 3: subtract pressure 
//Step 4: ???
//Step 5: Profit!
void fluid_project(fluid* f,int stepmax,int eastwall, int westwall, int northwall, int southwall) {
    int i,j,im,ip,jm,jp,step;
    real pim,pip,pjm,pjp;
    //Precalculate some coefficients for solving the pressure
    //equation since this is a good share of the processing time.
    real ci=f->hy*f->hy/(f->hx*f->hx+f->hy*f->hy)*0.5;
    real cj=f->hx*f->hx/(f->hx*f->hx+f->hy*f->hy)*0.5;
    real cd=-f->hx*f->hx*f->hy*f->hy/(f->hx*f->hx+f->hy*f->hy)*0.5;

    //This seems like a fine time and place to enforce some boundary conditions.
    //dv/dy=0 on south, north
    for(i=0;i<f->nx;i++) {
        f->v1[i] = f->v1[i+f->nx];
        f->v1[i+f->nx*(f->ny)] = f->v1[i+f->nx*(f->ny-1)];
    }
    //du/dx=0 on east,west
    for(j=0;j<f->ny;j++) {
        f->u1[(f->nx+1)*j] = f->u1[1+(f->nx+1)*j];
        f->u1[f->nx+(f->nx+1)*j]=f->u1[f->nx-1+(f->nx+1)*j];
    }
    //Set boundary velocities to zero or at least cause them to decay
    //very slightly so strong currents don't get started and keep going
    real Nfactor=0.99*northwall;
    real Sfactor=0.95*southwall;
    real Efactor=0.95*eastwall;
    real Wfactor=0.95*westwall;
    for(j=0;j<f->ny;j++) {
        f->u1[(f->nx+1)*j]*=Wfactor;
        f->v1[f->nx*j]*=Wfactor;
        f->v1[f->nx*(j+1)]*=Wfactor;
        f->T1[f->nx*j] = 0.0;
    }
    i=f->nx-1;
    for(j=0;j<f->ny;j++) {
        f->u1[i+1+(f->nx+1)*j]*=Efactor;
        f->v1[i+f->nx*j]*=Efactor;
        f->v1[i+f->nx*(j+1)]*=Efactor;
        f->T1[i+f->nx*j] = 0.0;
    }
    for(i=0;i<f->nx;i++) {
        f->v1[i]*=Sfactor;
        f->u1[i]*=Sfactor;
        f->u1[i+1]*=Sfactor;
        f->T1[i]=0.0;
    }
    j=f->ny-1;
    for(i=0;i<f->nx;i++) {
        f->v1[i+f->nx*(j+1)]*=Nfactor;
        f->u1[i+(f->nx+1)*j]*=Nfactor;
        f->u1[i+1+(f->nx+1)*j]*=Nfactor;
    }
    //Make the entire domain divergence-free by subtracting some
    //velocity from each of the four walls.  This has the effect
    //of an overall source/sink so that no net stuff is going in
    //or out.  It's not really a well-posed problem otherwise.
    //(But it's hard to tell the difference here, so I've skipped it.)
    if(0) {
        real totaldiv = fluid_calc_div(f);
        for(i=0;i<f->nx;i++) {
            f->v1[i] += totaldiv*f->hx*f->hy*0.125;
            f->v1[i+f->nx*f->ny] -= totaldiv*f->hx*f->hy*0.125;
        }
        for(j=0;j<f->ny;j++) {
            f->u1[(f->nx+1)*j] += totaldiv*f->hx*f->hy*0.125;
            f->u1[f->nx+(f->nx+1)*j] -= totaldiv*f->hx*f->hy*0.125;
        }
    }

    //Check to make sure it's divergence-free
    //fluid_calc_div(f);
    //totaldiv=0;
    //for(i=0;i<f->nx;i++) {
        //for(j=0;j<f->ny;j++) {
            //totaldiv+=f->div[i+f->nx*j];
        //}
    //}
    //fprintf(stderr,"total divergence = %15.10f\n",totaldiv);

    //Solve the pressure-poisson equation with SOR.  The overrelaxation
    //is about as high as it can going without being obvious.  Notice
    //there is no calculation of error since about twenty iterations
    //seems to be pretty good most of the time.
    real pnew, omega=1.87;
    for(step=0;step<stepmax;step++) {
        for(j=0;j<f->ny;j++) {
            for(i=0;i<f->nx;i++) {
                //Enforce dp/dn=0 while we're at it.
                im=(i-1<0)?(0):(i-1);
                ip=(i+1>f->nx-1)?(f->nx-1):(i+1);
                jm=(j-1<0)?(0):(j-1);
                jp=(j+1>f->ny-1)?(f->ny-1):(j+1);
                pim=f->p[im+f->nx*j];
                pip=f->p[ip+f->nx*j];
                pjm=f->p[i+f->nx*jm];
                pjp=f->p[i+f->nx*jp];
                //SOR G-S iteration
                pnew = ci*(pim+pip)+cj*(pjm+pjp)+cd*f->div[i+f->nx*j];
                f->p[i+f->nx*j] += (pnew-f->p[i+f->nx*j])*omega;
            }
        }
    }
    //Subtract the pressure gradient from velocity.  The result is a divergence-free
    //velocity field which goes back into the u and v arrays.
    i=f->nx;
    real oofx=1.0/f->hx;
    real oofy=1.0/f->hy;
    //Just copy the edges since the component of the pressure gradient is defined to be zero.
    for(j=0;j<f->ny;j++) {
        f->u[i+(f->nx+1)*j]=f->u1[i+(f->nx+1)*j];
    }
    for(j=0;j<f->ny;j++) {
        f->u[(f->nx+1)*j]=f->u1[(f->nx+1)*j];
    }
    for(j=0;j<f->ny;j++) {
        for(i=1;i<f->nx;i++) {
            f->u[i+(f->nx+1)*j]=f->u1[i+(f->nx+1)*j]-(f->p[i+f->nx*j]-f->p[i-1+f->nx*j])*oofx;
        }
    }
    for(i=0;i<f->nx;i++) {
        f->v[i]=f->v1[i];
    }
    j=f->ny;
    for(i=0;i<f->nx;i++) {
        f->v[i+f->nx*j]=f->v1[i+f->nx*j];
    }
    for(j=1;j<f->ny;j++) {
        for(i=0;i<f->nx;i++) {
            f->v[i+f->nx*j]=f->v1[i+f->nx*j]-(f->p[i+f->nx*j]-f->p[i+f->nx*(j-1)])*oofy;
        }
    }
    fluid_calc_vbar(f);
    fluid_calc_vorticity(f);
    for(i=0;i<f->nx*f->ny;i++) {
        f->T[i]=MIN(MAX(0,f->T1[i]),1.0);
        //f->T[i]=0.5+f->vort[i]*0.05;
    }
}

//Draw some OpenGL vectors.  Lines, really.
void fluid_draw(fluid* f,real scale) {
    int i,j;
    real ubar,vbar,x0,y0;
    glColor3f(1.0,1.0,1.0);
    glBegin( GL_LINES );
    for(j=0;j<f->ny;j++) {
        for(i=0;i<f->nx;i++) {
            ubar=0.5*(f->u[i+(f->nx+1)*j]+f->u[i+1+(f->nx+1)*j]);
            vbar=0.5*(f->v[i+f->nx*j]+f->v[i+f->nx*(j+1)]);
            x0=f->xmin+(i+0.5)*f->hx;
            y0=f->ymin+(j+0.5)*f->hy;
            glVertex2f(x0,y0);
            glVertex2f(x0+scale*ubar,y0+scale*vbar);
        }
    }
    glEnd();
}
