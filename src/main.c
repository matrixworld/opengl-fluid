#include "main.h"



//Lots of unorganized global variables.  Ugh.
int eastwall=1;
int westwall=1;
int northwall=1;
int southwall=0;
int cmode=0;
int rbuttondown=0;
int lbuttondown=0;
int mbuttondown=0;
int texallocd=0;
static int win;
int width; //of window
int height; //of window
int full=0; //fullscreen?
float mpx=0.0,mpx0=0.0; //variables for past and present window extents
float mpy=0.0,mpy0=0.0;
float xcen0,xcen=0.0;
float ycen0,ycen=0.0;
float ymin0,ymin,ymax0,ymax;
float xmin0,xmin,xmax0,xmax;
float dx0,dy0,dx,dy=2.1;
float mpx_new, mpy_new;
float mpx_old, mpy_old;

int buoyOn=1;
int vortOn=1;
int gni,gnj;
GLuint texture_id;
GLfloat *colormap;
GLfloat *Tcolors;
fluid* f;

void loadcolors(int type) {
    int c,i;
    FILE* colors;
    int r=0,g=0,b=0;
    switch(type) {
        case 0:
            c=0;
            colors = fopen("fire.dat","r");
            for(i=0;i<256;i++) {
                fscanf(colors,"%i %i %i",&r,&g,&b);
                colormap[c]=r/255.0;
                colormap[c+1]=g/255.0;
                colormap[c+2]=b/255.0;
                c+=3;
            }
            fclose(colors);
        break;
        case 1:
            c=0;
            for(i=0;i<256;i++) {
                colormap[c]=i/255.0;
                colormap[c+1]=i/255.0;
                colormap[c+2]=i/255.0;
                c+=3;
            }
        break;
    }
    
}
//Call this for each frame to turn the fluid into a texture mapped
//onto a single quad in the disp() function.
void specifyTexture() {
    if(texallocd) {
        glDeleteTextures(1,&texture_id);
        texallocd=0;
    }
    int i,j;
    memset(Tcolors,0,sizeof(float)*(gni+1)*(gnj+1));
    unsigned char colindex;
    int c=0;
    float t;
    if(cmode) {
        for(j=0;j<gnj;j++) {
            for(i=0;i<gni;i++) {
                colindex=(unsigned char)(255*f->T[i+f->nx*j]);
                Tcolors[c]=colormap[3*colindex];
                Tcolors[c+1]=colormap[3*colindex+1];
                Tcolors[c+2]=colormap[3*colindex+2];
                c+=3;
            }
        }
    } else {
        for(j=0;j<gnj;j++) {
            for(i=0;i<gni;i++) {
                //colindex=(unsigned char)(255*pow(fluid_interpolate_T(f,2.0*i/(gni-1)-1.0,2.0*j/(gnj-1)-1.0),0.25));
                t=1.0-f->T[i+f->nx*j];
                colindex=(unsigned char)(255*((1-t*t*t*t)+sqrt(1-t))*0.5);
                Tcolors[c]=colormap[3*colindex];
                Tcolors[c+1]=colormap[3*colindex+1];
                Tcolors[c+2]=colormap[3*colindex+2];
                c+=3;
            }
        }
    }
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, gni, gnj, 0, GL_RGB, GL_FLOAT, Tcolors);
    texallocd=1;
}

int main(int argc, char **argv) {
    int ni=256;
    int nj=256;
    f = alloc_fluid(ni,nj);
    f->xmin = -1.0;
    f->xmax = 1.0;
    f->ymin = -1.0;
    f->ymax = 1.0;
    f->dt = 0.1;
    fluid_calc_spacing(f);
    f->nu = 0.001;
    f->eps = 1.5e-6;
    f->buoy=1.0e-3;

    gni=256;
    gnj=256;
    colormap = malloc(sizeof(*colormap)*256*3);
    Tcolors = malloc(sizeof(*Tcolors)*(gni+1)*(gnj+1)*3);
    loadcolors(0);

    //INITIALIZE GLUT WITH SOME OPTIONS
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE );
    glutInitWindowSize(500,500);
    glutInitWindowPosition(100,100);
    win = glutCreateWindow("2D Window");

    //callback functions
    glutIdleFunc(disp);
    glutDisplayFunc(disp);    //called every time glut thinks a redraw is necessary
    glutReshapeFunc(reshape); //called when window is reshaped
    glutKeyboardFunc(keyb);   //called when keys pressed
    glutMouseFunc(mouseClick);//called on mouse click
    glutMotionFunc(mouseMove);//called when clicked mouse moves

    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial( GL_FRONT, GL_AMBIENT_AND_DIFFUSE );

    glGenTextures(1, &texture_id);
    specifyTexture();

    fprintf(stdout,"\n\n2D Fluid Demo\n\nby Ricky Reusser\n\n");
    fprintf(stdout,"Left click + drag to move fluid\n");
    fprintf(stdout,"Middle click to add fluid divergence\n");
    fprintf(stdout,"Right click to subtract divergence\n");
    fprintf(stdout,"Click edge of window to toggle walls\n");
    fprintf(stdout,"Window resizable (for the most part)\n");
    fprintf(stdout,"Keys:\n");
    fprintf(stdout,"\tb: buoyance on/off\n");
    fprintf(stdout,"\tv: vorticity confinement on/off\n");
    fprintf(stdout,"\t=/-: increase/decrease vorticity confinement force\n");
    fprintf(stdout,"\t]/[: increase/decrease buoyance force\n");
    fprintf(stdout,"\tz: zero out fluid temperature\n");
    fprintf(stdout,"\tspace: zero out velocities\n");
    fprintf(stdout,"\tC: image mode\n");
    fprintf(stdout,"\tc: color mode\n");
    fprintf(stdout,"\tf: toggle fullscreen\n");
    fprintf(stdout,"\tq: quit\n");
    fflush(stdout);

    glutMainLoop();
    //free_fluid(f);
    //free(colormap); 
    return 0;
}

void reshape(int w, int h) {
    width=w;
    height=h;

    ymin = ycen - dy/2.0; //recalculate the window extents
    ymax = ycen + dy/2.0;
    dx = dy*(float)w/h;
    xmin = xcen-dx/2.0;
    xmax = xcen+dx/2.0;
    f->xmin=xmin+0.1;
    f->xmax=xmax-0.1;
    fluid_calc_spacing(f);

    glViewport(0, 0, w, h); //tell OpenGL what kind of view to use
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(xmin,xmax,ymin,ymax); //2-D orthographic
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void disp(void) {
    //Do some stuff like add temperature for each frame if the mouse
    //button is still depressed.  Just a copy-paste from the mouse
    //callback functions
    if(lbuttondown) {
        int i,j,ii,ji;
        fluid_get_cell(f,mpx,mpy,&i,&j);
        int rad=5;
        int im=MIN(MAX(i-rad,0),f->nx-1);
        int ip=MIN(MAX(i+rad,0),f->nx-1);
        int jm=MIN(MAX(j-rad,0),f->ny-1);
        int jp=MIN(MAX(j+rad,0),f->ny-1);
        for(ji=jm;ji<jp;ji++) {
            for(ii=im;ii<ip;ii++) {
                f->T[ii+f->nx*ji] = MAX(f->T[ii+f->nx*ji],1.0-1.0*((i-ii)*(i-ii)+(j-ji)*(j-ji))/(float)rad/(float)rad);
            }
        }
    }
    if(mbuttondown) {
        int i,j,ii,ji;
        fluid_get_cell(f,mpx,mpy,&i,&j);
        int rad=8;
        int im=MIN(MAX(i-rad,0),f->nx-1);
        int ip=MIN(MAX(i+rad,0),f->nx-1);
        int jm=MIN(MAX(j-rad,0),f->ny-1);
        int jp=MIN(MAX(j+rad,0),f->ny-1);
        for(ji=jm;ji<jp;ji++) {
            for(ii=im;ii<ip;ii++) {
                f->T[ii+f->nx*ji] = MAX(f->T[ii+f->nx*ji],1.0-1.0*((i-ii)*(i-ii)+(j-ji)*(j-ji))/(float)rad/(float)rad);
            }
        }
        fluid_zero_fdiv(f);
        im=MIN(MAX(i-3,0),f->nx-1);
        ip=MIN(MAX(i+3,0),f->nx-1);
        jm=MIN(MAX(j-3,0),f->ny-1);
        jp=MIN(MAX(j+3,0),f->ny-1);
        for(ji=jm;ji<jp;ji++) {
            for(ii=im;ii<ip;ii++) {
                f->fdiv[ii+f->nx*ji] = -50.0;
            }
        }
    }
    if(rbuttondown) {
        int i,j,ii,ji;
        fluid_get_cell(f,mpx,mpy,&i,&j);
        int im=MIN(MAX(i-3,0),f->nx-1);
        int ip=MIN(MAX(i+3,0),f->nx-1);
        int jm=MIN(MAX(j-3,0),f->ny-1);
        int jp=MIN(MAX(j+3,0),f->ny-1);
        fluid_zero_fdiv(f);
        for(ji=jm;ji<jp;ji++) {
            for(ii=im;ii<ip;ii++) {
                f->fdiv[ii+f->nx*ji] = 100.0;
            }
        }
    }
    //Do the fluid calculations:
    if(vortOn) {
        fluid_calc_vbar(f);
        fluid_calc_vorticity(f);
    }
    if(vortOn||buoyOn) {
        fluid_calc_forces(f,vortOn,buoyOn);
        fluid_add_forces(f);
    }
    fluid_advect_u(f);
    fluid_advect_v(f);
    fluid_advect_T(f);
    fluid_calc_div(f);
    fluid_project(f,20,eastwall,westwall,northwall,southwall);
    //Done.


    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0, 0.0, -1.0);

    //Create a fluid texture
    specifyTexture();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture_id);

    //Map the texture onto a quad
    glBegin( GL_QUADS );
        glTexCoord2f(0.0,0.0); glVertex2f(f->xmin,f->ymin);
        glTexCoord2f(1.0,0.0); glVertex2f(f->xmax,f->ymin);
        glTexCoord2f(1.0,1.0); glVertex2f(f->xmax,f->ymax);
        glTexCoord2f(0.0,1.0); glVertex2f(f->xmin,f->ymax);
    glEnd();
    glDisable(GL_TEXTURE_2D);

    //Draw the boundaries
    glColor3f(1.0,1.0,1.0);
    glBegin( GL_LINES );
        if(!westwall)  { glVertex2f(f->xmin,f->ymin); glVertex2f(f->xmin,f->ymax); }
        if(!eastwall)  { glVertex2f(f->xmax,f->ymin); glVertex2f(f->xmax,f->ymax); }
        if(!northwall) { glVertex2f(f->xmin,f->ymax); glVertex2f(f->xmax,f->ymax); }
        if(!southwall) { glVertex2f(f->xmin,f->ymin); glVertex2f(f->xmax,f->ymin); }
    glEnd();

    glPopMatrix();
    glutSwapBuffers();
    glutPostRedisplay();
}



void keyb(unsigned char key, int x, int y) {
    switch(key) { //respond to keys
    case 'F': //full screen without window border
        glutReshapeWindow(1024,768);
        glutPositionWindow(0,-10);
        full=0;
    break;
    case 'f': //full screen with window border
        if(full==1) {
            glutReshapeWindow(500,500);
            glutPositionWindow(100,100);
        } else {
            glutFullScreen();
        }
        full=!full;
    break;
    case ' ':
        fluid_zero_velocity(f);
    break;
    case 'q': //quit
        free_fluid(f);
        free(colormap); 
        glutDestroyWindow(win);
        exit(0);
    break;
    case 'C':
        fluid_zero_velocity(f);
        buoyOn=0;
        cmode=1;
        fluid_read_image(f);
        loadcolors(1);
    break;
    case 'c':
        cmode=0;
        loadcolors(0);
    break;
    case 'z':
        fluid_zero_temperature(f);
    break;
    case 'v':
        vortOn=1-vortOn;
    break;
    case 'b':
        buoyOn=1-buoyOn;
    break;
    case ']':
        f->buoy *= 1.1;
    break;
    case '[':
        f->buoy /= 1.1;
    break;
    case '=':
        f->eps *= 1.1;
    break;
    case '-':
        f->eps /= 1.1;
    break;
    }
}

void mouseClick(int button, int state, int x, int y) {
    //function called only when mouse clicked or unclicked
    mpx = xmin+(xmax-xmin)*((float)x/width);
    mpy = ymin+(ymax-ymin)*(1.0-((float)y/height));
    mpx0 = mpx;   mpy0 = mpy;   xmin0 = xmin; ymin0 = ymin;
    xmax0 = xmax; ymax0 = ymax; xcen0 = xcen; ycen0 = ycen;
    dy0 = dy;
    mpx_old = mpx;
    mpy_old = mpy;

    if(mpx<f->xmin || mpx>f->xmax || mpy<f->ymin || mpy>f->ymax) {
        if(state==GLUT_DOWN) {
            if(mpx<f->xmin)
                westwall=1-westwall;
            if(mpx>f->xmax)
                eastwall=1-eastwall;
            if(mpy>f->ymax)
                northwall=1-northwall;
            if(mpy<f->ymin)
                southwall=1-southwall;
        }
    } else {
        int i,j,ii,ji,im,ip,jm,jp;
        fluid_get_cell(f,mpx,mpy,&i,&j);
        int rad;
        if(state==GLUT_DOWN) { //if button is down
            switch(button) {
            case GLUT_LEFT_BUTTON:
                lbuttondown=1;
                rad=5;
                im=MIN(MAX(i-rad,0),f->nx-1);
                ip=MIN(MAX(i+rad,0),f->nx-1);
                jm=MIN(MAX(j-rad,0),f->ny-1);
                jp=MIN(MAX(j+rad,0),f->ny-1);
                for(ji=jm;ji<jp;ji++) {
                    for(ii=im;ii<ip;ii++) {
                        f->T[ii+f->nx*ji] = MAX(f->T[ii+f->nx*ji],1.0-1.0*((i-ii)*(i-ii)+(j-ji)*(j-ji))/(float)rad/(float)rad);
                    }
                }
            break;
            case GLUT_MIDDLE_BUTTON:
                mbuttondown=1;
                rad=8;
                im=MIN(MAX(i-rad,0),f->nx-1);
                ip=MIN(MAX(i+rad,0),f->nx-1);
                jm=MIN(MAX(j-rad,0),f->ny-1);
                jp=MIN(MAX(j+rad,0),f->ny-1);
                for(ji=jm;ji<jp;ji++) {
                    for(ii=im;ii<ip;ii++) {
                        f->T[ii+f->nx*ji] = MAX(f->T[ii+f->nx*ji],1.0-1.0*((i-ii)*(i-ii)+(j-ji)*(j-ji))/(float)rad/(float)rad);
                    }
                }
                im=MIN(MAX(i-3,0),f->nx-1);
                ip=MIN(MAX(i+3,0),f->nx-1);
                jm=MIN(MAX(j-3,0),f->ny-1);
                jp=MIN(MAX(j+3,0),f->ny-1);
                for(ji=jm;ji<jp;ji++) {
                    for(ii=im;ii<ip;ii++) {
                        f->fdiv[ii+f->nx*ji] = -50.0;
                    }
                }
            break;
            case GLUT_RIGHT_BUTTON:
                rbuttondown=1;
                im=MIN(MAX(i-3,0),f->nx-1);
                ip=MIN(MAX(i+3,0),f->nx-1);
                jm=MIN(MAX(j-3,0),f->ny-1);
                jp=MIN(MAX(j+3,0),f->ny-1);
                for(ji=jm;ji<jp;ji++) {
                    for(ii=im;ii<ip;ii++) {
                        f->fdiv[ii+f->nx*ji] = 100.0;
                    }
                }
            break;
            }
        } else {
            fluid_zero_fdiv(f);
            rbuttondown=0;
            lbuttondown=0;
            mbuttondown=0;
        }
    }
}

void mouseMove(int x, int y) {
    mpx = xmin0+(xmax0-xmin0)*((float)x/width);
    mpy = ymin0+(ymax0-ymin0)*(1.0-((float)y/height));
    float xdiff=mpx-mpx_old;
    float ydiff=mpy-mpy_old;
    mpx_old=mpx;
    mpy_old=mpy;
    int i,j,ii,ji;
    if(lbuttondown|mbuttondown) {
        fluid_get_cell(f,mpx,mpy,&i,&j);
        int rad=5;
        int im=MIN(MAX(i-rad,0),f->nx-1);
        int ip=MIN(MAX(i+rad,0),f->nx-1);
        int jm=MIN(MAX(j-rad,0),f->ny-1);
        int jp=MIN(MAX(j+rad,0),f->ny-1);
        for(ji=jm;ji<jp;ji++) {
            for(ii=im;ii<ip;ii++) {
                f->u[ii+(f->nx+1)*ji] += xdiff/f->dt;
                f->v[ii+f->nx*ji] += ydiff/f->dt;
                f->T[ii+f->nx*ji] = MAX(f->T[ii+f->nx*ji],1.0-1.0*((i-ii)*(i-ii)+(j-ji)*(j-ji))/(float)rad/(float)rad);
            }
        }
    }
}
