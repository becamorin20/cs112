

#ifdef WIN32
#include <windows.h>

#endif

#include <math.h>
#include <GLUT/GLUT.h>
#include <iostream>
using namespace std;
// Defines maximum stack capacity.
#define STACK_CAP 16
// Defines pi for converting angles.
#define PI 3.14159265359

// Structure for the matrix stack, top specifies current top position on the stack, initially zero (which means one matrix in the stack)
struct matrix_stack
{
    GLdouble m[STACK_CAP][16];
    int top;
};

// Define a macro for retrieving current matrix from top of current stack.
#define current_matrix (current_stack->m[current_stack->top])

// Identity matrix constant.
const GLdouble identity[16] =
{1, 0, 0, 0,
 0, 1, 0, 0,
 0, 0, 1, 0,
 0, 0, 0, 1};

// the model view matrix stack.
struct matrix_stack model_view = {{{0}}, 0};
// the projection matrix stack.
struct matrix_stack projection = {{{0}}, 0};
// The current stack pointer that specifies the matrix mode.
struct matrix_stack *current_stack = &model_view;

// Multiplies current matrix with matrix b, put the result in current matrix.
// current = current * b
void matrix_multiply(const GLdouble *b)
{
    GLdouble *current = current_matrix;
    //current matrix defined...now multiply them
    //first row
    
    current[0] = (current[0] * b[0]) + (current[1] * b[4]) + (current[2] * b[8]) + (current[3] * b[12]);
    current[1] = (current[0] * b[1]) + (current[1] * b[5]) + (current[2] * b[9]) + (current[3] * b[13]);
    current[2] = (current[0] * b[2]) + (current[1] * b[6]) + (current[2] * b[10]) + (current[3] * b[14]);
    current[3] = (current[0] * b[3]) + (current[1] * b[7]) + (current[2] * b[11]) + (current[3] * b[15]);
    //second row
    current[4] = (current[4] * b[0]) + (current[5] * b[4]) + (current[6] * b[8]) + (current[7] * b[12]);
    current[5] = (current[4] * b[1]) + (current[5] * b[5]) + (current[6] * b[9]) + (current[7] * b[13]);
    current[6] = (current[4] * b[2]) + (current[5] * b[6]) + (current[6] * b[10]) + (current[7] * b[14]);
    current[7] = (current[4] * b[3]) + (current[5] * b[7]) + (current[6] * b[11]) + (current[7] * b[15]);
    //third row
    current[8] = (current[8] * b[0]) + (current[9] * b[4]) + (current[10] * b[8]) + (current[11] * b[12]);
    current[9] = (current[8] * b[1]) + (current[9] * b[5]) + (current[10] * b[9]) + (current[11] * b[13]);
    current[10] = (current[8] * b[2]) + (current[9] * b[6]) + (current[10] * b[10]) + (current[11] * b[14]);
    current[11] = (current[8] * b[3]) + (current[9] * b[7]) + (current[10] * b[11]) + (current[11] * b[15]);
    //fourth row
    current[12] = (current[12] * b[0]) + (current[13] * b[4]) + (current[14] * b[8]) + (current[15] * b[12]);
    current[13] = (current[12] * b[1]) + (current[13] * b[5]) + (current[14] * b[9]) + (current[15] * b[13]);
    current[14] = (current[12] * b[2]) + (current[13] * b[6]) + (current[14] * b[10]) + (current[15] * b[14]);
    current[15] = (current[12] * b[3]) + (current[13] * b[7]) + (current[14] * b[11]) + (current[15] * b[15]);
    
}

// Calculates cross product of b and c, put the result in a
// a = b x c
void cross_product(GLdouble *ax, GLdouble *ay, GLdouble *az,
    GLdouble bx, GLdouble by, GLdouble bz,
    GLdouble cx, GLdouble cy, GLdouble cz)
{
    *ax = (by * cz) - (bz * cy);
    *ay = (bz * cx) - (bx * cz);
    *az = (bx * cy) - (by * cx);
}

// Normalizes vector (x, y, z).
void normalize(GLdouble *x, GLdouble *y, GLdouble *z)
{
    
    GLdouble length = (GLdouble) sqrt((((*x)*(*x)) + ((*y)*(*y)) + ((*z)*(*z))));
    *x = *x/length;
    *y = *y/length;
    *z = *z/length;
}

// Switches matrix mode by changing the current stack pointer.
void I_my_glMatrixMode(GLenum mode)
{
    //if input mode is model view
    if (mode == GL_MODELVIEW)
    {
        current_stack = &model_view;
        return;
    }
    //if input mode is projection view
    else if (mode == GL_PROJECTION)
    {
        current_stack = &projection;
        return;
    }
}

// Overwrites current matrix with identity matrix.
void I_my_glLoadIdentity(void)
{
    GLdouble *current = current_matrix;
    for (int i = 0; i< STACK_CAP; i++)
    {
        current[i] = identity[i];
    }
}

// Pushes current matrix onto current stack, report error if the stack is already full.
void I_my_glPushMatrix(void)
{
    //get current matrix to push
    GLdouble *current = current_matrix;
    //push the matrix to the top of the matrix
    
    if (current_stack == &model_view)
    {
        if (model_view.top < STACK_CAP) //Stack cap = 16
        {
            for (int i = 0; i < STACK_CAP; i++)
            {
                //copy each value in current matrix into the model view matrix
                model_view.m[model_view.top + 1][i] = current[i];
            }
            model_view.top++;
        }
        else
        {
            cout << "Model View Stack Limit reached" << endl;
        }
    }
    else if (current_stack == &projection)
    {
        if (projection.top < STACK_CAP)
        {
            for (int i = 0; i < STACK_CAP; i++)
            {
                projection.m[projection.top + 1][i] = current[i];
            }
            projection.top++;
        }
        else
        {
            cout << "Projection View Stack limit reached"<< endl;
        }
    
    }
}

// Pops current matrix from current stack, report error if the stack has only one matrix left.
void I_my_glPopMatrix(void)
{
    if (current_stack == &model_view)
    {
        if (model_view.top != 0)//last value at 0
        {
            //remove model_view.top
            model_view.top = model_view.top - 1;
        }
        else
        {
            cout << "Model view stack has only 1 matrix left" << endl;
        }
    }
    else if (current_stack == &projection)
    {
        if (projection.top != 0)//last value is 0
        {
            //remove projection.top
            projection.top = projection.top - 1;
        }
        else
        {
            cout << "Projection stack has only 1 matrix left" << endl;
        }
    }
}

// Overwrites currentmatrix with m.
void I_my_glLoadMatrixf(const GLfloat *m)
{
    GLdouble *current = current_matrix;
    for (int i = 0; i < STACK_CAP; i++)
    {
        current[i] = m[i];
    }
}

void I_my_glLoadMatrixd(const GLdouble *m)
{
    GLdouble *current = current_matrix;
    for(int i = 0; i < STACK_CAP; i++)
    {
        current[i] = m[i];
    }
}

void I_my_glTranslated(GLdouble x, GLdouble y, GLdouble z)
{
    const GLdouble translating[16] =
    {1, 0, 0, x,
     0, 1, 0, y,
     0, 0, 1, z,
     0, 0, 0, 1};
    //transpose the translating matrix
    GLdouble transpose[STACK_CAP];
    for (int i = 0; i < STACK_CAP/4; i++)
    {    for (int j = 0; j < STACK_CAP/4; j++)
        {
            transpose[j*4+i] = translating[i*4+j];
        }
    }
    
    //multiply current with transpose to get the result
    matrix_multiply(transpose);
    
}

void I_my_glTranslatef(GLfloat x, GLfloat y, GLfloat z)
{
    I_my_glTranslated((GLdouble)x, (GLdouble)y, (GLdouble)z);
}

// Remember to normalize vector (x, y, z), and to convert angle from degree to radius before calling sin and cos.
void I_my_glRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z)
{
    // normalize vector x y z
    normalize(&x, &y, &z);
    
    //degree = radians * 180/pi
    GLdouble radian = angle * ((GLdouble)PI/ (GLdouble)180);
    GLdouble c = (GLdouble)cos(radian);
    GLdouble s = (GLdouble)sin(radian);
    GLdouble xx = x*x;
    GLdouble xy = x*y;
    GLdouble xz = x*z;
    GLdouble yy = y*y;
    GLdouble yz = y*z;
    GLdouble zz = z*z;
    
    
 
    //make a temp change
    GLdouble result_rotation[STACK_CAP];
    result_rotation[0] = (xx * (1 - c)) + c;
    result_rotation[1] = (xy * (1 - c)) - (z*s);
    result_rotation[2] = (xz * (1 - c)) + (y*s);
    result_rotation[3] = 0;
    
    result_rotation[4] = (xy * (1-c)) + (z*s);
    result_rotation[5] = (yy * (1-c)) + c;
    result_rotation[6] = (yz * (1-c)) - (x*s);
    result_rotation[7] = 0;
    
    result_rotation[8] = (xz * (1-c)) - (y*s);
    result_rotation[9] = (yz * (1-c)) + (x*s);
    result_rotation[10] = (zz * (1-c)) + c;
    result_rotation[11] = 0;
    
    result_rotation[12] = 0;
    result_rotation[13] = 0;
    result_rotation[14] = 0;
    result_rotation[15] = 1;
    
    
    //transpose the matrix
    GLdouble transpose[STACK_CAP];
    for (int i = 0; i < STACK_CAP/4; i++)
    {    for (int j = 0; j < STACK_CAP/4; j++)
        {
            transpose[j*4+i] = result_rotation[i*4+j];
        }
    }
    //change to this matrix
    matrix_multiply(transpose);
    
}

void I_my_glRotatef(GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
    I_my_glRotated((GLdouble)angle, (GLdouble)x, (GLdouble)y, (GLdouble)z);
}

void I_my_glScaled(GLdouble x, GLdouble y, GLdouble z)
{
    const GLdouble scale_matrix[16] =
    {   x,0,0,0,
        0,y,0,0,
        0,0,z,0,
        0,0,0,1};
    //need transpose?
    GLdouble transpose[STACK_CAP];
    for (int i = 0; i < STACK_CAP/4; i++)
    {
        for (int j = 0; j < STACK_CAP/4; j++)
        {
            transpose[j+4+i] = scale_matrix[i*4+j];
        }
    }
    matrix_multiply(transpose);
    
}

void I_my_glScalef(GLfloat x, GLfloat y, GLfloat z)
{
    I_my_glScaled((GLdouble)x, (GLdouble)y, (GLdouble)z);
}

// Copes current matrix to m.
void I_my_glGetMatrixf(GLfloat *m)
{
    GLdouble *current = current_matrix;
    
    if (current_stack == &model_view)
    {
        for(int i = 0; i < STACK_CAP; i++)
        {
            m[i] = (GLfloat) current[i];
        }
    }
    else if (current_stack == &projection)
    {
        for(int i = 0; i < STACK_CAP; i++)
        {
            m[i] = (GLfloat)current[i];
        }
    }
}

void I_my_glGetMatrixd(GLdouble *m)
{
    GLdouble *current = current_matrix;
    
    if (current_stack == &model_view)
    {
        for(int i = 0; i < STACK_CAP; i++)
        {
            m[i] = current[i];
        }
    }
    else if (current_stack == &projection)
    {
        for(int i = 0; i < STACK_CAP; i++)
        {
            m[i] = current[i];
        }
    }

}

// Remember to normalize vectors.
void I_my_gluLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ, 
    GLdouble centerX, GLdouble centerY, GLdouble centerZ, 
    GLdouble upX, GLdouble upY, GLdouble upZ)
{

    GLdouble F1 = centerX - eyeX;
    GLdouble F2 = centerY - eyeY;
    GLdouble F3 = centerZ - eyeZ;
    normalize(&F1, &F2, &F3);
    //normalize(&upX, &upY, &upZ);
    GLdouble s1,s2,s3;
    cross_product(&s1,&s2,&s3,F1,F2,F3,upX,upY,upZ);
    normalize(&s1, &s2, &s3);
    GLdouble u1,u2,u3;
    cross_product(&u1,&u2,&u3, s1,s2,s3,F1,F2,F3);
    normalize(&u1, &u2, &u3);
    //cross product to create f
    //GLdouble f1, f2, f3;//normalized
    //cross_product(&f1, &f2, &f3, F1, F2, F3, F1, F2, F3);
    
    //cross product to create UP''
    //GLdouble up1,up2,up3;//normalized
    //cross_product(&up1,&up2,&up3, upX,upY,upZ,upX,upY,upZ);
    
    
    //cross_product(&s1,&s2,&s3, f1,f2,f3,up1,up2,up3);
    
    //GLdouble sx,sy,sz;
    //cross_product(&sx,&sy,&sz, s1,s2,s3,s1,s2,s3);
    
    
    //cross_product(&u1,&u2,&u3,sx,sy,sz,f1,f2,f3);
    
    GLdouble M[16] =
    {
        s1,s2,s3,
        u1,u2,u3,
        -F1,-F2,-F3,
        0,0,0,1 };
    matrix_multiply(M);
    I_my_glTranslated(-eyeX, -eyeY, -eyeZ);
    
}

void I_my_glFrustum(GLdouble left, GLdouble right, GLdouble bottom,
    GLdouble top, GLdouble zNear, GLdouble zFar)
{
    GLdouble A = (right + left) / (right - left);
    GLdouble B = (top + bottom) / (bottom - top);
    GLdouble C = -(zFar + zNear) / (zFar - zNear);
    GLdouble D = (-2 * zFar* zNear) / (zFar - zNear);
    
    
    GLdouble Frust[16] =
    {
        (2*zNear/(right - left)), 0, A, 0,
        0, ((2 *zNear)/(bottom - top)), B, 0,
        0, 0, C, D,
        0, 0, -1, 0 };
    
    matrix_multiply(Frust);
    
    
}

// Based on the inputs, calculate left, right, bottom, top, and call I_my_glFrustum accordingly

// remember to convert fovy from degree to radius before calling tan
void I_my_gluPerspective(GLdouble fovy, GLdouble aspect, 
    GLdouble zNear, GLdouble zFar)
{
    
    GLdouble f = tan(fovy/((GLdouble)PI * (GLdouble) 360));
    GLdouble top, right;
    top = zNear * f;
    right = top * aspect;
    
    I_my_glFrustum(-right, right, -top, top, zNear, zFar);
    /*
    GLdouble generated[16] =
    {
        (f/aspect), 0, 0, 0,
        0 ,f ,0 ,0,
        0,0,((zFar + zNear)/(zNear - zFar)),((2 * zFar * zNear)/(zNear - zFar)),
        0 ,0 ,-1 ,0
    };
     */
    
    
}
