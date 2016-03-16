#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double Matrix[3][3];
typedef double Vector[3];

void getRotXMatrix(double angle, Matrix c){
    c[0][0] = 1.;
    c[0][1] = c[0][2] = c[1][0] = c[2][0] = 0.;
    c[1][1] = c[2][2] = cos(angle);
    c[2][1] = sin(angle);
    c[1][2] = -c[2][1];
}

void getRotYMatrix(double angle, Matrix c){
    c[1][1] = 1.;
    c[1][0] = c[1][2] = c[0][1] = c[2][1] = 0.;
    c[0][0] = c[2][2] = cos(angle);
    c[0][2] = sin(angle);
    c[2][0] = -c[0][2];
}

void getRotZMatrix(double angle, Matrix c){
    c[2][2] = 1.;
    c[2][0] = c[2][1] = c[0][2] = c[1][2] = 0.;
    c[0][0] = c[1][1] = cos(angle);
    c[1][0] = sin(angle);
    c[0][1] = -c[1][0];
}

void matrixAdd(const Matrix a, const Matrix b, Matrix c) {
    for (int y=0;y<3;y++)
        for (int x=0;x<3;x++)
            c[y][x] = a[y][x] + b[y][x];
}

void matrixMult(const Matrix a, const Matrix b, Matrix c) {
    for (int y=0;y<3;y++)
        for (int x=0;x<3;x++) {
            c[y][x] = 0.;
            for (int i=0;i<3;i++)
                c[y][x] += a[y][i] * b[i][x];
        }
}

void vectorAdd(const Vector a, const Vector b, Vector c) {
    for (int x=0;x<3;x++) c[x]=a[x] + b[x];
}

void vectorSub(const Vector a, const Vector b, Vector c) {
    for (int x=0;x<3;x++) c[x]=a[x] - b[x];
}

double vectorLen(const Vector a) {
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

double vectorDist(const Vector a, const Vector b) {
    Vector d;
    vectorSub(a,b,d);
    return vectorLen(d);
}

void vectorSMult(const double a, const Vector b, Vector c) {
    for (int x=0;x<3;x++) c[x]=a * b[x];
}

void vectorMult(const Vector a, const Vector b, double &c) {
    c=0;
    for (int x=0;x<3;x++) c += a[x] * b[x];
}

void vectorCross(const Vector a, const Vector b, Vector c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void matrixVectorMult(const Matrix a, const Vector b, Vector c) {
    for (int y=0;y<3;y++) {
        c[y] = 0.;
        for (int i=0;i<3;i++)
            c[y] += a[y][i] * b[i];
    }
}

void cartesianToSpherical(const double x, const double y, const double z, double &r, double &theta, double &phi) {
    r = sqrt(x*x+y*y+z*z);
    theta = asin(z/r);
    phi = atan2(y,x);
}

void sphericalToCartesian(const double r, const double theta, const double phi, double &x, double &y, double &z) {
    x = r * cos(theta) * cos(phi);
    y = r * cos(theta) * sin(phi);
    z = r * sin(theta);
}

void tripleTowerToCenter(const Vector a, const Vector b, const Vector c, const double la, const double lb, const double lc, Vector center)
{
    Vector oa, ob, oc;
    // calculate offset vectors to Vector a
    vectorSub(a,a,oa);
    vectorSub(b,a,ob);
    vectorSub(c,a,oc);
    // Rotate Vector b to positive x-Axis
    double r, theta, phi;
    cartesianToSpherical(ob[0],ob[1],ob[2],r,theta,phi);
    Matrix rotY, rotZ, rot, rotb;
    getRotYMatrix(theta,rotY);
    getRotZMatrix(-phi,rotZ);
    matrixMult(rotY, rotZ, rot);
    Vector rb, rc;
    matrixVectorMult(rot, ob, rb);
    matrixVectorMult(rot, oc, rc);
    fprintf (stdout, "(%f|%f|%f)(%f|%f|%f)\n" ,ob[0],ob[1], ob[2],r, theta, phi);
    fprintf (stdout, "(%f|%f|%f)\n" ,rb[0],rb[1], rb[2]);
    fprintf (stdout, "(%f|%f|%f)\n" ,rc[0],rc[1], rc[2]);
    // Rotate around x-Axis in a way that Vector c lies on the x-y-Plane
    double omega = atan2(rc[2],rc[1]);
    Matrix rotX;
    getRotXMatrix(-omega, rotX);
    Vector fc;
    matrixVectorMult(rotX, rc, fc);
    fprintf (stdout, "(%f|%f|%f)(%f)\n" ,fc[0],fc[1], fc[2],omega);
    // Calculate Center point
    Vector ot;
    ot[0]=(la*la-lb*lb+rb[0]*rb[0])/(2*rb[0]);
    ot[1]=(la*la-lc*lc+fc[0]*fc[0]+fc[1]*fc[1]-2*fc[0]*ot[0])/(2*fc[1]);
    ot[2]=-sqrt(la*la-ot[0]*ot[0]-ot[1]*ot[1]);
    fprintf(stdout, "(%f|%f|%f)\n",ot[0], ot[1], ot[2]);
    // Rotate Center around x-Axis
    // Rotate Center back using Spherical coordinate angles
    Vector os;
    getRotXMatrix(omega, rotX);
    getRotYMatrix(-theta,rotY);
    getRotZMatrix(phi,rotZ);
    matrixMult(rotZ, rotY, rot);
    matrixMult(rot, rotX, rotb);
    matrixVectorMult(rotb, ot, os);
    fprintf(stdout, "(%f|%f|%f) %f=%f %f=%f %f=%f\n",os[0], os[1], os[2], la, vectorLen(os), lb, vectorDist(os,ob), lc, vectorDist(os,oc));
    // Move to the right offset
    vectorAdd(os, a, center);
}

void CenterToTripleTower()
{

}
#endif // GEOMETRY_H
