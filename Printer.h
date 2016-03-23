#ifndef PRINTER_H
#define PRINTER_H

#include "Geometry.h"

class Tower {
public:
    Tower() {}

private:
    Vector origin;          // Origin of coordinate system
    double theta,phi;       // direction perpendicular to z-Axis and the carriage. Defines the Effector plane.
    double zstop;           // distance from origin to endstop
    double cofs;            // offset between carriage rod connection point and z-Axis in normal direction
    double pos;             // position of carriage relative to origin of coordinate system
};

class Connector;

class Effector {
public:
    Effector() {}

private:
    Connector *con[3];
    double  conofs[3];
    Vector  nozofs;         // offset of the nozzle
};

class Connector {
public:
    Connector() {}

private:
    Tower       *tow;
    Effector    *eff;
    double      len;
};

class Platform {
public:
    Platform() {}

private:
    Vector  origin;
    double  dia;
    Vector  normal;

};

#endif // PRINTER_H
