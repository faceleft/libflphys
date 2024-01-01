#ifndef LIBFLPHYS_H
#define LIBFLPHYS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

extern const double PHYS_G;                   //6.6743015151515151514e-11
extern const double PHYS_PI;                  //3.1415926535897932385
extern const double PHYS_AIR_DENSITY;         //1.225
extern const double PHYS_ACCEL_OF_FREE_FALL;  //9.80665
extern const double PHYS_BALL_DRAG_COEF;      //0.47

typedef struct {
   double x; // m | m/s | N
   double y; // m | m/s | N
} pvec_t;

typedef struct {
    pvec_t pos;    //m
    pvec_t mov;    //m/s
    double mass;   //kg
    double radius; //m
    double area;   //m^2
    double volume; //m^3
    pvec_t force;  //N
} pobj_t;

typedef struct {
    double density;           //ambient air density kg/m^3
    pvec_t accel_of_gravity;  //constant acceleration acting on all objects m/s^2
    pvec_t wind;              //ambien wind m/s
    bool is_gravity;          //inter-object gravity flag
    double time;              //total simulation time
    pobj_t * objects;         //pointer to objects array
    size_t objects_num;       //number of objects in array
} phys_t;

//enumeration returned by phys_run() and pobj_run()
typedef enum {
    OK = 0,        //success
    ERR_NULL_PTR,  //pointer to objects is null when their number is greater than 0
    ERR_ZERO_DIST, //2 objects are at the same point, with gravity turned on
    ERR_ZERO_MASS  //massless object
} pres_t;

//returns a new initialized vector using length, XOY angle and ZOY angle
pvec_t pvec_scs_create(double len, double xy_angle);

//calculates the length of the vector
double pvec_len(pvec_t vector);

//calculates the XOY angle of the vector
double pvec_xy_angle(pvec_t vector);

//returns a new initialized vector using position, movement, mass and radius
pobj_t pobj_create(pvec_t pos, pvec_t mov, double mass, double radius);

//calculates the displacement of an object with obj->force
pres_t pobj_run(pobj_t * obj, double time);

//sets the radius and recalculates area and volume
void pobj_set_radius(pobj_t * obj, double radius);

//sets the area and recalculates radius and volume
void pobj_set_area(pobj_t * obj, double area);

//sets the volume and recalculates radius and area
void pobj_set_volume(pobj_t * obj, double volume);

//returns a new initialized scene
phys_t phys_create(double density,
                   pvec_t accel_of_gravity,
                   pvec_t wind,
                   pobj_t objects[],
                   size_t objects_num,
                   bool is_gravity);

//steps times calculates the strength for all objects, and calls pobj_run(phys->obj+i, step_time)
pres_t phys_run(phys_t * phys, double step_time, uint64_t steps);

#endif
