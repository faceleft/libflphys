A simple library for calculating multi-body motion that can use air resistance, wind, mutual gravity, total gravity

# Installing 

```bash
make install
```
This can be used without installation, you just need to copy the `flphys.c` and `flphys.h` files into the desired project,
only the standard library (including `math.h`) has dependencies. The `-lm` compiler/linker flag may be required to use it.

### Uninstall 
```bash
make uninstall
```

# Test
Tests the functions of the library
```bash
make test
```
# Bench
Runs simulations with different numbers of objects, with 1 ms accuracy, and outputs the ratio of simulation time to real time.
```bash
make bench
```



# Provides

## Functions

### pvec_scs_create()

Alternatively, you can create a vector by its length and angles XOY and ZOY (radians), for two-dimensional simulations the angle ZOY must be equal to PI/2. Returns initialized vector.

```C
pvec_t pvec_scs_create(double len, double xy_angle, double zy_angle);
```

### pvec_len()

Gets the vector and returns the calculated length

```C
double pvec_len(pvec_t vector);
```

### pvec_xy_angle()

Gets the vector and returns the calculated XOY angle in radians

```C
double pvec_xy_angle(pvec_t vector);
```

### pvec_zy_angle()

Gets the vector and returns the calculated ZOY angle in radians

```C
double pvec_zy_angle(pvec_t vector);
```

### pobj_create()
Gets the required data, computes the others, and returns the initialized object.
* `pos` - initial position of the object 
* `mov` - initial motion of the object
* `mass` - object mass
* `radius` - object radius, all objects are spherical

```C
pobj_t pobj_create(pvec_t pos, pvec_t mov, double mass, double radius);
```

### pobj_run()
Сalculates the motion of the passed object along the parabola by `mov` and `force` fields, for `time` seconds. Returns `pres_t` structure (see below).  Usually used only inside the library
```C
pres_t pobj_run(pobj_t * obj, double time);
```
### pobj_set_radius()
Sets the radius of the object and recalculates its cross-sectional area and volume.
```C
void pobj_set_radius(pobj_t * obj, double radius);
```
### pobj_set_area()
Sets the cross-sectional area of the object and recalculates its radius and volume.
```C
void pobj_set_area(pobj_t * obj, double area);
```
### pobj_set_volume()
Sets the volume of the object and recalculates its cross-sectional area and radius.
```C
void pobj_set_volume(pobj_t * obj, double volume);
```
### phys_create()
Gets the required data, computes the others, and returns the initialized physics scene.
* `double density` - ambient air density, set it to 0 if there's a vacuum around. You can use the constant PHYS_AIR_DENSITY
* `accel_of_gravity` - Free fall acceleration of the environment. Set it to 0 if there is no common gravity. You can use the constant PHYS_ACCEL_OF_FREE_FALL.
* `wind` - air movement vector
* `objects[]` - pointer to objects array
* `objects_num` - number of objects in the array
* `is_gravity` - flag enabling gravity between objects, works slowly, *O(n^2)* where *n* is the number of objects.
```C
phys_t phys_create(double density,
                   pvec_t accel_of_gravity,
                   pvec_t wind,
                   pobj_t objects[],
                   size_t objects_num,
                   bool is_gravity);
```

### phys_run()
Basic function for calculations. Calculates the movement by steps. 
* `steps` - number of steps
* `step_time` - stepping time

During a step, all objects move at an equal acceleration, after which the forces are recomputed. Smaller step times improve the accuracy of the computation, but require more repeats, usually 1us - 1ms is sufficient.
```C
pres_t phys_run(phys_t * phys, double step_time, uint64_t steps);
```

## Constants
**Constants used for calculations**

Gravitational constant:
```C
extern const double PHYS_G;                   //6.6743015151515151514e-11 m^3/(kg*s^2)
```

Number of pi:
```C
extern const double PHYS_PI;                  //3.1415926535897932385
```

Normal air density:
```C
extern const double PHYS_AIR_DENSITY;         //1.225 kg/m^3
```

The normal acceleration of free fall:
```C
extern const double PHYS_ACCEL_OF_FREE_FALL;  //9.80665 m/s^2
```

Approximate shape resistance coefficient for a sphere:
```C
extern const double PHYS_BALL_DRAG_COEF;      //0.47
```

## Structures

### Vector

A three-dimensional vector that is used in all vector quantities, it is better to store values in SI units, you can create an instance of it: `(pvec-t){x, y, z}`

```C
typedef struct {
   double x; // m | m/s | N
   double y; // m | m/s | N
   double z; // m | m/s | N
} pvec_t;
```
### Object

The structure of an object that contains its position, velocity, mass, radius, cross-sectional area, volume, and the force acting on it. The object is considered to be a ball, the cross-sectional area and volume are calculated at creation using the appropriate formulas. The field `force` is needed as a buffer for calculations.

```C
typedef struct {
    pvec_t pos;    //m
    pvec_t mov;    //m/s
    double mass;   //kg
    double radius; //m
    double area;   //m^2
    double volume; //m^3
    pvec_t force;  //N
} pobj_t;
```
### Scene

The scene uniting the objects, through which all calculations are performed. Stores environment parameters, pointer to objects and their number. Also in `time` stores the total simulation time in seconds.

```C
typedef struct {
    double density;           //ambient air density kg/m^3
    pvec_t accel_of_gravity;  //constant acceleration acting on all objects m/s^2
    pvec_t wind;              //ambien wind m/s
    bool is_gravity;          //inter-object gravity flag
    double time;              //total simulation time
    pobj_t * objects;         //pointer to objects array
    size_t objects_num;       //number of objects in array
} phys_t;
```
### Result

The enumeration returned by the motion calculation functions
* `OK` - successfully
* `ERR_NULL_PTR` - pointer to objects is null when their number is greater than 0
* `ERR_ZERO_DIST` - 2 objects are at the same point, with gravity turned on
* `ERR_ZERO_MASS` - massless object
```C
typedef enum {
    OK = 0,
    ERR_NULL_PTR, 
    ERR_ZERO_DIST,
    ERR_ZERO_MASS
} pres_t;
```

## Example
```C
#include <stdio.h>
#include <flphys.h>

#define BALL_SPEED_MS 16.67 //60 km/h
#define BALL_ANGLE 45.0 //degrees

int main() {
    pobj_t volleyball = pobj_create(
        (pvec_t){0,1,0}, 
        pvec_scs_create(BALL_SPEED_MS, BALL_ANGLE*(PHYS_PI/180), PHYS_PI/2),
        0.27, 
        0.105
    );

    /*
    A volleyball at the point {0,1,0} (let's say one meter 
    above the ground) flying at a velocity of 16.67 m/s, 
    at an angle of 45 degrees to the horizon. 
    It has a mass of 270 grams and a radius of 10.5 cm.
    */

    phys_t scene = phys_create(
        PHYS_AIR_DENSITY, 
        (pvec_t){0, -PHYS_ACCEL_OF_FREE_FALL, 0}, 
        (pvec_t){0}, 
        &volleyball, 
        1, 
        false
    );

    /*
    Environment with air density PHYS_AIR_DENSITY (1.225 kg/m^3), 
    with constant gravity in PHYS_ACCEL_OF_FREE_FALL to negative 
    semi-major axis OY (9.8 m/s^2), no wind, with 1 'volleyball' 
    object, no inter-object gravity
    */
    
    puts("time, [pos.x, pos.y], speed");
    while(volleyball.pos.y >= 0) { //while above ground

        printf("%.2f, [%.2f, %.2f], %.2f\n", 
            scene.time, 
            volleyball.pos.x, 
            volleyball.pos.y, 
            pvec_len(volleyball.mov)
        );    
        phys_run(&scene, 0.0001, 1000);
    }
    
    /*
    While the ball is in the air, the time, position and speed of 
    the ball are displayed every tenth of a second of its flight.
    */
   
    return 0;
}
```