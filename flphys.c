#include "flphys.h"
#include <math.h>

const double PHYS_G = 6.6743015151515151514e-11;
const double PHYS_PI = 3.1415926535897932385;
const double PHYS_AIR_DENSITY = 1.225;
const double PHYS_ACCEL_OF_FREE_FALL = 9.80665;
const double PHYS_BALL_DRAG_COEF = 0.47;

pvec_t pvec_scs_create(double len, double xy_angle) {
    return (pvec_t) {
        .x = len * cos(xy_angle),
        .y = len * sin(xy_angle),
    };
}

double pvec_len(pvec_t vector) {
    return sqrt(vector.x*vector.x + vector.y*vector.y);
}

double pvec_xy_angle(pvec_t vector) {
    if(vector.x == 0) return PHYS_PI/2.0;
    return atan(vector.y/vector.x);
}

pobj_t pobj_create(pvec_t pos, pvec_t mov, double mass, double radius) {
    return (pobj_t) {
        .pos = pos,
        .mov = mov, 
        .mass = mass, 
        .radius = radius, 
        .area = (PHYS_PI * radius * radius), 
        .volume = ((4.0L/3.0L) * PHYS_PI * radius * radius * radius),
        .force = (pvec_t) {0}
    };
}

pres_t pobj_run(pobj_t * obj, double time) {
    if(obj->mass == 0) {
        return ERR_ZERO_MASS;
    };
    pvec_t acceleration = {.x = obj->force.x/obj->mass,
                           .y = obj->force.y/obj->mass,};

    obj->pos.x += (obj->mov.x + acceleration.x * time * 0.5L) * time;
    obj->pos.y += (obj->mov.y + acceleration.y * time * 0.5L) * time;

    obj->mov.x += acceleration.x * time;
    obj->mov.y += acceleration.y * time;
    return OK;
}

void pobj_set_radius(pobj_t * obj, double radius) {
    obj->radius = radius;
    obj->area = (PHYS_PI * radius * radius);
    obj->volume = ((4.0/3.0) * PHYS_PI * radius * radius * radius);
}

void pobj_set_area(pobj_t * obj, double area) {
    obj->radius = sqrt(area) / PHYS_PI;
    obj->area = area;
    obj->volume = area * (4.0/3.0) * obj->radius;
}

void pobj_set_volume(pobj_t * obj, double volume) {
    obj->radius =  pow(volume, 1.0/3.0) * (3.0/4.0) / PHYS_PI;
    obj->area = PHYS_PI * obj->radius * obj->radius;
    obj->volume = volume;
}

phys_t phys_create(double density,
                   pvec_t acceleration_of_gravity,
                   pvec_t wind,
                   pobj_t objects[],
                   size_t objects_num,
                   bool is_gravity) {
    return (phys_t) {
        .density = density,
        .accel_of_gravity = acceleration_of_gravity,
        .wind = wind,
        .objects = objects,
        .objects_num = objects_num,
        .is_gravity = is_gravity,
        .time = 0
    };

}

static pres_t objects_force_compute(const phys_t * phys, pobj_t * obj) {                              
    obj->force = (pvec_t){0};

    double real_mass = obj->mass - (obj->volume * phys->density);

    pvec_t real_obj_mov = {
        obj->mov.x - phys->wind.x,
        obj->mov.y - phys->wind.y,      
    };

    double real_obj_speed = pvec_len(real_obj_mov);
    
    if (real_obj_speed != 0) {
        double air_f = obj->area 
                        * phys->density
                        * real_obj_speed
                        * real_obj_speed
                        * 0.5
                        * PHYS_BALL_DRAG_COEF;
    
    
        double k = air_f/real_obj_speed;
        obj->force.x -= real_obj_mov.x * k;
        obj->force.y -= real_obj_mov.y * k;
    }

    obj->force.x += phys->accel_of_gravity.x * real_mass;
    obj->force.y += phys->accel_of_gravity.y * real_mass;
    
    if (phys->is_gravity && phys->objects_num>1) {
        for(size_t i = 0; i<phys->objects_num; i++) {
            pobj_t * other_obj = phys->objects + i;
            if (other_obj != obj) {
                
                pvec_t dist = {.x = obj->pos.x - other_obj->pos.x,
                               .y = obj->pos.y - other_obj->pos.y,};

                double dist_len = pvec_len(dist);
                if(dist_len == 0) return ERR_ZERO_DIST;

                double gravy_f = real_mass
                                   * (other_obj->mass - (other_obj->volume * phys->density))
                                   * PHYS_G
                                   / (dist_len*dist_len);

                double k = gravy_f/dist_len;
                obj->force.x -= dist.x * k;
                obj->force.y -= dist.y * k;
            }
        }
    }
    return OK;
}

pres_t phys_run(phys_t * phys, double step_time, uint64_t steps) {
    if(phys->objects_num == 0) return OK;
    if(phys->objects == NULL) return ERR_NULL_PTR;
    
    uint64_t counter = 0;
    pres_t err;

    while(counter++ < steps) {
        for (size_t i = 0; i < phys->objects_num; i++) {
            err = objects_force_compute(phys, phys->objects+i);
            if(err) return err;
        }
        for (size_t i = 0; i < phys->objects_num; i++) {
            err = pobj_run(phys->objects+i, step_time);
            if(err) return err;
        }
    }
    phys->time+=step_time*steps;
    return OK;
}
