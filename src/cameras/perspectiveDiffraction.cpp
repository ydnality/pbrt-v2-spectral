
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// cameras/perspective.cpp*
#include "stdafx.h"
#include "cameras/perspectiveDiffraction.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include <gsl/gsl_randist.h>
#include <iostream>

// PerspectiveCamera Method Definitions
PerspectiveDiffractionCamera:: PerspectiveDiffractionCamera(const AnimatedTransform &cam2world,
        const float screenWindow[4], float sopen, float sclose,
        float lensr, float focald, float fov, Film *f)
    : ProjectiveCamera(cam2world, Perspective(fov, 1e-2f, 1000.f),
                       screenWindow, sopen, sclose, lensr, focald, f) {
    // Compute differential changes in origin for perspective camera rays
    dxCamera = RasterToCamera(Point(1,0,0)) - RasterToCamera(Point(0,0,0));
    dyCamera = RasterToCamera(Point(0,1,0)) - RasterToCamera(Point(0,0,0));


    //code borrowed from http://www.gnu.org/software/gsl/manual/gsl-ref.html#Random-Number-Generation  18.13
   const gsl_rng_type * T;
 
   int i, n = 10;
 
   gsl_rng_env_setup();
 
   T = gsl_rng_default;
   r = gsl_rng_alloc (T);

}


float PerspectiveDiffractionCamera::GenerateRay(const CameraSample &sample,
                                     Ray *ray) const {
    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    *ray = Ray(Point(0,0,0), Normalize(Vector(Pcamera)), 0.f, INFINITY);
    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft = focalDistance / ray->d.z;
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point(lensU, lensV, 0.f);
        ray->d = Normalize(Pfocus - ray->o);


        // --------------------add effect of diffraction----------------
        double initx = 0;
        double inity = 0;       
        double * x = &initx;
        double * y = &inity;   
        double currentAperture = lensRadius * 2;

        //calculate min distance direction
       
        //calculate radius
        Point intersectPoint(lensU, lensV, 0.f);
        double radius = sqrt(intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y );

//cout << "radius: " << radius << "\n";

//test for spectral
std::cout << "rayWavelength: " << ray->wavelength << "\n";

        //calculate direction 
        Vector direction(intersectPoint.x, intersectPoint.y, 0);
        
        Vector orthoDirection(-intersectPoint.y, intersectPoint.x, 0);  
                   // double direction = atan (intersectPoint.y/intersectPoint.x);
        
        double a = currentAperture/2 - radius;
        double b = sqrt((currentAperture/2 * currentAperture/2) - radius * radius);

        a = a;
        b = b;
        double pi = 3.14159265359;
        double lambda = .000000550;  //550 nanometers for now
        
        double sigma_x = atan(1/(2 * a *.001 * 2 * pi/lambda)); 
        double sigma_y = atan(1/(2 *b * .001 * 2 * pi/lambda)); 

        //gsl_ran_bivariate_gaussian (const gsl_rng * r, double sigma_x, double sigma_y, double rho, double * x, double * y)             
       gsl_ran_bivariate_gaussian (r, sigma_x, sigma_y, 0, x, y);    //experiment for now
               
        //add r.v. in directions of direction and orthoDirection


        //calculate component of these vectors based on 2 random degrees
        direction = Normalize(direction);
        orthoDirection = Normalize(orthoDirection);


        float noiseA = (float)(*x);
        float noiseB = (float)(*y);

        
         //project the original ray onto the new bases         
         double projA = (ray->d.x * direction.x + ray->d.y * direction.y)/sqrt(direction.x * direction.x + direction.y * direction.y);
         double projB = (ray->d.x * orthoDirection.x + ray->d.y * orthoDirection.y)/sqrt(orthoDirection.x * orthoDirection.x + orthoDirection.y * orthoDirection.y);
         double projC = ray->d.z;

         double rA = sqrt(projA * projA + projC * projC);
         double rB = sqrt(projB * projB + projC * projC);
         double thetaA = acos(projA/rA);          
         double thetaB = acos(projB/rB);

         //add uncertainty
         //thetaA = thetaA + noiseA;
        // thetaB = thetaB + noiseB;
         
         //convert angles back into cartesian coordinates, but in a,b space
         double newProjA = cos(thetaA) * rA;

         ray->d.z = sin(thetaA) * rA;

         projC = ray->d.z;
         rB = sqrt(projB * projB + projC * projC);
       
         //rB = sqrt(ray->d.y * ray->d.y + ray->d.z * ray->d.z); // THIS LINE IS WRONG!!!!! 

        //need to recalculate thetaB after you modify z (this is a new addition, also questionable)          
        thetaB = acos(projB/rB);

        double newProjB = cos(thetaB) * rB;
        ray->d.z = sin(thetaB) * rB;
    
        //convert from ab space back to x,y space
        ray->d.x = direction.x * newProjA + orthoDirection.x * newProjB;
        ray->d.y = direction.y * newProjA + orthoDirection.y * newProjB;

        ray->d = Normalize(ray->d);       


        //------------end diffraction-----------------


    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    return 1.f;
}


float PerspectiveDiffractionCamera::GenerateRayDifferential(const CameraSample &sample,
                                                 RayDifferential *ray) const {
    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    Vector dir = Normalize(Vector(Pcamera.x, Pcamera.y, Pcamera.z));
//std::cout << "rayWavelength: " << ray->wavelength << "\n";
    float tempWavelength = ray->wavelength;   //Andy added
    *ray = RayDifferential(Point(0,0,0), dir, 0.f, INFINITY);  //this gets rid of the wavelength field, so we need to add it back in
    ray->wavelength = tempWavelength;   //Andy added
//std::cout << "rayWavelength: " << ray->wavelength << "\n";

    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft = focalDistance / ray->d.z;
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point(lensU, lensV, 0.f);
        ray->d = Normalize(Pfocus - ray->o);

        

        // --------------------add effect of diffraction----------------
        double initx = 0;
        double inity = 0;       
        double * x = &initx;
        double * y = &inity;   
        double currentAperture = lensRadius * 2;

        //calculate min distance direction
       
        //calculate radius
        Point intersectPoint(lensU, lensV, 0.f);
        double radius = sqrt(intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y );

//cout << "radius: " << radius << "\n";
//std::cout << "rayWavelength: " << ray->wavelength << "\n";

        //calculate direction 
        Vector direction(intersectPoint.x, intersectPoint.y, 0);
        
        Vector orthoDirection(-intersectPoint.y, intersectPoint.x, 0);  
                   // double direction = atan (intersectPoint.y/intersectPoint.x);
        
        double a = currentAperture/2 - radius;
        double b = sqrt((currentAperture/2 * currentAperture/2) - radius * radius);

        a = a;
        b = b;
        double pi = 3.14159265359;
        //double lambda = .000000550;  //550 nanometers for now
        double lambda = ray->wavelength * 1e-9;
        double sigma_x = atan(1/(2 * a *.001 * 2 * pi/lambda)); 
        double sigma_y = atan(1/(2 *b * .001 * 2 * pi/lambda)); 

        //gsl_ran_bivariate_gaussian (const gsl_rng * r, double sigma_x, double sigma_y, double rho, double * x, double * y)             
       gsl_ran_bivariate_gaussian (r, sigma_x, sigma_y, 0, x, y);    //experiment for now
               
        //add r.v. in directions of direction and orthoDirection


        //calculate component of these vectors based on 2 random degrees
        direction = Normalize(direction);
        orthoDirection = Normalize(orthoDirection);


        float noiseA = (float)(*x);
        float noiseB = (float)(*y);

        
         //project the original ray onto the new bases         
         double projA = (ray->d.x * direction.x + ray->d.y * direction.y)/sqrt(direction.x * direction.x + direction.y * direction.y);
         double projB = (ray->d.x * orthoDirection.x + ray->d.y * orthoDirection.y)/sqrt(orthoDirection.x * orthoDirection.x + orthoDirection.y * orthoDirection.y);
         double projC = ray->d.z;

         double rA = sqrt(projA * projA + projC * projC);
         double rB = sqrt(projB * projB + projC * projC);
         double thetaA = acos(projA/rA);          
         double thetaB = acos(projB/rB);

         //add uncertainty
         thetaA = thetaA + noiseA;
         thetaB = thetaB + noiseB;
         
         //convert angles back into cartesian coordinates, but in a,b space
         double newProjA = cos(thetaA) * rA;

         ray->d.z = sin(thetaA) * rA;

         projC = ray->d.z;
         rB = sqrt(projB * projB + projC * projC);
       
         //rB = sqrt(ray->d.y * ray->d.y + ray->d.z * ray->d.z); // THIS LINE IS WRONG!!!!! 

        //need to recalculate thetaB after you modify z (this is a new addition, also questionable)          
        thetaB = acos(projB/rB);

        double newProjB = cos(thetaB) * rB;
        ray->d.z = sin(thetaB) * rB;
    
        //convert from ab space back to x,y space
        ray->d.x = direction.x * newProjA + orthoDirection.x * newProjB;
        ray->d.y = direction.y * newProjA + orthoDirection.y * newProjB;

        ray->d = Normalize(ray->d);       


        //------------end diffraction-----------------









    }

    // Compute offset rays for _PerspectiveCamera_ ray differentials
    ray->rxOrigin = ray->ryOrigin = ray->o;
    ray->rxDirection = Normalize(Vector(Pcamera) + dxCamera);
    ray->ryDirection = Normalize(Vector(Pcamera) + dyCamera);
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    ray->hasDifferentials = true;
    return 1.f;
}


PerspectiveDiffractionCamera *CreatePerspectiveDiffractionCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
    // Extract common camera parameters from _ParamSet_
    float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        swap(shutterclose, shutteropen);
    }
    float lensradius = params.FindOneFloat("lensradius", 0.f);
    float focaldistance = params.FindOneFloat("focaldistance", 1e30f);
    float frame = params.FindOneFloat("frameaspectratio",
        float(film->xResolution)/float(film->yResolution));
    float screen[4];
    if (frame > 1.f) {
        screen[0] = -frame;
        screen[1] =  frame;
        screen[2] = -1.f;
        screen[3] =  1.f;
    }
    else {
        screen[0] = -1.f;
        screen[1] =  1.f;
        screen[2] = -1.f / frame;
        screen[3] =  1.f / frame;
    }
    int swi;
    const float *sw = params.FindFloat("screenwindow", &swi);
    if (sw && swi == 4)
        memcpy(screen, sw, 4*sizeof(float));
    float fov = params.FindOneFloat("fov", 90.);
    float halffov = params.FindOneFloat("halffov", -1.f);
    if (halffov > 0.f)
        // hack for structure synth, which exports half of the full fov
        fov = 2.f * halffov;
    return new PerspectiveDiffractionCamera(cam2world, screen, shutteropen,
        shutterclose, lensradius, focaldistance, fov, film);
}
