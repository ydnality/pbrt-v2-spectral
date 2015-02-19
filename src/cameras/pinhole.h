#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_PINHOLE_H
#define PBRT_CAMERAS_PINHOLE_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include <gsl/gsl_randist.h>

/*struct LensElement{
    float radius;
    float separation;
    float n;
    float aperture;
};*/

class pinholeCamera : public Camera {
public:
   pinholeCamera(const AnimatedTransform &cam2world,
                                 float filmdistance, 
                                 float filmdiag,
								 Film *f);
   ~pinholeCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
   float GenerateCameraRay(const CameraSample &sample, Ray *) const;
   float getSensorWidth();
private:
   float ShutterOpen;
   float ShutterClose;
   float filmDiag;
   float filmDistance;

   gsl_rng * r;
   Film * film;
};

pinholeCamera *CreatePinholeCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
