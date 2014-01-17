#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_IDEAL_H
#define PBRT_CAMERAS_IDEAL_H

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

class idealDiffractionCamera : public Camera {
public:
   idealDiffractionCamera(const AnimatedTransform &cam2world,
                                 float filmdistance, float aperture_diameter_,
								 float focallength,
                                 float filmdiag,
								 Film *f, 
                                 bool diffractFlag,
                                 bool chromaticFlag,
                                 float xOffset,
                                 float yOffset);
   ~idealDiffractionCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
   float GenerateCameraRay(const CameraSample &sample, Ray *) const;
    void runLensFlare(const Scene * scene, const Renderer * renderer) const;
   float getFStop();
   float getFocalLength();
   float getSensorWidth();
private:
   float ShutterOpen;
   float ShutterClose;
   float apertureDiameter;
   float filmDiag;
   float filmDistance;
   float xApertureOffset;
   float yApertureOffset;
   gsl_rng * r;
   Film * film;
   //vector<LensElement> lensEls;
   bool IntersectLensEl(const Ray &r, float *tHit, float radius, float dist,  Vector & normalVec) const;
   bool diffractionEnabled;
   bool chromaticAberrationEnabled;
   float fstop;
   float focalLength;

;
};

idealDiffractionCamera *CreateIdealDiffractionCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
