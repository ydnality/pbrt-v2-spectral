// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/pinhole.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/spectralImage.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"
#include "floatfile.h"
#include "camera.h"
#include "scene.h"
#include "lights/point.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <time.h>
#include "../vdb.h"
//#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

pinholeCamera *CreatePinholeCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {

	   float filmdistance = params.FindOneFloat("filmdistance", 70.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	   return new pinholeCamera(cam2world, filmdistance, filmdiag, film);
}




pinholeCamera::pinholeCamera(const AnimatedTransform &cam2world,
                                 float filmdistance,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, -1, -1, f),
								   ShutterOpen(-1),
								   ShutterClose(-1),
								   film(f)
{

    filmDiag = filmdiag;
    filmDistance = filmdistance;
	std::cout <<"filmDistance: " << filmDistance << "\n";
	std::cout <<"filmDiag: " << filmDiag << "\n";

}


pinholeCamera::~pinholeCamera()
{

}



float pinholeCamera::getSensorWidth()
{
    
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    //float width = filmDiag /sqrt((1.f + aspectRatio * aspectRatio));
    float width = filmDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));
    return width;
}

float pinholeCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens

  // GenerateRay() should return the weight of the generated ray

	//compute starting Point using sensor size and aspect ratio and film distance
    Point startingPoint;
	//std::cout << "imageX: " << sample.imageX << "imageY: " << sample.imageY << "\n";
    startingPoint.x = -(sample.imageX - film->xResolution/2.f)/(film->xResolution/2.f);  //film is the sensor.  range [-1, 1]
    startingPoint.y = (sample.imageY - film->yResolution/2.f)/(film->yResolution/2.f);
    startingPoint.z = -filmDistance;  

	//cout << "startingPoint.z: " << startingPoint.z;
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    //float width = filmDiag /sqrt((1.f + aspectRatio * aspectRatio));
    float width = filmDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));  
    float height = width/aspectRatio;
    startingPoint.x *= width/2.f;  //scales starting point to correspond to physical width
    startingPoint.y *= height/2.f;

	//assign initial ray origin and directionS from sensor to random point on aperture
    float tempWavelength = ray->wavelength; 
    *ray = Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
    ray->wavelength = tempWavelength;  //so that wavelength information is retained.

	//compute ray going from sensor to the center point of lens
	Ray * centerRay = new Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
	Point centerOfLens = Point(0,0,0);
	centerRay->o = startingPoint;
	centerRay->d = Normalize(centerOfLens - startingPoint);
	centerRay->o = centerOfLens;  //put origin at origin (center of lens)
	ray->o = centerRay->o;
	ray->d = centerRay->d;

    ray->time = Lerp(sample.time, ShutterOpen, ShutterClose);
    CameraToWorld(*ray, ray);
    ray->d = Normalize(ray->d); 
	delete centerRay;
    return 1.f;
}





