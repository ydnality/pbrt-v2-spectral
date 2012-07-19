// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
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
#include "../vdb.h"

using namespace std;

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);
	   // Realistic camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }

	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter_,
                                 const string &specfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   film(f)
{

	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.
    
    filmDiag = filmdiag;
    filmDistance = filmdistance;
    //parse the specfile 
    string fn = AbsolutePath(ResolveFilename(specfile));


    vector<float> vals;

    //check to see if there is valid input in the lens file.
    if (!ReadFloatFile(fn.c_str(), &vals)) {
        Warning("Unable to read lens file!");
        return;
    }

    //check to see if the number of floats is a multiple of 4
    if (vals.size() % 4 != 0)
    {
        Warning("Wrong number of float values in lens file!");
        return;
    }

    for (int i = 0; i < vals.size(); i+=4)
    {
        
        std::cout << vals[i] << "  " << vals[i+1] << "  " << vals[i+2] << "  " << vals[i+3] << "\n";
        LensElement currentLensEl;
        currentLensEl.radius = vals[i];
        currentLensEl.separation = vals[i+1];
        currentLensEl.n = vals[i+2];
        currentLensEl.aperture = vals[i+3];
        lensEls.push_back(currentLensEl);
    }   

    for (int i = lensEls.size()-1; i >=0 ; i--)
    {   
        cout << "radius: " << lensEls[i].radius << " ";
        cout << "separation: " << lensEls[i].separation << " ";
        cout << "n: " << lensEls[i].n << " ";
        cout << "aperture: " << lensEls[i].aperture << "\n";
    }
}


RealisticCamera::~RealisticCamera()
{

}

void RealisticCamera::runLensFlare(const Scene * scene, const Renderer * renderer) const
{
/*    std::cout << "in runLensFlare()!\n";
    
    for (int i = 0; i < scene->lights.size(); i++)
    {
        std::cout << "isPointLight?: " << scene->lights[i]->isPointLight() << "\n";
        if (!scene->lights[i]->isPointLight())
            continue;
        
        //if it is a point light, make lens flare!  shoot rays into the lens! 
        Point pointLight = ((const PointLight*)(scene->lights[i]))->getLightPos();
        std::cout << "pointLight position: " << pointLight.x << ", " << pointLight.y << ", " << pointLight.z << "\n";

        //now we want to shoot rays uniformly from the light source to the front of the lens
        
    }*/
    
}

bool RealisticCamera::IntersectLensEl(const Ray &r, float *tHit, float radius, float dist, Vector & normalVec) const{
    float phi;
    Point phit;
    // Transform _Ray_ to object space

    Transform shiftZ =  Translate(Vector(0,0,radius - dist));
    Ray ray = r;
    (shiftZ)(r, &ray);
    //(*WorldToObject)(r, &ray);


    if (radius < 0)
        radius = -radius;
    // Compute quadratic sphere coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
              ray.o.z*ray.o.z - radius*radius;
//std::cout << "ray.d: " << ray.d.x << " " << ray.d.y << " " << ray.d.z <<"\n";
//std:: cout << "A: " << A << " B: " << B << " C: " << C << "\n";
    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;
//std::cout <<"got here";
    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
//std::cout <<"got here1";
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    // Compute sphere hit position and $\phi$
    phit = ray(thit);
    if (phit.x == 0.f && phit.y == 0.f) phit.x = 1e-5f * radius;
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;
//std::cout << "got here2!";
    // Test sphere intersection against clipping parameters
    //we don't need this test for now -> assume lens is within clipping parameters
//    if ((zmin > -radius && phit.z < zmin) ||
//        (zmax <  radius && phit.z > zmax) || phi > phiMax) {
//        if (thit == t1) return false;
//        if (t1 > ray.maxt) return false;
//        thit = t1;
        // Compute sphere hit position and $\phi$
//        phit = ray(thit);
//        if (phit.x == 0.f && phit.y == 0.f) phit.x = 1e-5f * radius;
//        phi = atan2f(phit.y, phit.x);
//        if (phi < 0.) phi += 2.f*M_PI;
 //       if ((zmin > -radius && phit.z < zmin) ||
 //           (zmax <  radius && phit.z > zmax) || phi > phiMax)
 //           return false;
 //   }

/*
    // Find parametric representation of sphere hit
    float u = phi / phiMax;
    float theta = acosf(Clamp(phit.z / radius, -1.f, 1.f));
    float v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Compute sphere $\dpdu$ and $\dpdv$
    float zradius = sqrtf(phit.x*phit.x + phit.y*phit.y);
    float invzradius = 1.f / zradius;
    float cosphi = phit.x * invzradius;
    float sinphi = phit.y * invzradius;
    Vector dpdu(-phiMax * phit.y, phiMax * phit.x, 0);
    Vector dpdv = (thetaMax-thetaMin) *
        Vector(phit.z * cosphi, phit.z * sinphi,
               -radius * sinf(theta));

    // Compute sphere $\dndu$ and $\dndv$
    Vector d2Pduu = -phiMax * phiMax * Vector(phit.x, phit.y, 0);
    Vector d2Pduv = (thetaMax - thetaMin) * phit.z * phiMax *
                    Vector(-sinphi, cosphi, 0.);
    Vector d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
                    Vector(phit.x, phit.y, phit.z);

    // Compute coefficients for fundamental forms
    float E = Dot(dpdu, dpdu);
    float F = Dot(dpdu, dpdv);
    float G = Dot(dpdv, dpdv);
    Vector N = Normalize(Cross(dpdu, dpdv));
    float e = Dot(N, d2Pduu);
    float f = Dot(N, d2Pduv);
    float g = Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    float invEGF2 = 1.f / (E*G - F*F);
    Normal dndu = Normal((f*F - e*G) * invEGF2 * dpdu +
                         (e*F - f*E) * invEGF2 * dpdv);
    Normal dndv = Normal((g*F - f*G) * invEGF2 * dpdu +
                         (f*F - g*E) * invEGF2 * dpdv);*/

    // Initialize _DifferentialGeometry_ from parametric information
    //const Transform &o2w = *ObjectToWorld;
    //*dg = DifferentialGeometry(o2w(phit), o2w(dpdu), o2w(dpdv),
      //                         o2w(dndu), o2w(dndv), u, v, this);

    // Update _tHit_ for quadric intersection
    *tHit = thit;

    // Compute _rayEpsilon_ for quadric intersection
    //*rayEpsilon = 5e-4f * *tHit;
    normalVec = Normalize(Vector(ray.d.x * thit + ray.o.x, ray.d.y*thit + ray.o.y, ray.d.z*thit + ray.o.z));
    return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
  // YOUR CODE HERE -- make that ray!

  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens

  // GenerateRay() should return the weight of the generated ray

    Point startingPoint;
//std::cout << "imageX: " << sample.imageX << "imageY: " << sample.imageY << "\n";
    startingPoint.x = -(sample.imageX - film->xResolution/2.f)/(film->xResolution/2.f);
    startingPoint.y = (sample.imageY - film->yResolution/2.f)/(film->yResolution/2.f);
    startingPoint.z = -filmDistance;  
    
    
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    float width = filmDiag /sqrt((1.f + aspectRatio * aspectRatio));
    float height = width/aspectRatio;
    
    startingPoint.x *= width/2.f;
    startingPoint.y *= height/2.f;

    float lensU, lensV;
    float prevN = 1;

    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    
    //use geometry to solve for z intercept of disc on sphere
    float firstAperture = lensEls[lensEls.size()-1].aperture/2;
    float firstRadius = lensEls[lensEls.size()-1].radius;
    float zIntercept = (-firstRadius - sqrt(firstRadius * firstRadius - firstAperture * firstAperture));
//std::cout << " firstRadius: " << firstRadius;
//std::cout << " firstAperture: " << firstAperture;
//std::cout << " zIntercept: " << zIntercept;
    lensU *= firstAperture;
    lensV *= firstAperture;
    //vdb visualize sampling surface
    Point pointOnLens = Point(lensU, lensV, zIntercept);   //can we even assume that lens is a flat disk, maybe this has problems later?
   
    *ray = Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
    ray->o = startingPoint;    //initialize ray origin
    ray->d = Normalize(pointOnLens - ray->o);

    // vdb_color(.7, 0, 0);
    // vdb_point(pointOnLens.x, pointOnLens.y, pointOnLens.z); 
    

    float lensDistance = 0;
    for (int i = lensEls.size()-1; i>=0 ; i--)
    {

//if (i == 5)  //this doesn't run!!
//    std::cout << "i =0 ";

//std::cout << "width: " << width << " height: " << height << "\n";
//std::cout << "startX: " << startingPoint.x << " " << "startY: " << startingPoint.y << " ";

        //generate a uniform sample
        //float lensRadius = -40;   //dummy value for now!! (mm)
        float lensRadius = lensEls[i].radius;   //(mm)
        lensDistance += lensEls[i].separation;
        float currentAperture = lensEls[i].aperture;
        //float currentAperture = 10;

// std::cout << "lensDistance: " << lensDistance << "\n";        
        //lensDistance = 0;
        
        
        // Compute point on plane of focus
        //float focalDistance = 36.5;  //dummy value for now!
        //float ft = focalDistance / ray->d.z;
        //Point Pfocus = (*ray)(ft);

        ray->o = startingPoint;
//std::cout << "ray->o: " << ray->o.x << " " << ray->o.y << " " << ray->o.z 
//          << "ray->d: " << ray->d.x << " " << ray->d.y << " " << ray->d.z << "   ";
        float dummy = 0;
        float *tHit;
        tHit = &dummy;

        
        bool intersected = false;
        Vector normalVec(0,0,1);

        if (lensRadius ==0)
        {
            
            float tAperture = (lensDistance - ray->o.z)/(ray->d.z);
            float xAperture = ray->o.x + ray->d.x * tAperture;
            float yAperture = ray->o.y + ray->d.y * tAperture;
            if (xAperture * xAperture + yAperture + yAperture > currentAperture * .5)
                return 0.f;            
            
            //intersected = true;
            //lensRadius = 1000000;
            //normalVec = Vector(0,0,1);
        }
        else
        {
            intersected = IntersectLensEl(*ray, tHit, lensRadius, lensDistance, normalVec);   //these are just test values for now!
        
            if (intersected)
            {
                //paint the intersection point green
                Point intersectPoint(0,0,0);
                intersectPoint.x = *tHit * ray->d.x + ray->o.x;
                intersectPoint.y = *tHit * ray->d.y + ray->o.y;
                intersectPoint.z = *tHit * ray->d.z + ray->o.z;

                //aperture test
                if (intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y >= currentAperture * currentAperture / 4)
                {
                    return 0.f;
                }

                //vdb visualization
    //std::cout << "*tHit: " << *tHit << "\n";
                float currentT = 0;
                /*
                while (currentT < *tHit)
                {
                    Point currentPoint(0,0,0);
                    currentPoint.x = currentT * ray->d.x + ray->o.x;
                    currentPoint.y = currentT * ray->d.y  + ray->o.y;
                    currentPoint.z = currentT * ray->d.z  + ray->o.z;
                    
                   // if (currentT == 0)
                   //     vdb_color(.7, 0, 0);
                   // else
                   //     vdb_color(0, 0, .7);
                   // vdb_point(currentPoint.x, currentPoint.y, currentPoint.z);   
                    currentT += .5;     
                }*/


                //vdb_color(0, .7, 0);
               // vdb_point(intersectPoint.x, intersectPoint.y, intersectPoint.z); 

                //-------------- snell's law----------------
                float n1 = lensEls[i].n;  
                       
                float n2 = 1;
                if (i-1 >= 0)
                    n2 = lensEls[i-1].n;
                
                //lens overall aperture case
                if (n1 == 0)
                    n1 = 1;   
                if (n2 ==0)
                    n2 = 1;

                Vector s1 = ray->d;
                if (lensRadius >0)
                    normalVec = -normalVec;
              

                float radicand = 1 - (n1/n2) * (n1/n2) * Dot(Cross(normalVec, s1), Cross(normalVec, s1));
                if (radicand < 0)
                    return 0;   //reflection, no refraction - might want to change for lens flare!

                Vector s2 = n1/n2 * (Cross(normalVec, Cross(-1 * normalVec, s1))) - normalVec * sqrt(radicand);
                
                //ray->d = normalVec;
                //reinitialize ray->o
                
                ray->d = Normalize(s2);  //reassign the direction to the new direction

                //---------------end snell's law------------

                //update starting point to current point
                startingPoint.x = intersectPoint.x;
                startingPoint.y = intersectPoint.y;
                startingPoint.z = intersectPoint.z;            
            }
            else
            {     
                //did not intersect -> draw red rays
                float currentT = 0;
                /*while (currentT < 10)
                {
                    Point currentPoint(0,0,0);
                    currentPoint.x = currentT * ray->d.x + ray->o.x;
                    currentPoint.y = currentT * ray->d.y  + ray->o.y;
                    currentPoint.z = currentT * ray->d.z  + ray->o.z;
                    
                    //if (currentT == 0)
                    //    vdb_color(.7, 0, 0);
                    //else
                    //    vdb_color(1, 0, 0);
                    //vdb_point(currentPoint.x, currentPoint.y, currentPoint.z);   
                    currentT += .5;     
                }*/
                return 0.f;
            }
        }
    }

    ray->time = Lerp(sample.time, ShutterOpen, ShutterClose);
    CameraToWorld(*ray, ray);
    ray->d = Normalize(ray->d); 
    return 1.f;
    

//this is the code from perspective camera to start out with...

/*
// Generate raster and camera samples

    float lensRadius = .02;   //dummy values for now
    float focalDistance = .75;

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
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    return 1.f;
*/


}
