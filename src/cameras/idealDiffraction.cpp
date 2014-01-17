// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/idealDiffraction.h"
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

idealDiffractionCamera *CreateIdealDiffractionCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}

	   float focallength = params.FindOneFloat("focal_length", 50.0); 
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
       		cout << "filmdistance: " << filmdistance << "\n";
	   float apdiameter = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
       
       float xOffset = params.FindOneFloat("x_aperture_offset", 0);
       float yOffset = params.FindOneFloat("y_aperture_offset", 0);

       //flags for enabling diffraction - default yes
       float diffractionEnabled = params.FindOneFloat("diffractionEnabled", 1.0);
       bool diffractFlag = diffractionEnabled == 1.f;
       //flags for enabling chromatic aberration - default no
       float chromaticAberrationEnabled = params.FindOneFloat("chromaticAberrationEnabled", 0.0);
       bool chromaticFlag =  chromaticAberrationEnabled ==1.f;


	   return new idealDiffractionCamera(cam2world, filmdistance, apdiameter,focallength, 
	       filmdiag, film, diffractFlag, chromaticFlag, xOffset, yOffset);
}




idealDiffractionCamera::idealDiffractionCamera(const AnimatedTransform &cam2world,
                                 float filmdistance, float aperture_diameter_,
								 float focallength,
                                 float filmdiag,
								 Film *f, 
                                 bool diffractFlag,
                                 bool chromaticFlag,
                                 float xOffset,
                                 float yOffset)
                                 : Camera(cam2world, -1, -1, f),
								   ShutterOpen(-1),
								   ShutterClose(-1),
								   film(f)
{

	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.

    focalLength = focallength;
    filmDiag = filmdiag;
    filmDistance = filmdistance;

    //parse the specfile 
    //string fn = AbsolutePath(ResolveFilename(specfile));

    diffractionEnabled = diffractFlag;
    chromaticAberrationEnabled = chromaticFlag;
    xApertureOffset = xOffset;
    yApertureOffset = yOffset;
	apertureDiameter = aperture_diameter_;


	std::cout <<"focalLength: " << focalLength << "\n";
	std::cout <<"filmDistance: " << filmDistance << "\n";
	std::cout <<"filmDiag: " << filmDiag << "\n";
    std::cout <<"xApertureOffset: " << xApertureOffset << "\n";
    std::cout <<"yApertureOffset: " << yApertureOffset << "\n";
    std::cout <<"DiffractionEnabled: " << diffractionEnabled << "\n";
    std::cout <<"chromaticAberrationEnabled: " << chromaticAberrationEnabled << "\n";
/*
    vector<float> vals;


    //check to see if there is valid input in the lens file.
    if (!ReadFloatFile(fn.c_str(), &vals)) {
        Warning("Unable to read lens file!");
        return;
    }

    //check to see if the number of floats is a multiple of 4
    if ((vals.size()-1) % 4 != 0)
    {
        Warning("0Wrong number of float values in lens file!  Check file format!  Did you forget to specify the focal length?");
        return;
    }

    float focallength = vals[0];   //read the focal length - this is new
    focalLength = focallength;
    std::cout << "focalLength :" << focalLength << "\n";
    fstop = focalLength/aperture_diameter_;

    std::cout << "apertureDiameter: " << aperture_diameter_ << "\n";
    std::cout << "fstop: f/" << fstop << "\n";
    
    for (int i = 1; i < vals.size(); i+=4)
    {
        std::cout << vals[i] << "  " << vals[i+1] << "  " << vals[i+2] << "  " << vals[i+3] << "\n";
        LensElement currentLensEl;
        currentLensEl.radius = vals[i];
        currentLensEl.separation = vals[i+1];
        currentLensEl.n = vals[i+2];
        currentLensEl.aperture = vals[i+3];
        if (currentLensEl.radius ==0)
            currentLensEl.aperture = aperture_diameter_;
        lensEls.push_back(currentLensEl);
    }   

    for (int i = lensEls.size()-1; i >=0 ; i--)
    {   
        cout << "radius: " << lensEls[i].radius << " ";
        cout << "separation: " << lensEls[i].separation << " ";
        cout << "n: " << lensEls[i].n << " ";
        cout << "aperture: " << lensEls[i].aperture << "\n";
    }
*/

//code borrowed from http://www.gnu.org/software/gsl/manual/gsl-ref.html#Random-Number-Generation  18.13
       const gsl_rng_type * T;
     
       int i, n = 10;
     
       gsl_rng_env_setup();
     
       T = gsl_rng_default;
       r = gsl_rng_alloc (T);


    //TODO: TELL FILM THAT THIS IS DEPTH MAP PROCESSING


}


idealDiffractionCamera::~idealDiffractionCamera()
{

}


float idealDiffractionCamera::getFStop()
{
	return fstop;
}

float idealDiffractionCamera::getFocalLength()
{
	return focalLength;
}

void idealDiffractionCamera::runLensFlare(const Scene * scene, const Renderer * renderer) const
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

bool idealDiffractionCamera::IntersectLensEl(const Ray &r, float *tHit, float radius, float dist, Vector & normalVec) const{
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

float idealDiffractionCamera::getSensorWidth()
{
    
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    //float width = filmDiag /sqrt((1.f + aspectRatio * aspectRatio));
    float width = filmDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));
    return width;
}

float idealDiffractionCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
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

	//generate a uniform set of samples on a disk of size unit circle
    float lensU, lensV;
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    lensU *= apertureDiameter/2.f;  //scale by aperture radius
    lensV *= apertureDiameter/2.f;

	//assign initial ray origin and directionS from sensor to random point on aperture
    Point pointOnLens = Point(lensU, lensV, 0);  //lens is at z = 0
    float tempWavelength = ray->wavelength; 
    *ray = Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
    ray->o = startingPoint;    //initialize ray origin
    ray->d = Normalize(pointOnLens - ray->o);
    ray->wavelength = tempWavelength;  //so that wavelength information is retained.

    // vdb_color(.7, 0, 0);
    // vdb_point(pointOnLens.x, pointOnLens.y, pointOnLens.z); 

	//-------new ideal lens stuff
	//we are calculating the refraction angle of the rays according to an ideal thin lens using the thin lens equation
	//we assume that rays through the middle of the lens go straight forever.
	//the other rays 

	//compute ray going from sensor to the center point of lens
	Ray * centerRay = new Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
	Point centerOfLens = Point(0,0,0);
	centerRay->o = startingPoint;
	centerRay->d = Normalize(centerOfLens - startingPoint);
	centerRay->o = centerOfLens;  //put origin at origin (center of lens)

	//use thin lens equation to compute new point of focus
	//this is the distance to a point that is in focus.  all rays from a single point must converge at this distance.  
	float focalDistance = 0;	
	float denominator = 1.f/focalLength - 1.f/filmDistance; 
	if (abs(focalLength - filmDistance) < .0000000001f)   //protect against division by 0
		focalDistance = 1000000 * focalLength;
	else
		focalDistance = 1.f/denominator;
	
//cout << "focalLength: " << focalLength << "\n";
//cout << "filmDistance: " << filmDistance << "\n";
//cout << "focalDistance: " << focalDistance << "\n";

	//compute point of focus in 3 space (world coordinates).
	//first, compute t intersection with in focus plane - then position
	float ft = focalDistance / centerRay->d.z;   
	//Pfocus is the object point in world coordinates
    Point Pfocus = (*centerRay)(ft);  //PBRT notation assigns point on a ray, if given the desired time

	//intersectPoint is a point in world coordinates corresponding to where the ray hits the pupil plane
	Point intersectPoint;
	intersectPoint.x = lensU;
	intersectPoint.y = lensV;
	intersectPoint.z = 0;

	//set the ray origin and directionS from the lens into the scene
	//set the lens intersection point as the origin
	//then, set the new directionS of ray to converge to the focus point (Pfocus)
	ray->o = intersectPoint;
	ray->d = Normalize(Pfocus - intersectPoint);

    // --------------------add effect of diffraction----------------
    if (diffractionEnabled)
    {
		//variables that are returned by random gaussian function
        double inits = 0;
        double initl = 0;       
        double * noiseSPointer = &inits;
        double * noiseLPointer = &initl;   

        //calculate min distance directionS
       
        //calculate the distance of the intersect point to the center of the lens
        //Point intersectPoint(lensU, lensV, 0.f);
        double ipLength = sqrt(intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y );

//cout << "radius: " << radius << "\n";
//std::cout << "rayWavelength: " << ray->wavelength << "\n";

        //calculate directionS and orthogonal directionL
        Vector directionS(intersectPoint.x, intersectPoint.y, 0);
        Vector directionL(-intersectPoint.y, intersectPoint.x, 0);  
        directionS = Normalize(directionS);
        directionL = Normalize(directionL);
        
		double apertureRadius = apertureDiameter/2;
        double pointToEdgeS = apertureRadius- ipLength;   //this is 'a' from paper  //pointToEdgeS stands for point to edge short
        double pointToEdgeL = sqrt((apertureRadius* apertureRadius) - ipLength * ipLength);  //pointToEdgeS stands for point to edge long

        double pi = 3.14159265359;
        double lambda = ray->wavelength * 1e-9;  //this converts lambda to meters
        double sigmaS = atan(1/(2 * pointToEdgeS *.001 * 2 * pi/lambda));  //the .001 converts mm to m
        double sigmaL = atan(1/(2 * pointToEdgeL * .001 * 2 * pi/lambda)); 

		//this function regenerates a 2D gaussian sample and returns it in x and y 
        //gsl_ran_bivariate_gaussian (const gsl_rng * r, double sigma_x, double sigma_y, double rho, double * x, double * y)             
        gsl_ran_bivariate_gaussian (r, sigmaS, sigmaL, 0, noiseSPointer, noiseLPointer);    //experiment for now

        //calculate component of these vectors based on 2 random degrees

		//assign noise in the s and l directions according to data at these pointers
        float noiseS = (float)(*noiseSPointer);
        float noiseL = (float)(*noiseLPointer);

         //project the original ray (in world coordinates) onto a new set of basis vectors in the s and l directions         
         double projS = (ray->d.x * directionS.x + ray->d.y * directionS.y)/sqrt(directionS.x * directionS.x + directionS.y * directionS.y);
         double projL = (ray->d.x * directionL.x + ray->d.y * directionL.y)/sqrt(directionL.x * directionL.x + directionL.y * directionL.y);
         
		 //double projC = ray->d.z;

         //double rA = sqrt(projS * projS + projC * projC);
         //double rB = sqrt(projL * projL + projC * projC);
         //double thetaS = acos(projS/rA);          
         //double thetaL = acos(projL/rB);

		 //calculate angles (also in s and l coordinate system)  l
         
		 //double thetaS = atan(projS/ray->d.z);          
         //double thetaL = atan(projL/ray->d.z);

		

		double thetaA = atan(projS/ray->d.z);   //azimuth - this corresponds to sigmaS
		double thetaE = atan(projL/sqrt(projS*projS + ray->d.z * ray->d.z));   //elevation - this corresponds to sigmaL

         //add uncertainty
         thetaA = thetaA + noiseS;  
         thetaE = thetaE + noiseL;
         
         //convert angles back into cartesian coordinates, but in s,l space
         //double newprojS = tan(thetaS);

			 //old stuff... 
		     //ray->d.z = cos(thetaS) * rA;
		     //projC = ray->d.z;
		     //rB = sqrt(projL * projL + projC * projC);
		     //thetaL = asin(projL/rB);
		    //need to recalculate thetaL after you modify z (this is a new addition, also questionable)   
		    //double newprojL = sin(thetaL) * rB;
		    //ray->d.z = cos(thetaL) * rB;
		       
		//double newprojL = tan(thetaL);

		double newprojL = sin(thetaE);
		double smallH = cos(thetaE);   //smallH corresponds to the projection of the ray onto the s-z plane
		double newprojS = smallH * sin(thetaA);
		ray->d.z = smallH * cos(thetaA); 

        //convert from s-l space back to x-y space
        ray->d.x = (directionS.x * newprojS + directionL.x * newprojL)/sqrt(directionS.x * directionS.x + directionL.x * directionL.x);
        ray->d.y = (directionS.y * newprojS + directionL.y * newprojL)/sqrt(directionS.y * directionS.y + directionL.y * directionL.y);

		//old code /// ray->d.z = 1;   // why is this assigned as 1?


        ray->d = Normalize(ray->d);       
    //------------end diffraction-----------------        
    }

    ray->time = Lerp(sample.time, ShutterOpen, ShutterClose);
    CameraToWorld(*ray, ray);
    ray->d = Normalize(ray->d); 
	delete centerRay;
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





