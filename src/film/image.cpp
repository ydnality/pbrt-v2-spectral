
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


// film/image.cpp*
#include <sstream>
#include <iostream>
#include <fstream>
#include "stdafx.h"
#include "film/image.h"
#include "spectrum.h"
#include "parallel.h"
#include "imageio.h"
#include "floatfile.h"

// ImageFilm Method Definitions
ImageFilm::ImageFilm(int xres, int yres, Filter *filt, const float crop[4],
                     const string &fn, bool openWindow)
    : Film(xres, yres) {
    filter = filt;
    memcpy(cropWindow, crop, 4 * sizeof(float));
    filename = fn;
    // Compute film image extent
    xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
    xPixelCount = max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
    yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
    yPixelCount = max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);

    // Allocate film image storage
    pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);

    // Precompute filter weight table
#define FILTER_TABLE_SIZE 16
    filterTable = new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
    float *ftp = filterTable;
    for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
        float fy = ((float)y + .5f) *
                   filter->yWidth / FILTER_TABLE_SIZE;
        for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
            float fx = ((float)x + .5f) *
                       filter->xWidth / FILTER_TABLE_SIZE;
            *ftp++ = filter->Evaluate(fx, fy);
        }
    }

    // Possibly open window for image display
    if (openWindow || PbrtOptions.openWindow) {
        Warning("Support for opening image display window not available in this build.");
    }
    //Andy added
    nCMRows = nSpectralSamples;
    nCMCols = nSpectralSamples;
}

//Andy: this needs to be modified - modified
void ImageFilm::AddSample(const CameraSample &sample,
                          const Spectrum &L, const Ray &currentRay) {

    //std::cout << "got to AddSample()\n";
    // Compute sample's raster extent
    float dimageX = sample.imageX - 0.5f;
    float dimageY = sample.imageY - 0.5f;
    int x0 = Ceil2Int (dimageX - filter->xWidth);
    int x1 = Floor2Int(dimageX + filter->xWidth);
    int y0 = Ceil2Int (dimageY - filter->yWidth);
    int y1 = Floor2Int(dimageY + filter->yWidth);
    x0 = max(x0, xPixelStart);
    x1 = min(x1, xPixelStart + xPixelCount - 1);
    y0 = max(y0, yPixelStart);
    y1 = min(y1, yPixelStart + yPixelCount - 1);
    if ((x1-x0) < 0 || (y1-y0) < 0)
    {
        PBRT_SAMPLE_OUTSIDE_IMAGE_EXTENT(const_cast<CameraSample *>(&sample));
        return;
    }

    // Loop over filter support and add sample to pixel arrays
    //float xyz[3];
    //L.ToXYZ(xyz);
    
    //Andy: put in spectrum data here
    float origC[nSpectralSamples];
    L.GetOrigC(origC);


    // Precompute $x$ and $y$ filter table offsets
    int *ifx = ALLOCA(int, x1 - x0 + 1);
    for (int x = x0; x <= x1; ++x) {
        float fx = fabsf((x - dimageX) *
                         filter->invXWidth * FILTER_TABLE_SIZE);
        ifx[x-x0] = min(Floor2Int(fx), FILTER_TABLE_SIZE-1);
    }
    int *ify = ALLOCA(int, y1 - y0 + 1);
    for (int y = y0; y <= y1; ++y) {
        float fy = fabsf((y - dimageY) *
                         filter->invYWidth * FILTER_TABLE_SIZE);
        ify[y-y0] = min(Floor2Int(fy), FILTER_TABLE_SIZE-1);
    }
    bool syncNeeded = (filter->xWidth > 0.5f || filter->yWidth > 0.5f);
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            // Evaluate filter value at $(x,y)$ pixel
            int offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0];
            float filterWt = filterTable[offset];

            // Update pixel values with filtered sample contribution
            Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
            if (!syncNeeded) {
                /*pixel.Lxyz[0] += filterWt * xyz[0];
                pixel.Lxyz[1] += filterWt * xyz[1];
                pixel.Lxyz[2] += filterWt * xyz[2];*/
		        for (int i = 0; i < nSpectralSamples; i++)  //fixed some of the brackets... need to check
                {
			        pixel.c[i] += filterWt * origC[i];                    
                }
                pixel.weightSum += filterWt;
            }
            else {
                    // Safely update _Lxyz_ and _weightSum_ even with concurrency
		        for (int i = 0; i < nSpectralSamples; i++)
                {   
			        AtomicAdd(&pixel.c[i], filterWt * origC[i]);    
                }   
                AtomicAdd(&pixel.weightSum, filterWt);
            }
            pixel.Z += currentRay.maxt;    //add the depth map data member to a pixel - depth is inherently stored in Ray
        }
    }
}



//Andy: this needs to be modified to allow for more than 3 channel splatting.  Try to understand this more.  Do we need depth map processing here?
void ImageFilm::Splat(const CameraSample &sample, const Spectrum &L) {
    if (L.HasNaNs()) {
        Warning("ImageFilm ignoring splatted spectrum with NaN values");
        return;
    }
    
    //don't need this anymore
    //float xyz[3];
    //L.ToXYZ(xyz);

    int x = Floor2Int(sample.imageX), y = Floor2Int(sample.imageY);
    if (x < xPixelStart || x - xPixelStart >= xPixelCount ||
        y < yPixelStart || y - yPixelStart >= yPixelCount) return;
    Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);

	for (int i = 0; i < nSpectralSamples; i++)
	{
		AtomicAdd(&pixel.splatC[i], pixel.splatC[i]);
	}	
    
    //AtomicAdd(&pixel.splatXYZ[1], xyz[1]);
    //AtomicAdd(&pixel.splatXYZ[2], xyz[2]);
}


void ImageFilm::GetSampleExtent(int *xstart, int *xend,
                                int *ystart, int *yend) const {
    *xstart = Floor2Int(xPixelStart + 0.5f - filter->xWidth);
    *xend   = Floor2Int(xPixelStart + 0.5f + xPixelCount  +
                        filter->xWidth);

    *ystart = Floor2Int(yPixelStart + 0.5f - filter->yWidth);
    *yend   = Floor2Int(yPixelStart + 0.5f + yPixelCount +
                        filter->yWidth);
}


void ImageFilm::GetPixelExtent(int *xstart, int *xend,
                               int *ystart, int *yend) const {
    *xstart = xPixelStart;
    *xend   = xPixelStart + xPixelCount;
    *ystart = yPixelStart;
    *yend   = yPixelStart + yPixelCount;
}

//Andy added
void ImageFilm::ParseConversionMatrix(string filename){
    string fn = AbsolutePath(ResolveFilename(filename));
    nCMRows = nSpectralSamples;
    nCMCols = nSpectralSamples;
    
    //default wavespecify
    waveSpecify = new float[nCMRows];
    for (int i = 0; i < nCMCols; i++)
        waveSpecify[i] = sampledLambdaStart + (((sampledLambdaEnd - sampledLambdaStart)/nSpectralSamples) * (i - .5));  //TODO: double check int and float issues
    //populate the identity matrix
    float * identity = new float[nCMRows * nCMCols];
    for (int i =0; i < nCMRows * nCMCols; i++)
        identity[i] = 0.f;
    for (int i = 0; i < nCMRows; i++)
        identity[i * nCMCols + i] = 1.f;

    vector<float> vals;
    if (!ReadFloatFile(fn.c_str(), &vals)) {
        Warning("Unable to read conversion matrix file \"%s\".  Using identity matrix.", 
                fn.c_str());
        
        conversionMatrix = identity;
        return;
    }
    // TODO: consider list of wavelengths
    if (vals.size() < 2 || ((vals.size() - 2 - nCMRows) % nSpectralSamples != 0))  //changed to consider the CMRows amount of wave specifications
    {
        Warning("Incorrect conversion matrix file format (wrong number of matrix elements! Using identity matrix.");
        conversionMatrix = identity;
        return;
    }
    nCMRows = vals[0];
    nCMCols = vals[1];

    waveSpecify = new float[nCMRows];
    //read wavelength representation data
    for (int i = 2; i < 2  + nCMRows; i++)
    {
        waveSpecify[i-2] = vals[i];
    }    

    if (nCMRows * nCMCols != vals.size() -2 - nCMRows || nCMCols != nSpectralSamples)
    {
        Warning("Incorrect conversion matrix file format (wrong number of matrix elements! Using identity matrix.");
        conversionMatrix = identity;
        return;
    }    

    conversionMatrix = new float[nCMCols * nCMRows];

    for (int i = 2 + nSpectralSamples; i< vals.size(); i++)
    {
        conversionMatrix[i-2-nSpectralSamples] = vals[i];
        if (debugMode)
            std::cout << "conversionMatrix[" << i-2 << "]=" << conversionMatrix[i-2-nCMRows] << "\t";
    }    
    return;
}

//Andy changed
void ImageFilm::WriteImage(float splatScale) {

	//static const int sampledLambdaStart = 400;
	//static const int sampledLambdaEnd = 700;
	//static const int nSpectralSamples = 30;

    //Andy: will need to multiply by the conversion matrix for final output!

    // Convert image to RGB and compute final pixel values
    int nPix = xPixelCount * yPixelCount;
    //float *rgb = new float[3*nPix];  old code
    float *finalC = new float[nSpectralSamples * nPix];
    float *finalZ = new float[nPix * 3];

    int offset = 0;
    for (int x = 0; x < xPixelCount; ++x) {
        for (int y = 0; y < yPixelCount; ++y) {
        
            // Convert pixel XYZ color to RGB
            //XYZToRGB((*pixels)(x, y).Lxyz, &rgb[3*offset]);     //Andy: don't need this for now - want to keep it in spectrum
	        for (int ind = 0; ind < nSpectralSamples; ind++)
	        {
                //*pixels(x,y).c[ind] contains computed spectral intensity of the image
	            finalC[(y * xPixelCount + x)*nSpectralSamples + ind] = (*pixels)(x,y).c[ind];
	        }
            //depth map
            finalZ[(y*xPixelCount + x) * 3] = (*pixels)(x,y).Z;  //Andy: do we need to do the weighting business?
            finalZ[(y*xPixelCount + x) * 3 + 1] = (*pixels)(x,y).Z;  //do this 3 times so we have a grayscale .exr image
            finalZ[(y*xPixelCount + x) * 3 + 2] = (*pixels)(x,y).Z;  

            // Normalize pixel with weight sum
            // TODO: we may need to eliminate this "filtering" later
            float weightSum = (*pixels)(x, y).weightSum;  // Andy: will still need this, but need to make a for loop for all channels
            if (weightSum != 0.f) {
                float invWt = 1.f / weightSum;
		        //Andy: new for loop
		        for (int i = 0; i < nSpectralSamples; i++)
		        {
			        finalC[nSpectralSamples * offset + i] =  max(0.f, finalC[nSpectralSamples*offset + i  ]);  //TODO: investigate invWt

		        }
                finalZ[3 * offset ] = max(0.f, finalZ[3*offset] * invWt);
                finalZ[3 * offset + 1] = max(0.f, finalZ[3*offset + 1] * invWt);
                finalZ[3 * offset + 2] = max(0.f, finalZ[3*offset + 2] * invWt);
		        /*
                rgb[3*offset  ] = max(0.f, rgb[3*offset  ] * invWt);
                rgb[3*offset+1] = max(0.f, rgb[3*offset+1] * invWt);
                rgb[3*offset+2] = max(0.f, rgb[3*offset+2] * invWt);*/
            }

            //Add splat value at pixel
	        //Andy: will probably still need this as well, but once again need all channels
            //float splatRGB[3];
            //XYZToRGB((*pixels)(x, y).splatXYZ, splatRGB);
	        float * splatC = (*pixels)(x, y).splatC;
            /*rgb[3*offset  ] += splatScale * splatRGB[0];
            rgb[3*offset+1] += splatScale * splatRGB[1];
            rgb[3*offset+2] += splatScale * splatRGB[2];*/

	        //Andy: new
	        for (int i = 0; i < nSpectralSamples; i++)
	        {
	        	finalC[nSpectralSamples * offset + i] += splatScale * splatC[nSpectralSamples]; 
		        //std::cout << finalC[nSpectralSamples * offset + i] << "\t";
	        }
            ++offset;
        }
    }

    
    float *finalCMultiplied = new float[nCMRows * nPix];
    //Andy added: multiply by conversion matrix to produce proper output
    for (int y = 0; y < yPixelCount; ++y) {
        for (int x = 0; x < xPixelCount; ++x) {
            //for each pixel, multiply conversion matrix by existing data
            //these nested for loops are for matrix multiplication
            for (int row = 0; row < nCMRows; row++)
            {
                float tempSum = 0;
                for (int iter = 0; iter < nCMCols; iter++)
                {
                    tempSum += conversionMatrix[row * nCMCols + iter] * finalC[nSpectralSamples * (y * xPixelCount + x) + iter];    //matrix multiplication by column vector
                }
                finalCMultiplied[nCMRows * (x *  yPixelCount + y) + row] = tempSum;  //check this later
            }

        }
    }    

    //declare text file stream
    std::ofstream myfile;
    int lastPos = filename.find_last_of(".");
    string newFileName = filename.substr(0, lastPos) + ".dat";

    
    myfile.open (newFileName.c_str());


    //print out the dimensions of the image    
    myfile << xPixelCount << " " << yPixelCount << " " << nCMRows << "\n";
    //for (int i = 0; i < nCMRows; i++)
    //    myfile << waveSpecify[i] << " ";
    //myfile << "\n";

    myfile.close();

    //open file for binary writing purposes
    FILE * spectralImageBin;
	spectralImageBin = fopen(newFileName.c_str(), "a");


    //TODO: perhaps dump the conversion matrix here

    //Write Binary image
    for (int i = 0; i < nCMRows; i++)
    {
		for (int j = 0; j < nPix; j++)
		{ 
            double r = (double)finalCMultiplied[nCMRows * j + i];
            fwrite((void*)(&r), sizeof(r), 1, spectralImageBin);
            //myfile << finalCMultiplied[nCMRows * j + i];
            //if (j < nPix - 1)
            //    myfile << " ";
		}
        //myfile << "\n";
    }

    fclose(spectralImageBin);

    //::WriteImage(filename, rgb, NULL, xPixelCount, yPixelCount,
    //             xResolution, yResolution, xPixelStart, yPixelStart);

    //write .exr depth map!
    string newFileNameDepth = filename.substr(0, lastPos) + "_depth.exr";
    ::WriteImage(newFileNameDepth, finalZ, NULL, xPixelCount, yPixelCount,
                xResolution, yResolution, xPixelStart, yPixelStart);

    // Release temporary image memory
    delete[] finalC;
    delete[] finalCMultiplied;
    delete[] finalZ;
    //delete[] rgb;
}


void ImageFilm::UpdateDisplay(int x0, int y0, int x1, int y1,
    float splatScale) {
}

//Andy: added conversionmatrix file here
ImageFilm *CreateImageFilm(const ParamSet &params, Filter *filter) {
    string filename = params.FindOneString("filename", PbrtOptions.imageFile);
    if (filename == "")
#ifdef PBRT_HAS_OPENEXR
        filename = "pbrt.exr";
#else
        filename = "pbrt.tga";
#endif

    int xres = params.FindOneInt("xresolution", 640);
    int yres = params.FindOneInt("yresolution", 480);     

    if (PbrtOptions.quickRender) xres = max(1, xres / 4);    
    if (PbrtOptions.quickRender) yres = max(1, yres / 4);
    bool openwin = params.FindOneBool("display", false);
    float crop[4] = { 0, 1, 0, 1 };
    int cwi;
    const float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
        crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
        crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
        crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
    }

    
    //Andy changed to allow for parsing of conversionmatrixfile name
    ImageFilm * newFilm = new ImageFilm(xres, yres, filter, crop, filename, openwin);
    string conversionMatrixFilename = params.FindOneString("conversionmatrixfile", "");   //params contains all the parsed attributes of film... we now look for the "conversionmatrixfile" entry, which is a string
    std::cout<< "\nconversionMatrixFilename: " << conversionMatrixFilename << "\n\n";
    newFilm->ParseConversionMatrix(conversionMatrixFilename);
    
    return newFilm;
}


