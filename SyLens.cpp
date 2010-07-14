/*
*	This plugin uses the lens distortion model provided by Russ Anderson of the SynthEyes camera tracker
*	fame. It is made so that it's output is _identical_ to the Image Preparation tool of the SY camera tracker.
*	so you can use it INSTEAD of the footage Syntheyes prerenders.
*	It is largely based on tx_nukeLensDistortion by Matti Gruener of Trixter Film and Rising Sun fame. However,
*	For filtering we also use the standard filtering methods provided by Nuke.
*	
*	Written by Julik Tarkhanov in Amsterdam in 2010. me(at)julik.nl
*	with kind support by HecticElectric.
*	The beautiful Crimean landscape shot provided by Tim Parshikov, 2010.
*	
*	The code has some more comments than it's 3DE counterpart since we have to do some things that the other plugin
*	did not (like expanding the image output)
*/

extern "C" {

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <string.h>
}


#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <cstdlib> 

// Nuke-Data
#include "DDImage/Iop.h"
#include "DDImage/Row.h"
#include "DDImage/Tile.h"
#include "DDImage/Pixel.h"
#include "DDImage/Filter.h"
#include "DDImage/Knobs.h"
#include "DDImage/Knob.h"
#include "DDImage/Vector2.h"
#include "DDImage/DDMath.h"
#include "DDImage/Hash.h"

using namespace DD::Image;

static const char* const CLASS = "SyLens";
static const char* const HELP =  "This plugin undistorts footage according"
" to the lens distortion model used by Syntheyes";
static const char* const VERSION = "0.0.5";
static const char* const mode_names[] = { "undistort", "redistort", 0 };

// The number of discrete points we sample on the radius of the distortion.
// The rest is going to be interpolated
static const int STEPS = 1024;

class LutTuple
{
public:
	double r2;
	double f;
	LutTuple(double nr2, double nf2) {
		r2 = nr2;
		f = nf2;
	}
};

class SyLens : public Iop
{
	//Nuke statics
	
	const char* Class() const { return CLASS; }
	const char* node_help() const { return HELP; }
	
	
	static const Iop::Description description;
	
	Filter filter;
	
	enum { UNDIST, REDIST };
	
	// The original size of the plate that we distort
	unsigned int _plateWidth, _plateHeight;

	// The size of the filmback we are actually sampling from
	unsigned int _extWidth, _extHeight;
	
	// When we extend the image and get undistorted coordinates, we need to add these
	// values to the pixel offset
	int _paddingW, _paddingH;
	
	// Image aspect and NOT the pixel aspect Nuke furnishes us
	double _aspect;
	
	// Stuff driven by knobbz
	double kCoeff, kCubeCoeff, kUnCrop;
	bool kDbg, kTrimToFormat;
	int kMode;
	
	int _lastScanlineSize;
	
	// Used to store the output format. This needs to be kept between
	// calls to _request and cannot reside on the stack, so... 
	Format _outFormat;
	
	// Lookup table based on the radius, and we add a value outside
	// of the domain to help us interpolate
	LutTuple* _values[STEPS +1];
	
public:
	SyLens( Node *node ) : Iop ( node )
	{
		// IMPORTANT! these are the knob defaults. If someone has a plugin in
		// his script and did not change these values this is the values they will expect.
		// Therefore THESE values have to be here for the life of the plugin, through versions to come.
		// These are sane defaults - they showcase distortion and uncrop but do not use cubics.
		// cast in stone BEGIN {{
		kMode = UNDIST;
		kCoeff = -0.01826;
		kCubeCoeff = 0.0f;
		kUnCrop = 0.038f;
		kDbg = false;
		kTrimToFormat = false;
		// }} END
		
		_aspect = 1.33f;
		_lastScanlineSize = 0;
		_paddingW = 0;
		_paddingH = 0;
		
		// Initialize the lookup table with null ptrs
		for(int i = 0; i < STEPS + 1; i++) {
			_values[i] = NULL;
		}
	}
	
	~SyLens () { }
	
	void _computeAspects() {
		// Compute the aspect from the input format
		Format f = input0().format();
		
		double sc = 1.0f + (2.0f * kUnCrop);
		
		if(kMode == UNDIST) {
			_plateWidth = f.width();
			_plateHeight = f.height();
			_extWidth = ceil(float(_plateWidth) * sc);
			_extHeight = ceil(float(_plateHeight) * sc);
		} else {
			_plateWidth = floor( float(f.width()) / sc);
			_plateHeight = floor( float(f.height()) / sc);
			_extWidth = f.width();
			_extHeight = f.height();
		}
		
		_paddingW = (_extWidth - _plateWidth) / 2;
		_paddingH = (_extHeight - _plateHeight) / 2;
		_aspect = float(_plateWidth) / float(_plateHeight) *  f.pixel_aspect();
		
		if(kDbg) printf("SyLens: true plate window with uncrop will be %dx%d\n", _extWidth, _extHeight);
	}
	
	
	// Here we need to expand the image and the bounding box. This is the most important method in a plug like this so
	// pay attention
	void _validate(bool for_real)
	  {
		filter.initialize();
		input0().validate(for_real);
		copy_info();
		set_out_channels(Mask_All);
		
		info_.black_outside(true);
		
		_computeAspects();    
		
		if(kDbg) printf("SyLens: _validate info box to  %dx%d\n", _extWidth, _extHeight);
		
		// Time to define how big our output will be
		int ow, oh;
		
		if(kMode == UNDIST) {
			ow = _extWidth;
			oh = _extHeight;
		} else {
			ow = _plateWidth;
			oh = _plateHeight;
		}
		
		// Nudge outputs to power of 2, up for undistort and down for redistort
		if (kMode == UNDIST) {
			if(ow % 2 != 0) ow +=1;
			if(oh % 2 != 0) oh +=1;
		} else {
			if(ow % 2 != 0) ow -=1;
			if(oh % 2 != 0) oh -=1;
		}
		
		// Crucial. Define the format in the info_ - this is what Nuke uses
		// to know how big OUR output will be. We also pretty much NEED to store it
		// in an instance var because we cannot keep it on the stack (segfault!)
		_outFormat = Format(ow, oh, input0().format().pixel_aspect());
		info_.format( _outFormat );
		
		// For the case when we are working with a 8k by 4k plate with a SMALL CG pink elephant rrright in the left
		// corner we want to actually translate the bbox of the elephant to our distorted pipe downstream. So we need to
		// apply our SuperAlgorizm to the bbox as well and move the bbox downstream too.
		// Grab the bbox from the input first
		Info i = input0().info();
		
		if(kDbg) printf("SyLens: input info was %dx %dy %dr %dt\n", i.x(), i.y(), i.r(), i.t());
		
		Vector2 xy(i.x(), i.y());
		Vector2 tr(i.r(), i.t());
		
		// At this point we need a lookup table
		buildTable();
		
		// Here the distortion is INVERTED with relation to the pixel operation. With pixels, we need
		// to obtain the coordinate to sample FROM. However, here we need a coordinate to sample TO
		// since this is where our bbox corners are going to be in the coordinate plane of the output
		// format
		if(kMode == UNDIST) {
			undistortVectorIntoDest(xy);
			undistortVectorIntoDest(tr);
		} else {
			distortVectorIntoSource(xy);
			distortVectorIntoSource(tr);
		}
		
		Box obox(xy.x, xy.y, tr.x, tr.y);
		if(kTrimToFormat) obox.intersect(_outFormat);
		
		if(kDbg) printf("SyLens: output format will be %dx%d\n", _outFormat.width(), _outFormat.height());
		if(kDbg) printf("SyLens: output bbox is %dx%d to %dx%d\n", obox.x(), obox.y(), obox.r(), obox.t());
		
		info_.set(obox);
	}
	
	// Request the same source area upstream. We pad in the output during engine() call
	void _request(int x, int y, int r, int t, ChannelMask channels, int count)
	{
		if(kDbg) printf("SyLens: Received request %d %d %d %d\n", x, y, r, t);
		ChannelSet c1(channels); in_channels(0,c1);
		// Request the same part of the input plus padding times two. This is an opportunistic
		// cargo cult approximation but it usually allows us to grab just enough pixels to survive
		input0().request(x - (_paddingW * 2), y - (_paddingW * 2), r + (_paddingW * 2), t + (_paddingH * 2), channels, count);
	}

	// nuke
	void engine( int y, int x, int r, ChannelMask channels, Row& out );
	
	// knobs. There is really only one thing to pay attention to - be consistent and call your knobs
	// "in_snake_case_as_short_as_possible", labels are also lowercase normally
	void knobs( Knob_Callback f) {
		// For info on knob flags see Knob.h
		const int KNOB_ON_SEPARATE_LINE = 0x1000;
		
		Knob* _modeSel = Enumeration_knob(f, &kMode, mode_names, "mode", "Mode");
		_modeSel->label("mode");
		_modeSel->tooltip("Pick your poison");
		
		Knob* _kKnob = Float_knob( f, &kCoeff, "k" );
		_kKnob->label("k");
		_kKnob->tooltip("Set to the same distortion as applied by Syntheyes");
		
		Knob* _kCubeKnob = Float_knob( f, &kCubeCoeff, "kcube" );
		_kCubeKnob->label("cubic k");
		_kCubeKnob->tooltip("Set to the same cubic distortion as applied by Syntheyes");
		
		Knob* _kUncropKnob = Float_knob( f, &kUnCrop, "uncrop" );
		_kUncropKnob->label("uncrop expansion");
		_kUncropKnob->tooltip("Set to the same uncrop value as applied by Syntheyes, it will be the same on all sides");
		
		// Add the filter selection menu that comes from the filter obj itself
		filter.knobs( f );
		
		Knob* kTrimKnob = Bool_knob( f, &kTrimToFormat, "trim");
		kTrimKnob->label("trim bbox");
		kTrimKnob->tooltip("When checked, SyLens will reduce the bouding box to be within the destination format");
		kTrimKnob->set_flag(KNOB_ON_SEPARATE_LINE);
		
		Knob* kDbgKnob = Bool_knob( f, &kDbg, "debug");
		kDbgKnob->label("debug info");
		kDbgKnob->tooltip("When checked, SyLens will output various debug info to STDOUT");
		
		Divider(f, 0);
		Text_knob(f, (std::string("SyLens v.") + std::string(VERSION)).c_str());
	}
	
	// called whenever a knob is changed
	int knob_changed(Knob* k) {
		// Dont dereference unless...
		if (k == NULL) return 1;
		
		// Touching the crop knob changes our output bounds
		if ((*k).startsWith("uncrop")) {
			_validate(false);
		}
		
		// Touching the mode changes everything
		if ((*k).startsWith("mode")) {
			_validate(false);
		}
	}
	
	// Hashing for caches. We append our version to the cache hash, so that when you update
	// the plugin all the caches will(should?) be flushed automatically
	void append(Hash& hash) {
		int intv;
		memcpy(&intv, VERSION, strlen(VERSION));
		hash.append(intv);
	}
	
private:
	
	double toUv(double, int);
	double fromUv(double, int);
	void vecToUV(Vector2&, Vector2&, int, int);
	void vecFromUV(Vector2&, Vector2&, int, int);
	void distortVector(Vector2& uvVec, double k, double kcube);
	void distortVectorIntoSource(Vector2& vec);
	void undistortVectorIntoDest(Vector2& vec);
	void Remove(Vector2& vec);
	void buildTable();
	double lerp(double y1,double y2,double mu);
	double curp(double y0, double y1, double y2, double y3, double mu);
	double interpolatedF(float);
};

// Since we do not need channel selectors or masks, we can use our raw Iop
// directly instead of putting it into a NukeWrapper
static Iop* SyLensCreate( Node* node ) {
	return new SyLens(node);
}

// The second item is ignored because all a compsitor dreams of is writing fucking init.py
// every time he installs a plugin
const Iop::Description SyLens::description(CLASS, "Transform/SyLens", SyLensCreate);

// Syntheyes uses UV coordinates that start at the optical center of the image,
// and go -1,1. Nuke offers a UV option on Format that goes from 0 to 1, but it's not
// exactly what we want
double SyLens::toUv(double absValue, int absSide)
{
	double x = (absValue / (double)absSide) - 0.5;
    // .2 to -.3, y is reversed and coords are double
	return x * 2;
}

double SyLens::fromUv(double uvValue, int absSide) {
    // First, start from zero (-.1 becomes .4)
	double value_off_corner = (uvValue / 2) + 0.5f;
	return absSide * value_off_corner;
}

void SyLens::vecToUV(Vector2& absVec, Vector2& uvVec, int w, int h)
{
	// Nuke coords are 0,0 on lower left
	uvVec.x = toUv(absVec.x, w);
	uvVec.y = toUv(absVec.y, h);
}

void SyLens::vecFromUV(Vector2& absVec, Vector2& uvVec, int w, int h)
{
	// Nuke coords are 0,0 on lower left
	absVec.x = fromUv(uvVec.x, w);
	absVec.y = fromUv(uvVec.y, h);
}

// Build a table of distortion values from the 0 UV coordinate outwards to 1.
// We dd another one to helpo us interpolate since a) nobody guaranteed floats never going out
// of 1.0 b) we might use cubic or cosine interp later
void SyLens::buildTable() 
{
	double inc = 1.0f / STEPS;
	double u, v, r2, f;
	
	for (int i = 0; i < STEPS + 1; i++) {
		// Dealloc the old struct 
		if(_values[i] != NULL) delete _values[i];
		
		u = i * inc;
		r2 = _aspect * _aspect * u * u + u * u;
		
		// Skipping the square root munges up precision some
		if (fabs(kCubeCoeff) > 0.00001) {
			f = 1 + r2*(kCoeff + kCubeCoeff * sqrt(r2));
		} else {
			f = 1 + r2*(kCoeff);
		}
		_values[i] = new LutTuple(r2, f);
	}
}

// Place the UV coordinates of the distorted image in the undistorted image plane into uvVec
void SyLens::distortVector(Vector2& uvVec, double k, double kcube)
{
	
	// The radius on which computation is based
	float r2 = _aspect * _aspect * uvVec.x * uvVec.x + uvVec.y * uvVec.y;
	
	// The distortion multiplier derived from the radius.
	float f = interpolatedF(r2);
	
	uvVec.x = uvVec.x * f;
	uvVec.y = uvVec.y * f;
}

double SyLens::interpolatedF(float r2)
{
	float fLeft;
	float fRight;
	float rLeft;
	float rRight;
	
	// find the nearest values in the table. We stashed an extreme
	// value at the end of _values as well
	for (int i = 0; i < STEPS + 1; i++) {
		if(_values[i]->r2 < r2) {
			rLeft = _values[i]->r2;
			fLeft = _values[i]->f;
		}
		
		if(_values[i]->r2 > r2) {
			rRight = _values[i]->r2;
			fRight = _values[i]->f;
			break;
		}
	}
	
	// linearly interpolate between them (no not smart but does the jobb)
	return lerp(fLeft, fRight, r2 - rLeft);
	
}
// Do a linear intepolation
double SyLens::lerp(
   double y1,double y2,
   double mu)
{
   return(y1*(1-mu)+y2*mu);
}

// Do a cubic interpolation
double SyLens::curp(
   double y0,double y1,
   double y2,double y3,
   double mu)
{
   double a0,a1,a2,a3,mu2;

   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;

   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}
	
// Get a coordinate that we need to sample from the SOURCE distorted image to get at the absXY
// values in the RESULT
void SyLens::distortVectorIntoSource(Vector2& absXY) {
	Vector2 uvXY(0, 0);
	// The gritty bits - get coordinates of the distorted pixel in the coordinates of the
	// EXTENDED film back
	vecToUV(absXY, uvXY, _extWidth, _extHeight);
	distortVector(uvXY, kCoeff, kCubeCoeff);
	vecFromUV(absXY, uvXY, _extWidth, _extHeight);
	absXY.x -= (double)_paddingW;
	absXY.y -= (double)_paddingH;
}

// This is still a little wrongish but less wrong than before
void SyLens::undistortVectorIntoDest(Vector2& absXY) {
	absXY.x += (double)_paddingW;
	absXY.y += (double)_paddingH;
	Vector2 uvXY(0, 0);
	vecToUV(absXY, uvXY, _extWidth, _extHeight);
	Remove(uvXY);
	vecFromUV(absXY, uvXY, _extWidth, _extHeight);
}

// THIS IS SEMI-WRONG! Ported over from distort.szl, does not honor cubic distortion
void SyLens::Remove(Vector2& pt) {
	double rp, f;
	
	rp = _aspect * _aspect * pt.x * pt.x + pt.y * pt.y;
	
	// The distortion multiplier derived from the radius.
	f = interpolatedF(rp);
	
	pt.x = pt.x / f;
	pt.y = pt.y / f;
}


// The image processor that works by scanline. Y is the scanline offset, x is the pix,
// r is the length of the row. We are now effectively in the undistorted coordinates, mind you!
void SyLens::engine ( int y, int x, int r, ChannelMask channels, Row& out )
{
	if(r != _lastScanlineSize) {
		if(kDbg) printf("SyLens: Rendering scanline %d pix\n", r);
		_lastScanlineSize = r;
	}
	
	foreach(z, channels) out.writable(z);
	
	Pixel pixel(channels);
	const float sampleOff = 0.5f;
	
	Vector2 sampleFromXY(0.0f, 0.0f);
	for (; x < r; x++) {
		
		sampleFromXY = Vector2(x, y);
		if( kMode == UNDIST) {
			distortVectorIntoSource(sampleFromXY);
		} else {
			undistortVectorIntoDest(sampleFromXY);
		}
		
		// Sample from the input node at the coordinates
		// half a pixel has to be added here because sample() takes the first two
		// arguments as the center of the rectangle to sample. By not adding 0.5 we'd
		// have to deal with a slight offset which is *not* desired.
		input0().sample(
			sampleFromXY.x + sampleOff , sampleFromXY.y + sampleOff, 
			1.0f, 
			1.0f,
			&filter,
			pixel
		);
		
		// write the resulting pixel into the image
		foreach (z, channels)
		{
			((float*)out[z])[x] = pixel[z];
		}
	}
}
