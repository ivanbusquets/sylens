/*
	This plugin uses the lens distortion model provided by Russ Anderson of the SynthEyes camera tracker
	fame. It is made so that it's output is _identical_ to the Image Preparation tool of the SY camera tracker.
	so you can use it INSTEAD of the footage Syntheyes prerenders.
	
	It implements the algorithm described here http://www.ssontech.com/content/lensalg.htm
	
	It is largely based on tx_nukeLensDistortion by Matti Gruener of Trixter Film and Rising Sun Films.
	The code has however been simplified and some features not present in the original version have been added.
	
	Written by Julik Tarkhanov in Amsterdam in 2010-2011 with kind support by HecticElectric.
	I thank the users for their continued support and bug reports.
	For questions mail me(at)julik.nl
	
	The beautiful Crimean landscape shot used in the test script is provided by Tim Parshikov
	and Mikhail Mestezky, 2010.
	
	The code has some more comments than it's 3DE counterpart since we have to do some things that the other plugin
	did not
*/

// for our friend printf
extern "C" {
#include <stdio.h>
}

// For max/min on containers
#include <algorithm>

// For string concats
#include <sstream>

#include "DDImage/Iop.h"
#include "DDImage/Row.h"
#include "DDImage/Pixel.h"
#include "DDImage/Filter.h"
#include "DDImage/Knobs.h"
#include "SyDistorter.cpp"

using namespace DD::Image;

static const char* const CLASS = "SyLens";
static const char* const HELP =  "This plugin undistorts footage according "
	"to the lens distortion model used by Syntheyes. "
	"Contact me@julik.nl if you need help with the plugin.";

#include "VERSION.h"

static const char* const mode_names[] = { "undistort", "redistort", 0 };


class SyLens : public Iop
{
	//Nuke statics
	
	const char* Class() const { return CLASS; }
	const char* node_help() const { return HELP; }
	static const Iop::Description description;
	
	Filter filter;
	
	enum { UNDIST, REDIST };
	
	// The original size of the plate that we distort
	unsigned int plate_width_, plate_height_;

	// Image aspect and NOT the pixel aspect Nuke furnishes us
	double _aspect;
	
	// Movable centerpoint offsets
	double centerpoint_shift_u_, centerpoint_shift_v_;
	
	// Stuff driven by knobbz
	bool k_enable_debug_, k_trim_bbox_to_format_, k_only_format_output_;
	int kMode;
	
	// The distortion engine
	SyDistorter distorter;
	
public:
	SyLens( Node *node ) : Iop ( node )
	{
		kMode = UNDIST;
		_aspect = 1.33f;
	}
	
	void _computeAspects();
	void _validate(bool for_real);
	void _request(int x, int y, int r, int t, ChannelMask channels, int count);
	void engine( int y, int x, int r, ChannelMask channels, Row& out );
	void knobs( Knob_Callback f);
	int knob_changed(Knob* k);
	
	// Hashing for caches. We append our version to the cache hash, so that when you update
	// the plugin all the caches will be flushed automatically
	void append(Hash& hash) {
		hash.append(VERSION);
		distorter.append(hash);
		Iop::append(hash); // the super called he wants his pointers back
	}
	
	~SyLens () { 
	}
	
private:
	
	int round(double x);
	double toUv(double, int);
	double fromUv(double, int);
	void absolute_px_to_centered_uv(Vector2&, int, int);
	void centered_uv_to_absolute_px(Vector2&, int, int);
	void distort_px_into_source(Vector2& vec);
	void undistort_px_into_destination(Vector2& vec);
};

// Since we do not need channel selectors or masks, we can use our raw Iop
// directly instead of putting it into a NukeWrapper. Besides, the mask input on
// the NukeWrapper cannot be disabled even though the Foundry doco says it can
// (Foundry bug #12598)
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
	double x = (absValue / (double)absSide) - 0.5f;
	return x * 2;
}

double SyLens::fromUv(double uvValue, int absSide) {
	double value_off_corner = (uvValue / 2) + 0.5f;
	return absSide * value_off_corner;
}

void SyLens::absolute_px_to_centered_uv(Vector2& xy, int w, int h)
{
	// Nuke coords are 0,0 on lower left
	xy.x = toUv(xy.x, w);
	xy.y = toUv(xy.y, h);
}

void SyLens::centered_uv_to_absolute_px(Vector2& xy, int w, int h)
{
	// Nuke coords are 0,0 on lower left
	xy.x = fromUv(xy.x, w);
	xy.y = fromUv(xy.y, h);
}

// Get a coordinate that we need to sample from the SOURCE distorted image to get at the absXY
// values in the RESULT
void SyLens::distort_px_into_source(Vector2& absXY) {
	absolute_px_to_centered_uv(absXY, plate_width_, plate_height_);
	distorter.apply_disto(absXY);
	centered_uv_to_absolute_px(absXY, plate_width_, plate_height_);
}

// This is still a little wrongish but less wrong than before
void SyLens::undistort_px_into_destination(Vector2& absXY) {
	absolute_px_to_centered_uv(absXY, plate_width_, plate_height_);
	distorter.remove_disto(absXY);
	centered_uv_to_absolute_px(absXY, plate_width_, plate_height_);
}

// The image processor that works by scanline. Y is the scanline offset, x is the pix,
// r is the length of the row. We are now effectively in the undistorted coordinates, mind you!
void SyLens::engine ( int y, int x, int r, ChannelMask channels, Row& out )
{
	
	foreach(z, channels) out.writable(z);
	
	Pixel pixel(channels);
	const float sampleOff = 0.5f;
	
	Vector2 sampleFromXY(0.0f, 0.0f);
	for (; x < r; x++) {
		
		sampleFromXY = Vector2(x, y);
		if( kMode == UNDIST) {
			distort_px_into_source(sampleFromXY);
		} else {
			undistort_px_into_destination(sampleFromXY);
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


// knobs. There is really only one thing to pay attention to - be consistent and call your knobs
// "in_snake_case_as_short_as_possible", labels are also lowercase normally
void SyLens::knobs( Knob_Callback f) {
	// For info on knob flags see Knob.h
	const int KNOB_ON_SEPARATE_LINE = 0x1000;
	const int KNOB_HIDDEN = 0x0000000000040000;
	
	Knob* _modeSel = Enumeration_knob(f, &kMode, mode_names, "mode", "Mode");
	_modeSel->label("mode");
	_modeSel->tooltip("Pick your poison");
	
	distorter.knobs(f);
	filter.knobs(f);
	
	Divider(f, 0);
	
	std::ostringstream ver;
	ver << "SyLens v." << VERSION;
	Text_knob(f, ver.str().c_str());
}

// called whenever a knob is changed
int SyLens::knob_changed(Knob* k) {
	
	// Touching the crop knob changes our output bounds
	if (k->startsWith("uncrop")) {
		_validate(false);
	}
	
	// Touching the mode changes everything
	if (k->startsWith("mode")) {
		_validate(false);
	}
	
	return Iop::knob_changed(k); // Super knows better
}

// http://stackoverflow.com/questions/485525/round-for-float-in-c
int SyLens::round(double x) {
	return floor(x + 0.5);
}

// The algo works in image aspec, not the pixel aspect. We also have to take the uncrop factor
// into account.
void SyLens::_computeAspects() {
	// Compute the aspect from the input format
	Format f = input0().format();
	
	plate_width_ = round(f.width());
	plate_height_ = round(f.height());

	_aspect = float(plate_width_) / float(plate_height_) *  f.pixel_aspect();
	
	if(k_enable_debug_) printf("SyLens: true plate window with uncrop will be %dx%d\n", plate_width_, plate_height_);
}

	
// Here we need to expand the image and the bounding box. This is the most important method in a plug like this so
// pay attention
void SyLens::_validate(bool for_real)
{
	
	// Set the autolabel
	Knob* al = knob("label");
	std::string label;
	if (kMode == UNDIST) {
		label = "undistort";
	} else {
		label = "redistort";
	}
	
	al->set_text(label.c_str());
	
	// Bookkeeping boilerplate
	filter.initialize();
	input0().validate(for_real);
	copy_info();
	set_out_channels(Mask_All);
	
	// Do not blank away everything
	info_.black_outside(false);
	
	// We need to know our aspects so prep them here
	_computeAspects();
	
	distorter.set_aspect(_aspect);
	distorter.recompute_if_needed();
	
	if(k_enable_debug_) printf("SyLens: _validate info box to  %dx%d\n", plate_width_, plate_height_);
	
	// Time to define how big our output will be in terms of format. Format will always be the whole plate.
	// If we only use a bboxed piece of the image we will limit our request to that. But first of all we need to
	// compute the format of our output.
	int ow, oh;
	
	ow = plate_width_;
	oh = plate_height_;
	
	// Nudge outputs to power of 2, upwards
	if(ow % 2 != 0) ow +=1;
	if(oh % 2 != 0) oh +=1;
	
	// For the case when we are working with a 8k by 4k plate with a SMALL CG pink elephant rrright in the left
	// corner we want to actually translate the bbox of the elephant to our distorted pipe downstream. So we need to
	// apply our SuperAlgorizm to the bbox as well and move the bbox downstream too.
	// Grab the bbox from the input first
	Info inf = input0().info();
	Format f = input0().format();
	
	// Just distorting the four corners of the bbox is NOT enough. We also need to find out whether
	// the bbox intersects the centerlines. Since the distortion is the most extreme at the centerlines if
	// we just take the corners we might be chopping some image away. So to get a reliable bbox we need to check
	// our padding at 6 points - the 4 extremes and where the bbox crosses the middle of the coordinates
	int xMid = ow/2;
	int yMid = oh/2;

	std::vector<Vector2*> pointsOnBbox;
	
	// Add the standard two points - LR and TR
	pointsOnBbox.push_back(new Vector2(inf.x(), inf.y()));
	pointsOnBbox.push_back(new Vector2(inf.r(), inf.t()));
	
	// Add the TL and LR as well
	pointsOnBbox.push_back(new Vector2(inf.x(), inf.t()));
	pointsOnBbox.push_back(new Vector2(inf.r(), inf.y()));
	
	// If our box intersects the midplane on X add the points where the bbox crosses centerline
	if((inf.x() < xMid) && (inf.r() > xMid)) {
		// Find the two intersections and add them
		pointsOnBbox.push_back( new Vector2(xMid, inf.y()) );
		pointsOnBbox.push_back( new Vector2(xMid, inf.t()) );
	}
	
	// If our box intersects the midplane on Y add the points where the bbox crosses centerline
	if((inf.y() < yMid) && (inf.t() > yMid)) {
		pointsOnBbox.push_back( new Vector2(inf.x(), yMid) );
		pointsOnBbox.push_back( new Vector2(inf.r(), yMid) );
	}
	
	std::vector<int> xValues;
	std::vector<int> yValues;
	
	// Here the distortion is INVERTED with relation to the pixel operation. With pixels, we need
	// to obtain the coordinate to sample FROM. However, here we need a coordinate to sample TO
	// since this is where our bbox corners are going to be in the coordinate plane of the output
	// format.
	for(unsigned int i = 0; i < pointsOnBbox.size(); i++) {
		if(kMode == UNDIST) {
			undistort_px_into_destination(*pointsOnBbox[i]);
		} else {
			distort_px_into_source(*pointsOnBbox[i]);
		}
		xValues.push_back(pointsOnBbox[i]->x);
		yValues.push_back(pointsOnBbox[i]->y);
	}
	
	int minX, minY, maxX, maxY;
	
	// Formally speaking, we have to allocate an std::iterator first. But we wont.
	minX = *std::min_element(xValues.begin(), xValues.end());
	maxX = *std::max_element(xValues.begin(), xValues.end());
	minY = *std::min_element(yValues.begin(), yValues.end());
	maxY = *std::max_element(yValues.begin(), yValues.end());
	
	Box obox(minX, minY, maxX, maxY);
	
	// If trim is enabled we intersect our obox with the format so that there is no bounding box
	// outside the crop area. Thiis handy for redistorted material.
	if(k_trim_bbox_to_format_) obox.intersect(input0().format());
	
	if(k_enable_debug_) printf("SyLens: output bbox is %dx%d to %dx%d\n", obox.x(), obox.y(), obox.r(), obox.t());
	
	info_.set(obox);
}

void SyLens::_request(int x, int y, int r, int t, ChannelMask channels, int count)
{
	
	if(k_enable_debug_) printf("SyLens: Received request %d %d %d %d\n", x, y, r, t);
	ChannelSet c1(channels); in_channels(0,c1);
	
	Vector2 bl(x, y), br(r, y), tr(r, t), tl(x, t);
	
	if(kMode == UNDIST) {
		distort_px_into_source(bl);
		distort_px_into_source(br);
		distort_px_into_source(tl);
		distort_px_into_source(tr);
	} else {
		undistort_px_into_destination(bl);
		undistort_px_into_destination(br);
		undistort_px_into_destination(tl);
		undistort_px_into_destination(tr);
	}
	
	// Request the same part of the input distorted. However if rounding errors have taken place 
	// it is possible that in engine() we will need to sample from the pixels slightly outside of this area.
	// If we don't request it we will get black pixels in there, so we add a small margin on all sides
	// to give us a little cushion
	const unsigned int safetyPadding = 16;
	input0().request(
		round(bl.x - safetyPadding),  
		round(bl.y  - safetyPadding),
		round(tr.x + safetyPadding),
		round(tr.y + safetyPadding),
		channels, count
	);
}
