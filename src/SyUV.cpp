/*

This node implements the first part of the workflow described here
http://forums.thefoundry.co.uk/phpBB2/viewtopic.php?p=1915&sid=497b22394ab62d0bf8bd7a35fa229913

*/

static const char* const CLASS = "SyUV";
static const char* const HELP = "This node will apply the Syntheyes undistortion in the UV space of your geometry. "
	"Use ProjectUV or just a card to get a basic geo, and apply this node - your UVs will be undistorted. "
	"Contact me@julik.nl if you need help with the plugin.";

#include "VERSION.h"

#include "DDImage/ModifyGeo.h"
#include "DDImage/Scene.h"
#include "DDImage/Knob.h"
#include "DDImage/Knobs.h"
#include <sstream>
#include "SyDistorter.cpp"

using namespace DD::Image;

class SyUV : public ModifyGeo
{
private:
	
	const char* uv_attrib_name;
	
	// The distortion engine
	SyDistorter distorter;
	
public:

	static const Description description;
	
	const char* Class() const
	{
		return CLASS;
	}
	
	const char* node_help() const
	{
		return HELP;
	}

	SyUV(Node* node) : ModifyGeo(node)
	{
		uv_attrib_name = "uv";
	}
	
	void append(Hash& hash) {
		hash.append(VERSION);
		
		// Knobs that change the SyLens algo
		hash.append(distorter.compute_hash());
		hash.append(uv_attrib_name);
		ModifyGeo::append(hash); // the super called he wants his pointers back
	}
	
	void knobs(Knob_Callback f)
	{
		ModifyGeo::knobs(f);
		distorter.knobs_with_aspect(f);
		Divider(f, 0);
		std::ostringstream ver;
		ver << "SyUV v." << VERSION;
		Text_knob(f, ver.str().c_str());
	}

	void get_geometry_hash()
	{
		// Get all hashes up-to-date
		ModifyGeo::get_geometry_hash();
		
		// Knobs that change the SyLens algo:
		geo_hash[Group_Points].append(distorter.compute_hash());
		geo_hash[Group_Points].append(uv_attrib_name);
	}
	
	void _validate(bool for_real)
	{
		distorter.recompute_if_needed();
		return ModifyGeo::_validate(for_real);
	}
	
	// This is needed to preserve UVs which are already there
	void keep_uvs(int index, GeoInfo& info, GeometryList& out)
	{
		
		// get the original uv attribute used to restore untouched uv coordinate
		const AttribContext* context = info.get_attribcontext(uv_attrib_name);
		AttributePtr uv_original = context ? context->attribute : AttributePtr();

		if(!uv_original){
			Op::error( "Missing \"%s\" channel from geometry", uv_attrib_name );
			return;
		}

		DD::Image::GroupType t_group_type = context->group; // texture coordinate group type

		// we have two possibilities:
		// the uv coordinate are stored in Group_Points or in Group_Vertices way
		// sanity check
		assert(t_group_type == Group_Points || t_group_type == Group_Vertices);

		// create a buffer to write on it
		Attribute* uv = out.writable_attribute(index, t_group_type, uv_attrib_name, VECTOR4_ATTRIB);
		assert(uv);

		// copy all original texture coordinate if available
		if (uv_original){

			// sanity check
			assert(uv->size() == uv_original->size());

			for (unsigned i = 0; i < uv->size(); i++) {
				uv->vector4(i) = uv_original->vector4(i);
			}
		}
	}
	
	void modify_geometry(int obj, Scene& scene, GeometryList& out)
	{
		const char* uv_attrib_name = "uv";
		// Call the engine on all the caches:
		for (unsigned i = 0; i < out.objects(); i++) {
			GeoInfo& info = out[i];
			
			// Copy over old UV attributes
			keep_uvs(i, info, out);
			
			// TODO: investigate difference between vertex and point UVs
			
			// Create a point attribute
			Attribute* uv = out.writable_attribute(i, Group_Points, uv_attrib_name, VECTOR4_ATTRIB);
			if(!uv) return;
			
			for (unsigned p = 0; p < info.points(); p++) {
				distorter.distort_uv(uv->vector4(p));
			}
		}
	}
};

static Op* build(Node* node)
{
	return new SyUV(node);
}
const Op::Description SyUV::description(CLASS, build);

