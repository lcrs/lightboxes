// Correct IBL reflections for latlong maps
// lewis@lewissanders.com

vec3 adsk_getCameraPosition();
vec3 adsk_getVertexPosition();
vec3 adsk_getComputedNormal();
mat4 adsk_getIBLRotationMatrix(in int idx);
vec3 adsk_getAngularMapIBL(in int idx, in vec2 coords, float lod);

uniform int adskUID_iblidx;
uniform bool adskUID_correct;
uniform float adskUID_lod;

const float adskUID_PI = 3.14159265358979323846;

// This function is from the Action ubershader and seems wrong
vec2 adskUID_getCylinderMappingTexCoord( in vec3 transfRefl ) {
   vec3 normTransfRefl = normalize( transfRefl );
   
   // compute the rotation angle between x/z on the y axis
   const float PI = 3.14159265358979;
   const float factor = 1.0 / ( 2.0 * PI );
   const float offset = -0.25;
   float x = atan( normTransfRefl.z, normTransfRefl.x ) * factor + offset;

   // make sure x isn't out of the 0..1 range and scale y from -1..1 to 0..1
   return vec2( fract( x ), normTransfRefl.y * 0.5 + 0.5 );
}

// This one matches Houdini's reflection of the same environment
vec2 adskUID_latlong(vec3 v) {
	float lat = asin(v.y) / adskUID_PI + 0.5;
	float lon = atan(v.z, v.x) / (2.0*adskUID_PI) + 0.75;
	return fract(vec2(lon, lat));
}

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
	vec3 vertex = adsk_getVertexPosition();
	vec3 view = normalize(vertex - cam);
	vec3 n = normalize(adsk_getComputedNormal());
	vec3 refl = reflect(view, n);

	refl = (adsk_getIBLRotationMatrix(adskUID_iblidx) * vec4(refl, 0.0)).xyz;

	vec2 st;
	if(adskUID_correct) {
		st = adskUID_latlong(refl);
	} else {
		st = adskUID_getCylinderMappingTexCoord(refl);
	}
	vec3 sample = adsk_getAngularMapIBL(adskUID_iblidx, st, adskUID_lod);

	return vec4(sample, 1.0);
}
