#line 1
// Correct IBL reflection mapping
// lewis@lewissanders.com

const float adskUID_PI = 3.14159265358979323846;

// This function is from the Action ubershader and seems wrong
// convert a 3d vector into the 2d coord (uv) to index a cylinder map
vec2 adskUID_getCylinderMappingTexCoord( in vec3 transfRefl )
{
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
    vec3 p = adsk_getVertexPosition();
    vec3 v = normalize(cam - p);
    vec3 n = adsk_getNormal();

    return vec4(5.0 * adsk_getAngularMapIBL(0, adskUID_latlong(reflect(-v, n)), 0.0), 1.0);

}
