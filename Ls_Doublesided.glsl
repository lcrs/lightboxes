// Light front faces one colour, back faces another
// lewis@lewissaunders.com

vec3 adsk_getCameraPosition();
vec3 adsk_getNormal();
vec3 adsk_getVertexPosition();
vec3 adsk_getLightColour();
vec3 adsk_hsv2rgb(vec3 hsv);

uniform vec3 adskUID_cfront, adskUID_cback;

vec4 adskUID_lightbox(vec4 i) {
	vec3 v = normalize(adsk_getCameraPosition() - adsk_getVertexPosition());
	vec3 n = normalize(adsk_getNormal());
	float facing = dot(n, v);
	vec3 hsv;
	if(facing > 0.0) {
		hsv = adskUID_cfront;
	} else {
		hsv = adskUID_cback;
	}
	hsv.x /= 360.0;
	hsv.y /= 100.0;
	hsv.z /= 100.0;
	vec3 c = adsk_hsv2rgb(hsv);
	vec3 cmul = i.rgb * c;
	i.rgb = cmul * i.a;
	return i;
}
