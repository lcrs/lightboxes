// Tint faces which face the light
// lewis@lewissaunders.com

vec3 adsk_getLightPosition();
vec3 adsk_getNormal();
vec3 adsk_getVertexPosition();
vec3 adsk_getLightColour();
vec3 adsk_hsv2rgb(vec3 hsv);

uniform vec3 adskUID_col;

vec4 adskUID_lightbox(vec4 i) {
	vec3 l = normalize(adsk_getLightPosition() - adsk_getVertexPosition());
	vec3 n = normalize(adsk_getNormal());
	float facing = dot(n, l);
	facing = clamp(facing, 0.0, 1.0);
	
	vec3 hsv = adskUID_col;
	hsv.x /= 360.0;
	hsv.y /= 100.0;
	hsv.z /= 100.0;
	vec3 rgbcol = adsk_hsv2rgb(hsv);
	vec3 fulltint = i.rgb * rgbcol;
	i.rgb = mix(i.rgb, fulltint, facing * i.a * adsk_getLightColour());
	return i;
}
