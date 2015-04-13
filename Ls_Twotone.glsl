// Uses RGB curves to adjust colour of faces which point toward the light vs. those which point away
// lewis@lewissaunders.com

vec3 adsk_getCameraPosition();
vec3 adsk_getVertexPosition();
vec3 adsk_getComputedNormal();
vec3 adsk_getLightPosition();
vec3 adsk_getLightColour();
vec3 adskEvalDynCurves(ivec3 curve, vec3 x);
vec4 adsk_getBlendedValue(int blendType, vec4 srcColor, vec4 dstColor);

uniform bool adskUID_usecam;
uniform bool adskUID_mono;
uniform int adskUID_blend;
uniform ivec3 adskUID_curves;

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
	vec3 light = adsk_getLightPosition();
	vec3 p = adsk_getVertexPosition();
	vec3 n = normalize(adsk_getComputedNormal());
	vec3 v = normalize(cam - p);
	vec3 l = normalize(light - p);

	float lambert;
	if(adskUID_usecam) {
		lambert = dot(n, v);
	} else {
		lambert = dot(n, l);
	}

	vec3 col = adskEvalDynCurves(adskUID_curves, vec3(lambert));

	if(adskUID_mono) col = col.rrr;
	col *= adsk_getLightColour();

	col = adsk_getBlendedValue(adskUID_blend, vec4(i.rgb, 1.0), vec4(col, 1.0)).rgb;
	col = mix(i.rgb, col, i.a);

	return vec4(col, 1.0);
}
