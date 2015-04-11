vec3 adsk_getCameraPosition();
vec3 adsk_getVertexPosition();
vec3 adsk_getComputedNormal();
vec3 adsk_getLightPosition();
vec3 adsk_getLightColour();
vec3 adskEvalDynCurves(ivec3 curve, vec3 x);

uniform bool adskUID_usecam;
uniform ivec3 adskUID_curves;

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
	vec3 light = -adsk_getLightPosition();
	vec3 p = adsk_getVertexPosition();
	vec3 n = adsk_getComputedNormal();
	vec3 v = normalize(cam - p);
	vec3 l = normalize(light - p);

	float lambert;
	if(adskUID_usecam) {
		lambert = dot(n, v);
	} else {
		lambert = dot(n, l);
	}

	vec3 col = adskEvalDynCurves(adskUID_curves, vec3(lambert));
	col *= i.rgb * adsk_getLightColour();
	i.rgb = mix(i.rgb, col, i.a);
	return i;
}
