vec3 adsk_getCameraPosition();
vec3 adsk_getVertexPosition();
vec3 adsk_getComputedNormal();
vec3 adsk_getLightPosition();
vec3 adsk_getLightColour();
float adskEvalDynCurves(int curve, float x);

uniform bool adskUID_oppsitecam;
uniform int adskUID_wrap;

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
	vec3 light = adsk_getLightPosition();
	vec3 p = adsk_getVertexPosition();
	vec3 n = adsk_getComputedNormal();
	vec3 v = normalize(cam - p);
	vec3 l = normalize(light - p);

	float lambert;
	if(adskUID_oppsitecam) {
		lambert = dot(n, v);
	} else {
		lambert = dot(n, l);
	}

	float backlight = 1.0 - lambert;
	backlight = adskEvalDynCurves(adskUID_wrap, backlight);

	i.rgb += backlight * i.a * adsk_getLightColour();
	return i;
}
