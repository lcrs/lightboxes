vec3 adsk_getLightPosition();
vec3 adsk_getNormal();
vec3 adsk_getVertexPosition();

vec4 adskUID_lightbox(vec4 i) {
	vec3 l = normalize(-adsk_getLightPosition() - adsk_getVertexPosition());
	vec3 n = normalize(adsk_getNormal());
	float facing = dot(n, l);

	i.rgb = n;
	i.a = 1.0;

	return i;
}
