vec3 adsk_getDiffuseMapCoord();
vec4 adsk_getDiffuseMapValue(in vec2 texCoord);

uniform vec2 adskUID_offset;
uniform bool adskUID_showarea;

vec4 adskUID_lightbox(vec4 i) {
	vec3 uv = adsk_getDiffuseMapCoord();
	uv.xy /= max(uv.z, 0.0001);

	uv.xy -= adskUID_offset * i.a;
	i.rgb = adsk_getDiffuseMapValue(uv.xy).rgb;

	if(adskUID_showarea) i.rgb += i.a * vec3(0.5, -0.5, 0.5);

	i.a = 1.0;
	return i;
}
