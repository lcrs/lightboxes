// UVtweak - adjust UVs to warp diffuse texture in region of light
// lewis@lewissaunders.com
vec3 adsk_getDiffuseMapCoord();
vec4 adsk_getDiffuseMapValue(in vec2 texCoord);

uniform vec2 adskUID_offset;
uniform bool adskUID_showarea;

vec4 adskUID_lightbox(vec4 i) {
	vec3 uv = adsk_getDiffuseMapCoord();
	uv.xy /= max(uv.z, 0.0001);
	uv.xy -= adskUID_offset * i.a;
	vec3 new = adsk_getDiffuseMapValue(uv.xy).rgb;

	if(i.a < 0.000001) {
		// Outside area of influence.  Might as well try to return current lighting result...
		new = i.rgb;
	}

	if(adskUID_showarea) {
		// Tint
		new.rgb += i.a * vec3(0.5, -0.5, 0.5);
	}

	return vec4(new, 1.0);
}
