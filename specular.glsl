vec4 adsk_getMaterialSpecular();

vec4 adskUID_lightbox(vec4 i) {
	i.rgb = adsk_getMaterialSpecular().rgb;
	return i;
}
