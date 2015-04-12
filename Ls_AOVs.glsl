// Returns various arbitrary outputs
// lewis@lewissaunders.com

vec3 adsk_getVertexPosition();
vec3 adsk_getNormal();
vec3 adsk_getComputedNormal();
vec3 adsk_getBinormal();
vec3 adsk_getTangent();
vec3 adsk_getLightPosition();
vec3 adsk_getLightDirection();
vec3 adsk_getLightTangent();
vec3 adsk_getLightColour();
float adsk_getLightAlpha();
vec3 adsk_getComputedDiffuse();
vec3 adsk_getComputedSpecular();
vec4 adsk_getComputedDiffuseMapValue(in vec3 vertexPos);
vec3 adsk_getDiffuseMapCoord();
vec3 adsk_getEmissiveMapCoord();
vec3 adsk_getSpecularMapCoord();
vec3 adsk_getNormalMapCoord();
vec2 adsk_getReflectionMapCoord(in vec3 vrtPos, in vec3 normal);
vec3 adsk_getParallaxMapCoord();
vec4 adsk_getDiffuseMapValue(in vec2 texCoord);
vec4 adsk_getEmissiveMapValue(in vec2 texCoord);
vec4 adsk_getSpecularMapValue(in vec2 texCoord);
vec4 adsk_getNormalMapValue(in vec2 texCoord);
vec4 adsk_getReflectionMapValue( in vec2 texCoord);
vec4 adsk_getParallaxMapValue(in vec2 texCoord);
vec4 adsk_getUVMapValue(in vec2 texCoord);
vec4 adsk_getMaterialDiffuse();
vec4 adsk_getMaterialSpecular();
vec4 adsk_getMaterialAmbient();
float adsk_getShininess();
vec3 adsk_getComputedShadowCoefficient();

uniform int w;

vec2 perspectiveDivide(vec3 v) {
	return v.xy / v.z;
}

vec4 adskUID_lightbox(vec4 i) {
	if(w==0) return vec4(adsk_getVertexPosition(), 1.0);
	if(w==1) return vec4(adsk_getNormal(), 1.0);
	if(w==2) return vec4(adsk_getComputedNormal(), 1.0);
	if(w==3) return vec4(adsk_getBinormal(), 1.0);
	if(w==4) return vec4(adsk_getTangent(), 1.0);
	if(w==5) return vec4(adsk_getLightPosition(), 1.0);
	if(w==6) return vec4(adsk_getLightDirection(), 1.0);
	if(w==7) return vec4(adsk_getLightTangent(), 1.0);
	if(w==8) return vec4(adsk_getLightColour(), 1.0);
	if(w==9) return vec4(vec3(adsk_getLightAlpha()), 1.0);
	if(w==10) return vec4(adsk_getComputedDiffuse(), 1.0);
	if(w==11) return vec4(adsk_getComputedSpecular(), 1.0);
	if(w==12) return vec4(adsk_getComputedDiffuseMapValue(adsk_getVertexPosition()).rgb, 1.0);
	if(w==13) return vec4(adsk_getDiffuseMapCoord(), 1.0);
	if(w==14) return vec4(adsk_getEmissiveMapCoord(), 1.0);
	if(w==15) return vec4(adsk_getSpecularMapCoord(), 1.0);
	if(w==16) return vec4(adsk_getNormalMapCoord(), 1.0);
	if(w==17) return vec4(adsk_getReflectionMapCoord(adsk_getVertexPosition(), adsk_getComputedNormal()), 0.0, 1.0);
	if(w==18) return vec4(adsk_getParallaxMapCoord(), 1.0);
	if(w==19) return vec4(adsk_getDiffuseMapValue(perspectiveDivide(adsk_getDiffuseMapCoord())).rgb, 1.0);
	if(w==20) return vec4(adsk_getEmissiveMapValue(perspectiveDivide(adsk_getEmissiveMapCoord())).rgb, 1.0);
	if(w==21) return vec4(adsk_getSpecularMapValue(perspectiveDivide(adsk_getSpecularMapCoord())).rgb, 1.0);
	if(w==22) return vec4(adsk_getNormalMapValue(perspectiveDivide(adsk_getNormalMapCoord())).rgb, 1.0);
	if(w==23) return vec4(adsk_getReflectionMapValue(adsk_getReflectionMapCoord(adsk_getVertexPosition(), adsk_getComputedNormal())).rgb, 1.0);
	if(w==24) return vec4(adsk_getParallaxMapValue(perspectiveDivide(adsk_getParallaxMapCoord())).rgb, 1.0);
	if(w==25) return vec4(adsk_getUVMapValue(perspectiveDivide(adsk_getDiffuseMapCoord())).rgb, 1.0);
	if(w==26) return vec4(adsk_getMaterialDiffuse().rgb, 1.0);
	if(w==27) return vec4(adsk_getMaterialSpecular().rgb, 1.0);
	if(w==28) return vec4(adsk_getMaterialAmbient().rgb, 1.0);
	if(w==29) return vec4(vec3(adsk_getShininess()), 1.0);
	if(w==30) return vec4(adsk_getComputedShadowCoefficient().rgb, 1.0);
	if(w==31) return vec4(i.rgb, 1.0);
	if(w==32) return vec4(i.aaa, 1.0);
}
