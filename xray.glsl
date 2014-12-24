 vec3 adsk_getNormal();							 vec3 adsk_getVertexPosition();
 vec3 adsk_getBinormal();						 vec3 adsk_getTangent();
 vec3 adsk_getCameraPosition();					float adsk_getTime();

 bool adsk_isLightActive();						float adsk_getLightAlpha();
 bool adsk_isLightAdditive();					 bool adsk_isPointSpotLight();
 vec3 adsk_getLightPosition();					 bool adsk_isDirectionalLight();
 vec3 adsk_getLightColour();					 bool adsk_isAmbientLight();
 vec3 adsk_getLightDirection();					 bool adsk_isAreaRectangleLight();
 vec3 adsk_getLightTangent();					 bool adsk_isAreaEllipseLight();
float adsk_getLightDecayRate();					float adsk_getAreaLightWidth();
 bool adsk_isSpotlightFalloffParametric();		float adsk_getAreaLightHeight();
float adsk_getSpotlightParametricFalloffIn();	float adsk_getSpotlightSpread();
float adsk_getSpotlightParametricFalloffOut();	

 bool adsk_isLightboxRenderedFromDiffuse();		float adsk_getShininess();
 bool adsk_isLightboxRenderedBeforeLight();		 vec3 adsk_getComputedDiffuse();
 vec3 adsk_getComputedSpecular();

 vec4 adsk_getDiffuseMapValue(vec2 texCoord);	 vec4 adsk_getParallaxMapValue(vec2 texCoord);
 vec4 adsk_getEmissiveMapValue(vec2 texCoord);	 vec4 adsk_getUVMapValue(vec2 texCoord);
 vec4 adsk_getSpecularMapValue(vec2 texCoord);   vec4 adsk_getReflectanceMapValue(vec2 texCoord);
 vec4 adsk_getNormalMapValue(vec2 texCoord);
 
  int adsk_getNumberIBLs();						 vec3 adsk_getCubeMapIBL(int idx, vec3 coords, float lod);
 bool adsk_isCubeMapIBL(int idx);				 vec3 adsk_getAngularMapIBL(int idx, vec2 coords, float lod);
 bool adsk_isAmbientIBL(int idx);				float adsk_getIBLDiffuseOffset(int idx);

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
	vec3 p = adsk_getVertexPosition();
	vec3 n = adsk_getNormal();
	vec3 view = normalize(cam - p);
	i.rgb = vec3(1.0 - dot(view, n));
	return i;
}
