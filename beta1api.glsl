//------------------------------------------------------------------------------
// NOTE: all positions and directions are in camera space.  that is to say their
// local positions are multiplied by the modelView matrix.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// General
//------------------------------------------------------------------------------
bool adsk_isSceneLinear();
//------------------------------------------------------------------------------
// Lighting/Shading
//------------------------------------------------------------------------------
// returns unaltered interpolated normal ( before maps )
vec3 adsk_getNormal();
// returns computed normal once (normal,parallax,displacement) maps are applied
vec3 adsk_getComputedNormal();
// returns unaltered interpolated binormal and tangent
vec3 adsk_getBinormal();
vec3 adsk_getTangent();
// returns interpolated vertex position for current fragment. Displacement and
// position maps have been applied.
vec3 adsk_getVertexPosition();
// returns the eye/camera position
vec3 adsk_getCameraPosition();
// returns current frame/time number
float adsk_getTime();
// returns activity status of parent light.
bool adsk_isLightActive();
// returns true for additive, or false for single light
bool adsk_isLightAdditive();
// returns light position vector
vec3 adsk_getLightPosition();
// returns light colour modulated by the intensity // lightColour * intensity
vec3 adsk_getLightColour();
// returns light direction for directional lights. the direction is computed
// whether or not the light is in target mode or free mode.
vec3 adsk_getLightDirection();
vec3 adsk_getLightTangent();
// Light and shadow decay types values:
// const int NO_DECAY = 0;
// const int LINEAR_DECAY = 1;
// const int QUADRATIC_DECAY = 2;
// const int CUBIC_DECAY = 3;
// const int EXP_DECAY = 4;
// const int EXP2_DECAY = 5;
void adsk_getLightShadowDecayType( out int lightDecayType,
                                   out int shadowDecayType );
float adsk_getLightDecayRate();
bool adsk_isSpotlightFalloffParametric();
float adsk_getSpotlightParametricFalloffIn();
float adsk_getSpotlightParametricFalloffOut();
float adsk_getSpotlightSpread();
// returns pre-computed alpha calculated from spotlight cutoff, decay and GMask
float adsk_getLightAlpha();
bool adsk_isPointSpotLight();
bool adsk_isDirectionalLight();
bool adsk_isAmbientLight();
bool adsk_isAreaRectangleLight();
bool adsk_isAreaEllipseLight();
float adsk_getAreaLightWidth();
float adsk_getAreaLightHeight();
//------------------------------------------------------------------------------
// Lightbox
//------------------------------------------------------------------------------
// false its from current lit fragment
bool adsk_isLightboxRenderedFromDiffuse();
// false its post light
bool adsk_isLightboxRenderedBeforeLight();
// gets computed diffuse value for current fragment before shading once all maps have
// been processed.
vec3 adsk_getComputedDiffuse();
float adsk_getShininess();
// gets computed specular value for current fragment before shading once all maps have
// been processed.
vec3 adsk_getComputedSpecular();
//------------------------------------------------------------------------------
// MAP access
//------------------------------------------------------------------------------
vec4 adsk_getDiffuseMapValue( in vec2 texCoord );
// returns the compute diffuse map value once all maps have been treated.
vec4 adsk_getComputedDiffuseMapValue( in vec3 vertexPos );
vec4 adsk_getEmissiveMapValue( in vec2 texCoord );
vec4 adsk_getSpecularMapValue( in vec2 texCoord );
vec4 adsk_getNormalMapValue( in vec2 texCoord );
vec4 adsk_getReflectionMapValue( in vec2 texCoord );;
vec4 adsk_getUVMapValue( in vec2 texCoord );
vec4 adsk_getParallaxMapValue( in vec2 texCoord );
// homogeneous coordinates need to be safe divided by z
vec3 adsk_getDiffuseMapCoord();
vec3 adsk_getEmissiveMapCoord();
vec3 adsk_getSpecularMapCoord();
vec3 adsk_getNormalMapCoord();
vec2 adsk_getReflectionMapCoord( in vec3 vrtPos, in vec3 normal );
vec3 adsk_getParallaxMapCoord();
//------------------------------------------------------------------------------
// Transforms 
//------------------------------------------------------------------------------
mat4 adsk_getModelViewMatrix();
mat4 adsk_getModelViewInverseMatrix();
//------------------------------------------------------------------------------
// IBL
//------------------------------------------------------------------------------
// depending on platform IBL supports 0, 1 or 2 maps
int adsk_getNumberIBLs();
// else its angular map
bool adsk_isCubeMapIBL( in int idx );;
// lod is the level of details in the mipmap
vec3 adsk_getCubeMapIBL( in int idx, in vec3 coords, float lod );
vec3 adsk_getAngularMapIBL( in int idx, in vec2 coords, float lod );
bool adsk_isAmbientIBL( in int idx );
float adsk_getIBLDiffuseOffset( in int idx );
mat4 adsk_getIBLRotationMatrix( in int idx );
//------------------------------------------------------------------------------
// Material Information
//------------------------------------------------------------------------------ 
ec4 adsk_getMaterialDiffuse();
vec4 adsk_getMaterialSpecular();
vec4 adsk_getMaterialAmbient();
//------------------------------------------------------------------------------
// Colour Management
//------------------------------------------------------------------------------ 
float adskEvalDynCurves(in int dynCurveId, in float x);
vec2 adskEvalDynCurves(in ivec2 dynCurveId, in vec2 x);
vec3 adskEvalDynCurves(in ivec3 dynCurveId, in vec3 x);
vec4 adskEvalDynCurves(in ivec4 dynCurveId, in vec4 x);
// Conversions between scene linear and Cineon log colour space
vec3 adsk_scene2log( in vec3 src );
vec3 adsk_log2scene( in vec3 src );
// Conversions between RGB, HSV, YUV
vec3 adsk_rgb2hsv( in vec3 c );
vec3 adsk_hsv2rgb( in vec3 c );
vec3 adsk_rgb2yuv( in vec3 c );
vec3 adsk_yuv2rgb( in vec3 c );
vec3 adsk_getLuminanceWeights();
float adsk_getLuminance( in vec3 color );
// Parametric blending curves for blending effects for Shadows, Midtones, Highlights.
// The first argument is the pixel colour. Second argument is the X-axis location,
// where the curve is 50%. To get the midtone weight, you need to compute
// 1 - adsk_shadows - adsk_highlights, for the same pixel colour
float adsk_highlights( in float pixel, in float halfPoint );
float adsk_shadows( in float pixel, in float halfPoint );