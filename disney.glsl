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

uniform vec3 adskUID_baseColor;
uniform float adskUID_metallic;
uniform float adskUID_subsurface;
uniform float adskUID_specular;
uniform float adskUID_roughness;
uniform float adskUID_specularTint;
uniform float adskUID_anisotropic;
uniform float adskUID_sheen;
uniform float adskUID_sheenTint;
uniform float adskUID_clearcoat;
uniform float adskUID_clearcoatGloss;
const float adskUID_PI = 3.14159265358979323846;

float adskUID_sqr(float x) { return x*x; }

float adskUID_SchlickFresnel(float u)
{
    float m = clamp(1.0-u, 0.0, 1.0);
    float m2 = m*m;
    return m2*m2*m; // pow(m,5)
}

float adskUID_GTR1(float NdotH, float a)
{
    if (a >= 1.0) return 1.0/adskUID_PI;
    float a2 = a*a;
    float t = 1.0 + (a2-1.0)*NdotH*NdotH;
    return (a2-1.0) / (adskUID_PI*log(a2)*t);
}

float adskUID_GTR2(float NdotH, float a)
{
    float a2 = a*a;
    float t = 1.0 + (a2-1.0)*NdotH*NdotH;
    return a2 / (adskUID_PI * t*t);
}

float adskUID_GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay)
{
    return 1.0 / ( adskUID_PI * ax*ay * adskUID_sqr( adskUID_sqr(HdotX/ax) + adskUID_sqr(HdotY/ay) + NdotH*NdotH ));
}

float adskUID_smithG_GGX(float Ndotv, float alphaG)
{
    float a = alphaG*alphaG;
    float b = Ndotv*Ndotv;
    return 1.0/(Ndotv + sqrt(a + b - a*b));
}

vec3 adskUID_BRDF( vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y )
{
    float NdotL = dot(N,L);
    float NdotV = dot(N,V);
    if (NdotL < 0.0 || NdotV < 0.0) return vec3(0.0);

    vec3 H = normalize(L+V);
    float NdotH = dot(N,H);
    float LdotH = dot(L,H);

    vec3 Cdlin = adskUID_baseColor;
    float Cdlum = 0.3*Cdlin[0] + 0.6*Cdlin[1]  + 0.1*Cdlin[2]; // luminance approx.

    vec3 Ctint = Cdlum > 0.0 ? Cdlin/Cdlum : vec3(1.0); // normalize lum. to isolate hue+sat
    vec3 Cspec0 = mix(adskUID_specular*.08*mix(vec3(1.0), Ctint, adskUID_specularTint), Cdlin, adskUID_metallic);
    vec3 Csheen = mix(vec3(1.0), Ctint, adskUID_sheenTint);

    // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing
    // and mix in diffuse retro-reflection based on roughness
    float FL = adskUID_SchlickFresnel(NdotL), FV = adskUID_SchlickFresnel(NdotV);
    float Fd90 = 0.5 + 2.0 * LdotH*LdotH * adskUID_roughness;
    float Fd = mix(1.0, Fd90, FL) * mix(1.0, Fd90, FV);

    // Based on Hanrahan-Krueger brdf approximation of isotropic bssrdf
    // 1.25 scale is used to (roughly) preserve albedo
    // Fss90 used to "flatten" retroreflection based on roughness
    float Fss90 = LdotH*LdotH*adskUID_roughness;
    float Fss = mix(1.0, Fss90, FL) * mix(1.0, Fss90, FV);
    float ss = 1.25 * (Fss * (1.0 / (NdotL + NdotV) - 0.5) + 0.5);

    // specular
    float aspect = sqrt(1.0-adskUID_anisotropic*0.9);
    float ax = max(0.001, adskUID_sqr(adskUID_roughness)/aspect);
    float ay = max(0.001, adskUID_sqr(adskUID_roughness)*aspect);
    float Ds = adskUID_GTR2_aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
    float FH = adskUID_SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1.0), FH);
    float roughg = adskUID_sqr(adskUID_roughness*0.5+0.5);
    float Gs = adskUID_smithG_GGX(NdotL, roughg) * adskUID_smithG_GGX(NdotV, roughg);

    // sheen
    vec3 Fsheen = FH * adskUID_sheen * Csheen;

    // clearcoat (ior = 1.5 -> F0 = 0.04)
    float Dr = adskUID_GTR1(NdotH, mix(0.1,0.001,adskUID_clearcoatGloss));
    float Fr = mix(0.04, 1.0, FH);
    float Gr = adskUID_smithG_GGX(NdotL, 0.25) * adskUID_smithG_GGX(NdotV, 0.25);

    return ((1.0/adskUID_PI) * mix(Fd, ss, adskUID_subsurface)*Cdlin + Fsheen)
        * (1.0-adskUID_metallic)
        + Gs*Fs*Ds + 0.25*adskUID_clearcoat*Gr*Fr*Dr;
}

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
    vec3 light = adsk_getLightPosition();
	vec3 p = adsk_getVertexPosition();
    vec3 l = normalize(light - p);
    vec3 v = normalize(cam - p);
	vec3 n = adsk_getNormal();
    vec3 t = adsk_getTangent();
    vec3 b = adsk_getBinormal();

    i.rgb *= adskUID_BRDF(l, v, n, t, b);

	return i;
}
