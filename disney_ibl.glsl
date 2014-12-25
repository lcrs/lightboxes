#line 1
#extension GL_EXT_gpu_shader4 : enable
#define uint unsigned int

// TODO:
//  o rotation control for anisotropic spec
//  o use uv as tangent space basis when flame bug is fixed

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
uniform int adskUID_minsamples;
uniform int adskUID_maxsamples;
uniform float adskUID_variancelimit;
uniform bool adskUID_outputvariance;
uniform bool adskUID_outputsamples;
uniform int adskUID_method;
uniform float adskUID_lod;
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

uint adskUID_hash(uint x, uint y) {
    const uint M = 1664525u, C = 1013904223u;
    uint seed = (x * M + y + C) * M;
    seed ^= (seed >> 11u);
    seed ^= (seed << 7u) & 0x9d2c5680u;
    seed ^= (seed << 15u) & 0xefc60000u;
    seed ^= (seed >> 18u);
    return seed;
}

float adskUID_radicalInverse_VdC(uint bits, uint seed) {
    bits = (bits << 16u) | (bits >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    bits ^= seed;
    return float(bits) * 2.3283064365386963e-10; // divide by 0x100000000
 }

vec2 adskUID_hammersley(uint i, uint n, uint seed) {
    return fract(vec2(float(i)/float(n), adskUID_radicalInverse_VdC(i, seed)));
}

float adskUID_rand(vec2 co) {
    return fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec2 adskUID_halton(uint i) {
    uvec2 base = uvec2(7u, 3u);
    vec2 result = vec2(0.0);
    vec2 f = 1.0 / vec2(base);
    uvec2 index = uvec2(i, i);
    while (index.x > 0u || index.y > 0u) {
        result = result + f * vec2(index % base);
        index = index / base;
        f = f / vec2(base);
    }
    return result;
}

float adskUID_luma(vec3 c) {
    return dot(c, vec3(0.2126, 0.7152, 0.0722));
}

vec3 adskUID_hemi(vec2 uv) {
    float phi = uv.y * 2.0 * adskUID_PI;
    float costheta = 1.0 - uv.x;
    float sintheta = sqrt(1.0 - costheta * costheta);
    return normalize(vec3(cos(phi) * sintheta, sin(phi) * sintheta, costheta));
}

vec2 adskUID_cylinder(vec3 v) {   
   const float factor = 1.0 / (2.0 * adskUID_PI);
   const float offset = -0.25;
   float x = atan(v.z, v.x ) * factor + offset;
   return vec2(fract(x), v.y * 0.5 + 0.5);
}

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
    vec3 p = adsk_getVertexPosition();
    vec3 v = normalize(cam - p);
    vec3 n = adsk_getNormal();
    vec3 t = cross(n, vec3(0.0, 1.0, 0.0));
    vec3 b = cross(n, t);
    mat3 world2tangent = mat3(t, b, n);

    // Per-fragment seed to decorrelate sampling pattern
    uint seed = adskUID_hash(adskUID_hash(uint(gl_FragCoord.x*abs(p.x*1234.5)), uint(gl_FragCoord.y*abs(p.y*1234.5))), uint(adsk_getTime()));

    // Accumulator
    vec3 a = vec3(0.0);

    // Max samples we will ever take
    uint sn = uint(adskUID_maxsamples);

    // How many samples we actually take to get a good low variance
    uint bail = 0u;

    // State variables for continuous variance calculation
    float m = 0.0, s = 0.0, variance = 0.0;

    for(uint i = 0u; i < sn; i++) {
        // Generate position for this sample
        vec2 sample = vec2(0.0, 0.0);
        if(adskUID_method == 0) {
            sample = adskUID_hammersley(i, sn, seed);
        } else if(adskUID_method == 1) {
            sample = adskUID_hammersley(i, sn, seed);
            sample.x = fract(float(adskUID_hash(seed, i))/123456.7); // otherwise x increments regularly
        } else if(adskUID_method == 2) {
            sample = fract(vec2(adskUID_hash(i, seed), adskUID_hash(i+1u, seed)) * 2.3283064365386963e-10);
        } else if(adskUID_method == 3) {
            sample = vec2(adskUID_rand(float(i+10u) * p.xy * 0.01 * adsk_getTime()), adskUID_rand(float(i+20u) * p.yx * 0.01 * adsk_getTime()));
        } else if(adskUID_method == 4) {
            sample = adskUID_halton(i + seed/24u);
        }

        // Go!
        vec3 l = world2tangent * adskUID_hemi(sample);
        vec3 env = adsk_getAngularMapIBL(0, adskUID_cylinder(l), adskUID_lod);
        vec3 b = env * adskUID_BRDF(l, v, n, t, b) * dot(n, l);

        // Accumulate
        a += b;

        // Update variance
        float mk1 = m;
        float sk1 = s;
        float k = float(i);
        float xk = adskUID_luma(pow(b, vec3(1.0/2.4))); // Perceptual weighting, luma of gamma'd colour
        if(i == 0u) m = xk; else m = mk1 + (xk - mk1)/k;
        s = sk1 + (xk - mk1) * (xk - m);
        variance = s/(k-1.0);
        bail++;
        if((adskUID_method != 0) && (bail >= uint(adskUID_minsamples)) && (variance < (adskUID_variancelimit/100.0))) break;
    }

    i.rgb = a / float(bail);
    i.a = 1.0;

    if(adskUID_outputsamples) {
        if(bail == uint(adskUID_minsamples)) i.r = 1.0; else i.r = 0.0;
        if(bail == uint(adskUID_maxsamples)) i.b = 1.0; else i.b = 0.0;
        i.g = float(bail)/100.0;
    }
    if(adskUID_outputvariance) i.rgb = vec3(variance);

	return i;
}