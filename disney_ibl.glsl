#line 1
#extension GL_EXT_gpu_shader4 : enable
#define uint unsigned int

// TODO:
//  o rotation control for anisotropic spec
//  o use uv as tangent space basis when flame bug is fixed

uniform vec3 adskUID_baseColor;
uniform float adskUID_metallic;
uniform float adskUID_roughness;
uniform float adskUID_subsurface;
uniform float adskUID_specular;
uniform float adskUID_specularTint;
uniform float adskUID_anisotropic;
uniform float adskUID_sheen;
uniform float adskUID_sheenTint;
uniform float adskUID_clearcoat;
uniform float adskUID_clearcoatGloss;
uniform int adskUID_samples;
uniform bool adskUID_importance;
uniform int adskUID_method;
uniform float adskUID_lod;
const float adskUID_PI = 3.14159265358979323846;

float adskUID_luma(vec3 c) { 
    return dot(c, vec3(0.2126, 0.7152, 0.072));
}

float adskUID_Cdlum = adskUID_luma(adskUID_baseColor);
vec3 adskUID_Ctint = adskUID_baseColor/adskUID_Cdlum;
vec3 adskUID_CspecTint = mix(vec3(1.0), adskUID_Ctint, adskUID_specularTint);
vec3 adskUID_CspecFZ = mix(adskUID_specular*0.08*adskUID_CspecTint, adskUID_baseColor, adskUID_metallic);
vec3 adskUID_Csheen = mix(vec3(1.0), adskUID_Ctint, adskUID_sheenTint);
vec3 adskUID_diff_sheen = adskUID_Csheen * adskUID_sheen;
float adskUID_aspect = sqrt(1.0-adskUID_anisotropic*0.9);
float adskUID_roughSqr = adskUID_roughness * adskUID_roughness;
float adskUID_alphax = max(0.001, adskUID_roughSqr/adskUID_aspect);
float adskUID_alphay = max(0.001, adskUID_roughSqr*adskUID_aspect);
float adskUID_alphaG = max(0.001, adskUID_roughSqr);
float adskUID_coatAlpha = mix(0.1, 0.001, adskUID_clearcoatGloss);
float adskUID_coatAlphaG = 0.25;
vec3 adskUID_coatFZ = vec3(0.04);

float adskUID_sqr(float v) {
    return v*v;
}

float adskUID_schlick_f(float u) {
    float m = clamp(1.0-u, 0.0, 1.0);
    float m2 = m*m;
    return m2*m2*m; 
}

float adskUID_smith_g(float Ndotv, float alphaG) {
    if (Ndotv <= 0.0) {
        return 0.0;
    }
    float a = alphaG*alphaG;
    float b = Ndotv*Ndotv;
    return 1.0/(Ndotv + sqrt(a + b - a*b));
}

vec3 adskUID_diffuse(vec3 u, vec3 v, vec3 n) {
    vec3 nn = normalize(n);
    vec3 un = normalize(v);
    vec3 vn = normalize(n);
    vec3 h = normalize(un+vn);
    float udoth = dot(un,h);
    float Ndotu = dot(un,nn);
    float a_udoth = abs(udoth);
    float a_Ndotu = abs(Ndotu);
    float Ndotv = dot(vn,nn);
    if(Ndotv > 0.0) {
        float tmp = udoth*udoth*adskUID_roughness;

        float Fd90 = 0.5 + 2.0*tmp;
        float FV = adskUID_schlick_f(a_Ndotu);
        float FL = adskUID_schlick_f(Ndotv);
        float Fd = mix(1.0, Fd90, FL) * mix(1.0, Fd90, FV);

        float Fss90 = tmp;
        float Fss = mix(1.0, Fss90, FL) * mix(1.0, Fss90, FV);
        float ss = 1.25 * (Fss * (1.0/(a_Ndotu+Ndotv)-0.5) + 0.5);

        float FH = adskUID_schlick_f(a_udoth);
        vec3 Fsheen = FH * adskUID_diff_sheen;

        vec3 eval = mix(Fd, ss, adskUID_subsurface) * adskUID_baseColor + Fsheen;
        eval *= Ndotv;
        eval *= 2.0;
        eval *= 1.0 - adskUID_metallic;
        return eval;
    } else {
        return vec3(0.0);
    } 
}

vec3 adskUID_spec(vec3 u, vec3 v, vec3 n, vec3 t, vec3 b) {
    float alphaxSqr = adskUID_alphax*adskUID_alphax;
    float alphaySqr = adskUID_alphay*adskUID_alphay;
    float rho = 0.5*(adskUID_alphax*adskUID_alphay);
    vec3 refl = vec3(rho);
    vec3 nn = normalize(n);
    vec3 xn = normalize(t);
    vec3 yn = normalize(b);
    vec3 un = normalize(v);
    vec3 vn = normalize(u);
    vec3 h = normalize(un+vn);
    float Ndotu = dot(nn,un);
    float Ndotv = dot(nn,vn);
    float udoth = dot(un,h);
    float a_udoth = abs(udoth);
    float a_Ndotu = abs(Ndotu);
    vec3 F0 = adskUID_CspecFZ;

    float cosTheta = dot(h,nn);

    if(Ndotv <= 0.0 || cosTheta <= 0.0) {
        return vec3(0.0);
    } else {
        vec3 F = vec3(1.0);
        float G = 1.0;
        float D = 1.0;
        float tmp = adskUID_sqr(dot(h,xn))/alphaxSqr +
                    adskUID_sqr(dot(h,yn))/alphaySqr +
                    adskUID_sqr(cosTheta);
        D = 1.0 / (tmp*tmp);
        F = F0 + (1.0-F0) * adskUID_schlick_f(a_udoth);
        G = adskUID_smith_g(a_Ndotu, adskUID_alphaG) * adskUID_smith_g(Ndotv, adskUID_alphaG); 
        vec3 eval = F * D * G * Ndotv;
        eval *= 1.0 / adskUID_luma(refl);
        return eval;
    }
}

bool adskUID_sameHemisphere(vec3 w, vec3 wp) {
    return (w.z * wp.z) > 0.0;
}

vec4 adskUID_spec_importance(vec3 v, vec3 n, vec3 t, vec3 b, vec2 s) {
    float alphaxSqr = adskUID_alphax*adskUID_alphax;
    float alphaySqr = adskUID_alphay*adskUID_alphay;
    float rho = 0.5*(adskUID_alphax*adskUID_alphay);
    vec3 refl = vec3(rho);
    vec3 nn = normalize(n);
    vec3 xn = normalize(t);
    vec3 yn = normalize(b);
    mat3 from_rspace = mat3(xn, yn, nn);
    xn = normalize(from_rspace * xn);
    yn = normalize(from_rspace * yn);
    vec3 F0 = adskUID_CspecFZ;
    float pdf = 0.0;

    vec3 wo = normalize(from_rspace * v);

    vec3 wh = sqrt(s.y/(1.0-s.y)) * (adskUID_alphax*cos(6.28318530717958647692 *s.x)*xn +
                                     adskUID_alphax*sin(6.28318530717958647692 *s.x)*yn) + 
                                     vec3(0.0, 0.0, 1.0);
    wh = normalize(wh);

    if(!adskUID_sameHemisphere(wo,wh)) {
            wh = -wh;
    }

    float cosTheta = wh.z;

    if ( wh.z > 0.0 ) {
        vec3 wi = 2.0*dot(wo,wh)*wh - wo;

        if (wi.z <= 0.0 ) {
            pdf = 0.0;
            refl = vec3(0.0);
        } else {
            float Ndotu = wo.z;
            float Ndotv = wi.z;
            float udoth = dot(wo,wh);
            float a_Ndotu = abs(Ndotu);
            float a_udoth = abs(udoth);
            float cosThetaSqr = cosTheta*cosTheta;

            float pdf_h_to_wi = 1.0 / (4.0*a_udoth);

            float D = 1.0;
            float tmp = adskUID_sqr(dot(wh,xn))/alphaxSqr +
                        adskUID_sqr(dot(wh,yn))/alphaySqr +
                        adskUID_sqr(cosTheta);
            D = 1.0 / (tmp*tmp);

            pdf = D * abs(wh.z);
            pdf /= rho;
            pdf *= pdf_h_to_wi;
            vec3 F = F0 + (1.0-F0) * adskUID_schlick_f(a_udoth);
            float G = adskUID_smith_g(a_Ndotu, adskUID_alphaG) * adskUID_smith_g(Ndotv, adskUID_alphaG);
            refl *=  F * G * Ndotv;
            refl /= pdf_h_to_wi * abs(wh.z);
        }
        vec3 dir = normalize(wi * from_rspace);
        pdf = 1.0/max(0.1, pdf);
        return vec4(dir, pdf);
    }
}

vec3 adskUID_coat(vec3 u, vec3 v, vec3 ng) {
    float alphaSqr = adskUID_coatAlpha*adskUID_coatAlpha;
    float rho = log(adskUID_coatAlpha)/(-1.0+alphaSqr);
    vec3 refl = vec3(rho);
    vec3 n = normalize(ng);
    vec3 un = normalize(v);
    vec3 vn = normalize(u);
    vec3 h = normalize(un+vn);
    float Ndotu = dot(n,un);
    float Ndotv = dot(n,vn);
    float udoth = dot(un,h);
    float a_udoth = abs(udoth);
    float a_Ndotu = abs(Ndotu);

    float cosTheta = dot(h,n);

    if(Ndotv <= 0.0 || cosTheta <= 0.0) {
        return vec3(0.0);
    } else {
        vec3 F = vec3(1.0);
        float G = 1.0;
        float D = 1.0;
        float alphaSqrM1 = alphaSqr - 1.0;
        float cosThetaSqr = cosTheta*cosTheta;

        D = 1.0/(cosThetaSqr*alphaSqrM1+1.0);

        F = adskUID_coatFZ + (1.0-adskUID_coatFZ) * adskUID_schlick_f(a_udoth);
        G = adskUID_smith_g(a_Ndotu, adskUID_coatAlphaG) * adskUID_smith_g(Ndotv, adskUID_coatAlphaG); 
        vec3 eval = F * D * G * Ndotv;
        eval *= 0.25 / adskUID_luma(refl);
        return eval;
    }
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

    // Samples we will take
    uint sn = uint(adskUID_samples);

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
        vec4 importance;
        vec3 u;
        if(adskUID_importance) {
            importance = adskUID_spec_importance(v, n, t, b, sample);
            u = importance.xyz;
        } else {
            importance = vec4(0.0, 0.0, 0.0, 1.0);
            u = world2tangent * adskUID_hemi(sample);
        }
        vec3 env = adsk_getAngularMapIBL(0, adskUID_cylinder(u), adskUID_lod);
        vec3 light = vec3(0.0);
        //light += env * adskUID_diffuse(u, v, n);
        light += env * adskUID_spec(u, v, n, t, b);
        //light += env * adskUID_coat(u, v, n);
        light = light * importance.a;
        //light = min(vec3(2.0), light);
        a += light;
    }
    a /= float(sn);

    i.rgb += adsk_getComputedDiffuse() * i.a * adsk_getLightColour() * a;
	return i;
}
