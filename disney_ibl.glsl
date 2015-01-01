#line 1
#extension GL_EXT_gpu_shader4 : enable
#define uint unsigned int

// TODO:
//  o rotation control for anisotropic spec
//  o use uv as tangent space basis when flame bug is fixed
//  o get camera rotation via expression link for now? does env rotation matrix do this?
//  o adaptive degradation levels: on manip, degrade, normal, preview
//  o adaptive mipmap selection "filtered importantace sampling" from gpu gems 3?
//  o at least adapt mipmap LOD when roughness is turned up?
//  o cube maps vs latlongs... mipmaps per face defeat FIS... what's distortion metric for latlongs?
//  o can we use another texture for environment importance? how do we combine w/brdf importance?
//  o optimize, shaderanalyzer
//  o max color limit for firefly problem... when, where, how much
//  o spec/coat multipliers don't seem to do much?

uniform vec3 adskUID_baseColor;
uniform float adskUID_metallic;
uniform float adskUID_roughness;
uniform float adskUID_subsurface;
uniform float adskUID_specular;
uniform float adskUID_specularTint;
uniform float adskUID_sheen;
uniform float adskUID_sheenTint;
uniform float adskUID_clearcoat;
uniform float adskUID_clearcoatGloss;
uniform float adskUID_anisotropic;
uniform bool adskUID_diffuseon, adskUID_specon, adskUID_coaton;
uniform int adskUID_samples;
uniform bool adskUID_diffuseimportance, adskUID_specimportance, adskUID_coatimportance;
uniform float adskUID_diffuselod, adskUID_speclod, adskUID_coatlod;
uniform int adskUID_method;
uniform vec3 adskUID_camrot;
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

bool adskUID_sameHemisphere(vec3 a, vec3 b) {
    return (a.z * b.z) > 0.0;
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

vec3 adskUID_diffuse(vec3 u, vec3 v, vec3 n, vec3 t, vec3 b) {
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
        // eval *= 2.0;  // this should be active in theory but doesn't match houdini...
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
    vec3 h = normalize(u+v);
    float Ndotu = dot(n,v);
    float Ndotv = dot(n,u);
    float udoth = dot(v,h);
    float a_udoth = abs(udoth);
    float a_Ndotu = abs(Ndotu);
    float cosTheta = dot(h,n);

    if(Ndotv <= 0.0 || cosTheta <= 0.0) {
        return vec3(0.0);
    } else {
        float tmp = adskUID_sqr(dot(h,t))/alphaxSqr + adskUID_sqr(dot(h,b))/alphaySqr + adskUID_sqr(cosTheta);
        float D = 1.0 / (tmp*tmp);
        vec3 F = adskUID_CspecFZ + (1.0-adskUID_CspecFZ) * adskUID_schlick_f(a_udoth);
        float G = adskUID_smith_g(a_Ndotu, adskUID_alphaG) * adskUID_smith_g(Ndotv, adskUID_alphaG); 
        vec3 eval = F * D * G * Ndotv;
        eval *= 1.0 / adskUID_luma(refl);
        return eval;
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

void adskUID_makebasis(inout vec3 x, inout vec3 y, vec3 z) {
        if (abs(z.x) < 0.6)
                y = vec3(1.0, 0.0, 0.0);
        else if (abs(z.z) < 0.6)
                y = vec3(0.0, 1.0, 0.0);
        else
                y = vec3(0.0, 0.0, 1.0);
        x = normalize(cross(y,z));
        y = cross(z,x);
}

void adskUID_makebasis(inout vec3 x, inout vec3 y, vec3 z, vec3 u) {
    x = normalize(cross(z,u));
    if (length(x) == 0.0)
        adskUID_makebasis(x,y,z);
    else
        y = cross(z,x);
}

vec4 adskUID_diffuse_importance(vec3 v, vec3 n, vec3 t, vec3 b, vec2 s) {
    vec3 refl;
    vec3 dir = vec3(cos(s.x*6.28318530717958647692), sin(s.x*6.28318530717958647692), 0);
    dir *= sqrt(s.y);
    dir.z = sqrt(1.0-s.y);

    float pdf = 2.0*dir.z;

    vec3 framex = vec3(0.0);
    vec3 framey = vec3(0.0);
    adskUID_makebasis(framex, framey, n, v);

    dir = dir.x*framex + dir.y*framey + dir.z*n;

    vec3 un = normalize(v);
    vec3 vn = normalize(dir);

    if (dot(n,vn) > 0.0) {

        vec3 h = normalize(v+dir);
        float udoth = dot(un,h);
        float Ndotu = dot(un,n);
        float Ndotv = dot(vn,n);
        float a_udoth = abs(udoth);
        float a_Ndotu = abs(Ndotu);
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

        refl = vec3(0.5);
        refl *= mix(Fd, ss, adskUID_subsurface) * adskUID_baseColor + Fsheen;
    } else {
        pdf = 0.0;
        refl = vec3(0.0);
    }
    pdf = 1.0/max(0.01, pdf);
    return vec4(dir, pdf);
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
    xn = normalize(xn*from_rspace);
    yn = normalize(yn*from_rspace);
    vec3 F0 = adskUID_CspecFZ;
    float pdf = 0.0;

    vec3 wo = normalize(v*from_rspace);

    vec3 wh = sqrt(s.y/(1.0-s.y)) * (adskUID_alphax*cos(6.28318530717958647692 *s.x)*xn +
                                     adskUID_alphay*sin(6.28318530717958647692 *s.x)*yn) + 
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
        vec3 dir = normalize(from_rspace*wi);
        pdf = 1.0/max(0.01, pdf);
        return vec4(dir, pdf);
    }
}

vec3 adskUID_sphericalDirection(float sinTheta, float cosTheta, float phi) {
    return vec3(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);
}

vec4 adskUID_coat_importance(vec3 v, vec3 n, vec3 t, vec3 b, vec2 s) {
    float alphaSqr = adskUID_coatAlpha*adskUID_coatAlpha;
    float rho = log(adskUID_coatAlpha)/(-1.0+alphaSqr);
    vec3 refl = vec3(rho);
    vec3 framex = vec3(0.0);
    vec3 framey = vec3(0.0);
    adskUID_makebasis(framex, framey, n, v);
    mat3 from_rspace = mat3(framex, framey, n);
    vec3 F0 = adskUID_coatFZ;
    float pdf = 0.0;

    vec3 wo = normalize(v*from_rspace);

    float cosThetaSqr = 1.0/(1.0-alphaSqr);
    cosThetaSqr *= 1.0-pow(alphaSqr, s.x);
    float cosTheta = sqrt(cosThetaSqr);
    float sinTheta = sqrt(max(0.0,1.0-cosThetaSqr));
    float phi = s.y * 6.28318530717958647692;
    vec3 wh = adskUID_sphericalDirection(sinTheta, cosTheta, phi);

    if(!adskUID_sameHemisphere(wo,wh)) {
            wh = -wh;
    }

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
            float pdf_h_to_wi = 1.0 / (4.0*a_udoth);

            float D = 1.0;
            float alphaSqrM1 = alphaSqr - 1.0;
            D = 1.0 / (cosThetaSqr*alphaSqrM1+1.0);

            pdf = D * wh.z;
            pdf /= rho;
            pdf *= pdf_h_to_wi;
            vec3 F = F0 + (1.0-F0) * adskUID_schlick_f(a_udoth);
            float G = adskUID_smith_g(a_Ndotu, adskUID_coatAlphaG) * adskUID_smith_g(Ndotv, adskUID_coatAlphaG);
            refl *= F * G * Ndotv;
            refl /= pdf_h_to_wi * wh.z;
        }
        vec3 dir = normalize(from_rspace*wi);
        pdf = 1.0/max(0.01, pdf);
        return vec4(dir, pdf);
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

// Degrees to radians
float adskUID_deg2rad(float angle) {
    return(angle/(180.0/adskUID_PI));
}

// Rotates in ZXY order
vec3 adskUID_rotate(vec3 p, vec3 angles) {
    float x = adskUID_deg2rad(angles.x);
    float y = adskUID_deg2rad(angles.y);
    float z = adskUID_deg2rad(angles.z);
    mat3 rx = mat3(1.0, 0.0, 0.0, 0.0, cos(x), sin(x), 0.0, -sin(x), cos(x));
    mat3 ry = mat3(cos(y), 0.0, -sin(y), 0.0, 1.0, 0.0, sin(y), 0.0, cos(y));
    mat3 rz = mat3(cos(z), sin(z), 0.0, -sin(z), cos(z), 0.0, 0.0, 0.0, 1.0);
    mat3 r = ry * rx * rz;
    return(p * r);
}

vec2 adskUID_latlong(vec3 v) {
    // Bodge: rotate by opposite of camera rotation
    v = adskUID_rotate(v, adskUID_camrot);
    float lat = asin(v.y) / adskUID_PI + 0.5;
    float lon = atan(v.z, v.x) / (2.0*adskUID_PI) + 0.75;
    return fract(vec2(lon, lat));
}

vec4 adskUID_lightbox(vec4 i) {
	vec3 cam = adsk_getCameraPosition();
    vec3 p = adsk_getVertexPosition();
    vec3 v = normalize(cam - p);
    vec3 n = adsk_getNormal();
    vec3 t = normalize(cross(n, vec3(1.0, -1.0, 0.0))); // constant tangent basis - geo tangents currently broken in Flame
    vec3 b = cross(n, t);
    mat3 world2tangent = mat3(t, b, n);

    // Per-fragment seed to decorrelate sampling pattern
    uint seed = adskUID_hash(adskUID_hash(uint(gl_FragCoord.x*abs(p.x*1234.5)), uint(gl_FragCoord.y*abs(p.y*1234.5))), uint(adsk_getTime()));

    vec3 accum = vec3(0.0);
    for(uint i = 0u; i < uint(adskUID_samples); i++) {
        vec2 sample;
        vec3 light;
        vec3 returned = vec3(0.0);

        // Generate random position for this sample
        if(adskUID_method == 0) {
            sample = adskUID_hammersley(i, uint(adskUID_samples), seed);
        } else if(adskUID_method == 1) {
            sample = adskUID_hammersley(i, uint(adskUID_samples), seed);
            sample.x = fract(float(adskUID_hash(seed, i))/123456.7); // otherwise x increments regularly
        } else if(adskUID_method == 2) {
            sample = fract(vec2(adskUID_hash(i, seed), adskUID_hash(i+1u, seed)) * 2.3283064365386963e-10);
        } else if(adskUID_method == 3) {
            sample = vec2(adskUID_rand(float(i+10u) * p.xy * 0.01 * adsk_getTime()), adskUID_rand(float(i+20u) * p.yx * 0.01 * adsk_getTime()));
        } else if(adskUID_method == 4) {
            sample = adskUID_halton(i + seed/24u);
        }

        // Default uniform hemisphere sampling
        vec4 importance = vec4(world2tangent * adskUID_hemi(sample), 1.0);

        // For diffuse/spec/coat: warp to important direction, sample, multiply by BRDF and importance, accumulate
        if(adskUID_diffuseon) {
            if(adskUID_diffuseimportance) importance = adskUID_diffuse_importance(v, n, t, b, sample);
            light = adsk_getAngularMapIBL(0, adskUID_latlong(importance.xyz), adskUID_diffuselod);
            returned += light * adskUID_diffuse(importance.xyz, v, n, t, b) * importance.a;
        }
        if(adskUID_specon) {
            if(adskUID_specimportance) importance = adskUID_spec_importance(v, n, t, b, sample);
            light = adsk_getAngularMapIBL(0, adskUID_latlong(importance.xyz), adskUID_speclod);
            returned += light * adskUID_spec(importance.xyz, v, n, t, b) * importance.a;
        }
        if(adskUID_coaton) {
            if(adskUID_coatimportance) importance = adskUID_coat_importance(v, n, t, b, sample);
            light = adsk_getAngularMapIBL(0, adskUID_latlong(importance.xyz), adskUID_coatlod);
            returned += light * adskUID_coat(importance.xyz, v, n) * importance.a;
        }

        accum += returned;
    }
    accum /= float(uint(adskUID_samples));

    i.rgb += adsk_getComputedDiffuse() * i.a * adsk_getLightColour() * accum;
	return i;
}
