// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Re-jigged for Matchbox by Ivar
// Further re-jigged for Lightbox by Lewis
// TODO:
//  o trace area around origin instead of around camera?
//  o do something about slowness, stripy sampling
//  o combining with lighting is a bit weird, gets too hot through clouds


vec3 adsk_getVertexPosition();
vec3 adsk_getCameraPosition();
vec3 adsk_getLightPosition();
vec3 adsk_getLightDirection();
vec3 adsk_getLightTangent();
mat4 adsk_getModelViewMatrix();
mat4 adsk_getModelViewInverseMatrix();
vec3 adsk_getComputedDiffuse();
float adsk_getTime();

uniform vec2 adskUID_paramPos;
uniform float adskUID_paramSpeed;
uniform bool adskUID_paramProcedural;
uniform int adskUID_steps;
uniform float adskUID_far, adskUID_scale, adskUID_fuzz;

uniform float adskUID_Detail;
uniform float adskUID_Density;
uniform float adskUID_Volume;

// hash based 3d value noise
float adskUID_hash( float n )
{
    return fract(sin(n)*43758.5453);
}

float adskUID_noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*((adskUID_Detail+1.0)-adskUID_Detail*f);
    float n = p.x + p.y*57.0 + 113.0*p.z;
    return mix(mix(mix( adskUID_hash(n+  0.0), adskUID_hash(n+  1.0),f.x),
           mix( adskUID_hash(n+ 57.0), adskUID_hash(n+ 58.0),f.x),f.y),
           mix(mix( adskUID_hash(n+113.0), adskUID_hash(n+114.0),f.x),
           mix( adskUID_hash(n+170.0), adskUID_hash(n+171.0),f.x),f.y),f.z);
}


vec4 adskUID_map( in vec3 p )
{
	float d = adskUID_Volume - p.y;

	vec3 q = p - vec3(1.0,0.1,0.0)*(adsk_getTime()*adskUID_paramSpeed/1000.0);
	float f;
    f  = 0.5000*adskUID_noise( q ); q = q*2.02;
    f += 0.2500*adskUID_noise( q ); q = q*2.03;
    f += 0.1250*adskUID_noise( q ); q = q*2.01;
    f += 0.0625*adskUID_noise( q );

	d += adskUID_Density * f;

	d = clamp( d, 0.0, 1.0 );
	
	vec4 res = vec4( d );

	res.xyz = mix( 1.15*vec3(1.0,0.95,0.8), vec3(0.7,0.7,0.7), res.x );
	
	return res;
}


vec3 adskUID_sundir = vec3(-1.0,0.0,0.0);


vec4 adskUID_raymarch( in vec3 ro, in vec3 rd, float z )
{
	vec4 sum = vec4(0, 0, 0, 0);

	float t = 0.1;
    float stepsize = (adskUID_far/adskUID_scale) / float(adskUID_steps);
	for(int i=0; i<adskUID_steps; i++)
	{
		if( sum.a > 0.99 ) continue;

		vec3 pos = ro + t*rd;

		if(t > z) {
			// We've gone beyond our geo depth! Comp the sum so far over the geo
			vec3 bg = adsk_getComputedDiffuse();
			vec3 comp = bg * (1.0 - sum.a) + sum.rgb;

			return vec4(comp, 1.0);
		}

		vec4 col = adskUID_map( pos );
		
		float dif =  clamp((col.w - adskUID_map(pos+0.3*adskUID_sundir).w)/0.6, 0.0, 1.0 );
        vec3 lin = vec3(0.65,0.68,0.7)*1.35 + 0.45*vec3(0.7, 0.5, 0.3)*dif;
		col.xyz *= lin;
		
		col.a *= 0.35;
		col.rgb *= col.a;

		sum = sum + col*(1.0 - sum.a);	

        //t += max(0.1,0.025*t);
        t += stepsize;
	}

	sum.xyz /= (0.001+sum.w);

	return clamp( sum, 0.0, 1.0 );
}

float adskUID_rand(vec2 co) {
    return fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec4 adskUID_lightbox(vec4 i)
{
    vec3 camera = adsk_getCameraPosition();
    vec3 vertex = adsk_getVertexPosition();
    vec3 light = -adsk_getLightPosition();
    vec3 lightdir = adsk_getLightDirection();
    vec3 lighttan = adsk_getLightTangent();
    vec3 lightbitan = cross(lightdir, lighttan);
    mat3 lightbasis = mat3(lighttan, -lightbitan, lightdir);
    vec3 dir = normalize(vertex - camera);
	float z = length(vertex-camera);

    vec3 ro = camera - light;
    ro *= lightbasis;
	vec3 rd = dir * lightbasis;

	ro /= adskUID_scale;
	z /= adskUID_scale;

	float fuzz = adskUID_rand(vertex.xy/1000.0);
	float fuzz2 = adskUID_rand(vec2(fuzz, adsk_getTime()));
	float fuzz3 = adskUID_rand(vec2(fuzz, fuzz2));
	ro += (adskUID_fuzz/float(adskUID_steps)) * (vec3(fuzz, fuzz2, fuzz3) - vec3(0.5));

    vec4 res = adskUID_raymarch(ro, rd, z);

	float sun = clamp( dot(adskUID_sundir,rd), 0.0, 1.0 );
	vec3 col = vec3(0.6,0.71,0.75) - rd.y*0.2*vec3(1.0,0.5,1.0) + 0.15*0.5;
	col += 0.2*vec3(1.0,.6,0.1)*pow( sun, 8.0 );
	col *= 0.95;
	col = mix( col, res.xyz, res.w );
	col += 0.1*vec3(1.0,0.4,0.2)*pow( sun, 3.0 );

    return vec4( col, 1.0 );
}
