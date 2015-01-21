// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

vec3 adsk_getVertexPosition();
vec3 adsk_getCameraPosition();
vec3 adsk_getLightPosition();
mat4 adsk_getModelViewMatrix();
mat4 adsk_getModelViewInverseMatrix();
float adsk_getTime();

uniform vec2 adskUID_paramPos;
uniform float adskUID_paramSpeed;
uniform bool adskUID_paramProcedural;
uniform int adskUID_steps;

uniform float adskUID_Detail;
uniform float adskUID_Density;
uniform float adskUID_Volume;

// hash based 3d value noise
float hash( float n )
{
    return fract(sin(n)*43758.5453);
}

float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*((adskUID_Detail+1.0)-adskUID_Detail*f);
    float n = p.x + p.y*57.0 + 113.0*p.z;
    return mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
           mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
           mix(mix( hash(n+113.0), hash(n+114.0),f.x),
           mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
}


vec4 map( in vec3 p )
{
	float d = adskUID_Volume - p.y;

	vec3 q = p - vec3(1.0,0.1,0.0)*(adsk_getTime()*adskUID_paramSpeed/1000.0);
	float f;
    f  = 0.5000*noise( q ); q = q*2.02;
    f += 0.2500*noise( q ); q = q*2.03;
    f += 0.1250*noise( q ); q = q*2.01;
    f += 0.0625*noise( q );

	d += adskUID_Density * f;

	d = clamp( d, 0.0, 1.0 );
	
	vec4 res = vec4( d );

	res.xyz = mix( 1.15*vec3(1.0,0.95,0.8), vec3(0.7,0.7,0.7), res.x );
	
	return res;
}


vec3 sundir = vec3(-1.0,0.0,0.0);


vec4 raymarch( in vec3 ro, in vec3 rd )
{
	vec4 sum = vec4(0, 0, 0, 0);

	float t = 0.0;
	for(int i=0; i<adskUID_steps; i++)
	{
		if( sum.a > 0.99 ) continue;

		vec3 pos = ro + t*rd;
		vec4 col = map( pos );
		
		#if 1
		float dif =  clamp((col.w - map(pos+0.3*sundir).w)/0.6, 0.0, 1.0 );

        vec3 lin = vec3(0.65,0.68,0.7)*1.35 + 0.45*vec3(0.7, 0.5, 0.3)*dif;
		col.xyz *= lin;
		#endif
		
		col.a *= 0.35;
		col.rgb *= col.a;

		sum = sum + col*(1.0 - sum.a);	

        #if 0
		t += 0.1;
		#else
		t += max(0.1,0.025*t);
		#endif
	}

	sum.xyz /= (0.001+sum.w);

	return clamp( sum, 0.0, 1.0 );
}

vec4 adskUID_lightbox(vec4 i)
{
    vec4 camera = adsk_getModelViewInverseMatrix() * vec4(adsk_getCameraPosition(), 1.0);
    vec4 vertex = adsk_getModelViewInverseMatrix() * vec4(adsk_getVertexPosition(), 1.0);
    vec3 dir = normalize(vertex.xyz - camera.xyz);
    //camera.y += 200.0;
    //camera.z += 600.0;
    vec3 ro = camera.xyz;
	vec3 rd = dir;

    vec4 res = raymarch( ro, rd );

	float sun = clamp( dot(sundir,rd), 0.0, 1.0 );
	vec3 col = vec3(0.6,0.71,0.75) - rd.y*0.2*vec3(1.0,0.5,1.0) + 0.15*0.5;
	col += 0.2*vec3(1.0,.6,0.1)*pow( sun, 8.0 );
	col *= 0.95;
	col = mix( col, res.xyz, res.w );
	col += 0.1*vec3(1.0,0.4,0.2)*pow( sun, 3.0 );

    return vec4( col, 1.0 );
}
