#define MARCH_ITER 200
#define RAYAMP_MIN 0.01
#define REFLECT_MAX 10.0
#define REFLECT_PER_PATH 1
#define INIT_LEN 0.005
#define FOG_LENGTH 0.1
#define FOG_THRESHOLD 0.0001

// -----

#define MTL_AIR 1
#define MTL_FOG 2
#define MTL_FLOOR 3
#define MTL_MANDELBULB 4
#define MTL_SKY 5
#define MTL_LIGHT 6
#define MTL_TUBE 7
#define MTL_AIR_FROM_TUBE 8
#define MTL_WATER 9
#define MTL_SHIP_SHELL 10
#define MTL_SHIP_CORE 11

// ------

#define PI 3.14159265
#define V vec2(0.,1.)
#define saturate(i) clamp(i,0.,1.)
#define lofi(i,m) (floor((i)/(m))*(m))

// ------

#extension GL_EXT_draw_buffers : require
precision highp float;

uniform float time;
uniform vec2 resolution;
uniform bool reset;

uniform sampler2D textureRandom;
uniform sampler2D textureRandomStatic;

uniform sampler2D textureDrawBuffers0;
uniform sampler2D textureDrawBuffers1;
uniform sampler2D textureDrawBuffers2;
uniform sampler2D textureDrawBuffers3;

// ------

vec4 seed;
float random() { // weird prng
  const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
  const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
  const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
  const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);

  vec4 beta = floor(seed / q);
  vec4 p = a * (seed - beta * q) - beta * r;
  beta = (sign(-p) + vec4(1.0)) * vec4(0.5) * m;
  seed = (p + beta);

  return fract(dot(seed / m, vec4(1.0, -1.0, 1.0, -1.0)));
}

vec4 random4() {
  return vec4(
    random(),
    random(),
    random(),
    random()
  );
}

// ------

mat2 rotate2D( float _t ) {
  return mat2( cos( _t ), sin( _t ), -sin( _t ), cos( _t ) );
}

bool isUvValid( vec2 _v ) {
  return ( 0.0 <= _v.x ) && ( _v.x <= 1.0 ) && ( 0.0 <= _v.y ) && ( _v.y <= 1.0 );
}

float smin( float _a, float _b, float _k ) {
  float h = clamp( 0.5 + 0.5 * ( _b - _a ) / _k, 0.0, 1.0 );
  return mix( _b, _a, h ) - _k * h * ( 1.0 - h );
}

vec4 noise( vec2 _uv ) {
  vec4 sum = V.xxxx;
  for ( int i = 0; i < 6; i ++ ) {
    float mul = pow( 2.0, float( i ) );
    sum += texture2D( textureRandomStatic, _uv / 64.0 * mul ) / 2.0 / mul;
  }
  return sum;
}

// ------

vec2 randomCircle() {
  vec2 v = V.xx;
  for ( int i = 0; i < 99; i ++ ) {
    v = random4().xy * 2.0 - 1.0;
    if ( length( v ) < 1.0 ) { break; }
  }
  return v;
}

vec3 randomSphere() {
  vec3 v = V.xxx;
  for ( int i = 0; i < 99; i ++ ) {
    v = random4().xyz * 2.0 - 1.0;
    if ( length( v ) < 1.0 ) { break; }
  }
  v = normalize( v );
  return v;
}

vec3 randomHemisphere( in vec3 _normal ) {
  vec3 v = randomSphere();
  if ( dot( v, _normal ) < 0.0 ) { v = -v; }
  return v;
}

// Ref: https://pathtracing.wordpress.com/2011/03/03/cosine-weighted-hemisphere/
vec3 randomHemisphereCosWeighted( in vec3 _normal ) {
  float theta = acos( sqrt( 1.0 - random() ) );
  float phi = 2.0 * PI * random();

  vec3 sid = normalize( cross( V.xyx, _normal ) );
  vec3 top = normalize( cross( _normal, sid ) );



  return (
    sid * sin( theta ) * cos( phi )
    + top * sin( theta ) * sin( phi )
    + _normal * cos( theta )
  );
}

// ------

struct Camera {
  vec3 pos;
  vec3 dir;
  vec3 sid;
  vec3 top;
};

struct Ray {
  vec3 dir;
  vec3 ori;
  int mtl;
};

struct Map {
  float dist;
  int mtl;
  vec4 props;
};

struct March {
  Ray ray;
  Map map;
  float len;
  vec3 pos;
  vec3 normal;
  float edge;
};

struct Material {
  vec3 color;
  vec3 emissive;
  vec3 edgeEmissive;
  float reflective;
  float reflectiveRoughness;
  float refractive;
  float refractiveRoughness;
  float refractiveIndex;
};

// ------

Camera camInit( in vec3 _pos, in vec3 _tar ) {
  Camera cam;
  cam.pos = _pos;
  cam.dir = normalize( _tar - _pos );
  cam.sid = normalize( cross( cam.dir, V.xyx ) );
  cam.top = normalize( cross( cam.sid, cam.dir ) );

  return cam;
}

Ray rayInit( in vec3 _ori, in vec3 _dir, in int _mtl ) {
  Ray ray;
  ray.dir = _dir;
  ray.ori = _ori;
  ray.mtl = _mtl;
  return ray;
}

Ray rayFromCam( in vec2 _p, in Camera _cam ) {
  vec3 dir = normalize( _p.x * _cam.sid + _p.y * _cam.top + _cam.dir * 2.0 * ( 1.0 - length( _p.xy ) * 0.3 ) );
  return rayInit( _cam.pos, dir, 1 );
}

Map mapInit( in float _dist ) {
  Map map;
  map.dist = _dist;
  map.mtl = 1;
  return map;
}

March marchInit( in Ray _ray ) {
  March march;
  march.ray = _ray;
  march.len = INIT_LEN;
  march.pos = _ray.ori + _ray.dir * march.len;
  march.normal = V.xxy;
  march.edge = 0.0;
  return march;
}

Material mtlInit() {
  Material material;
  material.color = V.xxx;
  material.emissive = V.xxx;
  material.edgeEmissive = V.xxx;
  material.reflective = 0.0;
  material.reflectiveRoughness = 0.0;
  material.refractive = 0.0;
  material.refractiveRoughness = 0.0;
  material.refractiveIndex = 1.0;
  return material;
}

// ------

float sphere( in vec3 _p, in float _r ) {
  return length( _p ) - _r;
}

vec3 mandelbulb( vec3 _p ) {
	vec3 p = _p.xzy;
	vec3 z = p;
	vec3 dz = vec3( 0.0 );
	float power = 8.0;
	float r, theta, phi;
	float dr = 1.0;

	float t0 = 1.0;

	for( int i = 0; i < 7; i ++ ) {
		r = length( z );
		if( r > 2.0 ) { continue; }
		theta = atan( z.y / z.x );
		phi = asin( z.z / r );

		dr = pow( r, power - 1.0 ) * dr * power + 1.0;

		r = pow( r, power );
		theta = theta * power;
		phi = phi * power;

		z = r * vec3( cos( theta ) * cos( phi ), sin( theta ) * cos( phi ), sin( phi ) ) + p;

		t0 = min( t0, r );
	}
	return vec3( 0.5 * log( r ) * r / dr, t0, 0.0 );
}

float metaball( in vec3 _p ) {
  vec3 p = _p;
  float dist = 1E9;

  for ( int i = 0; i < 6; i ++ ) {
    vec3 translate = vec3(
      sin( 3.0 + time * mod( float( i + 3 ) * 3.0, 4.0 ) * PI * 2.0 + float( i ) ),
      sin( 1.0 + time * mod( float( i + 1 ) * 3.0, 5.0 ) * PI * 2.0 + float( i ) ),
      sin( 5.0 + time * mod( float( i + 1 ) * 4.0, 3.0 ) * PI * 2.0 + float( i ) )
    ) * 0.15;
    p = _p - translate;
    dist = smin( dist, length( p ) - 0.12, 0.1 );
  }

  return dist;
}

float box( in vec3 _p, in vec3 _size ) {
  vec3 d = abs( _p ) - _size;
  return min( max( d.x, max( d.y, d.z ) ), 0.0 ) + length( max( d, 0.0 ) );
}

Map distFunc( in vec3 _p, in int _mtl ) {
  Map map = mapInit( 1E9 );
  vec3 pp = _p;

  { // mandelbulb
    vec3 p = pp;
    vec3 res = mandelbulb( p / 9.0 );
    res.x = 9.0 * ( ( _mtl == MTL_MANDELBULB ) ? -res.x : res.x );

    if ( res.x < map.dist ) {
      map = mapInit( res.x );
      map.mtl = ( _mtl == MTL_MANDELBULB ) ? MTL_AIR : MTL_MANDELBULB;
      map.props.x = res.y;
    }
  }

  { // light
    vec3 p = pp;
    p.yz = rotate2D( -lofi( atan( p.z, p.y ) + PI / 12.0, PI / 6.0 ) ) * p.yz;
    p.xy = rotate2D( -lofi( atan( p.y, p.x ), PI / 6.0 ) - PI / 12.0 ) * p.xy;

    {
      float dist = box( p, vec3( 11.1, 0.1, 0.1 ) );
      dist *= ( _mtl == MTL_TUBE ) ? -1.0 : 1.0;

      if ( dist < map.dist ) {
        map = mapInit( dist );
        map.mtl = ( _mtl == MTL_TUBE ) ? MTL_AIR_FROM_TUBE : MTL_TUBE;
        map.props.x = 10.2 < p.x ? 1.0 : 0.0;
      }
    }

    if ( _mtl == MTL_TUBE ) {
      float dist = box( p, vec3( 11.0, 0.04, 0.04 ) );

      if ( dist < map.dist ) {
        map = mapInit( dist );
        map.mtl = MTL_LIGHT;
      }
    }
  }

  { // water
    vec3 p = pp;
    float dist = sphere( p, 8.0 );
    dist *= ( _mtl == MTL_WATER ) ? -1.0 : 1.0;

    if ( dist < map.dist ) {
      map = mapInit( dist );
      map.mtl = ( _mtl == MTL_WATER ) ? MTL_AIR : MTL_WATER;
    }
  }

  { // sky
    vec3 p = pp;
    float dist = -sphere( p, 50.0 );

    if ( dist < map.dist ) {
      map = mapInit( dist );
      map.mtl = MTL_SKY;
    }
  }

  { // ship
    vec3 p = pp;
    p.yz = rotate2D( time * PI ) * p.yz;
    p.y -= 10.4;
    p.yz = rotate2D( 0.2 ) * p.yz;

    {
      vec3 p2 = p;
      p2.xy = rotate2D( time * PI ) * p2.xy;
      p2.xy = rotate2D( -lofi( atan( p2.y, p2.x ) + PI / 12.0, PI / 6.0 ) ) * p2.xy;
      p2.x -= 0.5;
      float dist = box( p2, vec3( 0.1, 0.03, 0.4 ) );
      dist *= ( _mtl == MTL_SHIP_SHELL ) ? -1.0 : 1.0;

      if ( dist < map.dist ) {
        map = mapInit( dist );
        map.mtl = ( _mtl == MTL_SHIP_SHELL ) ? MTL_AIR : MTL_SHIP_SHELL;
      }
    }

    {
      vec3 p2 = p;
      p2.zx = rotate2D( -time * PI * 4.0 ) * p2.zx;
      float dist = max( length( p2.zx ) - 0.2, -( length( p2.zx ) - 0.1 ) );
      dist = max( dist, -( abs( p2.x ) - 0.04 ) );
      dist = max( dist, abs( p2.y ) - 0.01 );

      if ( dist < map.dist ) {
        map = mapInit( dist );
        map.mtl = MTL_SHIP_CORE;
      }
    }
  }

  return map;
}

vec3 normalFunc( in vec3 _p, in float _d, in int _mtl ) {
  vec2 d = V * _d;
  return normalize( vec3(
    distFunc( _p + d.yxx, _mtl ).dist - distFunc( _p - d.yxx, _mtl ).dist,
    distFunc( _p + d.xyx, _mtl ).dist - distFunc( _p - d.xyx, _mtl ).dist,
    distFunc( _p + d.xxy, _mtl ).dist - distFunc( _p - d.xxy, _mtl ).dist
  ) );
}

// ------

March march( in Ray _ray ) {
  Ray ray = _ray;
  March march = marchInit( ray );
  float fogLen = 0.0;
  bool fog = false;

  for ( int iMarch = 0; iMarch < MARCH_ITER; iMarch ++ ) {
    Map map = distFunc( march.pos, ray.mtl );
    map.dist *= 0.8;

    march.map = map;
    march.len += map.dist;
    march.pos = ray.ori + ray.dir * march.len;

    if ( ray.mtl == MTL_AIR ) {
      for ( int i = 0; i < 99; i ++ ) {
        if ( march.len < fogLen ) { break; }
        fogLen += random() * FOG_LENGTH;
        if ( random() < FOG_THRESHOLD ) {
          fog = true;
          break;
        }
      }
    }

    if ( fog || 1E3 < march.len || abs( map.dist ) < INIT_LEN * 0.01 ) { break; }
  }

  if ( fog ) {
    march.len = fogLen;
    march.pos = ray.ori + ray.dir * march.len;
    march.normal = randomHemisphere( -ray.dir );
    march.map.dist = 0.0;
    march.map.mtl = MTL_FOG;
  } else {
    march.normal = normalFunc( march.pos, 1E-4, ray.mtl );
    march.edge = 1.0 - smoothstep( 0.9, 0.98, dot( normalFunc( march.pos, 4E-4, ray.mtl ), march.normal ) );
  }

  return march;
}

// ------

Material getMtl( int _mtl, vec4 _props ) {
  Material mtl = mtlInit();

  if ( _mtl == MTL_AIR ) {
    mtl.color = vec3( 1.0 );
    mtl.refractive = 1.0;
    mtl.refractiveIndex = 1.0;

  } else if ( _mtl == MTL_FOG ) {
    mtl.color = vec3( 0.7 );

  } else if ( _mtl == MTL_LIGHT ) {
    mtl.emissive = 100.0 * vec3( 1.0, 0.7, 0.9 );

  } else if ( _mtl == MTL_TUBE ) {
    if ( _props.x == 0.0 ) {
      mtl.color = vec3( 0.3 );
      mtl.reflective = 0.4;
    } else {
      mtl.color = vec3( 0.94 );
      mtl.refractive = 0.8;
      mtl.reflective = 0.1;
    }

  } else if ( _mtl == MTL_AIR_FROM_TUBE ) {
    if ( _props.x == 0.0 ) {
      mtl.color = vec3( 0.01 );
    } else {
      mtl.color = vec3( 1.0 );
      mtl.refractive = 1.0;
      mtl.refractiveIndex = 1.0;
    }

  } else if ( _mtl == MTL_MANDELBULB ) {
    mtl.color = mix(
      vec3( 0.5, 0.9, 0.1 ),
      vec3( 0.9, 0.1, 0.3 ),
      _props.x
    );

  } else if ( _mtl == MTL_WATER ) {
    mtl.color = vec3( 0.87, 0.91, 0.98 );
    mtl.emissive = 0.3 * vec3( 0.3, 0.5, 0.8 );
    mtl.refractive = 0.97;
    mtl.refractiveIndex = 1.6;
    mtl.reflective = 0.03;

  } else if ( _mtl == MTL_SKY ) {
    mtl.color = vec3( 0.0 );
    mtl.emissive = vec3( 0.04, 0.07, 0.1 );

  } else if ( _mtl == MTL_SHIP_SHELL ) {
    mtl.color = vec3( 0.9, 0.4, 0.5 );
    mtl.reflective = 0.3;
    mtl.refractive = 0.4;
    mtl.refractiveIndex = 1.3;

  } else if ( _mtl == MTL_SHIP_CORE ) {
    mtl.emissive = 400.0 * vec3( 0.2, 0.6, 0.9 );

  }

  return mtl;
}

// ------

Ray shade( in March _march, inout vec3 colorAdd, inout vec3 colorMul ) {
  March march = _march;

  if ( abs( march.map.dist ) < 1E-2 ) {
    vec3 normal = march.normal;
    float edge = march.edge;

    int rayMtl = march.ray.mtl;
    Material material = getMtl( march.map.mtl, march.map.props );

    vec3 dir = V.xxx;
    float dice = random4().x;

    // colorAdd += colorMul * max( 0.0, dot( normal, -march.ray.dir ) ) * march.map.material.emissive;
    colorAdd += colorMul * material.emissive;
    colorAdd += colorMul * edge * material.edgeEmissive;

    colorMul *= material.color;
    if ( dice < material.reflective ) { // reflect
      vec3 ref = normalize( reflect(
        march.ray.dir,
        normal
      ) );
      vec3 dif = randomHemisphere( normal );
      dir = normalize( mix(
        ref,
        dif,
        material.reflectiveRoughness
      ) );
      colorMul *= max( dot( dir, ref ), 0.0 );

    } else if ( dice < material.reflective + material.refractive ) { // refract
      vec3 inc = normalize( march.ray.dir );
      float eta = getMtl( march.ray.mtl, V.xxxx ).refractiveIndex / material.refractiveIndex;

      vec3 ref = refract( inc, normal, eta );
      ref = ( ref == V.xxx )
      ? ( normalize( reflect(
        march.ray.dir,
        normal
      ) ) )
      : normalize( ref );

      vec3 dif = randomHemisphere( -normal );
      dir = normalize( mix(
        ref,
        dif,
        material.refractiveRoughness
      ) );
      colorMul *= max( dot( dir, ref ), 0.0 );

      rayMtl = march.map.mtl;

    } else { // diffuse
      dir = randomHemisphereCosWeighted( normal );
      colorMul *= 1.0;
    }

    Ray ray = rayInit( march.pos, dir, rayMtl );
    return ray;
  } else {
    colorMul *= 0.0;

    return rayInit( V.xxy, V.xxy, MTL_AIR );
  }
}

// ------

void main() {
  vec2 uv = gl_FragCoord.xy / resolution;
  seed = texture2D( textureRandom, gl_FragCoord.xy / resolution );

  vec4 tex0 = texture2D( textureDrawBuffers0, uv );
  vec4 tex1 = texture2D( textureDrawBuffers1, uv );
  vec4 tex2 = texture2D( textureDrawBuffers2, uv );
  vec4 tex3 = texture2D( textureDrawBuffers3, uv );

  vec3 colorAdd = abs( tex1.xyz ) - 1E-2;
  vec3 colorMul = abs( tex2.xyz ) - 1E-2;
  vec3 colorOut = tex3.xyz;
  int rayMtl = int( abs( tex2.w ) );
  float depth = ( tex1.x < 0.0 ? 0.0 : 1.0 ) + ( tex1.y < 0.0 ? 0.0 : 2.0 ) + ( tex1.z < 0.0 ? 0.0 : 4.0 ) + ( tex2.x < 0.0 ? 0.0 : 8.0 ) + ( tex2.y < 0.0 ? 0.0 : 16.0 ) + ( tex2.z < 0.0 ? 0.0 : 32.0 );
  float samples = abs( tex3.w );

  Ray ray;
  vec3 dir = vec3( tex0.w, tex1.w, 0.0 );
  dir.z = sqrt( 1.0 - dot( dir, dir ) ) * sign( tex2.w );
  ray = rayInit( tex0.xyz, dir, rayMtl );

  if ( reset ) {
    colorOut = V.xxx;
    colorAdd = V.xxx;
    samples = 0.0;
  }

  for ( int i = 0; i < REFLECT_PER_PATH; i ++ ) {

    if ( reset || REFLECT_MAX <= depth || length( colorMul ) < RAYAMP_MIN ) {
      samples += 1.0;
      depth = 0.0;

      colorOut = mix(
        colorOut,
        max( V.xxx, colorAdd ),
        1.0 / samples
      );

      // ------

      float camRot = -time * PI;
      vec3 camTar = vec3( 0.0, 10.4, 0.0 );
      Camera cam = camInit(
        vec3( 1.0, 12.5, 2.5 ),
        camTar
      );

      // dof
      vec2 dofCirc = randomCircle() * 0.01;
      cam.pos += dofCirc.x * cam.sid;
      cam.pos += dofCirc.y * cam.top;

      cam = camInit( cam.pos, camTar );

      cam.pos.yz = rotate2D( camRot ) * cam.pos.yz;
      cam.dir.yz = rotate2D( camRot ) * cam.dir.yz;
      cam.sid.yz = rotate2D( camRot ) * cam.sid.yz;
      cam.top.yz = rotate2D( camRot ) * cam.top.yz;

      // antialias
      vec2 pix = gl_FragCoord.xy + random4().xy - 0.5;

      vec2 p = ( pix * 2.0 - resolution ) / resolution.x;
      ray = rayFromCam( p, cam );

      colorAdd = V.xxx;
      colorMul = V.yyy;
    } else {
      depth += 1.0;
    }

    March m = march( ray );
    ray = shade( m, colorAdd, colorMul );

  }

  // ------

  vec3 depthBits1 = vec3(
    mod( depth, 2.0 ) < 1.0 ? -1.0 : 1.0,
    mod( depth / 2.0, 2.0 ) < 1.0 ? -1.0 : 1.0,
    mod( depth / 4.0, 2.0 ) < 1.0 ? -1.0 : 1.0
  );

  vec3 depthBits2 = vec3(
    mod( depth / 8.0, 2.0 ) < 1.0 ? -1.0 : 1.0,
    mod( depth / 16.0, 2.0 ) < 1.0 ? -1.0 : 1.0,
    mod( depth / 32.0, 2.0 ) < 1.0 ? -1.0 : 1.0
  );

  gl_FragData[ 0 ] = vec4( ray.ori, ray.dir.x );
  gl_FragData[ 1 ] = vec4( ( colorAdd + 1E-2 ) * depthBits1, ray.dir.y );
  gl_FragData[ 2 ] = vec4( ( colorMul + 1E-2 ) * depthBits2, float( ray.mtl ) * ( ( 0.0 < ray.dir.z ) ? 1.0 : -1.0 ) );
  gl_FragData[ 3 ] = vec4( colorOut, samples );
}
