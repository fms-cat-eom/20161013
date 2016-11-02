#define MARCH_ITER 200
#define RAYAMP_MIN 0.01
#define REFLECT_MAX 8.0
#define REFLECT_PER_PATH 3
#define INIT_LEN 0.01
#define SKY_COLOR vec3( 1.0, 1.0, 1.0 )

// -----

#define MTL_AIR 1
#define MTL_AIR_FROM_COOL 2
#define MTL_COOL 3
#define MTL_HOT 4
#define MTL_LIGHT 5

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
  bool anotherWorld;
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
  float anotherWorld;
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
  ray.anotherWorld = false;
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
  material.anotherWorld = 0.0;
  return material;
}

// ------

float sphere( in vec3 _p, in float _r ) {
  return length( _p ) - _r;
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

float slasher( vec3 _p, float _ratio ) {
  float phase = ( _p.x + _p.y );
  float slash = abs( 0.5 - ( phase - floor( phase ) ) ) * 2.0;
  return ( slash - _ratio ) / sqrt( 2.0 );
}

vec3 typeIfs( vec3 _p, vec3 _rot, vec3 _shift ) {
  vec3 pos = _p;

  vec3 shift = _shift;

  for ( int i = 0; i < 5; i ++ ) {
    float intensity = pow( 2.0, -float( i ) );

    pos.y -= 0.0;

    pos = abs( pos ) - shift * intensity;

    shift.yz = rotate2D( _rot.x ) * shift.yz;
    shift.zx = rotate2D( _rot.y ) * shift.zx;
    shift.xy = rotate2D( _rot.z ) * shift.xy;

    if ( pos.x < pos.y ) { pos.xy = pos.yx; }
    if ( pos.x < pos.z ) { pos.xz = pos.zx; }
    if ( pos.y < pos.z ) { pos.yz = pos.zy; }
  }

  return pos;
}

Map distFunc( in vec3 _p, in int _mtl ) {
  Map map = mapInit( 1E9 );

  { // cube
    vec3 p = _p;
    p.yz = rotate2D( 0.5 + time * PI * 2.0 ) * p.yz;
    p.zx = rotate2D( 0.5 ) * p.zx;
    p.xy = rotate2D( time * PI * 2.0 ) * p.xy;
    float dist = box( p, vec3( 0.4 ) );
    dist *= _mtl == MTL_COOL ? -1.0 : 1.0;

    if ( dist < map.dist ) {
      map = mapInit( dist );
      map.mtl = _mtl == MTL_COOL ? MTL_AIR_FROM_COOL : MTL_COOL;
    }
  }

  vec3 pc = _p;
  pc.xy = rotate2D( sin( time * PI * 2.0 ) * 0.12 ) * pc.xy;
  pc.xz = rotate2D( time * PI * 2.0 ) * pc.xz;

  { // ifs
    vec3 p = pc;
    float dist = -box( p, vec3( 3.0 ) );

    if ( dist < map.dist ) {
      map = mapInit( dist );
      map.mtl = MTL_HOT;
    }
  }

  { // light
    vec3 p = pc;
    p = abs( p );
    p -= vec3( 2.0 );
    float dist = sphere( p, 0.5 );

    if ( dist < map.dist ) {
      map = mapInit( dist );
      map.mtl = MTL_LIGHT;
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

  for ( int iMarch = 0; iMarch < MARCH_ITER; iMarch ++ ) {
    Map map = distFunc( march.pos, ray.mtl );
    map.dist *= 0.8;

    march.map = map;
    march.len += map.dist;
    march.pos = ray.ori + ray.dir * march.len;

    if ( 1E3 < march.len || abs( map.dist ) < INIT_LEN * 0.01 ) { break; }
  }

  march.normal = normalFunc( march.pos, 1E-4, ray.mtl );
  march.edge = 1.0 - smoothstep( 0.9, 0.98, dot( normalFunc( march.pos, 4E-4, ray.mtl ), march.normal ) );

  return march;
}

// ------

Material getMtl( int _mtl, vec4 _props ) {
  Material mtl = mtlInit();

  if ( _mtl == MTL_AIR ) {
    mtl.refractive = 1.0;
    mtl.refractiveIndex = 1.0;

  } else if ( _mtl == MTL_AIR_FROM_COOL ) {
    mtl.color = vec3( 1.0 );
    mtl.refractive = 1.0;
    mtl.refractiveRoughness = 0.03;
    mtl.refractiveIndex = 1.2;
    mtl.edgeEmissive = vec3( 0.1, 0.9, 0.3 ) * 400.0;

  } else if ( _mtl == MTL_COOL ) {
    mtl.color = vec3( 0.5, 0.8, 0.9 );
    mtl.refractive = 1.0;
    mtl.refractiveIndex = 1.0;
    mtl.edgeEmissive = vec3( 0.1, 0.9, 0.3 ) * 400.0;

  } else if ( _mtl == MTL_HOT ) {
    mtl.color = vec3( 0.7 );
    mtl.reflective = 0.8;

  } else if ( _mtl == MTL_LIGHT ) {
    mtl.emissive = vec3( 1.0 );

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
    bool rayAW = march.ray.anotherWorld;
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
      bool toAir = ( 0.0 < dot( normal, inc ) );
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
      dir = randomHemisphere( normal );
      colorMul *= max( dot( dir, normal ), 0.0 );
    }

    Ray ray = rayInit( march.pos, dir, rayMtl );
    ray.anotherWorld = rayAW;
    return ray;
  } else {
    colorAdd += colorMul * SKY_COLOR;
    colorMul *= 0.0;

    return rayInit( V.xxy, V.xxy, 1 );
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

  vec3 colorAdd = tex1.xyz;
  vec3 colorMul = abs( tex2.xyz ) - 1E-2;
  vec3 colorOut = tex3.xyz;
  int rayMtl = int( abs( tex2.w ) );
  bool rayAw = tex3.w < 0.0;
  float depth = ( tex2.x < 0.0 ? 0.0 : 1.0 ) + ( tex2.y < 0.0 ? 0.0 : 2.0 ) + ( tex2.z < 0.0 ? 0.0 : 4.0 );
  float samples = abs( tex3.w );

  Ray ray;
  vec3 dir = vec3( tex0.w, tex1.w, 0.0 );
  dir.z = sqrt( 1.0 - dot( dir, dir ) ) * sign( tex2.w );
  ray = rayInit( tex0.xyz, dir, rayMtl );
  ray.anotherWorld = rayAw;

  if ( reset ) {
    colorOut = V.xxx;
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

      float camRot = time * PI * 2.0;
      Camera cam = camInit(
        vec3( 0.0, 0.0, 2.0 ),
        vec3( 0.0, 0.0, 0.0 )
      );

      // dof
      vec2 dofCirc = randomCircle() * 0.01;
      cam.pos += dofCirc.x * cam.sid;
      cam.pos += dofCirc.y * cam.top;

      cam = camInit( cam.pos, vec3( 0.0, 0.0, 0.0 ) );

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

  vec3 depthBits = vec3(
    mod( depth, 2.0 ) < 1.0 ? 1.0 : -1.0,
    mod( depth / 2.0, 2.0 ) < 1.0 ? 1.0 : -1.0,
    mod( depth / 4.0, 2.0 ) < 1.0 ? 1.0 : -1.0
  );
  depthBits = V.yyy;

  gl_FragData[ 0 ] = vec4( ray.ori, ray.dir.x );
  gl_FragData[ 1 ] = vec4( colorAdd, ray.dir.y );
  gl_FragData[ 2 ] = vec4( ( colorMul + 1E-2 ) * depthBits, float( ray.mtl ) * ( ( 0.0 < ray.dir.z ) ? 1.0 : -1.0 ) );
  gl_FragData[ 3 ] = vec4( colorOut, samples * ( ray.anotherWorld ? -1.0 : 1.0 ) );
}
