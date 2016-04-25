attribute vec2 coord2d;
attribute vec3 v_color;

varying vec3 f_color;
varying vec2  uv;

void main(void) {
  //gl_Position = vec4(coord2d, 0.0, 1.0);
  
  gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * vec4(coord2d, 0.0, 1.0);
  f_color = v_color;
  
  // Texture coordinate for screen aligned (in correct range):
  uv = (vec2( gl_Position.x, gl_Position.y ) + vec2( 1.0 ) ) / vec2( 2.0 );
  
}
