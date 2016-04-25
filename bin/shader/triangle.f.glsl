uniform float fade;

varying vec3 f_color;
varying vec2 uv;
varying vec2 TexCoord;

// some const, tweak for best look
const float sampleStrength = 1.8;

void main(void) {

  vec4 resultColor;
  
  // original
  // j first one
  resultColor = vec4(f_color.x,f_color.y,f_color.z, fade);
  //resultColor += 0.0;
  
  gl_FragColor = 1.1*((resultColor-0.5)*1.0 + 0.5);
  


  // this one has no fade value - pass thru
  //gl_FragColor = vec4(f_color.x, f_color.y, f_color.z, 1.0);

}
