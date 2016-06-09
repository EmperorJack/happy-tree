#version 120

// Constant across both shaders
uniform sampler2D texture0;
uniform bool useTexture;
uniform bool useLighting;

// Values to pass to the fragment shader
varying vec2 vTextureCoord0;
varying vec3 n;
varying vec3 v;

void main() {
	vTextureCoord0 = gl_MultiTexCoord0.xy;

  v = (gl_ModelViewMatrix * gl_Vertex).xyz;
  n = normalize(gl_NormalMatrix * gl_Normal);

	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}

