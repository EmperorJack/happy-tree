#version 120

// Constant across both shaders
uniform sampler2D texture0;
uniform bool useTexture;
uniform bool useLighting;

// Values passed in from the vertex shader
varying vec2 vTextureCoord0;
varying vec3 n;
varying vec3 v;

// Method for determining if fragment within spotlight
bool withinSpotlight(int, vec3, float);

void main() {
	vec4 finalColor = vec4(0, 0, 0, 0);

  if(useLighting){

    int spotlightIndex = -1;

    // Compute the illumination of each light
    for (int lightIndex = 0; lightIndex < gl_MaxLights; lightIndex++) {
      vec3 color = vec3(0, 0, 0);

      vec3 l = normalize(gl_LightSource[lightIndex].position.xyz - v);
      vec3 e = normalize(-v);
      vec3 r = normalize(reflect(-l, n));

      float s_dot_n = max(dot(l, n), 0.0);

      // Ambient
      vec3 ambient = gl_LightSource[lightIndex].ambient.rgb *
                     gl_FrontMaterial.ambient.rgb;

      // Diffuse
      vec3 diffuse = gl_LightSource[lightIndex].diffuse.rgb *
                     gl_FrontMaterial.diffuse.rgb *
                     s_dot_n;

      // Specular
      vec3 specular = vec3(0.0, 0.0, 0.0);
      if (s_dot_n > 0.0) {
        specular = gl_LightSource[lightIndex].specular.rgb *
                   gl_FrontMaterial.specular.rgb *
                   pow(max(dot(r, normalize(-v)), 0.0), gl_FrontMaterial.shininess);
      }

      if (lightIndex == spotlightIndex) {
        if (withinSpotlight(lightIndex, l, s_dot_n)) {
          color = ambient + diffuse + specular;
        }
      } else {
        color = ambient + diffuse + specular;
      }

      finalColor += vec4(color, 1);
    }
  } else{
    finalColor = gl_FrontMaterial.ambient.rgba * gl_FrontMaterial.diffuse.rgba;
  }

  if (useTexture) {
    finalColor *= texture2D(texture0, vTextureCoord0).rgba;
  }

  gl_FragColor = finalColor;
}

bool withinSpotlight(int lightIndex, vec3 l, float s_dot_n) {
  if (s_dot_n > 0.0) {
    float angle = dot(normalize(gl_LightSource[lightIndex].spotDirection), normalize(-l));
    if (acos(angle) < radians(gl_LightSource[lightIndex].spotCutoff)) {
      return true;
    }
  }

  return false;
}
