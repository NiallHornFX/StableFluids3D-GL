#version 400 core 

in vec4 vpos; 
in vec4 gl_FragCoord; 

out vec4 frag_color; 

// 2D Samplers \\
// Unit Cube Textures -   // TU
uniform sampler2D cf_tex; // 0
uniform sampler2D cb_tex; // 1

// 3D Samplers \\
// 3D Texture Grids - 
uniform sampler2D dens_tex; // 2
uniform sampler2D vel_tex;  // 3


/* ----------------------------------------------------------------- */

void main()
{
	// Map from 0-N FragCoord_Space to 0-1 UV Space. 
	vec2 uv = (gl_FragCoord.xy - 0.0) / 512; 
	vec4 samp_cf = texture(cf_tex, uv); 
	vec4 samp_cb = texture(cb_tex, uv); 
	frag_color = vec4(mix(samp_cf.xyz, samp_cb.xyz, 0.5), 1.0); 
	//frag_color = vec4(uv.xy, 0.0, 1.0); 
	
}
