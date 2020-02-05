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
uniform sampler3D dens_tex; // 2
uniform sampler3D vel_tex;  // 3


/* ----------------------------------------------------------------- */

void main()
{
	// Map from 0-N FragCoord_Space to 0-1 UV Space. 
	vec2 uv = (gl_FragCoord.xy - 0.0) / 512; 
	
	// Sample Baked Cube Textures -
	vec4 samp_cf = texture(cf_tex, uv); 
	vec4 samp_cb = texture(cb_tex, uv); 
	
	//vec3 ray_P = vec3(samp_cf.xyz); 
	vec3 ray_P = vec3(uv, 0.0); 
	//vec3 ray_dir = normalize(samp_cb.xyz - samp_cf.xyz); // Dir From BackFaceCol - FrontFaceCol
	vec3 ray_dir = vec3(0.0, 0.0, -1.0); 

	float dens = 0.0; 
	int step_count = 100;
	float step_size = 0.1 / step_count; 
	
	for (int i = 0; i < step_count; i++)
	{
		dens += texture(dens_tex, ray_P).x;
		ray_P += ray_dir * step_size;
		if (dens > 0.99) {break;}
	}		

	
	//frag_color = vec4(mix(samp_cf.xyz, samp_cb.xyz, 0.5), 1.0); // Check Cube Faces Blended. 
	frag_color = vec4(dens, dens, dens, 1.0); 
	//frag_color = vec4(ray_dir, 1.0); // Check RayDir
	//frag_color = vec4(uv.xy, 0.0, 1.0); // Check ScreenSpace UV.
	
}
