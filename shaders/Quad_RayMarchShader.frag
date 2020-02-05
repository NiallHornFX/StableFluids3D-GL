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
	
	// Sample Baked Cube Textures for Entry/Exit Postions -
	vec4 samp_cf = texture(cf_tex, uv); 
	vec4 samp_cb = texture(cb_tex, uv); 
	
	// Ray Inital Postion From Sampled Cube Texture RGB-Local XYZ Postion. 
	//vec3 ray_P = vec3(samp_cb.xyz); 
	vec3 ray_P = vec3(samp_cf.xyz); 
	
	//vec3 ray_P = vec3(uv, 0.0); 
	//vec3 ray_dir = normalize(samp_cf.xyz - samp_cb.xyz); // Dir From FrontFace-BackFace
	vec3 ray_dir = normalize(samp_cb.xyz - samp_cf.xyz); // Dir From BackFaceCol - FrontFaceCol
	//vec3 ray_dir = vec3(0.0, 0.0, -1.0); 
	ray_P = clamp(ray_P, 0.0, 1.0); // Clamp Ray_P to 3D Texture Space 0-1,xyz. 

	// Map 3D Texture Space Cube Local Space Texture Sampled Locations.
	float dens = 0.0; 
	int step_count = 100;
	float step_size = 0.01; /// step_count; 
	
	for (int i = 0; i < step_count; i++)
	{
		dens += texture(dens_tex, ray_P).x;
		ray_P += ray_dir * step_size;
		if (dens > 0.99) {break;}
	}		

	
	//frag_color = vec4(mix(samp_cf.xyz, samp_cb.xyz, 0.5), 1.0); // Check Cube Faces Blended. 
	vec3 dens_vec = vec3(dens, dens, dens); 
	vec3 col_vec = mix(samp_cf.xyz, dens_vec, 0.5);
	frag_color = vec4(col_vec, 1.0); 
	//frag_color = vec4(dens, dens, dens, 1.0); 
	//frag_color = vec4(ray_dir, 1.0); // Check RayDir
	//frag_color = vec4(uv.xy, 0.0, 1.0); // Check ScreenSpace UV.
	
}
