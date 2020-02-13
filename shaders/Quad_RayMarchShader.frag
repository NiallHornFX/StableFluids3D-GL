#version 420 core 

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

uniform float mx; // 2D Mouse Input. 
uniform float my;

vec4 desat (vec3 color, float DD)
{
	vec3 Aa = vec3(0.3, 0.59, 0.11);
	vec3 Ab = vec3(dot(Aa, color));
	return vec4(mix(color, Ab, DD), 1.0);
}

/* ----------------------------------------------------------------- */
void main()
{
	// Map from 0-N FragCoord_Space to 0-1 UV Space. 
	vec2 uv = (gl_FragCoord.xy - 0.0) / 512; 
	
	// Sample Baked Rasterized Cube Textures for Entry/Exit Postions -
	vec3 samp_cf_start = texture(cf_tex, uv).rgb; 
	vec3 samp_cb_end = texture(cb_tex, uv).rgb; 
	
	vec3 ray = samp_cb_end - samp_cf_start;
	float rayL = length(ray); 
	vec3 ray_dir = normalize(ray); 
	
	// Backwards March.
	vec3 ray_P = samp_cb_end;
	
	// Map 3D Texture Space Cube Local Space Texture Sampled Locations.
	float dens = 0.0; vec3 vel = vec3(0.0, 0.0, 0.0); 
	int step_count = 100, l_step_count = 30; 
	vec3 step = (samp_cf_start - samp_cb_end) / (step_count - 1);
	
	for (int i = 0; i < step_count; i++)
	{
		dens += texture(dens_tex, ray_P).x;
		vel += texture(vel_tex, ray_P).xyz; 
		ray_P += step; 
		
		// Early Ray Termination 
		//if (dens > 0.99) { dens = 1.0; break;}
	}	

	vel /= step_count; vel *= clamp(dens, 0.0, 1.0); 
	dens /= float(step_count); dens *= 5.0;
	float dens_a = dens; //dens_f + dens_b;
	
	// Colour Using BlackBody Lookup of vel_l
	float vel_l = length(vel); 
	vec3 vel_col = (vec4(1,1./4.,1./16.,1) * exp(4.*vel_l - 1.)).xyz;
	vel_col *= dens_a * 35.0; 
	
	// Colour Density + Cube. 
	vec3 dens_vec = vec3(dens_a + vel.x, dens_a + vel.y, dens_a + vel.z); 
	vec3 cv_0 = mix(samp_cf_start.xyz, samp_cb_end.xyz, 0.5);
	cv_0 = desat(cv_0, 1.0).xyz; cv_0.xy += 0.1 * length(cv_0); 
	cv_0 *= 0.1; // Grid Opac. 
	
	// Dens - 
	vec3 cv_1 = mix(cv_0, dens_vec, 0.5); 
	// Vel BB - 
	//vec3 cv_1 = mix(cv_0, vel_col, 0.5);
	
	// Final FragColor. 
	frag_color = vec4(clamp(cv_1, 0.0, 1.0), 1.0); 
	

	//frag_color = vec4(cv_0, 1.0); // Cube Only Test. 
	//frag_color = vec4(dens_vec, 1.0); // Dens Only.  	
}
