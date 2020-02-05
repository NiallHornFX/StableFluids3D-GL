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

uniform float mx; // 2D Mouse Input. 
uniform float my;

vec4 desat (vec3 color, float DD)
{
	vec3 grayXfer = vec3(0.3, 0.59, 0.11);
	vec3 gray = vec3(dot(grayXfer, color));
	return vec4(mix(color, gray, DD), 1.0);
}

/* ----------------------------------------------------------------- */
// BIN SHADER
void main()
{
	// Map from 0-N FragCoord_Space to 0-1 UV Space. 
	vec2 uv = (gl_FragCoord.xy - 0.0) / 512; 
	
	// Sample Baked Cube Textures for Entry/Exit Postions -
	vec3 samp_cf_start = texture(cf_tex, uv).rgb; 
	vec3 samp_cb_end = texture(cb_tex, uv).rgb; 
	vec3 ray = samp_cb_end - samp_cf_start;
	float rayL = length(ray); 
	vec3 ray_dir = normalize(ray); 
	
	// Ray Inital Postion From Sampled Cube Texture RGB-Local XYZ Postion. 
	vec3 ray_P = samp_cf_start;
	
	// Clamp Ray_P to 3D Texture Space 0-1,xyz. 
	//ray_P_f = clamp(ray_P_f, 0.0, 1.0); 
	//ray_P_b = clamp(ray_P_b, 0.0, 1.0); 

	// Map 3D Texture Space Cube Local Space Texture Sampled Locations.
	float dens = 0.0; 
	int step_count = 100, l_step_count = 30; 
	float step_size = 0.01; /// step_count; 
	
	for (int i = 0; i < step_count; i++)
	{
		dens += texture(dens_tex, ray_P).x;
		ray_P += ray_dir * step_size;
		float l = 0.0; 
		vec3 lray_P = ray_P; 
		for (int j = 0; j < l_step_count; j++) // Basic Shadow Ray...
		{
			//vec3 lpos = vec3(0.9, 1.0, -0.2);
			//vec3 lpos = vec3(mx, my, -clamp(mx+my, 0.0, 1.0)); 
			vec3 lpos = vec3(mx, my, -0.8); 
			vec3 ldir = normalize(lpos - ray_P); 
			l += texture(dens_tex, lray_P).x * 0.25;
			lray_P += ldir * step_size; 
			if (l >= 0.99) {break;}
		}
		
		l /= l_step_count; l = clamp(l, 0.0, 1.0); // Hacky Clamping... Scaled into SDR.
		dens -= l * 1.0; // Subtract ShadowRay Accumlated Density. In Primary Ray Loop. 
		
		// Early Ray Termination - Creating Clipping. 
		//if ((length(ray_P) >= rayL)) {break;}
		//if (dens > 0.99) { dens = 1.0; break;}
	}	

	//dens /= (step_count / 10); //dens = clamp(dens, 0.0, 1.0); 
	dens *= 0.75f; // Oppose to div by stepcount. 
	float dens_a = dens; //dens_f + dens_b;
	
	vec3 dens_vec = vec3(dens_a, dens_a, dens_a); 
	vec3 cv_0 = mix(samp_cf_start.xyz, samp_cb_end.xyz, 0.5);
	cv_0 = desat(cv_0, 1.0).xyz; cv_0.xy += 0.8 * length(cv_0); cv_0 *= 0.1; // cv_0.x += 0.1;
	vec3 cv_1 = mix(cv_0, dens_vec, 0.5); 
	frag_color = vec4(clamp(cv_1, 0.0, 1.0), 1.0); 
	
	//frag_color = vec4(mix(samp_cf.xyz, samp_cb.xyz, 0.5), 1.0); // Check Cube Faces Blended. 
	//vec3 dens_vec = vec3(dens, dens, dens); 
	//vec3 add_both = (ray_P_f + ray_P_b) / 2.0f;
	//frag_color = vec4(add_both, 1.0); 
	//frag_color = vec4(dens_vec, 1.0);
	//frag_color = vec4(dens, dens, dens, 1.0); 
	//frag_color = vec4(ray_dir, 1.0); // Check RayDir
	//frag_color = vec4(uv.xy, 0.0, 1.0); // Check ScreenSpace UV.
	
}
