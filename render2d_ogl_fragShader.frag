#version 400 core 

in vec4 vpos; 
in vec4 gl_FragCoord; 

out vec4 frag_color; 

// Grid Texture Samplers - 

// Density Grid - 
uniform sampler2D d_tex; // Texture Unit 0
// Velocity Texture - 
uniform sampler2D v_u_tex; // Texture Unit 1
uniform sampler2D v_v_tex; // Texture Unit 2
// Collision Grid - 
uniform sampler2D c_tex; // Texture Unit 3 
// Vorticity Confinement Gradient Grid - 
uniform sampler2D vc_u_tex; // Texture Unit 4
uniform sampler2D vc_v_tex; // Texture Unit 5
// Colour RGB - 
uniform sampler2D img_rgb_tex; // Packed RGB Texture Unit 6.

// Debug Samplers - 
uniform sampler2D ppv_u_tex; // Texture Unit 7
uniform sampler2D ppv_v_tex; // Texture Unit 8

// Constant Uniforms - 
uniform int N_Size; // Grid Size + Edge Cells. Per Dimension (N).
uniform int Mode; // 0 = Render Density, 1 = Render Velocity. 
uniform int Step; // Current Solve Step. 

// Util Functions - 
float fit (float value, float min_a, float max_a, float min_b, float max_b)
{
	return min_b + (value - min_a)*(max_b - min_b) / (max_a - min_a);
}

float deg2rad (int deg)
{
	float PI = 3.14159265;
	return deg * PI / 180.0; 
}

// Noise Functions - Patrico Gonzalez. 
float rand(float n)
{
	return fract(sin(n) * 43758.5453123);
}

/*
float noise(float p)
{
	float fl = floor(p);
	float fc = fract(p);
	return mix(rand(fl), rand(fl + 1.0), fc);
}


float noise(vec2 n)
{
	const vec2 d = vec2(0.0, 1.0);
	vec2 b = floor(n), f = smoothstep(vec2(0.0), vec2(1.0), fract(n));
	return mix(mix(rand(b), rand(b + d.yx), f.x), mix(rand(b + d.xy), rand(b + d.yy), f.x), f.y);
}
*/

void main()
{
	// Map from 0-N FragCoord_Space to 0-1 UV Space. 
	vec2 uv = (gl_FragCoord.xy - 0.0) / N_Size; // N_Size; // (Currently NOT) With -0.5f Half Pixel Offset Subtracted off. 

	// vec4.x contains sampleted channel texel data (Single Float Channel Textures).
	vec4 dens = clamp(texture(d_tex, uv), 0.0, 1.0); 
	vec4 cold = clamp(texture(c_tex, uv), 0.0, 1.0); 
	vec4 vel_x = clamp(texture(v_u_tex, uv), 0.0, 1.0); 
	vec4 vel_y = clamp(texture(v_v_tex, uv), 0.0, 1.0); 
	vec4 vc_x = clamp(texture(vc_u_tex, uv), 0.0, 1.0);
	vec4 vc_y = clamp(texture(vc_v_tex, uv), 0.0, 1.0); 
	
	// Packed RGB Channels. // Flip in GL Stb Loaded TopLeft - GL Flipped to BotLeft...
	vec4 rgb = clamp(texture(img_rgb_tex, uv), 0.0, 1.0); 
	
	// Debug Textures Sampling - 
	vec4 pp_vel_x = clamp(texture(ppv_u_tex, uv), 0.0, 1.0); 
	vec4 pp_vel_y = clamp(texture(ppv_v_tex, uv), 0.0, 1.0); 
	
	
	// Final Velocity Colour - 
	vec3 vel_f = vec3(vel_x.x, vel_y.x, 0.0);
	// Final PreProjected Velocity Colour - 
	vec3 ppvel_f = vec3(pp_vel_x.x, pp_vel_y.x, 0.0);
	// Sub Final Projected, from Pre-Projected Vel Field 
	vec3 resid_f = vel_f - ppvel_f; 
	
	
	if (Mode == 0)
	{
		// Density + Collider Viz
		frag_color = vec4(clamp(dens.x + cold.x, 0.0, 1.0), dens.x, dens.x, 1.0); 
		
		// Test Shading Colour via Velocity. 
		vec2 vel = vec2(vel_x.x, vel_y.x); 
		float vel_rm = fit(length(vel), 0.0, 1.25, 0.0, 1.0);
		//float vel_rm = length(vel); 
		//float vel_dens = clamp((dens.x * vel_rm), 0.0, 1.0); 
		
		// Mult Up Remapped Vel Coeff - 
		vel_rm = min(vel_rm * 2.5, 1.0); 
		vec3 vel_col = mix(vec3(0.0, 0.0, 1.0), vec3(1.0, 1.0, 1.0), vel_rm); 
		vel_col *= min(dens.x, 1.0); // Only Where Density Is. 
		
		// Oscillate Colour - (Disco Mode)
	//	vel_col *= vec3(abs(sin(deg2rad(Step))), abs(cos(deg2rad(Step + 1234))), abs(sin(deg2rad(Step + 4312)))); 
		
		// Add Collider Viz Set FragCol. 
		vel_col.x += cold.x; // Add Collider Colour to Red.
		frag_color = vec4(vel_col, 1.0); 
		
		/* Velocity Colour Rendering - 
		frag_color = vec4(vel_dens + cold.x, vel_dens, vel_dens, 1.0); 
		//RGB Colour Rendering - 
		frag_color = vec4(clamp(rgb.x + cold.x, 0.0, 1.0), rgb.y, rgb.z, 1.0); 
		*/
		
		/////////
		
		// Pseduo 2D RayMarching. (UV ScreenSpace) \\  
		float sinx = fit(sin(deg2rad(Step) * 5.0), -1.0, 1.0, 0.0, 1.0); 
		float cosx = fit(cos(deg2rad(Step) * 3.0), -1.0, 1.0, 0.0, 1.0); 
		
		// Animated "Light" Pos. 
		vec2 light_pos = vec2(0.99, 0.99); // UV Space 
		//light_pos = vec2(sinx, cosx); 
		//light_pos = vec2(sinx, 0.99); 
		
		vec4 acc; // Accumulated Color // + Alpha
		vec2 ray_P = uv; // Inital RayPos at Fragment UV. 
		//vec2 dir = normalize(vec2(light_pos - uv)); // UV Space. 
		vec2 dir = normalize(vec2(uv - light_pos)); // Direction is Inversed in OGL. == LightPos - UVPos. 
		
		int max_step = 10; 
		float step_size = 0.01; 
		//float step_size = 0.15 / max_step; // Step_Size Based on StepCount. 
		
		for (int i = 0; i < max_step; i++)
		{
			// Iteration to 0-1 Coeff
			float mult = float(i) / float(max_step); 
			// Iteration 1-0 Coeff
			float inv_mult = 1 - mult; 
			inv_mult = fit(inv_mult, 0.0, 1.0, 0.1, 1.5); 
			
			/*
			// Stop Overshooting SamplePos -> LightPos (WIP)
			if (dot(normalize(ray_P), normalize(light_pos)) >= 0.99) 
			{
				//acc = vec4(1.0, 1.0, 1.0, 1.0); 
				break;
			}
			*/
			
			// Accumlation Of Colour w/ Reverse Absorbtion / Increase Brightness as Approach Light
			acc += clamp(texture(d_tex, clamp(ray_P, 0.0, 1.0)), 0.0, 1.0); 
			//mult = smoothstep(0.0, 1.0, mult);// mult = pow(mult, 0.75); 
			acc *= mult; 
			
			// Increment Ray Towards Light 
			//float noise_v = vec2(noise(ray_P.x, ray_P.y), noise(ray_P.x + 0.5, ray_P.y + 0.1)); 
			//noise_v *= 0.5; 
			
			ray_P += dir * step_size; 
			
			// Decrease StepSize Each Iter - 
			//step_size *= fit(inv_mult, 0.0, 1.0, 0.2, 1.0);  
		}
		
		// Oppose to Full Averaging, Lerp between Non Averaged, and Averaged Accumlation Colour. 
		float acc_A = acc.x;// max_step; 
		// Avoid Generic Averaged Blurrly Look !
		float acc_mix = mix(acc_A, acc.x, 0.5); // Mix Between Accumlated Average with Non Averaged. 
		
		//frag_color = vec4(acc.x + cold.x, acc.x, acc.x, 1.0); 
		
		frag_color = vec4(acc_A + cold.x, acc_A, acc_A, 1.0); 
		
		//frag_color = vec4(acc_mix + cold.x, acc_mix, acc_mix, 1.0); 
		//frag_color = vec4(uv.x, uv.y, 0.0, 1.0); 
		//frag_color = vec4(dir, 0.0, 1.0); 
		
	}
	else if (Mode == 1)
	{
		// Velocity + Collider Viz
		frag_color = vec4(clamp(vel_x.x + cold.x, 0.0, 1.0), vel_y.x, 0.0, 1.0); 
		
		// Vorticity + Collider Viz
	//	frag_color = vec4(clamp(vc_x.x + cold.x, 0.0, 1.0), vc_y.x, 0.0, 1.0); 
	
		// Residual Velocity Divergence Test Viz - 
	//	frag_color = vec4(resid_f, 1.0); 
		
		// Pre Projected Velocity Viz - 
	//	frag_color = vec4(pp_vel_x.x,pp_vel_y.x,0.0,1.0);
	}
	
}
