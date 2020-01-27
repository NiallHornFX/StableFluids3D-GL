#version 400 core 

in vec4 vpos; 
in vec4 gl_FragCoord; 

out vec4 frag_color; 

// Grid Texture Samplers - 

// Density Grid - 
uniform sampler3D d_tex; // Texture Unit 0
uniform sampler3D v_tex; // Texture Unit 1 (Packed x,y,z vel comps). 

// Constant Uniforms - 
// Util Uniforms - 
uniform int W_Size; // Window Size to Get GLFragCoords -> UV Space for Inital XY RayMarch Postions.
uniform int Mode; // 0 = Render Density, 1 = Render Velocity. 
uniform int Step; // Current Solve Step. 

// Volume Uniforms - 
uniform vec3 offset; // Sample Offset of Grid/3DTex. Pre Matrix Transformation Implementation. 

// Light Uniforms - 
struct pt_light 
{
	vec3 pos;
	float radius, strength;
	
} light_00;

// Render Uniforms - 


/* ----------------------------------------------------------------- */

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

/* ----------------------------------------------------------------- */

void main()
{
	// Map from 0-N FragCoord_Space to 0-1 UV Space. 
	vec2 uv = (gl_FragCoord.xy - 0.0) / W_Size; // N_Size; // (Currently NOT) With -0.5f Half Pixel Offset Subtracted off. 

	if (Mode == 0)
	{
		// Single Light (fn) - 
		light_00.pos = vec3(0.5, 1.25, 0.0); 
		light_00.radius = 1.0; light_00.strength = 1.0; 
		
		// Basic RayMarching Inital - 
		int max_steps = 50;
		float step_size = 1.0 / max_steps;
		//float step_size = 0.05;
		vec3 dir = vec3(0.0, 0.0, -1.0); 
		vec3 ray_P = vec3(uv, 0.0); 
		
		vec4 acc = vec4(0.0, 0.0, 0.0, 0.0); 
		//vec4 accv = vec4(0.0, 0.0, 0.0, 0.0); 
		int total_i = 0; 
		for (int i = 0; i < max_steps; i++)
		{
			// Primary Camera Ray (a)
			acc += texture(d_tex, (ray_P)); // ray_P + offset
			//accv += texture(v_tex, (ray_P)); 
			ray_P += dir * step_size; 
			total_i++;
			if (acc.x >= 1.0) {break;}
			
			// Secondary Shadow Ray (c) - 
		}
		
		frag_color = vec4(acc.x, acc.x, acc.x, 1.0); 
		
		//frag_color = vec4(uv, 0.0, 1.0);
		//frag_color = vec4(1.0, 0.0, 0.0, 1.0);
		
		float viz = float(total_i) / float(max_steps); // float(max_steps);
		// if (total_i <= 5)
		// {
			// frag_color = vec4(1.0, 0.0, 0.0, 1.0); 
		// }
		// Shade By Ray Step (Depth till >= 1.0 exit) Count. 
		frag_color = vec4((1.0-viz) * 5.0, (1.0-viz) * 5.0, (1.0-viz) * 5.0, 1.0); 
		
		
		// Vel Tex Tiling, Not due to UV/Fragment Sample pos, probs due to XYZ->RGB Packing Issue. 
		//vec4 v_test = texture(v_tex, vec3(uv, -0.1)); 
		//frag_color = vec4(vec3(accv.xyz), 1.0); 
		
	}
	else if (Mode == 1)
	{
		// Velocity + Collider Viz
		//frag_color = vec4(clamp(vel_x.x + cold.x, 0.0, 1.0), vel_y.x, 0.0, 1.0); 

		//frag_color = vec4(accv.x * 10.0, accv.y * 10.0 , accv.z * 10.0, 1.0); 
	}
	
}
