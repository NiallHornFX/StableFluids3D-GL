#version 400 core 

in vec4 vpos; 
in vec4 gl_FragCoord; 

out vec4 frag_color; 

// Unit Cube Textures -
uniform sampler2D cf_tex; // 0
uniform sampler2D cb_tex; // 1


/* ----------------------------------------------------------------- */

void main()
{
	// Map from 0-N FragCoord_Space to 0-1 UV Space. 
	vec2 uv = (gl_FragCoord.xy - 0.0) / 512; 
	vec4 samp = texture(cf_tex, uv); 
	frag_color = vec4(samp.xyz, 1.0); 
	//frag_color = vec4(uv.xy, 0.0, 1.0); 
	
}
