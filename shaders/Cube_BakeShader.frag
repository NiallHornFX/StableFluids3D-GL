#version 440 core

out vec4 frag_colour; 

in vec3 vtc; // Get Vertex Uniform TexCoord.


void main()
{
	// Colour Fragments Local Cube/Texture Coords - 
	
	frag_colour = vec4(vtc, 1.0); 
}