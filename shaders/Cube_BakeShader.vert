#version 420 core 

layout (location = 0) in vec3 vPos; // Pos
layout (location = 1) in vec3 vTC; // TexCoord Attrib

uniform mat4 model; 
uniform mat4 view; 
uniform mat4 persp; 

out vec3 vpos; 
out vec3 vtc; 

void main()
{
	gl_Position = persp * view * model * vec4(vPos, 1.0);
	// persp * view *

	// Marshall TexCoord to Fragment
	vtc = vTC; 
	vpos = vPos; 
}

