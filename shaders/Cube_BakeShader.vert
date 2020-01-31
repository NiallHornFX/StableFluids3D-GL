#version 400 core 

layout (location = 0) in vec3 vPos;

out vec4 vpos; 

void main()
{
	gl_Position = vec4(vPos.x, vPos.y, vPos.z, 1.0); 
	vpos = vec4(vPos.x , vPos.y , vPos.z , 1.0); 
}

