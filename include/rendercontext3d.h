#ifndef RENDERCONTEXT_H
#define RENDERCONTEXT_H
// Render Context Creation and Window Initalization. (Eventually) For Muliplte APIs. 

#define GLEW_STATIC // Static Linking GLEW and GLFW For now. 
#include <GLEW\glew.h>
#include <GLFW\glfw3.h>

// Abstract Base Class Interface for Each Specfic API To Implement Render Context Setup. 
class render_context
{
public:
	render_context(int sizex, int sizey); // Constructor For Inherited Initalization Calls Only. 
	virtual ~render_context() = default; 

protected:
	virtual void window_context() = 0; // Window and Context Setup. 
	virtual void extensions_setup() = 0; // Extensions Loading Setup. 

	int wind_x, wind_y;
};


class render_context_OGL : public render_context
{
public:
	render_context_OGL(int xx, int yy, short v_maj, short v_min);
	~render_context_OGL(); 

	// Window Ptr Getter. 
	GLFWwindow* get_window(); 

protected:
	// Override of Init/Setup MFs. 
	void window_context() override;
	void extensions_setup() override;

private:
	GLFWwindow *window = nullptr;
	const GLubyte *render_device = nullptr; 
	const GLubyte *version = nullptr; 
	short gl_ver_major, gl_ver_minor; 
};

//class render_context_DX11 : public render_context { };

#endif 