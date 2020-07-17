#include "GLRenderer.h"




GLRenderer::GLRenderer()
{
	_init_ready = false;
	_projMatrix = glm::perspective(1.2f, (float)800 / (float)600, 0.1f, 100.f);
	_viewMatrix = glm::lookAt(glm::vec3(0.0f, 0.0, 4.5f), glm::vec3(0.0f, 0.0f, 00.f), glm::vec3(0.0f, 1.0f, 0.0f));
	//_modelMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f)); 

	// Set up our green background color
	for (int i = 0; i < 4; i++ ) _clear_color[i] = 1.0f;
	for (int i = 0; i < 4; i++ ) _clear_depth[i] = 1.0f;

}


GLRenderer::~GLRenderer()
{

}

/*
Create the renderer instance
@param window_width - width of the window in pixel
@param window_height - height of the window in pixel
@param name - label for the window as string
*/
bool GLRenderer::create(int window_width, int window_height, string name)
{
	_window_width = window_width;
	_window_height = window_height;

	_projMatrix = glm::perspective(1.2f, (float)window_width / (float)window_height, 0.1f, 100.f);

	
	// Init the GLFW Window and glew
    _window = cs557::initWindow(window_width, window_height, name);
    cs557::initGlew();

	// coordinate system
	_cs.create(3.0);


	return true;
}


/*
Start the renderer
*/
bool GLRenderer::start(void)
{
	_running = true;
	draw_loop();
	return true;
}


/*
The main draw loop
*/
void GLRenderer::draw_loop(void)
{

	 // Enable depth test
    glEnable(GL_DEPTH_TEST); 
    glEnable(GL_BLEND); 
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

	 // Init the view matrix. 
    cs557::InitControlsViewMatrix(_viewMatrix);

	glm::mat4 cs_modelMatrix =  glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f)); 
	clock_t begin = clock();

	_init_ready = true;

    while (!glfwWindowShouldClose(_window))
    {
		  // Clear the entire buffer with our green color (sets the background to be green).
        glClearBufferfv(GL_COLOR, 0, _clear_color);
        glClearBufferfv(GL_DEPTH, 0, _clear_depth);

		glm::mat4 rotated_view =    cs557::GetCamera().getViewMatrix() * glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f)) ;

		// draw the coordinate frame
		_cs.draw(_projMatrix, rotated_view, cs_modelMatrix);

		/*
		Call all render callbacks
		*/
		for (auto f : _render_calls) {
			f(_projMatrix, rotated_view);
		}

		glfwSwapBuffers(_window);
        glfwPollEvents();

		//_sleep(0.03);
	}

}


bool GLRenderer::addRenderFcn(std::function<void(glm::mat4 pm, glm::mat4 vm)> function)
{
	_render_calls.push_back(function);

	return true;
}


/*
Set a view matrix 
@param vm - 4x4 view matrix
*/
bool  GLRenderer::setViewMatrix(glm::mat4 vm)
{
	_viewMatrix = vm;
	 // Init the view matrix. 
    cs557::InitControlsViewMatrix(_viewMatrix);

	return true;
}


/*
Add a keyboard callback of type
void name (int key, int action)
*/
bool  GLRenderer::addKeyboardCallback(std::function<void(int, int)> function)
{
	cs557::AddKeyboardCallbackPtr(function);

	return true;
}