#include "ModelRenderer.h"


ModelRenderer::ModelRenderer(int window_width, int window_height, int image_width, int image_height):
	_width(window_width), _height(window_height), _image_width(image_width),  _image_height(image_height)
{
	_obj_model = NULL;
	_fboHidden = -1;
	_color_texture_idx = -1;
	_depth_texture_idx = -1;

	_fboHiddenNormals = -1;
	_normal_texture_idx = -1;
	_normal_depth_texture_idx = -1;

	_output_file_id = 0;
	_save = false;
	_output_file_path = "out";
	_output_file_name = "model";
	_writer_enabled = true;
	
	_with_roi = true; 
	_with_mask = true;

	_verbose = false;

	_projectionMatrix = glm::perspective(1.2f, (float)800 / (float)600, 0.1f, 100.f);
	_projectionMatrix = glm::perspective( glm::radians(40.0f), (float)480 / (float)480, 0.1f, 100.f);
	_viewMatrix = glm::lookAt(glm::vec3(0.0f, 0.0, 1.3f), glm::vec3(0.0f, 0.0f, 00.f), glm::vec3(0.0f, 1.0f, 0.0f));
	_modelMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f));  //glm::rotate(3.1415f, glm::vec3(0.0f, 1.0f, 0.0f)) *  glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f));


	_light0.index = 0;
	_light0.pos = glm::vec3(0.0f, 6.0f, -16.0f);
	_light0.dir = glm::vec3(0.0f, 0.0f, 0.0f);
	_light0.k1 = 0.1;
	_light0.intensity = 1.0;
	_light1.pos = glm::vec3(0.0f, 3.0f, 3.0f);
	_light1.dir = glm::vec3(0.0f, 0.0f, 0.0f);
	_light1.k1 = 0.1;
	_light1.index = 1;
	_light1.intensity = 1.0;
	_light2.index = 2;
	_light2.pos = glm::vec3(0.0f, 6.0f, 10.0f);
	_light2.dir = glm::vec3(0.0f, 0.0f, 0.0f);
	_light2.color = glm::vec3(0.0f, 0.0f, 1.0f);
	_light2.intensity = 1.0;


	_mat0.diffuse_mat = glm::vec3(1.0, 0.0, 0.0);
	_mat0.diffuse_int = 0.8;
	_mat0.ambient_mat = glm::vec3(1.0, 0.0, 0.0);
	_mat0.ambient_int = 0.2;
	_mat0.specular_mat = glm::vec3(1.0, 1.0, 1.0);
	_mat0.specular_int = 0.2;
	_mat0.specular_s = 4.0;
	
#ifndef _DEVELOP
	_light0.with_error_check = false;
	_light1.with_error_check = false;
	_mat0.with_error_check = false;
#endif

	_data_rgb = (unsigned char*)malloc(_image_width * _image_height * 3 * sizeof(int));
	_data_depth = (unsigned char*)malloc(_image_width * _image_height * 1 * sizeof(float));
	_data_normals = (unsigned char*)malloc(_image_width * _image_height * 3 * sizeof(float));

	_writer = new ImageWriter();
}


ModelRenderer::~ModelRenderer()
{
	free(_data_rgb);
	free(_data_depth);
	free(_data_normals);

	delete _writer;
}


/*
Set the model to render
@param path_and_file - string containg the relative or absolute
						path to the model.
@return - true, if model was successfully loaded.
*/
bool ModelRenderer::setModel(string path_and_file)
{
	if (path_and_file.empty()) return false;
	
//#define _DEVELOP
//#ifdef _DEVELOP
//	// load shader
//	int program = cs557::LoadAndCreateShaderProgram("./shaders/image_renderer.vs", "./shaders/image_renderer.fs");
//#else
//	int program = cs557::CreateShaderProgram(glslshader::image_renderer_vs, glslshader::image_renderer_fs);
//#endif
	
	// create model
	_obj_model = new cs557::OBJModel();
	_obj_model->create(path_and_file, model_pgm);
	//_obj_model->setTextureParam(cs557::TextureMode::REPLACE);

	_light0.apply(model_pgm);
	_light1.apply(model_pgm);
	//_light2.apply(model_pgm);
	//_mat0.apply(program);
	//-----------------------
	// create a second object to render only normal vectors
	// load shader
//#ifdef _DEVELOP
//	int program_normals = cs557::LoadAndCreateShaderProgram("./shaders/normal_renderer.vs", "./shaders/normal_renderer.fs");
//#else
//	int program_normals = cs557::CreateShaderProgram(glslshader::normal_renderer_vs, glslshader::normal_renderer_fs);
//#endif
	int program_normals = cs557::LoadAndCreateShaderProgram("./shaders/normal_renderer.vs", "./shaders/normal_renderer.fs");
	_obj_model_normals = new cs557::OBJModel();
	_obj_model_normals->create(path_and_file, program_normals);

	_light0.apply(program_normals);
	_light1.apply(program_normals);
	//_light2.apply(program_normals);
	

	//CreatePrerendererScene();
	//CreateHelperContent();

}

bool ModelRenderer::setMaterial(string name)
{
	if (name.empty()) return false;
	////////////////////////////////////////////////////////////////////////////////////// 
	//Setup for Environment Mapping starts here.
	
	//cs557::LoadBMPFromFile("../PBR/reflactance_map_lake.bmp", &width, &height, &channels, &g_data);
	cs557::LoadBMPFromFile("../PBR/reflactance_map_city.bmp", &width, &height, &channels, &g_data);
	// Activate the texture unit and bind the texture. 
	// Note that this function binds the texture to texture unit GL_TEXTURE0 
	glActiveTexture(GL_TEXTURE4);
	glGenTextures(1, &texture_id0);
	// Set a texture as active texture.
	glBindTexture(GL_TEXTURE_2D, texture_id0);
	// Change the parameters of your texture units.
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	//**********************************************************************************************
	// Create a texture and load it to your graphics hardware. This texture is automatically associated
	// with texture 0 and the textuer variable "texture" / the active texture.
	if (channels == 3)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_BGR, GL_UNSIGNED_BYTE, g_data);
	else if (channels == 4)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_BGRA, GL_UNSIGNED_BYTE, g_data);

	glGenerateMipmap(GL_TEXTURE_2D);

	// Fetch the texture location and set the parameter to 0.
	// Note that 0 is the number of the texture unit GL_TEXTURE0.
	int texture_location = glGetUniformLocation(model_pgm, "gradient");
	glUniform1i(texture_location, 4);
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup for Environment mapping ends here
	
	if (name == "Fe") {
		cs557::BRDFLoader::ReadBMP(brdf_tex0, "../PBR/Iron-Scuffed/Iron-Scuffed-albedo.bmp",
			"../PBR/Iron-Scuffed/Iron-Scuffed-ao.bmp",
			"../PBR/Iron-Scuffed/Iron-Scuffed-rough.bmp",
			"../PBR/Iron-Scuffed/Iron-Scuffed-metal.bmp");
		brdf_tex0.lightColor = glm::vec3(73.47, 71.31, 74.79);
		brdf_tex0.F0 = glm::vec3(0.86, 0.87, 0.88);
		brdf_tex0.k1 = 0.0;
		brdf_tex0.k2 = 0.0;
		brdf_tex0.apply(model_pgm);
	}

	else if (name == "WOOD") {
		cs557::BRDFLoader::ReadBMP(brdf_tex0, "../PBR/bamboo-wood-semigloss/bamboo-wood-semigloss-albedo.bmp",
			"../PBR/bamboo-wood-semigloss/bamboo-wood-semigloss-metal.bmp",
			"../PBR/bamboo-wood-semigloss/bamboo-wood-semigloss-rough.bmp",
			"../PBR/bamboo-wood-semigloss/bamboo-wood-semigloss-ao.bmp");
		brdf_tex0.lightColor = glm::vec3(53.47, 51.31, 50.79);
		brdf_tex0.F0 = glm::vec3(0.25, 0.25, 0.25);
		brdf_tex0.k1 = 0.05;
		brdf_tex0.k2 = 0.02;
		brdf_tex0.apply(model_pgm);
	}

	else if (name == "Cu") {
		cs557::BRDFLoader::ReadBMP(brdf_tex0, "../PBR/Copper-scuffed/Copper-scuffed-albedo.bmp",
		"../PBR/Copper-scuffed/Copper-scuffed-metal.bmp",
		"../PBR/Copper-scuffed/Copper-scuffed-rough.bmp",
		"../PBR/Copper-scuffed/streaked-metal1-ao.bmp");
		brdf_tex0.lightColor = glm::vec3(73.47, 71.31, 74.79);
		brdf_tex0.F0 = glm::vec3(0.86, 0.87, 0.88);
		brdf_tex0.k1 = 0.0;
		brdf_tex0.k2 = 0.0;
		brdf_tex0.apply(model_pgm);
	}
	else if (name == "Painted") {
		cs557::BRDFLoader::ReadBMP(brdf_tex0, "../PBR/chipped-paint-metal/chipped-paint-metal-albedo.bmp",
			"../PBR/chipped-paint-metal/chipped-paint-metal-metal.bmp",
			"../PBR/chipped-paint-metal/chipped-paint-metal-rough.bmp",
			"../PBR/chipped-paint-metal/chipped-paint-metal-ao.bmp");
		brdf_tex0.lightColor = glm::vec3(73.47, 71.31, 74.79);
		brdf_tex0.F0 = glm::vec3(0.86, 0.87, 0.88);
		brdf_tex0.k1 = 0.0;
		brdf_tex0.k2 = 0.0;
		brdf_tex0.apply(model_pgm);
	}
	/*cs557::BRDFLoader::ReadBMP(brdf_tex0, "../PBR/layerd-rock1/layered-rock1-albedo.bmp",
	"../PBR/layerd-rock1/layered-rock1-Metalness.bmp",
	"../PBR/layerd-rock1/layered-rock1-rough.bmp",
	"../PBR/layerd-rock1/layered-rock1-ao.bmp");*/

	/*cs557::BRDFLoader::ReadBMP(brdf_tex0, "../PBR/worn_braided_carpet/worn-braided-carpet-albedo.bmp",
	"../PBR/worn_braided_carpet/worn-braided-carpet-Metallic.bmp",
	"../PBR/worn_braided_carpet/worn-braided-carpet-Roughness.bmp",
	"../PBR/worn_braided_carpet/worn-braided-carpet-ao.bmp");*/
	

	
}

/*
Draw the scene into an fbo to get textures
*/
bool ModelRenderer::drawFBO(void)
{
	if (_obj_model == NULL) return false;

	glBindFramebuffer(GL_FRAMEBUFFER, _fboHidden);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	// Set up our green background color
	static GLfloat clear_color[] = { 0.0f, 0.0f, 0.0f, 0.0f };//{0.6f, 0.7f, 1.0f, 1.0f};
	static GLfloat clear_depth[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	// Clear the entire buffer with our green color.
	glClearBufferfv(GL_COLOR, 0, clear_color);
	glClearBufferfv(GL_DEPTH, 0, clear_depth);

	// set the viewport. It must match the texture size.
	glViewport(0, 0,  _image_width, _image_height);
	brdf_tex0.activateMaps(_obj_model->getProgram());
	_obj_model->draw(_projectionMatrix, _viewMatrix, _modelMatrix);

	//-------------------------------------------------------------------------------------
	// get the data back 

	// read back the pixels
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, _image_width, _image_height, GL_BGR, GL_UNSIGNED_BYTE, _data_rgb);

	glReadBuffer(GL_DEPTH_ATTACHMENT);
	glReadPixels(0, 0, _image_width, _image_height, GL_DEPTH_COMPONENT, GL_FLOAT, _data_depth);

	// switch back to the regular output buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// set the viewport to window size
	glViewport(0, 0, _width, _height);

	// rgb image
	cv::Mat image(_image_height, _image_width, CV_8UC3, _data_rgb);
	cv::Mat dst, output_rgb;
	cv::flip(image, dst, 0);
	

	// depth image
	cv::Mat imaged(_image_height, _image_width, CV_32FC1, _data_depth);
	cv::Mat dst_depth, output_depth, normalized;
	cv::flip(imaged, dst_depth, 0);



	//-------------------------------------------------------------------------------------
	// Draw normals
	glBindFramebuffer(GL_FRAMEBUFFER, _fboHiddenNormals);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	static GLfloat clear_normals[] = { 0.0f, 0.0f, 0.0f, 0.0f };

	// Clear the entire buffer with our green color.
	glClearBufferfv(GL_COLOR, 0, clear_normals);
	glClearBufferfv(GL_DEPTH, 0, clear_depth);

	// set the viewport. It must match the texture size.
	glViewport(0, 0, _image_width, _image_height);
	_obj_model_normals->draw(_projectionMatrix, _viewMatrix, _modelMatrix);


	//-------------------------------------------------------------------------------------
	// get the data back 

	// read back the pixels
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, _image_width, _image_height, GL_BGR, GL_FLOAT, _data_normals);

	cv::Mat image_normals(_image_height, _image_width, CV_32FC3, _data_normals);
	cv::Mat dst_norm, output_norm;
	cv::flip(image_normals, dst_norm, 0);


	//-------------------------------------------------------------------------------------
	// region of interest extraction
	cv::Rect2f roi;
	if (_with_roi) {
		RoIDetect::Extract(dst, roi);
	}

	//-------------------------------------------------------------------------------------
	// Extract an image mask
	cv::Mat mask;
	if (_with_mask) {
		ImageMask::Extract(dst, mask);
	}


	// switch back to the regular output buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// set the viewport to window size
	glViewport(0, 0, _width, _height);


	if (_verbose) {

		// to normalize the depth image
		//cv::normalize(dst_norm, output_norm, 0, 255, cv::NORM_MINMAX, CV_8UC3);
		cv::resize(dst, output_rgb, cv::Size(512, 512));
		cv::resize(dst_depth, output_depth, cv::Size(512, 512));
		cv::resize(dst_norm, output_norm, cv::Size(512, 512));

		cv::imshow("RGB image (3 x uchar)", output_rgb);
		cv::imshow("Depth image (float)", output_depth);
		cv::imshow("Normal image (float)", output_norm);
		cv::waitKey(1);

		if (_verbose && _with_roi) 
			RoIDetect::RenderRoI(dst, roi);
	
	}

	

	if (_save && _writer_enabled && _writer){

		ImageWriter::IWData odata;
		odata.index = _output_file_id;
		odata.rgb = &dst;
		odata.normals = &dst_norm;
		odata.depth = &dst_depth;
		if (_with_mask) {
			odata.mask = &mask;
		}
		// this view matrix describes the object's pose in camera coordinates since the object is at 0,0,0
		odata.pose = _viewMatrix;
		odata.roi = roi;
		_writer->write(odata);

		//_writer->write(_output_file_id, dst, dst_norm, dst_depth, glm::inverse(_viewMatrix));
		_output_file_id++;
	}

	return true;
}

/*
Draw the current object in a window
*/
bool ModelRenderer::draw(void)
{
	if (_obj_model == NULL) return false;


	drawFBO();

	//brdf_tex0.activateMaps(_obj_model->getProgram());
	_obj_model->draw(_projectionMatrix, _viewMatrix, _modelMatrix);
    _coordinateSystem.draw(_projectionMatrix, _viewMatrix, _modelMatrixCoordSystem);
	_display.draw(_projectionMatrix, _viewMatrix, _display_m);
	brdf_tex0.activateMaps(_obj_model->getProgram());
	glActiveTexture(GL_TEXTURE4);                   // Environment texture 
	glBindTexture(GL_TEXTURE_2D, texture_id0);
	_obj_model->draw(_projectionMatrix, _viewMatrix, glm::translate(glm::mat4(1.0f), glm::vec3(4.5f, 2.0f, 0.0f)));
	//glUseProgram(_display.getProgram());
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, _depth_texture_idx);
}


/*
	Draw the current object into a window and into the fbo.
	Save the fbo to file
	*/
bool  ModelRenderer::draw_and_save(void)
{
	if (_obj_model == NULL) return false;
	_save = true;

	draw();
}


/*
	Create a scene for the prerenderer
	*/
void  ModelRenderer::CreatePrerendererScene(void)
{

	// This function is part of RenderToTexture.h
	CreateRenderToTexture(_image_width, _image_height, _fboHidden, _color_texture_idx, _depth_texture_idx);

	// Reset to the regular buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);


	CreateRenderToTexture32Bit(_image_width, _image_height, _fboHiddenNormals, _normal_texture_idx, _normal_depth_texture_idx);
	//cout << _fboHiddenNormals << " : " << _normal_texture_idx << " : " << _normal_depth_texture_idx << endl;
	
	// Reset to the regular buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glUseProgram(0);
}


/*
Creates some helper displays
*/
void ModelRenderer::CreateHelperContent(void)
{
	// Load the shader program
//#ifdef _DEVELOP
//	int shader = cs557::LoadAndCreateShaderProgram("./shaders/display.vs", "./shaders/display.fs");
//#else
//	int shader = cs557::CreateShaderProgram(glslshader::display_renderer_vs, glslshader::display_renderer_fs);
//#endif
	int shader = cs557::LoadAndCreateShaderProgram("./shaders/display.vs", "./shaders/display.fs");
	// create a plane
	_display.create(0.56, 0.45, shader);
	_display_m = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f));
	_display_p = glm::perspective(1.57f, (float)400 / (float)400, 0.01f, 10.f);
	_display_v = glm::lookAt(glm::vec3(0.0f, 0.0f, 1.f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));

	glUseProgram(shader); 
	/*glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, _color_texture_idx);
	int texture_location = glGetUniformLocation(shader, "tex");
	glUniform1i(texture_location, 0);
	glUniform1f(glGetUniformLocation(shader, "display_scale"), 1.0f);*/

	glUseProgram(0);


	// create a coordinate system
	_modelMatrixCoordSystem = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f)); 
    _coordinateSystem.create(25);
}


/*
Set the intrinsic parameters
@param  fx, fy, Field of view in pixels into the x and y direction.
@parm	px, py, the camera principle point
*/
void ModelRenderer::setIntrinsic(float fx, float fy, float px, float py)
{
	_projectionMatrix = glm::perspective(fx/fy, (float)fx / (float)fy, 0.1f, 10.f);
	
}

/*
Set the view matrix.
@param viewmatrix - a 4x4 view matrix.
*/
void ModelRenderer::setCameraMatrix(glm::mat4 viewmatrix)
{
	_viewMatrix = viewmatrix;

	// update the light position so that it stays a camera front light. 
	glm::mat4 inv  = glm::inverse(viewmatrix);
	_light0.pos = glm::vec3(inv[3][0], inv[3][1], inv[3][2]);
	_light0.apply(_obj_model->getProgram());
	_light0.apply(_obj_model_normals->getProgram());
}

/*
Set an output path for the images
@param path - relative or absolute output path
@param name - name template for the output images.
*/
void ModelRenderer::setOutputPath(string path, string name)
{
	_output_file_path = path;
	_output_file_name = name;

	if(_writer)
		_writer->setPathAndImageName(path, name);
}


/*
Disable and enable the file writer
*/
bool  ModelRenderer::enable_writer(bool enable)
{
	_writer_enabled = enable;

	return true;
}