#version 410 core     
#define MAX_LIGHTS 4

uniform mat4 projectionMatrix;                                    
uniform mat4 viewMatrix;                                           
uniform mat4 modelMatrix;  
                                                                   
in vec3 pass_Normal;     
in vec3 pass_Position;   
in vec4 pass_Color;   
in vec4 pass_Coordinates;
in vec2 pass_Texture; 
// The material parameters
uniform struct LightSource {
    vec3 position;
	vec3 direction;
	vec3 color;
	float intensity;
	float cutoff_in; // inner
	float cutoff_out; // outer
	float k1;    //attenuation
	float k2;    //attenuation
	bool used;
	int  type;  //0:point, 1:spot, 2:directional
} light[MAX_LIGHTS];

uniform struct BRDFTexMaterial {
    sampler2D albedoMap;
    sampler2D normalMap;
    sampler2D metallicMap;
    sampler2D roughnessMap;
    sampler2D aoMap;

    vec3 lightColor; // color of the light. 
    vec3 F0; 

    float k1; // linear attenuation
    float k2; // quadratic attenuation
} brdf_tex[0];

/*// The material parameters
uniform struct Material {
    vec3  diffColor;
 	float diffInt;
 	vec3  ambColor;
 	float ambInt;
 	vec3  specColor;
 	float specInt;
 	float shininess;
} mat[1];*/


out vec4 frag_out;
float calculateAttenuation(vec3 light_position, vec3 fragment_position, float k1, float k2)
{
    float distance    = length(light_position - fragment_position);
    float attenuation = 1.0 / (1.0 + k1 * distance +  k2 * (distance * distance));  
    return attenuation;
}



vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}


float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;
	
    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
	
    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float num   = NdotV;
    float denom = NdotV * (1.0 - k) + k;
	
    return num / denom;
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);
	
    return ggx1 * ggx2;
}
       

/*
Per-fragment light. 
Note that all vectors need to be in camera/eye-space. 
@param L - vector between light and fragment position.
@param E - eye point
@param s - the light source
@param m - the material 
*/
/*vec4 useLight(vec3 L, vec3 E, vec3 N, LightSource s, Material m)
{
	vec4 color = vec4(0.0,0.0,0.0,0.0);

	// diffuse light
	vec3 Idiff =  m.diffInt * m.diffColor *  max(dot(L, N), 0.0); 
	Idiff = clamp(Idiff, 0.0, 1.0); 

	// ambient light
	vec4 Iamb = vec4( m.ambColor, 1.0) * m.ambInt;

	// specular light
	vec3 R = reflect(L, N);
	vec3 Ispec =  m.specInt * m.specColor *  pow(max(dot(R, E), 0.0), m.shininess); 

	// calculate color                                     
	color = max(vec4( ( Idiff + Ispec) * s.color, 1.0), Iamb)  * s.intensity;  

	// attenuation 
	float dist = length(L);
	float Iatt = 1 / (1 + s.k1 * dist + s.k2*s.k2 * dist);

	if(s.type == 0)//pointlight
	{
		color = color * Iatt;
	}
	else if(s.type == 1)// spotlight
	{
		vec4 K = viewMatrix * vec4(normalize(s.direction - s.position), 0.0);
		float a = dot(-L,K.xyz); // angle between light dir and fragment position. 

		float c = smoothstep(1.0-s.cutoff_out, 1.0-s.cutoff_in, a); 
		Iatt = Iatt * c; // multiple with attenuation to maintain the distance effect. 

		color = color * Iatt; // attenutation
	}


	return color;
}*/

                       
void main(void)                                                  
{    
	// read the color values from the map
    vec3 albedo_color = pow(texture(  brdf_tex[0].albedoMap, pass_Texture).rgb, vec3(2.2));
    float metallic_value = texture(  brdf_tex[0].metallicMap, pass_Texture).r;
    float roughness_value = texture(  brdf_tex[0].roughnessMap, pass_Texture).r;
    float ao_value = texture(  brdf_tex[0].aoMap, pass_Texture).r;

    //vec3 roughness_color = texture(  brdf_tex[0].roughnessMap, pass_Texture).rgb;

	// eye position 
	vec3 E = normalize( vec3(viewMatrix[3][0], viewMatrix[3][1], viewMatrix[3][2]) );
    
	vec3 F0 = brdf_tex[0].F0;
    F0 = mix(F0, albedo_color, metallic_value);

    // camera view to fragment position
    vec3 V = normalize( -pass_Position);

     // normal vector
    vec3 N = normalize(pass_Normal);
	vec3 Lo = vec3(0.0);
    for (int i=0; i<MAX_LIGHTS; i++){

        if(light[i].used == false) continue;

        // Light in view space
        vec3 L_viewspace = (viewMatrix * vec4(light[i].position, 1.0)).rgb;

        // light to fragment position
        vec3 L = normalize(L_viewspace - pass_Position);

          // halfway -vector
        vec3 H = normalize(V + L);


        // calculate the radiance comming from the light source. 
        vec3  lightColor  = brdf_tex[0].lightColor;
        float cosTheta    = max(dot(N, L), 0.0);
        float attenuation = calculateAttenuation(pass_Position, L_viewspace, brdf_tex[0].k1, brdf_tex[0].k2);
        vec3  radiance    = lightColor * attenuation * cosTheta;
        

         // cook-torrance reflection brdf
        float NDF = DistributionGGX(N, H, roughness_value);        
        float G   = GeometrySmith(N, V, L, roughness_value);      
        vec3 F    = fresnelSchlick(max(dot(V, H), 0.0), F0); 

        vec3 numerator    = NDF * G * F;
        float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0);
        vec3 specular     = numerator / max(denominator, 0.001);  
     

        vec3 kS = F;
        vec3 kD = vec3(1.0) - kS;
        kD *= 1.0 - metallic_value;	    

         // add to outgoing radiance Lo
        float NdotL = max(dot(N, L), 0.0);  

        Lo += (kD * albedo_color / PI + specular) * radiance * NdotL; 

        //Lo = vec3(pass_View);
    }

    //ambient add
    //vec3 ambient = vec3(0.03) * brdf[0].albedo * brdf[0].ao;
    vec3 ambient = vec3(0.03) * albedo_color * ao_value;
    vec3 color = ambient + Lo;
	
    // gamma correction
    color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0/2.2));  

	// calc light.old program---------------------------
	//vec4 mixed = vec4(0.0,0.0,0.0,1.0);
	//for (int i=0; i<MAX_LIGHTS; i++){

		//if(light[i].used == false) continue;

		// light position if camera frame
		//vec4 L_c = viewMatrix * vec4(light[i].position, 1.0);
	
		// light to fragment 
		//vec3 L = normalize( L_c.xyz - pass_Position );
		//if(light[i].type == 2) L = light[i].direction;// direct light

		// checks whether the light was set.
		// Multiple lights blend adative
		//mixed += useLight( L,  E,  pass_Normal, light[i], mat[0]);
		
	//}

	//color = mixed;      
	// back to new program----------------------------
	frag_out = vec4(color, 1.0);  
	

	//------------------------------------------------
	// Get linear depth back

	// must match  the projection
	const float n = 0.01; // camera z near
    const float f = 10.0; // camera z far

	//   range without z/w -> [0, high value, linear but not scaled], with z/w -> range [-1,1]
	float current_depth = pass_Coordinates.z/ pass_Coordinates.w;
	          
	// Adjust the range to [0,1]
	//current_depth = 0.5 * current_depth + 0.5;  // range [0,1]

	// linear depth [0, depth]. requires range [-1, 1] as input.
	// Calculate the real depth
	current_depth  = (2.0 * f  * n) / (f + n - current_depth * (f - n)) ;  

	gl_FragDepth =  current_depth;

	                   
}                                                      