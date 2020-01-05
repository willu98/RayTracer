/*
*Name: William Pu
*215057060
*willu@my.yorku.ca
*/
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <list> 
#include <vector>
#include "glm/glm.hpp"
#include "glm/ext.hpp" 
#include <chrono> 
#include <algorithm> 

using namespace std;
using namespace std::chrono;
using namespace glm;

//coordinates of the near
float near;
float left1;
float right1;
float top;
float bottom;

//resolution of the image
int resolution[2];

//constant indentity matrix for the 4x4 
const glm::mat4 idMat = glm::mat4(1.0f);

//structure for a sphere
struct sphere {
	//name of the sphere
	string name;

	float position[3];
	float scaling[3];
	float color[3];

	//position, represented by a vec3
	glm::vec4 spherePos;

	//non-uniform scaling factor, rep by vec 3
	glm::vec4 sphereScaling;

	//color of the obj/sphere(rgb format), rep by vec3
	glm::vec3 sphereColor;

	//inversematrix of the sphere
	glm::mat4 sphereInverseMat;

	//coefficients for the sphere for the A.D.S model
	//0 represents K ambient
	//1 represents K diffuse
	//2 represents K specular
	//3 represents K reflective
	float kCoefficients[4];

	//exponent for the ads model
	int specularExp;

	//4x4 matrices for the sphere
	glm::mat4x4 transformMatrix = glm::mat4x4(1.0f);
	glm::mat4x4 invTransformMatrix = glm::mat4x4(1.0f);
};

//list of sphere, as there could be any number of spheres in the scene
vector <sphere> sphereObjects;

//structure for a light object
struct lightObj {
	//name of the light object
	string name;

	//light position for point source light
	glm::vec4 lightPos;

	//light intensity, rep by vec 3
	glm::vec3 lightIntensity;

};
//list of light objects as there could be any number of light objects int he scene
vector <lightObj> lightObjects;

//structure for rays
struct rays {

	//camera is situated at the origin, also where the rays start
	glm::vec4 origin = glm::vec4(0, 0, 0, 1);

	//vec3 representation of the direction of the ray
	glm::vec4 direction;

	//vec3 representation of the direction of the ray, applying the inverse transformation
	glm::vec4 directionTrans;

	//vec3 representation of the origin  of the ray applying M^-1
	glm::vec4 originTrans;

	//vec4 representing the intersections
	glm::vec4 closestIntersection;

	//boolean for wether or not a collision is found for the ray
	bool foundCollision = false;

	//float value for the smallest T value
	float smallestT = std::numeric_limits<float>::max();

	//depth of ray initailized to one, max 3 according to assighnment specifications
	int rayDepth = 0;

	//keeps track of which sphere the ray collided with in order to get the color
	int whichSpheres;
};

//background color
float backColor[3];

//scene's ambient intensity
float ambientIntensity[3];

//the file that is to be output
string outputFile;
char outputFileCh[20];

//char array for pixels
unsigned char* pixels;

void fileReader(string fileName);
vec4 detectRayCollision(rays& ray, sphere& sphere, int currentSphere);
glm::vec3 mainRayTracingMethod(rays& ray);
void save_image(int Width, int Height, char* fname, unsigned char* pixels);


int main(int argc, char* arg[])
{

	auto start = high_resolution_clock::now();

	//variables for height and width of the image plane
	float height;
	float width;
	//fileReader(arg[1]);
	fileReader("Project3TestsAndKeys/testBehind.txt");

	cout << "\nCURRENTLY: DONE READING FILE\n" << endl;

	cout << "CURRENTLY: GETTING THE TRANSFORMED MATRICES & PREPARING RAYTRACING" << endl;

	//getting the transformationmatrix for each of the spheres
	//and theuir respective inverse matrices
	for (int i = 0; i < sphereObjects.size(); i++) {
		glm::mat4 transMatrix = glm::translate(idMat, glm::vec3(sphereObjects[i].position[0], sphereObjects[i].position[1], sphereObjects[i].position[2])) * glm::scale(idMat, glm::vec3(sphereObjects[i].scaling[0], sphereObjects[i].scaling[1], sphereObjects[i].scaling[2]));
		sphereObjects[i].transformMatrix = transMatrix;
		//transMat = glm::transpose(transMat);
		sphereObjects[i].invTransformMatrix = glm::inverse(sphereObjects[i].transformMatrix);

	}


	//getting dimenions of the image
	width = abs(left1 - right1) / (float)2;
	height = abs(top - bottom) / (float)2;
	cout << "IMAGE PLANE WIDTH: " << width << " " << endl;
	cout << "IMAGE PLANE HEIGHT: " << height << " " << endl;

	//array size of 
	pixels = new unsigned char[resolution[0] * resolution[1] * 3];
	//represents the current index for when saving th ppm file
	int currentImgIndex = 0;


	cout << "STARTING CASTING RAYS" << endl;

	//iterating to each pixel and casting a ray in each pixel
	for (int i = resolution[0]; i > 0; i--)
	{
		for (int j = 0; j < resolution[1]; j++) {

			//creating a temp ray
			rays tempRay;

			//getting the coordinates of the pixels in the camera coordinate system using the formulas
			//Px = -w + w * ((2 * i)/resolution.y)
			//Py = -h + h * ((2 * j)/resolution.y)
			float pixelX = -width + width * ((2 * j) / (float)resolution[0]);
			float pixelY = -height + height * ((2 * i) / (float)resolution[1]);
			//cout << "PIXEL X: " << pixelX << endl;

			//by default the color found by the ray tracer is the background color
			glm::vec3 color = glm::vec3(backColor[0], backColor[1], backColor[2]);

			//direction vector for the ray
			tempRay.direction = glm::vec4(pixelX, pixelY, -near, 0);

			//getting the color for the pixel
			color = mainRayTracingMethod(tempRay);

			pixels[currentImgIndex] = glm::min(255.f, color.r * 255.f);
			pixels[currentImgIndex + 1] = glm::min(255.f, color.g * 255.f);
			pixels[currentImgIndex + 2] = glm::min(255.f, color.b * 255.f);
			currentImgIndex += 3;
		}
	}
	cout << "DONE CASTING RAYS" << endl;

	//saving the image
	save_image(resolution[0], resolution[1], outputFileCh, pixels);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Time taken by function: " << duration.count() / 1000 << " milliseconds" << endl;

	return 0;
}


//saves the image as a ppm file
//code provided for the assighnment
void save_image(int Width, int Height, char* fname, unsigned char* pixels) {
	FILE* fp;
	const int maxVal = 255;

	printf("Saving image %s: %d x %d\n", fname, Width, Height);
	fp = fopen(fname, "wb");
	if (!fp) {
		printf("Unable to open file '%s'\n", fname);
		return;
	}
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);

	for (int j = 0; j < Height; j++) {
		fwrite(&pixels[j * Width * 3], 3, Width, fp);
	}

	fclose(fp);
}



//detects if there is a collision between a ray and another object
//returns a vec4 for where a collision is detetcted
//ray param is the ray being used to check the collision
//sphere is the sphere used to check for collision with the ray
//sphereNum is the index of the sphere, of which sphere is currently being checked
vec4 detectRayCollision(rays& ray, sphere& sphere, int currentSphere) {
	//getting the values for M^-1(rayDir)
	//applying the inverse transformation to the direction of the ray
	ray.directionTrans = sphere.invTransformMatrix * ray.direction;

	//getting the values for M^-1(rayOrigin)
	//applying the inverse transformation to the origin of the ray
	ray.originTrans = sphere.invTransformMatrix * ray.origin;
	ray.originTrans[3] = 0;

	//getting a, b, c for the quadratic equation of form
	//-B/A +- sqrt(b^2-ac)/a
	//A = |dir|^2 = dir dot dir, by property of dot product
	float a = glm::dot(ray.directionTrans, ray.directionTrans);
	//B = Point dot dir
	float b = glm::dot(ray.directionTrans, ray.originTrans);
	//C = |origin|^2 = origin dot origin, by property of dot product
	float c = glm::dot(ray.originTrans, ray.originTrans) - 1.f;

	//calculating b^2-ac
	float discriminant = pow(b, 2) - (a * c);

	//value of t, the solution for the quadratic formula 
	float tInt1 = 0;
	float tInt2 = 0;

	//only one solution
	if (discriminant == 0) {
		//getting th
		tInt1 = -(b / a);
		//cout << "tInt1: " << tInt1 << endl;

		//checking for behind
		if (ray.closestIntersection[2] < -1)
		{
			if (!ray.foundCollision) {
				ray.foundCollision = false;
			}
			return ray.closestIntersection;
		}

		//the smallest t value for the sphere is being checked for, as there could be multiple collisions
		if (ray.smallestT > tInt1) {
			ray.smallestT = tInt1;
			ray.whichSpheres = currentSphere;
		}
		ray.foundCollision = true;

		//clostest intersection with the ray is found using the smaller of the two t's
		return ray.smallestT * ray.direction + ray.origin;
	}
	//two solutions/intersections
	else if (discriminant > 0) {
		//getting both t values
		tInt1 = -(b / a) + sqrt(discriminant) / a;
		tInt2 = -(b / a) - sqrt(discriminant) / a;

		//checking for behind
		if (ray.closestIntersection[2] > 1)
		{
			if (!ray.foundCollision) {
				ray.foundCollision = false;
			}
			return ray.closestIntersection;
		}

		//getting smallest of two t's
		float t = tInt1 < tInt2 ? tInt1 : tInt2;


		//do a comparison to find the smallest t for the ray
		//the smallest t value for the sphere is being checked for, as there could be multiple collisions
		if (ray.smallestT > t) {
			ray.smallestT = t;
			ray.whichSpheres = currentSphere;
		}

		ray.foundCollision = true;
		//clostest intersection with the ray is found using the smallest t
		return ray.smallestT * ray.direction + ray.origin;
	}
	if (!ray.foundCollision) {
		ray.foundCollision = false;
	}
	//if no collision is found then return false fro collision detection
	return ray.closestIntersection;
}


//main raytracing method, uses the ray
//takes a ray as a parameter
glm::vec3 mainRayTracingMethod(rays& ray) {

	ray.rayDepth++;

	//if the ray depth has reached at least 3 simply return black
	if (ray.rayDepth >= 3) {
		return glm::vec3(0.0f, 0.0f, 0.0f);
	}


	//checking for collisions with each object in the space
	for (int i = 0; i < sphereObjects.size(); i++) {

		//calling collision detection method
		//and checking for collision
		//getting the closest inmtersection for the sphere
		ray.closestIntersection = detectRayCollision(ray, sphereObjects[i], i);
	}

	//if no collision is found or the collision is behind the img plane 
	//then return the background color
	if (ray.foundCollision == false) {
		return glm::vec3(backColor[0], backColor[1], backColor[2]);
	}


	//the transpose of the inverse of the sphere's inverse matrix
	glm::mat4 transInvTransformMatrix = glm::transpose(sphereObjects[ray.whichSpheres].invTransformMatrix);

	//getting the collision point transformed onto the cononical sphere of radius one
	vec4 NVector = sphereObjects[ray.whichSpheres].invTransformMatrix * ray.closestIntersection;
	NVector = normalize(transInvTransformMatrix * NVector);
	NVector[3] = 0;

	//color vector for local illumination
	glm::vec3 cLocal = glm::vec3(0.0f, 0.0f, 0.0f);


	//iterating through each light object
	//to check to see if the rays are blocked by an object
	for (int i = 0; i < lightObjects.size(); i++) {

		//L vector is the vector from the point to the light 
		vec4 LVector = normalize(lightObjects[i].lightPos - ray.closestIntersection);


		//view vector 
		vec4 VVector = -ray.closestIntersection;
		VVector[3] = 0;
		VVector = normalize(VVector);

		//reflection of the light
		vec4 RVector = normalize(2.0f * (glm::max(dot(normalize(NVector), LVector), 0.0f)) * normalize(NVector) - LVector);

		//shadow ray var
		rays shadowRay;
		
		//direcional vector from the point to the light object
		//found by subtracting the point of collision from 1st object and the position of the light
		//then normalizing the vector to get the normalized direction vector
		shadowRay.direction = LVector;

		//getting the origin of the shadow ray and setting it to the intersection point 
		shadowRay.origin = ray.closestIntersection;

		//if a collision is not found then return color based on local illumination method
		//checking to see if the ray from the first collision point to the light has any collisions inbetween
		for (int k = 0; k < sphereObjects.size(); k++) {

			if (k != ray.whichSpheres) {
				shadowRay.closestIntersection = detectRayCollision(shadowRay, sphereObjects[k], k);
			}
		}
		if (!shadowRay.foundCollision) {

			//returns the local illumination
			//calulating Kdiffuse * (Norm dot shadowRay) * Sphere.red * Ired
			vec3 colorTemp;
			//kd = Light vector dot sphere normal vector
			float kD = 0;
			if (sphereObjects[ray.whichSpheres].kCoefficients[1] > 0) {
				kD = glm::max(dot(normalize(LVector), normalize(NVector)), 0.0f) * sphereObjects[ray.whichSpheres].kCoefficients[1];
			}

			//ks = (reflection vector dot view vector)^n
			float kS = 0;
			if (sphereObjects[ray.whichSpheres].kCoefficients[2] > 0) {
				kS = glm::max(dot(normalize(RVector), normalize(VVector)), 0.0f);
				kS = pow(kS, sphereObjects[ray.whichSpheres].specularExp) * sphereObjects[ray.whichSpheres].kCoefficients[2];
			}
			//cout << "DOT:" << glm::max(dot(normalize(LVector), normalize(NVector)), 0.0f) << endl;
			for (int n = 0; n < 3; n++) {
				colorTemp[n] = kD * lightObjects[i].lightIntensity[n] * sphereObjects[ray.whichSpheres].color[n] + kS * lightObjects[i].lightIntensity[n];
			}

			//get the color from light source i, using ads model
			cLocal += colorTemp;
		}
	}

	//calculating the ambient color
	glm::vec3 ambientColor = sphereObjects[ray.whichSpheres].kCoefficients[0] * glm::vec3(ambientIntensity[0] * sphereObjects[ray.whichSpheres].color[0], ambientIntensity[1] * sphereObjects[ray.whichSpheres].color[1], ambientIntensity[2] * sphereObjects[ray.whichSpheres].color[2]);

	vec3 colorReflect = vec3(0.0f, 0.0f, 0.0f);
	//calculating the reflection iff the reflection coeff > 0
	if (sphereObjects[ray.whichSpheres].kCoefficients[3] > 0) {
		//ray obj for the reflected ray
		rays reflectedRay;

		//setting reflected ray's origin
		reflectedRay.origin = ray.closestIntersection;

		//getting the directional vector for the reflected ray
		reflectedRay.direction = normalize(-2.f * (glm::max(dot(normalize(NVector), normalize(ray.direction)), 0.0f) * (normalize(NVector)) + normalize(ray.direction)));


		//setting ray depth for recursive call
		reflectedRay.rayDepth = ray.rayDepth;
		//recursive call to the reflection ray, iff the ray depth < 3
		colorReflect = sphereObjects[ray.whichSpheres].kCoefficients[3] * mainRayTracingMethod(reflectedRay);
	}

	//adding all the color values from the various lighting components
	return ambientColor + cLocal + colorReflect;
}


//this method does the file reading, and looks for key words to obtain the appropiate values
//fileName parameter is the name of the file that is to be read
void fileReader(string fileName) {
	ifstream inFile(fileName);
	string newLine;
	string keyWord;
	while (inFile)
	{
		//reading the file line by line
		getline(inFile, newLine);

		//if the current line is not an empty line, otherwise just continue to read to next line
		if ((int)newLine[0] != 0) {
			if (newLine.find(9) != std::string::npos && newLine.find(" ") > newLine.find(9)) {
				//getting the first word in the line, ie the key word
				keyWord = newLine.substr(0, newLine.find(9));;
				//removing keyword from the string
				newLine = newLine.substr(newLine.find(9) + 1, newLine.size());
			}

			if (newLine.find(" ") != std::string::npos && newLine.find(9) > newLine.find(" ")) {
				//getting the first word in the line, ie the key word
				keyWord = newLine.substr(0, newLine.find(" "));;
				//removing keyword from the string
				newLine = newLine.substr(newLine.find(" ") + 1, newLine.size());
			}

			//the split part of each string that will be used as the individual values
			string token;



			//if near is the keyword that is found
			if ((keyWord.compare("NEAR")) == 0) {

				//while there is white space at the front, continue to remove
				while (newLine[0] == ' ' || newLine[0] == 9) {
					//getting the substring without the white space at the front
					newLine = newLine.substr(1, newLine.size());
				}

				//getting the first number inline
				token = newLine.substr(0, newLine.find(" "));
				near = std::stof(token);
				cout << "NEAR: " << near << endl;
			}
			//if left key word is found
			else if ((keyWord.compare("LEFT")) == 0) {
				//while there is white space at the front, continue to remove
				while (newLine[0] == ' ' || newLine[0] == 9) {
					//getting the substring without the white space at the front
					newLine = newLine.substr(1, newLine.size());
				}

				//getting the first number inline
				token = newLine.substr(0, newLine.find(" "));
				left1 = std::stof(token);
				cout << "LEFT: " << left1 << endl;
			}
			//if right key word is found
			else if ((keyWord.compare("RIGHT")) == 0) {
				//while there is white space at the front, continue to remove
				while (newLine[0] == ' ' || newLine[0] == 9) {
					//getting the substring without the white space at the front
					newLine = newLine.substr(1, newLine.size());
				}

				//getting the first number inline
				token = newLine.substr(0, newLine.find(" "));
				right1 = std::stof(token);
				cout << "RIGHT: " << right1 << endl;
			}
			//if key word bottom is found
			else if ((keyWord.compare("BOTTOM")) == 0) {
				//while there is white space at the front, continue to remove
				while (newLine[0] == ' ' || newLine[0] == 9) {
					//getting the substring without the white space at the front
					newLine = newLine.substr(1, newLine.size());
				}

				//getting the first number inline
				token = newLine.substr(0, newLine.find(" "));
				bottom = std::stof(token);
				cout << "BOTTOM: " << bottom << endl;
			}
			//if key word top is found
			else if ((keyWord.compare("TOP")) == 0) {
				//while there is white space at the front, continue to remove
				while (newLine[0] == ' ' || newLine[0] == 9) {
					//getting the substring without the white space at the front
					newLine = newLine.substr(1, newLine.size());
				}

				//getting the first number inline
				token = newLine.substr(0, newLine.find(" "));
				top = std::stof(token);
				cout << "TOP: " << top << endl;
			}
			//if keyword res is found
			else if ((keyWord.compare("RES")) == 0) {
				for (int i = 0; i < 2; i++) {
					//while there is white space at the front, continue to remove
					while (newLine[0] == ' ' || newLine[0] == 9) {
						//getting the substring without the white space at the front
						newLine = newLine.substr(1, newLine.size());
					}

					//getting the first number inline
					if (newLine.find(" ") != std::string::npos) {
						token = newLine.substr(0, newLine.find(" "));
					}
					else if (newLine.find(9) != std::string::npos) {
						token = newLine.substr(0, newLine.find(9));
					}

					resolution[i] = std::stoi(token);
					if (i < 1) {
						//going to the next value
						//the character before the next value is a space
						if (newLine.find(" ") != std::string::npos) {
							//if space comes before tab
							if (newLine.find(" ") < newLine.find(9)) {
								//move past the space
								newLine = newLine.substr(newLine.find(" "), newLine.size());
							}
						}

						//if the tab is front of the next value
						if (newLine.find(9) != std::string::npos) {
							//if tab comes before space
							if (newLine.find(9) < newLine.find(" ")) {
								//move past the tab
								newLine = newLine.substr(newLine.find(9), newLine.size());
							}
						}
					}
					cout << "RESOULTION: " << resolution[i] << endl;
				}
			}
			//if keyword ambient is found
			else if ((keyWord.compare("AMBIENT")) == 0) {
				for (int i = 0; i < 3; i++) {
					//while there is white space at the front, continue to remove
					while (newLine[0] == ' ' || newLine[0] == 9) {
						//getting the substring without the white space at the front
						newLine = newLine.substr(1, newLine.size());
					}

					//getting the first number inline
					if (newLine.find(" ") != std::string::npos) {
						token = newLine.substr(0, newLine.find(" "));
					}
					else if (newLine.find(9) != std::string::npos) {
						token = newLine.substr(0, newLine.find(9));
					}

					ambientIntensity[i] = std::stof(token);
					if (i < 2) {
						//going to the next value
						//the character before the next value is a space
						if (newLine.find(" ") != std::string::npos) {
							//if space comes before tab
							if (newLine.find(" ") < newLine.find(9)) {
								//move past the space
								newLine = newLine.substr(newLine.find(" "), newLine.size());
							}
						}

						//if the tab is front of the next value
						if (newLine.find(9) != std::string::npos) {
							//if tab comes before space
							if (newLine.find(9) < newLine.find(" ")) {
								//move past the tab
								newLine = newLine.substr(newLine.find(9), newLine.size());
							}
						}
					}
					cout << "ambientIntensity: " << ambientIntensity[i] << endl;
				}
			}
			//if key word back is found
			else if ((keyWord.compare("BACK")) == 0) {
				for (int i = 0; i < 3; i++) {
					//while there is white space at the front, continue to remove
					while (newLine[0] == ' ' || newLine[0] == 9) {
						//getting the substring without the white space at the front
						newLine = newLine.substr(1, newLine.size());
					}

					//getting the first number inline
					if (newLine.find(" ") != std::string::npos) {
						token = newLine.substr(0, newLine.find(" "));
					}
					else if (newLine.find(9) != std::string::npos) {
						token = newLine.substr(0, newLine.find(9));
					}


					backColor[i] = std::stof(token);
					if (i == 2) {
						backColor[i] = std::stof(newLine);
					}
					if (i < 2) {
						//going to the next value
						//the character before the next value is a space
						if (newLine.find(" ") != std::string::npos) {
							//if space comes before tab
							if (newLine.find(" ") < newLine.find(9)) {
								//move past the space
								newLine = newLine.substr(newLine.find(" "), newLine.size());
							}
						}

						//if the tab is front of the next value
						if (newLine.find(9) != std::string::npos) {
							//if tab comes before space
							if (newLine.find(9) < newLine.find(" ")) {
								//move past the tab
								newLine = newLine.substr(newLine.find(9), newLine.size());
							}
						}
					}
					cout << "BACKGROUND COLOR: " << backColor[i] << endl;
				}
			}
			//IF SPHERE KEYWORD IS FOUND
			else if ((keyWord.compare("SPHERE")) == 0) {
				//temporary sphere object to push into the list
				sphere tempSphere;

				//GETTING POSITION, COLOR AND SCALING FACTORS OF THE SPHERE
				for (int i = 0; i < 15; i++) {
					//while there is white space at the front, continue to remove
					while (newLine[0] == ' ' || newLine[0] == 9) {
						//getting the substring without the white space at the front
						newLine = newLine.substr(1, newLine.size());
					}
					//SETTING the value number inline
					if (newLine.find(" ") != std::string::npos) {
						token = newLine.substr(0, newLine.find(" "));
					}
					else if (newLine.find(9) != std::string::npos) {
						token = newLine.substr(0, newLine.find(9));
					}

					if (i == 0) {
						//SETTING the name of the sphere
						tempSphere.name = token;
						cout << "SPHERE NAME: " << tempSphere.name << endl;
					}
					//SETTING THE (X,Y,Z) POSITIONS OF THE SPHERE
					else if (i < 4) {
						tempSphere.position[i - 1] = std::stof(token);
						cout << "SPHERE POSITION: " << tempSphere.position[i - 1] << endl;
					}
					//SETTING THE (X,Y,Z) SCALING FACTORS OF THE SPHERE
					else if (i < 7) {
						tempSphere.scaling[i - 4] = std::stof(token);
						cout << "SPHERE SCALING: " << tempSphere.scaling[i - 4] << endl;
					}
					//SETTING THE (R,G,B) VALUES OF THE COLOR
					else if (i < 10) {
						tempSphere.color[i - 7] = std::stof(token);
						cout << "SPHERE COLOR: " << tempSphere.color[i - 7] << endl;
					}
					//SETTING THE COEFFICIENTS FOR THE ADS MODEL
					else if (i < 14) {
						tempSphere.kCoefficients[i - 10] = std::stof(token);
						cout << "SPHERE COEFFICIENTS: " << tempSphere.kCoefficients[i - 10] << endl;
					}
					//SETTING THE SPECULAR EXPONENT VALUE
					else if (i == 14) {
						tempSphere.specularExp = std::stoi(newLine);
						cout << "SPHERE SPEC. EXP: " << tempSphere.specularExp << endl;
					}


					if (i < 14) {
						//going to the next value only if not at the last part of the line
						//the character before the next value is a space
						if (newLine.find(" ") != std::string::npos) {
							//if space comes before tab
							if (newLine.find(" ") < newLine.find(9)) {
								//move past the space
								newLine = newLine.substr(newLine.find(" "), newLine.size());
							}
						}

						//if the tab is front of the next value
						if (newLine.find(9) != std::string::npos) {
							//if tab comes before space
							if (newLine.find(9) < newLine.find(" ")) {
								//move past the tab
								newLine = newLine.substr(newLine.find(9), newLine.size());
							}
						}
					}
				}

				tempSphere.spherePos = glm::vec4(tempSphere.position[0], tempSphere.position[1], tempSphere.position[2], 1);
				tempSphere.sphereScaling = glm::vec4(tempSphere.scaling[0], tempSphere.scaling[1], tempSphere.scaling[2], 0);
				tempSphere.sphereColor = glm::vec3(tempSphere.color[0], tempSphere.color[1], tempSphere.color[2]);

				//adding a sphere object to the list
				sphereObjects.push_back(tempSphere);

			}
			//IF LIGHT KEYWORD IS FOUND
			else if ((keyWord.compare("LIGHT")) == 0) {
				//temporary sphere object to push into the list
				lightObj tempLight;

				//light position for point source light
				float posLight[3];
				//intensity of point source light
				float intensityLight[3];

				//GETTING POSITION, intensity and name FACTORS OF THE light
				for (int i = 0; i < 7; i++) {
					//while there is white space at the front, continue to remove
					while (newLine[0] == ' ' || newLine[0] == 9) {
						//getting the substring without the white space at the front
						newLine = newLine.substr(1, newLine.size());
					}
					//SETTING the value number inline
					if (newLine.find(" ") != std::string::npos) {
						token = newLine.substr(0, newLine.find(" "));
					}
					else if (newLine.find(9) != std::string::npos) {
						token = newLine.substr(0, newLine.find(9));
					}
					if (i == 0) {
						//SETTING the name of the light
						tempLight.name = token;
						cout << "LIGHT NAME: " << tempLight.name << endl;
					}
					//SETTING THE (X,Y,Z) POSITIONS OF THE LIGHT
					else if (i < 4) {
						posLight[i - 1] = std::stof(token);
						cout << "LIGHT POSITION: " << posLight[i - 1] << endl;
					}
					//SETTING THE (X,Y,Z) SCALING FACTORS OF THE INTENSITY
					else if (i < 7) {
						if (i != 6) {
							intensityLight[i - 4] = std::stof(token);
						}
						else
						{
							intensityLight[i - 4] = std::stof(newLine);
						}
						cout << "LIGHT INTENSITY: " << intensityLight[i - 4] << endl;
					}

					if (i < 6) {
						//going to the next value only if not at the last part of the line
						//the character before the next value is a space
						if (newLine.find(" ") != std::string::npos) {
							//if space comes before tab
							if (newLine.find(" ") < newLine.find(9)) {
								//move past the space
								newLine = newLine.substr(newLine.find(" "), newLine.size());
							}
						}
						//if the tab is front of the next value
						if (newLine.find(9) != std::string::npos) {
							//if tab comes before space
							if (newLine.find(9) < newLine.find(" ")) {
								//move past the tab
								newLine = newLine.substr(newLine.find(9), newLine.size());
							}
						}
					}
				}
				tempLight.lightIntensity = glm::vec3(intensityLight[0], intensityLight[1], intensityLight[2]);
				tempLight.lightPos = glm::vec4(posLight[0], posLight[1], posLight[2], 1);

				//adding a sphere object to the list
				lightObjects.push_back(tempLight);

			}

			//if the keyword found is OUTPUT
			else if ((keyWord.compare("OUTPUT ")) == 0 || (keyWord.compare("OUTPUT")) == 0) {

				//while there is white space at the front, continue to remove
				while (newLine[0] == ' ' || newLine[0] == 9) {
					//getting the substring without the white space at the front
					newLine = newLine.substr(1, newLine.size());
				}
				//SETTING the name of the output file
				outputFile = newLine;
				strcpy(outputFileCh, outputFile.c_str());
				cout << "OUTPUT FILE NAME: " << outputFile<< endl;
			}
			newLine = "";
			cout << "KEYWORD: " << keyWord << endl;
			cout << endl;

		}

	}
	inFile.close();
}



