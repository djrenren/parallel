#include <iostream>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <stdlib.h>
double inf;
using namespace std;
class Vector{
	public:
		double x,y,z;
		Vector(){z=0;y=0;z=0;}
		Vector(double ix, double iy, double iz){
			x = ix; y = iy; z = iz;
		}
		Vector(Vector* v1, Vector* v2){
			x = v2->x - v1->x;
			y = v2->y - v1->y;
			z = v2->z - v1->z;
		}
		Vector* normalize(){
			double invlen = 1.0/this->getMagnitude();
			x *= invlen; y *= invlen; z *= invlen;
			return this;
		}
		double getMagnitude(){
			return sqrt(x*x+y*y+z*z);
		}
		double dot(Vector* v){
			return x*v->x + y*v->y + z*v->z;
		}
		Vector operator*=(const double scalar){
			x*= scalar; y*= scalar; z*= scalar;
			return *this;
		}
		Vector operator+=(const Vector& v){
			x+= v.x; y+= v.y; z+= v.z;
			return *this;
		}
		Vector operator-=(const Vector& v){
			x-= v.x; y-= v.y; z-= v.z;
			return *this;
		}
		void print(){
			printf("(%f, %f, %f)",x,y,z);
		}
};

class Color{
	private:
		int r,g,b;
	public:
		Color(){r=0;g=0;b=0;}
		Color(Color* c){
			r=c->r; g=c->g; b=c->b;
		}
		Color(int ir, int ig, int ib){
			r = ir; g = ig; b = ib;
		}
		Color* add(int r, int g, int b){
			this->r = r+this->r > 255 ? 255 : r+this->r;
			this->g = g+this->g > 255 ? 255 : g+this->g;
			this->b = b+this->b > 255 ? 255 : b+this->b;
			return this;
		}
		Color operator+=(const Color c){
			return *(this->add(c.r,c.g,c.b));
		}
		Color operator*(const double s){
			return *(new Color(r*s,g*s,b*s));
		}
		void print(){
			printf("%d %d %d",r,g,b);
		}
};

class Object{
	protected:
		Vector* loc;
		Color* rgb;
	public:
		double reflectivity;
		Object(){
			this->rgb = new Color();
			this->loc = new Vector();
			this->reflectivity = 0;
		}
		Object(Vector* loc, Color* rgb, double ref){
			this->rgb = rgb; this->loc = loc;
			this->reflectivity = ref;
		}
		virtual Color* getColor(Vector*){
			return rgb;
		}
		virtual double getCollision(Vector*, Vector*) = 0;
		Color* getColor(Vector*, Object**, Vector*){
			return this->rgb;
		};
		virtual Vector* getNormal(Vector*) = 0;
};

class Sphere: public Object{
	double r;
	public:
		Sphere(Vector* loc, double r, Color* col, double ref)
		:Object(loc,col,ref){
			this->r = r;
		}
		double getCollision(Vector*, Vector*);
		Vector* getNormal(Vector*);
};

double Sphere::getCollision(Vector* eye, Vector* ray){
	double b = 2*(ray->x*(eye->x-loc->x) + ray->y*(eye->y-loc->y) + ray->z*(eye->z-loc->z));
	/*double c = (eye->x-loc->x)*(eye->x-loc->x) + 
				(eye->y-loc->y)*(eye->y-loc->y) +
				(eye->z-loc->z)*(eye->z-loc->z) - 
				r*r;
				*/
	double c = loc->x*loc->x + loc->y*loc->y + loc->z*loc->z + 
				eye->x*eye->x + eye->y*eye->y + eye->z*eye->z -
				r*r - 2*(eye->x*loc->x+eye->y*loc->y+eye->z*loc->z);
	double det = b*b-4*c;
	if(det < 0)
		return inf;
	double pos = (-b+sqrt(det))/2;
	double neg = (-b-sqrt(det))/2;
	return pos > neg ? neg : pos;
}

Vector* Sphere::getNormal(Vector* p){
	return (new Vector(loc,p))->normalize();
}
class Floor: public Object{
	private:
		Color* rgb2;
	public:
		Floor(double y,Color* col, Color* col2, double ref){
			this->loc = new Vector(0,y,0);
			this->rgb = col;
			this->rgb2 = col2;
			this->reflectivity = ref;
		};
		double getCollision(Vector*, Vector*);
		Vector* getNormal(Vector*);
		Color* getColor(Vector*);
};

Color* Floor::getColor(Vector* coord){
	if((int)fabs(ceil(coord->x/0.5))%2 != (int)fabs(floor(coord->z/0.5))%2 )
		return rgb;
	return rgb2;
}

double Floor::getCollision(Vector* eye, Vector* ray){
	double ret = (loc->y-eye->y)/ray->y;
	return ret > 0 ? ret : inf;
}

Vector* Floor::getNormal(Vector* col){
	return new Vector(0,1,0);
}





class Screen{
		double zindex;
		double width, height;
		Color* bg;
	public:
		Screen(double,double,double, Color*);
		void render(int,int,Vector*, Vector*, Object**, int);
		Color* castRay(Vector*, Vector*, Vector*, Object**, int, int);
		double getShadow(Vector*, Vector*, int,Object**, int);
};

Screen::Screen(double z, double w, double h, Color* bg){
	width = w;
	height = h;
	zindex = z;
	this->bg = bg;
}

Color* Screen::castRay(Vector* origin, Vector* dir, Vector* light, Object** objects, int numObj, int depth){
	double minb = inf;
	Object* closeObj = objects[0];
	int closei = 0;
	for(int o=0;o<numObj;o++){
		double b = objects[o]->getCollision(origin,dir);
		//cout << b << "\n";
		if(b > 0.001 && b < minb){
			minb = b;
			closeObj = objects[o];
			closei = o;
		}
	}
	//cout << minb << "\n";
	if(minb < inf){
		*dir*=minb;
		*dir+=*origin;
		Color* outCol = new Color();
		*outCol += *(closeObj->getColor(dir))*0.4;
		Vector* lv = new Vector(dir,light);
		double ldist = lv->getMagnitude();
		lv->normalize();
		Vector* normal = closeObj->getNormal(dir);
		double dot = normal->dot(lv);
		if(dot > 0 && ldist - getShadow(dir,lv,closei,objects,numObj) > 0){
			*outCol += *(closeObj->getColor(dir))*0.6*dot;
		}
		if(depth > 0 && closeObj->reflectivity > 0){
			origin = new Vector(new Vector(),dir);
			dir->normalize();
			*outCol = *outCol*(1-closeObj->reflectivity);
			*outCol += *(castRay(origin, (*dir -= (*normal *= dir->dot(normal)*2)).normalize(), light, objects, numObj, depth--))*closeObj->reflectivity;
		}
		return outCol;
	}
	return bg;
	

}

void Screen::render(int px, int py, Vector* eye, Vector* light, Object** objects, int numObj){
	Color* imgOut[py][px];
	double pix = width/px;
	double piy = height/py;
	
	int numshadow = 0;
	printf("P3\n%i %i\n255\n",px,py);
	#pragma omp parallel for
	for(int r=0; r<py; r++)
		for(int c=0; c<px; c++){
			Vector* dir = (new Vector(pix*c-eye->x, 
								height-piy*r-eye->y, 
								zindex-eye->z))->normalize();
			imgOut[r][c] = castRay(eye,dir, light, objects, numObj,1);
			//delete dir;
		}
	for(int r=0; r<py; r++)
		for(int c=0; c<px; c++){
			imgOut[r][c]->print();
			cout << ' ';
		}
}
double Screen::getShadow(Vector* sv, Vector* lv, int close, Object** obs, int numObj){

	for(int i=0; i < numObj; i++){
	double col = obs[i]->getCollision(sv,lv);
		if(col > 0 && col < inf){
			return inf;
		}
	}
	return 0;
}

int main(int argc, char** argv) {
	inf = numeric_limits<double>::infinity();
	Vector* eye = new Vector(0.5,0.5,-1);
	Vector* light = new Vector(0,1.0,0.5);
	Screen* screen = new Screen(0.0,1.5,1.0, new Color(0,0,0));
	Color* red = new Color(255,0,0);
	Color* blue = new Color(0,0,255);
	Color* green = new Color(0,255,0); 
	Color* brown = new Color(184,134,11);
	Color* indigo = new Color(75, 0, 130);
	Sphere* s1 = new Sphere(new Vector(0.5,0.166667,0.666667),0.166667,red,0.1);
	Sphere* s2 = new Sphere(new Vector(0.833333,0.866667,1.0),0.166667,blue,0);
	Sphere* s3 = new Sphere(new Vector(0.633333,0.333333,1.166667),0.333333,green,0.2);
	Floor* f = new Floor(0,brown,indigo,0.3);
	Object** obs = (Object**)malloc(4*sizeof(Object*));
	obs[0] = s1;
	obs[1] = s2;
	obs[2] = s3;
	obs[3] = f;
	screen->render(900,600,eye,light,obs,4);
}
