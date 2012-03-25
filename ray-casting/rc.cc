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
		Object(){
			this->rgb = new Color();
			this->loc = new Vector();
		}
		Object(Vector* loc, Color* rgb){
			this->rgb = rgb; this->loc = loc;
		}
		Color* getColor(Vector*){
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
		Sphere(Vector* loc, double r, Color* col)
		:Object(loc,col){
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
		Floor(double y,Color* col, Color* col2){
			this->loc = new Vector(0,y,0);
			this->rgb = col;
			this->rgb2 = col2;
		};
		double getCollision(Vector*, Vector*);
		Vector* getNormal(Vector*);
		Color* getColor(Vector*);
};

Color* Floor::getColor(Vector* coord){
	if((int)(coord->x/0.008) % 2 == 1 )
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
		double getShadow(Vector*, Vector*, int,Object**, int);
};

Screen::Screen(double z, double w, double h, Color* bg){
	width = w;
	height = h;
	zindex = z;
	this->bg = bg;
}

void Screen::render(int px, int py, Vector* eye, Vector* light, Object** objects, int numObj){
	Color* imgOut[py][px];
	double pix = width/px;
	double piy = height/py;
	double minb;
	Object* closeObj = objects[0];
	int closei = 0;
	int numshadow = 0;
	printf("P3\n%i %i\n255\n",px,py);
	#pragma omp parallel for private(minb,closei,closeObj)
	for(int r=0; r<py; r++)
		for(int c=0; c<px; c++){
			minb = inf;
			closei = 0;
			Vector* v = (new Vector(pix*c-eye->x, 
								height-piy*r-eye->y, 
								zindex-eye->z))->normalize();
			for(int o=0;o<numObj;o++){
				double b = objects[o]->getCollision(eye,v);
				if(b > 0 && b < minb){
					minb = b;
					closeObj = objects[o];
					closei = o;
				}
			}
			if(minb < inf){
				//cout << minb << "\n";
				//v->print();
				//cout << " * " << minb << " = ";
				*v*=minb;
				*v+=*eye;
				//v->print();
				//cout << "\n";
				Color* outCol = new Color();
				*outCol += *(closeObj->getColor(v))*0.4;
				Vector* lv = new Vector(v,light);
				double ldist = lv->getMagnitude();
				lv->normalize();
				//v->print();
				//cout << "---";
				//lv->print();
				double dot = closeObj->getNormal(v)->dot(lv);
				if(dot > 0 && ldist - getShadow(v,lv,closei,objects,numObj) > 0){
					*outCol += *(closeObj->getColor(v))*0.6*dot;
				}

				imgOut[r][c] = outCol;
			}
			else
				imgOut[r][c] = bg;
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
	Color* brown = new Color(139,119,101);
	Sphere* s1 = new Sphere(new Vector(0.5,0.166667,0.666667),0.166667,red);
	Sphere* s2 = new Sphere(new Vector(0.833333,0.866667,1.0),0.166667,blue);
	Sphere* s3 = new Sphere(new Vector(0.333333,0.333333,1.166667),0.333333,green);
	Floor* f = new Floor(0,brown,green);
	Object** obs = (Object**)malloc(4*sizeof(Object*));
	obs[0] = s1;
	obs[1] = s2;
	obs[2] = s3;
	obs[3] = f;
	screen->render(1500,1000,eye,light,obs,4);
}
