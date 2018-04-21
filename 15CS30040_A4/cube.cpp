#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#define MAX 255					//MAX pixel intensity
#define fine 100				//no of pixels in one unit of coordinate axes (e.g. between 0 and 1 on coordinate axes)
								//controls finesse of the image
#define SX 401					//frame size
#define CMAX (SX / (2 * fine))  //coordinate max fitting into frame
#define MAXROT 360				//number of images to generate
#define zEps 0


using namespace std;

double mind (double a, double b) {return (a < b ? a : b);}	//double min
double maxd (double a, double b) {return (a > b ? a : b);}	//double max
int abs(int a) {return (a >= 0 ? a : -a);}					//int abs
int minin (int a, int b) {return (a < b ? a : b);}			//int min
int maxin (int a, int b) {return (a > b ? a : b);}			//int max


//structure to operate on 3D points and vectors
struct D3point
{
	double x, y, z;	 //3D coordinates (double as vector components)

	D3point() {x = y = z = 0;}
	D3point(double a, double b, double c) {x = a; y = b; z = c;}
	D3point operator + (D3point const u) const {return (D3point(x + u.x, y + u.y, z + u.z));}
	D3point operator - (D3point const u) const {return (D3point(x - u.x, y - u.y, z - u.z));}
	D3point operator / (double const u) const {return (D3point(x/u, y/u, z/u));}			//scale
	D3point operator * (double const u) const {return (D3point(x*u, y*u, z*u));}			//scale
	double operator * (D3point const u) const {return (x * u.x + y * u.y + z * u.z);}		//dot
	double norm() const {return (sqrt(x*x + y*y + z*z));}
	D3point unify() const {return (*this / this->norm());}
};

D3point reject(0, 0, 2);

//attenuation characteristics of environment
struct fatt
{
	double d1, d2;  //refer atten function for meaning

	fatt() {d1 = d2 = 0;}
	double atten(D3point r) {double d = r.norm(); return (1 / (1 + d1*d + d2*d*d));}
} envir;


//illumination characteristics of object
struct illchar
{
	double kd;	  //diffused reflection coefficient
	double ks;	  //specular reflection coefficient
	double specexp; //specular reflection exponent

	illchar() {kd = ks = specexp = 0;}
	illchar(double a, double b, double c) {kd = a; ks = b; specexp = c;}
};


//structure for HSV color
struct HSV
{
	double h, s, v;
	HSV() {h = s = v = 0;}
	HSV(double a, double b, double c) {h = a; s = b; v = c;}
	HSV operator * (double u) const {return (HSV(h, s, v*u));}  //scale
};


//structure for RGB color
struct RGB
{
	int R, G, B;

	RGB() {R = G = B = 0;}
	RGB(int a, int b, int c) {R = a; G = b; B = c;}
	RGB(HSV);																		   //typecast HSV->RGB
	RGB operator + (RGB u) const {return (RGB(R + u.R, G + u.G, B + u.B));}
	RGB& operator += (RGB u) {*this = *this + u; return *this;}
	RGB operator / (double u) const {return (RGB(R/u, G/u, B/u));}					  //final scaling
	RGB operator * (RGB u) const {return (RGB(R * u.R, G * u.G, B * u.B) / 255);}	   //colour absorption
	RGB operator * (double u) const {return (RGB(HSV(*this)*u));}					   //HSV brightness scaling
	operator HSV() const;															   //typecast RGB->HSV
	RGB max(RGB u, RGB v) {return (RGB(maxin(u.R, v.R), maxin(u.B, v.B), maxin(u.B, v.B)));}
} amb, src;


//typecast RGB->HSV
RGB::operator HSV () const
{
	double h, s, v;
	int imin, imax, delta;

	imin = minin(R, minin(G, B));
	imax = maxin(R, maxin(G, B));
	delta = imax - imin;
	v = ((double)imax) / MAX;

	if (delta == 0)
	{
		s = 0;
		h = 0;
	}
	else
	{
		if (imax > 0)
			s = (double)delta / imax;
		else
			s = 0;

		if (imax == R)
			h = (double)(G - B) / delta;
		else if (imax == G)
			h = 2 + (double)(B - R) / delta;
		else
			h = 4 + (double)(R - G) / delta;
		h *= 60;

		if (h < 0.0)
			h += 360;
	}
	return HSV(h, s, v);
}


//typecast HSV->RGB
RGB::RGB (HSV in)
{
	if (in.s <= 0.0)
	{
		R = in.v * MAX;
		G = in.v * MAX;
		B = in.v * MAX;
	}
	else
	{
		double hh = in.h;
		if (hh >= 360.0)
			hh -= 360.0;
		hh /= 60;
		long i = (long)hh;
		double ff = hh - i;
		int p = in.v * (1.0 - in.s) * MAX;
		int q = in.v * (1.0 - (in.s * ff)) * MAX;
		int r = in.v * MAX;
		int t = in.v * (1.0 - (in.s * (1.0 - ff))) * MAX;
		switch(i)
		{
			case 0: R = r; G = t; B = p; break;
			case 1: R = q; G = r; B = p; break;
			case 2: R = p; G = r; B = t; break;
			case 3: R = p; G = q; B = r; break;
			case 4: R = t; G = p; B = r; break;
			default: R = r; G = p; B = q; break;
		}
	}
}

//true and easy scaling : premultiplication
RGB operator * (double u, RGB col) {return RGB(u * col.R, u * col.G, u * col.B);}

//real points in 2d view plane
struct px
{
	double x, y;
	px (double a = 0, double b = 0) {x = a; y = b;}
	px operator * (double u) const {return px(x * u, y * u);}   //scale
};


//discrete points on frame
struct fr
{
	int x, y;
	fr (int a = 0, int b = 0) {x = a; y = b;}
	fr (px v) {x = SX/2 + v.x * fine; y = (SX/2) - (v.y * fine);}
	operator px() const {return px((double)(x - (SX/2)) / fine, (double)((SX/2) - y) / fine);}
};


//structure for object faces
struct obj
{
	D3point vtx[4];	//face vertices
	D3point cen;	//face centre
	D3point nml;	//normal to face
	RGB amb;		//ambient light color of face
	RGB diff;		//diffuse light color of face
	RGB spec;		//specular light color of face
	double spep;	//specular reflection exponent
};


//structure for frame points
struct fpt
{
	RGB col;		//color of the point
	double z;		//Z-buffer

	fpt() {z = INFINITY;}
};


//print 3D point
ostream & operator << (ostream & s, D3point const & u) {s << u.x << " " << u.y << " " << u.z << "\n"; return s;}
//display RGB element
ostream & operator << (ostream & s, RGB const & u) {s << u.R << " " << u.G << " " << u.B << "\n"; return s;}
//display HSV element
ostream & operator << (ostream & s, HSV const & u) {s << u.h << " " << u.s << " " << u.v << "\n"; return s;}
//display px element
ostream & operator << (ostream & s, px const & u) {s << u.x << " " << u.y << "\n"; return s;}
//display fr element
ostream & operator << (ostream & s, fr const & u) {s << u.x << " " << u.y << "\n"; return s;}


//Rotates a 3D point about z-axis by thz, y-axis by thy, x-axis by thx (in that order i.e. z then y then x)
D3point rotator (D3point ori, double thx, double thy, double thz)
{
	D3point cha;
	cha.x = ori.x * cos(thy) * cos(thz) - ori.y * cos(thy) * sin(thz) + ori.z * sin(thy);
	cha.y = ori.x * (sin(thx) * sin(thy) * cos(thz) + cos(thx) * sin(thz)) +
			ori.y * (cos(thx) * cos (thz) - sin(thx) * sin(thy) * sin(thz)) -
			ori.z * sin(thx) * cos(thy);
	cha.z = ori.x * (sin(thx) * sin(thz) - cos(thx) * sin(thy) * cos(thz)) +
			ori.y * (cos(thx) * sin(thy) * sin(thz) + sin(thx) * cos (thz)) +
			ori.z * cos(thx) * cos(thy);
	return cha;
}

//mapping faces to vertices
int fv[6][4] = {{1, 5, 7, 3}, {0, 4, 6, 2}, {1, 0, 2, 3}, {5, 4, 6, 7}, {2, 6, 7, 3}, {1, 5, 4, 0}};


//finds px image of vertex wrt to viewer vi
//screen assumed at z = 0
//scales down by scale
px primage (D3point vi, D3point p, double scale = 1.0)
{
	double r = ((p.x * vi.z - vi.x * p.z) / (vi.z - p.z));
	double s = ((p.y * vi.z - vi.y * p.z) / (vi.z - p.z));
	return px(r, s) * scale;
}

//checks if a ray intersects a face, returns reject if not
//otherwise returns where the ray produced intersects the face
//obcheck tells whether we want the fact to obstruct the ray or not
D3point project (D3point vi, D3point p, obj * obc, bool obcheck = false)
{
	double test = (vi - p) * obc->nml;
	if (abs(test) <= 1e-5)
	{
		//Line parallel to plane.
		return reject;
	}
	double lamb = ((vi - obc->cen) * obc->nml) / test;

	if ((lamb <= 0) || (lamb >= 1 && obcheck))
	{
		//Line on one side of plane.
		return reject;
	}
	D3point ip = vi + (p - vi) * lamb, x, y, z, ycap, xcap;
	z = ip - obc->vtx[0];
	x = obc->vtx[1] - obc->vtx[0];
	y = obc->vtx[3] - obc->vtx[0];
	ycap = y.unify();
	xcap = x.unify();
	double alpha = z*xcap/x.norm(), beta = z*ycap/y.norm();

	if ((alpha <= zEps) || (alpha >= 1) || (beta <= zEps) || (beta >= 1))
	{
		//Line passes outside of face boundaries.
		return reject;
	}
	return ip;
}


//renders an object using Z-buffer and ray tracing
//illuminates the points
//if obstr is not NULL, it creates shadows of the obstruction as required during illumination
void render (fpt frame[SX][SX], D3point vi, obj obc[], int nobc, double scale, D3point si[], int nsrc, obj obstr[] = NULL, int nobs = 0)
{
	int i, j, k, l, m;
	double dot;
	RGB extra;
	px rv;
	D3point p, prime, hiccup, v, vcap, s, lcap, rcap, ncap;
	bool seeable;

	for (i = 0; i < SX; ++i)
	{
		for (j = 0; j < SX; ++j)
		{
			rv = px(fr(j, i)) * (1.0 / scale);
			p = D3point(rv.x, rv.y, 0);

			for (k = 0; k < nobc; ++k)
			{
				prime = project(vi, p, &obc[k]);
				if (prime.z < 1) //checking if not reject
				{
					if (-prime.z <= frame[i][j].z) //checking if this new point is ahead of previous point
					{
						ncap = obc[k].nml;
						frame[i][j].col = obc[k].amb;
						frame[i][j].z = -prime.z;
						v = vi - prime;
						vcap = v.unify();

						for (l = 0; l < nsrc; ++l)
						{
							seeable = true;
							
							for (m = 0; m < nobs; ++m)
							{
								hiccup = project (si[l], prime, &obstr[m], true);
								if (hiccup.z < 1)	//obstr[m] obstructs the ray
								{
									seeable = false;
									break;
								}
							}

							if (seeable)
							{
								s = si[l] - prime;
								lcap = s.unify();
								dot = lcap * ncap;
                                if (dot > 0)
                                {
                                    extra = dot * obc[k].diff;
                                    rcap = ncap * 2 * dot - lcap;
                                    dot = rcap * vcap;
                                    if (dot > 0)
                                        dot = pow(dot, obc[k].spep);
                                    else
                                        dot = 0;
                                    extra += dot * obc[k].spec;
                                    extra = envir.atten(s) * extra;
                                }
								else
                                    extra = RGB(0, 0, 0);
								frame[i][j].col += extra;
							}
						}
					}
				}
			}
		}
	}
}

void progressBar (double progress)
{
	int barWidth = 70;
	cout << "[";
	int pos = barWidth * progress;
	for (int i = 0; i <= barWidth; ++i)
	{
		if (i < pos) cout << "=";
		else if (i == pos) cout << ">";
		else cout << " ";
	}
	cout << "] " << int(progress * 100.0) << " %\r";
	cout.flush();
}


int main()
{
	//freopen("../input.txt", "r", stdin);

	D3point vi;				//viewer coords
	int nsrc;				//number of sources
	D3point * si;			//source coords

	double rd; 				//room dimension
	D3point rc; double rz;	//room centre
	D3point rv[8];			//room coords
	px r[8];				//room's pixel coords
	obj room[6];			//room's visible walls
	RGB dr[6], sr;			//room colors
	illchar ric;			//room illumination characteristics

	double cd;				//cube dimension
	D3point cc;				//cube centre
	D3point cv[8], mv[8];	//cube vertices
	int vis[6];				//cube face visiblibity
	obj cube[6];			//cube's faces
	RGB dc, sc;				//cube colors
	px c[8];				//cube's pixel coords
	illchar cic;			//cube illumination characteristics

	fpt fi[SX][SX];			//frame

	//local support
	string fn;
	ofstream fp;
	int i, j, k, l, colmax, rot;
	double x[2], zmax, zmin, dot, rotx, roty, rotz, scale;
	RGB curcol;
	D3point nv;
	unsigned char z = 0;

	//user input
	cout << "Please enter the coordinates of the viewer:\n";
	cin >> vi.x >> vi.y >> vi.z;
	cout << "Please enter color of ambient light:\n";
	cin >> amb.R >> amb.G >> amb.B;
	cout << "Please enter attenuation coefficients (d1 and d2):\n";
	cin >> envir.d1 >> envir.d2;

	cout << "\nPlease enter number of light sources:\n";
	cin >> nsrc;
	cout << "Please enter color of the light sources:\n";
	cin >> src.R >> src.G >> src.B;
	si = new D3point [nsrc];
	for (i = 0; i < nsrc; ++i)
	{
		cout << "\nLight source " << i+1 << " \n";
		cout << "\nPlease enter the coordinates of the light source:\n";
		cin >> si[i].x >> si[i].y >> si[i].z;
	}

	cout << "\nPlease enter size of the room:\n";
	cin >> rd;
	cout << "Please enter the coordinates of the centre of the room:\n";
	cin >> rc.x >> rc.y >> rz;
	cout << "Please enter diffuse color of room walls (back, left, right, top, bottom):\n";
	for (int i = 0; i < 5; ++i)
		cin >> dr[i].R >> dr[i].G >> dr[i].B;
	cout << "Please enter specular color of room:\n";
	cin >> sr.R >> sr.G >> sr.B;
	cout << "Please enter illumination characteristics of room (kd, ks, m or specular exponent):\n";
	cin >> ric.kd >> ric.ks >> ric.specexp;


	cout << "\nPlease enter size of the cube:\n";
	cin >> cd;
	cout << "Please enter the coordinates of the centre of the cube:\n";
	cin >> cc.x >> cc.y >> cc.z;
	cout << "Please enter diffuse color of cube:\n";
	cin >> dc.R >> dc.G >> dc.B;
	cout << "Please enter specular color of cube:\n";
	cin >> sc.R >> sc.G >> sc.B;
	cout << "Please enter illumination characteristics of cube (kd, ks, m or specular exponent):\n";
	cin >> cic.kd >> cic.ks >> cic.specexp;


	//centring everything to the room centre
	vi = vi - rc;
	for (i = 0; i < nsrc; ++i)
		si[i] = si[i] - rc;
	cc = cc - rc;

	//finding 3D coordinates of room vertices w.r.t. centre
	rd /= 2.0;
	x[0] = -rd;
	x[1] = rd;
	rc = D3point(0, 0, 0);

	i = 0;
	for (j = 0; j < 2; ++j)
	{
		for (k = 0; k < 2; ++k)
		{
			for (l = 0; l < 2; ++l)
			{
				rv[i] = D3point(x[j], x[k], mind(0, x[l] + rz));
				rc = rc + rv[i];
				i++;
			}
		}
	}
	rc = rc/8;

	//finding room face attributes
	for (i = 0; i < 5; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			room[i].vtx[j] = rv[fv[i+1][j]];
			room[i].cen = room[i].cen + room[i].vtx[j];
		}
		room[i].cen = room[i].cen/4;
		room[i].nml = rc - room[i].cen;
		room[i].nml = room[i].nml.unify();
		room[i].amb = ric.kd * (dr[i] * amb);
		room[i].diff = ric.kd * (dr[i] * src);
		room[i].spec = ric.ks * (sr * src);
		room[i].spep = ric.specexp;
	}

	//finding coordinates of projection of room vertices on image plane
	for (i = 0; i < 8; ++i)
		r[i] = primage (vi, rv[i]);

	//finding image scale
	zmax = maxd(r[0].x, r[0].y);
	zmin = mind(r[0].x, r[0].y);
	for (i = 1; i < 8; ++i)
	{
		zmax = maxd(zmax, maxd(r[i].x, r[j].x));
		zmin = mind(zmin, mind(r[i].x, r[j].x));
	}

	if (abs(zmin) > zmax)
		zmax = -zmin;

	//scaling the room
	if (zmax > CMAX)
		scale = (zmax > CMAX ? ((double) CMAX / zmax) : 1.0);

	//finding 3D coordinates of cube vertices w.r.t. centre
	cd /= 2.0;
	x[0] = -cd;
	x[1] = cd;

	i = 0;
	for (j = 0; j < 2; ++j)
		for (k = 0; k < 2; ++k)
			for (l = 0; l < 2; ++l)
				cv[i++] = D3point(x[j], x[k], x[l]);

	for (i = 0; i < 6; ++i)
	{
		cube[i].amb = cic.kd * (dc * amb);
		cube[i].diff = cic.kd * (dc * src);
		cube[i].spec = cic.ks * (sr * src);
		cube[i].spep = cic.specexp;
	}

	cout << "\nRendering Images...\nProgress...\n";
	for (rot = 0; rot < MAXROT; ++rot)
    {
		//rotating the cube
		rotx = (rot * M_PI) / 180;
		roty = 2 * rotx;
		rotz = 3 * rotx;
		for (i = 0; i < 8; ++i)
			mv[i] = rotator(cv[i], rotx, roty, rotz) + cc;

		//finding cube face attributes
		for (i = 0; i < 6; ++i) {
			cube[i].cen = D3point();
			for (j = 0; j < 4; ++j) {
				cube[i].vtx[j] = mv[fv[i][j]];
				cube[i].cen = cube[i].cen + cube[i].vtx[j];
			}
			cube[i].cen = cube[i].cen / 4;
			cube[i].nml = cube[i].cen - cc;
			cube[i].nml = cube[i].nml.unify();
		}

		//finding cube face visibility
		for (i = 0; i < 6; ++i) {
			//finding dot of normal with view vector
			nv = cube[i].cen - vi;
			dot = cube[i].nml * nv;

			//setting visibility
			if (dot < 0)
				vis[i] = 1;
			else
				vis[i] = 0;
		}

		//drawing the room
		render(fi, vi, room, 5, scale, si, nsrc, cube, 6);
		render(fi, vi, cube, 6, scale, si, nsrc);

		int oldcol = colmax = 0;
		for (i = 0; i < SX; ++i) {
			for (j = 0; j < SX; ++j) {
				curcol = fi[i][j].col;
				colmax = maxin(colmax, curcol.R);
				colmax = maxin(colmax, curcol.G);
				colmax = maxin(colmax, curcol.B);
				if (oldcol != colmax)
					oldcol = colmax;
			}
		}

		if (colmax > MAX)
			for (i = 0; i < SX; ++i)
				for (j = 0; j < SX; ++j)
					fi[i][j].col = (MAX * fi[i][j].col) / colmax;

		//creating file
		string u = to_string(rot + 1);
		if (u.size() == 1)
			u = "00" + u;
		else if (u.size() == 2)
			u = "0" + u;
		fn = "Images/cube_" + u + ".ppm";
		fp.open(fn.c_str(), ios::out);

		//file header
		fp << "P6\n" << SX << " " << SX << "\n" << MAX << "\n";

		//writing the pixel intensities to the file
		for (i = 0; i < SX; ++i)
		{
			for (j = 0; j < SX; ++j)
			{
				fp << (unsigned char) fi[i][j].col.R << (unsigned char) fi[i][j].col.G
				   << (unsigned char) fi[i][j].col.B;
				fi[i][j].z = INFINITY;
			}
		}

		//closing the file
		fp.close();
		progressBar((rot + 1) * 1.0 / MAXROT);
	}

	cout << "\nRendering complete.\n\n";

	return 0;
}