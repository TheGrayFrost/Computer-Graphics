# include <bits/stdc++.h>

# define MAX 255	// MAX pixel intensity

# define fine 100	// no of pixels in one unit of coordinate axes (e.g. between 0 and 1 on coordinate axes)
					// controls finesse of the image

# define bmargin 10 // border size

# define visi MAX	// visible edge color
# define invisi 63	// invisible edge color
# define bg 0		// background color
# define bcol 0	 	// border color

# define SX 460		// frame size


using namespace std;

int abs(int a)
{
	if (a >= 0)
		return a;
	return -a;
}

// takes frame of dimensions size x size
// and draws line joining (x1, y1) to (x2, y2)
// using color col
void drawline (unsigned char ** frame, int sizex, int sizey, double x1, double y1, double x2, double y2, unsigned char col)
{
	int wx1, wy1, wx2, wy2;
	int hx1, hx2, hy1, hy2;
	int d1, d2, flag, h, hn;

	// points in real 2D space to pixel coordinates
	hx1 = (sizex/2) + (x1 * fine);
	hy1 = (sizey/2) - (y1 * fine);
	hx2 = (sizex/2) + (x2 * fine);
	hy2 = (sizey/2) - (y2 * fine);

	// starting coordinates of line: 
	// (wx1, wy1) -> (wx2, wy2) s.t. wx1 <= wx2
	if (hx1 <= hx2)
	{
		wx1 = hx1;
		wy1 = hy1;
		wx2 = hx2;
		wy2 = hy2;
	}
	else
	{
		wx2 = hx1;
		wy2 = hy1;
		wx1 = hx2;
		wy1 = hy2;
	}

	// determining line orientation
	d1 = wx2 - wx1;
	d2 = wy1 - wy2;
	if ((d1 >= d2) && (d2 >= 0))
		flag = 0;
	else if (d2 > d1)
		flag = 1;	
	else if (-d1 > d2)
		flag = 2;
	else
		flag = 3;

	// drawing based on line orientation
	frame[wy1][wx1] = col;
    h = 0;
    //d2 = abs(d2);
	switch (flag)
	{
		case 0: while (wx1 < wx2)
				{
					wx1 += 1;
					hn = h - d1 + d2;
					h = h + d2;
					if (abs(h) > abs(hn))
					{
						wy1 -= 1;
						h = hn;
					}
					frame[wy1][wx1] = col; 
				}
                break;

        case 1: while (wy1 > wy2)
                {
                    wy1 -= 1;
                    hn = h - d1 + d2;
                    h = h - d1;
                    if (abs(h) > abs(hn))
                    {
                        wx1 += 1;
                        h = hn;
                    }
                    frame[wy1][wx1] = col; 
                }
                break;

        case 2: while (wy1 < wy2)
                {
                    wy1 += 1;
                    hn = h + d1 + d2;
                    h = h + d1;
                    if (abs(h) > abs(hn))
                    {
                        wx1 += 1;
                        h = hn;
                    }
                    frame[wy1][wx1] = col; 
                }
                break;
                
        case 3: while (wx1 < wx2)
                {
                    wx1 += 1;
                    hn = h + d1 + d2;
                    h = h + d2;
                    if (abs(h) > abs(hn))
                    {
                        wy1 += 1;
                        h = hn;
                    }
                    frame[wy1][wx1] = col; 
                }
                break;
	}	
}

// mapping edges to their endpoints
int ev[12][2] = {{1, 3}, {1, 0}, {0, 2}, {2, 3}, {7, 5}, {5, 4}, {4, 6}, {6, 7}, {7, 3}, {2, 6}, {1, 5}, {0, 4}};
// mapping faces to vertices
int fv[6][4] = {{1, 5, 7, 3}, {0, 4, 6, 2}, {1, 0, 2, 3}, {5, 4, 6, 7}, {2, 6, 7, 3}, {1, 5, 4, 0}};
// mapping faces to edges
int fe[6][4] = {{0, 10, 4, 8}, {2, 11, 6, 9}, {0, 1, 2, 3}, {4, 5, 6, 7}, {3, 9, 7, 8}, {1, 10, 5, 11}};

// structure to operate on 3D points and vectors
struct D3point
{
	double x, y, z;

	D3point()	{x = y = z = 0;}

	D3point(double a, double b, double c)
	{
		x = a; y = b; z = c;
	}

	D3point operator + (D3point u)
	{
		return (D3point(x + u.x, y + u.y, z + u.z));
	}

	D3point operator - (D3point u)
	{
		return (D3point(x - u.x, y - u.y, z - u.z));
	}

	D3point operator / (double u)
	{
		return (D3point(x/u, y/u, z/u));
	}

	double operator * (D3point u)
	{
		return (x * u.x + y * u.y + z * u.z);
	}
};

// Rotates a 3D point about z-axis by thz, y-axis by thy, x-axis by thx (in that order i.e. z then y then x)
D3point rotator (D3point ori, double thx, double thy, double thz)
{
	D3point cha;
	cha.x = ori.x * cos(thy) * cos(thz) - ori.y * cos(thy) * sin(thz) + ori.z * sin(thy);
	cha.y = ori.x * (sin(thx) * sin(thy) * cos(thz) + cos(thx) * sin(thz)) + ori.y * (cos(thx) * cos (thz) - sin(thx) * sin(thy) * sin(thz)) - ori.z * sin(thx) * cos(thy);
	cha.z = ori.x * (sin(thx) * sin(thz) - cos(thx) * sin(thy) * cos(thz)) + ori.y * (cos(thx) * sin(thy) * sin(thz) + sin(thx) * cos (thz)) + ori.z * cos(thx) * cos(thy);
	return cha;
}

int main()
{
	int cx, cy, cz; 	// cube centre
	int vx, vy, vz; 	// viewer coords
	int a;				// cube length
	D3point ch[8], v[8], vi, norm, nv;	// cube vertices
	double x[2], hx[8], hy[8], dot;		// image coords
	int i, j, k, l, rot;
	double rx, ry, rz;		// rotation angles	
	string fn;		// file created
	int vis[12];	// edge visiblibity
	double zmax, zmin;
	unsigned char ** fi;	// image frame

	// user input
	cout << "Please enter length of the cube:\n";
	cin >> a;
	cout << "Please enter the coordinates of the centre of the cube:\n";
	cin >> cx >> cy >> cz; 
	cout << "Please enter the coordinates of the camera:\n";
	cin >> vx >> vy >> vz;

	// finding coordinates of vertices of cube w.r.t. centre
	double r = a/2.0;
	x[0] = -r;
	x[1] = r;
	
	i = 0;
	for (j = 0; j < 2; ++j)
		for (k = 0; k < 2; ++k)
			for (l = 0; l < 2; ++l)
				v[i++] = D3point(x[j], x[k], x[l]);

	// view vector
	vi = D3point(cx-vx, cy-vy, cz-vz);


	// creating frame
	fi = new unsigned char * [SX];
	for (i = 0; i < SX; ++i)
		fi[i] = new unsigned char [SX];

	for (rot = 0; rot < 360; ++rot)
	{
		// rotating the cube
		rx = (rot * M_PI) / 180; ry = 2 * rx; rz = 3 * rx;
		for (i = 0; i < 8; ++i)
			ch[i] = rotator(v[i], rx, ry, rz);

		// initialising all edges as invisible
		for (i = 0; i < 12; ++i)
			vis[i] = 0;

		//cout << rot << "\n";

		// finding edge visibility
		for (i = 0; i < 6; ++i)
		{
			// finding normal of face i
			norm = D3point(0, 0, 0);
			for (j = 0; j < 4; ++j)
				norm = norm + ch[fv[i][j]];
			norm = norm / 4;
			nv = norm - vi;

			// finding dot vith view vector
			dot = norm * nv;

			// setting visibility
			if (dot < 0)
				for (j = 0; j < 4; ++j)
					vis[fe[i][j]] = 1;
		}
	

		// finding coordinates of projection of cube vertices on image plane
		for (i = 0; i < 8; ++i)
		{
		    hx[i] = ((ch[i].x * vi.z - vi.x * ch[i].z) / (vi.z - ch[i].z));
		    hy[i] = ((ch[i].y * vi.z - vi.y * ch[i].z) / (vi.z - ch[i].z));
		}
        

		// scaling image
		zmax = zmin = hx[0];
		for (i = 1; i < 8; ++i)
		{
			if (hx[i] > zmax)
				zmax = hx[i];
			if (hx[i] < zmin)
				zmin = hx[i];
			if (hy[i] > zmax)
				zmax = hy[i];
			if (hy[i] < zmin)
				zmin = hy[i];
		}
		if (abs(zmin) > zmax)
			zmax = -zmin;
		if (zmax > 2)
		{
			for (i = 0; i < 8; ++i)
			{
				hx[i] = (hx[i]*2)/zmax;
				hy[i] = (hy[i]*2)/zmax;
			}
		}

		// creating frame and border
		for (i = 0; i < SX; ++i)
			for (j = 0; j < SX; ++j)
				fi[i][j] = bcol;

		// painting background
		for (i = SX - bmargin; i >= bmargin; --i)
			for (j = SX - bmargin; j >= bmargin; --j)
				fi[i][j] = bg;

		// drawing invisible edges
		for (i = 0; i < 12; ++i)
			if (vis[i] == 0)
				drawline (fi, SX, SX, hx[ev[i][0]], hy[ev[i][0]], hx[ev[i][1]], hy[ev[i][1]], invisi);

		// drawing visible edges
		for (i = 0; i < 12; ++i)
			if (vis[i] == 1)
				drawline (fi, SX, SX, hx[ev[i][0]], hy[ev[i][0]], hx[ev[i][1]], hy[ev[i][1]], visi);

		// creating file
		ofstream fp;
		string u = to_string(rot+1);
		if (u.size() == 1)
			u = "00" + u;
		else if (u.size() == 2)
			u = "0" + u;
		fn = "Sample Images/cube_" + u + ".pgm";
		fp.open (fn.c_str(), ios::out); 

		// file header
		fp << "P5\n" << SX << " " << SX << "\n" << MAX << "\n";
			
		// writing the pixel intensities to the file
		for (i = 0; i < SX; ++i)
			for (j = 0; j < SX; ++j)
				fp << fi[i][j];

		// closing the file
		fp.close();
	}

	return 0;
}