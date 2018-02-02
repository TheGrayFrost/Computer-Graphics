# include <fstream>
# include <iostream> 
# include <math.h>

# define MAX 255  // MAX pixel intensity

# define fine 200 // no of pixels in one unit of coordinate axes (e.g. between 0 and 1 on coordinate axes)
                  // controls finesse of the image

# define bmargin 10 // border size

# define visi MAX // visible edge color
# define invis 63 // invisible edge color
# define bg 0     // background color
# define bcol 0   // border color

using namespace std;

// takes frame of dimensions size x size
// and draws line joining (x1, y1) to (x2, y2)
// using color col
void drawline (unsigned char ** frame, int sizex, int sizey, double x1, double y1, double x2, double y2, int col)
{
  double lf, hx, hy, ux, uy;
  int wx, wy, u;

  // length of line
  lf = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) * fine;

  // starting coordinates of line
  hx = (sizex/2) + (x1 * fine);
  hy = (sizey/2) - (y1 * fine);

  // measure of the length of the line
  ux = ((x2 - x1) * fine) / lf;
  uy = ((y2 - y1) * fine) / lf;

  // drawing the line on the drawn coordinate frame
  for (u = 0; u <= lf; ++u)
  {
    wx = hx;
    wy = hy;
    frame[wy][wx] = col;
    hx += ux;
    hy -= uy;
  }
}


int main()
{
  ofstream fp;
  int cx, cy, cz; // cube centre
  int vx, vy, vz; // viewer coords
  int a;          // cube length
  double x[2], y[2], z[2];  // cube vertices
  double hx[8], hy[8];      // image coords
  int vis[12];  // edge visiblibity

  // user input
  cout << "Please enter length of the cube:\n";
  cin >> a;
  cout << "Please enter the coordinates of the centre of the cube:\n";
  cin >> cx >> cy >> cz; 
  cout << "Please enter the coordinates of the camera:\n";
  cin >> vx >> vy >> vz;

  // finding coordinates of vertices of cube
  double r = a/2.0;
  x[0] = cx - r;
  x[1] = cx + r;
  y[0] = cy - r;
  y[1] = cy + r;
  z[0] = cz - r;
  z[1] = cz + r;
  
  // finding coordinates of projection of cube vertices on image plane
  int i = 0, j, k, l;
  for (j = 0; j < 2; ++j)
  {
    for (k = 0; k < 2; ++k)
    {
      for (l = 0; l < 2; ++l)
      {
        hx[i] = ((x[j] * vz - vx * z[l]) / (vz - z[l]));
        hy[i] = ((y[k] * vz - vy * z[l]) / (vz - z[l]));
        i++;
      }
    }
  }

  // finding image centre
  double xmax, xmin, ymax, ymin, cenx, ceny;
  xmax = xmin = hx[0];
  ymax = ymin = hy[0];
  for (i = 1; i < 8; ++i)
  {
    if (hx[i] > xmax)
      xmax = hx[i];
    if (hx[i] < xmin)
      xmin = hx[i];
    if (hy[i] > ymax)
      ymax = hy[i];
    if (hy[i] < ymin)
      ymin = hy[i];
  }
  cenx = (xmax + xmin) / 2;
  ceny = (ymax + ymin) / 2;

  // shifting origin to centre the image
  for (i = 0; i < 8; ++i)
  {
    hx[i] -= cenx;
    hy[i] -= ceny;
  }
  xmax -= cenx;
  xmin -= cenx;
  ymax -= ceny;
  ymin -= ceny;

  // mapping edges to their endpoints
  int ep[12][2] = {{1, 3}, {1, 0}, {0, 2}, {2, 3}, {7, 5}, {5, 4}, {4, 6}, {6, 7}, {7, 3}, {2, 6}, {1, 5}, {0, 4}};

  // setting all edges to invisible
  for (i = 0; i < 12; ++i)
    vis[i] = 0;

  // finding edge visibility
  // making front face visible
  vis[0] = vis[4] = vis[8] = vis[10] = 1;

  if (vx < x[0])  // left face visible
    vis[0] = vis[1] = vis[2] = vis[3] = 1;
  else if (vx > x[1]) // right face visible
    vis[4] = vis[5] = vis[6] = vis[7] = 1;

  if (vy < y[0])  // bottom face visible
    vis[1] = vis[10] = vis[5] = vis[11] = 1;
  else if (vy > y[1]) // top face visible
    vis[3] = vis[8] = vis[7] = vis[9] = 1;

  // frame size
  int sizex = (xmax - xmin) * fine + (2 * bmargin);
  int sizey = (ymax - ymin) * fine + (2 * bmargin);

  // creating frame and border
  unsigned char ** fi = new unsigned char * [sizey];
  for (i = 0; i < sizey; ++i)
  {
    fi[i] = new unsigned char [sizex];
    for (j = 0; j < sizex; ++j)
      fi[i][j] = bcol;
  }

  // painting background
  for (i = sizey - bmargin; i >= bmargin; --i)
    for (j = sizex - bmargin; j >= bmargin; --j)
      fi[i][j] = bg;

  // drawing invisible edges
  for (i = 0; i < 12; ++i)
    if (vis[i] == 0)
      drawline (fi, sizex, sizey, hx[ep[i][0]], hy[ep[i][0]], hx[ep[i][1]], hy[ep[i][1]], invis);

  // drawing visible edges
  for (i = 0; i < 12; ++i)
    if (vis[i] == 1)
      drawline (fi, sizex, sizey, hx[ep[i][0]], hy[ep[i][0]], hx[ep[i][1]], hy[ep[i][1]], visi);

  // creating file
  fp.open ("../Sample Images/cube.pgm", ios::out); 

  // file header
  fp << "P5\n" << sizex << " " << sizey << "\n" << MAX << "\n";
    
  // writing the pixel intensities to the file
  for (i = 0; i < sizey; ++i)
    for (j = 0; j < sizex; ++j)
      fp << fi[i][j];

  // closing the file
  fp.close();

  return 0;
}