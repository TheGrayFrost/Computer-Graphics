# include <fstream>
# include <iostream> 
# include <math.h>

# define MAX 255  // MAX pixel intensity

# define fine 100 // no of pixels in one unit of coordinate axes (e.g. between 0 and 1 on coordinate axes)
                  // controls finesse of the image

# define scale 10 // scales the circle if it has any co-ordinate more than scale
                  // e.g. for scale = 10, 
                  // circle of radius 500 centred at (100, -500) 
                  // will have extreme ends at (600, -500), (-400, -500), (100, 0), (100, -1000)
                  // -1000 being the largest ordinate, it would get rescaled so that 1000 becomes 10
                  // and thus would become a 5 radius circle centred at (1, -5)
                  // but a circle of radius 3 centred at (1, -2) would remain unchanged

using namespace std;

// absolute value
int abs(int u)
{
  if (u > 0)
    return u;
  return -u;
}

// maximum of absolute values
int absmax (int g, int h)
{
  int a = abs(g), b = abs(h);
  if (a > b)
    return a;
  return b;
}

int main()
{
  ofstream fp;
  int b, max, size, r, i, j, x, cx, cy, rad, ux, uy, wx, wy, nt; 
  double radf, circum, ang, px, py;

  // user input
  cout << "Please enter the coordinates of centre of circle:\n";
  cin >> cx >> cy; 
  cout << "Please radius of circle:\n";
  cin >> rad;
  cout << "Please enter border size (x " << fine << "px):\n";
  cin >> b;

  // creating file
  fp.open ("circle.pgm", ios::out); 
  
  // scaling largest absolute value coordinate to scale
  // only if the lagest ordinate is more than scale
  max = absmax(absmax(cx+rad, cx-rad), absmax(cy+rad, cy-rad));
  if (max > scale)
  {
    px = (double(cx)/max)*scale;
    py = (double(cy)/max)*scale;     
  }
  else
  {
  	px = cx;
    py = cy;
  }


  size = 2 * (scale + 2 + b) * fine;

  // file header
  fp << "P5\n" << size << " " << size << "\n" << MAX << "\n";
  
  // white background
  unsigned char ** fi = new unsigned char * [size];
  for (i = 0; i < size; ++i)
  {
    fi[i] = new unsigned char [size];
    for (j = 0; j < size; ++j)
      fi[i][j] = MAX;
  }

  // drawing axes
  // we are drawing axes 100px longer than the line 
  // say, u = scale + 1
  // x-axis: (-u, u)
  // y-axis: (-u, u)
  r = (b + 1) * fine;
  for (x = 0; x < size - (2 * r); ++x)
  {
    fi[size/2][r+x] = 150;
    fi[r+x][size/2] = 150;
  }

  // drawing borders at b * 100px distance from axes
  r = fine;
  for (x = 0; x < size - (2 * r); ++x)
  {
    fi[r][r+x] = 0;
    fi[r+x][r] = 0;
    fi[size-r][r+x] = 0;
    fi[r+x][size-r] = 0;
  }

  // finding coordinates of points on circle
  // and coloring them black
  nt = rad * fine;
  circum = 2 * M_PI * nt;
  ux = size/2 + px * fine;
  uy = size/2 - py * fine;
  for (i = 0; i <= circum; ++i)
  {
    ang = double(i)/nt;
    wx = ux + nt * cos(ang);
    wy = uy + nt * sin(ang);
    fi[wy][wx] = 0;
  }
  
  // writing the pixel intensities to the file
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      fp << fi[i][j];

  // closing the file
  fp.close();

  return 0;
}