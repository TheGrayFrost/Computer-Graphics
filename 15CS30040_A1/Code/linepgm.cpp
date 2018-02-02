# include <fstream>
# include <iostream> 
# include <math.h>

# define MAX 255  // MAX pixel intensity

# define fine 100 // no of pixels in one unit of coordinate axes (e.g. between 0 and 1 on coordinate axes)
                  // controls finesse of the image

# define scale 10 // scales the line if it has any co-ordinate more than scale
                  // e.g. for scale = 10, (100, -200) to (-100, 100) would get scaled to (5, -10) to (-5, 5)
                  // but (1, -2) to (-1, -1) would remain unchanged

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
  int b, max, size, r, i, j, x, cx, cy, hx, hy, px1, py1, px2, py2; 
  double x1, y1, x2, y2, px, py, lf, u;

  // user input
  cout << "Please enter the coordinates through which you want to draw line:\n";
  cin >> px1 >> py1 >> px2 >> py2; 
  cout << "Please enter border size (x " << fine << "px):\n";
  cin >> b;

  // creating file
  fp.open ("line.pgm", ios::out); 
  
  // scaling largest absolute value coordinate to scale
  // only if the lagest ordinate is more than scale
  max = absmax(absmax(px1, py1), absmax(px2, py2));
  if (max > scale)
  {
  	x1 = (double(px1) / max) * scale;
    y1 = (double(py1) / max) * scale;
    x2 = (double(px2) / max) * scale;
    y2 = (double(py2) / max) * scale;
  }
  else
  {
  	x1 = px1;
  	y1 = py1;
  	x2 = px2;
  	y2 = py2;
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

  // finding how fine the line should be
  // and a measure of the slope of the line
  lf = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) * fine;
  hx = cx = (size/2) + (x1 * fine);
  hy = cy = (size/2) - (y1 * fine);
  px = ((x2 - x1) * fine) / lf;
  py = ((y2 - y1) * fine) / lf;

  // drawing the line on the drawn coordinate frame
  for (u = 0; u <= lf; ++u)
  {
    fi[hy][hx] = 0;
    hx = cx + u * px;
    hy = cy - u * py;
  }
  
  // writing the pixel intensities to the file
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      fp << fi[i][j];

  // closing the file
  fp.close();

  return 0;
}