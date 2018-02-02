# include <fstream>
# include <iostream> 
# include <math.h>

# define MAX 255  // MAX pixel intensity

# define fine 100 // no of pixels in one unit of coordinate axes (e.g. between 0 and 1 on coordinate axes)
                  // controls finesse of the image

# define scale 10 // scales the line if it has any co-ordinate more than scale
                  // e.g. for scale = 10, (100, -200) to (-100, 100) would get scaled to (5, -10) to (-5, 5)
                  // but (1, -2) to (-1, -1) would remain unchanged

# define axisthick 11    // axis thickness

# define borderthick 11  // border thickness

# define axisr 106  // axis red
# define axisg 90  // axis blue
# define axisb 205  // axis green

# define borderr 0  // border red
# define borderg 0  // border blue
# define borderb 0  // border green


using namespace std;


// structure to store pixel color
struct pix
{
  unsigned char r;
  unsigned char g;
  unsigned char b;
};


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
  int b, max, size, r, i, j, x, cx, cy, hx, hy, wx, wy, px1, py1, px2, py2, bgR, bgG, bgB, lR, lG, lB, t, progt; 
  double x1, y1, x2, y2, px, py, lf, u, l, progtx, progty;

  // user input
  cout << "Please enter the coordinates through which you want to draw line:\n";
  cin >> px1 >> py1 >> px2 >> py2; 
  cout << "Please enter border size (x " << fine << "px):\n";
  cin >> b;
  cout << "Please enter background color (in RGB format, with max intensity " << MAX << "):\n";
  cin >> bgR >> bgG >> bgB;
  cout << "Please enter line color (in RGB format, with max intensity " << MAX << "):\n";
  cin >> lR >> lG >> lB;
  cout << "Please enter line thickness (in px):\n";
  cin >> t;

  // creating file
  fp.open ("line.ppm", ios::out); 
  
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
  fp << "P6\n" << size << " " << size << "\n" << MAX << "\n";
  
  // white background
  pix ** fi = new pix * [size];
  for (i = 0; i < size; ++i)
  {
    fi[i] = new pix [size];
    for (j = 0; j < size; ++j)
    {
      fi[i][j].r = bgR;
      fi[i][j].g = bgG;
      fi[i][j].b = bgB;
    }
  }

  // drawing axes axisthick px wide and axiscolor colored
  // we are drawing axes 100px longer than the line 
  // say, u = scale + 1
  // x-axis: (-u, u)
  // y-axis: (-u, u)
  r = (b + 1) * fine;
  progt = (axisthick-1)/2;
  for (x = 0; x < size - (2 * r); ++x)
  {
    for (i = -progt; i <= progt; ++i)
    {
      fi[(size/2)+i][r+x].r = axisr;
      fi[(size/2)+i][r+x].g = axisg;
      fi[(size/2)+i][r+x].b = axisb;

      fi[r+x][(size/2)+i].r = axisr;
      fi[r+x][(size/2)+i].g = axisg;
      fi[r+x][(size/2)+i].b = axisb;    
    }   
  }

  // drawing borders borderthick px wide and bordercolor colored 
  // at b * 100px distance from axes
  r = fine;
  progt = (borderthick-1)/2;
  for (x = 0; x < size - (2 * r); ++x)
  {
    for (i = -progt; i <= progt; ++i)
    {
      fi[r+i][r+x].r = borderr;
      fi[r+i][r+x].g = borderg;
      fi[r+i][r+x].b = borderb;

      fi[r+x][r+i].r = borderr;
      fi[r+x][r+i].g = borderg;
      fi[r+x][r+i].b = borderb;

      fi[size-r+i][r+x].r = borderr;
      fi[size-r+i][r+x].g = borderg;
      fi[size-r+i][r+x].b = borderb;

      fi[r+x][size-r+i].r = borderr;
      fi[r+x][size-r+i].g = borderg;
      fi[r+x][size-r+i].b = borderb;
    }
  }

  // finding how fine the line should be
  // and a measure of the slope of the line
  l = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
  lf = l * fine;
  hx = cx = (size/2) + (x1 * fine);
  hy = cy = (size/2) - (y1 * fine);
  px = ((x2 - x1) * fine) / lf;
  py = ((y2 - y1) * fine) / lf;

  // drawing the line on the drawn coordinate frame
  // with t thickness
  progtx = ((y2 - y1)/l)/2;
  progty = ((x1 - x2)/l)/2;
  for (u = 0; u <= lf; ++u)
  {
    for (i = -t; i <= t; ++i)
    {
      wx = hx + i * progtx;
      wy = hy - i * progty;

      fi[wy][wx].r = lR;
      fi[wy][wx].g = lG;
      fi[wy][wx].b = lB;
    }
    hx = cx + u * px;
    hy = cy - u * py;
  }
  
  // writing the pixel intensities to the file
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      fp << fi[i][j].r << fi[i][j].g << fi[i][j].b;

  // closing the file
  fp.close();

  return 0;
}
