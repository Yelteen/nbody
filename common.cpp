#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  quadtree constructor :
//
QuadTreeNode::QuadTreeNode(QuadTreeNode* parent, double x, double y, 
                           double width, double height)
{
  parent   = parent;
  x        = x;
  y        = y;
  width    = width;
  height   = height;
  m        = 0.0;
  external = true;
  NW       = NULL;
  NE       = NULL;
  SW       = NULL;
  SE       = NULL;
  wn       = width/2;        // new quadrant width
  hn       = height/2;       // new quadrant height
  xmid     = x + wn;         // x-midpoint of this quadrant
  ymid     = y + hn;         // y-midpoint of this quadrant
}

//
// insert a new particle into the quadtree :
//
void QuadTreeNode::insert(particle_t* p)
{
  // if this quadtree is external :
  if (external)
  {
    // if this quadtree is empty, put the particle in it :
    if (this->p == NULL)
    {
      this->p  = p;
      this->m += mass;
    }
    // otherwise we need to subdivide and re-insert the particles :
    else
    {
      // subdivide this quadtee :
      NW = new QuadTreeNode(this, x,    y,    wn, hn);
      NE = new QuadTreeNode(this, xmid, y,    wn, hn);
      SW = new QuadTreeNode(this, x,    ymid, wn, hn);
      SE = new QuadTreeNode(this, xmid, ymid, wn, hn);
      
      // it is no longer external :
      external = false;
      
      // re-insert the particle in the correct quadrant :
      this->insert(this->p);
      this->insert(p);
    }
  }
  // else we insert the particles in the appropriate quadrant :
  else
  {
    bool pltx = p->x < xmid;
    bool plty = p->y < ymid;
    
    if (pltx and plty)
    {
      NW->insert(p);
    }
    else if (not pltx and plty)
    {
      NE->insert(p);
    }
    else if (pltx and not plty)
    {
      SW->insert(p);
    }
    else if (not pltx and not plty)
    {
      SE->insert(p);
    }
  }
}

//
//  update a quadtree's center of mass :
//
void QuadTreeNode::computeCOM()
{
  // if this is an external node :
  if (external)
  {
    // and there is a particle in it :
    if (this->p !=NULL)
    {
      this->m     = mass;
      this->com_x = this->p->x;
      this->com_y = this->p->y;
    }
    // otherwise mass is zero :
    else
    {
      this->m     = 0.0;
      this->com_x = 0.0;
      this->com_y = 0.0;
    }
  }
  // otherwise recurse on each quadrant :
  else
  {
    // compute center of mass for each quadrant :
    NW->computeCOM();
    NE->computeCOM();
    SW->computeCOM();
    SE->computeCOM();
    
    // calculate the cetner of mass for this quadrant :
    this->m      = NW->m + NE->m + SW->m + SE->m;
    this->com_x  =   NW->m * NW->com_x
                   + NE->m * NE->com_x
                   + SW->m * SW->com_x
                   + SE->m * SE->com_x;
    this->com_x /= this->m;
    this->com_y  =   NW->m * NW->com_y
                   + NE->m * NE->com_y
                   + SW->m * SW->com_y
                   + SE->m * SE->com_y;
    this->com_y /= this->m;
  }
}

//
// compute the force from this quadrant on a particle :
// 
void QuadTreeNode::computeF(particle_t* p,double* dmin,double* davg,int* navg)
{
  // if this is an external quadtree node :
  if (external)
  {
    // if the quadrant is not empty and the particles being compared 
    // are not the same :
    if (this->p != NULL and this->p != p)
    {
      double dx = this->p->x - p->x;
      double dy = this->p->y - p->y;
      double r2 = dx*dx + dy*dy;
      
      // update the minimum distance between particles :
      if (r2/(cutoff*cutoff) < *dmin * (*dmin))
      {
        *dmin  = sqrt(r2)/cutoff;
      }
      (*davg) += sqrt(r2)/cutoff;
      (*navg) ++;
      
      r2       = fmax( r2, min_r*min_r );
      double r = sqrt( r2 );
      
      //
      //  very simple short-range repulsive force
      //
      double coef = ( 1 - cutoff / r ) / r2 / mass;
      p->ax += coef * dx;
      p->ay += coef * dy;
    }
  }
  // otherwise evaluate the distance to the center of mass :
  else
  {
    double dx = x - p->x;
    double dy = y - p->y;
    double r  = sqrt( dx*dx + dy*dy + 1e-16);
    
    // if the distance is within tolerance, treat quadtree as a single body :
    if (width / r < 0.5)
    {
      // update the minimum distance between particles :
      if (r/cutoff < *dmin)
      {
        *dmin  = r/cutoff;
      }
      (*davg) += r/cutoff;
      (*navg) ++;
      
      r = fmax( r, min_r );
      
      //
      //  very simple short-range repulsive force
      //
      double coef = ( 1 - cutoff / r ) / (r*r*mass);
      p->ax += coef * dx;
      p->ay += coef * dy;
    }
    // otherwise recurse on each quadrant :
    else
    {
      NW->computeF(p, dmin, davg, navg);
      NE->computeF(p, dmin, davg, navg);
      SW->computeF(p, dmin, davg, navg);
      SE->computeF(p, dmin, davg, navg);
    }
  }
}


//
//  timer
//
double read_timer( )
{
  static bool initialized = false;
  static struct timeval start;
  struct timeval end;
  if( !initialized )
  {
    gettimeofday( &start, NULL );
    initialized = true;
  }
  gettimeofday( &end, NULL );
  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
  size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
  srand48( time( NULL ) );
      
  int sx = (int)ceil(sqrt((double)n));
  int sy = (n+sx-1)/sx;
  
  int *shuffle = (int*)malloc( n * sizeof(int) );
  for( int i = 0; i < n; i++ )
  {
    shuffle[i] = i;
  }
  
  for( int i = 0; i < n; i++ ) 
  {
    //
    //  make sure particles are not spatially sorted
    //
    int j = lrand48()%(n-i);
    int k = shuffle[j];
    shuffle[j] = shuffle[n-i-1];
    
    //
    //  distribute particles evenly to ensure proper spacing
    //
    p[i].x = size*(1.+(k%sx))/(1+sx);
    p[i].y = size*(1.+(k/sx))/(1+sy);

    //
    //  assign random velocities within a bound
    //
    p[i].vx = drand48()*2-1;
    p[i].vy = drand48()*2-1;
  }
  free( shuffle );
}

//
//  interact two particles
//
void apply_force(particle_t &particle, particle_t &neighbor,
                 double* dmin, double* davg, int* navg)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  
  if( r2 > cutoff*cutoff )
  {
    return;
  }
  if (r2 != 0)
  {
    if (r2/(cutoff*cutoff) < *dmin * (*dmin))
    {
      *dmin  = sqrt(r2)/cutoff;
    }
    (*davg) += sqrt(r2)/cutoff;
    (*navg) ++;
  }
  
  r2       = fmax( r2, min_r*min_r );
  double r = sqrt( r2 );
  
  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
  //
  //  slightly simplified Velocity Verlet integration
  //  conserves energy better than explicit Euler method
  //
  p.vx += p.ax * dt;
  p.vy += p.ay * dt;
  p.x  += p.vx * dt;
  p.y  += p.vy * dt;

  //
  //  bounce from walls
  //
  while( p.x < 0 || p.x > size )
  {
    p.x  = p.x < 0 ? -p.x : 2*size-p.x;
    p.vx = -p.vx;
  }
  while( p.y < 0 || p.y > size )
  {
    p.y  = p.y < 0 ? -p.y : 2*size-p.y;
    p.vy = -p.vy;
  }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
  static bool first = true;
  if( first )
  {
    fprintf( f, "%d %g\n", n, size );
    first = false;
  }
  for( int i = 0; i < n; i++ )
  {
    fprintf( f, "%g %g\n", p[i].x, p[i].y );
  }
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
  for( int i = 1; i < argc; i++ )
    if( strcmp( argv[i], option ) == 0 )
      return i;
  return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
  int iplace = find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return atoi( argv[iplace+1] );
  return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
  int iplace = find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return argv[iplace+1];
  return default_value;
}



