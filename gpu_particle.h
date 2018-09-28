#ifndef PARTICLE_H
#define PARTICLE_H

#include "types.h"

#include <stdio.h>
#include <string>
#include <string.h>
#include <iostream>
#include <math.h>
#include "run_control.h"

using namespace std;

string FortranExpWrite(double t)
{
      char str[100],res[100];
      char *prev,*next,*dot;
      string res_str;
      
      sprintf(str,"%11.4E",t);
      
      dot = strstr(str,".");
      
      prev = dot - 1;
      next = dot + 1;

      sprintf(res,"0.%c%s",*prev,next);
      
      res_str = res;
      
      return res_str;
}

class GParticle
{
public:  
   double x,y,z,pu,pv,pw,m,q_m;
#ifdef DEBUG_PLASMA
   double3 next_x; 
#endif
   int jmp;
#ifdef ATTRIBUTES_CHECK
   double *d_ctrlParticles;

   void SetControlSystem(int j,double *c){jmp = j;d_ctrlParticles = c;}
#endif

   GParticle(){}
   GParticle(double x1,double y1, double z1,double u1,double v1,double w1,double m1,double q_m1): x(x1), y(y1), z(z1), pu(u1), pv(v1), pw(w1), m(m1), q_m(q_m1) {}
  ~GParticle(){}

#ifdef ATTRIBUTES_CHECK
   void SetControlSystem(int j,double *c){jmp = j;d_ctrlParticles = c;}
#endif

   virtual void Move(double3 E,double3 H,double tau,double *p_control,int jmp_control)
{
        double bx,by,bz,tau1,u,v,w,ps,su,sv,sw,s1,s2,s3,s4,s5,s6,s;
	double sx,sy,sz,x1,y1,z1,pu1,pv1,pw1;
	
	p_control[ParticleAttributePosition(jmp_control,fortran_number,m,1)] = x;

	tau1=q_m*tau*0.5;
    
	pu += tau1*E.x;
	pv += tau1*E.y;
	pw += tau1*E.z;
	ps = tau1 * pow((pu * pu + pv * pv + pw * pw) * 1. + 1.,-0.5);
	bx = ps * H.x;
	by = ps * H.y;
	bz = ps * H.z;
	su = pu + pv * bz - pw * by;
	sv = pv + pw * bx - pu * bz;
	sw = pw + pu * by - pv * bx;
	s1 = bx * bx;
	s2 = by * by;
	s3 = bz * bz;
	s4 = bx * by;
	s5 = by * bz;
	s6 = bz * bx;
	s = s1 + 1. + s2 + s3;
	
	sx = tau1*E.x;
	sy = tau1*E.y;
	sz = tau1*E.z;

	pu1 = ((s1 + 1.) * su + (s4 + bz) * sv + (s6 - by) * sw) / s;
	pv1 = ((s4 - bz) * su + (s2 + 1.) * sv + (s5 + bx) * sw) / s;
	pw1 = ((s6 + by) * su + (s5 - bx) * sv + (s3 + 1.) * sw) / s;
	
	pu = pu1 + sx;
	pv = pv1 + sy;
	pw = pw1 + sz;
	ps = pu * pu + pv * pv + pw * pw;
	ps = pow((pu * pu + pv * pv + pw * pw) + 1.0,-0.5);
	
	u = ps * pu;
	v = ps * pv;
	w = ps * pw;
	x1 = x + tau * u;
	y1 = y + tau * v;
	z1 = z + tau * w;
	
	x = x1;
	y = y1;
	z = z1;
	
}
   
   
   virtual void Collide(double sect){}
   
   double3 GetX(){double3 d3x; d3x.x = x; d3x.y = y; d3x.z = z; return d3x;}
   double3 GetV(){double3 d3x; d3x.x = x; d3x.y = y; d3x.z = z; return d3x;}
   void    SetX(double3 x1){x = x1.x;y = x1.y;z = x1.z;}
   double  GetMass(){return m;}
   double  GetQ2M(){return q_m;}
   
#ifdef DEBUG_PLASMA
   double3 GetXnext(){return next_x;}
   void    SetXnext(double3 x1){next_x.x = x1.x;next_x.y = x1.y;next_x.z = x1.z;}
   int checkParticle(){return (
                               (fabs(x -next_x.x) < TOLERANCE) &&
                               (fabs(y -next_x.y) < TOLERANCE) &&
                               (fabs(z -next_x.z) < TOLERANCE) 
				 
			      );
     
                      }

#endif
   
void Print(FILE* f, int num)
{
     char num_str[20]; 
     sprintf(num_str,"num %05d",num);
     
     string print_str = num_str;
     
     print_str += " x "    + FortranExpWrite(x);
     print_str += " y "    + FortranExpWrite(y);
     print_str += " z "    + FortranExpWrite(z); 
     print_str += " px "   + FortranExpWrite(pu);
     print_str += " py "   + FortranExpWrite(pv);
     print_str += " pz "   + FortranExpWrite(pw);
     print_str += " mass " + FortranExpWrite(m);
     print_str += " q/m"   + FortranExpWrite(q_m);

     fprintf(f,"%s\n",print_str.c_str());
   }
};

#endif
