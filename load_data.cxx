/*
 * read_particle.cxx
 *
 *  Created on: Jun 9, 2018
 *      Author: snytav
 */


#include <string>


#include "load_data.h"

std::string getBinaryFileName(int nt)
	  {
		  char part_name[100];
		  std::string s;

		  sprintf(part_name,"mumu000%08d.dat",nt);

		  s = part_name;

		  return s;
	  }

