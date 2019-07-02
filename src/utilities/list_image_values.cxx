/*
    Copyright (C) 2002- 2011, Hammersmith Imanet Ltd
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the Lesser GNU General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    Lesser GNU General Public License for more details.

    See STIR/LICENSE.txt for details
*/
/*!
  \file
  \ingroup utilities
  \brief Utility to extract image values along profiles (or another box-shape)
  \author Sanida Mustafovic
  \author Kris Thielemans

  \par Usage
  \code
   list_image_values output_file_name input_image \\
       min_plane max_plane  min_col max_col min_row max_row
  \endcode
  Indices need to be in the STIR convention (plane starts from 0, col,row are centred around 0)

*/

#include "stir/DiscretisedDensity.h"
#include "stir/IO/read_from_file.h"
#include <fstream>
#include <iomanip>



USING_NAMESPACE_STIR

int main(int argc, char *argv[])
{ 
  
  if (argc!= 9)
  {
    std::cerr << "Usage:\n" << argv[0] << " output_profile_name input_image min_plane max_plane  min_col max_col min_row max_row\n";    
    return(EXIT_FAILURE);
  }

  const char * const output_filename = argv[1];
  const char * const input_filename = argv[2];

    
  shared_ptr<DiscretisedDensity<3,float> > 
    input_image_sptr(read_from_file<DiscretisedDensity<3,float> >(input_filename));
  const DiscretisedDensity<3,float>& input_image =
    *input_image_sptr;
  
  const int min_plane_num = atoi(argv[3]);
  const int max_plane_num = atoi(argv[4]);    
  const int min_column_num = atoi(argv[5]);
  const int max_column_num = atoi(argv[6]);    
  const int min_row_num   = atoi(argv[7]);
  const int max_row_num   = atoi(argv[8]);  
    
  std::ofstream  profile_file(output_filename);

  using std::setw;
  profile_file << setw(8) << "plane" << setw(8) << "row" 
	       << setw(8) << "column" << setw(10) << "value" <<'\n';
   
  for (int plane = min_plane_num;plane <=max_plane_num;plane++) 
    for (int row = min_row_num;row <=max_row_num;row++)
      for (int column = min_column_num;column<=max_column_num;column++)
      {
	profile_file << setw(8) << plane << setw(8) << row 
		     << setw(8) << column
		     << setw(10) << input_image[plane][row][column] << '\n';
      }
  
  return EXIT_SUCCESS;
}

