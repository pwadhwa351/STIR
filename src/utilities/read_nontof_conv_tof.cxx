/*
 Copyright (C) 2013, University College London
 Copyright (C) 2018, University of Leeds
 This file is part of STIR.
 This file is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 This file is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 See STIR/LICENSE.txt for details
 */
/*!
 \file
 \ingroup utilities
 \brief  This program reads non time-of-flight (ToF) sinogram as STIR interfile and then saves it as a ToF sinogram such that the data in each ToF bin is i
 same and equal to non-ToF-data/num_ToF_bins (such that adding all ToF_bins is equal to the non-ToF sinogram).
 as STIR interfile.
 \author Palak Wadhwa
 */

#include "stir/ProjDataInterfile.h"
#include "stir/ExamInfo.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/ProjDataFromHDF5.h"
#include "stir/IO/read_data.h"
#include "stir/Succeeded.h"
#include "stir/NumericType.h"
#include "stir/IO/HDF5Wrapper.h"
#include "stir/IndexRange3D.h"
#include "stir/DetectionPosition.h"
#include "stir/DetectionPositionPair.h"
#include "stir/shared_ptr.h"
#include "stir/RelatedViewgrams.h"
#include "stir/ViewSegmentNumbers.h"
#include "stir/IndexRange2D.h"
#include "stir/IndexRange.h"
#include "stir/Bin.h"
#include "stir/display.h"
#include "stir/IO/read_data.h"
#include "stir/IO/InterfileHeader.h"
#include "stir/ByteOrder.h"
#include "stir/is_null_ptr.h"
#include "stir/modulo.h"
#include <algorithm>
#include <fstream>
#include <cctype>
#ifndef STIR_NO_NAMESPACES
using std::ofstream;
using std::fstream;
using std::ios;
#endif

#define NUMARG 4

int main(int argc,char **argv)
{
    using namespace stir;

    static const char * const options[]={
        "argv[1]  output_filename\n"
        "argv[2]  non_ToF_projdata\n"
        "argv[3]  template_ToF_projdata"
    };
    if (argc!=NUMARG){
        std::cerr << "\n\nConvert the non-ToF sinogram into ToF sinogram such that each ToF bin is same and equal to non-ToF-data/num_ToF_bins.\n\n";
        std::cerr << "Not enough arguments !!! ..\n";
        for (int i=1;i<NUMARG;i++) std::cerr << options[i-1];
        exit(EXIT_FAILURE);
    }

    const std::string output_filename(argv[1]);
    const std::string non_ToF_projdata(argv[2]);
    shared_ptr<ProjDataInfo> template_projdata_info_sptr =
            ProjData::read_from_file(argv[3])->get_proj_data_info_sptr();

    shared_ptr<ProjData> template_projdata_ptr =
            ProjData::read_from_file(argv[3]);

    ProjDataInterfile proj_data(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename, std::ios::out);
    shared_ptr<ProjDataInfo> proj_data_non_tof_info_sptr = ProjData::read_from_file(argv[2]);
    shared_ptr<ProjData> proj_data_non_tof_ptr = ProjData::read_from_file(argv[2]);

    Bin bin;

    for (bin.segment_num() = proj_data.get_min_segment_num();
     bin.segment_num() <= proj_data.get_max_segment_num();
     ++ bin.segment_num())
      {

    for (bin.axial_pos_num() = proj_data.get_min_axial_pos_num(bin.segment_num());
         bin.axial_pos_num() <= proj_data.get_max_axial_pos_num(bin.segment_num());
         ++bin.axial_pos_num())
      {
        for (bin.timing_pos_num() = proj_data.get_min_tof_pos_num();
        bin.timing_pos_num() <= proj_data.get_max_tof_pos_num();
        ++ bin.timing_pos_num())
        {
        Sinogram<float> sinogram =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,bin.timing_pos_num());
        Sinogram<float> non_ToF_sinogram = proj_data_non_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
            for (bin.view_num() = proj_data.get_min_view_num();
             bin.view_num() <= proj_data.get_max_view_num();
             ++ bin.view_num())
              {

            for (bin.tangential_pos_num() = proj_data.get_min_tangential_pos_num();
                 bin.tangential_pos_num() <= proj_data.get_max_tangential_pos_num();
                 ++bin.tangential_pos_num())
              {
                /*(*segment_ptr)[bin.axial_pos_num()]*/
                sinogram[bin.view_num()][bin.tangential_pos_num()] =
                non_ToF_sinogram[bin.view_num()][bin.tangential_pos_num()]/template_projdata_info_sptr->get_max_tof_pos_num();
                }
              }



        proj_data.set_sinogram(sinogram);
        }
      }
        }

    return EXIT_SUCCESS;
}
