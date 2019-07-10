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
 \brief  This program reads time-of-flight (ToF) sinogram as STIR interfile and then saves it as a non ToF sinogram such that the ToF sinogram can be in extracted as non ToF 
  sinograms each for every ToF bin. This is further helpful to read the nonToF sinogram and scale the scatter estimate using a factor of counts(ToF bin i)/total counts(in ToF sino)
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

#define NUMARG 30

int main(int argc,char **argv)
{
    using namespace stir;

    static const char * const options[]={
        "argv[1]  output_filename_tof0\n" 
	"argv[2]  output_filename_tof1\n"
        "argv[3]  output_filename_tof2\n"
	"argv[4]  output_filename_tof3\n"
	"argv[5]  output_filename_tof4\n"
	"argv[6]  output_filename_tof5\n"
	"argv[7]  output_filename_tof6\n"
	"argv[8]  output_filename_tof7\n"
	"argv[9]  output_filename_tof8\n"
	"argv[10] output_filename_tof9\n"
	"argv[11] output_filename_tof10\n"
	"argv[12] output_filename_tof11\n"
	"argv[13] output_filename_tof12\n"
	"argv[14] output_filename_tof13\n"
	"argv[15] output_filename_tof-1\n"
	"argv[16] output_filename_tof-2\n"
        "argv[17] output_filename_tof-3\n"
        "argv[18] output_filename_tof-4\n"
        "argv[19] output_filename_tof-5\n"
        "argv[20] output_filename_tof-6\n"
        "argv[21] output_filename_tof-7\n"
        "argv[22] output_filename_tof-8\n"
        "argv[23] output_filename_tof-9\n"
        "argv[24] output_filename_tof-10\n"
        "argv[25] output_filename_tof-11\n"
        "argv[26] output_filename_tof-12\n"
        "argv[27] output_filename_tof-13\n"
	"argv[28] ToF_projdata\n"
        "argv[29] template_non_ToF_projdata"
    };
    if (argc!=NUMARG){
        std::cerr << "\n\nConvert the ToF sinogram into its 27 nonToF counterparts sinogram.\n\n";
        std::cerr << "Not enough arguments !!! ..\n";
        for (int i=1;i<NUMARG;i++) std::cerr << options[i-1];
        exit(EXIT_FAILURE);
    }

    const std::string output_filename_tof0(argv[1]);
    const std::string output_filename_tof1(argv[2]);
    const std::string output_filename_tof2(argv[3]);
    const std::string output_filename_tof3(argv[4]);
    const std::string output_filename_tof4(argv[5]);
    const std::string output_filename_tof5(argv[6]);
    const std::string output_filename_tof6(argv[7]);
    const std::string output_filename_tof7(argv[8]);
    const std::string output_filename_tof8(argv[9]);
    const std::string output_filename_tof9(argv[10]);
    const std::string output_filename_tof10(argv[11]);
    const std::string output_filename_tof11(argv[12]);
    const std::string output_filename_tof12(argv[13]);
    const std::string output_filename_tof13(argv[14]);
    const std::string output_filename_tofneg1(argv[15]);
    const std::string output_filename_tofneg2(argv[16]);
    const std::string output_filename_tofneg3(argv[17]);
    const std::string output_filename_tofneg4(argv[18]);
    const std::string output_filename_tofneg5(argv[19]);
    const std::string output_filename_tofneg6(argv[20]);
    const std::string output_filename_tofneg7(argv[21]);
    const std::string output_filename_tofneg8(argv[22]);
    const std::string output_filename_tofneg9(argv[23]);
    const std::string output_filename_tofneg10(argv[24]);
    const std::string output_filename_tofneg11(argv[25]);
    const std::string output_filename_tofneg12(argv[26]);
    const std::string output_filename_tofneg13(argv[27]);
    const std::string ToF_projdata(argv[28]);
    shared_ptr<ProjDataInfo> template_projdata_info_sptr =
            ProjData::read_from_file(argv[29])->get_proj_data_info_sptr();

    shared_ptr<ProjData> template_projdata_ptr =
            ProjData::read_from_file(argv[29]);

    ProjDataInterfile proj_data_tof0(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof0, std::ios::out);
    ProjDataInterfile proj_data_tof1(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof1, std::ios::out);
    ProjDataInterfile proj_data_tof2(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof2, std::ios::out);
    ProjDataInterfile proj_data_tof3(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof3, std::ios::out);
    ProjDataInterfile proj_data_tof4(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof4, std::ios::out);
    ProjDataInterfile proj_data_tof5(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof5, std::ios::out);
    ProjDataInterfile proj_data_tof6(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof6, std::ios::out);
    ProjDataInterfile proj_data_tof7(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof7, std::ios::out);
    ProjDataInterfile proj_data_tof8(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof8, std::ios::out);
    ProjDataInterfile proj_data_tof9(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof9, std::ios::out);
    ProjDataInterfile proj_data_tof10(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof10, std::ios::out);
    ProjDataInterfile proj_data_tof11(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof11, std::ios::out);
    ProjDataInterfile proj_data_tof12(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof12, std::ios::out);
    ProjDataInterfile proj_data_tof13(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tof13, std::ios::out);
    ProjDataInterfile proj_data_tofneg1(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg1, std::ios::out);
    ProjDataInterfile proj_data_tofneg2(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg2, std::ios::out);
    ProjDataInterfile proj_data_tofneg3(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg3, std::ios::out);
    ProjDataInterfile proj_data_tofneg4(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg4, std::ios::out);
    ProjDataInterfile proj_data_tofneg5(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg5, std::ios::out);
    ProjDataInterfile proj_data_tofneg6(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg6, std::ios::out);
    ProjDataInterfile proj_data_tofneg7(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg7, std::ios::out);
    ProjDataInterfile proj_data_tofneg8(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg8, std::ios::out);
    ProjDataInterfile proj_data_tofneg9(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg9, std::ios::out);
    ProjDataInterfile proj_data_tofneg10(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg10, std::ios::out);
    ProjDataInterfile proj_data_tofneg11(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg11, std::ios::out);
    ProjDataInterfile proj_data_tofneg12(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg12, std::ios::out);
    ProjDataInterfile proj_data_tofneg13(template_projdata_ptr->get_exam_info_sptr(), template_projdata_ptr->get_proj_data_info_sptr(),
                                output_filename_tofneg13, std::ios::out);

    shared_ptr<ProjDataInfo> proj_data_tof_info_sptr = ProjData::read_from_file(argv[28])->get_proj_data_info_sptr();
    shared_ptr<ProjData> proj_data_tof_ptr = ProjData::read_from_file(argv[28]);

    Bin bin;

    for (bin.segment_num() = proj_data_tof0.get_min_segment_num();
     bin.segment_num() <= proj_data_tof0.get_max_segment_num();
     ++ bin.segment_num())
      {

    for (bin.axial_pos_num() = proj_data_tof0.get_min_axial_pos_num(bin.segment_num());
         bin.axial_pos_num() <= proj_data_tof0.get_max_axial_pos_num(bin.segment_num());
         ++bin.axial_pos_num())
      {
       // for (bin.timing_pos_num() = proj_data.get_min_tof_pos_num();
        //bin.timing_pos_num() <= proj_data.get_max_tof_pos_num();
      //  ++ bin.timing_pos_num())
       // {
        Sinogram<float> non_ToF_sinogram_tof0 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
	Sinogram<float> non_ToF_sinogram_tof1 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof2 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof3 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
	Sinogram<float> non_ToF_sinogram_tof4 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof5 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof6 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof7 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
	Sinogram<float> non_ToF_sinogram_tof8 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof9 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof10 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof11 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof12 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tof13 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg1 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg2 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
	Sinogram<float> non_ToF_sinogram_tofneg3 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg4 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg5 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg6 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg7 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg8 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg9 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg10 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg11 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg12 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
        Sinogram<float> non_ToF_sinogram_tofneg13 =
          template_projdata_info_sptr->get_empty_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);

	Sinogram<float> sinogram_tof0 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,0);
	Sinogram<float> sinogram_tof1 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,1);
        Sinogram<float> sinogram_tof2 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,2);
        Sinogram<float> sinogram_tof3 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,3);
	Sinogram<float> sinogram_tof4 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,4);
        Sinogram<float> sinogram_tof5 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,5);
        Sinogram<float> sinogram_tof6 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,6);
        Sinogram<float> sinogram_tof7 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,7);
	Sinogram<float> sinogram_tof8 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,8);
        Sinogram<float> sinogram_tof9 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,9);
        Sinogram<float> sinogram_tof10 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,10);
        Sinogram<float> sinogram_tof11 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,11);
        Sinogram<float> sinogram_tof12 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,12);
        Sinogram<float> sinogram_tof13 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,13);
        Sinogram<float> sinogram_tofneg1 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-1);
        Sinogram<float> sinogram_tofneg2 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-2);
	Sinogram<float> sinogram_tofneg3 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-3);
        Sinogram<float> sinogram_tofneg4 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-4);
        Sinogram<float> sinogram_tofneg5 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-5);
        Sinogram<float> sinogram_tofneg6 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-6);
        Sinogram<float> sinogram_tofneg7 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-7);
        Sinogram<float> sinogram_tofneg8 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-8);
        Sinogram<float> sinogram_tofneg9 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-9);
        Sinogram<float> sinogram_tofneg10 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-10);
        Sinogram<float> sinogram_tofneg11 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-11);
        Sinogram<float> sinogram_tofneg12 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-12);
        Sinogram<float> sinogram_tofneg13 = proj_data_tof_ptr->get_sinogram(bin.axial_pos_num(),bin.segment_num(),false,-13);
   
	     for (bin.view_num() = proj_data_tof0.get_min_view_num();
             bin.view_num() <= proj_data_tof0.get_max_view_num();
             ++ bin.view_num())
              {

            for (bin.tangential_pos_num() = proj_data_tof0.get_min_tangential_pos_num();
                 bin.tangential_pos_num() <= proj_data_tof0.get_max_tangential_pos_num();
                 ++bin.tangential_pos_num())
              {
                /*(*segment_ptr)[bin.axial_pos_num()]*/
                non_ToF_sinogram_tof0[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof0[bin.view_num()][bin.tangential_pos_num()];
		non_ToF_sinogram_tof1[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof1[bin.view_num()][bin.tangential_pos_num()];
 		non_ToF_sinogram_tof2[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof2[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof3[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof3[bin.view_num()][bin.tangential_pos_num()];
		non_ToF_sinogram_tof4[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof4[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof5[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof5[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof6[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof6[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof7[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof7[bin.view_num()][bin.tangential_pos_num()];
		non_ToF_sinogram_tof8[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof8[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof9[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof9[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof10[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof10[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof11[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof11[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof12[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof12[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tof13[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tof13[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg1[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg1[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg2[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg2[bin.view_num()][bin.tangential_pos_num()];
		non_ToF_sinogram_tofneg3[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg3[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg4[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg4[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg5[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg5[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg6[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg6[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg7[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg7[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg8[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg8[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg9[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg9[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg10[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg10[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg11[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg11[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg12[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg12[bin.view_num()][bin.tangential_pos_num()];
                non_ToF_sinogram_tofneg13[bin.view_num()][bin.tangential_pos_num()] +=
                sinogram_tofneg13[bin.view_num()][bin.tangential_pos_num()];
                }
              }



        proj_data_tof0.set_sinogram(non_ToF_sinogram_tof0);
	proj_data_tof1.set_sinogram(non_ToF_sinogram_tof1);
	proj_data_tof2.set_sinogram(non_ToF_sinogram_tof2);
        proj_data_tof3.set_sinogram(non_ToF_sinogram_tof3);
	proj_data_tof4.set_sinogram(non_ToF_sinogram_tof4);
        proj_data_tof5.set_sinogram(non_ToF_sinogram_tof5);
        proj_data_tof6.set_sinogram(non_ToF_sinogram_tof6);
        proj_data_tof7.set_sinogram(non_ToF_sinogram_tof7);
	proj_data_tof8.set_sinogram(non_ToF_sinogram_tof8);
        proj_data_tof9.set_sinogram(non_ToF_sinogram_tof9);
        proj_data_tof10.set_sinogram(non_ToF_sinogram_tof10);
        proj_data_tof11.set_sinogram(non_ToF_sinogram_tof11);
        proj_data_tof12.set_sinogram(non_ToF_sinogram_tof12);
        proj_data_tof13.set_sinogram(non_ToF_sinogram_tof13);
        proj_data_tofneg1.set_sinogram(non_ToF_sinogram_tofneg1);
        proj_data_tofneg2.set_sinogram(non_ToF_sinogram_tofneg2);
	proj_data_tofneg3.set_sinogram(non_ToF_sinogram_tofneg3);
        proj_data_tofneg4.set_sinogram(non_ToF_sinogram_tofneg4);
        proj_data_tofneg5.set_sinogram(non_ToF_sinogram_tofneg5);
        proj_data_tofneg6.set_sinogram(non_ToF_sinogram_tofneg6);
        proj_data_tofneg7.set_sinogram(non_ToF_sinogram_tofneg7);
        proj_data_tofneg8.set_sinogram(non_ToF_sinogram_tofneg8);
        proj_data_tofneg9.set_sinogram(non_ToF_sinogram_tofneg9);
        proj_data_tofneg10.set_sinogram(non_ToF_sinogram_tofneg10);
        proj_data_tofneg11.set_sinogram(non_ToF_sinogram_tofneg11);
        proj_data_tofneg12.set_sinogram(non_ToF_sinogram_tofneg12);
        proj_data_tofneg13.set_sinogram(non_ToF_sinogram_tofneg13);
       // }
      }
        }

    return EXIT_SUCCESS;
}
